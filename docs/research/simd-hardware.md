# flashvec · SIMD 与硬件利用率调研（Zen 2 / Ryzen 7 5700U）

调研日期 2026-10-17 · 目标机器 AMD Ryzen 7 5700U（Lucienne, Zen 2, 8C/16T）· g++ 12.4 · Linux 6.17
本机实验源码：`/tmp/fvsimd/{cpuinfo.c,bench3.cpp,bench4.cpp,bench5.cpp,asmcheck.cpp,dispatch.cpp}`（临时文件，未入库）
所有「实测」数字都在本机跑出；所有链接均已 HTTP 200 核验，未核验的标「未验证」。

## 0 · 本机事实（命令与输出，原文粘贴）

    $ grep -m1 'model name' /proc/cpuinfo
    model name    : AMD Ryzen 7 5700U with Radeon Graphics
    $ lscpu -C
    NAME ONE-SIZE ALL-SIZE WAYS TYPE   LEVEL SETS PHY-LINE COHERENCY-SIZE
    L1d       32K     256K    8 Data      1   64        1             64
    L2       512K       4M    8 Unified   2 1024        1             64
    L3         4M       8M   16 Unified   3 4096        1             64
    $ cc -O2 cpuinfo.c -o cpuinfo && ./cpuinfo
    arch=amd
    sse2=16 sse4.2=256 avx=512 avx2=1024 fma=16384 f16c=131072 bmi2=131072
    avx512f=0 avx512f_var=0
    $ perf stat -e cycles,instructions true
    access to performance monitoring and observability operations ...
    perf_event_paranoid setting is 4        # 无 CAP_PERFMON；sudo -n 需要密码

三条硬结论：① 本机 **无 AVX-512**（`/proc/cpuinfo` 里 `grep -c avx512` = 0，`__builtin_cpu_supports("avx512f")` = 0）；
② `__builtin_cpu_supports` 返回**位掩码而非 0/1**（见上：sse2=16），判断必须写 `!= 0`；
③ `-march=native` 在本机解析为 **znver2**（`gcc -march=native -Q --help=target`）；`perf stat` 被 `perf_event_paranoid=4` 挡住，
所以第 5 节的周期数一律由 rdtsc + cpufreq 采样换算，不来自性能计数器。

## 1 · 结论速览

1. **当前 API 形态让自动向量化结构性失效**：循环体是一次不可内联的 `extern "C"` 调用时，GCC 直接拒绝向量化
   （`missed: statement clobbers memory: *_5 = fv_vec2_add(...)`）。（对 flashvec：先解决「不可内联」，再谈 intrinsics。）
2. **按值 ABI 本身几乎免费**：`fv_vec2`（8B）走 1 个 XMM 寄存器，`fv_vec3`（12B）走 2 个；库里 `fv_vec2_add` 就是
   `addps %xmm1,%xmm0; ret`。（别动结构体布局和按值签名，动「调用边界」。）
3. **贵的是跨 TU 调用边界**：同一个 vec2 加法，库调用 1.66 ns/元素、内联 0.57、256-bit 打包 0.21（同一次运行内）。
   边界 ≈ 4.2 cycle/元素，是「打包收益」（≈1.3 cycle/元素）的 **3 倍**。（先批量化 API；否则 SIMD 收益被调用开销吃掉。）
4. **vec2 装进 `__m256` 不需要任何数据重排**：add/sub/scale/lerp 都是逐分量运算，`_mm256_load_ps` 直接把连续 4 个
   `fv_vec2` 当 4 个元素用；只有 dot/len/normalize 需要 4→1 归约。（对 flashvec：`_n` 批量接口先覆盖 add/sub/scale/lerp。）
5. **水平归约要避开 `haddps`**：打包 dot 用 `hadd` 是 0.374 ns/元素，用 `permute_ps+add_ps+shuffle_ps` 只要 0.217 ns
   （1.7×）。GLM 的 `_mm_hadd_ps` 写法（`glm/simd/geometric.h#L33`）在 Zen 2 上不划算。
6. **要打满 FMA 需要 ≥10 条独立链**（延迟 5 × 吞吐 2）：8 个 ymm 累加器实测只有 1.4 FMA/周期，12 个才到 ~2.0。
   （对 flashvec：批量块至少 12×4 = 48 个 vec2，或者干脆把块做得比 L1 小得多。）
7. **`-fno-math-errno` 是性价比最高的一步**：默认 `sqrtf` 带 errno 语义 → 生成分支 + 库调用（`fv_vec2_len` 14 条指令、
   `fv_vec3_normalize` **48 条**）；加该开关后分别降到 10 / **22** 条，且让含 sqrt 的循环可向量化（实测
   `missed: ... sqrtf` → `optimized: loop vectorized using 32 byte vectors`）。（只加在库的编译开关上，不改头文件语义。）
8. **不要全局 `-march=native`，用 `__attribute__((target("avx2,fma")))` 单函数开**：实测不加 `-mavx2` 时 intrinsics
   直接编译失败（`target specific option mismatch`），而 target 属性在 `-O2` 基线下可编译；派发用
   `__builtin_cpu_supports`（实测编译成 1 次内存 load + test + jmp，无 CPUID）。

## 2 · 参考库解剖

| 库 | 怎么用 SIMD | 可借鉴 | 代价 |
|---|---|---|---|
| GLM | 模板 `vec<L,T,Q>` 用宏切 ABI/对齐；显式 SIMD 只在 load/store 与几何函数（`_mm_dp_ps`/`_mm_hadd_ps`），vec4 四则运算源码是标量，靠编译器 | 「SIMD 只用在值得的地方」；4 分量 padding 让 vec3 也能塞进一个 XMM | 头文件重、编译慢；按值 ABI 随平台变；vec4 强绑定 |
| Eigen | 每个 arch 一张表定义 packet（`Packet4f = __m128`），packet 之上再套表达式模板；自带对齐分配器 | 「先定义 packet + traits，再在其上写算法」的分层，是 flashvec 扩到 mat4/quat 后最该抄的骨架 | 表达式模板=编译期爆炸、报错难读；对齐/别名要求传导到用户结构体 |
| DirectXMath | `XMVECTOR = __m128`，调用约定 `XM_CALLCONV __vectorcall`，强制 16B 对齐 | 「向量类型就是一个寄存器」的接口哲学；`XM_ALIGNED_DATA` 的显式对齐 | `__vectorcall` 是 MSVC ABI，GCC/Linux 忽略它；强制 4 分量浪费 25% 带宽；C++ only |
| MathFu | `Vector<T,N>` 模板 + SIMD-aware new/delete，用 4 分量内部表示承载 SIMD | 「按类型给对齐/分配器」的思路 | 依赖分配器与 C++ 模板，与 flashvec「无分配 + 纯 C」契约冲突 |
| ultraviolet | 用多一个 w 分量把 vec2/vec3 撑到 16B 对齐（**未验证**：`github.com/kwgoodman/ultraviolet` 与 jsdelivr 均 404，npm 同名包是无关的 0.0.1 占位包，本机无法取得原文） | 若属实，「padding 到对齐」是 AoS 打包最省事的手法 | 内存 +25%~100%；无来源可引用，不建议照搬 |
| google/highway | 可移植 batch + 编译期 target 表 + `HWY_EXPORT`/`HWY_DYNAMIC_DISPATCH` 运行期派发；用 `__attribute__((target(...)))` 而不是全局 `-march` | 「一份源码 N 个 target」的派发机制，正是 flashvec 若要同时上 AVX2/NEON 的正确答案 | 引入重依赖、宏与模板复杂度极高、编译时间显著上升 |
| xsimd | `batch<T,N>` 类型 + 对齐类型/分配器，load/store 断言检查对齐 | batch 的对齐契约（`requires_aligned_allocator`）与 load/store API | C++ only；收益依赖编译器把 batch 正确映射到寄存器 |

关键文件：GLM [simd/platform.h](https://github.com/g-truc/glm/blob/master/glm/simd/platform.h)、
[detail/type_vec4.inl#L1018](https://github.com/g-truc/glm/blob/master/glm/detail/type_vec4.inl#L1018)（`#if GLM_ARCH & GLM_ARCH_SSE2_BIT` 门控块，`grep -c _mm_` 全文只有 8 处）、
[simd/geometric.h#L30](https://github.com/g-truc/glm/blob/master/glm/simd/geometric.h#L30)（`_mm_dp_ps`，L33 `_mm_hadd_ps`）；
Eigen [arch/SSE/PacketMath.h#L46](https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/Core/arch/SSE/PacketMath.h#L46)（`using Packet4f = __m128;`）；
DirectXMath [Inc/DirectXMath.h#L28](https://github.com/microsoft/DirectXMath/blob/main/Inc/DirectXMath.h#L28)（`#define XM_CALLCONV __vectorcall`）、[#L402](https://github.com/microsoft/DirectXMath/blob/main/Inc/DirectXMath.h#L402)（`using XMVECTOR = __m128;`）、
[Inc/DirectXMathVector.inl](https://github.com/microsoft/DirectXMath/blob/main/Inc/DirectXMathVector.inl)；
MathFu [include/mathfu/vector.h#L529](https://github.com/google/mathfu/blob/master/include/mathfu/vector.h#L529)（`MATHFU_DEFINE_CLASS_SIMD_AWARE_NEW_DELETE`）；
highway [hwy/highway.h#L368](https://github.com/google/highway/blob/master/hwy/highway.h#L368)（`HWY_DYNAMIC_DISPATCH`）、[hwy/base.h#L222](https://github.com/google/highway/blob/master/hwy/base.h#L222)（`__attribute__((target(targets_str)))`）、[hwy/ops/x86_256-inl.h](https://github.com/google/highway/blob/master/hwy/ops/x86_256-inl.h)；
xsimd [include/xsimd/types/xsimd_batch.hpp#L215](https://github.com/xtensor-stack/xsimd/blob/master/include/xsimd/types/xsimd_batch.hpp#L215)（`class batch`）。

## 3 · 编译器自动向量化的真实条件

**按值聚合体的 ABI**（本机反汇编实测，`-O2`）：

    <fv_vec2_add>:  addps %xmm1,%xmm0 ; ret                    # 8B 结构体 = 1 个 XMM
    <fv_vec3_add>:  addss %xmm3,%xmm1 ; addps %xmm2,%xmm0; ret  # 12B = 2 个 XMM，第二个只用低 4 字节
    call site:      vmovq 0x0(%rbp),%xmm1 ; vmovq (%r12),%xmm0 ; call fv_vec2_add

即按值传参**不产生栈流量**（符合 SysV x86-64 [psABI](https://gitlab.com/x86-psABIs/x86-64-ABI) 的 eightbyte 分类）。
真正的杀手是**调用本身**：`-fopt-info-vec-missed` 给出

    asmcheck.cpp:16:43: missed: statement clobbers memory: *_5 = fv_vec2_add (MEM[...], MEM[...]);
    asmcheck.cpp:19:52: missed: statement clobbers memory: _14 = fv_vec2_dot (MEM[...], MEM[...]);

不透明调用被当成「可能写内存」，循环直接放弃向量化 —— 这三条是 flashvec 现状的数学底座：
**想向量化，必须让编译器看见循环体**（头文件内联、LTO、或库内自带 `_n` 循环）。

**别名分析 / `__restrict`**：同一份内联代码，GCC 的报告是
`optimized: loop vectorized using 32 byte vectors` 外加 `loop versioned for vectorization because of possible aliasing`
（若指针重叠，运行期退回标量版本）；不加 restrict 的 `sum_copy` 则没有任何向量化记录。→ 批量接口请一律写
`fv_vec2 * __restrict`。

**-O3 与 -O2**：GCC 12 在 `-O2 -march=native` 下只报了 2 条向量化记录，`-O3 -march=native` 下同一文件报了 12 条
（-O2 用 very-cheap 代价模型，大量放弃）。→ 库若要走自动向量化，至少要 `-O3`。

**-ffast-math / -funsafe-math-optimizations 的语义风险**：会一并打开 `-fno-math-errno`、`-fno-trapping-math`、
`-ffinite-math-only`（假定无 NaN/Inf），允许重结合、向量化归约；物理引擎若依赖 NaN 传播做调试、或用
`x != x` 判 NaN，会静默失效。**廉价子集是 `-fno-math-errno`**：只承诺「不设置 errno」，实测已足以让
`1.0f/sqrtf(d)` 的循环从 `missed: statement clobbers memory: _32 = sqrtf (d_23)` 变成
`optimized: loop vectorized using 32 byte vectors`，并把 `fv_vec3_normalize` 从 48 条指令压到 22 条。
（flashvec 里 sqrt 的参数恒为平方和 ≥0，除了 errno 之外语义零变化。）见 [GCC Optimize Options](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)。

**-march=native 的可移植性代价**：二进制绑定 znver2，换机器 SIGILL；且它会把整 TU 的基线抬高。实测在默认
`-O2`（无 `-mavx2`）下写 intrinsics **直接编译失败**：

    error: inlining failed in call to 'always_inline' '__m256 _mm256_add_ps(...)': target specific option mismatch

正解见第 4 节（单函数 target 属性）。

**`__builtin_assume_aligned`**：给编译器对齐承诺后可消除版本化/对齐检查（本机 `asum_aligned` 报 32B 向量化）。
代价：它不做运行期检查，承诺不成立就是 UB/段错误，只能用在由 API 契约保证对齐的参数上。

**`always_inline`**：`static inline __attribute__((always_inline))` 是让「库调用」消失的最便宜手段，但要求实现进头文件，
与 flashvec 当前「纯 C 头 + 独立 .cpp 实现」冲突。折中：**加 `_n` 批量接口，循环放进 .cpp** ——调用边界从
「每元素一次」降到「每数组一次」，不需要动头文件的 ABI。

**LTO**：实测有效。把 `asmcheck.cpp` 与 `src/flashvec.cpp` 用 `-flto -O3 -march=native` 一起链接后，
原本「clobbers memory」的循环变成：

    asmcheck.cpp:16:19: optimized: loop vectorized using 32 byte vectors   # call_vec2_add
    asmcheck.cpp:19:30: optimized: loop vectorized using 32 byte vectors   # call_vec2_dot

即 **LTO 能救回现有 API**，代价是构建链路必须全程带 `-flto`（预编译的静态库交出去就失效）、链接时间与调试复杂度上升。

## 4 · 手写 intrinsics 的工程细节

- **头文件 guard**：intrinsics 依赖 `-mavx2`/target 属性，别写进纯 C 的 `flashvec.h`（该头必须保持 `-x c` 可编译）。
  建议新增 `include/flashvec_simd.h`（C++/内部用）并把 `#include <immintrin.h>` 放在 `.cpp` 里，用
  `#if defined(__AVX2__) || defined(FLASHVEC_TARGET_AVX2)` 之类的自有宏门控。
- **运行时派发三种做法**（均在本机验证）：
  1. `__builtin_cpu_supports("avx2")` 分支 —— 反汇编是 `mov __cpu_model(%rip),%rax` + `testb $0x4,0xd(%rax)` + `jmp`，
     即「一次全局 load + test」，**不是**每次 CPUID，开销可忽略；
  2. `__attribute__((target("avx2,fma")))` 单函数 —— `-O2` 基线下即可编译，未选中的机器永不执行它；
  3. `__attribute__((target_clones("avx2","default")))` —— 生成 `scale8.avx2`/`scale8.default` + ifunc 符号，调用点全走 PLT，
     首次调用做 resolver；静态链接或非 glibc 平台可能不可用。见 [GCC Function Attributes](https://gcc.gnu.org/onlinedocs/gcc/Common-Function-Attributes.html)。
- **对齐**：`alignas(32)` 只对 `_mm256_load_ps/store_ps`（要求对齐）有意义；`fv_vec2` 只有 4 字节对齐，
  把 `fv_vec2*` 直接当 `__m256*` 用是 UB（虽然 `loadu` 在 L1 命中的代价很小，本机未单独测量）。
  批量接口请把对齐写成契约：`alignas(32) float buf[]` 或全程 `loadu`。
- **水平操作**：`vhaddps` 实测 1.10 ns/次（≈4.4 cycle，@~4.1GHz），比 `vpermilps+vaddps` 组合慢；打包 dot 的两种写法
  同一次运行内：`mul+2*hadd` = 0.374 ns/元素，`mul+permute_ps(p,0xB1)+add_ps+shuffle_ps` = 0.217 ns/元素。
  4 路 vec2 点积用后者，不要用 `hadd`。
- **AoS↔SoA 转换成本**：`vgatherdps` 实测 2.49–2.57 ns/元素（≈10 cycle/条指令，8 lane），比顺序 load 贵一个数量级，
  **不要**用它做 AoS→SoA 重排；`vpermps` 约 4.4 cycle/条。结论：vec2/vec4 这类「一个向量=一条记录」的 AoS，
  逐元素阶段**零转换**直接打包；只有归约阶段值得转置，而转置本身往往比省下的归约更贵。
- **x86 ↔ ARM NEON 可迁移性**：SSE 的 `_mm_*` 有 [sse2neon](https://github.com/DLTcollab/sse2neon)（头文件级翻译）与
  [SIMDe](https://github.com/simd-everywhere/simde)（运行时模拟）可用；AVX2 的 256-bit 在 NEON 上是两组 128-bit 寄存器，
  **不存在一一对应**。若将来真要跨架构，正确做法是像 highway 那样写一层可移植 batch，而不是散落 `_mm256_*`。
  参考 [ARM NEON intrinsics](https://developer.arm.com/architectures/instruction-sets/simd-isas/neon/intrinsics)。

## 5 · Zen 2 硬数字（本机核对 + 文献对照）

**方法学警告（必须先读）**：`perf stat` 在本机不可用（`perf_event_paranoid=4` 且无 sudo），改用 `__rdtsc` + 采样
`scaling_cur_freq`（取中位数）；且本 APU 虽 `/sys/.../boost` 显示 0、governor=performance，实测频率仍在
1.64–3.97 GHz 之间漂移 —— 同一操作在不同程序/不同运行时段的绝对值差过 **2×**（例：`vpermps` 一次测 0.52 ns，
  另一次 1.11 ns）。**因此本文只使用「同一次运行内的比值」**，绝对周期数由两条延迟基准反推有效主频
（`vaddps` 依赖链 0.735 ns/次 = 3 cycle ⇒ ≈4.08 GHz；`vfmadd` 链 1.224 ns = 5 cycle ⇒ ≈4.08 GHz，二者自洽）。

**(a) 延迟与端口（同一次运行，min-of-6，taskset -c 2）**

| 项目 | 实测 | 换算 | 文献对照 |
|---|---|---|---|
| `fv_vec2_add` 库调用 依赖链 | 2.710 ns | ≈10.7 cycle | — |
| 内联 vec2 加法 依赖链 | 0.802 ns | ≈3.2 cycle | — |
| `fv_vec3_add` 库调用 依赖链 | 3.160 ns | ≈12.4 cycle | — |
| `vaddps` ymm 依赖链 | 0.735 ns | **3 cycle**（延迟） | Agner/uops.info：Zen 2 ADDPS 延迟 3 |
| `vfmadd213ps` ymm 依赖链 | 1.224 ns | **5 cycle**（延迟） | Zen 2 FMA 延迟 5 |
| `vaddps` ymm 吞吐（12 累加器） | 0.130 ns/条 | ≈1.9 条/cycle | 2 条/cycle（2 个 256-bit ADD 管） |
| `vfmadd` ymm 吞吐（12 累加器） | 0.137 ns/条 | ≈2.0 条/cycle | 2 条/cycle |
| `vfmadd` ymm 吞吐（**8** 累加器） | 0.170 ns/条 | ≈1.4 条/cycle | 依赖受限：5×2=10 条链才满 |
| `vsqrtps` ymm 吞吐 | 1.354 ns/条 | ≈5.5 cycle/8 lane | ≈0.69 cycle/元素 |
| `vpermps` 吞吐 | 0.52–1.11 ns/条 | 2.2–4.4 cycle | 变长延迟/跨 lane 洗牌贵 |
| `vhaddps` 吞吐 | 1.10 ns/条 | ≈4.4 cycle | 水平加法=多 µop |
| `vgatherdps` 吞吐 | 2.49 ns/条 | ≈10 cycle/8 lane | Zen 2 gather 以「每元素 1 cycle 量级」慢 |

文献：Agner Fog [Instruction tables](https://www.agner.org/optimize/instruction_tables.pdf)、[Microarchitecture](https://www.agner.org/optimize/microarchitecture.pdf)、
[uops.info](https://uops.info/table.html)（本机只能证明这些文档页可达，**未能逐格核对**；§0 的缓存参数则是本机直接读出）。
AMD 官方 Zen 2 SOG（55723）与 WikiChip 在本机**网络不可达**（curl 返回 000），故未引用。

**(b) 数组工作负载（同一次运行，输入 24KB 全在 L1 内，ns/元素）**

| 运算 | 现状（库调用） | 内联（编译器可见） | 256-bit 打包 |
|---|---|---|---|
| vec2 加法 | 1.657 | 0.565 | **0.214** |
| vec2 点积 | 1.809 | — | 0.217（shuffle 版）/ 0.374（hadd 版） |
| vec2 缩放 | 1.446 | — | — |
| vec2 归一化 | 3.253 | — | — |
| vec3 加法 | 1.849 | 1.323 | — |
| vec3 点积 | 2.395 | — | — |
| vec3 叉积 | 12.954 | — | — |

（输入 384KB、越出 L1 后：库调用 1.657→1.737 ns、打包 0.214→0.230 ns —— 说明这组规模下还没到内存带宽墙，
瓶颈是调用边界与指令数，不是缓存。若把数组做到 L2 以上，打包版会先撞 store 带宽。）

**回答题目里的具体问题：vec2 只有 8 字节、一个 `__m256` 能装 4 个 —— 值传递 ABI 开销 vs 打包收益，谁先到瓶颈？**

**是调用边界先到瓶颈，不是 ABI，也不是寄存器宽度。** 三条证据：

1. **ABI 本身近乎免费**：`fv_vec2` 走 1 个 XMM（`addps; ret`），vec3 走 2 个，都不落栈；把 `fv_vec2` 换成两个
   `float` 参数也不会更快（唯一可见的成本是 `fv_vec2_dot` 里 GCC 用 `vmovq/shr/vmovd` 绕 GPR 取 lane 1，
   `-O2` 12 条 / `-O3 -march=native` 13 条指令 —— 这是可以只重写 dot 就修掉的小问题，与 ABI 形态无关）。
2. **边界比打包收益大 3 倍**：库调用 vs 内联 = 1.657−0.565 = **1.09 ns/元素（≈4.2 cycle）**；
   内联 vs 打包 = 0.565−0.214 = **0.35 ns/元素（≈1.3 cycle）**。即：先吃掉调用边界，才谈得上收打包的钱。
3. **算术根本不在瓶颈上**：vec2 加法只有 2 FLOP，本机 FP32 峰值 ≈ 2 条 FMA/cycle × 8 lane × 2 = 32 FLOP/cycle
   （≈130 GFLOP/s @4.08GHz）。现状 API 的吞吐 ≈ 2/1.657ns ≈ **1.2 GFLOP/s（峰值的 1%）**；打包版 ≈9 GFLOP/s，
   而且已经被 store 端口（每元素 8 字节）和循环开销挡住，不是被 FLOPs 挡住。

**因此**：`_n` 批量接口（`fv_vec2_add_n(fv_vec2 *restrict c, const fv_vec2 *restrict a, ...)`）应当先落地；
块大小按「≥12 条独立链 × 4 元素 = 48 个 vec2」取，且块内保持 AoS（不需要重排）。
值传递的 ABI 保留——它不是问题。

## 6 · 反例与陷阱

- **什么情况下 SIMD 会更慢**：
  ① 一次调用只算 1–2 个向量：打包权重（4 元素/指令）还没摊平调用与循环开销（实测库调用 1.66 ns vs 打包 0.21 ns/元素，
     但那是「批量」才成立的数字）；
  ② 归约密集：dot/len/normalize 每 4 个元素要多花 1 次 hadd/permute（实测 hadd 版 dot 比 shuffle 版慢 70%）；
  ③ 需要 AoS→SoA 重排：gather 每条指令 ≈10 cycle，重排成本远高于省下的算术；
  ④ 短依赖链：独立链 <10 条时 FMA 端口打不满（实测 8 条只有 1.4/cycle，12 条 2.0/cycle），此时向量化只增加
     寄存器压力和代码尺寸；
  ⑤ 数据不连续 / 跨步访问：打包反而放大未对齐惩罚（本机未量化，标**未验证**）。
- **不该做的**：
  ① 为了 SIMD 把 `fv_vec2` / `fv_vec3` 改成 4 分量 padding（内存 +33%~100%，破坏现有 ABI 与测试）；
  ② 全局 `-march=native` 或全局 `-mavx2`（二进制绑死机器；intrinsics 走 target 属性）；
  ③ 在没有批量接口之前手写 intrinsics（收益会被调用边界吃掉，见 §5）；
  ④ 把 flashvec 改成 SoA（会破坏「值语义、按值传、无分配」的 C API 契约，且物理引擎的调用点大多按粒子批处理）；
  ⑤ 全局 `-ffast-math` 做物理（NaN 传播、`x != x` 判定、归约重排全部改变语义）；
  ⑥ 用 `-mfma` 重写 `fv_vec2_dot` 之类的 2 元素运算（一条 `vfmadd` 只处理 1 个 vec2 的 2 个 lane，
     收益小于函数汇总开销）。
- **什么时候该做**：调用是批量的（n ≥ 数十）、运算是逐元素无归约的、数据在内存里连续且对齐可控 ——
  对 flashvec 就是 `add/sub/scale/lerp` 的 `_n` 版本，以及物理引擎里「一次处理一整段粒子数组」的积分/碰撞阶段。
- **可验证但本机没做的**：跨域（整数↔浮点）转发惩罚、store-forwarding 失败惩罚、AVX2 降频幅度。
  本机 `perf` 不可用（`perf_event_paranoid=4`，见 §0），要测它们需要 `sudo sysctl kernel.perf_event_paranoid=1`
  或 CAP_PERFMON；机制说明见 [kernel perf-security](https://www.kernel.org/doc/html/latest/admin-guide/perf-security.html)。
  在这些数字被实测之前，**不要**拿网上的「AVX2 会降频」传闻当设计依据。
