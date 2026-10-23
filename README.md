# flashvec

极小的向量基础库：二维 / 三维向量，按值传递，不分配内存，无全局状态。
标量 API 与批量 API 并存：默认路径保证**逐位可复现**，批量路径在缓存驻留规模上快 2–3 倍（`dot` 这类归约快 20 倍）。

## 为什么是这个形状

「用 C++ 写 C」：头文件 `include/flashvec.h` 是纯 C 接口（带 `extern "C"` 保护），C 和 C++ 都能直接包含；
实现用 C++ 编译器编译，但代码只使用 C 语言子集。这样以后想引入 C++ 特性（重载、模板、RAII、`constexpr`）时，
既不用换编译器，也不用改对外 ABI。

**批量入口默认是头文件里的 `static inline`**：跨 TU 调用边界被实测证明是一阶瓶颈（见
`docs/optimization-report.md` §3.1），而链接 `libflashvec.a` 的调用方无法把这层边界消掉，
所以实现只写一份、放在头文件里，`src/flashvec_batch.cpp` 只保留 `extern` 同名符号供
`FV_USE_EXTERN` 形态与已编译的旧消费方使用。头文件因此会在**支持的目标**上条件包含
`<immintrin.h>`（`FV_I_INTRINSICS && !FV_USE_EXTERN`）；用纯 C 编译器（`cc -std=c11`）、
以及 `FV_USE_EXTERN` 形态下，它仍然只是一个纯声明头文件。

## 目录

    include/flashvec.h        公开接口（纯 C，C11）；默认形态下批量入口的**实现也在这里**
    src/flashvec.cpp          标量实现（C++ 编译，C 风格写法）：fv_vec2 / fv_vec3
    src/flashvec_vec4.cpp     标量实现：fv_vec4（单独 TU，理由见头文件「四维」契约）
    src/flashvec_batch.cpp    extern 回退形态（FV_USE_EXTERN / FV_BUILD）：33 个批量入口的库内定义
    tests/test_flashvec.cpp   原有单元测试（30 条断言）
    tests/test_c_api.c        用 C 编译器编译，验证头文件是纯 C
    tests/test_p0.cpp         P0 新增函数（min/max/clamp/madd/aabb2）测试
    tests/test_batch.cpp      批量入口测试：位级一致性、尾长、掩码/索引语义
    tests/test_batch_c.c      批量入口的纯 C 调用验证
    tests/test_header_c.c     严格编译档（C11 -pedantic-errors）+ 同一 TU 双重包含
    tests/test_parity.c       默认（inline）形态 vs FV_USE_EXTERN 形态的逐字节对拍
    tests/test_stream.c       _stream（NT store）与 _mask_inplace 的行为回归
    tests/test_vec4.c         四维标量与 _x2/_soa 位级一致性
    tests/test_vec34_parity.c 三维 SoA / 四维 _x2 的形态对拍
    tests/consumer_probe.c    消费者默认选项下的内联效果普查
    bench/                    三档实现的性能基准（含内存带宽探针）
    bench/ab/                 同进程 A/B 夹具：A/A 地板、RFO 记账、五条验收判据
    verify/                   独立验收：位级对拍 + 性能对比 + 结构检查（第三方编写）
    docs/                     优化报告、落地记录、验收判据、调研笔记

## 构建与测试

    make            # 生成 build/libflashvec.a（3 个 TU）
    make test       # 原有测试（30 条断言 + C 兼容性）
    make test-p0    # P0 新增函数测试
    make test-batch # 批量入口测试（位级一致性等）
    make test-all   # 上面三个一起跑

    # 头文件 inline 形态（重构后）的验收目标
    make test-header         # 严格编译档（C11 -pedantic-errors / C++17 -Wpedantic）+ 双重包含
    make test-inline         # 默认（inline）形态 vs FV_USE_EXTERN 形态逐字节对拍
    make test-fma            # 消费者带 -march=native -mfma -ffp-contract=fast，仍必须位级一致
    make test-inline-effects # 内联效果普查（哪些入口必须完整内联、哪些只许每批一次调用）
    make test-stream         # _stream（NT store）与 _mask_inplace 的行为回归
    make test-vec4           # 四维标量 + _x2/_soa
    make test-vec34-inline   # 三维/四维批量入口的形态对拍
    make test-vec34-fma      # 同上，带敌意 FMA 选项
    make test-refactor       # 上面全部 + 导出符号清单 + 反汇编普查
    make clean

默认编译选项：`-O3 -fno-math-errno -ffp-contract=off`。

- `-fno-math-errno`：放弃「数学函数失败时设置 errno」。本库 `sqrt` 的参数恒 >= 0，errno 从不被设置，
  但实测这能让含 `sqrtf` 的循环从「拒绝向量化」变成「向量化」，数值逐位不变。
- `-ffp-contract=off`：禁止把 `a + b*c` 收缩成 FMA。逐元素 API 承诺逐位可复现，而 FMA 收缩会让
  同一表达式在不同编译单元 / 不同目标下差 1 ulp（实测：默认的 `fast` 下会出现 `vfmadd`，
  且会让已有 API 的输出位模式发生变化）。源码里另有 `pragma GCC optimize ("fp-contract=off")` 做双保险。
- 不使用 `-march=native` / `-mavx2` / `-ffast-math`：默认产物保持 x86-64 基线。AVX2 只在
  **头文件里**用函数级 `__attribute__((target("avx2")))` + 运行时
  `__builtin_cpu_supports("avx2") != 0` 派发启用（实测 32 处派发点全在 `include/flashvec.h`，
  `src/*.cpp` 里一个都没有）；非 AVX2 机器与所有尾长走标量回退，结果相同。
- 头文件里另有 `#pragma GCC optimize("fp-contract=off")` 挡在可能被收缩成 FMA 的入口外面：
  `target("avx2")` 管不住调用方命令行带来的 `-mfma`（实测 `-march=native -mfma -ffp-contract=fast`
  下带 target 属性的内核照样出现 `vfmadd`），消费者场景里 pragma 是唯一的防线。代价是这几个
  入口在消费者的默认选项下不内联（`make test-inline-effects` 把这个契约钉成可复现证据）。

在自己的 C 或 C++ 项目里使用：把 `include/` 加到头文件路径，链接 `libflashvec.a`。

    #include "flashvec.h"

    fv_vec2 a = fv_vec2_make(3, 4);
    fv_vec2 n = fv_vec2_normalize(a);      /* (0.6, 0.8) */
    fv_real  d = fv_vec2_len(a);           /* 5 */

## API

**标量（逐元素）**：`fv_vec2` 的 `make / add / sub / neg / scale / lerp / dot / len_sqr / len / dist_sqr /
dist / normalize / project / near`；`fv_vec3` 同上，另有 `cross`；`fv_vec4` 同样 18 个（`make / add /
sub / neg / scale / lerp / dot / len_sqr / len / dist_sqr / dist / normalize / project / near /
min / max / clamp / madd`，没有 `cross`）。

**P0 新增**：`fv_vec2` / `fv_vec3` / `fv_vec4` 的 `min / max / clamp / madd`；`fv_aabb2_merge / fv_aabb2_overlaps`。
NaN 语义显式定义为 `(a<b) ? a : b` 风格（**不是** IEEE `fminf` 的 NaN 吞噬语义），已写进头文件。

**批量入口**（共 33 个符号，`include/flashvec.h` 里有完整契约注释；默认形态是头文件里的 `static inline`，
`FV_USE_EXTERN` 形态转发到 `libflashvec.a` 里的同名符号）：

- AoS 4 路交错（一个 256-bit 寄存器正好装 4 个 8 B 的 `fv_vec2`，调用方数据布局**零改动**）：
  `fv_vec2_add_x4 / sub_x4 / scale_x4 / lerp_x4 / madd_x4`
- **AoS 2 路交错（`fv_vec4` 专用）**：`fv_vec4_add_x2 / sub_x2 / scale_x2 / lerp_x2 / madd_x2`。
  名字里的数字就是「一个 256-bit 寄存器装几个元素」：8 B 的 `fv_vec2` 装 4 个（`_x4`），
  16 B 的 `fv_vec4` 装 2 个（`_x2`）。两者都是**整装**，没有 lane 空转。
- 稀疏变体（休眠体跳过 / broad-phase 稀疏更新用）：`fv_vec2_add_x4_mask`（掩码外写回 `a[i]` 原值，
  不是「不写」）、`fv_vec2_add_x4_mask_inplace`（就地加，掩码外的**值**保持不变；按掩码密度在
  稠密/稀疏两条内核之间二选一）、`fv_vec2_add_x4_idx`、`fv_vec2_scale_x4_idx`。
- SoA 分量数组（归约型运算的自然布局）：`fv_vec2_add_soa / scale_soa / dot_soa / len_soa / normalize_soa`
- **三维只走 SoA**：`fv_vec3_add_soa / scale_soa / dot_soa / len_soa / normalize_soa`。
  `fv_vec3` 是 12 B，一个 ymm 装 2 个还空 8 B：AoS 打包要么空转 **1/4 lane**（浪费 25% 带宽），
  要么先 shuffle 拆包再打包（省下的宽度又赔回去）。所以三维不给 AoS 批量入口，也**不给 `fv_vec3`
  加 padding**（那会破坏 ABI、并让现有结构体直接变大 33%）；改为「同一分量的 8 个元素装一条 ymm」，
  lane 利用率 100%。
- **四维的 SoA 批量**：`fv_vec4_dot_soa / len_soa / normalize_soa`。
- **NT store（`_stream`）6 个**：`fv_vec2_add_x4_stream / sub_x4_stream / scale_x4_stream /
  add_soa_stream / scale_soa_stream / normalize_soa_stream`。结果用 `vmovntps` 写出，不经过 cache、
  不产生 RFO 的读 —— 由调用方显式选择，只在工作集明显大于 L3 时划算（物理原因见
  `docs/optimization-report.md` §2.1 的字节账）。

### 批量入口的语义承诺（重要）

- **逐元素型批量函数与对应标量 API 逐位一致**（memcmp 浮点位模式）：`add/sub/scale/lerp/madd` 的
  `_x4` / `_soa` / `_mask` / `_idx`、`len_soa`、`normalize_soa`。为此这些路径**禁止 FMA 与倒数近似**，
  只用 `vaddps/vsubps/vmulps/vdivps/vsqrtps`，归约顺序也与标量一致。
- **唯一例外是 `dot_soa`**：它用 4 个累加器 + 末端一次归约，结果依赖归约顺序，**不保证逐位一致**；
  实测相对误差约 1e-6 量级（n < 8 时逐位一致）。需要逐位复现时请用标量 `fv_vec2_dot`。
- 长度 `n` 完全任意（含 0 与 `NULL`）：主循环处理 `n / 4`（或 `n / 8`），其余走标量尾巴。
- `_x4` / `_soa` 的指针须满足 `__restrict` 前提（输入输出不重叠）；`_idx` 变体目前是标量循环，
  **不**使用 `vgatherdps`（Zen 2 的 gather 以微码实现，实测约 10 cycle/条，不划算）。
- **分段不变性（split invariance）**：逐元素型批量入口（`add/sub/scale/lerp/madd/len/normalize`
  及其 `_x4` / `_soa` / `_mask` / `_idx` 变体）具备该性质 —— 同一数组切成任意多段分别调用，
  输出与一次调用**逐位相同**。**归约型入口不具备**：`dot_soa`（含 `fv_vec3_dot_soa` /
  `fv_vec4_dot_soa`）的归约树随切分改变，实测 n=8192 时同一份输入「整段 vs 切 2 段」相对差
  可达 **2.4e-5**（n=1000 时约 3.9e-6），12/12 组语料全部不一致；同一次实验里 `add_x4` 作为
  对照组 12/12 逐位一致。这是归约的**设计必然**，不是缺陷：需要「结果与切分无关」的调用方
  必须用标量 `fv_vec2_dot` 顺序累加，或自己固定切分方式（分块并行时尤其注意）。
- 上一条容易与 `dot_soa` 的**误差**混为一谈，两者必须分开读：头文件里「n ≤ 1000 相对误差 ≤ 4e-8」
  讲的是**与标量路径的差**；「整段 vs 切段 2.4e-5」讲的是**与另一种切分的差**。前者不能当后者的
  上界用。

### 性能（本机实测，AMD Ryzen 7 5700U，8C/16T，AVX2/FMA，无 AVX-512）

N = 8192（每数组 32 KB，L2 驻留）相对**库调用式标量循环**的加速比：

| 运算 | `_x4` | `_soa` |
|---|---|---|
| add | ~2.1x | ~1.9x |
| scale | ~2.6x | ~2.7x |
| lerp | ~2.5x | — |
| dot | — | **~20x** |
| len | — | ~2.7x |
| normalize | — | ~3.4x |
| collision 微内核（法线归一化 + 速度反射） | — | ~5.1x |
| add_mask（5% / 50% / 95% 命中） | ~1.7x / ~3.8x / ~3.3x | — |
| add_idx（5% / 50% 命中） | ~1.1x / ~1.4x | — |

必须一起读的限定：

1. **规模决定一切**。本机 L3 是 **4 MiB/CCX、两个 CCX 合计 8 MiB**（`/sys/devices/system/cpu/
   cpu0/cache/index3` 报 `level=3 size=4096K shared_cpu_list=0-7`）——**单线程只吃得到自己那 4 MiB**；
   8 MiB 是芯片总量，不等于一个线程能用的量（旧文档按 8 MB 算，把拐点算高了一倍）。按 add 的
   24 B/元素算，拐点是 **4 MiB ÷ 24 B ≈ 1.75×10⁵ 个活跃元素**。实测正好卡在这条线上：同进程
   A/B 夹具在 2^17（工作集 3 MB）测得 38–57 GB/s，到 2^18（6 MB）掉到 5.6–8.8 GB/s。
   N = 1<<20（工作集 24 MB）时三档全部收敛到 ~1.0–1.5x，因为都撞在本机实测 6–7 GB/s 的
   流带宽墙上。**报加速比必须同时报工作集大小。**
2. **收益主要来自「跨 TU 调用边界消失」，不是 intrinsics**。同一份数学放进头文件内联后，
   `-O3 -march=native` 的自动向量化几乎追平手写 AVX2（实测 2.15x vs 2.10x）；
   而**使用预编译静态库的调用方无法内联**——这就是批量入口存在的理由。
3. 这台笔记本变频（1.64–3.97 GHz）且为共享机，**同一份代码跨运行可差 20–25%**，
   不同测试夹具的绝对值不可直接比较，只应看同一次运行内的比值。

`_stream`（NT store）与 `_mask_inplace` 的加硬实测（`bench/ab` **同进程交错 A/B + 同轮 A/A 地板**，
AB/BA × 两套独立构建；比值只可同次运行内比，括号里是该档的读数离散度）：

| 场景 | 实测 | 为什么是这个方向 |
|---|---|---|
| `_stream` @ 24 MiB 工作集（n=1M） | add **1.350** / sub **1.293** / scale **1.313**（A/A 地板 1.08–1.12，配对 IQR 30–35%） | 远超 L3，写流是一次性的：普通 store 要为写 miss 多读 8 B/元素（add 真实 32 B vs 名义 24 B），NT store 省掉它 |
| `_stream` @ 192 MiB 工作集（n=8M） | add **1.422** / sub **1.519** / add_soa **1.582**（地板 1.03–1.08，IQR 28–35%） | 同上（DRAM 档）；`add_soa` 自己的上限是 **32/24 = 1.333**（它读 4 条数组，不是 2 条） |
| `_stream` @ 96 MiB（n=4M） | 6 个入口**全部 INCONCLUSIVE**（A/A 地板 1.34–1.61，IQR 66–98%） | 这一档测不准：同一份代码的两侧比值就能到 1.5 量级，**不得当判据** |
| `_stream` @ 3 MiB（n=131072，压在 L3 边界上） | 6 个入口**全部 INCONCLUSIVE（黑名单档）** | L3 边界 + 4K 别名：该档 A/A 比实测 0.50–2.22，可测、可记录，但不得用来判定 |
| `_stream` @ 工作集在 L2 内（n=8192） | add **0.217** / sub **0.301** / scale **0.130** / add_soa **0.229** / scale_soa **0.138** / normalize_soa **0.506** → **0.13–0.51x**（慢 2–8 倍） | 反例：数据本来就在 cache 里，NT store 只多付清行的代价（6 个入口全部判 FAIL(变慢)） |
| `_mask_inplace` @ 1%/5%/20%/50%/95% × 4 档规模 | 20 个点**全部 INCONCLUSIVE**，增益 **-0.10…+0.00**（同轮 A/A 地板 1.06–1.56） | 就地写的目标数组本来就要先读，**不产生 RFO**，省不掉流量；价值在形态、不在速度（旧的「量级 1.5–2x」已作废） |

> 三条读法：① 上限是**字节账的期望值**，不是物理保证 —— add 系 **32/24 = 1.333**、scale 与
> normalize **24/16 = 1.500**，一个元素搬几个 float 只能逐入口算；② 4M 档配对 IQR 66–98%、1M/8M
> 档 20–35%，所以**单点超过上限 20% 以内一律视为噪声**（`add_soa` 1.582 对 1.333 超出 19%，
> 正在噪声带内），不得声称「实测超过了理论上限」；③ 每个数都带同轮 A/A 地板，判定表在
> `docs/research/raw/ab_campaign_verdict.txt`（每次 `make run` 覆盖写，引用请连同时刻），判据见
> `docs/acceptance-protocol.md`（含 §9「走不通的路」、§10「已知局限」）。

**还没测到的，不要引用成数字**：三维 SoA 批量、四维 `_x2`/`_soa` 批量相对**标量**的加速比 —— 待测。
本轮 campaign 只测了它们的**形态差**（库内定义 vs 头文件内联）：13 个入口 × 5 档规模 = 65 个点，
5 PASS / 8 FAIL(变慢) / 52 INCONCLUSIVE，同一份源码的读数在 **0.97–1.93** 之间散开 —— 这组数字
说明的是形态，不是 SIMD 收益。这台机器长期有外部负载，**单点比值不可信**，判据必须建立在逐规模
中位比、全点中位比与同轮 A/A 地板上，见 `docs/acceptance-protocol.md`。

完整的优化分析、证据与被否决的方案见 `docs/optimization-report.md`；
一条优化要「被接受」需要同时满足哪五条，见 `docs/acceptance-protocol.md`。

## 还没做（留给以后）

- `_fast` 变体（例如 `rsqrtps` + 1 次牛顿迭代的 normalize，约 2 ulp 误差、约 10% 更快）：
  必须作为**新函数**提供，不能改默认语义。
- 近似三角函数（`include/flashvec_math.h`）：`atan` 多项式约比 glibc 快 2.7 倍，
  但必须写明输入域（参数归约才是 float 的真瓶颈）。
- 矩阵（`fv_mat3` / `fv_mat4`）、四元数、角度/旋转。
- 就地运算的其余形态（`add_to` 之类）：目前测量显示收益来自批量入口，不来自减少拷贝。
  唯一的例外 `_mask_inplace` 已经做了 —— 它存在是因为稀疏更新需要「只改命中元素的值」这个
  **语义**，不是为了少写一个参数（`_mask` 的契约是掩码外写回 `a[i]`，`_mask_inplace` 是掩码外
  值不变，两者不是同一个契约）。

已经做完、旧版这里说要留到以后的：三维批量入口（走 SoA）与四维向量（`fv_vec4` + `_x2` AoS 批量
+ `_soa` 归约批量）——见上面的 API 节。

详细的调研笔记（SIMD/硬件、算法、数据布局、实测）在 `docs/research/`。
