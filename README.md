# flashvec

极小的向量基础库：二维 / 三维向量，按值传递，不分配内存，无全局状态。
标量 API 与批量 API 并存：默认路径保证**逐位可复现**，批量路径在缓存驻留规模上快 2–3 倍（`dot` 这类归约快 20 倍）。

## 为什么是这个形状

「用 C++ 写 C」：头文件 `include/flashvec.h` 是纯 C 接口（带 `extern "C"` 保护），C 和 C++ 都能直接包含；
实现用 C++ 编译器编译，但代码只使用 C 语言子集。这样以后想引入 C++ 特性（重载、模板、RAII、`constexpr`）时，
既不用换编译器，也不用改对外 ABI。

## 目录

    include/flashvec.h        公开接口（纯 C，C11）
    src/flashvec.cpp          标量实现（C++ 编译，C 风格写法）
    src/flashvec_batch.cpp    批量入口 + AVX2 内核 + 运行时派发（P1）
    tests/test_flashvec.cpp   原有单元测试（30 条断言）
    tests/test_c_api.c        用 C 编译器编译，验证头文件是纯 C
    tests/test_p0.cpp         P0 新增函数（min/max/clamp/madd/aabb2）测试
    tests/test_batch.cpp      批量入口测试：位级一致性、尾长、掩码/索引语义
    tests/test_batch_c.c      批量入口的纯 C 调用验证
    bench/                    三档实现的性能基准（含内存带宽探针）
    verify/                   独立验收：位级对拍 + 性能对比（第三方编写）
    docs/                     优化报告与调研笔记

## 构建与测试

    make            # 生成 build/libflashvec.a
    make test       # 原有测试（30 条断言 + C 兼容性）
    make test-p0    # P0 新增函数测试
    make test-batch # 批量入口测试（位级一致性等）
    make test-all   # 上面三个一起跑
    make clean

默认编译选项：`-O3 -fno-math-errno -ffp-contract=off`。

- `-fno-math-errno`：放弃「数学函数失败时设置 errno」。本库 `sqrt` 的参数恒 >= 0，errno 从不被设置，
  但实测这能让含 `sqrtf` 的循环从「拒绝向量化」变成「向量化」，数值逐位不变。
- `-ffp-contract=off`：禁止把 `a + b*c` 收缩成 FMA。逐元素 API 承诺逐位可复现，而 FMA 收缩会让
  同一表达式在不同编译单元 / 不同目标下差 1 ulp（实测：默认的 `fast` 下会出现 `vfmadd`，
  且会让已有 API 的输出位模式发生变化）。源码里另有 `pragma GCC optimize ("fp-contract=off")` 做双保险。
- 不使用 `-march=native` / `-mavx2` / `-ffast-math`：默认产物保持 x86-64 基线，AVX2 只在
  `src/flashvec_batch.cpp` 里用函数级 `__attribute__((target("avx2,fma")))` + 运行时
  `__builtin_cpu_supports("avx2") != 0` 派发启用；非 AVX2 机器走标量回退，结果相同。

在自己的 C 或 C++ 项目里使用：把 `include/` 加到头文件路径，链接 `libflashvec.a`。

    #include "flashvec.h"

    fv_vec2 a = fv_vec2_make(3, 4);
    fv_vec2 n = fv_vec2_normalize(a);      /* (0.6, 0.8) */
    fv_real  d = fv_vec2_len(a);           /* 5 */

## API

**标量（逐元素）**：`fv_vec2` 的 `make / add / sub / neg / scale / lerp / dot / len_sqr / len / dist_sqr /
dist / normalize / project / near`；`fv_vec3` 同上，另有 `cross`。

**P0 新增**：`fv_vec2` / `fv_vec3` 的 `min / max / clamp / madd`；`fv_aabb2_merge / fv_aabb2_overlaps`。
NaN 语义显式定义为 `(a<b) ? a : b` 风格（**不是** IEEE `fminf` 的 NaN 吞噬语义），已写进头文件。

**P1 批量入口**（`include/flashvec.h` 里有完整契约注释）：

- AoS 4 路交错（一个 256-bit 寄存器正好装 4 个 `fv_vec2`，调用方数据布局**零改动**）：
  `fv_vec2_add_x4 / sub_x4 / scale_x4 / lerp_x4 / madd_x4`
- 稀疏变体（休眠体跳过 / broad-phase 稀疏更新用）：`fv_vec2_add_x4_mask`、`fv_vec2_add_x4_idx`、
  `fv_vec2_scale_x4_idx`。掩码语义是「掩码外写回 `a[i]` 原值」，不是「不写」。
- SoA 分量数组（归约型运算的自然布局）：`fv_vec2_add_soa / scale_soa / dot_soa / len_soa / normalize_soa`

### 批量入口的语义承诺（重要）

- **逐元素型批量函数与对应标量 API 逐位一致**（memcmp 浮点位模式）：`add/sub/scale/lerp/madd` 的
  `_x4` / `_soa` / `_mask` / `_idx`、`len_soa`、`normalize_soa`。为此这些路径**禁止 FMA 与倒数近似**，
  只用 `vaddps/vsubps/vmulps/vdivps/vsqrtps`，归约顺序也与标量一致。
- **唯一例外是 `dot_soa`**：它用 4 个累加器 + 末端一次归约，结果依赖归约顺序，**不保证逐位一致**；
  实测相对误差约 1e-6 量级（n < 8 时逐位一致）。需要逐位复现时请用标量 `fv_vec2_dot`。
- 长度 `n` 完全任意（含 0 与 `NULL`）：主循环处理 `n / 4`（或 `n / 8`），其余走标量尾巴。
- `_x4` / `_soa` 的指针须满足 `__restrict` 前提（输入输出不重叠）；`_idx` 变体目前是标量循环，
  **不**使用 `vgatherdps`（Zen 2 的 gather 以微码实现，实测约 10 cycle/条，不划算）。

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

1. **规模决定一切**。N = 1<<20（工作集 24 MB，超过 8 MB L3）时三档全部收敛到 ~1.0–1.5x，
   因为都撞在本机实测 6–7 GB/s 的流带宽墙上。**报加速比必须同时报工作集大小。**
2. **收益主要来自「跨 TU 调用边界消失」，不是 intrinsics**。同一份数学放进头文件内联后，
   `-O3 -march=native` 的自动向量化几乎追平手写 AVX2（实测 2.15x vs 2.10x）；
   而**使用预编译静态库的调用方无法内联**——这就是批量入口存在的理由。
3. 这台笔记本变频（1.64–3.97 GHz）且为共享机，**同一份代码跨运行可差 20–25%**，
   不同测试夹具的绝对值不可直接比较，只应看同一次运行内的比值。

完整的优化分析、证据与被否决的方案见 `docs/optimization-report.md`。

## 还没做（留给以后）

- 三维的批量入口（`fv_vec3` 12 字节，4 个装不满一个 256-bit 寄存器，需要单独设计布局；
  届时建议新增 16 字节的 `fv_vec4` 而不是给 `fv_vec3` 加 padding）。
- `_fast` 变体（例如 `rsqrtps` + 1 次牛顿迭代的 normalize，约 2 ulp 误差、约 10% 更快）：
  必须作为**新函数**提供，不能改默认语义。
- 近似三角函数（`include/flashvec_math.h`）：`atan` 多项式约比 glibc 快 2.7 倍，
  但必须写明输入域（参数归约才是 float 的真瓶颈）。
- 四维向量、矩阵、四元数、角度/旋转。
- 就地运算（`add_to` 之类）：目前测量显示收益来自批量入口，不来自减少拷贝。

详细的调研笔记（SIMD/硬件、算法、数据布局、实测）在 `docs/research/`。
