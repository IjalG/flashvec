# flashvec

> [English](../README.md) | 简体中文

一个很小的向量基础库，面向 2D / 3D / 4D 实时计算。

- 全部按值传递，库内不分配内存，没有全局状态。
- 头文件是纯 C 接口，C 和 C++ 都能直接包含。
- 提供标量和批量两套 API。标量用于正确性基线，批量用于热路径。
- 默认路径逐位可复现：同样的输入永远给出同样的位模式。

## 快速开始

```c
#include "flashvec.h"

fv_vec2 a = fv_vec2_make(3, 4);
fv_vec2 n = fv_vec2_normalize(a);   /* (0.6, 0.8) */
fv_real d = fv_vec2_len(a);         /* 5 */

/* 批量：一次处理整个数组 */
fv_vec2 pos[N], vel[N];
fv_vec2_add_x4(pos, vel, pos, N);
```

在自己的项目里使用：把 `include/` 加进头文件搜索路径。默认形态完全由头文件提供实现，不需要链接库。想用预编译库的形态，就定义 `FV_USE_EXTERN` 并链接 `libflashvec.a`。

## 类型

| 类型 | 大小 | 说明 |
|---|---|---|
| `fv_vec2` | 8 B | |
| `fv_vec3` | 12 B | 不加 padding，保持用户既有布局 |
| `fv_vec4` | 16 B | `static_assert` 锁死大小 |
| `fv_aabb2` | 16 B | 轴对齐包围盒 |
| `fv_quat` | 16 B | 旋转；单位元 (0,0,0,1)；与 `fv_vec4` 布局兼容，但是不同类型 |
| `fv_mat3` | 36 B | 行主序 3x3，与 `fv_mat4` 同约定（与 GLM 相反），9 float 紧凑无 padding |

## 标量 API

每个向量类型各有 18 个核心函数（`fv_vec3` 额外多一个 `cross`）：

`make` `add` `sub` `neg` `scale` `lerp` `min` `max` `clamp` `madd`
`dot` `len_sqr` `len` `dist_sqr` `dist` `normalize` `project` `near`

在此之上，三个类型都有 `addmul`（`x + v*s`），只有 `fv_vec2` 有 `maddv`（第三操作数是向量的 `a*b + c`）。

包围盒另有 `fv_aabb2_merge` 和 `fv_aabb2_overlaps`。

默认形态把每个标量函数体写成 header 里的 `static inline`，与批量内核共用同一个 `fp-contract=off` pragma；`FV_USE_EXTERN` 切到库内符号，两种形态逐位一致（`make test-scalar-parity`）。混合标量循环实测（`make bench-scalar`）：默认形态比 `FV_USE_EXTERN` 快 **1.70–1.82x**；`-mfma -ffp-contract=fast` 的消费者仍然 **0 条 `vfmadd`**。默认形态是 `static inline`；需要稳定的库符号 / 函数地址 / ABI 边界时用 `FV_USE_EXTERN`。`normalize` 即使内联也保留标量 `sqrt`/`div` 与零分支，追平它要用 `_fast` 族。

`min` / `max` / `clamp` 的 NaN 语义是 `(a < b) ? a : b`，不是 IEEE 的 `fminf`。这一点写进了头文件。

## 批量 API

52 个向量入口，按内存布局分组；另有 4 个矩阵变换入口。名字里的数字表示一个 256-bit 寄存器装几个元素：`fv_vec2` 装 4 个，叫 `_x4`；`fv_vec4` 装 2 个，叫 `_x2`。两者都能装满，没有 lane 浪费。

其中 19 个是列表里已有运算的 `_stream` 变体，只有 store 指令不同；见下面「流式写入」一节。

矩阵变换 `fv_mat4_transform_points` 与 `fv_mat4_transform_dirs` 各带一个 `_stream` 形态，见下面的「矩阵变换」一节。

### AoS 交错

调用方数组原样传入，不需要改布局。

| 类型 | 函数 |
|---|---|
| `fv_vec2` | `add_x4` `sub_x4` `scale_x4` `lerp_x4` `madd_x4` |
| `fv_vec4` | `add_x2` `sub_x2` `scale_x2` `lerp_x2` `madd_x2` |

### SoA 分量数组

把 x、y、z 分成独立数组。归约型运算只有这种布局才高效。

| 类型 | 函数 |
|---|---|
| `fv_vec2` | `add_soa` `scale_soa` `dot_soa` `len_soa` `normalize_soa` |
| `fv_vec3` | `add_soa` `scale_soa` `dot_soa` `len_soa` `normalize_soa` |
| `fv_vec4` | `dot_soa` `len_soa` `normalize_soa` |

三维只有 SoA 形态。`fv_vec3` 是 12 B，一个寄存器装两个会空出 8 B；要填满就得 shuffle 拆包再打包，省下的宽度又赔回去。
给结构体加 padding 会破坏 ABI，也让它变大 33%。
所以三维走「同一分量的 8 个元素装一条寄存器」，lane 利用率 100%。

### 融合、就地、三操作数

物理步进反复需要同样这三个表达式，而每个都写成两次库调用，就多付一趟访存。

```c
/* out = normalize(a * s + b)，一趟；与 scale → add → normalize 逐位相同 */
void fv_vec2_madd_normalize_soa(const float* ax, const float* ay,
                                const float* bx, const float* by,
                                fv_real s, float* ox, float* oy, size_t n);

/* x += v * s，就地 */
void fv_vec2_addmul_x4 (fv_vec2* x, const fv_vec2* v, fv_real s, size_t n);
void fv_vec2_addmul_soa(float* x, float* y, const float* vx, const float* vy,
                        fv_real s, size_t n);
void fv_vec3_addmul_soa(float* x, float* y, float* z,
                        const float* vx, const float* vy, const float* vz,
                        fv_real s, size_t n);
void fv_vec4_addmul_x2 (fv_vec4* x, const fv_vec4* v, fv_real s, size_t n);

/* o = a * b + c；三向量操作数，而 madd 的第三参数是标量 */
void fv_vec2_maddv_x4(const fv_vec2* a, const fv_vec2* b, const fv_vec2* c, fv_vec2* o, size_t n);
```

对应的逐值形态是 `fv_vec2_addmul`、`fv_vec3_addmul`、`fv_vec4_addmul` 和 `fv_vec2_maddv`；
上面每个批量入口都与自己的逐值形态逐位相同。

就地写不产生 RFO，因为这次 store 落在刚刚读进来的那条 cache line 上。所以这一族是唯一上限超过 2 倍的：
两趟 56 B/元素，就地 24 B/元素。

`x` 既读又写，`v` 不得与它重叠。`n == 0` 时不碰任何指针。

融合管线与 `maddv` 都写独立输出，所以两者都有 `_stream` 版本。就地 `addmul` 一族故意没有：
那次 store 已经落在核自己持有的行上。

### 稀疏更新

| 函数 | 语义 |
|---|---|
| `add_x4_mask` | 掩码非零的元素做加法，掩码为零的元素写回 `a[i]` 原值 |
| `add_x4_mask_inplace` | 就地做加法。掩码为零的元素值不变；不保证没有被物理写回 |
| `add_x4_idx` | 只处理 `idx` 指定的下标，可重复、可乱序 |
| `scale_x4_idx` | 同上，做缩放 |

`_idx` 是标量循环，没有用 `vgatherdps`。Zen 2 上 gather 由微码实现，约 10 周期一条，不划算。

### 流式写入

一共 21 个 `_stream` 变体：19 个向量入口，加两个矩阵变换。向量那边覆盖 `_x4` 与 `_x2` 形态的
`add` / `sub` / `scale` / `lerp` / `madd`，`fv_vec2` 的 SoA `add` / `scale` / `normalize`（含融合管线），
`fv_vec3` 的 `add` / `scale` / `normalize`，`fv_vec4` 的 `normalize_soa`，以及 `maddv_x4`。

用非临时存储写出结果，绕过 cache，也省掉 RFO 那一次读。**只在工作集不小于 L3（本机 4 MiB per CCX）、
且结果不会马上被读回时使用。** 低于这条线时它不是略慢，而是**慢 2–8 倍**：L2 内 n=8192 实测，
与普通 store 的比值是 **0.13–0.37**。同一轮在 L3 之外读到 **1.48–1.60x**。

这个选择必须由调用方做，库不会自动切换。

## 物理原语（P6 / P7）

向量层之上的领域类型：旋转 `fv_quat`、3x3 线性代数 `fv_mat3`，以及既有 `fv_mat4` 上的通用 4x4 运算。四元数不是 `fv_vec4`：乘法是不可交换的复合，零不是旋转，批量布局是四个分量数组。

| 类型 | 大小 | 说明 |
|---|---|---|
| `fv_quat` | 16 B | `x/y/z/w`，单位元 `(0,0,0,1)`；与 `fv_vec4` 布局兼容但是不同类型 |
| `fv_mat3` | 36 B | 行主序 `m[row*3+col]`，与 `fv_mat4` 同约定（与 GLM 相反），9 float 紧凑无 padding |

标量 API：`fv_quat_make / add / sub / neg / scale / dot / len_sqr / len / normalize / conjugate / inverse / mul / rotate_vec3 / nlerp / slerp / from_axis_angle / from_euler_xyz`；`fv_mat3_from_rows / transform_vec3 / mul / transpose / from_quat / from_mat3 / determinant / inverse / inverse_checked`；`fv_mat4_mul / fv_mat4_inverse / fv_mat4_inverse_checked`；`fv_atan_approx`。

`fv_quat_mul(a,b)` 是 `a*b`，先作用 `b` 再作用 `a`；`fv_quat_normalize` 对长度 <= `FV_EPS` 返回单位元（不是零），这是与向量 normalize 唯一的有意差异；`rotate_vec3` 与 `from_quat` 假定 `q` 是单位四元数，不自动归一化。

`fv_quat_inverse` 是一般逆 `conjugate(q) / |q|²`（单位四元数下等于共轭）；长度 <= `FV_EPS` 返回**零**四元数而不是单位元——逆不存在，零是信号。`fv_mat3_inverse` 是伴随矩阵 / 行列式；当 `|det| <= FV_EPS * max|entry|³`（相对阈值，所以 `1e-3 * I` 仍可逆）时视为奇异，按值版本返回零矩阵。`_checked` 版本返回 1 / 0，失败时写入零值；`fv_mat3_inverse_checked` 允许 `out == m`。

`fv_quat_nlerp` / `fv_quat_slerp` 走最短弧（dot 为负时先翻转第二个四元数），结果归一化，`t` 不截断；slerp 在 dot 接近 1 时回退到 nlerp。`fv_quat_from_axis_angle` 内部归一化轴（零轴返回单位元），`fv_quat_from_euler_xyz` 的顺序是先 X、再 Y、再 Z（`q = qz*qy*qx`）。这些入口用到 `sin/cos`，是精度有界、不可跨 libm 逐位复现的。

`fv_quat_from_mat3` 用 Shepperd 的“最大对角分支”法，假定输入是旋转矩阵；`fv_mat4_mul` / `fv_mat4_inverse` 是**通用** 4x4 运算（不只仿射），逆用伴随矩阵 / 行列式闭式解，奇异阈值同样是相对的（`|det| <= FV_EPS * max|entry|⁴`），按值版本返回零矩阵，`fv_mat4_inverse_checked` 返回 1 / 0 且允许 `out == m`。

批量（SoA，每 ymm 8 个元素）：`fv_quat_mul_soa`、`fv_quat_normalize_soa`、`fv_quat_rotate_vec3_soa`、`fv_mat3_transform_vec3_soa`、`fv_mat3_mul_soa`、`fv_mat3_inverse_soa`。前六个对任意 n 与尾长都与逐值 API 逐位一致；mat3 乘 / 逆现在有独立 AVX2 内核：乘把 9 个 B 列常驻寄存器、每次算 3 个输出行；逆每 8 个元素只做一次 `vdivps(1/det)`（而不是 9 次标量除法），奇异 lane 用 blend 写零。尾长仍走标量，逐位不变。本地 `make bench-mat3`（P0 后的默认内联标量基线）：乘 8192 4.9x、1M 3.1x；逆 9.7x、2.7x。相对旧的 `FV_USE_EXTERN` 标量形态，同一批内核测得 13.3x / 3.5x 与 14.7x / 3.2x。这些是本地基线，不是 GLM 报告的比值。`fv_vec2/vec3/vec4_normalize_soa_fast` 与 `fv_quat_normalize_soa_fast` 是独立的近似族（`rsqrtps` + 一次牛顿迭代，约 1–2 ulp，标量回退用精确 `1/sqrtf`），不承诺与精确入口逐位一致。

实现：quat 乘 / 旋转、mat3 变换是 `target("avx2")` 里的纯 C，GCC 自动向量化成 256 位且零 shuffle（实测 40 / 40 / 21 条 AVX2 指令，0 个 `vpermps`/`vshufps`）；quat normalize 是编译器唯一没向量化的（保留 `sqrtss`/`divss` 与分支），因此单独写了 AVX2 内核（29 条 ymm）。`_fast` 族用 `rsqrtps` + 一次牛顿迭代；`fv_vec2_madd_normalize_soa_fast` 是融合近似路径，8192 档比精确融合快 1.07x，1M 档在 A/A 噪声内。精确融合入口另有**可选 AVX-512 backend**：`target("avx512f")` 的 16-lane 内核，运行时由 `__builtin_cpu_supports("avx512f")` 选择；AVX2 路径不变，`FV_I_NO_AVX512` 可移除。本机没有 AVX-512，只验证了 codegen（存在 `vrsqrt14ps` / `zmm`），运行时路径仍需 AVX-512 硬件验证。`fv_atan_approx` 是 `[-1,1]` 上的 5 项多项式（最大绝对误差约 1.2e-5，域外返回 NaN）。

初步测量（`make bench-physics`，同进程交错 A/B，11 轮配对，`taskset -c 0`，两次运行中位数；A/A 地板 8192 约 0.99–1.00、1M 约 0.93–1.07）：

| 入口 | n = 8192 | n = 1M |
|---|---|---|
| `fv_quat_mul_soa` | 18.0–19.5x | 2.39–2.51x |
| `fv_quat_normalize_soa` | 41.9–42.3x | 4.56–5.13x |
| `fv_quat_rotate_vec3_soa` | 25.8x | 3.39–3.57x |
| `fv_mat3_transform_vec3_soa` | 24.7–25.1x | 3.33–3.80x |

这些表测于 P0 之前：基线是 **`FV_USE_EXTERN`** 标量形态（跨 TU 调用）。P0 之后默认标量形态是 header 内联，批量相对默认标量的比值会更小；表仍对旧/extern 基线有效。这些是初步数据，尚未走五条验收标准。

## 语义承诺

### 逐位一致

下面这些入口保证与对应标量 API 的浮点位模式完全相同，可以直接 `memcmp`：

- `add` / `sub` / `scale` / `lerp` / `madd` 的 `_x4`、`_x2`、`_soa`、`_mask`、`_idx` 变体
- `len_soa`
- `normalize_soa`
- 融合管线 `madd_normalize_soa` 与 `scale_soa` → `add_soa` → `normalize_soa` 三次调用
- `addmul` / `maddv` 与对应逐值形态
- `fv_quat_mul_soa` / `normalize_soa` / `rotate_vec3_soa`，`fv_mat3_transform_vec3_soa` / `mul_soa` / `inverse_soa` 与对应逐值形态
- `_stream` 系列与对应非 stream 版本

实现上只用 `vaddps`、`vsubps`、`vmulps`、`vdivps`、`vsqrtps`，不用 FMA，不用倒数近似，归约顺序与标量一致。`_fast` 入口与 `fv_atan_approx` 是**独立的近似族**，不在上面的逐位承诺内，不能混用。

### 唯一的例外：`dot_soa`

`dot_soa` 用 4 个累加器再归约一次，结果依赖归约顺序，不保证逐位一致。n < 8 时它退化成顺序累加，这时与标量路径一致。需要逐位复现请用标量 `dot`。

头文件里两组误差数字要分开读：

| 说法 | 含义 |
|---|---|
| n ≤ 1000 相对误差 ≤ 4e-8 | 与**标量路径**的差 |
| 整段 vs 切段 2.4e-5 | 与**另一种切分**的差 |

前者不能当作后者的上界。

### 分段不变性

逐元素型入口（所有 `add` / `sub` / `scale` / `lerp` / `madd` / `len` / `normalize` 的批量变体）具备这个性质：同一数组切成任意多段分别调用，输出与一次调用逐位相同。

归约型入口不具备。`dot_soa` 的归约树随切分改变，实测 n=8192 时「整段」与「切两段再相加」相对差 2.4e-5，12 组语料全部不一致。同一实验里 `add_x4` 作为对照 12 组全部逐位一致。

这是归约的设计必然，不是缺陷。分块并行时如果需要结果与切分无关，请用标量 `dot` 顺序累加，或自己固定切分方式。

### 尾长与指针

长度 `n` 任意，包括 0 和 `NULL`。主循环处理整块，剩下的走标量尾巴。

`_x4` 与 `_soa` 的指针需满足 `__restrict` 前提，输入输出不得重叠。

### 矩阵变换

```c
typedef struct { float m[16]; } fv_mat4;   /* 行主序：m[row*4 + col] */

fv_mat4 fv_mat4_identity(void);

void fv_mat4_transform_points(const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n);
void fv_mat4_transform_dirs  (const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n);
```

`transform_points` 把 `in[i].w` 当 1，平移生效；`transform_dirs` 当 0，平移被忽略，`out[i].w` 为 0。

**这是行主序，与 GLM、OpenGL 相反。** GLM 的矩阵传入前必须先转置。

这是全库唯一一个手写 SIMD 真的赢过编译器的地方。编译器确实会向量化朴素的 AoS 循环，但它生成 440 条指令、14 次 `vpermps`；手写内核是 107 条指令、3 次。在 L1/L2 驻留的规模上，它比**用本库默认编译档编的标量循环**快 **1.81 倍**（n=8192 时 7 次全部落在 1.81–1.84），若把标量基线改用 `-march=native` 编，则是 2.25 倍。引用这个数时请连编译档一起给：基线的编译档能把倍数挪动 30%。超过约 131072 个元素后它收敛到 1.0x，因为每元素 48 字节，和别的运算一样撞 DRAM 墙。

没有 SoA 版本。实测它是 1.01x——编译器已经把 SoA 做到最优（探针在 `bench/tfprobe/probe.cpp`，任何人都能重跑）。

#### 逐位一致唯一够不到的地方

如果某条算术指令的**两个操作数同时是 NaN 且位模式不同**，留下哪一个由这条指令的操作数顺序决定。这个顺序由编译器挑，而它**在不同构建里挑得不一样**：实测同一个入口，调用方用 `-march=native` 编译时与链接库内实现时，会给出不同的 NaN payload。于是标量尾巴、向量内核、以及同一内核换个编法，都可能对同一个元素给出不同（但都合法）的 NaN 位模式——**结果可能随 `n` 和你的编译选项一起变**。

除这一格之外，对所有 `n` 和所有尾长都逐位一致。语料里若最多只有一种 NaN 位模式，这一格也成立。

## 性能

测试机：AMD Ryzen 7 5700U，8 核 16 线程，AVX2 + FMA，无 AVX-512。

### 缓存边界

L3 是 **4 MiB per CCX**，两个 CCX 合计 8 MiB。单线程只能用到自己那 4 MiB。
这一点可以从 `/sys/devices/system/cpu/cpu0/cache/index3` 读到：`size=4096K`, `shared_cpu_list=0-7`。

按 `add` 的 24 B/元素算，拐点在 4 MiB ÷ 24 B ≈ **1.75×10⁵ 个活跃元素**。实测吻合：2^17（工作集 3 MB）是 38–57 GB/s，2^18（6 MB）掉到 5.6–8.8 GB/s。

### 机器的上限

下面的数字都是这台机器能搬的量的一个分数，所以先给上限。来自 `bench/membw.c`，单线程、taskset 绑核、5 次取最小：

| 工作集 | triad | triad16（含 RFO） | read | memcpy |
|---|---|---|---|---|
| 96 MB | 7.88 GB/s | 10.51 GB/s | 13.28 GB/s | 12.93 GB/s |
| 48 MB | 9.16 GB/s | 12.21 GB/s | 14.31 GB/s | 14.31 GB/s |

`triad` 按逻辑字节算，12 B/元素。`triad16` 是同一次测量把写那一次的 RFO 也算进去，16 B/元素。

逐元素基线：`fv_vec2_add_x4` 在 96 MB 工作集上是 **3.12–3.33 ns/元素**，约 **9.6–10.3 GB/s**；同工作集下机器 triad16 上限是 10.51 GB/s，即 **91–98%**。
这个规模上没有 10 倍可拿，剩下的就是流量本身。

### 相对标量调用的加速比

下表用历史基线：**`FV_USE_EXTERN` 形态的逐元素库调用**（跨 TU、不可内联）。P0 之后默认标量形态已是 header 内联，同 TU 的循环也会自动向量化；表保留为旧基线记录，不再代表默认形态的比值。`make bench-scalar` 直接测两种标量形态（混合负载 1.70–1.82x）。

N = 8192，每数组 32 KB，L2 驻留：

| 运算 | `_x4` | `_soa` |
|---|---|---|
| add | ~2.1x | ~1.9x |
| scale | ~2.6x | ~2.7x |
| lerp | ~2.5x | — |
| dot | — | ~20x |
| len | — | ~2.7x |
| normalize | — | ~3.4x |
| 碰撞微内核 | — | ~5.1x |
| add_mask（5% / 50% / 95%） | ~1.7x / ~3.8x / ~3.3x | — |
| add_idx（5% / 50%） | ~1.1x / ~1.4x | — |

### `_stream` 的实测

数据来自 `bench/ab` 的同进程交错 A/B，带同轮 A/A 地板，AB/BA 换序 × 两套独立构建。比值只在同一次运行内可比。

| 场景 | 结果 |
|---|---|
| 24 MiB 工作集（n=1M） | add 1.350 / sub 1.293 / scale 1.313 |
| 192 MiB 工作集（n=8M） | add 1.422 / sub 1.519 / add_soa 1.582 |
| 96 MiB（n=4M） | 全部无法判定，该档 A/A 地板 1.34–1.61 |
| 3 MiB（L3 边界） | 全部无法判定，该档 A/A 比值 0.50–2.22 |
| L2 内（n=8192） | add 0.217 / sub 0.301 / scale 0.130 / add_soa 0.229 / scale_soa 0.138 / normalize_soa 0.506 |

最后一行是反例：数据本来就在 cache 里，非临时存储只会多付清行的代价。六个入口在这一档全部被判为变慢，区间 0.13–0.51x。

关于上限：字节账给出的期望值，`add` 系是 32/24 = 1.333，`scale` 与 `normalize` 是 24/16 = 1.500。一个元素搬几个 float 只能逐入口算。
实测值分布较散，单点超过上限 20% 以内都算噪声。

### 融合与就地入口的实测

同一副夹具、同样的条件：单线程、taskset 绑核、默认档 `-O3 -fno-math-errno -ffp-contract=off`（**不带 `-march`**）、
同进程交替 A/B、每行都配同轮 A/A 地板。

| 项目 | n=8192 | n=1M | n=8M | A/A 地板 |
|---|---|---|---|---|
| `_stream` vs 普通 store | — | 1.48–1.60x | 1.54–1.60x | 0.97–1.07 |
| 融合管线 `normalize(a*s+b)` vs 三次调用 | 1.73x | 2.43–2.55x | 2.53–2.63x | — |
| 同上，再叠 `_stream` | 1.68x | **3.53–3.68x** | **3.72–3.78x** | — |
| 就地 `x += v*s` vs `scale` → `add` 两趟 | 1.06–2.18x | 2.05–2.17x | 1.92–2.03x | 0.97–1.06 |
| `maddv` vs `mul` → `add` | 1.60–1.67x | 1.61–1.79x | 1.59–1.76x | ~1.00 |
| `maddv_stream` vs `mul` → `add` | 1.16–1.37x | 2.16–2.23x | **2.34–2.41x** | ~1.00 |

`_stream` 那一行来自另一轮。上一节那次 campaign 在同样两档读到 1.29–1.58x、L2 内读到 0.13–0.51x，这里是 1.48–1.60x 与 0.13–0.37；两者在本机跨轮散布之内是同一个结论。

有三件事这张表没说：

- **`_stream` 在小规模上更慢，而且不是略慢：慢 2–8 倍**（n=8192 时 0.13–0.37）。NT store 绕过 cache，
  在 L2 内既没有 RFO 可省，又失去了缓存命中。只在工作集不小于 L3（4 MiB per CCX）时用它。
- **`maddv_stream` 在 n=8192 是「相对更快、绝对更慢」**：相对两趟 workaround 是 1.16–1.37x，
  但绝对耗时 0.99 ns/元素，非 stream 融合版是 0.58。L2 内应该用非 stream 版。
- **字节账是上限，超过上限 20% 以内属于噪声。** `add` / `sub` / `lerp` / `madd` 的上限是 32/24 = 1.33x，
  `scale` / `normalize` / 矩阵变换是 1.50x，就地 `addmul` 是 2.33x（56/24），`maddv` 非 stream 1.60x、stream 2.00x。
  不要凭单点声称「超过理论上限」。`_stream` 那一行读到 1.48–1.60x，略高于 `add` 系的 1.33x 上限；
  这个超出是真实的，也有原因——NT store 还避免了输出污染 cache，这一层字节账没有计入。

### `_mask_inplace` 的实测

5 种密度 × 4 档规模共 20 个点，全部无法判定，增益在 −0.10 到 +0.00 之间。

原因在流量账上。就地写之前目标数组已经被读过，对应的 cache line 本来就在 cache 里，普通 store 这时不产生 RFO。所以「就地」省掉的是 0 B，不是 8 B/元素。它的价值是形态：省掉一块输出缓冲，不是速度。

### 三维 / 四维绝对加速比

`make bench-vec34`，同进程交错 A/B，11 轮配对，`taskset -c 0`，单次运行；基线是逐元素库调用：

| 入口 | n = 8192 | n = 1M |
|---|---|---|
| `fv_vec3_add_soa` | 8.8x | 2.16x |
| `fv_vec3_scale_soa` | 13.9x | 1.43x |
| `fv_vec3_normalize_soa` | 10.1x | 2.11x |
| `fv_vec4_add_x2` | 35.2x | 2.53x |
| `fv_vec4_scale_x2` | 25.3x | 2.05x |
| `fv_vec4_normalize_soa` | 38.8x | 4.34x |

### `_fast` normalize

`make bench-fast`，同样协议；A = `_fast`，B = 精确 SoA，speedup = exact/fast：

| 入口 | n = 8192 | n = 1M |
|---|---|---|
| `fv_vec2_normalize_soa_fast` | 1.125x | 1.03x |
| `fv_vec3_normalize_soa_fast` | 1.054x | 1.04x |

8192 档的 A/A 地板为 0.995–0.997；1M 档为 0.89–1.08，所以 1M 的 1.03–1.04x 落在噪声带内，只能读作“没有变慢”，不能当作确定的加速。

### 读这些数字时要注意

- 加速比必须连同工作集大小一起报。N = 1<<20（24 MB）时三档实现收敛到 1.0–1.5x，因为都受限于本机 6–7 GB/s 的流带宽；历史上批量相对标量的比值更大，因为那个基线是 `FV_USE_EXTERN` 跨 TU 调用。P0 之后默认标量已内联，引用时必须说明基线。
- 历史收益主要来自消除跨 TU 调用边界，不是 intrinsics。同一份数学放进头文件内联后，`-O3 -march=native` 的自动向量化几乎追平手写 AVX2；P0 把这条结论用到了标量 API 自身。
- 这台机器长期有外部负载，同一份代码跨运行可差 20–25%。不同夹具的绝对值不要直接比，只看同一次运行内的比值。
- 三维 / 四维批量的绝对加速比现在由 `make bench-vec34` 测量（见上表）；引用时带上 harness 与工作集，不要当成单点精确倍数。
- AoS↔SoA 拆装由 `make bench-pack` 测量：往返成本在 n=8192 约 0.28–0.32 ns/元素、n=1M 约 4.2–4.4 ns/元素。AoS 输入下，`pack + add_soa + unpack` 在 8192 比 `fv_vec2_add_x4` 慢 2.5 倍、1M 慢 3.1 倍；normalize 在 8192 仍快 3.9 倍，但 1M 慢 2.3 倍。SoA 入口适合数据本来就是 SoA，或数据驻留 cache 且运算足够贵的情况。

判定标准见 `docs/acceptance-protocol.md`。原始判定表在 `docs/research/raw/ab_campaign_verdict.txt`，每次 `make run` 会覆盖，引用时请连同时刻。

### 小批量与流式写入的可移植性

- 批量 API 面向 **n >= 4**。外部 GLM 报告实测：n=1 约打平，n=2 约 0.59x，n=4 起 2x 以上；小于 4 个元素时调用/守卫开销就是全部成本，建议直接用标量 API。
- `_stream` 的符号不跨机可移植：只有工作集超过 L3 且机器受 RFO 限制时才赚；在带宽充裕、read-bound 的机器上会亏（报告机器 0.47–1.03x）。库不会自动切换，先用 `bench/membw.c` 实测。
- 编译档对 GLM 自动向量化的影响远大于 flashvec 的手写内核：报比值必须固定 `-O2` 或 `-O3 -march=native`，不要跨档比较单点。

## 构建与测试

```bash
make                         # build/libflashvec.a
make test-all                # 全量套件：1.97M 项检查（tests/count_checks.sh 可复算）

make test-header             # 严格编译档 + 双重包含
make test-inline             # 默认形态 vs FV_USE_EXTERN 形态逐字节对拍
make test-fma                # 消费者带 -march=native -mfma 时仍须位级一致
make test-inline-effects     # 内联效果普查
make test-stream             # _stream 与 _mask_inplace 的行为回归
make test-physics            # 四元数 + mat3，全形态（C/C++ inline/extern/pure-C/FMA）
make test-quat               # 四元数标量语义 + SoA 批量逐位一致
make test-mat3               # 3x3 标量语义 + SoA 批量逐位一致
make test-physics-c          # 纯 C11 翻译单元下的同一契约
make test-mat4-ops           # 通用 4x4 乘法 / 求逆
make test-fast               # _fast normalize 精度 + fv_atan_approx
make bench-physics           # 逐元素调用 vs SoA 批量，同进程 A/B
make bench-fast              # 精确 vs _fast normalize
make bench-vec34             # 三维 / 四维相对标量调用的绝对加速比
make bench-pack              # AoS↔SoA 拆装成本与 SoA 判断规则
make bench-mat3              # mat3 AVX2 乘/逆 vs 逐元素标量调用
make bench-scalar            # 标量 header 内联 vs FV_USE_EXTERN 库调用
make test-vec4               # 四维标量与批量
make test-vec34-inline       # 三维/四维批量形态对拍
make test-refactor           # 以上全部 + 导出符号清单 + 反汇编普查
```

### 编译选项

默认是 `-O3 -fno-math-errno -ffp-contract=off`。

`-fno-math-errno` 放弃「数学函数失败时设置 errno」。本库 `sqrt` 的参数恒非负，errno 不会被设置；加上它之后含 `sqrtf` 的循环才能向量化，数值逐位不变。

`-ffp-contract=off` 禁止把 `a + b*c` 收缩成 FMA。FMA 收缩会让同一表达式在不同编译单元或不同目标下差 1 ulp。
头文件里另有 `#pragma GCC optimize("fp-contract=off")` 作为第二道保险。

不使用 `-march=native`、`-mavx2`、`-ffast-math`。默认产物保持 x86-64 基线。
`-ffast-math` 不只是不必要：实测它在任何入口上都没有加速，还让 `normalize_soa` 在 L1/L2 规模上慢 15–17%，并改变 32.5% 的输出位模式。见 `docs/zh-CN/acceptance-protocol.md` §9 第 7 条。
AVX2 只在头文件里用函数级 `target("avx2")` 属性加运行时 `__builtin_cpu_supports` 派发启用。
非 AVX2 机器与所有尾长都走标量路径，结果相同。

## 设计取舍

### 用 C++ 编译器写 C

头文件是纯 C 接口，带 `extern "C"` 保护。实现用 C++ 编译器编译，但只使用 C 语言子集。这样将来想引入 C++ 特性时，不用换编译器，也不用改对外 ABI。

### 批量入口放在头文件里

跨 TU 调用边界曾是实测的一阶瓶颈，链接静态库的调用方无法消除它；把函数体搬进 header 就是修法。
所以实现只写一份，放在头文件里；`src/flashvec_batch.cpp` 保留同名 `extern` 符号，供 `FV_USE_EXTERN` 形态和已经编译好的旧调用方使用。

代价是头文件会在支持的目标上条件包含 `<immintrin.h>`。用纯 C 编译器或者定义 `FV_USE_EXTERN` 时，它仍然只是一个纯声明头文件。

### 头文件里的 pragma 会挡住内联

`target("avx2")` 管不住调用方命令行带来的 `-mfma`。实测在 `-march=native -mfma -ffp-contract=fast` 下，带 target 属性的内核照样会出现 `vfmadd`。
所以有收缩风险的入口必须用 pragma 罩住，而 GCC 不会把带不同优化选项的函数内联进调用方。

结果是这部分入口（`lerp`、`madd`、`dot`、`len`、`normalize`）在消费者的默认选项下每批多一次 TU 内调用，其余八个入口可以完整内联。
`make test-inline-effects` 把这个契约钉成了可复现的断言。

## 范围

库现在覆盖了向量类型及其标量 / 批量 API、四元数、3x3 与通用 4x4 矩阵、四元数插值与角度构造、`_fast` 近似 normalize 族、`fv_atan_approx`，以及融合 madd-normalize 的可选 AVX-512 路径。本文档不再保留待办清单。

以下内容是**有意不做**（不是待办）：线程 / 调度等引擎设施、分配器或 arena 所有权（库从不分配）、双精度批量族、以及每个入口的 `_fast` 版本（目前只有 normalize 与融合 madd-normalize）。这些会改变值语义 / 无分配契约，而不是补全它。

更多分析见 `docs/optimization-report.md`。哪些方向已经实测过、结论是什么，见 `docs/acceptance-protocol.md` 的 §9。