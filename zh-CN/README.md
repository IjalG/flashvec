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

## 标量 API

每个向量类型各有 18 个函数（`fv_vec3` 额外多一个 `cross`）：

`make` `add` `sub` `neg` `scale` `lerp` `min` `max` `clamp` `madd`
`dot` `len_sqr` `len` `dist_sqr` `dist` `normalize` `project` `near`

包围盒另有 `fv_aabb2_merge` 和 `fv_aabb2_overlaps`。

`min` / `max` / `clamp` 的 NaN 语义是 `(a < b) ? a : b`，不是 IEEE 的 `fminf`。这一点写进了头文件。

## 批量 API

共 33 个入口，按内存布局分五类。名字里的数字表示一个 256-bit 寄存器装几个元素：`fv_vec2` 装 4 个，叫 `_x4`；`fv_vec4` 装 2 个，叫 `_x2`。两者都能装满，没有 lane 浪费。

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

### 稀疏更新

| 函数 | 语义 |
|---|---|
| `add_x4_mask` | 掩码非零的元素做加法，掩码为零的元素写回 `a[i]` 原值 |
| `add_x4_mask_inplace` | 就地做加法。掩码为零的元素值不变；不保证没有被物理写回 |
| `add_x4_idx` | 只处理 `idx` 指定的下标，可重复、可乱序 |
| `scale_x4_idx` | 同上，做缩放 |

`_idx` 是标量循环，没有用 `vgatherdps`。Zen 2 上 gather 由微码实现，约 10 周期一条，不划算。

### 流式写入

`add_x4_stream` `sub_x4_stream` `scale_x4_stream` `add_soa_stream` `scale_soa_stream` `normalize_soa_stream`

用非临时存储写出结果，绕过 cache，也省掉 RFO 那一次读。结果不马上被读、且工作集远大于 L3 时才划算。这个选择必须由调用方做，库不会自动切换。

## 语义承诺

### 逐位一致

下面这些入口保证与对应标量 API 的浮点位模式完全相同，可以直接 `memcmp`：

- `add` / `sub` / `scale` / `lerp` / `madd` 的 `_x4`、`_x2`、`_soa`、`_mask`、`_idx` 变体
- `len_soa`
- `normalize_soa`
- `_stream` 系列与对应非 stream 版本

实现上只用 `vaddps`、`vsubps`、`vmulps`、`vdivps`、`vsqrtps`，不用 FMA，不用倒数近似，归约顺序与标量一致。

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

## 性能

测试机：AMD Ryzen 7 5700U，8 核 16 线程，AVX2 + FMA，无 AVX-512。

### 缓存边界

L3 是 **4 MiB per CCX**，两个 CCX 合计 8 MiB。单线程只能用到自己那 4 MiB。
这一点可以从 `/sys/devices/system/cpu/cpu0/cache/index3` 读到：`size=4096K`, `shared_cpu_list=0-7`。

按 `add` 的 24 B/元素算，拐点在 4 MiB ÷ 24 B ≈ **1.75×10⁵ 个活跃元素**。实测吻合：2^17（工作集 3 MB）是 38–57 GB/s，2^18（6 MB）掉到 5.6–8.8 GB/s。

### 相对标量调用的加速比

基线是**逐元素调用库里的标量函数**，那是跨 TU 调用、无法内联的写法。如果换成手写 C 循环并在同一编译单元内编译，GCC 会自动向量化，差距会小得多（`dot` 只剩 4 倍左右）。

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

### `_mask_inplace` 的实测

5 种密度 × 4 档规模共 20 个点，全部无法判定，增益在 −0.10 到 +0.00 之间。

原因在流量账上。就地写之前目标数组已经被读过，对应的 cache line 本来就在 cache 里，普通 store 这时不产生 RFO。所以「就地」省掉的是 0 B，不是 8 B/元素。它的价值是形态：省掉一块输出缓冲，不是速度。

### 读这些数字时要注意

- 加速比必须连同工作集大小一起报。N = 1<<20（24 MB）时三档全部收敛到 1.0–1.5x，因为都受限于本机 6–7 GB/s 的流带宽。
- 收益主要来自消除跨 TU 调用边界，不是 intrinsics。同一份数学放进头文件内联后，`-O3 -march=native` 的自动向量化几乎追平手写 AVX2。
- 这台机器长期有外部负载，同一份代码跨运行可差 20–25%。不同夹具的绝对值不要直接比，只看同一次运行内的比值。
- 三维与四维批量相对标量的加速比尚未测量，不要引用成数字。本轮只测了它们的形态差：65 个点里 5 通过、8 变慢、52 无法判定，读数散布在 0.97–1.93 之间。

判定标准见 `docs/acceptance-protocol.md`。原始判定表在 `docs/research/raw/ab_campaign_verdict.txt`，每次 `make run` 会覆盖，引用时请连同时刻。

## 构建与测试

```bash
make                         # build/libflashvec.a
make test-all                # 30 + 57 + 1098 条断言 + 两个 C 兼容性测试

make test-header             # 严格编译档 + 双重包含
make test-inline             # 默认形态 vs FV_USE_EXTERN 形态逐字节对拍
make test-fma                # 消费者带 -march=native -mfma 时仍须位级一致
make test-inline-effects     # 内联效果普查
make test-stream             # _stream 与 _mask_inplace 的行为回归
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

跨 TU 调用边界是实测的一阶瓶颈，而链接静态库的调用方无法消除它。
所以实现只写一份，放在头文件里；`src/flashvec_batch.cpp` 保留同名 `extern` 符号，供 `FV_USE_EXTERN` 形态和已经编译好的旧调用方使用。

代价是头文件会在支持的目标上条件包含 `<immintrin.h>`。用纯 C 编译器或者定义 `FV_USE_EXTERN` 时，它仍然只是一个纯声明头文件。

### 头文件里的 pragma 会挡住内联

`target("avx2")` 管不住调用方命令行带来的 `-mfma`。实测在 `-march=native -mfma -ffp-contract=fast` 下，带 target 属性的内核照样会出现 `vfmadd`。
所以有收缩风险的入口必须用 pragma 罩住，而 GCC 不会把带不同优化选项的函数内联进调用方。

结果是这部分入口（`lerp`、`madd`、`dot`、`len`、`normalize`）在消费者的默认选项下每批多一次 TU 内调用，其余八个入口可以完整内联。
`make test-inline-effects` 把这个契约钉成了可复现的断言。

## 还没做

- `_fast` 变体，例如 `rsqrtps` 加一次牛顿迭代的 normalize，约 2 ulp 误差，约快 10%。必须作为新函数提供，不能改默认语义。
- 近似三角函数。`atan` 多项式约比 glibc 快 2.7 倍，但必须写明输入域。
- 矩阵、四元数、角度与旋转。
- 三维与四维批量的绝对加速比测量。

已经做完、旧版这里写「留给以后」的：三维批量入口（走 SoA）、四维向量（`fv_vec4` 加 `_x2` 与 `_soa` 批量）。见上面的 API 节。

更多分析见 `docs/optimization-report.md`。哪些方向已经实测否决，见 `docs/acceptance-protocol.md` 的 §9。