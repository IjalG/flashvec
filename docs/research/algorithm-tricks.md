# flashvec 算法层优化调研（逼近 / 常数技巧 / 数值稳定性 / 物理仿真接口）

**测量环境**：AMD Ryzen 7 5700U（8C/16T，`/proc/cpuinfo` 确认 avx2 + fma，**无 avx512f**），g++ 12.4.0，
Linux 6.17.0-19-generic。

除标注「shipped」的行外，计时均为 `-std=c++17 -O2/-O3 -march=znver2`（各表内注明）。
数据取自仓库**现有** `build/libflashvec.a` 或 `src/flashvec.cpp`，未修改任何既有源码。

精度数字一律对 double 参考值取最大误差；计时为单进程 best-of-3。**笔记本变频：同一份手写代码两次测得 4.479 vs 5.583 ns/el（差 25%），故 <20% 的差异不作结论。**

---

## 1. 结论速览

**1 · API 形态比任何算术技巧都值钱。** 这里的形态指可否内联、是否批量。
同一份 AoS 数组、同一编译选项下：`fv_vec2_normalize` 按现有 API 调用是 3.18 ns/vec2，手写 branch-free SoA 加 `rsqrtps`+1NR 批量是 0.215 ns/vec2。
收益 **~15x**（N=2048，L1 驻留）。成本中：新增一批 `_soa`/`_batch` 入口，现有标量函数一行不改。

**2 · 不要手写 Quake rsqrt。** 实测 NR0 3.40e-2 → NR1 1.75e-3 → NR2 4.63e-6；而 `rsqrtss`+1NR 是 1.85e-7（3.1 ulp）。
两次牛顿迭代仍差 25 倍，属负收益。成本 0，直接放弃。

**3 · GCC 在 `-ffast-math` 下已经自动把 `1.0f/sqrtf(x)` 编成 `vrsqrtss`+1NR**（objdump 实测）。
但 `1.0f/x` 在任何选项下都是 `vdivss`，GCC 不生成 `rcpps`。
收益批量 1.8x、延迟 1.28x。成本低：新增 `fv_*_fast` 变体，不动默认语义。

**4 · 朴素 float 累加是精度黑洞，但多累加器同时赢速度与精度。** 1e6 次 `x+=v*dt` 朴素 float 相对误差 2.42e-3。
1e6 个数求和时 8 路累加器 **0.338 ns/add 且误差 1.0e-3**，朴素是 0.926 ns/add、误差 9.9e-3。
收益速度 2.7x **且** 精度 10x。成本低：引擎侧累加循环改写法，库侧只需保证 `dot/len_sqr` 可内联。

**5 · 换积分器比换精度类型值钱。** 谐振子 dt=1/60、6e4 步，半隐式欧拉最大能量漂移 float 5.522e-2 / double 5.521e-2，
换 double 只改善 0.02%；速度 Verlet 是 2.744e-3，**20x**。收益 20x。成本 0：库无关，属引擎选型。

**6 · atan 用 7 系数多项式加 π/2 恒等式。** 最大误差 1.16e-6 rad，3.71 ns/call；glibc `atanf` 是 10.15 ns/call。
收益 2.7x。成本低：一张系数表加一个除法，放在 `math` 扩展头里。

**7 · SIMD 下多项式反超查表。** sin 8-wide on [-π,π]：poly7+FMA 1.663 ns/val（err 1.19e-7），
gather-LUT4096 2.042 ns/val（err 5.36e-7），glibc `sinf` 4.164 ns/val。
收益 2.5x vs libm。成本中：需要参数归约；**角度必须保持 [-π,π]**，见 §2 的归约误差。

**8 · 零分支化有真收益但依赖数据。** `near` 型逐分量比较在 50% 命中率合成数据上从 5.349 降到 2.803 ns/call，即 1.9x。
但这是合成最坏情况，不能外推。收益 1.9x（最坏）。成本低：新增 `fv_vec2_near_bl` 或改默认实现。

---

## 2. 逼近与近似

### 2.1 快速平方根倒数：Quake trick 的历史与现状

Quake III 的原始实现含著名的 `0x5f3759df` 与 "what the fuck?" 注释，见
<https://raw.githubusercontent.com/id-Software/Quake-III-Arena/master/code/game/q_math.c>，第 557–565 行，已下载核对。

数值分析见 Lomont 2003 <https://www.lomont.org/papers/2003/InvSqrt.pdf>。

本机实测，20 万样本，2^-20..2^20，指标为相对误差：

| 方法 | 最大相对误差 | 折算 ulp |
|---|---|---|
| Quake `0x5f3759df`，0 次牛顿 | 3.396e-02 | 5.7e5 |
| Quake `0x5f3759df`，1 次牛顿 | 1.750e-03 | 2.9e4 |
| Quake `0x5f3759df`，2 次牛顿 | 4.630e-06 | 78 |
| `rsqrtss` 原始 | 2.577e-04 | 4323 |
| `rsqrtss` + 1 次牛顿 | **1.848e-07** | **3.1** |
| `rsqrtss` + 2 次牛顿 | 1.198e-07 | 2.0 |
| `1.0f/sqrtf(x)` | 8.875e-08 | 1.5 |
| `rcpss` 原始 / +1 次牛顿 / `1.0f/x` | 3.142e-04 / 1.705e-07 / 5.953e-08 | 5271 / 2.9 / 1.0 |

硬件近似精度符合 Intel 文档约束：SDM 条目（RSQRTPS）写 "|Relative Error| ≤ 1.5 * 2^-12"，等于 3.66e-4，
见 <https://www.felixcloutier.com/x86/rsqrtps>（已抓取正文核对）。实测 2.58e-4 在界内。

两个魔数 `0x5f3759df` / `0x5f375a86` 在本机 40 万样本上**看不到可测差别**：
3.438e-2 vs 3.437e-2，1 次牛顿后 1.752e-3 vs 1.751e-3。常量之争在最大误差指标上没有意义。

真正的结论是：**在 AVX2 上这两条路都被 `rsqrtss`+1NR 严格支配**，后者还少一整条整数/位运算依赖链。

关于 `-ffast-math` 的精度边界：本机 objdump 实测（`-O2`/`-O3 -march=znver2`），
`1.0f/sqrtf(x)` 被 GCC 编成 `vrsqrtss` 加一次牛顿（`vmulss`/`vfmadd213ss` 序列）。

单独开 `-freciprocal-math` **不触发**，仍是 `vsqrtss`+`vdivss`。
也就是说 `-ffast-math` 下的 `sqrtf` 依然精确，但 `1.0f/sqrtf(x)` 会静默掉到 ~3 ulp。`1.0f/x` 在所有测试选项下都是 `vdivss`。

微架构参考：<https://www.agner.org/optimize/instruction_tables.pdf>、<https://www.agner.org/optimize/optimizing_cpp.pdf>。
编译选项语义以 GCC 官方文档为准：<https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>。

### 2.2 三角函数与 atan 的多项式近似

真实做法参考 musl 的 `atanf`，即 atanhi/atanlo 分段表加按 |x| 区间选择，
见 <https://git.musl-libc.org/cgit/musl/tree/src/math/atanf.c>。

`sinf` 是象限归约加 `__sindf` 多项式，见 <https://git.musl-libc.org/cgit/musl/tree/src/math/sinf.c>
与 <https://git.musl-libc.org/cgit/musl/tree/src/math/__sindf.c>。

Cephes 的 atan/atan2/atanh 文件清单见 <https://netlib.org/cephes/>。

自测用 numpy 2.4.4 Chebyshev 拟合转幂基，**不是 Remez 极小极大，故误差略高于最优**。exact-arithmetic 拟合误差 1.156e-6，float 求值几乎不退化。

- **atan(x)，7 系数**：`t=x²` 的 6 次，|x|≤1，最大误差 1.161e-06 rad（float 求值）。
  配合 π/2 恒等式 `atan(x)=π/2−atan(1/x)`：|x|≤10 为 1.160e-06，|x|≤1000 为 1.051e-06 rad。
  **3.71 ns/call**，glibc `atanf` 10.15 ns/call，同机 `-O2 -march=znver2`。

- **对照：查表 4096 项加线性插值**，|x|≤1 时 7.283e-08 rad，比多项式准 16 倍，但只在标量路径可用（见下）。

- **sin，7 系数**：`t=x²` 的 6 次，在拟合域 [-π/4,π/4] 最大误差 5.916e-08。
  配 `n=round(x/π)` 归约加奇偶符号修正后，[-π,π] 为 1.192e-07。

**参数归约是 float 的真正瓶颈。** 在 [-1000π,1000π] 上用 float 归约误差 8.742e-05，用 double 归约 1.192e-07，相差 700 倍。
物理引擎把角度保持 [-π,π]（或改用旋转矩阵/四元数）就是这个原因。

**SIMD 8-wide 实测（-O3 -march=znver2，同一进程）**

| sin 实现（[-π,π]） | ns/value | 最大误差 |
|---|---|---|
| 多项式（poly7 + FMA，零内存访问） | **1.663** | 1.192e-07 |
| 查表 gather（`_mm256_i32gather_ps`，LUT4096 + 线性插值） | 2.042 | 5.364e-07 |
| glibc `sinf`（标量） | 4.164 | 参考 |

结论：**在 AVX2 上「查表换计算」这个竞赛/老游戏引擎套路会失效。** 随机索引需要 `vgatherdps`，恰恰是 Zen 2 的弱项；
多项式既能进 FMA 流水又不需要内存。LUT 只在「标量、无 SIMD、表小到常驻 L1」时才划算，1024 项等于 4 KB。

**→ 对 flashvec 的具体动作**：新建独立头 `include/flashvec_math.h`，不改 `flashvec.h`，
提供 `fv_atan_approx` / `fv_sin_approx` / `fv_cos_approx`。

文档写明最大误差量级（atan ~1.2e-6 rad，sin ~1.2e-7）与输入域约束（sin 输入须在 [-π,π] 或由调用方先归约）。
**不要**把近似三角函数塞进 `flashvec.h` 的向量接口。

---

## 3. 数值稳定性与累加：为什么「换个算法省掉精度问题」

补偿求和与精确点积的理论来源均已核对元数据：Kahan 1965 <https://doi.org/10.1145/363707.363723>；
Ogita–Rump–Oishi《Accurate Sum and Dot Product》<https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf>，
SIAM J. Sci. Comput. 2005，DOI 经 Crossref API 核对。

Higham《Accuracy and Stability of Numerical Algorithms》<https://doi.org/10.1137/1.9780898718027>，成对求和章节。

注：`doi.org` 对 curl 返回 403（反爬），浏览器可正常打开；我通过 `api.crossref.org` 核对了标题/期刊/年份。

### 3.1 实测：1e6 个 0.1f 求和（long double 参考 = 104857.6015625）

| 方法 | 结果 | 相对误差 | ns/add |
|---|---|---|---|
| 朴素 float | 105891.844 | 9.863e-03 | 0.926 |
| Kahan float | 104857.602 | 0（打印精度内） | 3.247 |
| Neumaier float | 104850.609 | 6.668e-05 | 1.116 |
| **8 路累加器（pairwise）** | 104748.945 | **1.036e-03** | **0.338** |
| double 累加 | 104857.602 | 0 | 0.873 |

- 朴素 float 累加 1e6 项丢掉约 **1%**。

- **8 路累加器同时赢速度与精度**：0.338 vs 0.926 ns/add，误差 10x 改善。原因是收益来自打断循环依赖链，而多路求和恰好也让误差变成「小部分和再合并」。这条值得在引擎侧直接用。

- Neumaier 在本例比 Kahan 差，因为补偿量 `c` 自身在 float 里继续累积误差。
  但教科书反例 `{1, 1e30, 1, -1e30}`（精确值 2）实测：朴素 **0**、Kahan **0**、Neumaier **2.000**，
  **Kahan 在大项存在时会丢失小项，Neumaier 不会**。选型建议：默认 Neumaier，其分支实测只比朴素慢 20%；
  需要极值时再上 Ogita–Rump–Oishi 的 Error-Free Transformation / Dot2。

- 点积（1e6 项，随机 [-1,1]）：朴素 5.219e-05，Kahan 5.076e-08，相差 1000x。

### 3.2 长时间积分：float vs double 的真实差距

**位置累加**：`x += v*dt`，v=1，dt=1/60，1e6 步，精确值 16666.667536。

| 累加方式 | 结果 | 相对误差 |
|---|---|---|
| float 朴素 | 16626.292969 | **2.422e-03** |
| float + Kahan | 16666.667969 | 2.597e-08 |
| double | 16666.667536 | 0 |

**谐振子积分器**：dt=1/60，w=2π，每 1000 步采样一次的能量最大相对漂移。

| 积分器 | 6e4 步 float | 6e4 步 double | 6e5 步 float | 6e5 步 double |
|---|---|---|---|---|
| 显式欧拉 | 1.13e71（发散） | 2.93e279（发散） | 发散 | 发散 |
| 半隐式欧拉 | 5.522e-02 | 5.521e-02 | 5.526e-02 | 5.525e-02 |
| 速度 Verlet | **2.744e-03** | 2.740e-03 | 2.767e-03 | 2.742e-03 |

**这是本节最重要的一句话**：辛积分器的漂移是算法性的、有界的，步数 ×10 后基本不变；
把 `fv_real` 从 float 换成 double 只改变 0.02%~0.09%，换积分器却改变 20 倍。

**所以「精度不够」的第一反应应该是换算法或换累加结构，而不是换标量类型。**

半隐式欧拉优于显式欧拉的工程共识见 <https://gafferongames.com/post/integration_basics/>。
积分器与稳定性的完整讨论见 Box2D 作者 GDC 笔记
<https://box2d.org/files/ErinCatto_NumericalIntegration_GDC2009.pdf>、
<https://gafferongames.com/post/fix_your_timestep/>。

**→ 对 flashvec 的具体动作**：

1. 保持 `fv_real=float` + `FV_EPS=1e-6` 不动。现有 30 条断言全绿，且实测精度损失与算法误差相比微不足道。

2. 新增 `fv_vec2_len_sqr_batch` / `fv_vec2_dot_batch`（SoA，8 路部分和）。这是全库唯一「既快又准」的改造点。

3. 新增 `fv_vec2_sum_kahan`/`fv_vec2_sum_neumaier`，或 `fv_kahan2` 累加器结构体 `{float s, c;}` 加 `fv_kahan2_add/result`。
   因为引擎的积分/求解循环需要它，而库现在的 API 里没有任何累加器形态。

4. **需要产品决策的量纲问题**：`fv_vec2_project` 用 `fv_abs(denom) <= FV_EPS` 判断「onto 接近零」。
   其中 `denom=len_sqr(onto)` 是**平方量**，而 `fv_vec2_normalize` 用 `len <= FV_EPS`，是线性量。
   同一个 `FV_EPS` 在 project 里等价于长度阈值 `sqrt(1e-6)=1e-3`，比 normalize 松 1000 倍。

   头文件注释却把 FV_EPS 描述为统一的「足够接近零」阈值。

---

## 4. 竞赛编程的常数优化套路：逐条迁移评估

**拆循环与预处理（先算不变量）**，判定 ✅ 高价值。
`fv_vec2_project` 里 `len_sqr(onto)` 对整批调用是常量：提供 `fv_vec2_project_inv(v, onto, inv_denom)`。
实测同类改造（消掉函数边界加冗余零检查）在 normalize 上是 3.18 → 2.51 ns/vec2，即 1.27x。

代价：新增 2–3 个 `_pre` 变体，API 面积变大。

**位运算替代取模/除法**，判定 ❌ 不适用。
向量库不做索引/哈希，这条留给引擎的 broad-phase（见 §5）。
唯一相关的是 `fv_abs`，其降级指令本轮**未核实**，不做结论。

**查表换计算**，判定 ⚠️ 反了。
标量划算：LUT4096 sin 3.27 vs glibc 5.80 ns/call，`-O2`。SIMD 反而更慢：2.042 vs 多项式 1.663 ns/val。表选项只在标量路径考虑。

**不变量外提**，判定 ✅ 中价值。同上第 1 条；另外 `dot` 在归一化前算的 `len_sqr` 不要算两遍，目前 `fv_vec2_normalize` 内部只算一次，是对的。

**手工展开 / 端口压力**，判定 ✅ 高价值。
真正的收益来自**打断依赖链**而非少几条指令：8 路累加器 0.338 ns/add，朴素 0.926。
深度展开会增加寄存器压力，SoA+AVX2 下 8 路是自然上限（`ymm` 8×float）。

**无分支 min/max / 条件移动**，判定 ✅ 有条件。
`near` 型比较在 50% 命中率合成数据上从 5.349 降到 2.803 ns/call，即 1.9x。**但这是合成最坏情况**，真实数据里距离比较常常高度可预测，收益会缩水。
另：x86 `minps/maxps` 与 C `fminf` 的 NaN 语义不同，本轮**未验证** SDM 的 MINPS 条目，采用前须查。

**SoA ＝ 把维度当循环**，判定 ✅✅ 最高价值。
这是全部调研里收益最大的一条：L1 驻留数据下 SoA 批量 normalize 0.215 ns/vec2，对比 AoS 按 API 调用的 3.18，约 **15x**。
同一批数据在 128 MB 工作集下两者收敛到 5.3 vs 7.0 ns/vec2，内存带宽主导，只剩 1.3x。

**先把问题规模降下来**，判定 ✅ 零成本。
比较距离用 `fv_vec2_dist_sqr`（已存在）而不是 `dist`；`project` 用平方和比较。

Box2D 同样提供 `b2DistanceSquared`，见
<https://github.com/erincatto/box2d/blob/v2.4.1/include/box2d/b2_math.h> 第 467 行，已核对。

**现状风险：`fv_vec2_len` 与 `fv_vec2_dist` 是最容易被误用的两个函数。**

**调用形态的完整实测（-O3 -march=znver2，同一份数据，ns/vec2，best-of-3）**

| 形态 | N=2048（64 KB，L1/L2） | N=64K（2 MB） | N=4M（128 MB） |
|---|---|---|---|
| shipped `fv_vec2_add`（无 LTO） | 1.624 | 1.697 | 8.745 |
| shipped `fv_vec2_normalize`（无 LTO） | 3.182 | 3.170 | 7.041 |
| shipped `fv_vec2_dot`（无 LTO） | 2.946 | 2.955 | 3.617 |
| 手写标量内联 `sqrtf`+`div` | 2.505 | 2.504 | 5.649 |
| **手写 SoA `rsqrtps`+1NR（8-wide，无分支）** | **0.215** | **0.230** | 5.376 |

- 结构性事实（objdump 硬证据）：链接仓库现有 `build/libflashvec.a` 的消费端里有 **5 处 `call fv_vec2_*`**；开 `-flto` 后 **0 处**。
  库是静态库且头文件里没有 inline 定义，所以**消费端不内联就必然一函数一调用**。

- LTO 还能让无分支的可归约循环吃到自动向量化：`len_sqr` 类循环实测 2.778 → 1.307 ns/el，2.1x。
  但本轮的 LTO 计时表受「循环不变式外提」污染（我的测试外层重复同一批数据），
  **故 LTO 的绝对数字不作结论，只保留「5→0 调用」这一结构性事实**。

- 精度代价：在 4M 个真实 vec2 上，`rsqrtps`+1NR 路径与库现有 normalize 路径的**逐分量最大偏差 1.192e-07**，约 2 ulp of 1.0。

**→ 对 flashvec 的具体动作**：

1. 新增 `fv_vec2_normalize_soa(const float* x, const float* y, float* ox, float* oy, int n)`，`fv_vec3` 同理，y/z 用第三个数组。
   AVX2 8-wide，`rsqrtps`+1NR+FMA；**必须无分支**，零长度元素用比较加 blend 归零，不能 early-return，否则向量化失败。
   收益 15x（缓存驻留）/1.3x（内存带宽），代价约 40 行。

2. `fv_vec2_normalize_fast`（标量，`_mm_rsqrt_ss`+`1NR`）：只在调用方显式要求时用。
   默认 `fv_vec2_normalize` **不改**，其误差 3 ulp vs 1.5 ulp，且现有测试断言按当前语义写的。

3. 提供 `fv_vec2_near_bl`（max 形式，无分支）而不是改默认 `near`。

4. `Makefile` **只追加**目标：`-O3 -march=znver2` 的 `libflashvec_fast`（或 `CXXFLAGS_FAST`）加 `-flto`。
   并在消费端文档里写清「要吃到顺序无关的向量化，消费端必须 `-flto` 或改用 `_soa` 入口」。
   现有默认目标（`-O2`、无 `-march`）实测产物里 **`vfmadd` = 0**，即通用 x86-64 无 FMA，等于白白放弃一半乘法吞吐。

---

## 5. 物理仿真特有的算法级优化：向量库要留的口子

这些不是向量库该实现的东西，但**库的形态会决定引擎能不能用它们**。

1. **固定时间步 + 子步（substep）。** dt 恒定是数值稳定性的前提，见 <https://gafferongames.com/post/fix_your_timestep/>。
   库要留的口子：批量积分入口需要接收**标量 dt**。标量广播进 SIMD 是免费操作，`_mm256_set1_ps` 即可；
   不要设计成「每元素一个 dt」。

2. **休眠刚体跳过积分。** Box2D 用每体 `m_sleepTime` 加岛级判定，见
   <https://github.com/erincatto/box2d/blob/v2.4.1/include/box2d/b2_body.h> 第 330/469/638 行已核对；
   岛求解见 <https://github.com/erincatto/box2d/blob/v2.4.1/src/dynamics/b2_island.cpp>。

   **关键的接口约束**：SoA 批量函数如果只能处理「连续的前 n 个」，休眠体就无法跳过。
   必须提供**索引数组或位掩码**形式：`fv_vec2_add_gather_soa(dst_x, dst_y, src_x, src_y, const int* idx, int n)`。
   索引形式天然可 SIMD 化；掩码形式在 AVX2 上需要 `blendvps`。

3. **Broad-phase。** SAP 的经典讲义 <https://www.cs.cmu.edu/~baraff/sigcourse/notesd2.pdf>；
   动态 AABB 树 <https://github.com/erincatto/box2d/blob/v2.4.1/src/collision/b2_dynamic_tree.cpp>，
   GDC 讲解 <https://box2d.org/files/ErinCatto_DynamicBVH_Full.pdf>；
   Bullet 的 `btDbvtBroadphase`
   <https://github.com/bulletphysics/bullet3/blob/master/src/BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>；
   空间哈希 <https://matthias-research.github.io/pages/publications/tetraederCollision.pdf>。

   **当前最明显的缺口：库里没有 `min`/`max`/`clamp`**，而 AABB merge、SAP 排序键、树更新全靠它们。
   应新增 `fv_vec2_min/max/clamp`、`fv_aabb2_merge/overlaps`；AABB 只作为两个 vec2 的参数对，不引入新结构体。

4. **约束求解的顺序与迭代次数。** Gauss–Seidel 类求解器对顺序敏感，迭代次数换精度，见
   <https://box2d.org/files/ErinCatto_IterativeDynamics_GDC2005.pdf>、
   <https://box2d.org/files/ErinCatto_SoftConstraints_GDC2011.pdf>、
   <https://box2d.org/files/ErinCatto_UnderstandingConstraints_GDC2014.pdf>。

   库要留的口子：`fv_vec2_madd(dst, v, s)`（saxpy，冲量累加 `v += Jᵀλ` 的最小原语）加批量形式；
   以及**不要**把批量 API 设计成全原地更新，求解需要「先收集所有冲量再统一应用」的两遍结构。

5. **对称性与守恒（动量/能量修正）。** PBD/PXPBD 把稳定性放在位置层投影上，而不是靠精度，
   见 <https://matthias-research.github.io/pages/publications/posBasedDyn.pdf>。
   库要留的口子：位置/速度可分离的 SoA 布局，以及不假设「积分与修正必须在一次调用内完成」。

6. **对称性的免费红利。** 力成对累加（Newton 第三定律）时，`fv_vec3_add/sub/neg` 是最热的三条路径。
   它们目前每次都是独立函数调用，证据是 §4 的 5 处 `call`。这是引擎侧最该被批量化的一处。

**→ 对 flashvec 的具体动作**：新增 `fv_vec2_madd`、`fv_vec2_min/max/clamp`、`fv_aabb2_merge/overlaps`、
索引/掩码形式的 SoA 批量入口。这些全部是「新函数」，不动既有符号与语义。

---

## 6. 明确劝退清单（在这个规模上不划算）

1. **为 3 维向量做编译期全展开的模板表达式。** 实测证明函数边界只值 1.27x，3.18 → 2.51 ns/vec2，而真正的 15x 来自 SoA 批量。
   模板会摧毁「纯 C 头 + `extern "C"`」的设计前提（`flashvec.h` 第 13–15 行），并让 C 消费者（`tests/test_c_api.c`）无法使用。

2. **过早引入 SIMD 复数。** 现在连 `fv_vec3` 都还没有批量入口；复数既非第一阶段需求（物理底座），收益面也窄。

3. **为了「省一次开方」引入精度不确定性。** 应该做成 `_fast` 变体加 `_soa` 批量（精度实测 1.19e-07 偏差），
   而不是改 `fv_vec2_normalize` 的默认行为。默认路径保持 1.5 ulp 可预测，是物理引擎调试期的基本盘。

4. **把 `fv_real` 换成 double 求精度。** 实测能量漂移只改善 0.02%（半隐式欧拉），却让 SoA 批量吞吐减半（**未实测**，属推断）。

5. **把 `-ffast-math` 写进库的默认编译选项。** GCC 文档明确它包含 `-ffinite-math-only`，即假设无 NaN/Inf，
   见 <https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>。

   而 `fv_vec2_normalize`/`project` 的现有契约（「零向量返回 (0,0)，不会放大出 NaN」）正依赖 NaN/比较语义可预测。
   要快就给 `_fast` 变体，不要改全局开关。

6. **给 `fv_vec2_near` 加 SIMD。** 单次只处理 2 个 float，SIMD 无意义；正确做法是提供一个批量 `near_mask`（SoA + 掩码输出）。

7. **在向量库里塞表（sin/atan LUT）或 BVH。** 职责越界。表属于 `flashvec_math`（可选头），BVH 属于引擎。

8. **继续放大 `FV_EPS` 的量纲混用**（§3.2 第 4 条）。这是「省一次比较」的反面教材：
   为了少算一次 `sqrt` 而在平方量上复用线性阈值，会在 `project` 的零向量边界上引入 1000 倍宽的静默区间。

---

### 附：本轮证据的可复现命令

```bash
# 调用形态 / 是否内联（结构性证据）
g++ -std=c++17 -O3 -march=znver2 -Iinclude /tmp/fvbench/consumer3.cpp build/libflashvec.a -o /tmp/c3 && objdump -d --no-show-raw-insn /tmp/c3 | grep -c 'call.*<fv_vec2_'   # -> 5
g++ -std=c++17 -O3 -march=znver2 -flto -Iinclude -c src/flashvec.cpp -o /tmp/fv_lto.o \
  && g++ -std=c++17 -O3 -march=znver2 -flto -Iinclude /tmp/fvbench/consumer3.cpp /tmp/fv_lto.o -o /tmp/c3lto && objdump -d --no-show-raw-insn /tmp/c3lto | grep -c 'call.*<fv_vec2_'   # -> 0
# -ffast-math 下 1.0f/sqrtf 的降级
g++ -std=c++17 -O2 -ffast-math -mavx2 -mfma -c t.cpp -o t.o && objdump -d --no-show-raw-insn t.o | grep -E 'vrsqrtss|vsqrtss|vdivss'
```
