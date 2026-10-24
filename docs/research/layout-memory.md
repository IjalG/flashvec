# flashvec · 数据布局 / 内存层次 / 批处理 API 形态

> 调研员 C2。本笔记只做**阅读 + 判断**，不新增实验。数字全部来自 [measured.md](measured.md)（三档实测）、
> [simd-hardware.md](simd-hardware.md)（SIMD/ABI）、[algorithm-tricks.md](algorithm-tricks.md)（休眠体约束）与本机 [raw/](raw/)。
> 无标记 = 有实测依据；标「**未验证**」= 本组未测、未逐字核对上游源码/文档。

> 本机：AMD Ryzen 7 5700U（Zen 2, 8C/16T）· L1d 32 KiB/核 · L2 512 KiB/核 ·
> **L3 4 MiB/CCX × 2 CCX = 8 MiB 芯片总量，单线程可见 4 MiB** · 行 64 B · AVX2+FMA，无 AVX-512。

---

## 0. 结论速览（7 条，每条带改造量）

1. **先加批量入口，不要先改结构体。** L2 驻留（N=8192）实测 add x6.54、scale x13.23、dot x8.98、
   len x9.63、normalize x8.41、collision x14.51；而 N=1<<20 时三档全贴在 5.3–6.8 ns/elem（measured.md §0 条 1–2、§3.1）。
   **改造量**：新增 1 个编译单元 + 6–10 个符号 ≈ 150–250 行，现有标量 API 一行不改。

2. **`fv_vec2` 保持 8 B，不要对齐到 16 B。** 8 B 走 1 个 xmm，指令就是 `addps;ret`。
   且**一个 ymm 正好装 4 个 fv_vec2**，连续 4 个元素无需任何重排（simd-hardware.md §1 条 2/4）。
   padding 到 16 B 会让 x4 入口一个寄存器只装 2 个。**改造量**：0，现状即最优。

3. **AoS 4 路交错入口（`_x4`）优先级高于纯 SoA 入口。**
   它把「D 实测 8–14x」与「物理引擎调用点本来就是 AoS 粒子数组」这对矛盾直接消解：调用方**零改动**，
   库里用 `_mm256_loadu_ps` 把连续 4 个 fv_vec2 当 8 个 float 用。

   **改造量**：add/sub/scale/lerp 各 1 个函数 ≈ 60 行；dot/len/normalize 需额外一次「偶数列打包」，成本更高。

4. **纯 SoA 入口留给归约型运算（dot/len/normalize），且必须承认调用方成本未知。**
   AoS↔SoA 拆装成本**本组未实测**（measured.md §1.2 自述），引擎若内部是 AoS，转换会抵消一部分收益。
   **改造量**：库里 +5 个符号；调用方需拆装代码（0 行～重写数据组织）。

5. **每个批量入口必须配「索引 / 掩码」变体，否则休眠体跳过与 broad-phase 稀疏更新用不了。**
   「只处理连续前 n 个」是硬约束：Box2D 用 per-body `m_sleepTime` 加岛级判定跳过积分，见
   [b2_body.h](https://github.com/erincatto/box2d/blob/v2.4.1/include/box2d/b2_body.h) /
   [b2_island.cpp](https://github.com/erincatto/box2d/blob/v2.4.1/src/dynamics/b2_island.cpp)。

   其中行号已由 algorithm-tricks.md §5 条 2 核对。
   **改造量**：每个运算 +1–2 个变体，≈ +40% 符号数。

6. **对象池 / arena 现在不做。** 库契约是「值语义、按值传、无分配」（`include/flashvec.h` L8）。
   加池会同时破坏 C 消费者、无全局状态、跨 TU 的 malloc 配对三件事。**改造量**：0。

7. **ABI 不是瓶颈，不要拿「减少拷贝」当改签名的理由。**
   8 B 走 1 个 xmm、12 B 走 2 个、16 B 仍走寄存器（measured.md §1.5、§3.6 B1/B7 实测）。真正的代价在**跨 TU 调用边界**。

---

## 1. AoS / SoA / AoSoA 的真实权衡

### 1.1 三种布局在本机上的具体含义

| | AoS（现状） | SoA（分量数组） | AoSoA（分块交错，4/8 个一组） |
|---|---|---|---|
| 内存形态 | `fv_vec2 p[n]`，步长 8 B | `float px[n], py[n]` | `struct { float x[4], y[4]; } blk[m]` |
| SIMD 打包 | 4 个 vec2/ymm，**零重排** | 8 个分量/ymm，**零重排** | 4 或 8 个 vec2/ymm，零重排 |
| 缓存利用率 | 顺序访问满分；只取 x 时也拉进 y（浪费 50%） | 只取所需分量，理论最优 | 折中：块内顺序、块间跳 |
| 可读性 / 调试 | 最好：`p[i].x` 一处可见 | 差：同一粒子散在 2–4 个数组，gdb 里要跳 | 最差：两层索引 |
| C 生态兼容 | 最好（POD、可 memcpy/memcmp） | 好（就是 float 数组） | 好，但块尺寸要进 ABI |
| 迁移成本 | 0 | 调用方重写数据组织 | 调用方重写 + 块尺寸约定 |

**AoSoA 的前提是「块小到进 L1、大到装满一个 ymm」**：4 个 vec2 = 32 B，8 个 = 64 B，正好 1 行。

[Agner Fog《Optimizing software in C++》](https://www.agner.org/optimize/optimizing_cpp.pdf) 系统讨论过 AoS/SoA 取舍，
measured.md §7 已验证该链接返回 200。

[Richard Fabian《Data-Oriented Design》](https://www.dataorienteddesign.com/dodbook/) 讲布局与「按消费者拆轴」，
**未验证**：未核对章节号与可达性。

### 1.2 真实案例（谁真的这么干了）

- **ECS（EnTT），按「组件类型」拆轴的极端 SoA。** 每组件类型一个 packed array，实体只持索引
  （<https://github.com/skypjack/entt>，**未验证**行号）。教训不是「SoA 更快」，而是**遍历谁就把谁排连续**。

- **Bullet，AoS 求解器加 4 分量 SIMD 向量。** `btSequentialImpulseConstraintSolver` 持有
  `btSolverBody`/`btSolverConstraint` 的**连续数组**（AoS），SIMD 藏在 `btVector3`（SSE 4-wide）里。
  求解器源码
  <https://github.com/bulletphysics/bullet3/blob/master/src/BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.cpp>

  **未验证**未逐行核对。它选的是「**AoS 数组 + 向量类型本身 4-wide**」，正是本笔记给 flashvec 的路线。

- **PhysX 4，按消费者拆轴的折中。** 求解体拆成 `PxSolverBody` 与 `PxSolverBodyData` 两条并行数组
  （[PhysX](https://github.com/NVIDIA-Omniverse/PhysX)，**未验证**未核对文件与版本）。

  最有价值的一条是：**不是全局 SoA，而是「热循环里被谁读，就为谁单开一条数组」**，对 flashvec 的形态就是「AoS 主数组 + 可选 pos/vel 影子数组」。

- **simdjson，逐阶段 SoA，阶段间换布局。** Stage 1（结构扫描）对 64 B 块并行产出多个独立位掩码/结构索引流，Stage 2 再按需重放，
  **同一份字节流，不同阶段用不同布局**（论文 <https://arxiv.org/abs/1902.08318>、<https://github.com/simdjson/simdjson>，**未验证**可达性）。

  教训：布局不是全局属性，是**阶段本地**属性。

- **DirectXMath 坚持 AoS 的理由。** `using XMVECTOR = __m128`，见
  <https://github.com/microsoft/DirectXMath/blob/main/Inc/DirectXMath.h#L402>，simd-hardware.md 已引用该行。
  「向量类型 = 一个寄存器」，强制 4 分量（vec3 也吃 16 B），配 MSVC `__vectorcall`。

  其收益是「一个 C++ 变量就是一个 SIMD 寄存器」、代码与数学公式 1:1；代价是 25% 带宽、纯 C++、绑死 4 分量
  （[编程指南](https://learn.microsoft.com/en-us/windows/win32/dxmath/pg-xnamath-programming)，**未验证**未逐字核对布局表述）。
  → flashvec 是**纯 C 的 8 B vec2 库**，结论恰好相反：不要沿用它的 AoS 理由，那理由建立在「类型即寄存器」上。

### 1.3 关键张力：库的契约 vs 实测的 8–14x

- **契约侧**：值语义、按值传、无分配、纯 C 头（`flashvec.h` L1–9）。用户写 `for (i) o[i] = fv_vec2_add(a[i], b[i]);`，
  每元素一次跨 TU 调用，GCC 明确拒绝向量化（simd-hardware.md §1 条 1：`missed: statement clobbers memory`）。

- **收益侧**：D / A 实测说 SoA 批量能拿 8–14x（L2 驻留），但**调用点本来就是 AoS 粒子数组**。

- **解法（本笔记的核心主张）：不要二选一，把「AoS 4 路交错」做成一档，SoA 做第二档。**

```c
/* 一档：AoS 原样传入，一个 ymm 装 4 个 fv_vec2，零重排、调用方零改动 */
void  fv_vec2_add_x4  (const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o, size_t n);
void  fv_vec2_scale_x4(const fv_vec2 *a, fv_real s, fv_vec2 *o, size_t n);
/* 二档：归约型运算，分量天然在连续数组里 */
float fv_vec2_dot_soa (const float *ax, const float *ay,
                       const float *bx, const float *by, size_t n);
```

`_x4` 的语义与现有 API **完全一致**，就是同一个 fv_vec2 加法、只是 4 个一批，所以它不是新契约，是**新入口**。
调用方的 `particles[]` 不需要拆分。
**矛盾的解法就是：布局不改，只把「每元素一次 call」换成「每 4 个元素一次 call」。**

### 1.4 「什么时候该切」的判据

| 判据 | 阈值 | 依据 |
|---|---|---|
| **活跃元素数** | ≤ **1.75×10`5`** 时 SIMD 收益能兑现，该值对应工作集 ≈ 单线程可见的 4 MiB L3；≥ 1×10`6` 回归带宽墙 | measured.md §0 条 1、§2.1 粗算 |
| **同构运算密度** | 同一批元素在一个热循环里被 ≥2 个同构运算处理才值得拆布局；只有 1 个运算时拆装成本必亏 | 定性；collide 微内核是典型「≥3 个」 |
| **是否被多个消费者读** | 位置同时被积分器 / broad-phase / 渲染读 → 值得按消费者拆轴（PhysX 式） | §1.2 |
| **稀疏度** | 休眠/剔除比例 > ~30% 时「连续前 n 个」的收益被空转消耗掉 → 必须先有索引/掩码变体 | algorithm-tricks.md §5 条 2 |

第一行的阈值已按 sysfs 实测更正：旧值 3×10`5` 是按芯片总量 8 MiB 算的，高了一倍。

---

## 2. 内存层次与本机数字

### 2.1 为什么 N=1<<20 时「代码怎么写都一样」

本机 `lscpu` 原文（[raw/env.txt](raw/env.txt)）：`L1d cache: 256 KiB (8 instances)`（32 KiB/核）、
`L2 cache: 4 MiB (8 instances)`（512 KiB/核）、`L3 cache: 8 MiB (2 instances)`。

**「2 instances」就是 2 个 CCX，每个 4 MiB。** 这不只是 lscpu 的计数口径：
`/sys/devices/system/cpu/cpu0/cache/index3` 直接给出 `level=3, size=4096K, shared_cpu_list=0-7`，
即同一 CCX 的 4 个物理核 / 8 个线程共享 4 MiB。

- **8 MiB 是芯片总量，单线程看得见的是 4 MiB。** 旧版本笔记里「L3 8 MiB 共享」这个表述是错的，它把总量当成了单线程可用量。

- 24 MB 工作集仍远大于任一 L3 实例，所以本节「三档都受带宽上限约束」的结论不变；但**拐点元素数必须按 4 MiB 算**，即 1.75×10`5`（见下条），旧值 3×10`5` 高了一倍。

`bench/membw.c` 实测（[raw/membw.txt](raw/membw.txt)，taskset -c 2，每项 5 次取最小）：

```
size         working_MB    triad_GBs  triad16_GBs    write_GBs     read_GBs   memcpy_GBs
N=1048576          12.0         6.91         9.21        12.33        12.22        17.19
N=4194304          48.0         6.32         8.43         4.24         8.77         9.71
N=8388608          96.0         6.07         8.09         3.78         7.63         8.91
```

**粗算，本笔记推导，非实测**：add 每元素读 2×8 B、写 1×8 B，合计 24 B 逻辑访存；N=1<<20 工作集 3×8 MB = 24 MB，
按 triad 6.9 GB/s 得 **≈3.5 ns/elem 的理论下限**，而实测 5.3–6.8 ns/elem，即有效 3.5–4.5 GB/s。
所以 measured.md §0 条 1 说，在 1<<20 这个规模上，自动向量化和手写 AVX2 都不产生收益。

三点必须说清：

- **拐点在 1.75×10`5` 元素量级**：单线程可见的 4 MiB ÷ 24 B/元素 = 174 762；若误按芯片总量 8 MiB 会算成 3×10`5`，高了一倍。
  这条有直接实测：同进程 A/B 夹具在 **2^17（工作集 3 MB）38–57 GB/s，2^18（6 MB）5.6–8.8 GB/s**，断崖正好压在这条线上。
  引擎每帧只处理几千～几万活跃体时正好落在「SIMD 收益兑现」一侧；超过百万回到带宽墙，此时唯一有效手段是**减少访存量**（降精度、跳过休眠体、分块），不是换指令。

- **membw 两次独立运行有 10–15% 漂移**：raw/membw.txt 的 write=12.33 对比 measured.md §3.3 粘贴的 14.05，triad 6.91 对比 7.01。
  本笔记引用 **raw/** 的值；10% 以内的差距在这台机器上不该当结论。

- **measured.md 的 8–14x 是 L2 驻留数字，不是 L1 峰值**：N=8192 的 SoA 形态下 4 输入 + 2 输出分量数组约 192 KB，
  超出 32 KB L1d、落在 512 KB L2 内。这个区间恰好是真实引擎的活跃集大小，所以结论**可迁移**。

### 2.2 结构体大小的影响（实测过的部分）

| 类型 | sizeof | 一条 64 B 行装 | ABI | 判断 |
|---|---|---|---|---|
| `fv_vec2` | 8 B | **8 个，整** | 1 个 xmm（`addps;ret`） | **最优**：4 个进一个 ymm，零重排 |
| `fv_vec3` | 12 B | 5 个 + 4 B 残渣 | **2 个 xmm，第二个只用 4 B**（§3.6 B7） | 步长不规则、会跨行；但**改布局 = 破坏 ABI**，不要改 |
| `fv_vec4` | 16 B | 4 个，整 | 1 个 ymm 装 2 个 | 只用于新功能（mat4/quat），不要替换 vec3 |

- `fv_vec3` 若 padding 到 16 B，每记录多 4 B 带宽，即 +33%，换来「1 个 xmm、无跨行」，正是 DirectXMath 的取舍（§1.2）。
  flashvec 保持 12 B：measured.md §1.5 的实测结论是「按值 ABI 本身几乎免费」，没有非改不可的理由，
  而改布局会让所有已编译使用方 ABI 失效。

- **三维扩展的正确做法**：新增 `fv_vec4`，而不是让 `fv_vec3` 变形。

- 唯一硬建议：**`_x4` 入口内部用 `_mm256_loadu_ps`，不要要求 32 B 对齐。**
  要求对齐会把「调用方零改动」这唯一优势也会失去；Zen 2 的 unaligned load 只要不跨页就与对齐同速，**未验证**：本组未做微基准。

### 2.3 false sharing（多线程求解器里的表现 + 判据）

- **表现**：两线程写同一 64 B 行内**不同**元素时，行在核间 ping-pong，吞吐可掉到单线程的几分之一。
  写共享的加速结构（broad-phase 计数、每线程累加和、`std::atomic` 计数）是重灾区。

- **判据**：① 线程分区边界**按 64 B 对齐**、每块 ≥ 64 B ⇒ 消失；② 按 `i % T` 轮转分配元素**必然**发生，同一行里混着不同线程的元素；
  ③ 每线程累加和用 `alignas(64)` 独立槽位或线程局部缓冲加一次归约。

- **对休眠体跳过的额外影响**：稀疏写让访问变成随机行命中、行利用率下降。这也是「掩码变体」优于「索引变体」的原因之一：掩码整行读改写，行利用率保持。

---

## 3. 批处理 API 的三条路线

### 3.1 对照表

三条路线各自的签名：

- (a) `for (i) o[i] = fv_vec2_add(a[i], b[i]);`

- (b) `fv_vec2_add_soa(const float*ax,const float*ay,const float*bx,const float*by,float*ox,float*oy,size_t n)`

- (c) `fv_vec2_add_x4(const fv_vec2*a,const fv_vec2*b,fv_vec2*o,size_t n)`

| 维度 | (a) 纯标量自组织 | (b) SoA 分量数组入口 | (c) AoS 4 路交错入口 |
|---|---|---|---|
| ABI | 不变（8 B → 1 xmm） | 6 指针 + n，全走整数寄存器；无结构体分类问题 | 3 指针 + n，与现有头文件风格一致 |
| 可向量化性 | **0**：GCC 拒绝，报 `clobbers memory`，每元素一次 call | **6.5–13.2x**（实测 N=8192）；归约型最自然 | add/sub/scale/lerp **≈满速** |
| | | | 1 条 `add_ps` 处理 4 个 vec2；dot/len/normalize 需一次偶数列打包（`vpermd`/`permutevar8x32_ps`），**未实测** |
| 调用方改造成本 | 大（自己拆分量数组），收益 0 | **高**：拆分/装配 | **0**：现有 AoS 粒子数组直接传 |
| | | AoS↔SoA 转换成本**未实测**（measured.md §1.2 自述） | |
| `n%4`/`n%8` 尾巴 | 无 | 8 元素主循环 + 0–7 个标量收尾 | 4 元素主循环 + 0–3 个标量收尾 |
| | | 收尾可用 `_mm256_maskload_ps` | 收尾可用 128-bit `_mm_loadu_ps` 加标量 |
| 休眠体/稀疏更新 | 调用方自由，但退化为标量 | 索引变体（gather）、掩码变体（`blendvps`） | 同上，且**掩码变体更划算**（整行读改写，行利用率不塌） |
| 尾巴期正确性风险 | 无 | 高：6 个指针的别名/重叠要写清 | 中：输入输出可同址，in-place 天然成立 |

### 3.2 休眠体跳过是硬约束（B 与 D 都提到的一点）

「只处理连续前 n 个」的接口**无法**支持两件事：① 休眠体跳过（Box2D 的 `m_sleepTime` 加岛级判定，§0 条 5）；
② broad-phase 稀疏更新（只处理 `moved[]` 里的索引）。所以每个批量入口都要成对给出：

```c
/* 索引形式：真的跳过，代价是随机访存 */
void fv_vec2_add_x4_idx (const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o,
                         const uint32_t *idx, size_t m);
/* 掩码形式：整行读改写，掩码外写回原值 */
void fv_vec2_add_x4_mask(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o,
                         const uint8_t *mask, size_t n);
```

- **掩码语义必须写进头文件**：mask 位为 0 的 lane 输出 = 输入原值（`_mm256_blendv_ps`），**不是**「跳过存储」。ymm 只能整条写回，跳过存储做不到，且会留半行脏数据。

- **不要在 Zen 2 上用 `_mm256_i32gather_ps` 然后假设它快**：Zen 2 的 gather 是微码循环，通常不如标量逐元素。**未验证**：本组未实测，属文献推断。

- **优先级**：先做 (c) 的 add/sub/scale/lerp（调用方 0 成本）→ 再做 (b) 的 dot/len/normalize（归约型，SoA 更自然）→
  最后才劝调用方迁到 SoA。**不要反过来：先改别人的数据布局，再给库加函数。**

---

## 4. 分配与生命周期

- **16 B 结构体按值返回走寄存器，已实测，不重复验证**（measured.md §1.5、§3.6 B1/B7）。
  推论：**不要**为「避免结构体拷贝」把签名改成 `void add(const fv_vec2*, const fv_vec2*, fv_vec2*)`。
  measured.md §4.6 明确说档1 慢的是调用开销加标量化，不是 ABI 拷贝。

- **库不分配是契约，不是偶然**（`flashvec.h` L8）。加内存池会同时破坏三件事：
  C 消费者（`tests/test_c_api.c` 用纯 C 编译）、「无全局状态」承诺、跨 TU 的 malloc/free 配对（`extern "C"` 边界上最易出事）。

- **对象池 / arena 什么时候才值得**：三条缺一不可，且全部属于**调用方**。
  ① profile 显示分配器占帧时间 > 5%，**未验证**：本组无帧级 profile，这条是判据不是结论；
  ② 分配次数与活跃体数同阶，即每帧每体一次 malloc；
  ③ 需要 **64 B 对齐**分量数组时用 `aligned_alloc`/`posix_memalign`，C11 标准函数，而不是自建池。

  参考 <https://en.cppreference.com/w/c/memory/aligned_alloc>，**未验证**可达性。

- **明确结论：不要在没有性能数据前给 flashvec 加内存池。**
  给一个 8 B POD 库加池，等于把「值语义」换成「生命周期管理」：收益上限是一次 malloc 的 ~50 ns，风险是所有权 bug 与 C/C++ 混编裂缝。

---

## 5. 劝退清单

1. **不要为了对齐把 `fv_vec2` 填充到 16 B。** 8 B 正好 4 个进一个 ymm，是这台机器上的最优；填充后 x4 入口每个寄存器只装 2 个，收益减半。

2. **不要给单元素 API 套 AVX2 分派**（`__builtin_cpu_supports` 分支）。单元素调用永远被 call 开销主导（measured.md §1.6）。

3. **不要在没有 AoS↔SoA 转换成本数据前，把库改成纯 SoA 契约。** 转换成本本组未实测，是当前最大的未知数。

4. **不要用 `_mm256_i32gather_ps` 实现索引变体并宣称加速。** Zen 2 的 gather 是微码循环，未实测，属文献推断。

5. **不要引入模板表达式来「解决」布局问题。** 它摧毁「纯 C 头 + `extern C`」前提（algorithm-tricks.md §6 条 1 已劝退），而收益在批量入口，不在模板。

6. **不要给 `fv_vec3` 加 padding 或改 `fv_vec2`/`fv_vec3` 的字段顺序。** 12 B→16 B 是 +33% 带宽且 ABI 破坏；新需求用新类型 `fv_vec4`。

---

## 6. 引用的既有证据与本笔记新增来源

**仓库内既有证据（直接依据）**

- [measured.md](measured.md) §0 条 1–8、§1.2（批量 API 形态与 x4 折中建议）、§1.5（ABI）、§3.1（三档主表）、
  §3.3（带宽）、§3.6（B1/B7 反汇编）、§4.6（档1 慢的原因）、§7（已验证 200 的链接）。

- [simd-hardware.md](simd-hardware.md) §1 条 1–8 与 §2。**本笔记不重复其 SIMD 细节，只交叉引用。**

- [algorithm-tricks.md](algorithm-tricks.md) §1 条 1 与 §4 的「SoA ＝ 把维度当循环」（约 15x）、§5 条 2（Box2D 休眠体加索引/掩码约束）、
  §5 末段、§6 条 1、§6 条 4、§6 条 6。

- [raw/env.txt](raw/env.txt)、[raw/membw.txt](raw/membw.txt)、[raw/evidence.txt](raw/evidence.txt)；
  `bench/membw.c`、`bench/bench.cpp`、`Makefile` 的 `bench*` 目标。三者**只读，未修改**。

**新增外部来源**（**均未在本机验证可达性**，除注明者）

1. Agner Fog《Optimizing software in C++》<https://www.agner.org/optimize/optimizing_cpp.pdf>，measured.md §7 已验证 200。
2. System V AMD64 ABI <https://gitlab.com/x86-psABIs/x86-64-ABI>，measured.md §7 已验证 200。

3. simdjson <https://arxiv.org/abs/1902.08318> · <https://github.com/simdjson/simdjson>
4. EnTT <https://github.com/skypjack/entt>
5. Bullet 顺序冲量求解器
   <https://github.com/bulletphysics/bullet3/blob/master/src/BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.cpp>
6. PhysX <https://github.com/NVIDIA-Omniverse/PhysX>
7. DirectXMath <https://github.com/microsoft/DirectXMath/blob/main/Inc/DirectXMath.h#L402>（simd-hardware.md 已引用）
   · <https://learn.microsoft.com/en-us/windows/win32/dxmath/pg-xnamath-programming>

8. Mike Acton《Data-Oriented Design and C++》(CppCon 2014) <https://www.youtube.com/watch?v=rX0ItVEVjHc>
9. Richard Fabian《Data-Oriented Design》<https://www.dataorienteddesign.com/dodbook/>
10. `aligned_alloc` <https://en.cppreference.com/w/c/memory/aligned_alloc>
