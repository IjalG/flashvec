# flashvec P0 / P1 落地记录与验收

> [English](../implementation-log.md) | 简体中文

> 对应 `docs/optimization-report.md` 的 §5 改造路线。
> 这一段记录「实际做了什么、证据是什么、还差什么」，不重复报告里的分析。

## 1. 交付范围

| 阶段 | 内容 | 状态 |
|---|---|---|
| P0 | Makefile 加 `-O3 -fno-math-errno -ffp-contract=off` | 完成 |
| P0 | `fv_vec2/vec3` 的 `min / max / clamp / madd` | 完成 |
| P0 | `fv_aabb2` + `fv_aabb2_merge / overlaps` | 完成 |
| P0 | 修 `FV_EPS` 量纲 bug：新增 `FV_EPS_SQR`，`project` 改用平方阈值 | 完成 |
| P1 | AoS 4 路交错：`add / sub / scale / lerp / madd` 的 `_x4` | 完成 |
| P1 | 稀疏变体：`add_x4_mask`、`add_x4_idx`、`scale_x4_idx` | 完成 |
| P1 | SoA：`add / scale / dot / len / normalize` 的 `_soa` | 完成 |
| P1 | AVX2 内核 + 运行时派发；非 AVX2 机器走标量回退 | 完成 |
| 重构 | **批量入口搬进头文件 `static inline`**（默认形态）；库内保留同名 `extern` 符号供 `FV_USE_EXTERN` 与旧消费方 | 完成 |
| P1+ | NT store（`_stream`）6 个入口 + `_mask_inplace` 1 个 | 完成 |
| P1+ | `fv_vec4` 类型（标量 18 个）+ `_x2` AoS 批量 5 个 + `_soa` 归约批量 3 个 | 完成 |
| P1+ | `fv_vec3` 的 SoA 批量 5 个（三维不加 padding、只走 SoA） | 完成 |
| — | 文档：头文件契约注释、`README.md`、本记录、`docs/acceptance-protocol.md` | 完成 |

最后一行的 `docs/acceptance-protocol.md` 指验收判据那份。

改动落在这些文件上：

| 文件 | 说明 |
|---|---|
| `src/flashvec_batch.cpp` | 新；批量入口与 AVX2 内核 |
| `src/flashvec.cpp` | P0 函数 + `project` 修正 |
| `include/flashvec.h` | 新声明 + 契约注释 |
| `tests/test_p0.cpp` | |
| `tests/test_batch.cpp` | |
| `tests/test_batch_c.c` | 新 |
| `Makefile` | 编译选项 + 新目标，既有 `all/test/clean` 行为未变 |
| `verify/` | 第三方独立验收 |

## 2. 验收证据（四路独立，互相印证）

### 2.1 正确性

| 来源 | 结果 |
|---|---|
| 实现者测试 | 全 0 失败 |
| 第三方验收（`verify/parity.cpp`） | 12/12 PASS，三档 exit 0，ASan/UBSan 无报告 |
| 调度者独立复算，不复用上述代码 | 423/423 通过 |
| 负对照 | 分别 15 处 / 26 处 FAIL |

四路的测试内容如下。

- 实现者测试跑 `test_flashvec` / `test_p0` / `test_batch` / `test_batch_c`，30 / 57 / 1098 条断言，外加 C 批量接口。
- 第三方验收：12 个逐元素批量函数，对 `-O2`、`-O3 -march=native`、ASan+UBSan 三档。
- 第三方验收的覆盖面：n = 0..17 全尾长；mask 四种模式；idx 顺序 / 逆序 / 重复 / 乱序；越界写 canary。
- 调度者独立复算：423 条检查，覆盖 x4 五件套、mask 语义、逆序 idx、SoA add/len/normalize、project 阈值。
- 负对照：故意把 lerp 内核改成 FMA、把 mask blend 选反，分别出现 15 处 / 26 处 FAIL。
- 负对照的意义：证明测试确实能捕获错误，且 AVX2 内核真的在跑。

### 2.2 语义与安全

- **逐位一致**：`add/sub/scale/lerp/madd` 的 `_x4`/`_soa`/`_mask`/`_idx` 变体、`len_soa`、`normalize_soa` 与标量 API 位模式完全相同。
- 实现手段：这些路径只用 `vaddps/vsubps/vmulps/vdivps/vsqrtps`，**不用 FMA、不用 rsqrt/rcp 近似**，归约顺序也一致。
- 库内 `objdump` 普查：`vfmadd` = 0、`rcpps/rsqrtps` = 0、`haddps` = 0。

- **唯一例外 `dot_soa`**：4 累加器 + 末端一次归约，结果依赖归约顺序。
- 实测：n ≤ 1000 相对误差 ≤ 4e-8。
- n = 8192 时单次点积绝对误差 ≤ 4.8e-7、相对误差 ≤ 6.2e-8；接近抵消时可达 ~3e-6。
- 头文件已按此改写，不再把语料相关值当上界。

- **掩码语义**：`mask[i] == 0` 时显式写回 `a[i]` 原值，不是「不写」。
- 验收含 canary 检查：`_idx` 未引用的位置被改写次数 = 0。

- **纯 C 可用**：`tests/test_batch_c.c` 与 `tests/test_header_c.c` 用 `cc -std=c11 -pedantic-errors` 编译并运行通过。
- **注意：重构后口径已变。** 批量入口的实现现在就在头文件里。
- 所以 x86 + GCC/Clang 的纯 C 编译单元**会**展开 `<immintrin.h>`，实测 `cc -std=c11 -E | grep -c immintrin` = 52。
- 「头文件里没有 intrinsics」这句话只在两种形态下仍成立：`FV_USE_EXTERN`（纯声明，包含集与重构前逐字相同）与 `FV_I_NO_INTRINSICS`（纯 C 回退）。
- `verify/checks.sh` 的 [1][2] 段就是在检查这两种形态。

- **既有行为未被改动**：旧 30 条断言一字未改。
- 四个构建（旧选项/新源码、新选项、`-march=native`）下，既有标量 API 输出的位模式哈希相同。
- 唯一变化是 `project` 的阈值修正：|onto| = 1e-4 由错误的「判零」变成正确投影，1e-7 仍判零。

### 2.3 性能（本机 AMD Ryzen 7 5700U，AVX2/FMA，无 AVX-512）

N = 8192，每数组 32 KB，L2 驻留。下表是相对**库调用式标量循环**的加速比。

数据来自第三方 `verify/perf.cpp`，取 5 次中位数，校验和位级相等。

| 运算 | 标量 API | `_x4` | `_soa` |
|---|---|---|---|
| add | 2.601 ns/elem | 1.240（2.10x） | 1.369（1.90x） |
| scale | 2.776 | 1.079（2.57x） | 1.027（2.70x） |
| lerp | 3.002 | 1.196（2.51x） | — |
| dot | 3.087 | — | 0.152（**20.33x**） |
| len | 3.046 | — | 1.149（2.65x） |
| normalize | 4.382 | — | 1.306（3.35x） |
| collision 微内核 | 8.481 | — | 1.670（5.08x） |
| add_mask（50% 命中） | 4.648 | 1.223（3.80x） | — |
| add_idx（50% 命中） | 2.248 | 1.670（1.35x） | — |

N = 1<<20 时工作集 24 MB，远超本机单线程可见的 4 MiB L3。

芯片总量 8 MiB，由两个 CCX 各 4 MiB 组成，拐点因此在 1.75×10⁵ 元素而不是 3×10⁵。

这一档 `add` 13.977 → `_x4` 9.722（1.44x）/ `_soa` 9.455（1.48x）。

`scale` 反而 0.92–0.93x，回到内存带宽上限，与报告 §2 的预测一致。

## 3. 必须与数字一起读的三条限定

1. **收益主要来自「跨 TU 调用边界消失」，不是 intrinsics。**
   - 同一份数学放进头文件内联后，`-O3 -march=native` 的自动向量化几乎追平手写 AVX2，add 上 `scalar_flat` 2.15x vs `_x4` 2.10x。
   - 但**链接预编译静态库的调用方无法内联**，对他们 `_x4` 的 2.10x 是实打实的。
   - 而 `normalize`（3.35x）、`dot`（20.33x）、`collision`（5.08x）是内联标量循环也拿不到的，对应标量档只有 1.13x / 3.02x / 1.36x。

2. **规模决定一切**：N = 1<<20 时三档收敛到 0.9–1.5x。报加速比必须同时报工作集大小。

3. **绝对值跨夹具不可比**：这台笔记本变频（1.64–3.97 GHz）且为共享机，同一份代码跨运行可差 20–25%。
   - 三套夹具（实现者 / 第三方 / 调度者）给出的绝对值互不相同，例如 `add_x4` 在 0.35–1.24 ns/elem 之间。
   - 只有**同一次运行内的比值**有意义。

## 4. 明确没做的（以及为什么）

- **`_idx` 变体没有 AVX2 内核**：是标量循环，实测只有 1.08–1.35x。原因是 Zen 2 的 `vgatherdps` 约 10 cycle/条（报告 §3.4）。

- **因此不要把 broad-phase 的收益按 `_mask` 的 1.7–3.8x 来估，要按 1.1–1.35x 估。**

- **尾巴走标量**（`_x4` 0–3 个、`_soa` 0–7 个），没做掩码加载：尾长没有性能意义。

- ~~三维批量入口未做~~ → **本轮已做，且走 SoA**：`fv_vec3_add_soa / scale_soa / dot_soa / len_soa / normalize_soa`。

- `fv_vec3` 是 12 B，装进 ymm 会空转 1/4 lane，所以既不按 AoS 打包、也**不给** `fv_vec3` 加 padding（破坏 ABI）。

- 三维改为「同一分量的 8 个元素装一条 ymm」。

- ~~四维向量未做~~ → **本轮已做**：16 B 的 `fv_vec4` 标量 18 个 + `_x2` AoS 批量 5 个 + `_soa` 归约批量 3 个。

- `_x2` 的一条 ymm 正好装 2 个元素；它的数字含义与 `_x4` 一致：一个 256-bit 寄存器装几个元素。

- **`_fast` 变体（rsqrt + 牛顿迭代）未做**：它改变最后一位，必须作为新函数提供；且实测只快约 10%。

- **运行时派发只在有 AVX2 的机器上验证过**：强制走标量回退的等价测试通过，把 11 处派发改成 `if (0)` 后 30+57+1098 全绿。

- 但**真实非 AVX2 硬件上未运行过**。这是当前唯一缺少硬件级证据的路径。

## 5. 复现

~~~bash
make clean && make test-all      # 30 + 57 + 1098 条断言 + 两个 C 兼容性测试
cd verify && make parity         # 第三方位级一致性（O2 / O3-native / ASan 三档）
cd verify && make perf           # 第三方性能对比（含校验和一致性断言）
cd verify && make check          # 结构检查：头文件纯净度、编译选项、派发、指令普查
cd ../bench && make run          # 三档实现 vs 标量基准（另一个夹具，看趋势）
~~~

## 6. 待决策（与报告 §6 的对应关系）

| 项 | 现状 |
|---|---|
| `FV_EPS` 量纲 | 已按「拆成 `FV_EPS` + `FV_EPS_SQR`」实施（报告 §6 决策 1 的推荐项） |
| 批量 API 是否引入 | 已引入 `_x4` + `_soa`（报告 §6 决策 2 的推荐项） |
| `-march=native` | 未引入（决策 3 的推荐项），默认构建保持 x86-64 基线 + 函数级 `target` 属性派发 |
| `-fno-math-errno` 写进默认 | 已写入（决策 4 的推荐项） |
| 是否需要 `fv_vec4` / 三维批量 | **已决定并实施**，见下 |
| 是否需要 `_fast` 变体 | 未决，等有真实性能数据（仍未做） |
| `_stream`（NT store） | 已实施 6 个入口，见下 |
| `_mask_inplace` | 已实施 1 个入口，见下 |

- `fv_vec4` 与三维批量：三维走 SoA（5 个入口），四维走 `_x2`（5 个）+ `_soa`（3 个）+ 标量 18 个。

- `_stream`：已实施 6 个入口，由调用方显式选择。

- `_stream` 的读数：L2 内实测 **0.13–0.51x**，即慢 2–8 倍，只在工作集远超 L3 时划算。

- `_stream` 的上限：L3 之外 PASS 点 1.29–1.58x，上限按入口形状 1.333 / 1.500（`docs/acceptance-protocol.md` §9）。

- `_mask_inplace`：已实施 1 个入口，就地稀疏加，掩码外的**值**不变，按掩码密度二选一内核。

- `_mask_inplace` 的判定：**全档 INCONCLUSIVE**，相对 `_mask` 无实测吞吐优势。

- 它的价值在形态：就地写不产生 RFO，省不掉流量。

## 7. 重构后的架构事实与验收判据（本轮新增）

1. **实现的唯一一份在头文件里。**
   - 批量入口现在是 `include/flashvec.h` 里的 `static inline`。
   - 它含 AVX2 内核、`target("avx2")` 属性，以及 `__builtin_cpu_supports("avx2") != 0` 派发。
   - 实测计数：头文件里 32 处派发点、7 处 `target(`；`src/*.cpp` 里**各 0 处**。
   - `src/flashvec_batch.cpp` 只剩 33 个 `extern` 转发符号，用于二进制兼容。
   - 为什么这么改：跨 TU 调用边界是一阶瓶颈（§3 第 1 条的实测），而链接预编译静态库的调用方无法把这层边界消掉。

2. **三种形态必须给同一个结果**：默认（头文件 inline）、`FV_USE_EXTERN`（库内实现）、`FV_I_NO_INTRINSICS` / 纯 C 回退。
   - `make test-inline` / `test-vec34-inline` 直接对拍主输出流。

3. **头文件现在条件包含 `<immintrin.h>`**，条件是 `FV_I_INTRINSICS && !FV_USE_EXTERN`。
   - 对可能被收缩成 FMA 的入口加了 `#pragma GCC optimize("fp-contract=off")`。
   - `target("avx2")` 挡不住调用方命令行带来的 `-mfma`：实测 `-march=native -mfma -ffp-contract=fast` 下照样出现 `vfmadd`。
   - `verify/checks.sh` 的 [8] 段就是在最严的一档下做产物级指令普查：7 项禁用助记符 = 0，`vmovntps` > 0 作为 NT store 正对照。

4. **验收判据已经独立成文**：见 `docs/acceptance-protocol.md`。
   - 一条优化被接受必须**同时**满足五条：
     - 位级对拍。
     - 目标规模同进程中位比 ≥ 1.15x。
     - 同向复现。
     - 增益超过该点 A/A 地板。
     - 其他规模不回退超过 5%。
   - 「单点比值不可信」这条有硬证据：304 个大 N 测量点里可分辨的只有约 30–35%。
   - 四次运行的可分辨点数分别为 90 / 92 / 93 / 107，分母恒为 304。这个数随运行漂移，引用请带轮次。
   - 写成区间不削弱这条证据，反而是它更硬的地方：同一台机器、同一份代码，「有多少点可分辨」本身就在 90→107 之间漂。
   - 那么「某个 `(op, n)` 的单点比值」当然更不可复现。

5. **归约型入口不具备分段不变性**，这是本轮新增的契约事实，必须一起读。
   - `fv_vec2_dot_soa` / `fv_vec3_dot_soa` / `fv_vec4_dot_soa` 把同一数组切成多段分别调用再相加，结果与整段一次调用**不同**。
   - 实测 n=1000 相对差 ≈3.9e-6、n=8192 ≈2.4e-5，12/12 组语料全部不一致。
   - 同一次实验里 `fv_vec2_add_x4` 的 12/12 组**逐位一致**，这是对照组，证明实验有效。
   - 这是归约树随切分改变的**设计必然**，不是缺陷。但它与「与标量路径的误差 ≤ 4e-8」是两件事，两句话不能混用。
   - 需要「结果与切分无关」的调用方请用标量 `fv_vec2_dot` 顺序累加，或固定切分方式。
   - 逐元素型批量入口（`add/sub/scale/lerp/madd/len/normalize` 及其 `_x4`/`_soa`/`_mask`/`_idx` 变体）具备分段不变性。

6. **缓存口径已更正。**
   - 本机 L3 是 **4 MiB/CCX × 2 CCX = 8 MiB 芯片总量，单线程可见 4 MiB**，读自 `/sys/.../index3`：`size=4096K shared_cpu_list=0-7`。
   - 旧文档写「L3 8 MiB 共享」并据此把拐点算成 `3×10⁵` 元素，是错的。
   - 真实拐点约 **1.75×10⁵** 元素（@24 B/元素），实测断崖落在 2^17（38–57 GB/s）与 2^18（5.6–8.8 GB/s）之间。
