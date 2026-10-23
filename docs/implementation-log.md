# flashvec P0 / P1 落地记录与验收

> 对应 `docs/optimization-report.md` 的 §5 改造路线。
> 这一段记录「实际做了什么、证据是什么、还差什么」，不重复报告里的分析。

## 1. 交付范围

| 阶段 | 内容 | 状态 |
|---|---|---|
| P0 | Makefile 加 `-O3 -fno-math-errno -ffp-contract=off` | 完成 |
| P0 | `fv_vec2/vec3` 的 `min / max / clamp / madd` | 完成 |
| P0 | `fv_aabb2` + `fv_aabb2_merge / overlaps` | 完成 |
| P0 | 修 `FV_EPS` 量纲 bug（新增 `FV_EPS_SQR`，`project` 改用平方阈值） | 完成 |
| P1 | AoS 4 路交错：`add / sub / scale / lerp / madd` 的 `_x4` | 完成 |
| P1 | 稀疏变体：`add_x4_mask`、`add_x4_idx`、`scale_x4_idx` | 完成 |
| P1 | SoA：`add / scale / dot / len / normalize` 的 `_soa` | 完成 |
| P1 | AVX2 内核 + 运行时派发（非 AVX2 机器走标量回退） | 完成 |
| — | 文档：头文件契约注释、`README.md`、本记录 | 完成 |

文件：`src/flashvec_batch.cpp`（新，批量入口与 AVX2 内核）、`src/flashvec.cpp`（P0 函数 + project 修正）、
`include/flashvec.h`（新声明 + 契约注释）、`tests/test_p0.cpp`、`tests/test_batch.cpp`、`tests/test_batch_c.c`（新）、
`Makefile`（编译选项 + 新目标，既有 `all/test/clean` 行为未变）、`verify/`（第三方独立验收）。

## 2. 验收证据（四路独立，互相印证）

### 2.1 正确性

| 来源 | 内容 | 结果 |
|---|---|---|
| 实现者测试 | `test_flashvec` / `test_p0` / `test_batch` / `test_batch_c` | 30 / 57 / 1098 条断言 + C 批量接口，全 0 失败 |
| 第三方验收（`verify/parity.cpp`） | 12 个逐元素批量函数，对 `-O2`、`-O3 -march=native`、ASan+UBSan 三档；n = 0..17 全尾长；mask 四种模式；idx 顺序/逆序/重复/乱序；越界写 canary | 12/12 PASS，三档 exit 0，ASan/UBSan 无报告 |
| 调度者独立复算（不复用上述代码） | 423 条检查：x4 五件套、mask 语义、逆序 idx、SoA add/len/normalize、project 阈值 | 423/423 通过 |
| 负对照 | 故意把 lerp 内核改成 FMA / 把 mask blend 选反 | 分别 15 处 / 26 处 FAIL —— 证明测试确实能咬人，且 AVX2 内核真的在跑 |

### 2.2 语义与安全

- **逐位一致**：`add/sub/scale/lerp/madd` 的 `_x4`/`_soa`/`_mask`/`_idx`、`len_soa`、`normalize_soa` 与标量 API 位模式完全相同。
  实现手段：这些路径只用 `vaddps/vsubps/vmulps/vdivps/vsqrtps`，**不用 FMA、不用 rsqrt/rcp 近似**，归约顺序也一致；
  库内 `objdump` 普查：`vfmadd` = 0、`rcpps/rsqrtps` = 0、`haddps` = 0。
- **唯一例外 `dot_soa`**：4 累加器 + 末端一次归约，结果依赖归约顺序。实测：n ≤ 1000 相对误差 ≤ 4e-8；
  n = 8192 单次点积绝对误差 ≤ 4.8e-7、相对误差 ≤ 6.2e-8；接近抵消时可达 ~3e-6。头文件已按此改写（不再把语料相关值当上界）。
- **掩码语义**：`mask[i] == 0` 时显式写回 `a[i]` 原值（不是「不写」）。验收含 canary 检查：`_idx` 未引用的位置被改写次数 = 0。
- **纯 C 可用**：`tests/test_batch_c.c` 用 `cc -std=c11 -pedantic-errors` 编译并链接运行通过；头文件里没有任何 intrinsics。
- **既有行为未被改动**：旧 30 条断言一字未改；四个构建（旧选项/新源码、新选项、`-march=native`）下既有标量 API 输出位模式哈希相同，
  唯一变化是 `project` 的阈值修正（|onto| = 1e-4 由错误的「判零」变成正确投影，1e-7 仍判零）。

### 2.3 性能（本机 AMD Ryzen 7 5700U，AVX2/FMA，无 AVX-512）

N = 8192（每数组 32 KB，L2 驻留），相对**库调用式标量循环**的加速比（第三方 `verify/perf.cpp`，median of 5，校验和位级相等）：

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

N = 1<<20（工作集 24 MB > 8 MB L3）：`add` 13.977 → `_x4` 9.722（1.44x）/ `_soa` 9.455（1.48x），
`scale` 反而 0.92–0.93x —— **回到带宽墙，与报告 §2 的预测一致**。

## 3. 必须与数字一起读的三条限定

1. **收益主要来自「跨 TU 调用边界消失」，不是 intrinsics。** 同一份数学放进头文件内联后，`-O3 -march=native` 的自动向量化
   几乎追平手写 AVX2（add：`scalar_flat` 2.15x vs `_x4` 2.10x）。但**链接预编译静态库的调用方无法内联**——
   对他们 `_x4` 的 2.10x 是实打实的；而 `normalize`（3.35x）、`dot`（20.33x）、`collision`（5.08x）是内联标量循环也拿不到的
   （对应标量档只有 1.13x / 3.02x / 1.36x）。
2. **规模决定一切**：N = 1<<20 时三档收敛到 0.9–1.5x。报加速比必须同时报工作集大小。
3. **绝对值跨夹具不可比**：这台笔记本变频（1.64–3.97 GHz）且为共享机，同一份代码跨运行可差 20–25%。
   三套夹具（实现者 / 第三方 / 调度者）给出的绝对值互不相同（例如 `add_x4` 在 0.35–1.24 ns/elem 之间），
   只有**同一次运行内的比值**有意义。

## 4. 明确没做的（以及为什么）

- **`_idx` 变体没有 AVX2 内核**：是标量循环，实测只有 1.08–1.35x。原因是 Zen 2 的 `vgatherdps` 约 10 cycle/条（报告 §3.4）。
  **因此不要把 broad-phase 的收益按 `_mask` 的 1.7–3.8x 来估，要按 1.1–1.35x 估。**
- **尾巴走标量**（`_x4` 0–3 个、`_soa` 0–7 个），没做掩码加载：尾长没有性能意义。
- **三维批量入口未做**：`fv_vec3` 只有 12 字节，4 个装不进一个 256-bit 寄存器；要做得先决定是否引入 16 字节的 `fv_vec4`。
- **`_fast` 变体（rsqrt + 牛顿迭代）未做**：它改变最后一位，必须作为新函数提供；且实测只快约 10%。
- **运行时派发只在有 AVX2 的机器上验证过**：强制走标量回退的等价测试通过（把 11 处派发改成 `if (0)` 后 30+57+1098 全绿），
  但**真实非 AVX2 硬件上未运行过**。这是当前唯一缺少硬件级证据的路径。

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
| `-march=native` | 未引入；默认构建保持 x86-64 基线 + 函数级 `target` 属性派发（决策 3 的推荐项） |
| `-fno-math-errno` 写进默认 | 已写入（决策 4 的推荐项） |
| 是否需要 `fv_vec4` / 三维批量 | 未决——取决于物理引擎是否真的会用 3D |
| 是否需要 `_fast` 变体 | 未决——等有真实性能数据 |
