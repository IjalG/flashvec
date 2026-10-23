# flashvec 批量 SIMD 改动 —— 独立验证与基准（验收者 F）

独立第三方验收：不读实现者的自述，全部数字由本目录的夹具现场跑出。
**只拥有 verify/**，未修改 include/ src/ tests/ bench/ 根 Makefile docs/。

> **有效性边界（2026-10 重构后补记，请先读这段）**
> 库在本文写作之后经历了三次大改：批量入口的实现搬进 `include/flashvec.h` 的 `static inline`、
> 新增 `_stream` / `_mask_inplace` / `fv_vec4` / vec3·vec4 的 SoA 批量、以及 L3 尺寸的更正。
> 因此：
> * **§6 的结构检查已经按新架构重写**（`checks.sh`，11 条规则，逐条写明旧定义→新定义），
>   本文末尾「结构检查覆盖…」那一节是它的规则清单。
> * **§1 的位级对拍、§2 的性能数字来自重构前的快照**（`build/parity_*`、`build/perf_*` 的时间戳
>   是 10-17 / 10-23，早于这次重构）。它们证明的是当时那份源码；要在当前工作树上重新取得证据，
>   必须重跑 `make parity` / `make perf`。**在那之前，这两个数字表不可当作当前代码的验收证据。**
> * 机器口径也已更正：L3 是 **4 MiB/CCX × 2 CCX = 8 MiB 芯片总量，单线程可见 4 MiB**（见下），
>   本文 §2 里凡是按「8 MB L3」推的拐点都要按 4 MiB 重算。

* 机器：AMD Ryzen 7 5700U（8 核 16 线程，AVX2 + FMA，L2 512 KiB/核，L3 4 MiB/CCX x2）
* 编译器：g++ (Ubuntu 12.4.0-2ubuntu1~24.04.1) 12.4.0，内核 6.17.0-19-generic
* 被验证版本（重构前）：src/flashvec.cpp（P0）+ src/flashvec_batch.cpp（P1）+ include/flashvec.h。
  重构后：实现的唯一一份在 `include/flashvec.h`（`static inline`），`src/flashvec_batch.cpp`
  只剩 33 个 `extern` 转发符号，另新增 `src/flashvec_vec4.cpp`。
* 本目录产物：parity.cpp（位级一致性）、perf.cpp（基准）、c_api_check.c（纯 C ABI）、checks.sh（结构检查）、
  ref_strict.{h,cpp}（我自己的标量参考，单独 TU、-ffp-contract=off）、Makefile（独立构建）、build/*.txt（原始输出留档）

---

## 1. 位级一致性（parity）：全部 PASS

比较方式：每个逐元素批量函数都与 **(a) 库自己的标量 API** 与 **(b) 我独立写的标量参考**逐元素比较 IEEE 位模式。
覆盖 **全部尾长 n = 0..17**、多组随机/特殊语料（含 0、-0.0、次正规 1e-45、1e-30、FV_EPS 邻域
9.9999994e-7 / 1e-6 / 1.0000001e-6、±1e30、±3.4e38、FLT_MIN）、mask 的四种模式、idx 的顺序/逆序/
大量重复/随机乱序，并对每个函数的输出缓冲区做越界写 canary（o[n..) 不允许被写）。

    [逐元素批量路径 vs 标量路径 位级一致性]   （-O2 默认档与 -O3 -march=native 档结果相同）
      函数                     判定   samples  vs库标量 vs独立参考 备注
      fv_vec2_add_x4             PASS     1224     PASS     PASS     NaN对=2
      fv_vec2_sub_x4             PASS     1224     PASS     PASS
      fv_vec2_scale_x4           PASS     11016    PASS     PASS
      fv_vec2_lerp_x4            PASS     11016    PASS     PASS
      fv_vec2_madd_x4            PASS     11016    PASS     PASS
      fv_vec2_add_x4_mask        PASS     1530     PASS     PASS     NaN对=8
      fv_vec2_add_x4_idx         PASS     914      PASS     PASS
      fv_vec2_scale_x4_idx       PASS     914      PASS     PASS
      fv_vec2_add_soa            PASS     1224     PASS     PASS
      fv_vec2_scale_soa          PASS     11016    PASS     PASS
      fv_vec2_len_soa            PASS     612      PASS     PASS
      fv_vec2_normalize_soa      PASS     1224     PASS     PASS

    [P0 检查]
      FV_EPS / FV_EPS_SQR 取值                 PASS  FV_EPS=9.99999997e-07 FV_EPS_SQR=9.99999996e-13 位=FV_EPS_SQR:0x2b8ccc
      fv_vec2 min/max/madd（100 组）           PASS
      fv_vec2_clamp 夹取（75 组）              PASS
      fv_vec3 min/max/clamp/madd（1000 组）    PASS
      fv_aabb2_merge / overlaps                PASS  A∩B=1 A∩A=1 A∩C=0 A∩D=0
      比较器自检（NaN payload / ±0 / 1ulp）    PASS
      fv_aabb2_overlaps 边界相触               INFO  角点相触=0（规格未定义，仅记录）
      min/max 的 ±0 行为                       INFO  min((-0,0),(0,-0))=(0x00000000,0x80000000)

三档全绿、退出码 0：**-O2 无 -march=native**（默认构建档）、**-O3 -march=native**、**ASan+UBSan**（-O1）。
原始输出：verify/build/parity_O2_plain.txt、parity_O3_native.txt、parity_asan.txt。

**NaN 处理方式（比较策略，夹具里已自检）**：默认按 IEEE 位模式 memcmp；唯一例外是**两侧都是 NaN**
（payload 可不同）判为相等并单独计数（上表 NaN对=2/8）。±0.0 的差异按位如实判为不一致并单独分类
SIGNED_ZERO；本次没有出现。含 NaN/Inf 的输入单独统计（edge），未计入主判定 —— 本次 edge 分歧为 0。
canary：**无任何越界写**；_idx 变体对未被 idx 引用的 o[i] 也**没有改写**（UNTOUCHED=0，与头文件 line 167 一致）。
ASan+UBSan 无报错（我先修掉了自己夹具里的 new Stat() 泄漏，之后 exit=0）。

**dot_soa（规格允许重排归约顺序 -> 不做位级判定，只统计误差）**

    fv_vec2_dot_soa            PASS*    144 cases   逐位相等 143/144
                              max|abs err| = 4.77e-07，max rel err = 6.22e-08（语料含 1e30 量级）
    dot_soa 头文件断言(n<8 纯顺序累加)  PASS  n<8 逐位一致 64/64（头文件 line 149 的断言成立）
    逐 n 逐位相等次数: 0..7 各 8/8；8:7/8；9..17 各 8/8

---

## 2. 基准（perf）：ns/元素、加速比、校验和

计时：clock_gettime(CLOCK_MONOTONIC)，每项 5 次取样取中位数，取样前先 warm-up 并校准内层重复次数
（目标约 6 ms/取样）。**校验和 = 单次调用的输出累加**（与内层重复次数无关）；同一次运行内，除 dot_soa
（后者报相对误差）外，各档校验和必须位级相等，否则该轮判无效 —— 本次全部位级相等。

口径提醒：这台笔记本变频（本次起始平均 **3125.2 MHz**）且与实现者的构建/基准并行争用 CPU，
**跨夹具的绝对值不可比**（实现者测到 add_x4 = 0.346 ns/元素，我测到 1.240）。可比的是同一次运行内的比值。

原始输出：verify/build/perf_O3_native.txt（修订后夹具，-O3 -march=native）—— 逐行照抄：

    起始平均 CPU 频率: 3125.2 MHz

    === N=8192 (L2 驻留, 工作集 AoS a/b/o=0.2 MB) ===
      op           impl              ns/elem    speedup           checksum cycle/el
      add          scalar_api          2.601      1.00x          -206.8001     8.13
      add          scalar_flat         1.207      2.15x          -206.8001     3.77
      add          x4                  1.240      2.10x          -206.8001     3.87
      add          soa                 1.369      1.90x          -206.8001     4.28
      scale        scalar_api          2.776      1.00x          7.7606516     8.68
      scale        scalar_flat         1.154      2.41x          7.7606516     3.61
      scale        x4                  1.079      2.57x          7.7606516     3.37
      scale        soa                 1.027      2.70x          7.7606516     3.21
      lerp         scalar_api          3.002      1.00x          -75.61176     9.38
      lerp         scalar_flat         1.633      1.84x          -75.61176     5.10
      lerp         x4                  1.196      2.51x          -75.61176     3.74
      lerp         soa(无此入口)          -          -                  -
      dot          scalar_api          3.087      1.00x          19.535582     9.65
      dot          scalar_flat         1.024      3.02x          19.535582     3.20
      dot          x4(无此入口)          -          -                  -
      dot          soa                 0.152     20.33x          19.535637     0.47
        [校验和] X dot/soa cs=19.535636901855469 ref=19.535581588745117 rel=2.83e-06 -> 触到我设的 1e-6 容差线
      len          scalar_api          3.046      1.00x          12501.319     9.52
      len          scalar_flat         1.649      1.85x          12501.319     5.15
      len          x4(无此入口)          -          -                  -
      len          soa                 1.149      2.65x          12501.319     3.59
      normalize    scalar_api          4.382      1.00x         -54.362158    13.69
      normalize    scalar_flat         3.882      1.13x         -54.362158    12.13
      normalize    x4(无此入口)          -          -                  -
      normalize    soa                 1.306      3.35x         -54.362158     4.08
      collision    scalar_api          8.481      1.00x            6.44906    26.51
      collision    scalar_flat         6.216      1.36x            6.44906    19.43
      collision    soa_hybrid          1.670      5.08x            6.44906     5.22
      collision    x4(无此入口)          -          -                  -
      add_mask     scalar(5%)          2.083      1.00x         -13.461817     6.51
      add_mask     x4_mask(5%)         1.227      1.70x         -13.461817     3.83
      add_mask     scalar(50%)         4.648      1.00x         -200.68404    14.53
      add_mask     x4_mask(50%)        1.223      3.80x         -200.68404     3.82
      add_mask     scalar(95%)         3.971      1.00x         -138.84498    12.41
      add_mask     x4_mask(95%)        1.210      3.28x         -138.84498     3.78
      add_idx      scalar(5%)          1.021      1.00x          15.822212     3.19
      add_idx      x4_idx(5%)          0.950      1.08x          15.822212     2.97
      add_idx      scalar(50%)         2.248      1.00x         -56.604818     7.02
      add_idx      x4_idx(50%)         1.670      1.35x         -56.604818     5.22

    === N=1048576 (撞带宽墙, 工作集 AoS a/b/o=24.0 MB) ===
      op           impl              ns/elem    speedup           checksum cycle/el
      add          scalar_api         13.977      1.00x          603.99968    43.68
      add          scalar_flat        13.904      1.01x          603.99968    43.45
      add          x4                  9.722      1.44x          603.99968    30.38
      add          soa                 9.455      1.48x          603.99968    29.55
      scale        scalar_api          7.772      1.00x           1274.378    24.29
      scale        scalar_flat         8.285      0.94x           1274.378    25.89
      scale        x4                  8.354      0.93x           1274.378    26.11
      scale        soa                 8.447      0.92x           1274.378    26.40
      lerp         scalar_api         12.920      1.00x          544.78358    40.38
      lerp         scalar_flat        10.500      1.23x          544.78358    32.81
      lerp         x4                 12.319      1.05x          544.78358    38.50
      lerp         soa(无此入口)          -          -                  -
      dot          scalar_api          5.268      1.00x          3677.5696    16.46

**N=1<<20 这一档的表在这里被截断**：这一档的后半段（dot/soa 之后的 len、normalize、collision、稀疏）
没有取得可信数据 —— 原因是我发现一个**早先版本基准进程**在 26 分钟后仍未结束、占着一整颗核（已 kill），
它同时污染了这一档的数字；我按调度者要求停止重测，**不补估**缺行。N=8192 一档完整、可用。

读数要点：
* 撞带宽墙（N=1<<20）后 add/scale/lerp 全部回落到约 0.9-1.5x，与「三档在带宽墙前收敛」的预期一致。
* soa 档在 dot 上 20.33x、normalize 3.35x、len 2.65x —— 归约型/带分支的运算是 SoA 的主场。
* _mask 档在 50% 稀疏度 3.80x；_idx 档只有 1.08-1.35x（见风险点 2）。

---

## 3. 结论：收益主要来自「跨 TU 调用边界消失」，不是 intrinsics

* add @N=8192：scalar_api 2.601 -> scalar_flat（同 TU 内联、-O3 自动向量化）**1.207** -> x4 **1.240**。
  内联标量循环与库的 AVX2 打包路径**齐平**（2.15x vs 2.10x），intrinsics 本身没有额外优势。
* 真正的差距在 scalar_api -> scalar_flat 这一段：2.15x 几乎全部是「每元素一次跨 TU 调用」的代价。
  头文件 line 128 声称边界成本约 4.2 cycle/元素、是打包收益的 3 倍 —— 我的独立测量支持这个判断
  （8.13 -> 3.77 cycle/元素）。
* **限定条件**：scalar_flat 只是「调用方自己拥有这段循环、且与本 TU 一起被 -O3 -march=native 编译」时的上限。
  **链接预编译静态库的调用方做不到**（GCC 拒绝跨 TU 向量化：missed: statement clobbers memory），
  他们只能拿到 scalar_api 的 2.601；对他们而言 x4 的 2.10x 是实打实的。这正是批量入口存在的理由。
* 批量入口另有内联标量循环也拿不到的收益：normalize 3.35x（标量带 per-element 分支，GCC 只做到 1.13x）、
  dot 20.33x（4 累加器 + 末端归约，标量只能 3.02x）、collision 微内核 5.08x（标量 1.36x）。

---

## 4. 复现命令

    cd /home/hjl/桌面/codes/QWCIF/flashvec/verify

    make parity   # 位级一致性：-O2 / -O3 -march=native / ASan+UBSan 三档，原始输出 build/parity_*.txt
    make perf     # 基准三档（-O3 -march=native / 发布默认 -O3 / 旧 -O2），原始输出 build/perf_*.txt
                  #   用 taskset -c 2 + stdbuf -o0 运行，避免迁移与「崩溃丢失缓冲」
    make check    # 纯 C11 头文件+ABI 检查（-pedantic-errors）、既有 30 条断言（我重新编译）、checks.sh 结构检查
    make extra    # 实现者 E 自己的测试套件（我只编译运行，不修改）
    make clean

    # 单跑某一份
    ./build/parity_O3_native          # 或 parity_O2_plain / parity_asan
    taskset -c 2 stdbuf -o0 ./build/perf_O3_native

关键原始输出（都在 verify/build/）：parity_O2_plain.txt、parity_O3_native.txt、parity_asan.txt、
perf_O3_native.txt、c_api_check.txt、test_flashvec_F.txt、test_c_api_F.txt、checks.txt。

### 结构检查（checks.sh，只读仓库文件）：11 条规则，逐条都有 PASS/FAIL，任何一条 FAIL → 退出码非 0

重构后重写过。每条规则在脚本里都写着「旧定义 → 新定义 + 为什么必须改」。要点：

| # | 现在检查什么 | 相对旧规则变了什么 |
|---|---|---|
| 1 | 纯 C（`cc -std=c11 -pedantic-errors`）可编译；`immintrin` 展开行数 **> 0**（正对照）；头文件无 C++ 语法 | **判据反转**：旧规则要求头文件里 `immintrin` 行数为 0，重构后实现就在头文件里，为 0 反而说明封装坏了 |
| 2 | `FV_USE_EXTERN` 形态的**包含集**里没有 `math.h` / `string.h` / `immintrin.h` / `stdlib.h`，且可被纯 C `-pedantic-errors` 编译 | **判据换形式**：旧规则是「头文件不得含 target 属性/builtin（应为纯声明）」，现在把等价事实变成可判定的包含集比较 |
| 3 | `make -C .. -n -B all` 打出的每条 `-c src/…` 命令都不含 `-march=native/-mavx/-ffast-math/-mfma`；`FMA_FLAGS` 只被 `*_cpp_fma.o` 引用 | **判据落到实际命令**：旧规则 grep 全文件，会把「敌意选项对拍夹具」判死；现在直接问 make 默认目标执行什么（更严） |
| 4 | `__builtin_cpu_supports` 与 `target(` 扫**头文件**；32 处派发全部写成 `!= 0`；并报告 `src/*.cpp` 里为 0 | **扫描目标迁移**：旧规则扫 `src/*.cpp`，重构后那里一处都没有 |
| 5 | `include/` + `src/` 里 FMA / rcpps / rsqrtps / **maskload / maskstore（VMASKMOVPS）** 内建名 = 0；注释里提到助记符不算 | **范围扩到实现真正所在处**，并新增 VMASKMOVPS 禁用项 |
| 6 | 头文件写明归约顺序免责 + 掩码语义（原样保留，新增 PASS 判据） | 不变 |
| 7 | 机器/编译器事实（INFO，不设判据）；额外打印 `index2/index3` 的 `size` 与 `shared_cpu_list` | 新增缓存层级实测，防止文档与机器再次漂移 |
| 8 | **产物级指令普查**：脚本自己按 `-O3 -march=native -mfma -ffp-contract=fast` 编 4 个对象（3 个库 TU + 纯 C 的 `test_stream.c`），`vfmadd/vfmsub/vfnmadd/vrcpps/vrsqrtps/haddps/vmaskmovps` **合计必须为 0**，`vmovntps` **必须 > 0** | **对象不再用 build/ 里的过期产物**；配置从 `-O3 -march=native` 升级到最狠的一档（并加 `-mfma`），禁用项从 6 项扩到 7 项，新增 NT store 正对照 |
| 9 | 不带 `-march=native` 的默认档里：`fv_i_*_avx2` 局部符号存在、有 AVX2 指令、对象引用 `__cpu_model` 且在 `0x400`（avx2 位）上判断 | **判据改成消费者视角**：内核是头文件 static inline，所以「内核还在」对每个消费 TU 都成立；派发证据从「call 关系」换成「`__cpu_model` + avx2 位」 |
| 10 | 默认档实测摘要（只引已有产物；本脚本**不做性能测量**），并识别 26 字节的已知垃圾文件 | 不变（新增垃圾文件识别） |
| 11 | **新增**：`fv_i_maskload_add4_avx2` / `fv_i_maskstore_add4_avx2` 在整仓出现次数必须为 0 | 全新规则：这两个名字绑定被实测否决的 VMASKMOVPS 稀疏策略，不许复活 |

运行：`bash verify/checks.sh`（同时写 `build/checks.txt`；退出码 0 = 全部 PASS）。

---

## 5. 我最担心的 3 个风险点

1. **非 AVX2 机器的回退路径我没能真正跑过。** 本机有 AVX2，我无法强制 __builtin_cpu_supports("avx2")==0，
   所以整段标量回退逻辑只被「尾长 n%4」间接覆盖（0..17 全覆盖）。我能给的只是结构性证据（派发写法、
   内核是文件内静态符号、默认档对象里内核存在）。若 E 的测试也没有真正覆盖回退分支，那么这条路径是
   全库唯一没有运行时证据的路径 —— 建议用 qemu-x86_64 或 CI 上的非 AVX2 机器补一次。
2. **_idx 变体实际上没有 AVX2 内核（等价于标量循环），收益只有 1.08-1.35x。** 指令普查里没有
   fv_*_idx_avx2 符号，fv_vec2_add_x4_idx 只有 1 条 vaddps。头文件只写了「代价是随机访存」，
   没有说它在 SIMD 意义上并不比调用方自己的标量 gather 循环快。把 broad-phase 建在 _idx 上的收益预期
   要按约 1.1-1.35x 估，而不是按 _mask 的 1.7-3.8x 估；另外 _idx 语义（不写未引用位置）虽然我验证通过，
   但它是「保留 o 旧值」的唯一入口，调用方极易与 _mask 的「写回 a[i]」混淆。
3. **dot_soa 没有真正的误差上界，头文件里的数字是语料相关的。** 头文件 line 146 写「误差量级约 1e-6 相对」，
   并引用 tests/test_batch.cpp 的「<= 4e-8」；而我用小 n（0..17）、含 1e30 量级与阈值附近的对抗语料实测到
   **max rel err = 6.22e-08**，在 perf 的 ±2 随机语料上 checksum 层面看到 **rel = 2.83e-06**（超出 4e-8，
   也触到我夹具的 1e-6 容差线）。相对误差在有抵消的语料上会被放大，所以「<=4e-8」不能当上界用 ——
   建议 docs 给出「绝对误差 + 语料条件」而不是裸的相对误差，或明确写「与语料相关，实测 4e-8 ~ 3e-6」。

（另有两个未定义的边界行为，只作记录、不判失败：fv_aabb2_overlaps 角点相触返回 0（exclusive）；
min/max 是 (a<b)?a:b 而非 IEEE fminf，min((NaN,1),(0,0))=(0,0)，头文件已写明该契约。）

---

## 6. 有没有发现实现者 E 未提及的问题？

**有 3 条**（都不是「实现写错了」，而是「头文件/回报里没有说清，或说了但我测出反例」）：

1. dot_soa 的误差数字：头文件/回报引用 4e-8，我的对抗语料实测 6.22e-08，perf 语料上 2.83e-06（见风险点 3）。
2. _idx 变体**没有 AVX2 内核**，实测只有 1.08-1.35x（见风险点 2）。E 只提到随机访存代价。
3. 不带 -march=native 的默认构建里 AVX2 内核**确实存在**、派发也确实是 != 0（我用 nm 找到 11 个
   fv_*_avx2 静态符号并在 objdump 里统计了指令）—— 这条是**支持** E 的，但 E 的回报里没有给这个经验证据，
   我和调度者一开始都误判成「默认档是纯标量」，值得记一笔。

位级一致性、掩码语义、idx 不越界、P0 小工具、纯 C 头文件、既有 30 条断言、FMA/倒数近似的禁用 —— 这 7 项
我独立复验**全部通过**，与 E 的回报一致。**没有任何一条位级不一致或内存语义问题被我测出来。**

---

## 7. 已知夹具局限（如实列出）

* 绝对值受频率调节与并行负载影响大（本次起始平均 3125.2 MHz）；只有同一次运行内的比值可比。
  跨夹具对比实现者的 0.346 ns/元素与我的 1.240 无意义。
* N=1<<20 的后半段缺失（见第 2 节）；我不补估、不外推。
* verify/build/perf_O2_plain.txt 与 perf_O3_plain.txt（各 26 字节，内容为 /bin/sh: 3: IN: not found）是 23:08
  那次 Makefile 单字符变量转义 bug（$P 被 make 吃掉，只留 IN）留下的过期垃圾，**不是有效测量**；
  这两个文件会在下一次 make perf 时被正常覆盖。（该 bug 与库本身无关。）
* x4 覆盖 add/sub/scale/lerp/madd，soa 覆盖 add/scale/dot/len/normalize —— 冻结接口本身没有
  dot_x4 / len_x4 / normalize_x4 / lerp_soa / sub_soa，因此这些格子是「无此入口」而不是我没测。
* collision 微内核的 soa_hybrid 只有归一化阶段向量化：反射需要「逐元素点积 + 逐元素缩放」，
  冻结接口给不出来（scale_soa 只接受单个标量 s），反射阶段只能标量循环。这条对 P2 的布局选型有用。
