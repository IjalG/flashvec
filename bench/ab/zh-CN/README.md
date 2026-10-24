# bench/ab —— 同进程 A/B 交错性能与吞吐夹具

> [English](../README.md) | 简体中文

这套夹具只回答一个问题：**候选实现比基线快多少**。

答案必须在这台机器上回答得**可证伪**。这台机器在变频，资源被共享，还被别的负载占去 30-50%。

一句话结论：**它测得出来的，它一定标出来；它测不出来的，它拒绝下结论。**

> **本夹具的分辨率下限是硬约束**：304 个大 N 点里只有 **85-107 个**可分辨，约占 28-35%。
>
> 这个统计量**本身随运行漂移**。六次不同轮次分别读到 **85 / 90 / 92 / 93 / 95 / 107**，分母恒为 304，跨度约 26%。
>
> 所以引用必须带轮次。写成区间不是含糊，而是它唯一正确的写法。
>
> 最后一次实测：182/182 门禁通过，全点中位比 1.0109，95/304 可分辨。
>
> 个别 (op, size) 上**完全相同的两份代码**能稳定差 15-20%。这就是 A/A 地板，见 §10。
>
> 因此「同运行内 median >= 1.15x」**单独不足以判真**。判据见 §11。

---

## 1. 为什么不能「跑两个二进制再比」

README「性能」第 3 条给了限定：同一份代码跨运行可差 20-25%，频率、邻居负载、页映射都变过。

所以比法只有一种：**同一个进程、同一份数据、同一段时间内，A/B/A/B… 交替跑**，再按档取中位数相除。

这就是 ab_bench 的全部设计出发点。

## 2. 同进程双实现是怎么做到的

| 角色 | 源码 | 编译方式 |
|---|---|---|
| baseline | snapshot/，即 git show HEAD: 抓到的快照 | -include 自动生成的 build/rename_a.h，入口变成 ab_a_* |
| candidate | 默认同一份快照，用于自检；或任意源码树 | -include build/rename_b.h，变成 ab_b_* |
| 标量锚点 | HEAD 快照的 src/flashvec.cpp | 不改名，循环在**另一个 TU** 里调用它 |

snapshot/ 由 git show HEAD: 抓到的 include/flashvec.h、src/flashvec_batch.cpp 与 src/flashvec.cpp 组成。

candidate 想换源码树，用 CAND_HDR / CAND_SRC / CAND_INC 指定。全库只有一份标量实现，所以锚点不需要改名。

**为什么不直接编工作树 src/**：另一个代理正在改它。baseline 永远来自 git show HEAD:。

make snapshot 在 HEAD 变了但内容没变时不动 mtime，所以不会多花一次重编译。

**为什么计时循环放在每个变体 TU 内部**，而不是由驱动用函数指针调用：

include/flashvec.h 正在被重构成头文件内联形态，默认就是 static inline。

内联形态下批量入口是 TU 私有的 static，外部根本链接不到；库形态下才是外部符号。

把循环放进变体 TU 后，两种形态都如实可测：

* 库形态：循环里是一次跨函数 call，也就是链接静态库的调用方拿到的形态。
* 内联形态：循环里被内联，也就是重构想要的形态。

驱动只认 ab_X_loop / ab_X_abi / ab_X_opflags 三类符号，nm 可查。

它完全不知道批量入口叫什么、有几个，所以**头文件怎么重构都不影响夹具**。

切换形态只需 make ... CAND_WITH_SRC=0。0 表示批量入口在头文件里，不编 .cpp。

**符号改名表是自动生成的**，脚本是 mkrename.sh。它从**本侧头文件**里抓出所有 `fv_xxx(` 形态的名字，生成 `#define fv_xxx ab_[ab]_xxx`。

好处是头文件新增任何入口都自动被隔离进命名空间，**不需要手工维护改名表**。旧的手写 ab_rename_a/b.h 已删除。

**踩过的坑**。四条都踩过，写在这里省下后来者两小时：

1. -D 宏改名必须用 -include 强制包含在源文件**第一行之前**。否则头文件里的声明和 .cpp 里的定义只改名一边，到链接期才发现。
2. `extern "C"` 的函数放进 namespace 不会改名，符号还是 C 名。所以「用命名空间隔离两套实现」这条路行不通。
3. 两套输出缓冲必须**连续分配，相距恰好 n*8 字节**，而 n*8 是 4096 的倍数。否则 A 与 B 落在不同的 store/load 4K 别名关系里，L1/L2 驻留档能凭空差 20%，这是实测值。
4. objcopy --prefix-symbols 会把 sqrtf / memcpy 这类未定义引用一起改名，结果是链接失败。只能逐符号 --redefine-sym，比宏改名麻烦，也容易漏。

## 3. 三道门禁（不过就不出数字）

1. **与独立参考实现逐位比较**。参考实现在 ab_bench.cpp 里手写：纯标量、纯顺序、不调用任何批量入口。
   比较级别是 **memcmp**，覆盖掩码「写回 a[i]」、idx 未引用位置与管线中间量。
   唯一例外是 fv_vec2_dot_soa。头文件明确不承诺逐位一致，按 |diff| <= max(1e-5*|结果|, 1e-6*Σ|逐项积|) 判等。
   **这个界与 n 和语料无关**。裸相对误差在结果接近抵消时会被放大，实测 n=2^18 就到 7.4e-5，第一版因此误报。
2. **A 与 B 的输出也必须逐位相等**，例外同上。计时结束后再比一遍。
   就地算子除外：它在同一个缓冲上反复累加，位级结论由计时前的 check_op 给出。
3. **输入缓冲的校验和在计时前后不变**，防止被计时的实现写坏 a/b/mask/idx。
   就地算子写的是自己的第一个输入。驱动把**它的副本**（O.o）当那个可写数组传进去，所以这条门禁对就地算子同样成立。
4. **--selfcheck**：baseline vs baseline 的比值必须落在 **[0.85, 1.15]**。不落就醒目报警，退出码 3。
   **这比任何加速比都重要。**

任何一条硬失败都会打印**函数名 + 首个不一致的元素**，**不输出任何计时/吞吐结果**，退出码 4。

make negctl 就是拿这个当断言用的，见 §7。

`--gate-only` 只跑门禁、不出计时，退出码 0/4。verify 用它做全尾长对拍。

## 4. 三个测量区制与口径

### 大 N 吞吐扫描（默认 2^12 … 2^24，共 8 档）

* 报 **ns/元素**、**名义 GB/s** 与 **真实 GB/s** 三套。名义 GB/s =（读字节+写字节）/ 时间；真实 GB/s 走 RFO 记账，见 §12。
* **B/元素逐函数人工核对**，不是「数组大小×常数」。这张表写在注册表 ab_ops.h 里，驱动只读表，表结构见 §9。
  实际数值用 `./build/ab_bench --list-ops` 查。
* **上限参考**（同规模、同一计时器）：memcpy（8R+8W）、memset（0R+8W）、纯读求和（8R）、triad（16R+8W，与 AoS 批量入口同形态）。
  批量函数以**占同规模 memcpy 吞吐的百分比**报出。读多写少的算子（dot / 纯读）**可以超过 100%**，因为读带宽本来就高于拷贝带宽。这不是错。
* **warm / cold 都报**。warm 是已驻留缓冲区反复跑；cold 是 madvise(MADV_DONTNEED) 后单次调用，**含 minor page fault**，页被内核清零。
  idx 数组永不丢弃。丢了之后所有倒排写都落到下标 0，随机写变成热点写，测出来是假数。
* **缓存层级按实测标**：L1d 32 KB、L2 512 KB、**L3 4 MB/CCX**。依据是 /sys/.../index3 的 size=4096K 与 shared_cpu_list=0-7。
  项目文档里「8 MB L3」是错的，一个线程只吃得到自己那 4 MB。扫描线里 2^17（3 MB，带内）与 2^18（6 MB，带外）正好压在这条边界两侧。

### 小 N 调用开销区制（默认 n = 4/8/16/32/64）

* 同时报 **ns/调用** 与 **ns/元素**。这一区制的瓶颈是单次调用开销，只报 ns/元素会把它摊薄到看不见。
* 覆盖哪几个算子由注册表的 `AB_OPF_SMALL` 标志决定，不是写死的下标。
* 对照组是 scalar_loop：逐元素调标量 API，且**循环与实现分处两个 TU**（ab_scalar_loop.cpp / ab_scalar_lib.cpp）。
  只有这样 GCC 才不能内联、不能向量化，才是「链接预编译静态库的调用方」的真实形态。
* 输出 A/scalar 倍数，直接回答「小 N 下批量入口还值不值」。

### 统计

* 每点**最少 5 次**取样取中位数，一直取到两侧变异系数都 < 5%（--cov）或到上限 --trials-max 为止。实际取样次数写进 TSV。
* **主判据 ratio = median(每轮 A/B)**：A 与 B 在同一段时间里配对，低频漂移一起消掉。
  同时给 ratio_medAB = median(A)/median(B) 作对照，以及 pair[min,max]、CoV、min/max。
* 默认**固定「先 A 后 B」**（--order=ab）。实测 2^20、9 次取样：ab → ratio 0.995 / CoV 9.7%，ba → 1.030 / 34.5%，alt → 1.305 / 24.3%。
  固定顺序去掉了「每个配对处在两种不同进入状态」这个混杂因素。ab 与 ba 都落在 1.0 附近，说明顺序本身没有系统性偏置。要交叉核对用 --order=ba。
* **计时基准默认用本线程 CPU 时间**（CLOCK_THREAD_CPUTIME_ID），不是墙钟：被别的进程抢占的时间不算进被测代码。
  它仍然包含邻居抢 L3 造成的真实停顿。--clock=mono 可切回墙钟自己核对。

## 5. 出口

    docs/research/raw/ab_baseline_head.txt         控制台全文（含矩阵 + RFO 记账 + 自检）
    docs/research/raw/ab_baseline_head.tsv         长表：一个 (phase,op,n) 一行，47 列
                                                   （32 基础 + RFO 记账 8 + 同轮 A/A 5 + mode/ref_op 2）
    docs/research/raw/ab_baseline_head_matrix.tsv  矩阵：phase × metric × op × 规模
                                                   （metric: gbs / real_gbs / pct_memcpy / ratio）
    docs/research/raw/ab_baseline_head_smallN.tsv  小 N 区制：n × op
    docs/research/raw/ab_aa_floor.tsv              A/A 地板表（见 §10）
    docs/research/raw/ab_campaign_verdict.txt      一轮覆盖所有新入口的判定表（campaign.sh）

make -C bench/ab run-worktree 会把 candidate 换成当前工作树，输出到 ab_vs_worktree.*。

> **raw 目录是每次 `make run` 覆盖的产物，不是不可变存档。**
>
> 每次重跑都会重写 `ab_baseline_head.{txt,tsv}`、`_matrix.tsv`、`_smallN.tsv` 与 `ab_aa_floor.tsv`。
>
> 文档里引用的数字（可分辨点数、逐档中位比、某点比值）会随之漂移。
>
> 实测同一台机器上「可分辨点数」在四次运行里分别是 **90 / 92 / 93 / 107**，分母都是 304；全点中位比在 1.0063 / 1.0070 / 1.0072 / 1.0079 之间。
>
> 引用这些数字时**必须连同时刻一起写**。`ab_baseline_head.txt` 的 header 有 label/rev，文件 mtime 就是那一轮的时刻。否则会把两个不同轮次的数字混着引用。
>
> 需要长期对比的结论请抄进 `docs/research/` 下的分析文档，不要直接引用 raw。

## 6. 命令速查

    make -C bench/ab run               自检 + 基线快照（8 档 + 小 N + A/A 地板表）
    make -C bench/ab run-worktree      candidate = 当前工作树
    make -C bench/ab verify OP=… N=…   同向证据五条判定（见 §11）
    make -C bench/ab aa-floor          单独重跑自检并落地板表
    sh bench/ab/campaign.sh            一轮覆盖 _stream/_mask_inplace/vec3/vec4（自带安静窗口判定）
    make bench-ab / bench-verify …     顶层 Makefile 的转发目标（只新增，不改既有目标）
    make -C bench/ab negctl            负对照
    make -C bench/ab build             只编译
    make -C bench/ab clean             只删 bench/ab/build

常用旋钮，都可从命令行覆盖：

    SIZES=4096,65536,1048576      大 N 扫描线
    SMALL=8,32                    小 N 取值
    TRIALS_MIN=5 TRIALS_MAX=12    自适应取样
    TARGET_MS=8 SMALL_MS=20       每取样目标时长
    CPU=6                         绑核（cpu0-7 同一 CCX）
    CAND_INC / CAND_HDR / CAND_SRC / CAND_WITH_SRC   换 candidate 源码树与形态
    AB_OPS_EXTRA=1                打开实验区（_stream / _mask_inplace）
    SELFCHECK=0                   关掉 [0.85,1.15] 带宽判定

## 7. 负对照（证明夹具会咬人）

    make -C bench/ab negctl

从 HEAD 快照复制一棵**私有**树，只改一处内核，编译成 candidate。断言三条：退出码 4、输出指名函数、没有任何计时结果。

三个 case：

1. _mm256_blendv_ps 的掩码选择方向选反 → 必须点名 fv_vec2_add_x4_mask
2. fv_scale4_avx2 里 mul 改 add → 必须点名 fv_vec2_scale_x4
3. normalize_soa 的 _CMP_NLE_UQ 取反 → 必须点名 fv_vec2_normalize_soa。
   这条**只因为**数据里每 512 个元素埋了一个精确零向量才会被抓到。零向量探针是刻意留的，不要删。

全程不碰仓库的 src/、include/，所以没有「跑完还原」这个动作：工作树从未被改过。

## 8. 已知局限（如实写）

* **本机长期被外部负载占用**，loadavg 常驻 4-9 / 16 核。L3 与 DRAM 档的取样变异系数常超 5%。
  这类点会被标 CoV! 并计入「不可分辨」，**照报数字但不进判据**。自检结论只对「可分辨」的子集成立，要用 1.15x 判据请优先看 L1/L2 档。
* cold 测量会把页清零，之后的 warm 计时跑在零数据上。对本库的批量入口无影响，因为它没有数据相关分支，尾巴长度也为 0。
  但如果将来加了带值相关分支的内核，就必须重灌数据。
* 纯读探针在小规模上是 load-port / ALU 受限的，不是带宽受限。它在 DRAM 档才代表带宽。
* 频率未固定：governor=performance，观测 1.5-3.9 GHz 随负载浮动。**绝对值不可跨会话比较**，只能看同一次运行内的比值。
* **代码放置**：同一个头文件的两种形态（库内定义 / 头文件内联）本身就会改变内联与对齐，L1/L2 档能差 5-10%。
  这不是「优化」，判据里靠「两套独立构建」（verify）与 A/A 地板一起把这类一次性因素关掉。
* 三维 `fv_vec3_*_soa` 与四维 `fv_vec4_*_x2` **还没有进注册表**。
  它们的元素宽度是 12 B / 16 B，与 fv_vec2 的 8 B 不同。Buf 的 n 语义与分配要一起改，见 §9 的「新布局」一条。
  驱动本身不用改，这就是表驱动的意义。

---

## 9. op 注册表：怎么加一个 op

唯一真源是 `bench/ab/ab_ops.h`。**加一个 op 只改两处**：

1. 在 `AB_LIB_OPS(X)` 里加一行，实验区则是 `AB_LIB_OPS_EXTRA(X)`；
2. 在 `ab_variant.cpp` 里加一个 `abt_<ID>` thunk，也就是那段调用代码。

其余全部自动派生，**驱动逻辑一行都不改**：

| 派生物 | 由谁展开 |
|---|---|
| `AB_OP_<ID>` 枚举（两侧 TU 与驱动共用同一编号） | ab_abi.h |
| `fv_xxx -> ab_a_xxx / ab_b_xxx` 改名表 | Makefile + mkrename.sh（从**头文件**抓） |
| OPS 表（名字/字节数/校验策略/参考实现/标志/RFO 类） | ab_bench.cpp |
| 变体 TU 的分派表 + `ab_X_opflags`（NT 与否） | ab_variant.cpp |
| 校验和门禁、小 N 覆盖面、矩阵、TSV、RFO 记账、verify | 全部读表 |

行格式是 10 列。NOTE 里不能有逗号，外部脚本按列解析：

```c
X(ADD_X4, "fv_vec2_add_x4", VEC2, 8, EXACT, ref_add_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,
  "a/b 读 16；o 写 8 触发 RFO 读 8 → 真实 32 B/e")
// ID    枚举名 AB_OP_<ID>（同时决定 thunk 名 abt_<ID>）
// NAME  报告/TSV 里的名字
// KIND  输出形态 VEC2 / SOA2 / SOA3 / SOA4 / SOA4V / VEC4 / FLOAT / NONE
// ELOP  **元素宽度**（字节）：vec2 8 / vec3 12 / vec4 16。同一档扫描里所有 op 处理相同
//       字节量（元素数 = n*8/ELOP），所以不同宽度 op 的比值天然按字节归一。
// POL   EXACT（逐位）/ REL（dot_soa 专用容差）/ NOCHECK（上限参考）
// REF   驱动内的独立标量参考实现
// FLAGS AB_OPF_PIPE / SMALL / INPLACE / NO_BASELINE（HEAD 里没有这个入口）；上限参考用 AB_MF_*
// RD/WR 名义每元素读/写字节（逐函数人工核对）
// RFO   AB_RFO_MISS / HIT / MIXED / LIBC / NT  —— 真实流量记账的分类
// NOTE  人读依据（会原样打印在"RFO 记账"表里）
```

vec3/vec4 的 13 个入口已经在实验区里，`AB_OPS_EXTRA=1` 时生效。

每行都带 ELOP 与 `AB_OPF_NO_BASELINE`，并在 `AB_INTRA_PAIRS(X)` 里登记「跟谁比」。

上限参考不经过变体 TU，走另一张表 `AB_MEM_OPS(X)`：8 列，无枚举、无改名。

**加一个 _stream 入口的完整例子**，这正是当前工作树里的形态：

```c
/* 1) ab_ops.h 里加一行（默认 182 项门禁不变，所以要放实验区） */
X(ADD_X4_STREAM, "fv_vec2_add_x4_stream", VEC2, EXACT, ref_add_x4, 0u, 16, 8, AB_RFO_NT,
  "NT store：16 读 + 8 写 = 24 B/e（普通版 32 B/e 的那 8 B RFO 没了）")

/* 2) ab_variant.cpp 里加 thunk（名字必须是 abt_ + ID） */
static void abt_ADD_X4_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM   /* Makefile 自动 grep 本侧头文件 */
    fv_vec2_add_x4_stream((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#else               /* baseline 头文件还没有它：退化成普通 store，于是这一行就是
                     * "普通 store vs NT store" 的对照（32 vs 24 B/e） */
    fv_vec2_add_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#endif
}
```

然后：

    make -C bench/ab AB_OPS_EXTRA=1 build      # 编译
    ./bench/ab/build/ab_bench --list-ops       # 确认在表里、字节数对
    make -C bench/ab verify OP=fv_vec2_add_x4_stream N=1048576

**试验区的开关**：`AB_OPS_EXTRA=0` 是默认值，此时 `AB_LIB_OPS_EXTRA` 展开为空，所以默认门禁仍是历史上那 182 项，即 15 op × 8 档 + 6 op × 5 小 N。

新入口稳定后想让它进默认矩阵，把 ab_ops.h 里那一行从 EXTRA 表挪到 AB_LIB_OPS 即可。

**新布局（vec3 / vec4）**：元素宽度不是 8 B 时，除上面两处外还要扩一次 `AbArgs` + `Buf`，加上 az/bz/oz 这类指针字段与元素宽度规则。

那次改动是一次性的，与 op 数量无关。之后所有同布局的新 op 仍然只加表和 thunk。

## 9.5 这台机器的外部负载是**常驻**的 —— 绝对 loadavg 不能当测量门

**证据**，同一时刻采集：

    $ cut -d' ' -f1-3 /proc/loadavg     ->  3.34 3.87 4.69
    $ pgrep -c cc1plus                  ->  0
    $ pgrep -c 'cc1 '                   ->  0

`loadavg` 常驻 3–9。本轮收尾 campaign 期间从 **3.03 涨到 5.41**，而**全过程没有任何 cc1plus / cc1 进程**。

也就是说这些负载**不是我们的编译器**造成的，是机器上的常驻外部进程。

**结论句**：**这台机器上绝对 loadavg 不能当测量门，必须用 A/A 地板控制。**

由此三条操作规则：

1. **「安静窗口」的判据是「无我方编译进程」**，不是 loadavg：
   `pgrep -f 'cc1plus|cc1 |make -C bench/ab|g\+\+ .*flashvec'` 连续 60 s 为空。
   campaign.sh 用的就是这个。本轮等到干净窗口只用了 5 s；旧判据 `loadavg < 1.5` 永远打不开。
2. **loadavg 只记录、不拒绝**：记下起止值，只有**涨幅 > 1.0** 时才提示「这段测量可能被咬」。
3. **用测量控制测量**：污染对每个 (op,size) 点的影响由**同轮交错的 A/A 地板**直接量化，规则是「地板 > 15% → 该点 INCONCLUSIVE」。
   本轮 campaign 的 loadavg 涨了 2.38，但没有贴任何全局「污染」标签，而是有 **38 个点**因为自身地板 > 15% 被判 INCONCLUSIVE。
   一个全局标签会把干净的点一起废掉，逐点地板不会。

## 10. A/A 地板表

`docs/research/raw/ab_aa_floor.tsv` 存的是每个 (phase, op, n) 的**自检比值**，也就是「两侧跑同一份代码」时这台机器能给到的分辨率：

    phase  op  n  ratio  pair_iqr_pct  trials  resolvable  floor_dev  in_band  rfo_flags

* `resolvable=0` 的点连地板都测不准（配对 IQR >= --cov），**不要拿它判**；
* `floor_dev = |ratio-1|`：任何「候选快 r 倍」的结论，只有 `(r-1) > floor_dev` 才可能为真；
* `in_band`：比值是否落在 [0.85, 1.15]，也就是自带的那条自检带；
* 本表只在**两侧同一份源码**的构建里成立。默认 `make run` 就是这种情况；run-worktree 不是，所以它不写这张表。

## 10.5 INTRA 模式：没有基线的新入口怎么测

`fv_vec4_*_x2` / `fv_vec3_*_soa` 这些入口在 **git HEAD 里根本不存在**。

拿它们跟「不存在」比，得到的是「从无到有」，不是「快了多少」。所以注册表给它们打上 `AB_OPF_NO_BASELINE`，并且：

* 普通 A/B 模式**拒绝**把它们算进来；--op 点名时自动切 INTRA。
* `--intra` 模式在 **candidate 内部**比：新入口 vs 注册表 `AB_INTRA_PAIRS` 里指定的同字节量参照入口。
  配对是 vec4 ↔ vec2_x4、vec3 ↔ vec2_soa、vec4_soa ↔ vec2_soa。
* **按字节归一**：同一档扫描里所有 op 处理相同字节量，元素数 = n*8/元素宽度。所以 ratio = t(ref)/t(new) 就是「同一份工作量谁更快」。
* 报告里同时给 ns/元素 与 GB/s 两套，并**显式标注** `INTRA（candidate 内部对照，不是相对 HEAD 的改进）`。

元素宽度写在注册表的 **ELOP 列**，取值 8 / 12 / 16，不是从输出形态推的。原因是 dot/len 的输出是 float/NONE，元素宽度仍可能是 12/16 B。

## 10.6 结论三档：PASS / FAIL / INCONCLUSIVE

**FAIL 会被误读成「这个优化无效」，而实际往往是「这台机器测不出来」。**两者必须分开：

| 档 | 触发条件 |
|---|---|
| FAIL | R1 位级对拍失败；R5 其它规模回退 > 5%；中位比显著低于 1 且超过地板（测出来的变慢） |
| INCONCLUSIVE | 增益 ≤ A/A 地板；A/A 地板 > 15%（该点不可分辨）；同向但未达 1.15；黑名单档（n=131072） |
| PASS | R1 过 + 四个比值都 ≥ 1.15 + 最小增益 > A/A p90 地板 + 地板 ≤ 15% + R5 无回退 |

## 10.7 逐规模判据：地板制，不是常数带

`--selfcheck` 里「每一档的比值中位数」原来必须落在写死的 `[0.90, 1.10]`。

实测 4 次运行里有 2 次因此 exit 3，而**校验和门禁 182/182 每次都过**。失败的是这个常数带，不是夹具。改成：

1. 判据 = **该档实测 A/A 地板 p90**。它由同一批运行里的 `--aa-each-round` 测出，属同一分布。
   即 `|该档中位比 - 1| <= 该档地板 p90`。常数带只作为「没有同轮 A/A 数据时」的退化路径。
2. **n=131072 从逐规模判据里剔除**（`--tier-exclude=131072`，可改）。它是 L3 边界，A/A 地板本身在 0.50–2.22 之间，拿它判等价于掷骰子。
   仍然测、仍然记录、**仍然参与门禁**。
3. **兜底**：某档地板本身 > 20% → 该档标 INCONCLUSIVE，不 PASS 也不 FAIL。因为「测不准」和「测出来偏了」是两件事。

## 10.8 负结果：`_mask_inplace` 没有可测的吞吐优势（5 密度 × 4 规模全 INCONCLUSIVE）

campaign 实测 `fv_vec2_add_x4_mask_inplace`：密度 1% / 5% / 20% / 50% / 95%，规模 8192 / 65536 / 1M / 4M。

**20 个点全部 INCONCLUSIVE**，增益 **-0.10 … +0.00**，同轮 A/A 地板 1.06–1.56。

**原因在流量账上**。就地写是 `a[i] += b[i]`，目标数组在写之前**已经被读过**，对应的 cache line 本来就在 L1/L2 里。

普通 store 在这种情况下**不产生 RFO**。所以「就地」相对「分离输出」少掉的不是 8 B/元素，而是 0。

注册表里这一行的 RFO 类别就是 `AB_RFO_HIT`，真实流量 25 B/e = 名义 25 B/e；分离版是 `AB_RFO_MISS`，真实 33 B/e。

稀疏度只改变「写多少元素」，不改变「每个被写元素是否要 RFO」，所以密度扫描也扫不出差别。

**对夹具的启示**，这才是重点：一个看起来天经地义的优化假设——「稀疏就地写更快」——在流量账上根本不成立。

**先建夹具、再动手优化**，否则这个入口会带着「应该更快」的直觉被写进契约，再被后续的人当成事实引用。

夹具的价值不是「证明优化有效」，而是**在写代码之前就把无效的优化假设挑出来**。

## 11. 五条验收判据（verify 子命令）

`make -C bench/ab verify OP=<op> N=<规模>` 一次跑完并打印判定表：

| 规则 | 判据 | 数据来源 |
|---|---|---|
| R1 | 位级对拍（全尾长 n=N-7…N+3，两套构建） | `--gate-only` 两遍 + 四次计时运行的门禁 |
| R2 | 目标规模同进程中位比 >= 1.15x | ratio(AB, build#1) |
| R3 | **同向复现**：AB/BA（换顺序）× 两次独立构建（换代码放置）四个比值全部 >= 1.15x 且同向 | 4 次运行 |
| R4 | 增益 (r-1) **超过**该点 A/A 地板 | **同轮交错的 A/A**，见下 |
| R5 | 其它规模不回退超过 5% | 一次多规模扫描；L2/L3 的预期性回退用 `EXPECT_REGRESS=65536,262144` 显式声明 |

R4 的地板这样测：每个 A/B 轮后面紧跟一对 A/A，同一进程、同一数据、同一浸入状态。

跨 2 构建 × 2 顺序 × ROUNDS 轮，≥16 个样本，取 p50/p90/max，判据用 p90。

FAIL 时会打印**具体是哪一条、差多少**，例如 `ratio=0.999 < 1.15（同进程同规模没有可测增益）`，不给红叉了事。

旋钮：`MIN_GAIN`（默认 1.15）、`MAX_REGRESS`（默认 0.95）、`CHECK_SIZES`、`EXPECT_REGRESS`。

**两次构建怎么算「独立」**：同一份源码，只改代码放置（`-falign-functions=64 -falign-loops=32 -falign-jumps=32`），语义零变化。

这一条专门用来关掉「代码放置」这类一次性因素，它是 A/A 地板的主要来源。

**为什么 A/A 必须同轮交错（R4 的修法）**：旧版拿「另一组运行的地板」比「本组运行的增益」，两者不是同一个分布。

现在每个进程里每轮先测 A/B 再测一对 A/A（`--aa-each-round`），地板与增益共享同一段频率漂移、同一份数据、同一个浸入状态。

实测效果：同一个无差异 op 上，旧版 A/A 能给出 1.208 并「通过」R4；新版把该点的地板 p90 定成 0.0805、增益 -0.2393，直接判 INCONCLUSIVE，因为增益在地板以下。

**n=131072 是黑名单档**：实测该档 A/A 比在 0.50–2.22 之间（L3 边界 + 4K 别名），可以测、可以记录，**不得当判据**。

verify 与 campaign 都按 BLACKLIST 处理。

## 12. RFO 记账：名义流量 vs 真实流量

GB/s 原来用的是**名义字节**（rd+wr）。但普通 store 在写 miss 时会先把整条 cache line **读进来**。

这是 RFO，read-for-ownership，真实流量是「写一次 + 读一次」：

    real_wr = wr * 2   普通 store（写分配 miss）
    real_wr = wr        NT store（_mm256_stream_ps：不经过 cache、不产生 RFO）
    real_wr = wr       就地更新（写目标之前已经被读过，行本来就在 cache 里）
    real_wr = wr       libc 自己选（memcpy/memset：标注 libc?，不硬填）

真实 GB/s = (rd + real_wr) / ns_per_elem，与名义 GB/s 并列输出。矩阵里 `gbs` 与 `real_gbs` 各占一个 metric。

**这条的物理意义**：NT store 的真实收益不是「每个字节更快」，而是**少走 25% 的字节**。同一份 add 算子名义都是 24 B/e：

| | 名义 B/e | 真实 B/e | 真实 GB/s | 名义/真实 |
|---|---|---|---|---|
| 普通 store | 24 | 32 | 同一值 | 0.75 |
| NT store | 24 | 24 | 同一值 | 1.00 |

所以两者的**真实** GB/s 应该基本相等，都在 DRAM 带宽上。而**名义** GB/s 差 1.33 倍。

如果名义差 1.33 而真实也差 1.33，那说明增益不是流量导致的。这就是判假阳性的一把尺子。

反过来说：NT store 只在**工作集远超 L3、写流是干净的一次性写**时才有这 25%。工作集落在 L2/L3 里时它比普通版慢。

**实测 0.13–0.51x，即慢 2–8 倍**。n=8192 时六个入口全部判 FAIL（变慢）：

| 入口 | 比值 |
|---|---|
| add | 0.217 |
| sub | 0.301 |
| scale | 0.130 |
| add_soa | 0.229 |
| scale_soa | 0.138 |
| normalize | 0.506 |

这组数字来自 campaign 判定表 `docs/research/raw/ab_campaign_verdict.txt`。

头文件契约里旧的 0.40x 是单个入口的单点测量，低估了这个惩罚的幅度。

不确定的点**标注而不是硬填**：

| 标注 | 含义 |
|---|---|
| `mixed?` | 管线中间数组可能被 L3 吸收，2x 是上界 |
| `libc?` | memcpy/memset 自己选 NT |
| `miss*` | 本侧没有 NT 入口、退化成普通 store |

`--list-ops` 会把每个 op 的 rd/wr/真实字节数与分类打出来。

## 13. 这套判据还有什么漏洞（自己先说出来）

1. **A/A 地板只度量 A 侧的噪声。** 地板来自 A 与 A 的比值，而 A/B 比值的抖动里还含 B 侧——另一个形态、另一份目标文件布局。
   两侧噪声不同时地板会低估 A/B 的真实抖动，所以它是**下界**，不是估计值。
2. **地板本身在轮次之间会翻倍。** 实测 65536 档三轮的 A/A 地板分别是 0.130 / 0.169 / 0.218，同一档能在「可判」与「INCONCLUSIVE」之间翻转。
   因此「某一档可判」这个结论本身就不稳定；只有**逐点 PASS 且该点地板小**才算硬结论。
3. **R5 只看 5 个规模**（默认），可能漏掉中间某个规模的坑。
4. **libc 与 mixed 类 op 的真实流量是区间而不是数**，本表只给端点并标 ?。任何基于 real_gbs 的「精确」结论都不成立。
5. **A/A 地板依赖构建形态**：基线是库内定义形态、候选是头文件内联形态时，两侧「同一份代码」其实也不是同一个目标文件布局，地板里混进了形态差异。
   verify 的两次构建只换了代码放置，没换形态，所以这一点**没有**被覆盖。
6. **verify 的 R1 只覆盖这个 op**，不是整个库的门禁，那归顶层 make test-all。它回答的是「这个改动没有破坏这个 op 的位级语义」。
