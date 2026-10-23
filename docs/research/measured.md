# flashvec 三档实现实测基准（AMD Ryzen 7 5700U / Zen 2 / AVX2+FMA）

> 目标：把"标量 AoS API / 编译器自动向量化 / 手写 AVX2"三档放在**同一份负载**下实测，给出每档的 ns/op、
> 相对档 1 的加速比，以及**为什么是这个数**（带指令级证据）。全部数字由 `make bench` 在本机真实跑出，
> 原始日志在 [raw/](raw/)；本文所有代码块都是原始输出的原样粘贴。

---

## 0. 结论速览（先看这 8 条）

1. **N = 1<<20（每数组 4 MB，工作集 24 MB > 8 MB L3）时，三档没有差别。** add 档1 = 5.33、档2 = 5.48、
   档3 = 6.77 ns/elem（x0.97 / x0.79）；只有 collision（x1.34–1.41）和 normalize（x1.06–1.14）有稳定收益。
   原因不是"SIMD 没用"，而是三档都撞在这台机器的**内存带宽墙**上：实测 triad 只有 5.3–7.0 GB/s、
   memcpy 9.1–17.5 GB/s（§3.3）。→ **在 1<<20 这个规模上，自动向量化和手写 AVX2 都白写。**
2. **N = 8192（每数组 32 KB，L1/L2 驻留）时差距立刻出来**，档3 手写 AVX2 相对档1：add **x6.54**、
   scale **x13.23**、dot_elem **x8.98**、len **x9.63**、normalize **x8.41**、collision **x14.51**。
3. **但档2（朴素 for + 自动向量化）几乎追平档3**：add x6.30 vs x6.54、scale x13.46 vs x13.23、
   dot_elem x9.46 vs x8.98。也就是说：**只要数据是 SoA、运算里没有条件分支/跨 lane 归约，
   手写 AVX2 相对 `-O3 -march=native` 没有净收益（差 0–8%，方向还不稳定）。**
4. **档2 的向量化被两件事卡住，代价完全不同**（GCC 原文裁决见 §3.5）：
   - `len`（循环里有 `sqrtf` 的 errno 分支）：GCC 报 `missed: not vectorized: control flow in loop.`。
     加 `-fno-math-errno`（**不改变任何数值结果**）后：x1.25 → **x8.16**（0.2395 ns/elem）。
   - `normalize`/`collision`（循环里有**条件执行的除法**，即零向量保护）：加 `-fno-math-errno` 仍不向量化；
     只有 `-ffast-math` 能解锁：x1.14 → x9.84、x2.01 → x13.53。
   - 手写 AVX2 用 `_mm256_cmp_ps` + `_mm256_and_ps` 写成**无分支掩码**，在不放宽任何 FP 语义的前提下
     拿到 x8.41 / x14.51，并且与标量路径**逐位一致**（§3.4，`bits_equal=1`）。
5. **dot 的瓶颈是累加器依赖链，不是水平归约指令。** 档3 单累加器 0.3213 ns/elem（x8.34）→ 4 个独立累加器
   **0.1428**（x18.78）；而被认为"慢"的 `vhaddps`（每 8 元素归约一次）实测 0.2492，
   **比单累加器的 shuffle 版还快 22%**。→ 批量 dot 的正确写法是「4 累加器 + 末端一次 shuffle 归约」。
6. **档1 的代价可以逐指令解释**：主循环每元素 `vmovq → vmovq → call → vmovq`，参数确实走 xmm 寄存器
   （**没有栈拷贝，ABI 不背锅**）；`dot` 还把标量累加器经 `ebp` 来回搬（§3.6）。6–13 倍差距来自
   「每元素一次函数调用 + 完全没有 SIMD」。
7. **`-ffast-math` 的代价是可测量的**：它把 8 个开关全开（§2.2 实测展开），换来档2 在
   len/normalize/collision 上的 7–14 倍；数值上 tier1/tier2/tier3 与严格 IEEE 参考的差 ≤1.79e-07 绝对、
   ≤4.95e-05 相对（出现在结果接近 0 的抵消处）。rsqrt+NR 路径相对误差 ≤2.62e-07，8192 个元素里 3146 个位模式不同。
8. **perf stat 在本机拿不到数据**：`kernel.perf_event_paranoid = 4`、无 CAP_PERFMON、`sudo` 需要密码，
   `perf stat` 退出码 255（原文见 §5）。本报告用「GCC 向量化裁决 + 反汇编 + 分级内存带宽实测」替代。
   一条命令补齐：`sudo sysctl -w kernel.perf_event_paranoid=1 && make bench-perf`。

---

## 1. 对 flashvec 的具体建议（每条都写清代价）

### 1.1 零语义代价、当天可做：批量入口所在编译单元加 `-fno-math-errno`
- 现状：`Makefile` 的 `CXXFLAGS ?= -std=c++17 -O2 ...`，无 `-O3`、无 `-march`、无 `-fno-math-errno`。
- 实测：同一条 len 循环，不加 → x1.25（1.56 ns/elem）；加 `-fno-math-errno` → **x8.16**（0.2395 ns/elem），且结果逐位不变。
- 代价：只放弃"数学函数失败时设置 `errno`"。对 `sqrtf` 的非负输入而言 errno 从不被设置，风险≈0。

### 1.2 批量 API 形态（收益最大的一条：让编译器有得向量化）
现在只有按值小结构体 API，用户写 `for (i) o[i] = fv_vec2_add(a[i], b[i]);` 时因为跨编译单元不能内联，
**每个元素一次 call**（§3.6 证据），且 AoS 的 8 字节元素对 AVX2 不友好。建议在 `include/flashvec.h` 新增：

```c
/* 批量入口：输入输出都是分离分量（SoA）。内部无分支、无 FV_EPS 保护。 */
void   fv_vec2_add_soa(const float *ax, const float *ay, const float *bx, const float *by,
                       float *ox, float *oy, size_t n);
void   fv_vec2_scale_soa(const float *ax, const float *ay, float s, float *ox, float *oy, size_t n);
float  fv_vec2_dot_soa(const float *ax, const float *ay, const float *bx, const float *by, size_t n);
void   fv_vec2_len_soa(const float *ax, const float *ay, float *ox, size_t n);
/* 归一化：零向量保护用掩码而不是分支；调用者显式接受"长度 <= eps 时输出 0"的语义。 */
void   fv_vec2_normalize_soa(const float *ax, const float *ay, float *ox, float *oy, size_t n);
```

- 收益（实测 N=8192）：add x6.5、scale x13.2、dot x18.8（4 累加器）、len x9.6、normalize x8.4。
- 代价：ABI 表扩大；调用者必须把数据拆成分量数组——对引擎是真成本：若内部是 AoS，
  拆/装（`vpermps`/`vunpck`）会吃掉一部分收益（**本报告未实测 AoS↔SoA 转换成本，标为未验证**）。
- 折中（未实测，建议下一步做）：保留 AoS，加 4 路交错入口 `void fv_vec2_add_x4(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o, size_t n)`，
  一次吃 4 个 `fv_vec2`（正好一个 ymm）；代价是要处理 `n % 4` 尾巴。

### 1.3 `-march=native` 的部署代价（必须说清）
- 本机展开为 `-march=znver2 -mtune=znver2`（§2.2），产出的二进制在**没有 AVX2 的机器上会 SIGILL**。
- 若要发预编译库：用 `__attribute__((target_clones("default","avx2")))` 或运行时 `__builtin_cpu_supports("avx2")` 分派。
  代价：代码体积、维护面、一次分派分支（**本报告只验证了能编译，未实测分派开销**）。
- 本次基准确认：**AVX-512 在本机完全禁用**（`-mavx512f/bw/dq/vl...` 全 [disabled]），按 16 float/次 的手写内核在这台机器上不会更快。

### 1.4 dot / len / normalize 的内核写法（有实测支撑）
| 运算 | 推荐写法 | 实测依据 |
|---|---|---|
| dot（归约） | **4 个 ymm 累加器**，末端一次 extract+shuffle 归约 | 单累加器 0.3213 → 4 累加器 0.1428 ns/elem（x2.25） |
| dot（长数组） | 同上；别每 8 元素归约一次 | 每 8 元素 hadd 归约 0.2492，比 4 累加器慢 74% |
| len | 批量写 `_mm256_sqrt_ps`，掩码清零零向量 | x9.63；与标量路径逐位一致 |
| normalize | 批量写 `_mm256_div_ps(1, l)` + `_mm256_and_ps` 掩码 | x8.41，**逐位一致** |
| normalize（快路） | 若接受 ~1e-7 相对误差，用 `rsqrtps`+1 次 NR | 0.3082 vs 0.3395 ns/elem，但 3146/8192 元素位不同 |

### 1.5 三维 / 矩阵扩展的 ABI 提醒（实测）
- `fv_vec2` = 8 字节，SysV 下**一个 xmm 寄存器**传参（§3.6 B1）。
- `fv_vec3` = 12 字节，实测**两个 xmm 寄存器**传参（§3.6 B7），第二个只用了 4 字节、浪费 4 字节。
  → 扩三维时建议直接定义 `typedef struct { float x, y, z, w; } fv_vec4;`（16 字节，1 个 ymm 装 2 个向量）；
  代价：多 4 字节内存、要与 `fv_vec3` 互转。

### 1.6 不要做的事
- 不要给 `fv_vec2_add` 这类单元素函数套 `if (__builtin_cpu_supports)` 做 AVX2 分派：单元素调用永远被 call 开销主导（§4.6）。
- 不要用 `-Ofast` 而不说明放宽了什么；本报告用 `-O3 -march=native -ffast-math` 并给出实测展开（§2.2）。
- 不要指望「把 `sqrtf` 换成 `_mm_sqrt_ss`」就能让档1 变快：档1 的瓶颈是调用与标量化。

---

## 2. 测试环境、口径与负载

### 2.1 机器（原始采集见 [raw/env.txt](raw/env.txt)）
| 项 | 值 |
|---|---|
| CPU | AMD Ryzen 7 5700U with Radeon Graphics（Zen 2，8 核 16 线程） |
| 指令集 | AVX2 + FMA 可用；**AVX-512 全部 disabled**（`-march=native` → `znver2`） |
| L1d / L2 / L3 | 256 KiB（8×32 KiB） / 4 MiB（8×512 KiB） / 8 MiB（2×4 MiB） |
| 编译器 | g++ (Ubuntu 12.4.0-2ubuntu1~24.04.1) 12.4.0 |
| perf | perf 6.17.13 存在，但被 `perf_event_paranoid=4` 挡住（§5） |
| 调速器 | `performance`；实测运行中 scaling_cur_freq 在 1.70–4.22 GHz 间变动 |
| 机器负载 | **共享机器**：全程 loadavg 3.0–5.7（同机还有其他编译/测试进程）；已绑核 cpu2 |

### 2.2 编译选项与实测展开
| cfg | 选项 |
|---|---|
| `O3_native` | `-O3 -march=native` |
| `O3_native_ffastmath` | `-O3 -march=native -ffast-math` |
| `O3_native_nomatherrno` | `-O3 -march=native -fno-math-errno`（严格语义安全，用来隔离原因） |
| `O2_shipped` | `-O2`（仓库当前 Makefile 的档1 编译选项） |

`-march=native` 实测展开（`g++ -Q --help=target`）：

```
-march=   znver2        -mtune=  znver2
-mavx2    [enabled]     -mfma    [enabled]
-mavx512f [disabled]    -mavx512bw [disabled]   -mavx512dq [disabled]  -mavx512vl [disabled]
```

**`-ffast-math` 到底放宽了什么**（`g++ -O3 -ffast-math -Q --help=optimizers` 原样输出，不是抄文档）：

```
-fassociative-math         [enabled]   ← 允许重排浮点加法链（归约变向量累加器的前提）
-ffinite-math-only         [enabled]   ← 假定没有 NaN/Inf
-freciprocal-math          [enabled]   ← 允许把 x/y 变成 x*(1/y)
-funsafe-math-optimizations [enabled]
-fmath-errno               [disabled]  ← 不再为数学函数维护 errno
-frounding-math            [disabled]  ← 假定默认舍入模式
-fsigned-zeros             [disabled]  ← 不再区分 +0.0 / -0.0
-ftrapping-math            [disabled]  ← 允许推测执行可能触发 FP 异常的运算
-ffp-contract=fast                    ← 默认就是 fast：允许把 a*b+c 收缩成 FMA
```

### 2.3 测量口径
- 负载：6 个运算（add / scale / dot / dot_elem / len / normalize）+ 1 个物理微内核（collision）；
  N ∈ {8192, 65536, 1048576}（SoA；档1 用 AoS 的 `fv_vec2` 数组）。
- 计时：`clock_gettime(CLOCK_MONOTONIC)`；每用例重复 reps 次取平均（N=1<<20 时 reps=100，按任务书要求），
  重复 5 次（trials），**报告值是 5 次中位数**，同时把 5 次原始值和最小值都写进日志（`trials_ns=` / `min_ns=`）。
- 防优化：每次重复之间插 `asm volatile("" ::: "memory")` 屏障；计时区**外**单独跑校验和（读输出数组），防死代码消除。
- 绑核：`sched_setaffinity` 到 cpu2（日志首行 `# pinned cpu=2 running_on=2`）。
- 精度对照：`bench/ref.cpp` 单独用 `-O2 -fno-fast-math -ffp-contract=off -fno-tree-vectorize` 编译作为严格 IEEE 参考；
  每用例报 `max_abs/max_rel/n_diff/n_nan`，并额外做逐位比较（`IDENT bits_equal=`）。

### 2.4 噪声与可信度（必须先说）
- 本机是共享环境（同机还有其他进程，loadavg 3–5.7）。同一份代码两次运行的 N=1<<20 数据可差 20%+：
  例如档3 add 第一次 5.28 ns/elem、第二次 6.77 ns/elem。
- **因此：N=1<<20 的表只用于说明"三档撞在带宽墙上"这一定性结论**；定量结论请用 N=8192 / 65536 两张表（两次偏差 <8%）。
- 每行日志都带 `trials_ns=`，抖动可自查：例如档1 collision@1M 的 5 次是
  `9679309 / 9289511 / 10827762 / 10868032 / 10019166` ns（中位数 10019166、最小值 9289511，差 8%）。

---

## 3. 实测结果

### 3.1 主表：三档 × 三规模（ns/elem，5 次中位数；加速比相对档1@O3_native）

**N = 8192（每数组 32 KB，L1/L2 驻留）** —— 最能说明"三个档位本身"的差别：

| 运算 | 档1 标量 AoS | 档2 自动向量化 | 档2 +fast-math | 档2 +no-math-errno | 档3 手写 AVX2 | 档3 变体 | 档3 GB/s |
|---|---|---|---|---|---|---|---|
| add | 1.5162 | 0.2405 (**x6.30**) | 0.2370 (x6.40) | 0.2392 (x6.34) | 0.2319 (x6.54) | – | 103.5 |
| scale | 2.2139 | 0.1645 (**x13.46**) | 0.1672 (x13.24) | 0.1674 (x13.23) | 0.1674 (x13.23) | – | 95.6 |
| dot | 2.6811 | 0.7607 (x3.52) | 0.1448 (x18.52) | 0.7646 (x3.51) | 0.3213 (x8.34) | 4累加 0.1428 (**x18.78**) | 49.8 |
| dot_elem | 1.8714 | 0.1979 (**x9.46**) | 0.2136 (x8.76) | 0.1994 (x9.39) | 0.2083 (x8.98) | – | 96.0 |
| len | 1.9551 | 1.5597 (x1.25) | 0.2750 (x7.11) | 0.2395 (**x8.16**) | 0.2031 (**x9.63**) | – | 98.5 |
| normalize | 2.8541 | 2.4967 (x1.14) | 0.2901 (x9.84) | 2.4844 (x1.15) | 0.3395 (x8.41) | rsqrt 0.3082 (x9.26) | 47.1 |
| collision | 7.9584 | 3.9682 (x2.01) | 0.5882 (x13.53) | 3.8351 (x2.08) | 0.5483 (**x14.51**) | rsqrt 0.5552 (x14.33) | 58.4 |

**N = 65536（每数组 256 KB，L2/L3 交界）**：

| 运算 | 档1 | 档2 | 档2 +fast-math | 档2 +no-math-errno | 档3 | 档3 变体 | 档3 GB/s |
|---|---|---|---|---|---|---|---|
| add | 1.4260 | 0.2923 (x4.88) | 0.3266 (x4.37) | 0.3146 (x4.53) | 0.2992 (x4.77) | – | 80.2 |
| scale | 2.1721 | 0.2206 (x9.85) | 0.2197 (x9.89) | 0.2263 (x9.60) | 0.2228 (x9.75) | – | 71.8 |
| dot | 2.6521 | 0.7721 (x3.43) | 0.1768 (x15.00) | 0.8172 (x3.25) | 0.3385 (x7.83) | 4累加 0.1836 (x14.44) | 47.3 |
| dot_elem | 1.8067 | 0.2728 (x6.62) | 0.2943 (x6.14) | 0.2742 (x6.59) | 0.2703 (x6.68) | – | 74.0 |
| len | 1.8037 | 1.5900 (x1.13) | 0.3126 (x5.77) | 0.2657 (x6.79) | 0.2294 (**x7.86**) | – | 87.2 |
| normalize | 2.7905 | 2.5382 (x1.10) | 0.3021 (x9.24) | 2.5544 (x1.09) | 0.3751 (x7.44) | rsqrt 0.3691 | 42.7 |
| collision | 7.9121 | 4.7021 (x1.68) | 0.6592 (x12.00) | 4.5492 (x1.74) | 0.6354 (x12.45) | rsqrt 0.6898 | 50.4 |

**N = 1048576（每数组 4 MB，工作集 24 MB > L3 8 MB）** —— 三档撞在内存墙上：

| 运算 | 档1 | 档2 | 档2 +fast-math | 档2 +no-math-errno | 档3 | 档3 GB/s |
|---|---|---|---|---|---|---|
| add | 5.3275 | 5.4777 (x0.97) | 5.4674 (x0.97) | 6.2021 (x0.86) | 6.7682 (x0.79) | 3.55 |
| scale | 4.1889 | 4.3238 (x0.97) | 4.0990 (x1.02) | 4.1111 (x1.02) | 4.0064 (x1.05) | 3.99 |
| dot | 3.5304 | 2.5921 (x1.36) | 2.5095 (x1.41) | 2.5752 (x1.37) | 2.5957 (x1.36) | 6.16 |
| dot_elem | 5.3121 | 5.3256 (x1.00) | 5.3644 (x0.99) | 5.7004 (x0.93) | 5.3940 (x0.98) | 3.71 |
| len | 4.2208 | 4.0837 (x1.03) | 3.9290 (x1.07) | 4.0235 (x1.05) | 4.0430 (x1.04) | 4.95 |
| normalize | 4.5879 | 4.3124 (x1.06) | 4.0814 (x1.12) | 4.2673 (x1.08) | 4.0104 (**x1.14**) | 3.99 |
| collision | 9.5550 | 7.1059 (x1.34) | 6.9829 (x1.37) | 7.2374 (x1.32) | 6.7836 (**x1.41**) | 4.72 |

> 档1 的 `dot_elem` 与档2/档3 的语义统一为"逐元素点积写成标量数组"（`ox[i]=oy[i]=a_i·b_i`）；
> 库没有逐分量乘，所以档1 这一行是"每元素一次 `fv_vec2_dot` 调用"。

### 3.2 原始输出（`make bench` 产物，原样粘贴）

```
$ make bench            # = make -C bench run：编译 4 组 cfg，跑 3 规模 × 3 档 × 7 运算，落盘 docs/research/raw/
$ make bench-matrix     # 由原始日志生成下面的紧凑矩阵
$ make bench-table      # 由原始日志生成宽表
$ make bench-checks     # 只跑校验段（DIFF/IDENT）
$ make bench-vec        # GCC 向量化报告
$ make bench-asm        # 反汇编归档
$ make bench-perf       # perf stat（需先放开 perf_event_paranoid）
```

`docs/research/raw/matrix.txt`（原样）：

```
-- N=1048576 --
add       t1=5.3275   t2=5.4777   t2ff=5.4674   t2noerr=6.2021   t3=6.7682   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=3.55  | x0.97 x0.97 x0.86 x0.79 -
scale     t1=4.1889   t2=4.3238   t2ff=4.0990   t2noerr=4.1111   t3=4.0064   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=3.99  | x0.97 x1.02 x1.02 x1.05 -
dot       t1=3.5304   t2=2.5921   t2ff=2.5095   t2noerr=2.5752   t3=2.5957   t3hadd=2.4619   t3acc4=2.4110   t3rsqrt=-        t3_gbps=6.16  | x1.36 x1.41 x1.37 x1.36 -
dot_elem  t1=5.3121   t2=5.3256   t2ff=5.3644   t2noerr=5.7004   t3=5.3940   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=3.71  | x1.00 x0.99 x0.93 x0.98 -
len       t1=4.2208   t2=4.0837   t2ff=3.9290   t2noerr=4.0235   t3=4.0430   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=4.95  | x1.03 x1.07 x1.05 x1.04 -
normalize t1=4.5879   t2=4.3124   t2ff=4.0814   t2noerr=4.2673   t3=4.0104   t3hadd=-        t3acc4=-        t3rsqrt=4.0088   t3_gbps=3.99  | x1.06 x1.12 x1.08 x1.14 -
collision t1=9.5550   t2=7.1059   t2ff=6.9829   t2noerr=7.2374   t3=6.7836   t3hadd=-        t3acc4=-        t3rsqrt=6.7604   t3_gbps=4.72  | x1.34 x1.37 x1.32 x1.41 -

-- N=65536 --
add       t1=1.4260   t2=0.2923   t2ff=0.3266   t2noerr=0.3146   t3=0.2992   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=80.22 | x4.88 x4.37 x4.53 x4.77 -
scale     t1=2.1721   t2=0.2206   t2ff=0.2197   t2noerr=0.2263   t3=0.2228   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=71.83 | x9.85 x9.89 x9.60 x9.75 -
dot       t1=2.6521   t2=0.7721   t2ff=0.1768   t2noerr=0.8172   t3=0.3385   t3hadd=0.2707   t3acc4=0.1836   t3rsqrt=-        t3_gbps=47.27 | x3.43 x15.00 x3.25 x7.83 -
dot_elem  t1=1.8067   t2=0.2728   t2ff=0.2943   t2noerr=0.2742   t3=0.2703   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=74.00 | x6.62 x6.14 x6.59 x6.68 -
len       t1=1.8037   t2=1.5900   t2ff=0.3126   t2noerr=0.2657   t3=0.2294   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=87.17 | x1.13 x5.77 x6.79 x7.86 -
normalize t1=2.7905   t2=2.5382   t2ff=0.3021   t2noerr=2.5544   t3=0.3751   t3hadd=-        t3acc4=-        t3rsqrt=0.3691   t3_gbps=42.65 | x1.10 x9.24 x1.09 x7.44 -
collision t1=7.9121   t2=4.7021   t2ff=0.6592   t2noerr=4.5492   t3=0.6354   t3hadd=-        t3acc4=-        t3rsqrt=0.6898   t3_gbps=50.36 | x1.68 x12.00 x1.74 x12.45 -

-- N=8192 --
add       t1=1.5162   t2=0.2405   t2ff=0.2370   t2noerr=0.2392   t3=0.2319   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=103.48 | x6.30 x6.40 x6.34 x6.54 -
scale     t1=2.2139   t2=0.1645   t2ff=0.1672   t2noerr=0.1674   t3=0.1674   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=95.57 | x13.46 x13.24 x13.23 x13.23 -
dot       t1=2.6811   t2=0.7607   t2ff=0.1448   t2noerr=0.7646   t3=0.3213   t3hadd=0.2492   t3acc4=0.1428   t3rsqrt=-        t3_gbps=49.79 | x3.52 x18.52 x3.51 x8.34 -
dot_elem  t1=1.8714   t2=0.1979   t2ff=0.2136   t2noerr=0.1994   t3=0.2083   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=96.03 | x9.46 x8.76 x9.39 x8.98 -
len       t1=1.9551   t2=1.5597   t2ff=0.2750   t2noerr=0.2395   t3=0.2031   t3hadd=-        t3acc4=-        t3rsqrt=-        t3_gbps=98.46 | x1.25 x7.11 x8.16 x9.63 -
normalize t1=2.8541   t2=2.4967   t2ff=0.2901   t2noerr=2.4844   t3=0.3395   t3hadd=-        t3acc4=-        t3rsqrt=0.3082   t3_gbps=47.13 | x1.14 x9.84 x1.15 x8.41 -
collision t1=7.9584   t2=3.9682   t2ff=0.5882   t2noerr=3.8351   t3=0.5483   t3hadd=-        t3acc4=-        t3rsqrt=0.5552   t3_gbps=58.36 | x2.01 x13.53 x2.08 x14.51 -
```

单行原始日志的样子（含 5 次 trial 全部数值，节选自 [raw/bench_O3_native.txt](raw/bench_O3_native.txt)）：

```
# pinned cpu=2 running_on=2
# flashvec bench: cfg=O3_native
# n_full=1048576 n_mid=65536 n_small=8192 reps=100 reps_mid=400 reps_small=4000 trials=5 checks=8192
# avx2=yes fma=1
# compiled_features: avx2=1 fma=1 fast_math=0
# tsc_hz=1796550553 (1.797 GHz) cpu_khz_before=1696973
CASE cfg=O3_native tier=tier3 op=dot n=8192 reps=4000 med_ns=2632.4 min_ns=2523.0 ns_elem=0.3213 ns_elem_min=0.3080 gbps=49.79 chk=-151.689 trials_ns= 2756.7 2632.4 2637.9 2523.0 2581.2
CASE cfg=O3_native tier=tier3 op=normalize n=8192 reps=4000 med_ns=2780.9 min_ns=2758.0 ns_elem=0.3395 ns_elem_min=0.3367 gbps=47.13 chk=714.896 trials_ns= 2914.4 2852.2 2774.2 2758.0 2780.9
CASE cfg=O3_native tier=tier1 op=collision n=1048576 reps=100 med_ns=10019166.0 min_ns=9289511.0 ns_elem=9.5550 ns_elem_min=8.8592 gbps=3.35 chk=-2800.42 trials_ns= 9679309.1 9289511.0 10827762.7 10868032.8 10019166.0
CASE cfg=O3_native_ffastmath tier=tier2 op=len n=8192 reps=4000 med_ns=2252.9 min_ns=2162.1 ns_elem=0.2750 ns_elem_min=0.2639 gbps=72.73 chk=31466 trials_ns= 2162.1 2347.0 2252.9 2250.4 2470.7
CASE cfg=O3_native_ffastmath tier=tier2 op=normalize n=8192 reps=4000 med_ns=2376.8 min_ns=2366.4 ns_elem=0.2901 ns_elem_min=0.2889 gbps=55.15 chk=714.896 trials_ns= 2381.7 2370.2 2506.4 2366.4 2376.8
CASE cfg=O3_native_nomatherrno tier=tier2 op=len n=8192 reps=4000 med_ns=1962.2 min_ns=1871.1 ns_elem=0.2395 ns_elem_min=0.2284 gbps=83.50 chk=31466 trials_ns= 1871.1 2019.6 1962.2 1964.3 1953.7
```

### 3.3 内存带宽实测（解释 N=1<<20 为什么三档一样）

`bench/membw.c`（taskset -c 2，每项 5 次取最小）原始输出：

```
# membw: 每项 5 次取最小值，CPU 绑定由外部 taskset 决定
size         working_MB    triad_GBs  triad16_GBs    write_GBs     read_GBs   memcpy_GBs
N=1048576          12.0         7.01         9.35        14.05        16.68        17.52
N=4194304          48.0         6.10         8.13         4.76         8.57        10.03
N=8388608          96.0         5.99         7.98         3.99         8.13         9.14
```

- `triad` = `c[i]=a[i]+b[i]`（12 B/元素逻辑访存）；`triad16` 把写分配的隐性读也算进去（16 B/元素）。
- 对照 §3.1：N=1<<20 时档3 各运算实测 3.55–6.16 GB/s，正好挡在 triad 的 6–7 GB/s 附近。
  **这就是三档数字贴在一起的直接原因。**
- 这台机器绝对带宽本身偏低（12 MB 工作集 memcpy 17.5 GB/s，96 MB 时 9.1 GB/s），同机 loadavg 3–5。
  **原因未验证**（可能是移动 APU 功耗限制 + 共享内存控制器 + 同机负载），但三档条件相同，结论不受影响。

### 3.4 精度与逐位一致性（原始输出，[raw/ident_check.txt](raw/ident_check.txt)）

```
DIFF cfg=checks tier=tier1 op=add n=8192 max_abs=0.000e+00 max_rel=0.000e+00 n_diff=0 n_nan=0
IDENT cfg=checks tier=tier2 op=add n=8192 bits_equal=1 n_diff=0
IDENT cfg=checks tier=tier3 op=add n=8192 bits_equal=1 n_diff=0
DIFF cfg=checks tier=tier1 op=dot n=8192 scalar=-30.3377914 ref=-30.3377762 rel=5.030e-07
DIFF cfg=checks tier=tier2 op=dot n=8192 scalar=-30.3377914 ref=-30.3377762 rel=5.030e-07
DIFF cfg=checks tier=tier3 op=dot n=8192 scalar=-30.3378315 ref=-30.3377762 rel=1.823e-06
DIFF cfg=checks tier=tier3hadd op=dot n=8192 scalar=-30.3378277 ref=-30.3377762 rel=1.698e-06
DIFF cfg=checks tier=tier3_4acc op=dot n=8192 scalar=-30.3378124 ref=-30.3377762 rel=1.195e-06
DIFF cfg=checks tier=tier1 op=len n=8192 max_abs=1.192e-07 max_rel=1.192e-07 n_diff=690 n_nan=0
DIFF cfg=checks tier=tier3 op=len n=8192 max_abs=1.192e-07 max_rel=1.192e-07 n_diff=690 n_nan=0
IDENT cfg=checks tier=tier3 op=len n=8192 bits_equal=1 n_diff=0
DIFF cfg=checks tier=tier1 op=normalize n=8192 max_abs=1.788e-07 max_rel=1.913e-07 n_diff=572 n_nan=0
DIFF cfg=checks tier=tier2 op=normalize n=8192 max_abs=1.788e-07 max_rel=1.913e-07 n_diff=572 n_nan=0
IDENT cfg=checks tier=tier2 op=normalize n=8192 bits_equal=1 n_diff=0
DIFF cfg=checks tier=tier3 op=normalize n=8192 max_abs=1.788e-07 max_rel=1.913e-07 n_diff=572 n_nan=0
IDENT cfg=checks tier=tier3 op=normalize n=8192 bits_equal=1 n_diff=0
DIFF cfg=checks tier=tier3_rsqrt op=normalize n=8192 max_abs=2.384e-07 max_rel=2.624e-07 n_diff=3342 n_nan=0
IDENT cfg=checks tier=tier3_rsqrt op=normalize n=8192 bits_equal=0 n_diff=3146
DIFF cfg=checks tier=tier1 op=collision n=8192 max_abs=1.431e-06 max_rel=8.403e-07 n_diff=2105 n_nan=0
DIFF cfg=checks tier=tier2 op=collision n=8192 max_abs=1.431e-06 max_rel=8.131e-07 n_diff=5810 n_nan=0
IDENT cfg=checks tier=tier2 op=collision n=8192 bits_equal=0 n_diff=4948
DIFF cfg=checks tier=tier3 op=collision n=8192 max_abs=1.431e-06 max_rel=8.131e-07 n_diff=5810 n_nan=0
IDENT cfg=checks tier=tier3 op=collision n=8192 bits_equal=0 n_diff=4948
DIFF cfg=checks tier=tier3_rsqrt op=collision n=8192 max_abs=1.669e-06 max_rel=9.563e-07 n_diff=6615 n_nan=0
IDENT cfg=checks tier=tier3_rsqrt op=collision n=8192 bits_equal=0 n_diff=6118
```

怎么读：

- **add/scale/dot_elem/len/normalize：档2、档3 与档1 逐位一致**（`bits_equal=1`）。手写 AVX2（sqrt+div 路径）
  在拿到 x6–x13 的同时**没有付出任何精度代价**。
- tier1/tier2/tier3 与严格 IEEE 参考的差（len 1.19e-07、normalize 1.79e-07）**三者完全相同**，来源是 `-O3` 默认的
  `-ffp-contract=fast`：`a.x*a.x + a.y*a.y` 被收缩成 FMA，而参考实现用 `-ffp-contract=off` 编，所以差 1 ULP。
- **collision 暴露了一个真实语义差异**：档1（跨函数、逐算子）与档2/档3（同一表达式、可 FMA 收缩）在 8192 个元素里
  有 4948 个位模式不同（最大 1.43e-06 绝对 / 8.1e-07 相对 ≈ 7 ULP）。这是"函数式 API 天然阻止 FMA 收缩"的直接后果——
  对物理引擎既是**可复现性优势**也是**性能代价**。
- rsqrt+NR（档3_rsqrt）把 normalize 从 0.3395 降到 0.3082 ns/elem（快 9%），代价是 3146/8192 个元素位模式不同。

### 3.5 编译器的向量化裁决（GCC 原文，[raw/vect_verdicts.txt](raw/vect_verdicts.txt)）

```
### 编译选项：-O3 -march=native
   ../bench.cpp:236:28: note:  LOOP VECTORIZED    (t2_add)
   ../bench.cpp:244:28: note:  LOOP VECTORIZED    (t2_scale)
   ../bench.cpp:253:28: note:  LOOP VECTORIZED    (t2_dot)
   ../bench.cpp:261:28: note:  LOOP VECTORIZED    (t2_dot_elem)
   ../bench.cpp:270:28: missed:   not vectorized: control flow in loop.   (t2_len)
   ../bench.cpp:279:28: missed:   not vectorized: control flow in loop.   (t2_normalize)
   ../bench.cpp:289:28: missed:   not vectorized: control flow in loop.   (t2_collision)

### 编译选项：-O3 -march=native -fno-math-errno
   ../bench.cpp:270:28: note:  LOOP VECTORIZED      ← len 解锁
   ../bench.cpp:279:28: missed:   not vectorized: control flow in loop.   ← normalize 仍不行
   ../bench.cpp:289:28: missed:   not vectorized: control flow in loop.   ← collision 仍不行

### 编译选项：-O3 -march=native -ffast-math
   ../bench.cpp:270:28: note:  LOOP VECTORIZED
   ../bench.cpp:279:28: note:  LOOP VECTORIZED      ← normalize 解锁（代价见 §2.2）
   ../bench.cpp:289:28: note:  LOOP VECTORIZED      ← collision 解锁
```

- 生成方式：`g++ ... -fdump-tree-vect-details -c bench/bench.cpp`，取 `*.vect` 里对应函数的 LOOP VECTORIZED / missed 行。
- `t2_len/t2_normalize/t2_collision` 里都有条件表达式，但只有 len 被 `-fno-math-errno` 解锁：
  它循环里没有除法，唯一的控制流来自 `sqrtf` 的 errno 分支。normalize/collision 里还有一个**条件执行的 `1.0f/l`**，
  GCC 在严格 FP 语义下不能把它变成无条件向量除法（不允许推测执行可能除零的运算），所以必须 `-ffast-math`
  （含关闭 `-ftrapping-math`）。手写 AVX2 用掩码天然绕开——**这是"手写 SIMD"在本库里唯一站得住脚的理由**。

### 3.6 指令级证据（反汇编节选，完整版 [raw/evidence.txt](raw/evidence.txt)）

B1. 档1 `t1_add` 主循环（-O3 -march=native）：

```
    23d4:	vmovq  (%r14,%rbp,8),%xmm1        ← b[i] 直接进 xmm1
    23da:	vmovq  0x0(%r13,%rbp,8),%xmm0    ← a[i] 直接进 xmm0（无栈拷贝）
    23e1:	call   79e0 <fv_vec2_add>          ← 每个元素一次调用
    23e6:	vmovq  %xmm0,(%r15,%rbp,8)        ← 返回值 xmm0 直接写回
    23ec:	cmp    %ebx,%r12d
    23ef:	jne    23d0 <_ZL6t1_addj+0x30>
```

B2. 档1 `t1_dot` 主循环：累加器被挤到通用寄存器 `ebp`：

```
    2494:	vmovq  0x0(%r13,%rax,8),%xmm0
    249b:	vmovq  (%r14,%rax,8),%xmm1
    24a1:	call   7a50 <fv_vec2_dot>
    24a6:	vmovd  %ebp,%xmm3        ← 累加器从 ebp 搬到 xmm
    24aa:	vaddss %xmm0,%xmm3,%xmm2
    24ae:	vmovd  %xmm2,%ebp        ← 又搬回去（寄存器压力）
    24b2:	cmp    %ebx,%r12d
    24b5:	jne    2490 <_ZL6t1_dotj+0x20>
```

B3. 档2 `t2_dot`（无 fast-math）：乘法向量化了，但每 8 个元素做一次 8 路串行水平归约：

```
    21f0:	vmovaps (%rsi,%rax,1),%ymm4
    21f5:	vmovaps (%r8,%rax,1),%ymm5
    21fb:	vmulps (%rcx,%rax,1),%ymm4,%ymm1
    2200:	vfmadd231ps (%rdi,%rax,1),%ymm5,%ymm1     ← 乘法用了 FMA
    220a:	vshufps $0x55,%xmm1,%xmm1,%xmm3           ┐
    220f:	vshufps $0xff,%xmm1,%xmm1,%xmm2           │ 每轮 11 条标量/洗牌指令
    2214:	vaddss %xmm1,%xmm0,%xmm0                  │ 构成串行依赖链
    2218:	vaddss %xmm0,%xmm3,%xmm3                  │ → 0.7607 ns/elem
    221c:	vunpckhps %xmm1,%xmm1,%xmm0               │（档3 4累加器是 0.1428）
    2220:	vextractf128 $0x1,%ymm1,%xmm1              │
    2226:	vaddss %xmm3,%xmm0,%xmm0                  │
    222a:	vaddss %xmm2,%xmm0,%xmm0                  ┘
```

B4. 档2 `t2_len`（无 fast-math）：`call sqrtf@plt` 出现在循环里（errno 语义的代价）：

```
    27ed:	vmulss %xmm1,%xmm1,%xmm1
    27fc:	vsqrtss %xmm0,%xmm0,%xmm0
    2859:	call   11f0 <sqrtf@plt>
```

B5/B6. 档3 单累加器 vs 4 累加器：

```
t3_dot  (单累加器):  25e0: vfmadd231ps (%r8,%rax,4),%ymm2,%ymm0   ← 两条 FMA 依赖同一个 ymm0
                     25e6: vfmadd231ps (%rcx,%rax,4),%ymm3,%ymm0   ← 串行链 = 2×FMA 延迟 / 8 元素
t3_dot4 (4 累加器):  266d: vfmadd231ps (%r8,%rax,4),%ymm4,%ymm1   ← 4 条互不依赖的 FMA
                     2673: vfmadd231ps (%rcx,%rax,4),%ymm5,%ymm2
                     2679: vfmadd231ps 0x20(%r8,%rax,4),%ymm6,%ymm0
                     2680: vfmadd231ps 0x20(%rcx,%rax,4),%ymm7,%ymm3
```

B7. `fv_vec3` 的传参（探针 `__attribute__((noinline)) fv_vec3 probe3add(fv_vec3 a, fv_vec3 b){ return fv_vec3_add(b, a); }`）：

```
   0:	endbr64
   4:	movaps %xmm1,%xmm4
   7:	movq   %xmm2,-0x40(%rsp)
   d:	movdqa %xmm0,%xmm2
  11:	movaps %xmm3,%xmm1
  14:	movq   -0x40(%rsp),%xmm0
  1a:	movaps %xmm4,%xmm3
  1d:	jmp    22 <...>
```
→ `fv_vec3`（12 字节）走 `xmm0..xmm3`，每个结构体 2 个寄存器，第二个只用 4 字节。

---

## 4. 逐项分析（为什么是这个数）

### 4.1 add / scale：纯粹是「每元素一次 call」 vs SIMD
- N=8192：档1 1.516 / 2.214 ns/elem；档2 0.240 / 0.165；档3 0.232 / 0.167。
- 档1 循环体（B1）是 `vmovq/vmovq/call/vmovq` + 2 条循环簿记 ≈ 10 条指令/元素，其中 `call`+`ret` 至少 2–4 周期，
  且取值/回写无法跨元素流水；档2/档3 是 8 元素一条 `vaddps`/`vmulps`。
- scale 的加速比（x13）比 add（x6）大：档2/档3 的 scale 每元素只搬 16 B，add 要搬 24 B；在 L1 带宽 ~100 GB/s 的天花板下
  add 先撞上限（档3 实测 add 103.5 GB/s、scale 95.6 GB/s）。`dot_elem`（20 B/元素）介于两者之间。

### 4.2 dot：决定因素是**累加器依赖链**，不是 hadd
| 写法 | ns/elem (N=8192) | 相对档1 |
|---|---|---|
| 档1（每元素一次 `fv_vec2_dot` 调用） | 2.6811 | x1.00 |
| 档2 严格 FP（向量化但每轮 8 路串行归约，见 B3） | 0.7607 | x3.52 |
| 档2 + `-fno-math-errno` | 0.7646 | x3.51 |
| 档2 + `-ffast-math`（GCC 自动多累加器 + 末端归约） | 0.1448 | x18.52 |
| 档3 单累加器 + 末端 shuffle 归约 | 0.3213 | x8.34 |
| 档3 每 8 元素 hadd 归约（"naive SIMD"） | 0.2492 | x10.76 |
| 档3 **4 累加器** + 末端 shuffle 归约 | **0.1428** | **x18.78** |

- 结论一：**"每轮做一次水平归约"比"单累加器 + 末端归约"更快**，因为它把 2 条串行 FMA 链换成 1 条独立的标量加链。
  "hadd 慢"在本例中不成立——真正慢的是依赖链。
- 结论二：手写 4 累加器只比 GCC 自动向量化快 1.4%——**单靠 dot 不值得手写 SIMD**。
- 结论三：N=1<<20 时 dot 只快 x1.36–1.41（4 条流、16 B/元素，撞 6.16 GB/s 带宽墙）。

### 4.3 len：errno 分支的代价 = 6.5 倍
- 档2 严格 FP：1.5597（x1.25，几乎等于档1 的 1.9551——因为 `call sqrtf@plt` 和档1 的 `call fv_vec2_len` 是同一种代价）。
- 档2 + `-fno-math-errno`：0.2395（x8.16，便宜 6.5 倍，**结果逐位不变**）。
- 档3 手写：0.2031（x9.63），相对档2 只多 18%（`vsqrtps` 吞吐成了瓶颈）。

### 4.4 normalize：条件除法是分水岭
- 档2 严格 FP / +no-math-errno 都在 x1.1 上下（2.4967 / 2.4844）——除法没被向量化。
- 档2 + `-ffast-math` 0.2901（x9.84）；档3（`_mm256_div_ps`+掩码）0.3395（x8.41）；档3 rsqrt 0.3082（x9.26，代价 2.6e-07 相对误差）。
- 即：**精确除法版（档3）比 fast-math 版慢 17%，但逐位一致**；rsqrt 版又快 9%，代价可测。
  物理引擎的取舍就在这里，flashvec 应把两条路都给出来（如 `fv_vec2_normalize_soa` 与 `fv_vec2_normalize_soa_fast`）。

### 4.5 collision 微内核（法线归一化 + 速度反射）
- N=8192：档1 7.9584 → 档2 严格 3.9682（x2.01）→ 档2 +fast-math 0.5882（x13.53）→ 档3 0.5483（x14.51）。
  档2 严格只快 2 倍，因为它**没有被向量化**（§3.5 同一句 `not vectorized: control flow in loop.`），只是省掉了函数调用。
- N=1<<20：档1 9.5550 → 档3 6.7836（x1.41）；32 B/元素、实测 4.72 GB/s，是最接近内存上限的运算之一。
- **最值得手写 SIMD 的运算**：同时命中「条件除法」和「多算子融合」两个痛点，档3 相对档2 严格 FP 快 7.2 倍。

### 4.6 档1 为什么慢：调用开销 + 标量化，不是 ABI 拷贝
- 逐指令看（B1/B2）：参数与返回值都在 xmm 寄存器，**没有任何栈上的结构体拷贝**（`fv_vec2` 是 8 字节、SSE 类）。
- 真正开销：每元素一次 `call`（加被调用方 push/pop/ret）+ 无法跨元素 SIMD/流水；`dot` 还多出累加器在 `ebp` 与 xmm 之间来回搬的 3 条指令。
- 实测印证：N=8192 时档1 的 `add`（1.5162）与 `dot_elem`（1.8714）同量级，而两者库内计算量差别很大
  → **耗时由调用与循环开销决定，而不是数学运算**（1.5 ns ≈ 5.6 周期/元素，正是 call + 4 条访存/算术语的量级）。

### 4.7 `-ffast-math` 值不值
- 值：档2 在 len/normalize/collision 上从 x1.1–x2.0 跳到 x7.1–x13.5；dot 从 x3.5 跳到 x18.5。
- 不值：add/scale/dot_elem 上无额外收益（x6.40 vs x6.30 之类，在噪声内）。
- 代价：见 §2.2 的实测展开（重排、无 NaN/Inf 假设、除法变倒数、忽略 ±0、关闭 trapping/rounding）。
  对物理引擎意味着 NaN 不再被"零向量保护"挡住、除法精度可能降到 ~1e-7 相对（实测 rsqrt 路径 2.62e-07、3146/8192 元素位不同）。
- 建议：**库本身不要开 `-ffast-math`**；把"精确"与"近似"做成两个函数名，让调用者选。

---

## 5. perf stat：为什么没有 IPC / cache-miss / branch-miss 数字

```
$ cat /proc/sys/kernel/perf_event_paranoid
4
$ perf stat -e cycles,instructions ./bench/build/bench_perf tier2 4
Error:
Access to performance monitoring and observability operations is limited.
... (perf_event_paranoid setting is 4 ...)
$ echo $?
255
```

- `perf_event_paranoid = 4`（≥3 即拒绝一切非特权 `perf_event_open`），当前用户无 `CAP_PERFMON`，`sudo -n` 报"需要密码"。
  这是环境限制，不是基准没写。
- 补齐只需一条命令（`bench/Makefile` 已备好目标，跑的是缩短重复次数的 `bench_perf`）：

```
$ sudo sysctl -w kernel.perf_event_paranoid=1
$ make bench-perf
# 等价于：taskset -c 2 perf stat -e cycles,instructions,cache-references,cache-misses,\
#          branch-instructions,branch-misses build/bench_perf all 2
```

- 替代证据三类：① GCC 向量化裁决（§3.5）说明哪条循环被向量化/为什么没有；② 反汇编（§3.6）说明实际指令与依赖链；
  ③ 分级内存带宽（§3.3）+ 三档规模扫描（§3.1）等价于 cache-miss 曲线——加速比从 N=8192 的 x13 掉到 N=1<<20 的 x1.05，
  曲线本身说明瓶颈从计算转移到了内存。

---

## 6. 复现步骤与产物清单

```
$ cd flashvec
$ make test                 # 现有测试仍全绿（30 checks, 0 failures）
$ make bench                # 4 组 cfg × 3 规模 × 3 档 × 7 运算；原始日志 -> docs/research/raw/
$ make bench-matrix         # 紧凑矩阵（本文 §3.2）
$ make bench-table          # 宽表（每个 cfg 一列）
$ make bench-checks         # DIFF/IDENT 精度与逐位一致性（本文 §3.4）
$ make bench-vec            # GCC 向量化报告
$ make bench-asm            # 反汇编归档 -> bench/build/disasm_*.txt
$ make bench-perf           # 需先 sudo sysctl -w kernel.perf_event_paranoid=1
```

| 路径 | 内容 |
|---|---|
| `bench/bench.cpp` | 三档实现 + 计时夹具（屏障、校验和、逐位比较） |
| `bench/ref.cpp` / `bench/ref.h` | 严格 IEEE 参考实现（独立编译单元，固定 `-ffp-contract=off -fno-tree-vectorize`） |
| `bench/membw.c` | 分级内存带宽探针 |
| `bench/Makefile` | 构建运行规则（4 组 cfg × 4 个二进制 + checks/vec/asm/perf 目标） |
| `bench/matrix.awk` / `bench/table.awk` | 由原始日志生成两张表 |
| `docs/research/raw/bench_*.txt` | 每行一个 CASE/DIFF/IDENT 的原始日志（含 5 次 trial） |
| `docs/research/raw/{matrix,table}.{txt,md}` | 由日志生成的两张表 |
| `docs/research/raw/{membw,env,evidence,vect_verdicts,ident_check}.txt` | 带宽、环境、反汇编节选、向量化裁决、精度对照 |

本次未改动 `include/`、`src/`、`tests/` 的任何既有内容；根 `Makefile` 只在末尾**追加**了 `bench*` 转发目标
（前 39 行与原文件逐字节相同，`make test` 仍全绿）。

---

## 7. 参考资料（链接均在 2026-10-17 用 curl -L -w '%{http_code}' 验证返回 200）

1. GCC 优化选项文档（`-ffast-math` 展开、`-fno-math-errno`、`-ffp-contract`、`-ftree-vectorize`）：<https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>
2. GCC 浮点语义 wiki（各组开关含义与风险）：<https://gcc.gnu.org/wiki/FloatingPointMath>
3. GCC 向量扩展与 `-fopt-info`：<https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html>
4. cppreference：`math_errhandling` / errno 语义（解释 `sqrtf` 为何带 errno 分支）：<https://en.cppreference.com/w/c/numeric/math/math_errhandling>
5. Agner Fog《Instruction tables》（含 Zen 2 指令延迟/吞吐）：<https://www.agner.org/optimize/instruction_tables.pdf>
6. Agner Fog《Optimizing software in C++》（AoS/SoA 与调用开销）：<https://www.agner.org/optimize/optimizing_cpp.pdf>
7. uops.info 指令数据库（含 Zen 2 列）：<https://uops.info/table.html> ；`VFMADD213PS` 条目：<https://uops.info/html-instr/VFMADD213PS_YMM_YMM_YMM.html>
8. Intel Intrinsics Guide（`_mm256_hadd_ps/_mm256_div_ps/_mm256_sqrt_ps/_mm256_cmp_ps` 语义）：<https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html>
9. System V AMD64 ABI（结构体分类与寄存器传参，`fv_vec2/fv_vec3` 传参依据）：<https://gitlab.com/x86-psABIs/x86-64-ABI>
10. 参考实现对照（**本报告未从中引用任何数字**）：Box2D <https://github.com/erincatto/box2d> 、JoltPhysics <https://github.com/jrouwe/JoltPhysics>

未能验证的来源（**不作为依据**）：AMD 官方 Zen 2 优化指南与本机型产品页在当前网络下返回 000/404，
因此本报告关于 Zen 2 微架构的论断只用**本机实测**（§3.3/§3.5/§3.6）支撑，不引用 AMD 文档数字。
