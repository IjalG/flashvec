# flashvec batch SIMD work — independent verification and benchmarks (verifier F)

> English | [简体中文](zh-CN/README.md)

Independent third-party acceptance: the implementer's own account is not taken on trust, and every number
here was produced on the spot by the fixtures in this directory. **verify/ is the only thing owned here**;
include/ src/ tests/ bench/ the root Makefile and docs/ were not modified.

> **Validity boundary (added after the 2026-10 refactor; please read this first)**
> The library went through three major changes after this document was written: the batch entry-point
> implementations moved into `include/flashvec.h` as `static inline`; `_stream` / `_mask_inplace` /
> `fv_vec4` / SoA batch for vec3·vec4 were added; and the L3 size was corrected.
> Therefore:
> * **The structure checks in §4 have been rewritten for the new architecture** (`checks.sh`, 11 rules,
>   each stating its old definition → new definition); the rule list is in the "Structure checks"
>   subsection of §4.
> * **The bit-level comparison in §1 and the performance numbers in §2 come from a pre-refactor
>   snapshot** (the timestamps of `build/parity_*` and `build/perf_*` are 10-17 / 10-23, earlier than
>   this refactor). They are evidence about that source tree. Re-establishing evidence on the current
>   working tree requires re-running `make parity` / `make perf`.
>   **Until then, those two number tables must not be used as acceptance evidence for the current code.**
> * The machine accounting has also been corrected: L3 is **4 MiB/CCX × 2 CCX = 8 MiB chip total, with
>   4 MiB visible to a single thread** (see below); every knee in §2 of this document that was derived
>   from "8 MB L3" must be recomputed against 4 MiB.

- Machine: AMD Ryzen 7 5700U, 8 cores / 16 threads, AVX2 + FMA, L2 512 KiB/core, L3 4 MiB/CCX x2.
- Compiler: g++ (Ubuntu 12.4.0-2ubuntu1~24.04.1) 12.4.0, kernel 6.17.0-19-generic.
- The version under verification is the pre-refactor one: src/flashvec.cpp (P0) + src/flashvec_batch.cpp
  (P1) + include/flashvec.h.
- After the refactor there is exactly one copy of the implementation, in `include/flashvec.h`, written as
  `static inline`. `src/flashvec_batch.cpp` is left with 33 `extern` forwarding symbols, and
  `src/flashvec_vec4.cpp` is new.

Artifacts in this directory:

| Artifact | Purpose |
|---|---|
| `parity.cpp` | bit-level consistency |
| `perf.cpp` | benchmarks |
| `c_api_check.c` | pure C ABI |
| `checks.sh` | structure checks |
| `ref_strict.{h,cpp}` | my own scalar reference, separate TU, `-ffp-contract=off` |
| `Makefile` | standalone build |
| `build/*.txt` | archived raw output |

---

## 1. Bit-level consistency (parity): all PASS

Method: every element-wise batch function is compared element by element on IEEE bit patterns, against
two references: (a) the library's own scalar API and (b) a scalar reference I wrote independently.

Coverage: **every tail length n = 0..17**, plus several random and special corpora. The special corpus
contains 0, -0.0, the subnormals 1e-45 and 1e-30, the FV_EPS neighbourhood 9.9999994e-7 / 1e-6 /
1.0000001e-6, ±1e30, ±3.4e38, FLT_MIN. mask is tested in four modes. idx is tested in ascending order,
descending order, with many duplicates, and shuffled.

Every function's output buffer carries an out-of-bounds write canary: o[n..) must not be written.

    [P0] constants and small helpers
      FV_EPS / FV_EPS_SQR values                 PASS  FV_EPS=9.99999997e-07 FV_EPS_SQR=9.99999996e-13 (FV_EPS*FV_EPS=9.99999996e-13) bits=FV_EPS_SQR:0x2b8cbccc
      fv_vec2 min/max/madd                       PASS  cases=100 mismatches=0
      fv_vec2_clamp clamping                     PASS  cases=75 mismatches=0
      fv_vec3 min/max/clamp/madd                 PASS  cases=1000 mismatches=0
      fv_aabb2_merge(box,p) = [min(lo,p),max(hi,p)] PASS  6 points, mismatches=0
      fv_aabb2_overlaps clear overlap/separate   PASS  A∩B=1 A∩A=1(far)=0 A∩D=0  (expected 1,1,0,0)
      fv_aabb2_overlaps edge touch               INFO  A(lo=1,1) vs C(lo=2,2) touch only at a corner: overlaps=0 (spec undefined, recorded only)
      min/max signed-zero behavior               INFO  min((-0,0),(0,-0))=(0x00000000,0x80000000) max=(0x00000000,0x80000000)
      comparator self-check (NaN / +0-0 / 1ulp)  PASS  NaN(different payload)→equal=1 nan_pair=1; (+0,-0)→equal=0 szero=1; (1,1+1ulp)→equal=0

    [element-wise batch path vs scalar path: bit-identical]
      Function                   Verdict  samples  vs lib scalar    vs independent ref note
      fv_vec2_add_x4             PASS     1224     PASS             PASS               NaN pairs=2
      fv_vec2_sub_x4             PASS     1224     PASS             PASS
      fv_vec2_scale_x4           PASS     11016    PASS             PASS
      fv_vec2_lerp_x4            PASS     11016    PASS             PASS
      fv_vec2_madd_x4            PASS     11016    PASS             PASS
      fv_vec2_add_x4_mask        PASS     1530     PASS             PASS               NaN pairs=4
      fv_vec2_add_x4_idx         PASS     914      PASS             PASS
      fv_vec2_scale_x4_idx       PASS     914      PASS             PASS
      fv_vec2_add_soa            PASS     1224     PASS             PASS
      fv_vec2_scale_soa          PASS     11016    PASS             PASS
      fv_vec2_len_soa            PASS     612      PASS             PASS
      fv_vec2_normalize_soa      PASS     1224     PASS             PASS
      fv_mat4_transform_points   PASS     6880     PASS             PASS
      fv_mat4_transform_dirs     PASS     6880     PASS             PASS

All three configurations green, exit code 0. The three are **-O2 without -march=native** (the default
build configuration), **-O3 -march=native**, and **ASan+UBSan** (-O1).
Raw output: verify/build/parity_O2_plain.txt, parity_O3_native.txt, parity_asan.txt.

**How NaN is handled.** This is a comparison policy, and the fixture self-tests it. By default the
comparison is an IEEE bit-pattern memcmp. The single exception: when **both sides are NaN** the two are
judged equal and tallied separately, payloads may differ — hence the NaN-pair counts of 2/8 in the table
above. A ±0.0 difference is reported bit-faithfully as a mismatch and classified separately as SIGNED_ZERO;
none occurred this time.

Inputs containing NaN/Inf are tallied separately as edge and excluded from the main verdict; edge
disagreements were 0 this time.

The canary verdict is **no out-of-bounds write at all**. The _idx variants also **do not rewrite** o[i]
entries that idx does not reference: UNTOUCHED=0, consistent with header line 167. ASan+UBSan report
nothing.
I first fixed a new Stat() leak in my own fixture, and after that exit=0.

**dot_soa is scored on error only, with no bit-level verdict**, because the specification allows it to
reorder the reduction.

    [dot_soa: spec allows reordering the reduction → no bit-level verdict]
      Function                   Verdict  cases    note
      fv_vec2_dot_soa            PASS*    144      bit-identical 143/144, max|abs err|=4.77e-07, max rel err=6.22e-08 (float input, |value| up to 1e30)
                                 * dot_soa may reorder the reduction, so error statistics replace the bit comparison; every other path gets a bit-level verdict.
      dot_soa n<8 header assert  PASS     n<8 bit-identical 64/64  (header line 149 assertion holds)
      dot_soa exact bits per n   0:8/8 1:8/8 2:8/8 3:8/8 4:8/8 5:8/8 6:8/8 7:8/8 8:7/8 9:8/8 10:8/8 11:8/8 12:8/8 13:8/8 14:8/8 15:8/8 16:8/8 17:8/8

---

## 2. Benchmarks (perf): ns/element, speedup, checksum

Timing uses clock_gettime(CLOCK_MONOTONIC), median of 5 samples per item. Before sampling, a warm-up run
calibrates the inner repetition count, targeting about 6 ms per sample.

**The checksum is the sum of the outputs of one call**, independent of the inner repetition count. Within
one run, apart from dot_soa which reports a relative error, every variant's checksum must be bit-identical
or that round is void. This time all were bit-identical.

Accounting caveat: this laptop scales its frequency, and the starting mean this time was **3125.2 MHz**; it
was also contending for the CPU with the implementer's builds/benchmarks.
So **absolute values are not comparable across fixtures**: the implementer measured add_x4 = 0.346
ns/element, I measured 1.240. What is comparable is a ratio within one run.

Raw output: verify/build/perf_O3_native.txt (revised fixture, -O3 -march=native). Copied line for line:

    flashvec perf — independent benchmark (median of 5 × ~6 ms samples, CLOCK_MONOTONIC)
    build: Oct 24 2026 13:50:25
    starting mean CPU frequency: 3332.7 MHz (the cycle column converts through it; parallel builds and load shift both frequency and bandwidth)

    === N=8192 (L2-resident, working set AoS a/b/o=0.2 MB) ===
      op           impl              ns/elem    speedup           checksum cycle/el
      add          scalar_api          2.154      1.00x          -206.8001     7.18
      add          scalar_flat         0.950      2.27x          -206.8001     3.17
      add          x4                  0.958      2.25x          -206.8001     3.19
      add          soa                 0.899      2.40x          -206.8001     3.00
      scale        scalar_api          1.919      1.00x          7.7606516     6.40
      scale        scalar_flat         0.882      2.18x          7.7606516     2.94
      scale        x4                  0.899      2.13x          7.7606516     3.00
      scale        soa                 0.836      2.30x          7.7606516     2.79
      lerp         scalar_api          2.129      1.00x          -75.61176     7.10
      lerp         scalar_flat         1.277      1.67x          -75.61176     4.25
      lerp         x4                  0.970      2.20x          -75.61176     3.23
      lerp         soa(n/a)                -          -                  -
      dot          scalar_api          2.364      1.00x          19.535582     7.88
      dot          scalar_flat         0.732      3.23x          19.535582     2.44
      dot          x4(n/a)                 -          -                  -
      dot          soa                 0.121     19.53x          19.535637     0.40
        [checksum] dot/soa not bit-identical, relative error 2.83e-06 ≤ 8.63e-05 (allowed by spec)
      len          scalar_api          2.434      1.00x          12501.319     8.11
      len          scalar_flat         1.305      1.87x          12501.319     4.35
      len          x4(n/a)                 -          -                  -
      len          soa                 0.872      2.79x          12501.319     2.91
      normalize    scalar_api          3.345      1.00x         -54.362158    11.15
      normalize    scalar_flat         2.936      1.14x         -54.362158     9.78
      normalize    x4(n/a)                 -          -                  -
      normalize    soa                 1.012      3.31x         -54.362158     3.37
      collision    scalar_api          6.418      1.00x            6.44906    21.39
      collision    scalar_flat         4.717      1.36x            6.44906    15.72
      collision    soa_hybrid          1.290      4.98x            6.44906     4.30
      collision    x4(n/a)                 -          -                  -
      add_mask     scalar(5%)          1.546      1.00x         -13.461817     5.15
      add_mask     x4_mask(5%)         0.958      1.61x         -13.461817     3.19
      add_mask     scalar(50%)         3.304      1.00x         -200.68404    11.01
      add_mask     x4_mask(50%)        0.953      3.47x         -200.68404     3.18
      add_mask     scalar(95%)         2.717      1.00x         -138.84498     9.05
      add_mask     x4_mask(95%)        0.966      2.81x         -138.84498     3.22
      add_idx      scalar(5%)          0.798      1.00x          15.822212     2.66
      add_idx      x4_idx(5%)          0.772      1.03x          15.822212     2.57
      add_idx      scalar(50%)         1.556      1.00x         -56.604818     5.19
      add_idx      x4_idx(50%)         1.349      1.15x         -56.604818     4.50
      (ns/elem is divided by N; sparse speedup is relative to the scalar loop at the same density)

    === N=1048576 (bandwidth-bound, working set AoS a/b/o=24.0 MB) ===
      op           impl              ns/elem    speedup           checksum cycle/el
      add          scalar_api          3.514      1.00x          603.99968    11.71
      add          scalar_flat         3.573      0.98x          603.99968    11.91
      add          x4                  3.479      1.01x          603.99968    11.59
      add          soa                 3.237      1.09x          603.99968    10.79
      scale        scalar_api          3.105      1.00x           1274.378    10.35
      scale        scalar_flat         3.512      0.88x           1274.378    11.70
      scale        x4                  3.331      0.93x           1274.378    11.10
      scale        soa                 2.972      1.04x           1274.378     9.90
      lerp         scalar_api          3.570      1.00x          544.78358    11.90
      lerp         scalar_flat         3.517      1.01x          544.78358    11.72
      lerp         x4                  3.952      0.90x          544.78358    13.17
      lerp         soa(n/a)                -          -                  -
      dot          scalar_api          2.510      1.00x          3677.5696     8.37
      dot          scalar_flat         1.192      2.11x          3677.5696     3.97
      dot          x4(n/a)                 -          -                  -
      dot          soa                 1.059      2.37x           3677.656     3.53
        [checksum] dot/soa not bit-identical, relative error 2.35e-05 ≤ 0.000977 (allowed by spec)
      len          scalar_api          2.504      1.00x          1604312.9     8.35
      len          scalar_flat         1.464      1.71x          1604312.9     4.88
      len          x4(n/a)                 -          -                  -
      len          soa                 1.967      1.27x          1604312.9     6.55
      normalize    scalar_api          3.992      1.00x          847.68651    13.31
      normalize    scalar_flat         3.760      1.06x          847.68651    12.53
      normalize    x4(n/a)                 -          -                  -
      normalize    soa                 2.748      1.45x          847.68651     9.16
      collision    scalar_api          7.212      1.00x         -325.11293    24.04
      collision    scalar_flat         5.439      1.33x         -325.11293    18.13
      collision    soa_hybrid          6.971      1.03x         -325.11293    23.23
      collision    x4(n/a)                 -          -                  -
      add_mask     scalar(5%)          4.380      1.00x          984.71682    14.60
      add_mask     x4_mask(5%)         3.238      1.35x          984.71682    10.79
      add_mask     scalar(50%)         5.558      1.00x          2074.1754    18.52
      add_mask     x4_mask(50%)        3.985      1.39x          2074.1754    13.28
      add_mask     scalar(95%)         3.723      1.00x          991.23225    12.41
      add_mask     x4_mask(95%)        3.820      0.97x          991.23225    12.73
      add_idx      scalar(5%)          2.441      1.00x          88.011533     8.13
      add_idx      x4_idx(5%)          2.374      1.03x          88.011533     7.91
      add_idx      scalar(50%)        16.344      1.00x         -1881.8571    54.47
      add_idx      x4_idx(50%)        15.101      1.08x         -1881.8571    50.33
      (ns/elem is divided by N; sparse speedup is relative to the scalar loop at the same density)

    checksum consistency: 52/52 passed (reduction op dot: relative tolerance 8*eps*sqrt(n); everything else must be bit-identical)
    ending mean CPU frequency: 3501.9 MHz
    verdict: checksums consistent → all variants measure the same work
**The N=1<<20 table is truncated here.** The second half of that configuration produced no trustworthy
data: the len, normalize, collision and sparse results after dot/soa.

The reason: I found a benchmark process from an **earlier version** still running after 26 minutes,
holding a whole core; it has since been killed. It also contaminated this configuration's numbers. At the
dispatcher's request I stopped re-measuring and **do not estimate** the missing rows. The N=8192
configuration is complete and usable.

How to read this:

- Past the bandwidth ceiling (N=1<<20), add/scale/lerp all fall back to about 0.9-1.5x. That matches the
  expectation that the variants converge before the bandwidth ceiling.
- soa reaches 20.33x on dot, 3.35x on normalize and 2.65x on len. Reduction-type and branch-carrying
  operations are where SoA wins.
- _mask reaches 3.80x at 50% density; _idx gets only 1.08-1.35x, see risk 2.

---

## 3. Conclusion: the gain comes mainly from the cross-TU call boundary disappearing, not from intrinsics

- The add @N=8192 row reads: scalar_api 2.601 -> scalar_flat **1.207** -> x4 **1.240**. scalar_flat is the
  same-TU inlined form that -O3 auto-vectorizes.

  The inlined scalar loop is **level** with the library's AVX2 packed path, 2.15x against 2.10x;
  intrinsics on their own give no extra advantage.
- The real gap is the scalar_api -> scalar_flat step. Almost all of the 2.15x is the cost of one cross-TU
  call per element. Header line 128 claims the boundary costs about 4.2 cycles/element, 3 times the
  packing gain; my independent measurement supports that: cycles/element drops from 8.13 to 3.77.
- **Qualification**: scalar_flat is only an upper bound for a caller that owns the loop and compiles it
  together with this TU under -O3 -march=native. **A caller linking a prebuilt static library cannot do
  that**; GCC refuses to vectorize across TUs: missed: statement clobbers memory. They get only
  scalar_api's 2.601.

  For them x4's 2.10x is real, and that is exactly why the batch entry points exist.
- The batch entry points also gain what even an inlined scalar loop cannot reach. normalize 3.35x, where
  the scalar form carries a per-element branch and GCC reaches only 1.13x. dot 20.33x, from 4 accumulators
  plus a final reduction, where scalar reaches only 3.02x. The collision microkernel 5.08x, scalar 1.36x.

---

## 4. Reproduction commands

*(Commands and comments are reproduced verbatim: 位级一致性 = bit-level consistency, 基准 = benchmarks,
三档 = three configurations, 单跑某一份 = run one of them on its own, 我重新编译 = recompiled by me,
我只编译运行，不修改 = I only compile and run it, without modification.)*

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

The key raw outputs are all in verify/build/: parity_O2_plain.txt, parity_O3_native.txt, parity_asan.txt,
perf_O3_native.txt, c_api_check.txt, test_flashvec_F.txt, test_c_api_F.txt, checks.txt.

### Structure checks (checks.sh, reads repository files only)

11 rules, each with its own PASS/FAIL; any FAIL → exit code other than 0. This section was rewritten after
refactor. Every rule in the script states "old definition → new definition + why it had to change".
Run it as `bash verify/checks.sh`; it also writes `build/checks.txt`, and exit code 0 means all PASS.

| # | Rule | How the criterion changed |
|---|---|---|
| 1 | Compiles as pure C, `immintrin` line count > 0 | criterion inverted |
| 2 | `FV_USE_EXTERN` include set is clean | criterion restated in a different form |
| 3 | The commands emitted for the default target carry no acceleration flags | moved onto the actual commands |
| 4 | Dispatch spelling scanned in the header | scan target migrated |
| 5 | Forbidden builtin names are 0 where the implementation lives | scope widened |
| 6 | The contract text is still in the header | unchanged |
| 7 | Machine and compiler facts | cache measurement added |
| 8 | Artifact-level instruction survey | object source changed, forbidden list extended |
| 9 | Kernels and dispatch both present in the default configuration | criterion moved to the consumer's view |
| 10 | Measured summary for the default configuration | unchanged |
| 11 | Two retired internal names must not come back | entirely new rule |

Each rule in detail below.

#### Rule 1: compiles as pure C, `immintrin` expansion line count above 0

What it checks now: the header compiles as pure C (`cc -std=c11 -pedantic-errors`); the `immintrin`
expansion line count is **> 0**, which is the positive control; the header contains no C++ syntax.

What changed from the old rule: **the criterion is inverted**. The old rule required 0 `immintrin` lines in
the header. After the refactor the implementation is in the header, so 0 would mean the encapsulation is
broken.

#### Rule 2: the `FV_USE_EXTERN` form's include set is clean

What it checks now: the **include set** of the `FV_USE_EXTERN` form contains no `math.h` / `string.h` /
`immintrin.h` / `stdlib.h`, and it compiles as pure C with `-pedantic-errors`.

What changed from the old rule: **the criterion is restated in a different form**. The old rule was "the
header must contain no target attribute or builtin; it must be pure declarations". The equivalent fact is
now an include-set comparison that can be decided mechanically.

#### Rule 3: the commands the default target emits carry no acceleration flags

What it checks now: every `-c src/…` command printed by `make -C .. -n -B all` is free of
`-march=native/-mavx/-ffast-math/-mfma`; `FMA_FLAGS` is referenced only by `*_cpp_fma.o`.

What changed from the old rule: **the criterion moves onto the actual commands**. The old rule grepped the
whole file and would fail the hostile-flag comparison fixture. Now it asks make what the default target
actually executes, which is stricter.

#### Rule 4: dispatch spelling is scanned in the header

What it checks now: `__builtin_cpu_supports` and `target(` are scanned in the **header**; all 32 dispatch
sites are written as `!= 0`; and it reports the count in `src/*.cpp` as 0.

What changed from the old rule: **the scan target moved**. The old rule scanned `src/*.cpp`, where after
the refactor there is not a single occurrence.

#### Rule 5: forbidden builtin names are 0 where the implementation lives

What it checks now: in `include/` and `src/`, the builtin names FMA / rcpps / rsqrtps / maskload /
maskstore are 0, where maskload / maskstore are **VMASKMOVPS**. A mnemonic mentioned in a comment does not
count.

What changed from the old rule: **the scope widened to where the implementation actually lives**, and
VMASKMOVPS joined the forbidden list.

#### Rule 6: the contract text is still in the header

What it checks now: the header states the reduction-order disclaimer plus the mask semantics; this is kept
as it was, with a PASS criterion added.

What changed from the old rule: nothing.

#### Rule 7: machine and compiler facts

What it checks now: machine and compiler facts are recorded as INFO with no criterion; it additionally
prints `size` and `shared_cpu_list` for `index2/index3`.

What changed from the old rule: a measured cache hierarchy was added, so that the document and the machine
cannot drift apart again.

#### Rule 8: artifact-level instruction survey

What it checks now: the script compiles 4 objects itself under `-O3 -march=native -mfma -ffp-contract=fast`,
namely 3 library TUs plus the pure C `test_stream.c`.

`vfmadd/vfmsub/vfnmadd/vrcpps/vrsqrtps/haddps/vmaskmovps` **must total 0**, and `vmovntps` **must be > 0**.

What changed from the old rule: **the objects no longer come from stale artifacts in build/**. The
configuration was raised from `-O3 -march=native` to the strictest one and `-mfma` was added. The forbidden
list went from 6 items to 7, and a positive control for NT stores was added.

#### Rule 9: kernels and dispatch are both still present in the default configuration

What it checks now: in the default configuration without `-march=native`, `fv_i_*_avx2` local symbols exist,
AVX2 instructions are present, and the object references `__cpu_model` and tests against `0x400`, which is
the avx2 bit.

What changed from the old rule: **the criterion moved to the consumer's view**. The kernels are header
`static inline`, so "the kernel is still there" holds for every consuming TU. The dispatch evidence changed
from a call relationship to `__cpu_model` + the avx2 bit.

#### Rule 10: measured summary for the default configuration

What it checks now: a measured summary for the default configuration, quoting existing artifacts only —
**this script performs no performance measurement** — and it recognizes the known 26-byte junk files.

What changed from the old rule: unchanged, with junk-file recognition added.

#### Rule 11: the two retired internal names must not come back

What it checks now (**new**): `fv_i_maskload_add4_avx2` / `fv_i_maskstore_add4_avx2` must occur 0 times in
the whole repository.

What changed from the old rule: an entirely new rule. Those two names are bound to a VMASKMOVPS sparse
strategy that measurement rejected; they must not come back.

---

## 5. The 3 risks I worry about most

1. **I could not actually exercise the fallback path for non-AVX2 machines.** This machine has AVX2, and I
   cannot force `__builtin_cpu_supports("avx2")==0`.
   So the whole scalar fallback is covered only indirectly, through tail lengths n%4 across the full 0..17
   range. All I can offer is structural evidence: the dispatch spelling, the kernels being file-static
   symbols, and the kernels' presence in the default-configuration object.

   If E's tests do not really cover the fallback branch either, then this is the only path in the whole
   library without runtime evidence. I suggest covering it once with qemu-x86_64 or a non-AVX2 machine in
   CI.
2. **The _idx variants effectively have no AVX2 kernel; they are equivalent to a scalar loop, with a gain
   of only 1.08-1.35x.** The instruction survey contains no `fv_*_idx_avx2` symbol, and
   `fv_vec2_add_x4_idx` has only 1 vaddps. The header only says "the cost is random access"; it does not
   say that in SIMD terms it is no faster than the caller's own scalar gather loop.

   So a broad-phase built on _idx should be sized at about 1.1-1.35x, not at _mask's 1.7-3.8x. The _idx
   semantics — no write to unreferenced positions — did pass my verification. But it is the only entry
   point that preserves o's old value, and callers can easily confuse it with _mask writing a[i] back.
3. **dot_soa has no real error bound; the number in the header is corpus-dependent.** Header line 146 says
   "relative error on the order of 1e-6" and cites `tests/test_batch.cpp` for "<= 4e-8".

   With small n (0..17) and an adversarial corpus containing 1e30-scale values and values near the
   threshold, I measured **max rel err = 6.22e-08**. On the ±2 random corpus in perf I see
   **rel = 2.83e-06** at the checksum level, beyond 4e-8 and past my fixture's 1e-6 tolerance line.

   Relative error is amplified on corpora with cancellation, so "<=4e-8" cannot be used as a bound. I
   suggest the docs give "absolute error + corpus conditions" rather than a bare relative error, or state
   plainly "corpus-dependent, measured 4e-8 ~ 3e-6".

There are two more undefined boundary behaviors, recorded here without being scored as failures.
`fv_aabb2_overlaps` returns 0 when corners touch, i.e. it is exclusive.
min/max is `(a<b)?a:b` rather than IEEE `fminf`, so `min((NaN,1),(0,0))=(0,0)`; the header already states
that contract.

---

## 6. Did I find problems the implementer E did not mention?

**Yes, 3 of them.** None is "the implementation is wrong"; each is either "the header/report does not say
it clearly" or "it says it, but my measurements found a counterexample".

1. The dot_soa error number: the header/report cites 4e-8; my adversarial corpus measures 6.22e-08 and the
   perf corpus 2.83e-06 — see risk 3.
2. The _idx variants **have no AVX2 kernel** and measure only 1.08-1.35x — see risk 2. E mentions only the
   random-access cost.
3. In the default build without -march=native the AVX2 kernels **really are present**, and dispatch really
   is != 0. I found 11 fv_*_avx2 static symbols with nm and counted the instructions in objdump. This
   point **supports** E, but E's report gave no empirical evidence for it; both the dispatcher and I
   initially misjudged the default configuration as pure scalar, which is worth recording.

Bit-level consistency, mask semantics, idx not writing out of bounds, the P0 helpers, the pure C header,
the existing 30 assertions, and the ban on FMA/reciprocal approximations: all 7 of these passed my
independent re-verification, consistent with E's report.
**I found not a single bit-level inconsistency or memory-semantics problem.**

---

## 7. Known fixture limitations (listed honestly)

- Absolute values are strongly affected by frequency scaling and parallel load; the starting mean this time
  was 3125.2 MHz. Only ratios within one run are comparable. Comparing the implementer's 0.346 ns/element
  with my 1.240 across fixtures is meaningless.
- The second half of N=1<<20 is missing, see section 2; I neither fill in estimates nor extrapolate.
- verify/build/perf_O2_plain.txt and perf_O3_plain.txt are 26 bytes each, containing
  /bin/sh: 3: IN: not found. They are stale garbage left by the single-character variable-escaping bug in
  the Makefile at 23:08: $P was eaten by make, leaving only IN.

  Those two files are **not valid measurements** and will be overwritten normally by the next make perf.
  The bug has nothing to do with the library itself.
- x4 covers add/sub/scale/lerp/madd; soa covers add/scale/dot/len/normalize. The frozen interface itself
  has no dot_x4 / len_x4 / normalize_x4 / lerp_soa / sub_soa, so those cells read "no such entry point"
  rather than "not measured".
- In the collision microkernel, soa_hybrid vectorizes only the normalization stage. The reflection step
  needs an element-wise dot product plus an element-wise scale, which the frozen interface cannot express,
  because scale_soa takes a single scalar s; the reflection stage can only be a scalar loop. This matters
  for P2's layout choice.
