# flashvec Optimization Report

> English | [简体中文](zh-CN/optimization-report.md)

> Audience: the author of flashvec.
> Goal: make "where to optimize, in what order, what each step is worth, and what it costs" clear, and
> explicitly mark which parts are measured and which are inferred.
> How this was produced: four parallel research tracks (SIMD hardware / ideas borrowed from other fields /
> data layout / three-tier measurement) plus the dispatcher's own adjudication experiments and synthesis.
> Source material lives in `docs/research/`: `simd-hardware.md`, `algorithm-tricks.md`,
> `layout-memory.md`, `measured.md`.
> `raw/` holds the raw logs.
> Machine: AMD Ryzen 7 5700U (Zen 2, 8C/16T), AVX2+FMA, no AVX-512, g++ 12.4, Linux 6.17.
> Cache: L1d 32 KiB/core, L2 512 KiB/core, **L3 4 MiB/CCX × 2 CCX = 8 MiB per chip, 4 MiB visible to one
> thread**, 64 B lines.

---

> **Post-P0 note (2026-10-31):** this report's scalar API baseline is the pre-P0 `.cpp` form (a
> cross-TU call). P0 moved the scalar bodies into the header as `static inline` under the same
> `fp-contract=off` pragma, so the default scalar form now inlines and auto-vectorizes; the
> `FV_USE_EXTERN` form retains the old baseline. The historical ratios here remain valid for that
> extern baseline, but batch-vs-default-scalar ratios are smaller.


## 0. One-Page Summary

**Five conclusions** (ordered by "how much do I get back if I change this"):

| # | Conclusion | Strength of evidence | Size of change |
|---|---|---|---|
| 1 | Scale decides everything | Measured on this machine (three tiers × three sizes + my own independent rerun) | 0 (awareness) |
| 2 | The call boundary is a first-order bottleneck | Three independent lines of evidence + my own recomputation | Medium (batch entry points) |
| 3 | `-fno-math-errno` is the largest single win at zero semantic cost | Measured on this machine (A / D both directions + my recomputation) | Small (one compiler flag) |
| 4 | What blocks auto-vectorization is the "conditionally executed division", not sqrt | Measured on this machine (D, including GCC's verbatim verdict and bits_equal=1) | Medium (a few dozen lines of intrinsics) |
| 5 | For accuracy problems, change the algorithm first, not the scalar type | Measured on this machine (B) | 0 (engine selection) |

The five, expanded:

1. **Scale decides everything**: the same code is 6–14x when L1/L2-resident and 0.8–1.4x on a 24 MB
   working set, where it is bound by memory bandwidth. A speedup quoted without a working-set size means
   nothing.
2. **The call boundary is a first-order bottleneck** — not lane width, not the ABI, not FLOPs: library call
   1.657 / inlined 0.565 / 256-bit packed 0.214 ns·element⁻¹.
   The boundary cost is ≈4.2 cycles, 3x the ≈1.3 cycle packing win. Today the code reaches about
   1% of FP peak.
3. **`-fno-math-errno` is the largest single win at zero semantic cost**: `fv_vec3_normalize` goes from
   40 to 26 instructions (my recomputation), loops containing `sqrtf` move from "vectorization refused" to
   "vectorized at 32 bytes", and the numbers are bit-identical.
4. **What blocks auto-vectorization is not sqrt but the "conditionally executed division"** — the branch of
   the zero-vector guard. A hand-written masked SIMD kernel gets 8.4x / 14.5x without relaxing any FP
   semantics, and is bit-identical to the scalar path.
5. **For accuracy problems, change the algorithm, not the scalar type**: for the harmonic oscillator's
   energy drift, semi-implicit Euler gives 5.522e-2 in float vs 5.521e-2 in double (switching to double
   buys 0.02%), while velocity Verlet gives 2.744e-3 (20x).

**One methodological rule**, see §4.8: **a benchmark that does not verify its results produces very
convincing wrong conclusions.** Two of the three paths in my first adjudication experiment were
miscalculated — off by exactly 4x — and the "performance numbers" looked perfectly reasonable.

**One-line recommendation**: give the library **batch entry points** first, AoS 4-way interleaved `_x4`
before anything else, and turn on `-fno-math-errno` at the same time.

**Do not** start with hand-written intrinsics, do not switch to double, do not introduce templates or a pure
SoA contract.

**Added after this report was written**: §2.2 records the measurements for the fused, in-place and
three-operand entry points. They confirm the inference §2 ends with — at large sizes the lever is
traffic, and fusing three passes into one is worth 2.43-2.55x at n=1M and 2.53-2.63x at n=8M.

---

## 1. How This Report Came About (Read This First for Credibility)

### 1.1 Division of Labor

| Role | Content | Artifact |
|---|---|---|
| Research A | SIMD / compiler vectorization / Zen 2 hard numbers / reference-library dissection | `docs/research/simd-hardware.md` |
| Research B | Approximation algorithms, numerical stability, porting competitive-programming constant-factor tricks, physics-simulation interface needs | `docs/research/algorithm-tricks.md` |
| Research C2 | AoS/SoA/AoSoA, memory hierarchy, batch-API forms | `docs/research/layout-memory.md` |
| Implementation D | Three-tier implementation (scalar / auto-vectorized / hand-written AVX2) × three sizes, measured | `docs/research/measured.md` + `bench/` |
| Me | Planning, acceptance, adjudication experiments, synthesis and judgement | This report |

Process disclosure: the first "data layout" researcher delivered nothing in 80 minutes, so I stopped it and
swapped in C2; C2 delivered in about 15 minutes.

For C2's output I trust no quantitative conclusion that has not been cross-confirmed.

### 1.2 Three Hard Limits of Measurement on This Machine (Know These Before Reading Any Number)

1. **perf is unavailable**: `kernel.perf_event_paranoid = 4`, no CAP_PERFMON, sudo needs a password. So
   there are no IPC / cache-miss / branch-miss counters.
   - The report substitutes "GCC's verbatim vectorization verdicts + disassembly + tiered memory-bandwidth
     measurements" as evidence.
   - To fix: `sudo sysctl -w kernel.perf_event_paranoid=1`, then `make bench-perf`.
2. **Laptop frequency scaling**: the measured clock drifts between 1.64 and 3.97 GHz, and two measurements
   of the same code can differ by 25%. Therefore **no difference below 20% is ever treated as a
   conclusion**.
   - Absolute ns figures are not comparable across runs; only ratios within one run are used.
3. **Shared machine**: loadavg 3–6 during measurement. Every key conclusion has at least two independent
   sources cross-confirming it.

### 1.3 What I Recomputed Independently (Not Second-Hand)

| Recomputation | Verdict |
|---|---|
| Effect of `-fno-math-errno` on `fv_vec3_normalize` | Order of magnitude holds, accounting differs; the report uses my numbers |
| Library calls block vectorization | Holds |
| `__builtin_cpu_supports` return values | Holds |
| Reproducibility of D's three-tier benchmarks | Holds |
| Whether hadd reduction is fast or slow (A and D disagree) | **D's claim not accepted**, see §3.4 |

Item by item:

- `-fno-math-errno` on `fv_vec3_normalize`: the subagent reported 48 → 22 instructions; I measured
  **40 → 26 instructions** (objdump count).
- Library calls block vectorization: the subagent cited `missed: statement clobbers memory`; I reproduced
  the identical message with `-fopt-info-vec-missed`.
- `__builtin_cpu_supports` return values: the subagent said it returns a bitmask and you must write
  `!= 0`; I measured `sse2=16 avx2=1024 avx512f=0`.
- Reproducibility of D's three-tier benchmarks: the subagent reported x6–14 (N=8192) / x0.8–1.4
  (N=1<<20); I reran `make bench` myself and got the same shape.
- Whether hadd reduction is fast or slow: A said slow, D said 22% faster; I wrote a decisive experiment,
  and the verdict is **hadd is about 2x slower**.

### 1.4 Not Yet Verified (Do Not Use as a Basis)

- AoS↔SoA **conversion cost** has not been measured. This is the biggest unknown right now, and it
  directly determines whether a SoA entry point is worth it.
- Cross-domain (integer↔float) forwarding penalties, store-forwarding failure penalties, AVX2 downclocking:
  not measured — perf is unavailable here. **Until you measure it, do not use the internet's "AVX2 downclocks" as a
  design basis.**
- Zen 2's `vgatherdps` is slow: A measured about 10 cycles/instruction (holds); but "the index variant is
  always worse than the mask variant" is an inference.
- The NaN semantics difference between `minps`/`maxps` and `fminf`, and the fallback instruction for
  `fv_abs`: not verified.
- "double would halve SoA batch throughput": inference, not measured.
- Internal drift of the evidence: the same bandwidth data pasted in `raw/membw.txt` and `measured.md`
  §3.3 differs by 10–15%, write 12.33 vs 14.05.
- `raw/` is authoritative. **A gap under 10% should not be treated as a conclusion on this machine.**

---

## 2. First Conclusion: Scale Decides Everything

D's three-tier measurements, in ns/element, median of 5 runs. Tier 1 = the existing scalar API, tier 2 = a
naive for loop + auto-vectorization, tier 3 = hand-written AVX2.

After rerunning `make bench` myself, the shape is identical.

| Operation | Tier 1 | Tier 2 | Tier 3 | Tier 1 | Tier 2 | Tier 3 |
|---|---|---|---|---|---|---|
| (N) | 8192 | 8192 | 8192 | 1048576 | 1048576 | 1048576 |
| add | 1.52 | 0.24 (**x6.30**) | 0.23 (**x6.54**) | 5.33 | 5.48 (x0.97) | 6.77 (x0.79) |
| scale | 2.21 | 0.16 (**x13.46**) | 0.17 (**x13.23**) | 4.19 | 4.32 (x0.97) | 4.01 (x1.05) |
| dot | 2.68 | 0.76 (x3.52) | 0.32 (x8.34); 4 accumulators **0.14 (x18.78)** | 3.53 | 2.59 (x1.36) | 2.60 (x1.36) |
| len | 1.96 | 1.56 (x1.25); **0.24 (x8.16)** with `-fno-math-errno` | 0.20 (**x9.63**) | — | — | — |
| normalize | 2.85 | 2.50 (x1.14) | 0.34 (**x8.41**) | — | — | — |
| collision | 7.96 | 3.97 (x2.01) | 0.55 (**x14.51**) | 9.56 | 7.11 (x1.34) | 6.78 (x1.41) |

`collision` is a batch microkernel for "normalize the normal + reflect the velocity", representing the
shape of a real physics loop.

**Why all three tiers flatten out on large arrays**: `bench/membw.c` measures this machine's three-stream
bandwidth at only 6.0–6.9 GB/s, with memcpy at 8.9–17.2 GB/s.

A vec2 add touches 24 B/element of logical memory traffic (2 reads, 1 write); at 6.9 GB/s the theoretical
floor is ≈**3.5 ns/element**.

Measured tier 1 is 5.3 and tier 3 is 6.8 — all three tiers sit right against that bandwidth ceiling.

**The knee must be computed from "L3 visible to one thread", not from the per-chip total.** An earlier
version got this wrong here; what was wrong is explained below.

This machine's L3 is **4 MiB/CCX × 2 CCX = 8 MiB per chip**.

`/sys/devices/system/cpu/cpu0/cache/index3` reports `level=3, size=4096K, shared_cpu_list=0-7`.

A thread runs on one CCX and **can only reach 4 MiB**.

So the knee is **4 MiB ÷ 24 B ≈ 1.75×10⁵ active elements**. The earlier version wrote 3×10⁵, computed from
8 MiB, i.e. twice as high.

That knee has direct measurement behind it: an in-process A/B harness measured **38–57 GB/s at 2^17
(3 MB working set) and 5.6–8.8 GB/s at 2^18 (6 MB)** — the cliff sits exactly on that line.

Using the real traffic below (32 B/element), the knee moves further forward to about 1.3×10⁵.

**One correction that has to be stated**: the 8–14x at N=8192 is an **L2-resident figure for a working set
of about 192 KB**, not an L1 peak.

That size happens to land in the range of a real engine's "active bodies per frame", so the conclusion
transfers; but it is not a promise that "any speed is fine".

**Inference (the most useful one for you)**:

- A physics engine processes **a few thousand to a few tens of thousands** of active bodies per frame → it
  falls on the "SIMD pays off" side, and batch entry points are worth it.
- Once you reach **millions**, changing instructions is useless; **the only lever is reducing memory
  traffic**: skip sleeping bodies, lower precision, tile the work, and do not store quantities you do not
  need. This is the same thing as B's "sleeping-body skipping" and C2's "mask variant", stated three ways.

### 2.1 Byte Accounting: Nominal Traffic vs Real Traffic with RFO (Where the `_stream` Win Physically Comes From)

Everything above uses **nominal** bytes (reads + writes). But on a write miss, an ordinary store makes the
CPU **fetch the whole cache line first** — that is RFO (read-for-ownership). Real traffic is "one write plus
one read".

So the per-element accounting needs two columns:

| Entry point | Nominal B/element | Real B/element |
|---|---|---|
| `add_x4` / `sub_x4` / `lerp_x4` / `madd_x4` / `add_soa` | 24 (16 read + 8 write) | **≈32** (+8 RFO) |
| `scale_x4` / `scale_soa` | 16 (8 read + 8 write) | ≈24 (+8 RFO) |
| `dot_soa` (vec2, and the same-named entry points for vec3/vec4) | **16**, four `float` arrays ax/ay/bx/by × 4 B | 16 (read-only, no writes → no RFO) |
| `normalize_soa` | **16** (8 read + 8 write) | ≈24 (+8 RFO) |
| `*_stream` (NT store) | 24 (same as the ordinary version) | **24** (no RFO) |

Notes per entry point:

- `add_x4` / `sub_x4` / `lerp_x4` / `madd_x4` / `add_soa`: the output array gets read one extra
  time, and that is exactly the 8 B an NT store saves.
- `scale_x4` / `scale_soa`: same, just with one fewer read stream.
- `dot_soa`: **accounting reminder: it is 16 B, not 32 B**. Each element reads only 4 floats; it was once
  recorded as 32 B, which was wrong.
- `normalize_soa`: **accounting reminder: 8 read + 8 write = 16 B, not 16+16**. It was once recorded as
  16+16, which was wrong.
- `*_stream`: the real win is not "faster per byte" but **25% fewer bytes moved**.

Two criteria follow, and `bench/ab`'s reading rules are written from them:

1. For the same add operator, the **real** GB/s of the ordinary and `_stream` versions should be
   essentially equal, both pinned at DRAM bandwidth, while the **nominal** GB/s differ by 1.33x.
   If both nominal and real differ by 1.33x, the gain is not coming from traffic. This is the ruler for
   spotting false positives.
2. An NT store only gets that 25% when the **working set is far larger than L3 and the write stream is a
   clean one-shot write**. When the working set fits in L2/L3 it is **2–8x slower** than the ordinary
   version; measured inside L2 it is **0.13–0.51x**, and all 6 entry points are judged FAIL (slower).
   The ceiling is counted per entry-point shape: add family **32/24 = 1.333**, scale and normalize
   **24/16 = 1.500**.
   Measured PASS points are 1.29–1.58x, and paired 4M-tier IQRs are 66–98%; **a single point exceeding the
   ceiling by less than 20% is noise**. See `docs/acceptance-protocol.md` §9.

### 2.2 Fusion: the Lever That Is Actually Left at Large Sizes

§2 ends with an inference: once the working set is in the millions, changing instructions buys
nothing and the only lever is memory traffic. This section is that inference measured, and it is
where the library's fused, in-place and three-operand entries come from.

The entries, with the full contract beside each declaration in `include/flashvec.h`:

| Entry | Expression |
|---|---|
| `fv_vec2_madd_normalize_soa` | `out = normalize(a*s + b)`, one pass |
| `fv_vec2_addmul_x4` / `_soa`, `fv_vec3_addmul_soa`, `fv_vec4_addmul_x2` | `x += v*s`, in place |
| `fv_vec2_maddv_x4` | `o = a*b + c` — three vector operands; `madd`'s third operand is a scalar |

Byte accounting decides what each one can be worth, and it is counted per entry shape:

| Entry | Real B/element | The two-call workaround | Ceiling |
|---|---|---|---|
| fused `normalize(a*s+b)` | 32 (read a and b 16 + write 8 + RFO 8) | 80 (three passes) | 80/32 = 2.5x |
| in-place `addmul` | 24 (read v 8 + read x 8 + write 8; the store hits the line the read owns, so no RFO) | 56 (`scale` 24 + `add` 32) | 56/24 = 2.33x |
| `maddv` | 40 (a, b, c 24 + write 8 + RFO 8) | 64 (two passes, each paying its own RFO) | 64/40 = 1.60x |
| `maddv` with `_stream` | 32 (the same, RFO dropped) | 64 | 64/32 = 2.00x |

The two workarounds do not carry the same baseline: quoting 56 B for both would understate `maddv`.

The denominator first, because every ratio below is a fraction of it. `bench/membw.c`, one thread,
taskset-pinned, best of 5:

| Working set | triad | triad16 (with RFO) | read | memcpy |
|---|---|---|---|---|
| 96 MB | 7.88 GB/s | 10.51 GB/s | 13.28 GB/s | 12.93 GB/s |
| 48 MB | 9.16 GB/s | 12.21 GB/s | 14.31 GB/s | 14.31 GB/s |

`triad` counts logical bytes, 12 B/element; `triad16` is the same run with the write's
read-for-ownership allocated in, 16 B/element. These readings are higher than the 6.0–6.9 GB/s
three-stream figure quoted in §2 — both come from `bench/membw.c`, and this section quotes the
newer run.

The element-wise baseline: `fv_vec2_add_x4` on a 96 MB working set runs at **3.33 ns/element**,
about **9.6-10.3 GB/s** against a machine triad16 ceiling of 10.51 GB/s at the same 96 MB working set -- **91-98%** of it. That is the number worth
planning against: at that size a batch entry point is already near the wall, and no amount of
instruction work moves it.

Measured, with the harness this report has used throughout — single thread, taskset-pinned, the
default tier `-O3 -fno-math-errno -ffp-contract=off` **without `-march`**, alternating A/B inside
one process:

| Change | n=8192 | n=1M | n=8M | A/A floor |
|---|---|---|---|---|
| `_stream` vs an ordinary store | — | 1.48–1.60x | 1.54–1.60x | 0.97–1.07 |
| fused `normalize(a*s+b)` vs three calls | 1.73x | 2.43–2.55x | 2.53–2.63x | — |
| the same, plus `_stream` | 1.68x | **3.53–3.68x** | **3.72–3.78x** | — |
| in-place `x += v*s` vs `scale` then `add` | 1.06–2.18x | 2.05–2.17x | 1.92–2.03x | 0.97–1.06 |
| `maddv` vs `mul` then `add` | 1.60–1.67x | 1.61–1.79x | 1.59–1.76x | ~1.00 |
| `maddv_stream` vs `mul` then `add` | 1.16–1.37x | 2.16–2.23x | **2.34–2.41x** | ~1.00 |

The floor belongs beside every row because each row is a ratio between two different bodies of code:
when the floor itself reads 0.97–1.07, a 1.06x row carries no information.

Reproduction: the `addmul` and `maddv` rows are `make bench-addmul`; the fused row is a separate
round on the same harness; the bit-level identity claim is `make test-maddnorm`, which compares
against the stepwise calls with `memcmp`.

Three limits, and they are the part that gets lost:

1. **`_stream` is slower at small sizes — not slightly, 2–8 times** (0.13–0.37 at n=8192). An NT
   store bypasses the cache, so inside L2 it has no RFO to save and loses the cache hit as well.
   Use it only when the working set is at least L3, 4 MiB per CCX.
2. **`maddv_stream` at n=8192 is faster relative and slower absolute**: 1.16–1.37x against the
   two-pass workaround, but 0.99 ns/element against 0.58 for the non-stream form. Inside L2 the
   non-stream form is the one to call; that ratio only says "still faster than two passes".
3. **The byte account is a ceiling, and 20% above it is noise.** `add` / `sub` / `lerp` /
   `madd` top out at 32/24 = 1.33x, `scale` / `normalize` / the matrix transform at 1.50x,
   in-place `addmul` at 2.33x, `maddv` at 1.60x and `maddv_stream` at 2.00x. The `_stream`
   row reads 1.48–1.60x, above the 1.33x ceiling for the `add` family; that excess is real — an
   NT store also keeps the output from polluting the cache — and the account above does not count
   it. Do not turn a single point into "above the theoretical limit".

What this means for the plan: §2's inference said that at millions of elements the only lever left
is reducing memory traffic. That lever now has two concrete forms. Fuse the element-wise passes.
Stop writing the output into the cache when it is not read back. The first works at every size and
is the larger win; the second only above L3. Neither reaches a reduction: changing the order of a
`dot` changes its result, which is the §3.4 point again.

---

## 3. The Five Main Threads, Expanded

### 3.1 The Call Boundary Is a First-Order Bottleneck (Three Independent Lines of Evidence + My Recomputation)

The same AoS data (24 KB, entirely in L1), within a single run:

| Form | ns/element |
|---|---|
| Existing API (`fv_vec2_add`, cross-TU call) | 1.657 |
| Same code placed in a header and inlined | 0.565 |
| 256-bit packed (4 vec2 per operation) | **0.214** |

- **Boundary cost = 1.657 − 0.565 = 1.09 ns/element ≈ 4.2 cycles**; **packing win = 0.565 − 0.214 =
  0.35 ns/element ≈ 1.3 cycles**. The boundary is 3x the packing win.
- Structural evidence, not an estimate: in a consumer that links the existing `libflashvec.a`, a for loop
  containing the call compiles to **5 calls to external `fv_vec2_xxx`**; with `-flto` it compiles to
  **0**.
- GCC refuses vectorization outright: `missed: statement clobbers memory`. **An opaque `extern "C"`
  call is treated as "may write memory"**.
- **The problem is not the ABI**: `fv_vec2` (8 B) goes in 1 xmm register, `fv_vec3` (12 B) in 2, neither
  spills to the stack; the body of `fv_vec2_add` is exactly `addps %xmm1,%xmm0; ret`.
- **Do not change the signature to "reduce copying".**
- Current throughput ≈ 2 FLOP ÷ 1.657 ns ≈ **1.2 GFLOP/s** against this machine's FP32 peak of about
  130 GFLOP/s, i.e. **only 1% of peak**. The bottleneck is still far from "arithmetic".

**Action**: batch entry points, and **AoS 4-way interleaved (`_x4`) before pure SoA**. One ymm holds
exactly 4 `fv_vec2`, with zero reshuffling for consecutive elements, and the caller does not change a
single line — see §4.6.

### 3.2 -fno-math-errno: Zero Semantic Cost, Doable Today

- By default `sqrtf` carries errno semantics → a branch plus a library call. My recomputation:
  `fv_vec3_normalize` goes from **40 to 26** instructions; A reported 48→22, the same direction.
- D's end-to-end effect: the same len loop goes **x1.25 → x8.16** (1.56 → 0.2395 ns/element), with
  **bit-identical numbers** (DIFF all 0).
- Why not `-ffast-math`: it also turns on `-ffinite-math-only` and reassociation.
- `-ffinite-math-only` assumes no NaN/Inf, while flashvec's contract (a zero vector normalizes to a zero
  vector, no NaNs emitted) depends on those semantics.
- **If you want speed, ship a `_fast` variant; do not touch the global switch.**
- Note: under `-ffast-math`, GCC compiles `1.0f/sqrtf(x)` into `vrsqrtss` plus 1 Newton step (B
  measured this with objdump), i.e. it silently drops to about 3 ulp.
- Whereas `1.0f/x` is `vdivss` under any option; GCC does not generate `rcpps`.

### 3.3 What Blocks Auto-Vectorization Is the "Conditionally Executed Division", Not sqrt

The verbatim GCC verdict D captured: under `-O3 -march=native`, len / normalize / collision all report
`missed: not vectorized: control flow in loop`.

| Operation | `-fno-math-errno` alone (zero numerical cost) | Requires `-ffast-math` to unlock | Hand-written masked SIMD |
|---|---|---|---|
| len | Yes: x1.25 → **x8.16** | — | x9.63 |
| normalize | No (still not vectorized) | Yes: x1.14 → x9.84 | **x8.41** |
| collision | No | Yes: x2.01 → x13.53 | **x14.51** |

**Conclusion**: hand-written AVX2 uses `_mm256_cmp_ps` + `_mm256_and_ps` to express the "zero-vector
guard" as a branchless mask.

It reaches 8.4x / 14.5x without relaxing any FP semantics, and is **bit-identical** to the scalar path
(IDENT bits_equal=1).

This answers "should flashvec have a SIMD implementation branch": **yes, but only inside batch entry
points**. A single-element `fv_vec2_normalize` will always be dominated by call overhead.

### 3.4 Writing a Reduction: 4 Accumulators + One Final Reduction (Plus: How I Adjudicated a Contradictory Conclusion)

A says "hadd reduction is slow"; D says "hadd is 22% faster than the single-accumulator shuffle version".
The two conclusions have opposite signs.

I wrote a decisive experiment: all three paths **verify correctness first**, then time; `taskset -c 2`,
200 runs taking the minimum.

| dot reduction form (8-wide) | n=8192 | n=262144 |
|---|---|---|
| hadd×2 reduction every 8 elements | 0.2507 | 0.1359 |
| shuffle+add reduction every 8 elements | 0.1284 | 0.1068 |
| **4 independent accumulators + one final reduction** | **0.0868** | **0.0912** |

All three paths give 2080.0 on the same reference input, which is the correct value.

- **Verdict: hadd really is slow (about 2x)**, agreeing with A. D's 22% difference falls inside the 25%
  frequency-scaling noise band D itself declared, and it has the opposite sign, so **it is not accepted**.
- The more important layer: **the real lever is reduction *frequency*, not the reduction *instruction***.
  Reducing every 8 elements (0.13–0.25) versus reducing once at the end (0.087) is a 1.5–2.9x gap. That
  matters far more than the hadd/shuffle argument.
- D's own finding on the same point: single accumulator 0.3213 → 4 accumulators 0.1428 (**x2.25**).
- A's "you need ≥10 independent chains to saturate FMA" and B's "8-way accumulator" are **the same trick**,
  and it wins on both throughput and numerical stability (§3.5).

### 3.5 Accuracy: Changing the Algorithm > Changing the Scalar Type

| Scenario | float | double | Conclusion |
|---|---|---|---|
| Harmonic oscillator, 6e4 steps, semi-implicit Euler (max energy drift) | 5.522e-2 | 5.521e-2 | double improves it by **0.02%** |
| Same oscillator, **velocity Verlet** | **2.744e-3** | 2.740e-3 | changing the integrator improves it **20x** |
| x += v*dt, 1e6 steps (relative error) | 2.422e-3 | 0 (double accumulation) | changing the accumulation structure (Kahan) = 2.597e-8 |
| Summing 1e6 values of 0.1f | 9.863e-3 (0.926 ns/add) | — | **8-way accumulator: 1.036e-3 at only 0.338 ns/add** (both faster and more accurate) |

- A symplectic integrator's drift is **algorithmic and bounded**; multiply the step count by 10 and it
  barely changes. float→double does not buy that.
- Kahan vs Neumaier: the textbook counterexample {1, 1e30, 1, -1e30} has exact value 2.
- The three answers: naive = 0, **Kahan = 0, the small term swallowed by the large one**, Neumaier = 2.000,
  correct. Default to Neumaier; it is only about 20% slower than naive.
- **Action**: keep `fv_real` as float. Add accumulator forms (`fv_kahan2` / `fv_neumaier2`) and batch
  reduction entry points on the engine side; do not change default semantics.

---

## 4. Ideas Borrowed from Elsewhere (My Part)

Each entry gives: source → what it originally solved → what form it takes in flashvec → cost.

**4.1 Columnar databases / vectorized execution engines (vector-at-a-time)**

They found long ago that a row-at-a-time processor's bottleneck is **per-tuple dispatch overhead**, not the
arithmetic itself. The fix is to process a batch at a time and amortize the dispatch.

This is completely isomorphic to flashvec. Our "per-tuple dispatch" is "one call per vec2", and the measured
4.2-cycle boundary cost is 3x the 1.3-cycle packing win. **Lesson: change the processing granularity before
you change the instructions.**

Databases offer an even closer analogy: their "vectors" are columns (SoA), while the caller holds rows
(AoS). DuckDB simply performs one layout conversion between the operator and the storage layer, and openly
admits the conversion cost. So we should not pretend conversion is free either — see §1.4.

**4.2 The compiler's "ordering of optimization levers"**

A compiler does not do register allocation before inlining; it orders passes by payoff magnitude. Sorting
the data from the three inputs by "how much does changing this get me back":

~~~
batch entry points (eliminate cross-TU calls + packing)   ≈ 6–14x    ← one change in call shape
-fno-math-errno (zero semantic cost)                      ≈ 1.2–8x   ← one compiler flag
reduction frequency (4 accumulators / final reduction)    ≈ 1.2–2.2x ← loop structure
mask instead of branch (hand-written SIMD only)           ≈ 8–14x    ← a few dozen lines of intrinsics,
                                                                       and only on top of the three above
~~~

**Counter-example**: starting with hand-written intrinsics means micro-optimizing at 1% utilization, and
the call boundary eats most of the gain. The boundary cost is 3x the packing win.

**4.3 Three transferable and three non-transferable ideas from competitive programming (Luogu / CF)**

- Transferable: **hoisting invariants**. In `fv_vec2_project(v, onto)`, `len_sqr(onto)` is constant
  across an entire batch of calls; provide `fv_vec2_project_inv(v, onto, inv_denom)`.
- Transferable: **shrink the problem's scale first**. Compare distances with `dist_sqr` (already present)
  rather than `dist` — exactly the reason Box2D provides `b2DistanceSquared`.
  - **Current risk: `fv_vec2_len` and `fv_vec2_dist` are the two functions in this library most likely
    to be misused.**
- Transferable: **the real purpose of manual unrolling is to break dependency chains**, not to write fewer
  lines: 8-way accumulator 0.338 vs naive 0.926 ns/add.

- Not transferable: **trading computation for table lookups** reverses under SIMD.
  - sin 8-wide: polynomial+FMA 1.663 ns/value (error 1.19e-7) is **faster than** gather-LUT 4096 at 2.042
    (5.36e-7), and also faster than glibc `sinf` at 4.164.
  - Reason: random indices need `vgatherdps`, which is a Zen 2 weak spot at about 10 cycles/instruction.
  - Table lookup only pays off with "scalar + a table that stays resident in L1".
- Not transferable: **replacing modulo/division with bit operations**. A vector library has no index
  arithmetic; that belongs to the engine's broad-phase (hashing / sort keys).
- Conditional: **branchless min/max** has a real payoff.
  - near-type comparisons on synthetic data at 50% hit rate: 5.349→2.803 ns/call, 1.9x.
  - But that is worst-case synthetic data; and the NaN semantics difference between `minps`/`fminf` is
    **unverified**, so consult the SDM before adopting it.

**4.4 Two tricks from ML kernels**

- **Multiple accumulators to break dependency chains** (the standard GEMM technique): here it is the same
  thing as the numerical-analysis "multiple partial sums" — one trick that buys throughput and accuracy at
  once. **Changes where both dimensions improve in the same direction are always done first.**
- **The tail must be designed explicitly**: the n%4 and n%8 tails, masked loads, and the "write the original
  value back outside the mask" semantics all have to be written into the header. Otherwise callers will
  assume a mask means "do not write" and produce dirty data.

**4.5 Data layout is not a global property, it is a per-stage property (borrowed from simdjson)**

simdjson uses different layouts in different stages: structural scanning uses a bitmask stream, and the
replay stage expands as needed.

Likewise, flashvec should not declare "we are AoS" or "we are switching to SoA": **storage is AoS; inside
hot loops you materialize a batch in whatever layout the reader wants** — that is enough.

The PhysX practice C2 summarized is closer to engineering reality: not global SoA, but "for each hot loop,
open a separate array for whoever reads it".

**4.6 The ECS/engine practice of "splitting axes per consumer" + One Concrete Design That Resolves the Contradiction**

Resolve the contradiction between "D measured 8–14x for SoA batches" and "physics-engine call sites already
hold AoS particle arrays" without choosing one:

~~~c
/* Tier one: pass AoS as-is. One ymm in the library holds 4 fv_vec2, zero reshuffling, zero caller changes */
void fv_vec2_add_x4  (const fv_vec2 * __restrict a, const fv_vec2 * __restrict b,
                      fv_vec2 * __restrict o, size_t n);
/* Tier two: reduction-style ops (dot / len / normalize); components are contiguous, so SoA is natural */
float fv_vec2_dot_soa(const float *ax, const float *ay,
                      const float *bx, const float *by, size_t n);
/* Every batch entry point needs a sparse variant, or sleeping-body skipping / broad-phase cannot use it */
void fv_vec2_add_x4_idx (..., const uint32_t *idx, size_t m);  /* true skip; the price is random memory access */
void fv_vec2_add_x4_mask(..., const uint8_t *mask, size_t n);  /* whole-line RMW; original kept outside the mask */
~~~

**4.7 Turning "reproducibility vs performance" into an explicit contract (D's finding, raised here to a principle)**

D measured: for the same collision logic, the functional-API version and the inlined version differ in their
**upper bits on 4948 / 8192 elements** (about 7 ULP). The cause is that inlining lets the compiler do FMA
contraction, while the cross-TU call does not.

This is not a bug, it is **an unwritten contract**. Recommend stating it explicitly in the header:

> The element-wise API guarantees bit-for-bit reproducibility (same input, same output); batch / SIMD entry
> points do not guarantee bit-identical results, only a bound on error.

"Two runs of the same scene give different results" during physics-engine debugging is an expensive bug;
this contract removes it.

**4.8 One meta-lesson: a benchmark must verify its results**

I made the mistake twice myself. The first time, two of three paths gave results **4x too large**: after
hadd×2 the vector lanes are a repeated fill of [S0,S0,S0,S0,S1,S1,S1,S1].

I counted them 4 times with a generic 8-lane reduction. After fixing that, they were still off by 2x.

And the **performance numbers of both versions looked perfectly "reasonable"**. Had I not verified the sum,
I would have written a wrong hadd conclusion into the report.

Three corollaries should be fixed in place:

1. Every microbenchmark must **print a result checksum at the same time**, and compare it against a
   reference implementation.
2. Any conclusion whose difference is smaller than the noise band (20–25% on this machine) does not hold,
   even if the direction is "intuitive".
3. Any performance claim must report **working-set size / cache level / compiler options** at the same time,
   or the number is meaningless (this is the whole point of §2).

---

## 5. Roadmap (Staged, with Acceptance Criteria)

### P0 · Zero Semantic Cost, Doable Today

**Add `-fno-math-errno` to the Makefile, optionally also `-O3`**

- Expected: len-type loops x1.25 → x8.16, results bit-identical.
- Cost: give up "math functions set errno". In this library the sqrt argument is always ≥0.
- Acceptance: `make test` all green + bit-for-bit comparison of old and new artifacts.

**Add `fv_vec2_min/max/clamp`, `fv_vec2_madd`, `fv_aabb2_merge/overlaps`**

- Expected: close the gap for broad-phase and impulse accumulation. The current library has **no
  min/max/clamp at all**.
- Cost: a dozen or so new symbols.
- Acceptance: unit tests covering NaN / boundaries.

**Fix the `FV_EPS` dimensional bug** (§6, decision 1)

- Expected: eliminate `project`'s 1000x loose threshold.
- Cost: one line.
- Acceptance: a new assertion checking the behavior of `project(v, tiny onto)`.

### P1 · Batch Entry Points (The Biggest Single Step)

**AoS 4-way interleaved: `fv_vec2_add/sub/scale/lerp_x4`**

- Form: AoS, zero reshuffling, `__restrict`.
- Expected: x6.3–13.2 at N=8192. Hand-written intrinsics are only 0–8% ahead of auto-vectorization.
- Cost: ≈60 lines + tail handling.
- Acceptance: **bit-identical** to the scalar path; n%4 boundary tests.

**Reduction-style SoA: `fv_vec2_dot/len/normalize_soa`**

- Form: branchless mask + 4 accumulators.
- Expected: normalize x8.41, collision x14.51, dot x18.8.
- Cost: ≈150–250 lines.
- Acceptance: bit-identical to the scalar path; zero vectors produce no NaN.

**Index / mask variants**

- Expected: the precondition for sleeping-body skipping.
- Cost: +40% symbol count.
- Acceptance: the original value is written back outside the mask, not "not written".

**Important**: do `_x4` first, and **do not** turn around and require the caller (a future physics engine)
to change its data layout first.

### P2 · Optional Fast Paths (Once There Is Performance Data)

- `fv_vec2_normalize_fast` (`rsqrtps` + 1 NR): 0.3082 vs 0.3395 ns/element (about 10%), error ≈2 ulp
  (max deviation 1.192e-07 over 4M samples).
- `include/flashvec_math.h` (a separate header): atan with a 7-term polynomial 3.71 vs glibc 10.15 ns/call
  (error 1.16e-6 rad); sin poly7 1.663 ns/value.
- **The input domain must be stated**: sin requires [-π, π]. float argument reduction has error 8.742e-5
  over [-1000π, 1000π], while double reduction gets 1.192e-7.
- Runtime dispatch: a single function with `__attribute__((target("avx2,fma")))` +
  `__builtin_cpu_supports`.
- `__builtin_cpu_supports` returns a bitmask, so you must write `!= 0`; the disassembly is one global
  load plus a test, negligible overhead.

### P3 · Not Recommended / Defer

- Global `-march=native`: binds the binary to znver2, and another machine gets SIGILL. If you ship prebuilt
  libraries, use target attributes / multiversion dispatch.
- Full-pipeline LTO: it removes the existing API's call boundary (5 calls → 0), but it requires the
  **consumer** to also use `-flto`; hand out a prebuilt static library and it stops working.
- AoSoA / pure SoA contract / expression templates / memory pool / SIMD complex numbers library-wide:
  see §7.

---

## 6. Four Decisions That Need Your Call

1. **How to fix the `FV_EPS` dimensional bug.** `fv_vec2_project` uses
   `fv_abs(len_sqr(onto)) <= FV_EPS` to decide "onto is near zero", but `len_sqr` is a **squared
   quantity**.
   - The equivalent length threshold is `sqrt(1e-6) = 1e-3`, **1000x looser** than `normalize`'s linear
     threshold, yet the header comment describes `FV_EPS` as a uniform "close enough to zero" threshold.
   - This is a real bug B found and I confirmed on review; it is also my own oversight when writing the
     code.
   - Options: (a) align with `normalize` (`len_sqr <= FV_EPS*FV_EPS`); (b) split into two constants,
     `FV_EPS` and `FV_EPS_SQR`; (c) keep it as is and fix the comment. **I lean toward (b)**, which is
     hardest to get wrong again.
2. **Whether to introduce `_x4` / `_soa` APIs**: trade API surface for 6–14x. I lean toward "yes, and add
   only `_x4` and the reduction-style `_soa`", because zero caller changes is the precondition for this
   library actually being adopted.
3. **Whether to go `-march=native` / ship prebuilt libraries**: I lean toward "keep the library on a
   generic x86-64 baseline + target-attribute dispatch for key functions", avoiding "fast on your machine,
   SIGILL on someone else's".
4. **Whether to put `-fno-math-errno` into the default compile options**: I lean toward "yes". Zero
   numerical cost, verifiable, just compare bit for bit.

---

## 7. What We Explicitly Will Not Do (Three "Not Recommended" Lists Merged)

| Not recommended | Reason |
|---|---|
| Pad `fv_vec2` to 16 B for SIMD | 8 B fits exactly 4 per ymm, which is optimal on this machine; padded, `_x4` holds only 2 per register and the win is halved |
| Change `fv_real` to double for accuracy | Energy drift improves by only 0.02%, while batch throughput may halve (inference) |
| Put `-ffast-math` in the default options | It assumes no NaN/Inf, while the current contract relies on predictable NaN and comparison semantics |
| Wrap the single-element API in AVX2 dispatch | Single-element calls are always dominated by call overhead |
| Hand-write intrinsics before batch entry points land | The gain gets eaten by the call boundary, which costs 3x the packing win |
| Introduce expression templates or a pure SoA contract | It destroys the pure-C-header premise, and the gain is in batch entry points, not in templates |
| Add a memory pool / arena to the library | The ceiling is about 50 ns of one malloc; the risk is ownership bugs and C/C++ interop cracks |
| Hand-write Quake rsqrt | `rsqrtss` + 1 NR strictly dominates it, see below |
| Use approximate trig functions in a physics simulation without stating the input domain | Argument reduction is float's real bottleneck |
| Make the vector library "batch-only usable" | The scalar API is the anchor for debugging and the correctness baseline |

Three expanded:

- Expression templates / a pure SoA contract would destroy the premise of "pure C header + `extern "C"` +
  value semantics + no allocation", and the gain is in batch entry points, not in templates.
- A memory pool / arena has a ceiling of about 50 ns of one malloc, and the risk is ownership bugs and
  C/C++ interop cracks.
- Hand-written Quake rsqrt: `rsqrtss` + 1 NR strictly dominates it, 1.848e-7 vs 1.750e-3, and even after
  two Newton iterations it is still 25x off.
  - The two magic numbers `0x5f3759df` / `0x5f375a86` show no measurable difference over 400k samples.

---

## 8. Appendix

### 8.1 Material Index

| File | Content |
|---|---|
| `docs/research/simd-hardware.md` | Reference-library dissection, auto-vectorization conditions, intrinsics engineering details, Zen 2 hard numbers, counter-example list |
| `docs/research/algorithm-tricks.md` | Approximation and approximation methods, numerical stability, assessment of porting competitive-programming tricks, physics-simulation interface needs, not-recommended list |
| `docs/research/layout-memory.md` | AoS/SoA/AoSoA trade-offs with real cases, memory hierarchy, three batch-API routes compared, allocation and lifetime |
| `docs/research/measured.md` | Three-tier implementation × three sizes measured, item-by-item analysis, reproduction steps |
| `docs/research/raw/` | Raw logs: per-CASE output for the four build-option groups, bandwidth probes, disassembly excerpts, verbatim vectorization verdicts, accuracy comparisons |

> Note: the files in the table above are **overwritten** by every `make bench` / `make run`; they are not
> immutable archives.
> When citing them, cite the moment too: the file mtime is that round's timestamp, and the header line
> carries label/rev.
> On the same machine, the "distinguishable point count" over four runs was 90 / 92 / 93 / 107 (denominator
> always 304), with an all-point median ratio of 1.0063–1.0079; these numbers drift on every rerun.
> Conclusions that need long-term comparison should be copied into an analysis document, not taken directly
> from the numbers in raw.

### 8.2 Reproduction Commands

~~~bash
make test          # unit tests (30 assertions) + C compatibility test
make bench         # rerun four build-option groups × three tiers × three sizes; raw output to docs/research/raw/
make bench-matrix  # generate the comparison matrix from the raw logs
make bench-vec     # GCC auto-vectorization report (-fopt-info-vec-*)
make bench-checks  # accuracy / bit-identical checks (DIFF / IDENT)
make bench-asm     # disassembly archive
sudo sysctl -w kernel.perf_event_paranoid=1 && make bench-perf   # requires privileges; perf is currently disabled
~~~

### 8.3 My Own Adjudication Experiment

Three dot reduction paths (hadd every 8 elements / shuffle every 8 elements / 4 accumulators + final
reduction), **verify the sum before timing**:

~~~bash
g++ -O3 -march=native -std=c++17 /tmp/haddtest2.cpp -o /tmp/haddtest2 && taskset -c 2 /tmp/haddtest2
# all three paths give 2080.0 on the reference input (the correct value); timings are in §3.4
~~~

The source currently exists only in the session record. If you want it checked into `bench/` as a
permanent regression, tell me and I will add a make target.

---

*End of report. Every measurement in the four research notes comes with its original command and output;
this report only cross-confirms, adjudicates, and makes trade-offs among them.*
