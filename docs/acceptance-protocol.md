# flashvec Acceptance Protocol

> English | [简体中文](zh-CN/acceptance-protocol.md)

> This document answers one question: **what makes an optimization count as "accepted"?**
> It separates "looks fast" from "is fast", and sets out the arithmetic, reproducibility, and semantics criteria.
> Data sources: `bench/ab/` (same-process A/B harness), `verify/` (third-party bit-level comparison and structural
> checks), `docs/research/raw/` (raw logs).
> Only verified numbers are cited here; anything not yet measured is marked "TBD".

> **`docs/research/raw/` is a product that every `make run` overwrites; it is not an immutable archive.**
> Every re-run rewrites `ab_baseline_head.{txt,tsv}`, `_matrix.tsv`, `_smallN.tsv`, and `ab_aa_floor.tsv`.
> On the same machine, the measured "resolvable point count" was 90 / 92 / 93 / 107 across four runs
> (denominator 304), and the all-points median ratio ranged from 1.0063 to 1.0079.
> When you cite it, cite the moment with it: the header line carries label/rev, and the file mtime is that
> round's timestamp. Otherwise you end up quoting numbers from two different rounds side by side.
> Conclusions you need to compare over time belong in an analysis document under `docs/research/`; do not depend
> on the files in raw.

Test machine: AMD Ryzen 7 5700U, Zen 2, 8C/16T. L1d 32 KiB per core, L2 512 KiB per core,
**L3 4 MiB/CCX × 2 CCX = 8 MiB total on the chip, 4 MiB visible to a single thread**, 64 B lines.
AVX2+FMA, no AVX-512.

---

## 0. One-page summary: five conditions that must all hold

An optimization counts as accepted only if it satisfies all five conditions below **simultaneously**. "Optimization"
here means changes such as a new batch entry point, different instructions, a different reduction formulation, or an
added NT store.

| # | Criterion | What failing it means |
|---|---|---|
| 1 | **Bit-level comparison passes**: all forms × all tails × three compile settings; each of the three is defined in §1 | The numerical contract has changed. Element-wise entry points promise the same bit pattern; any single difference is a contract break, not "floating-point noise" |
| 2 | **Same-process median ratio ≥ 1.15x at the target size** | No measurable gain at the size where it was measured. "Same process" is the point: cross-process and cross-run differences on this machine reach 20–25% |
| 3 | **Same-direction reproduction**: swap the measurement order (AB→BA), or rebuild the binary and re-run; the direction agrees and is still ≥ 1.15x | The direction of one measurement may be noise or an artifact of code placement — alignment, a change in cross-TU inlining |
| 4 | The **A/A floor** at that (op, size) point is < the claimed gain; A/A floor = the ratio deviation of two pieces of **identical** code | The gain is below this machine's resolution for the same code, so it cannot be separated from measurement noise |
| 5 | **No other size regresses by more than 5%**; an expected regression must be **explicitly labelled**, e.g. the NT store at the L2/L3 sizes | Trading a win on one large array for a general slowdown on small arrays is a net-negative optimization |

The three parts of criterion 1 are expanded below; each is defined in §1.
"All forms" means the three forms: library-internal definition, header inline, and the pure-C fallback.
"All tails" means `n = 0..17` plus a number of sizes that are not multiples of 4 or 8.
"Three compile settings" means the default setting, `-O3 -march=native`, and
`-march=native -mfma -ffp-contract=fast`.

Beyond the five there is one methodological point, expanded in §4: **a single-point ratio cannot be trusted.** A
verdict must rest on the **per-size median ratio** and the **all-points median ratio**, never on one reading of one
(op, n).

The automated entry point is `make -C bench/ab verify OP=<op> N=<size>`. It runs everything in one pass and prints a
verdict per criterion; on FAIL it says which criterion failed and by how much. Details in `bench/ab/README.md` §11.

The companion structural check is `bash verify/checks.sh`. It checks header purity, dispatch placement, banned
builtin names, a binary-level instruction survey, and that retired internal names do not come back.

---

## 1. R1 Bit-level comparison: why "all forms × all tails × all three settings"

- **All forms**: one implementation must give one result in all three forms. The three forms are the library-internal
  `extern` definition (`FV_USE_EXTERN` / `FV_BUILD`), the header `static inline` (default), and the pure-C fallback
  (`-DFV_I_NO_INTRINSICS`, which goes scalar). A form difference changes inlining and instruction selection; FMA
  contraction shows up here most easily.
- **All tails**: the main loop is `n/4` and `n/8` (`n/2` for vec4), with the tail going scalar. `n = 0..17` covers
  every remainder. Also test the boundary where `n` is just above an integer multiple of the batch width, plus
  `n = 0` and `ptr = NULL`.
- **Three compile settings**: the default setting (`-O3 -fno-math-errno -ffp-contract=off`), `-O3 -march=native`,
  and the strictest, **`-march=native -mfma -ffp-contract=fast`** — the consumer bringing their own FMA flags.

  The header holds contraction-risky entry points down with `#pragma GCC optimize("fp-contract=off")`. This setting
  is the only evidence for that pragma: remove it and `vfmadd` appears silently.

The equality criterion: IEEE bit patterns via `memcmp`. The sole exception is when **both sides are NaN**, which are
judged equal and tallied separately, since the payloads may differ. A ±0.0 difference is reported as unequal, bit for
bit.

Every output buffer needs an out-of-bounds-write canary: nothing in `o[n..)` may be written. `_idx` variants must also
check that any `o[i]` not referenced by idx is left untouched.

**The single exemption: `dot_soa`, including the vec3/vec4 versions.** It uses 4 accumulators plus one final
reduction at the end to break the dependency chain, so its reduction order differs from the element-wise path; it is
therefore not judged at the bit level, only its error is tallied. Its **contract must be written down**, see §5.

---

## 2. R2/R3 Ratios and reproduction: why swap the order and swap the binary

- **Same process**: A and B must alternate inside one process (A/B/A/B…). Otherwise frequency migration, external
  load, and page-allocation differences leak into the ratio. This machine clocks between 1.64–3.97 GHz, and the same
  code can differ by 20–25% across runs.
- **Swapping the order (AB → BA)**: cancels systematic bias such as "whichever runs second benefits from a warm cache
  or a recovered clock".
- **Swapping the binary**: identical source, only code placement changes
  (`-falign-functions=64 -falign-loops=32 -falign-jumps=32`); semantics are unchanged. This one isolates one-off
  factors such as code placement. Measured, it moves the L1/L2 sizes by 5–10%, and that is not an optimization.

---

## 3. R4 A/A floor: how small a difference this machine can actually resolve

The A/A floor is the ratio deviation produced by **two copies of exactly the same code** alternating in one process at
one size; the ideal value is 1.000. It gives the lower bound on this machine's resolution for "the same code". Any
claim that a candidate is r times faster can only be true if **(r − 1) > floor_dev**.

One measured set on this machine, taken from **a single round** of `docs/research/raw/ab_aa_floor.tsv`: 304 large-N
measurement points, warm + cold:

- Distribution of `floor_dev = |ratio − 1|`: **29/304 points ≥ 15%**, **17/304 points ≥ 20%**.
- The worst point is 0.583, at `add_x4` @ 131072 in the warm setting. Restricting to the "resolvable" points, the
  maximum is 0.312.
- In other words, **a 15–20% floor is normal on this machine, and the extreme points are worse still**. An
  optimization claiming 1.15x cannot be told apart from noise at a point without an A/A floor.

> Those numbers are 29 / 17 / a worst point of 0.583. Like the resolvable-point count in §4, they are readings from
> **a single round** of raw, not constants. Re-run the same code on the same machine and the resolvable count drifts
> between 90 / 92 / 93 / 107 with the denominator fixed at 304; see §4 and the raw-drift note at the top.
> **The conclusion does not depend on which round it is**: "the floor is normally 15–20%, so a single point is
> undecidable" held in all four readings.

Two companion reading rules:

1. **Do not judge with points where `resolvable = 0`.** The paired IQR of that size is ≥ 5%; external load has
   already contaminated it.
2. **Prefer the L1/L2 sizes**: sampling coefficients of variation in the L3 and DRAM sizes routinely exceed 5% and can
   be contaminated by external load.

---

## 4. Why a single-point ratio cannot be trusted

This is a shared machine at loadavg 4–9 on 16 cores, and sampling coefficients of variation in the L3 and DRAM sizes
routinely exceed 5%. One measured set: about **30–35%** of points are resolvable, with 90 / 92 / 93 / 107 across four
runs and a fixed denominator of 304.

That number drifts from run to run, so cite the round when you cite it. "Resolvable" means paired IQR < 5.0%, and the
denominator 304 is the count of large-N measurement points. The same row also reports, separately, the group of cold
points that miss the bar; see the "sampling stability" line in `docs/research/raw/ab_baseline_head.txt`.

**The four readings are four `make run` invocations; raw keeps only the last round.** The specific numbers from
earlier rounds can only be read from the record in this table — which is exactly why citing raw requires the
timestamp.

So:

- Judge by the **per-size median ratio** — one ratio per size n, then the median — and by the **all-points median
  ratio**, the median over the whole sweep line.
- One reading of a single (op, n) can only **raise a question**; it cannot settle one.
- Reporting a number means reporting four things with it: working-set size, cache tier, compile flags, and whether it
  is resolvable (CoV / paired IQR). A bare "speedup" figure means nothing.

---

## 5. Semantic acceptance: split invariance

Test suites do not cover this class of property automatically, but callers always hit it: block-parallel processing,
streaming, accumulating chunk by chunk.

**Element-wise batch entry points have split invariance**: `add / sub / scale / lerp / madd / len / normalize` and
their `_x4` / `_soa` / `_mask` / `_idx` variants. Split the same array into any number of segments and call them
separately: the output is **bit-identical** to a single call.

**Reduction entry points do not**: for `fv_vec2_dot_soa` / `fv_vec3_dot_soa` / `fv_vec4_dot_soa` the reduction tree
changes with the split. One measured case on identical input:

| Case | One call over the whole range | Split into 2 segments, called separately, then summed | Relative difference |
|---|---|---|---|
| `fv_vec2_dot_soa`, n = 1000 | −64.4546509 | −64.4543991 | **≈3.9e-6** |
| `fv_vec2_dot_soa`, n = 8192 | 79.6643372 | 79.6624146 | **≈2.4e-5** |
| Control `fv_vec2_add_x4` (same corpus, same split) | — | — | **bit-identical** |

Across 12 corpora: `dot_soa` is **12/12 non-identical**, `add_x4` is **12/12 bit-identical**. The latter is the
control, and it shows the experiment is sound.

Conclusion and requirements:

- This is a **design consequence of reduction, not a defect**: changing the reduction order always changes rounding.
- But it **must be written into the contract**, or callers will take it for floating-point noise. Callers who need
  "the result is independent of the split" must use the scalar `fv_vec2_dot` with sequential accumulation, or fix
  the split themselves.
- The magnitude is not small: 2.4e-5 is roughly 400 times the "relative error ≤ 4e-8 for n ≤ 1000" figure. The two
  statements are about **two different things**. "The difference from the scalar path" is the error number in the
  header; "the difference from another split" is this section. The former is not an upper bound on the latter.
- Therefore, **when accepting an optimization, the bit-level exemption for `dot_soa` is not a precision exemption.**
  Its error must be reported as "corpus condition + absolute error", never as a bare relative error, and the fact
  that "a split call will change it" must be stated as well.

---

## 6. Byte accounting: nominal traffic vs real traffic including RFO

Computing GB/s from **nominal** bytes (read + write) misreads the NT-store gain of `_stream`. On a write miss an
ordinary store must **first read in** the entire cache line; that step is RFO, read-for-ownership. Real traffic is
"one write plus one read".

| Entry point | Nominal B/element | Real B/element | Notes |
|---|---|---|---|
| `add_x4` / `sub_x4` / `lerp_x4` / `madd_x4` / `add_soa` | 24 (16 read + 8 write) | **≈32** (+8 RFO) | the output array is read one extra time |
| `scale_x4` / `scale_soa` | 16 (8 read + 8 write) | ≈24 (+8 RFO) | as above, one read stream fewer |
| `dot_soa` (vec2/vec3/vec4) | **16** (four `float` arrays ax/ay/bx/by × 4 B) | 16 (pure read; no write, no RFO) | **not 32 B**: only 4 floats are read per element |
| `normalize_soa` | **16** (8 read + 8 write) | ≈24 (+8 RFO) | **not 16+16** |
| `*_stream` (NT store) | 24 (same as the ordinary version) | **24** (no RFO) | the real gain is **25% fewer bytes moved**, not faster bytes |

Two criteria:

1. For the same operator, the **real** GB/s of the ordinary version and of the `_stream` version should be roughly
   equal — both sit at the DRAM bandwidth limit — while the **nominal** GB/s differ by 1.33 times. If both nominal
   and real differ by 1.33 times, the gain is not coming from traffic. This is the ruler for spotting false
   positives.
2. NT stores only deliver that 25% of bytes when the **working set is far beyond L3 and the write stream is a
   clean, one-shot write**. When the working set lands in L2/L3, they are **2–8 times** slower than the ordinary
   version: measured inside L2 at **0.13–0.51x**, with all 6 entry points judged FAIL (slower), see §9 item 2.

This is why the **knee** must be recomputed with real traffic. For `add`, nominal 24 B/element gives
4 MiB ÷ 24 B ≈ 1.75×10⁵ elements. On real 32 B/element the knee comes earlier still, at about 1.3×10⁵.

---

## 7. Cache accounting (a precondition for reporting numbers)

L3 on this machine is **4 MiB/CCX × 2 CCX = 8 MiB total on the chip**.
`/sys/devices/system/cpu/cpu0/cache/index3` reports `level=3, size=4096K, shared_cpu_list=0-7`.
**A single thread only ever gets its own 4 MiB.**

- Saying "L3 is 8 MiB shared" is wrong: 8 MiB is the total on the chip, not what one thread can use.
- Compute the knee from 4 MiB: `1.75×10⁵` active elements @ 24 B/element. Older documentation used 8 MiB and got
  `3×10⁵`, twice too high.
- There is direct measurement: the same-process A/B harness gives **38–57 GB/s at 2^17 (working set 3 MB)** and
  **5.6–8.8 GB/s at 2^18 (6 MB)** — the cliff sits exactly on this line.
- This repository's consistency check prints the measured `index3` value every time, in section [7] of
  `bash verify/checks.sh`, so that documentation and machine cannot drift apart again.

---

## 8. Where no trustworthy data exists (do not cite it as numbers)

- Speedup **relative to scalar** for the three-component SoA batch (`fv_vec3_*_soa`) and the four-component
  `_x2` / `_soa` batch: **TBD**.

  This round measured only their **form difference**, i.e. library-internal definition vs header inline, over
  13 entry points × 5 sizes = 65 points: 5 PASS / 8 FAIL (slower) / 52 INCONCLUSIVE, readings 0.97–1.93.
  Those numbers may only be read as form, never cited as speedup; see §9 item 5.
- `AoS↔SoA pack/unpack cost`: still unmeasured, and it directly decides whether the SoA entry points are worth it.
- The scalar fallback path on non-AVX2 hardware: this machine has AVX2, so `__builtin_cpu_supports("avx2") == 0`
  cannot be forced. Only structural evidence exists, namely section [9] of `verify/checks.sh`: the kernel is still
  present in the baseline setting, and dispatch goes through the avx2 bit of `__cpu_model`. It has never run on a
  real non-AVX2 machine.

---

## 9. Dead ends (recorded so that nobody redoes them)

This section records directions **already rejected by measurement**. It is not "not tried yet" but "tried, and the
answer was no". Read it before redoing anything, and give the physical reason why "this time is different". If you
cannot give one, do not re-run; a re-run only buys another noisy reading. Every item carries the location of its
evidence.

1. **The in-place write of `_mask_inplace` brings no throughput gain; every size is INCONCLUSIVE.** 5 densities
   (1% / 5% / 20% / 50% / 95%) × 4 sizes = 20 points, gains **-0.10…+0.00**, same-round A/A floor **1.06–1.56**.

   Its 4 size points in the campaign verdict table are all INCONCLUSIVE too: -0.41 gain and 1.48 floor at the 4M
   size. **The reason lies in the byte accounting.** The target array of an in-place write already has to be read
   before it is written, so the corresponding cache lines are in cache, and an ordinary store **produces no RFO** in
   that situation.

   In-place saves **0 B**, not 8 B/element. Its value is therefore **form** — one output buffer fewer, one
   allocation fewer — not speed; do not use it to buy throughput.
   Evidence: `docs/research/raw/ab_campaign_verdict.txt`, `bench/ab/README.md` §10.8.

2. **NT stores inside L2/L3 run at 0.13–0.51x, 2–8 times slower.**

   The 6 `_stream` entry points measured at n=8192, where the working set is 128–192 KiB and lands inside L2:
   `add` 0.217 / `sub` 0.301 / `scale` 0.130 / `add_soa` 0.229 / `scale_soa` 0.138 / `normalize_soa` 0.506 —
   **6/6 judged FAIL (slower)**. Only a working set **far beyond L3** with a clean, one-shot write stream gets that
   25% of bytes, see §6.

   And the choice must be made **explicitly by the caller**; the library never switches on its own. Nor is it
   "always effective beyond L3": the 1M/8M sizes PASS (1.29–1.58x), the 4M size is 6/6 INCONCLUSIVE (floor > 15%),
   and n=131072 is a blacklisted size. Evidence: as in §1, and the `_stream` contract in `include/flashvec.h`.

3. **VMASKMOVPS (`maskload` / `maskstore`) is a microcode-level slow instruction on Zen 2; this path has been
   removed from the library.** The two internal names `fv_i_maskload_add4_avx2` / `fv_i_maskstore_add4_avx2` are
   marked **retired**.

   Rule [11] of `bash verify/checks.sh` requires **0 occurrences** across the whole repository; reviving them is
   caught by the gate on the spot. The reason that has nothing to do with numbers matters more. With them, "whether
   a missed lane is physically written" becomes a property that **varies with mask density**, and no caller can
   depend on a store guarantee that shifts with the input distribution.

   Early probes put them at 1.2–3.8 times slower than a plain `storeu` at the same density; this round **did not
   re-measure**, so read that as an order of magnitude only. Evidence: `verify/checks.sh` [11], the
   `include/flashvec.h` contract.

4. **Split invariance does not hold for reduction entry points.** This is a design consequence, not a defect, but it
   **must be written into the contract**: for `fv_vec2_dot_soa` / `fv_vec3_dot_soa` / `fv_vec4_dot_soa` the
   reduction tree changes with the split.

   At n=8192 the relative difference between "the whole range" and "two segments summed" is **2.4e-5**; at n=1000
   it is about **3.9e-6**; 12/12 are non-identical. The control `add_x4` is 12/12 bit-identical on the same corpus.
   Callers who need "the result is independent of the split" must use the scalar `fv_vec2_dot` with sequential
   accumulation.

   Evidence: the table in §5, and the contract in `include/flashvec.h`.

5. **The "inline form" is not faster in itself.** Moving an implementation from "library-internal definition" to
   "header inline" changes only code placement; the algorithm does not change at all.
   This round measured 13 three- and four-component entry points × 5 sizes = 65 points in intra mode, where intra
   mode means A and B differ in this one variable. The result: **5 PASS / 8 FAIL (slower) / 52 INCONCLUSIVE**.

   In the smallest size alone, n=8192, readings for **the same source** run the whole way from **0.97 to 1.93**.
   Reading 0.97 as 3% slower shows that "inlining can also read as slower"; reading 1.93 as 93% faster shows that
   "changing only the form can read out an almost 2-fold false gain".

   The other face of the same phenomenon is the same-round A/A floor: for that **identical** code at the 4M size,
   the ratio across the two sides reaches **1.34–1.61** — changing only the object-file layout reads as roughly
   1.5-fold. So any "it is faster because it is inlined" conclusion must clear R3 (same-direction reproduction) and
   R4 (same-round A/A floor). A single-point ratio can only raise a question; it cannot serve as evidence.

   Evidence: the intra rows of `docs/research/raw/ab_campaign_verdict.txt`.

6. **`_idx` software prefetch does nothing.** The measured configuration was three index distributions × three
   sizes {65536, 1M, 8M} × three prefetch distances {4, 8, 16} × two hints {T0, NTA}.
   **Random distribution: 0 PASS / 30 points, all INCONCLUSIVE**, from the harness, full table of 30 rows.

   **All three distributions combined: 0 PASS / 2 resolvable regressions / 52 INCONCLUSIVE**, from a purpose-built
   same-process A/B driver, 54 rows. R1 bit-level comparison passed 8/8: macro on × four forms × two hints, with
   the 175 296-byte main output stream byte-for-byte identical. **Prefetch changed no semantics; it simply has no
   performance.**

   Resolvable points mean paired IQR < 5% and same-round floor ≤ 15%; the actual readings at those points follow.
   - The 8M **random** size is the most latency-limited of the sizes: gains **+0.3% … +1.9%**, and -1.8% as well;
     same-round floor 1.5%–4.6%. Nothing resolvable separates D=4/8/16 or T0/NTA.
   - 65536 **ordered ascending** is **0.606x**, **blocked clustering** is **0.700x**, the latter grouping every
     64 indices into one 4 KB window. That is a 39% / 30% regression — 6–8 times the R5 criterion (anything above
     5% is rejected), consistent in direction under R3, and coming from two independent builds, D=8 and D=16.

     `_MM_HINT_NTA` is worse still at 65536: **0.12–0.41x**.
   - The 65536 random size looks like it has a 1.11–1.24 "gain", but the **same-round A/A floor is itself
     1.48–1.93**. Across three repeats the ratios are 1.115 / 1.150 / 1.157, with floor p90 deviations of
     0.48 / 0.93 / 0.57. That size is one of the hardest points to judge on this machine; under R4 it must be
     judged INCONCLUSIVE.

   **The reason is read out of the data**:
   - Under a random distribution the bottleneck is not "not enough lookahead" but **the parallelism of outstanding
     misses**. At the 8M size, 55–60 ns per element, i.e. 20 B/element ≈ 0.34 GB/s ≈ one DRAM latency ÷ 1.5
     concurrent accesses. Prefetch requests compete for the same miss slots as demand misses, that is, the L1 fill
     buffer.

     It issues the same accesses **earlier** without adding parallelism, so the gain is zero. That also explains why
     sweeping distance and hint turns up nothing: changing the lookahead does not change the number of slots.
   - Under the ordered and blocked distributions the working set is in L1/L2, which the hardware prefetcher already
     covers. There the loop is **issue-limited**, not latency-limited, and 1–2 extra prefetch µops per element
     convert straight into a 25–40% throughput loss. NTA is worse: the "non-temporal" path is harmful for lines
     that are **used repeatedly**, which stay in L1 yet get demoted or pushed out.
   - So the physical ceiling of this path is the stretch that turns latency-bound into MLP-bound, and that stretch
     has already been taken by the out-of-order window: measured MLP ≈ 1.5–2. Going further requires changing the
     **data layout** — batching or sorting the indices first, then walking them sequentially — which is
     preprocessing that belongs to the gather performed by the caller. That is not something a batch entry point can
     do.

   Evidence and reproduction. The harness is `bench/ab`: A comes from the same tree with `FV_I_IDX_PREFETCH=0`, B
   from the same tree with `=1`. Compile `var_a.o` and `var_b.o` from the same `ab_variant.cpp` using
   `bench/ab/mkrename.sh` and reuse the driver, with `--aa-each-round`.

   The data is in three files:
   - `bench/ab/build/idxpf/sweep_random.tsv`, 30 points.
   - `bench/ab/build/idxpf/chk_{1,2,3}.tsv`, three repeats of the 65536 random size.
   - `pf_D8T0_fv_vec2_{add,scale}_x4_idx.tsv`, the pair from the coarse sweep.

   The one-off driver behind the 54 points of the three-distribution run was deleted together with the
   implementation. That driver was two `tests/idxpf_side.c` wrapper TUs plus `tests/idxpf_dist.c`, running an ABBA
   scheme with a homogeneous A/A on every round.

   The data remains in `bench/ab/build/idxpf/dist/dist_sweep.tsv` and `dist/run_{8_0,16_0,8_1}.log`.

   The implementation was **deleted wholesale** as well: the three macros `FV_I_IDX_PREFETCH` / `FV_I_IDX_PF_DIST` /
   `FV_I_IDX_PF_HINT`, plus two prefetch branches in the `_idx` loops. The header is back to its original state;
   `grep -c IDX_PREFETCH include/flashvec.h` = 0.

   Re-checked with `make test-header` / `make test-inline` / `bash verify/checks.sh`: 11 rules, FAIL=0.

7. **`-ffast-math` buys nothing, costs accuracy, and makes one entry point measurably slower.**

   The library is already at the hardware limit for its workloads, so there is nothing left for fast-math to
   exploit. Measured here with two shared objects loaded side by side into one process and run alternately,
   median of 7, 5 entry points × 4 sizes:

   | Entry point | size | default | `-ffast-math` | ratio |
   |---|---|---|---|---|
   | `add_x4`, `scale_x4`, `dot_soa`, `len_soa` | all | — | — | 1.00 ± 0.05, i.e. noise |
   | `normalize_soa` | 8192 | 0.35 ns/e | 0.42 ns/e | **0.83–0.85x**, stable over four runs |
   | every entry point | 8M | — | — | ≈1.00, DRAM-bound |

   The instruction census explains it. fast-math replaces every `vdivps` with `vrcpps` plus a Newton step:

       default       vdivps=5  vrcpps=0  vrsqrtps=0
       -ffast-math   vdivps=0  vrcpps=5  vrsqrtps=2
       both          vfmadd=0                    (the pragma still holds)

   `normalize` runs 8 independent lanes, so `vdivps` throughput is already pipelined away. Trading one divide
   for three or four instructions spends issue slots without shortening the dependency chain.

   Accuracy, on 2^20 random vectors:

   | Entry point | elements with a different bit pattern | max relative error |
   |---|---|---|
   | `len_soa` | 0 / 1048576 | 0 |
   | `normalize_soa` | 340525 / 1048576 (32.5%) | 2.8e-7 |
   | `dot_soa` | — | 1.1e-7 |

   Reproduction:

       c++ -std=c++17 -O3 -fno-math-errno -ffp-contract=off -fPIC -shared -Iinclude \
           src/flashvec.cpp src/flashvec_batch.cpp src/flashvec_vec4.cpp -o /tmp/libbase.so
       c++ -std=c++17 -O3 -ffast-math               -fPIC -shared -Iinclude \
           src/flashvec.cpp src/flashvec_batch.cpp src/flashvec_vec4.cpp -o /tmp/libfast.so
       objdump -d /tmp/libbase.so | grep -c vrcpps   # 0
       objdump -d /tmp/libfast.so | grep -c vrcpps   # 5

   Then `dlopen` both and time the same entry point alternately. This is the exact scenario a consumer hits:
   because the batch entry points are `static inline`, defining `-ffast-math` compiles the header's
   implementation under it.

**How to use this**: if a new proposal is structurally the same as one of the items above, first say which layer's
assumption changed — byte accounting, microarchitecture, or semantics. If you cannot say, go back to the five criteria
in §0 and work the "build the harness first, then touch code" way.
The negative result in `bench/ab` (§1) was caught exactly by building the harness first.

---

## 10. Known limitations (which way the current criteria are biased; unfixed)

§9 records "approaches that were tried and rejected"; this section records **problems with the criteria themselves**.
They have not been fixed, so each one limits what a PASS above actually proves. The criteria are implemented per
`bench/ab/README.md`: §10.7 the per-size criterion, §11 the acceptance criteria, §13 the six holes the author lists.

- **The per-size criterion is loose (known, unfixed).** The self-check's per-size rule compares "the median ratio of
  all points in that size" against "the p90 of the single-point A/A floor", i.e.
  `|size median ratio − 1| ≤ that size's floor p90`. The standard error of a median over 15 points is about 1/√15
  of a single point's, i.e. ≈1/4.

  The comment in `bench/ab/ab_bench.cpp` says exactly that: "the median of 15 samples has a standard error of about
  1/4 of the single-point noise". So the rule lets through deviations larger than the true per-size noise: it
  measures a number that has already averaged most of the noise away against a ruler graduated from a single
  reading, and so it is necessarily loose. Why it was not tightened by √N: tightening needs an extra assumption
  about **independence between points**.

  The noise of 15 operators is not independent and identically distributed; they share the same stretch of frequency
  drift and the same external load. And the per-size rule is only a **health check**, not evidence — the real
  criteria are R2/R3/R4 point by point. **Do not take "per-size PASS" as evidence that "this optimization works".**
- **The A/A floor only measures noise on the A side**, i.e. A vs A, while the noise in an A/B ratio also contains the
  B side. When the two sides differ in noise, the floor underestimates jitter. That is, the floor is a **lower
  bound**: a small floor does not mean the point is decidable, only that "A against itself" is resolvable.
- **The INTRA-mode floor is struck against a reference op**, i.e. ref vs ref, not new vs new, which raises the
  problem of representing one op's noise with another's. A new entry point and the reference op may differ in memory
  access pattern, tail length, and register pressure, so the floor's magnitude cannot be assumed to carry over.
- **The per-size floor itself doubles between rounds.** Measured at the 65536 size: 0.130 / 0.169 / 0.218. The same
  size can flip between "decidable" and "INCONCLUSIVE".

  In the round that the current raw left behind, the n=65536 warm floor is **0.218**, which has already tripped the
  "floor > 20% → this size is INCONCLUSIVE" fallback, while the same size was only 0.130 / 0.169 in the other two
  rounds. So only "**a per-point PASS whose floor is small**" counts as a hard conclusion.

  A "whole-size PASS" may become INCONCLUSIVE the next time you re-run; it describes one measurement, not a
  property.

One sentence to close: all four items above say the same thing. **A PASS at the size or aggregate level is the
harness doing a health check; only the per-point R2/R3/R4 is evidence.** Any inference of the form "the whole size
passed, so this optimization works" lies outside that boundary.
