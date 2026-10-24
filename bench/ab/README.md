# bench/ab — In-Process Interleaved A/B Performance and Throughput Harness

> English | [简体中文](zh-CN/README.md)

This harness answers one question and nothing else: **how much faster is the candidate implementation than the baseline**.

The answer has to be **falsifiable** on this machine. The machine shifts frequency, resources are shared, and other load eats 30-50% of it.

The one-line verdict: **what it can measure, it labels; what it cannot measure, it refuses to conclude.**

> **The resolution floor of this harness is a hard constraint**: of 304 large-N points, only **85-107** are resolvable, roughly 28-35%.
>
> That statistic **drifts from run to run by itself**. Six separate rounds read **85 / 90 / 92 / 93 / 95 / 107**, with a constant denominator of 304 and a spread of about 26%.
>
> So a citation has to carry its round. Writing it as a range is not hedging; it is the only correct way to write it.
>
> Most recent measurement: 182/182 gates pass, all-point median ratio 1.0109, 95/304 resolvable.
>
> On individual (op, size) points, **two identical copies of the code** sit a stable 15-20% apart. That is the A/A floor; see §10.
>
> So "median >= 1.15x within one run" is **not enough on its own** to call a result real. The criteria are in §11.

---

## 1. Why "run two binaries and compare" does not work

Item 3 under "Performance" in the README states the bound: the same code varies 20-25% across runs, with frequency, neighbor load, and page mapping all changing.

So there is exactly one way to compare: **in one process, over one data set, inside one time window, alternate A/B/A/B…**, then divide per-tier medians.

That one requirement is the entire design premise of ab_bench.

## 2. How two implementations coexist in one process

| Role | Source | How it is compiled |
|---|---|---|
| baseline | snapshot/, the snapshot captured by git show HEAD: | -include the generated build/rename_a.h; entry points become ab_a_* |
| candidate | by default the same snapshot, for the self-check; or any source tree | -include build/rename_b.h; they become ab_b_* |
| scalar anchor | src/flashvec.cpp from the HEAD snapshot | not renamed; the loop calls it from **another TU** |

snapshot/ consists of include/flashvec.h, src/flashvec_batch.cpp, and src/flashvec.cpp as captured by git show HEAD:.

To point the candidate at a different source tree, use CAND_HDR / CAND_SRC / CAND_INC. The library has only one scalar implementation, so the anchor needs no renaming.

**Why not compile the working-tree src/ directly**: another agent is editing it. The baseline always comes from git show HEAD:.

make snapshot leaves mtime alone when HEAD moved but the content did not, so it never costs an extra rebuild.

**Why the timing loop sits inside each variant TU**, instead of the driver calling through a function pointer:

include/flashvec.h is being refactored into the header-inline form, static inline by default.

In the inline form the batch entry points are TU-private static and nothing outside can link to them; only the library form gives them external symbols.

With the loop inside the variant TU, both forms are measured as they actually are:

* Library form: the loop contains one cross-function call — exactly what a caller linking the static library gets.
* Inline form: the body is inlined into the loop — exactly what the refactor is aiming for.

The driver knows only three symbol families, ab_X_loop / ab_X_abi / ab_X_opflags, all visible to nm.

It has no idea what the batch entry points are called or how many exist, so **no header refactor can reach the harness**.

Switching form takes only make … CAND_WITH_SRC=0. 0 means the batch entry points live in the header and no .cpp is compiled.

**The symbol rename table is generated**, by mkrename.sh. It pulls every name of the form `fv_xxx(` out of **this side's header** and emits `#define fv_xxx ab_[ab]_xxx`.

The payoff: any entry point added to the header is isolated into its own namespace automatically, with **no hand-maintained rename table**. The old hand-written ab_rename_a/b.h files are gone.

**Pitfalls already hit.** All four were hit for real; they are written down here to save the next person two hours:

1. A -D macro rename must be force-included with -include **before the first line** of the source file. Otherwise only one side gets renamed — the header declarations or the .cpp definitions — and you find out at link time.
2. Putting a function declared `extern "C"` inside a namespace does not rename it; the symbol keeps its C name. So "isolate the two implementations with namespaces" is a dead end.
3. The two output buffers must be **allocated back to back, exactly n*8 bytes apart**, where n*8 is a multiple of 4096. Otherwise A and B land in different store/load 4K aliasing relationships, and the L1/L2-resident tier differs by 20% out of nowhere — a measured value.
4. objcopy --prefix-symbols renames undefined references such as sqrtf / memcpy along with everything else, and the link then fails. The only way left is per-symbol --redefine-sym, which is more trouble than the macro rename and easy to miss cases in.

## 3. Three gates (fail one, get no numbers)

1. **Bit-for-bit comparison against an independent reference implementation.** The reference is hand-written in ab_bench.cpp: pure scalar, pure sequential, calling no batch entry point at all.
   The comparison is **memcmp**, covering the mask's "write back to a[i]", the positions idx does not reference, and pipeline intermediates.
   The only exception is fv_vec2_dot_soa. The header explicitly does not promise bit-identical results, so equality is judged by |diff| <= max(1e-5*|result|, 1e-6*Σ|termwise product|).
   **This bound is independent of n and of the corpus.** A bare relative error gets inflated when the result nearly cancels; at n=2^18 it already reaches 7.4e-5, which is why the first version produced a false report.
2. **The outputs of A and B must also be bit-identical**, same exception as above. Compared once more after timing ends.
   In-place operators are excluded: they accumulate repeatedly over the same buffer, so the bit-level verdict comes from check_op before timing.
3. **The input buffer checksums are unchanged before and after timing**, which catches a timed implementation that corrupts a/b/mask/idx.
   An in-place operator writes to its own first input. The driver passes **a copy of it** (O.o) as the writable array, so this gate holds for in-place operators too.
4. **--selfcheck**: the baseline-vs-baseline ratio must land inside **[0.85, 1.15]**. If it does not, the run raises a loud alarm and exits with code 3.
   **This matters more than any speedup number.**

Any hard failure prints **the function name and the first mismatching element**, prints **no timing or throughput result at all**, and exits with code 4.

make negctl uses exactly that as its assertion; see §7.

--gate-only runs the gates only, emits no timing, and exits 0/4. verify uses it for the all-tail-length cross-check.

## 4. Three measurement regimes and their accounting

### Large-N throughput sweep (2^12 … 2^24 by default, 8 tiers)

* Reports three sets: **ns/element**, **nominal GB/s**, and **real GB/s**. Nominal GB/s = (bytes read + bytes written) / time; real GB/s follows RFO accounting, see §12.
* **B/element is checked by hand, per function**, and is not derived as "array size × constant". The table lives in the registry ab_ops.h and the driver only reads it; the table layout is in §9.
  Look up the actual values with `./build/ab_bench --list-ops`.
* **Upper-bound references** (same size, same timer): memcpy (8R+8W), memset (0R+8W), pure-read sum (8R), triad (16R+8W, the same form as the AoS batch entry point).
  Batch functions are reported **as a percentage of memcpy throughput at the same size**. Read-heavy operators (dot / pure read) **can exceed 100%**, because read bandwidth is inherently higher than copy bandwidth. That is not an error.
* **Both warm and cold are reported**. Warm runs repeatedly over buffers already resident; cold is a single call after madvise(MADV_DONTNEED), **including minor page faults**, with the pages zeroed by the kernel.
  The idx array is never dropped. Once it is dropped, every scattered write lands on index 0, random writes turn into hot-spot writes, and the numbers are fake.
* **Cache levels are labelled from measurement**: L1d 32 KB, L2 512 KB, **L3 4 MB/CCX**. The evidence is size=4096K and shared_cpu_list=0-7 under /sys/.../index3.
  The "8 MB L3" in the project docs is wrong; one thread only ever gets its own 4 MB. In the sweep line, 2^17 (3 MB, in-band) and 2^18 (6 MB, out-of-band) sit exactly on either side of that boundary.

### Small-N call-overhead regime (n = 4/8/16/32/64 by default)

* Reports both **ns/call** and **ns/element**. The bottleneck in this regime is per-call overhead, and reporting ns/element alone dilutes it out of sight.
* Which operators are covered is decided by the registry's `AB_OPF_SMALL` flag, not by hard-coded indices.
* The control is scalar_loop: call the scalar API element by element, with **the loop and the implementation in two separate TUs** (ab_scalar_loop.cpp / ab_scalar_lib.cpp).
  Only then can GCC neither inline nor vectorize, which is the true form of "a caller linking a precompiled static library".
* Emits the A/scalar multiple, answering directly whether the batch entry point is still worth it at small N.

### Statistics

* **At least 5 samples** per point, taking the median, and continuing until both sides' coefficient of variation is < 5% (--cov) or the --trials-max cap is hit. The actual sample count is written into the TSV.
* **Primary criterion is ratio = median(per-round A/B)**: A and B are paired inside the same time window, so low-frequency drift cancels out.
  ratio_medAB = median(A)/median(B) is reported alongside as a cross-check, together with pair[min,max], CoV, and min/max.
* By default the order is **fixed to A then B** (--order=ab). Measured at 2^20 with 9 samples: ab → ratio 0.995 / CoV 9.7%, ba → 1.030 / 34.5%, alt → 1.305 / 24.3%.
  A fixed order removes the confounder of "each pair sitting in two different entry states". Both ab and ba land near 1.0, which says the order itself carries no systematic bias. Use --order=ba to cross-check.
* **The default clock is this thread's CPU time** (CLOCK_THREAD_CPUTIME_ID), not wall time: time preempted by other processes is not charged to the measured code.
  It still includes the real stalls caused by neighbors stealing L3. --clock=mono switches back to wall time if you want to check for yourself.

## 5. Outputs

    docs/research/raw/ab_baseline_head.txt         full console output (matrix + RFO accounting + self-check)
    docs/research/raw/ab_baseline_head.tsv         long table: one row per (phase,op,n), 47 columns
                                                   (32 base + 8 RFO accounting + 5 same-round A/A + 2 mode/ref_op)
    docs/research/raw/ab_baseline_head_matrix.tsv  matrix: phase × metric × op × size
                                                   (metric: gbs / real_gbs / pct_memcpy / ratio)
    docs/research/raw/ab_baseline_head_smallN.tsv  small-N regime: n × op
    docs/research/raw/ab_aa_floor.tsv              A/A floor table (see §10)
    docs/research/raw/ab_campaign_verdict.txt      verdict table for one round across all new entry points (campaign.sh)

make -C bench/ab run-worktree swaps the candidate for the current working tree and writes to ab_vs_worktree.*.

> **The raw directory is an artifact that every `make run` overwrites; it is not an immutable archive.**
>
> Every rerun rewrites `ab_baseline_head.{txt,tsv}`, `_matrix.tsv`, `_smallN.tsv`, and `ab_aa_floor.tsv`.
>
> The numbers the docs quote (resolvable point count, per-tier median ratios, one point's ratio) drift along with it.
>
> On the same machine, the "resolvable point count" measured **90 / 92 / 93 / 107** across four runs, all with a denominator of 304; the all-point median ratio moved between 1.0063 / 1.0070 / 1.0072 / 1.0079.
>
> Quoting these numbers **requires quoting the moment with them**. The header of `ab_baseline_head.txt` carries label/rev, and the file mtime is that round's timestamp. Otherwise you end up quoting numbers from two different rounds in the same breath.
>
> Conclusions you want to compare over time belong in an analysis document under `docs/research/`; do not cite raw directly.

## 6. Command reference

    make -C bench/ab run               self-check + baseline snapshot (8 tiers + small N + A/A floor table)
    make -C bench/ab run-worktree      candidate = current working tree
    make -C bench/ab verify OP=… N=…   five verdict rules from same-direction evidence (see §11)
    make -C bench/ab aa-floor          rerun the self-check alone and land the floor table
    sh bench/ab/campaign.sh            one round over _stream/_mask_inplace/vec3/vec4 (own quiet-window check)
    make bench-ab / bench-verify …     top-level Makefile forwarding targets (additions only, existing ones unchanged)
    make -C bench/ab negctl            negative control
    make -C bench/ab build             compile only
    make -C bench/ab clean             delete bench/ab/build only

Knobs in common use, all overridable from the command line:

    SIZES=4096,65536,1048576      large-N sweep line
    SMALL=8,32                    small-N values
    TRIALS_MIN=5 TRIALS_MAX=12    adaptive sampling
    TARGET_MS=8 SMALL_MS=20       target duration per sample
    CPU=6                         core pinning (cpu0-7 share one CCX)
    CAND_INC / CAND_HDR / CAND_SRC / CAND_WITH_SRC   swap the candidate source tree and its form
    AB_OPS_EXTRA=1                enable the experimental area (_stream / _mask_inplace)
    SELFCHECK=0                   disable the [0.85,1.15] band check

## 7. Negative control (proof that the harness bites)

    make -C bench/ab negctl

Copies a **private** tree from the HEAD snapshot, changes one kernel, and compiles it as the candidate. Three assertions: exit code 4, the output names the function, and no timing result appears at all.

Three cases:

1. The mask select direction of _mm256_blendv_ps inverted → must name fv_vec2_add_x4_mask
2. mul changed to add inside fv_scale4_avx2 → must name fv_vec2_scale_x4
3. The _CMP_NLE_UQ comparison in normalize_soa negated → must name fv_vec2_normalize_soa.
   This one is caught **only because** every 512 elements of the data hides an exact zero vector. The zero-vector probe is deliberate; do not remove it.

src/ and include/ in the repository are never touched, so there is no "restore afterwards" step: the working tree is never modified.

## 8. Known limitations (stated honestly)

* **This machine is permanently occupied by external load**, with loadavg sitting at 4-9 on 16 cores. Sample CoV in the L3 and DRAM tiers routinely exceeds 5%.
  Such points are tagged CoV! and counted as unresolvable: **the numbers are still reported but never enter the criteria**. Self-check conclusions hold only for the resolvable subset, so if you plan to use the 1.15x criterion, look at the L1/L2 tiers first.
* A cold measurement zeroes the pages, so the warm timing that follows runs on zeroed data. This does not affect this library's batch entry points, which have no data-dependent branches and a tail length of 0.
  If a kernel with data-dependent branches is ever added, the data has to be refilled.
* At small sizes the pure-read probe is load-port / ALU bound, not bandwidth bound. It only represents bandwidth in the DRAM tier.
* Frequency is not pinned: governor=performance, and 1.5-3.9 GHz is observed, floating with load. **Absolute values are not comparable across sessions**; only ratios within one run mean anything.
* **Code placement**: the two forms of the same header (defined in the library / inlined in the header) change inlining and alignment by themselves, shifting the L1/L2 tiers by 5-10%.
  This is not "optimization". The criteria turn such one-off factors off by combining two independent builds (verify) with the A/A floor.
* The three-dimensional `fv_vec3_*_soa` and four-dimensional `fv_vec4_*_x2` entry points **are not in the registry yet**.
  Their element width is 12 B / 16 B, unlike the 8 B of fv_vec2. The n semantics and allocation of Buf have to change with them; see the "new layout" item in §9.
  The driver itself needs no change, and that is what table-driven buys.

---

## 9. The op registry: how to add an op

The single source of truth is `bench/ab/ab_ops.h`. **Adding an op touches two places only**:

1. add one line inside `AB_LIB_OPS(X)`, or `AB_LIB_OPS_EXTRA(X)` for the experimental area;
2. add one `abt_<ID>` thunk in `ab_variant.cpp` — that is the snippet that makes the call.

Everything else is derived automatically, and **not one line of driver logic changes**:

| Derived artifact | Expanded by |
|---|---|
| `AB_OP_<ID>` enum (both side TUs and the driver share one numbering) | ab_abi.h |
| `fv_xxx -> ab_a_xxx / ab_b_xxx` rename table | Makefile + mkrename.sh (read from the **header**) |
| OPS table (name / byte count / check policy / reference implementation / flags / RFO class) | ab_bench.cpp |
| variant-TU dispatch table + `ab_X_opflags` (NT or not) | ab_variant.cpp |
| checksum gate, small-N coverage, matrix, TSV, RFO accounting, verify | all read the table |

The row format has 10 columns. A NOTE must contain no comma; external scripts parse by column:

```c
X(ADD_X4, "fv_vec2_add_x4", VEC2, 8, EXACT, ref_add_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,
  "a/b read 16; o write 8 triggers an 8 B RFO read → real 32 B/e")
// ID    enum name AB_OP_<ID> (also fixes the thunk name abt_<ID>)
// NAME  the name used in reports and the TSV
// KIND  output form: VEC2 / SOA2 / SOA3 / SOA4 / SOA4V / VEC4 / FLOAT / NONE
// ELOP  **element width** (bytes): vec2 8 / vec3 12 / vec4 16. Within one tier sweep every op moves
//       the same byte count (element count = n*8/ELOP), so ratios across op widths are already byte-normalized.
// POL   EXACT (bit-identical) / REL (tolerance for dot_soa only) / NOCHECK (upper-bound reference)
// REF   independent scalar reference implementation inside the driver
// FLAGS AB_OPF_PIPE / SMALL / INPLACE / NO_BASELINE (absent from HEAD); upper-bound references use AB_MF_*
// RD/WR nominal read/write bytes per element (hand-checked per function)
// RFO   AB_RFO_MISS / HIT / MIXED / LIBC / NT  — the real-traffic accounting class
// NOTE  human-readable rationale (printed verbatim in the "RFO accounting" table)
```

The 13 vec3/vec4 entry points are already in the experimental area and take effect with `AB_OPS_EXTRA=1`.

Every one of those rows carries ELOP plus `AB_OPF_NO_BASELINE`, and registers its comparison partner in `AB_INTRA_PAIRS(X)`.

Upper-bound references bypass the variant TU and go through a separate table, `AB_MEM_OPS(X)`: 8 columns, no enum, no renaming.

**A completely worked example of adding a _stream entry point**, in exactly the form the current working tree has:

```c
/* 1) one line in ab_ops.h (the default 182-item gate must not change, so it goes in the experimental area) */
X(ADD_X4_STREAM, "fv_vec2_add_x4_stream", VEC2, EXACT, ref_add_x4, 0u, 16, 8, AB_RFO_NT,
  "NT store: 16 read + 8 write = 24 B/e (the 8 B RFO of the 32 B/e normal version is gone)")

/* 2) one thunk in ab_variant.cpp (its name must be abt_ plus the ID) */
static void abt_ADD_X4_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM   /* the Makefile greps this side's header automatically */
    fv_vec2_add_x4_stream((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#else               /* baseline header does not have it yet: it degrades to a normal store, so this line becomes
                     * "normal store vs NT store" comparison (32 vs 24 B/e) */
    fv_vec2_add_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#endif
}
```

Then:

    make -C bench/ab AB_OPS_EXTRA=1 build      # compile
    ./bench/ab/build/ab_bench --list-ops       # confirm it is in the table and the byte counts are right
    make -C bench/ab verify OP=fv_vec2_add_x4_stream N=1048576

**The experimental-area switch**: `AB_OPS_EXTRA=0` is the default value, which expands `AB_LIB_OPS_EXTRA` to nothing, so the default gate is still the historical 182 items — 15 ops × 8 tiers + 6 ops × 5 small N.

Once a new entry point is stable and you want it in the default matrix, move its line in ab_ops.h from the EXTRA table into AB_LIB_OPS.

**The new layout (vec3 / vec4)**: when the element width is not 8 B, one more change is needed beyond those two places — extend `AbArgs` + `Buf` with pointer fields such as az/bz/oz and an element-width rule.

That change is one-off and independent of how many ops exist. After it, every new op with the same layout still costs only a table row and a thunk.

## 9.5 External load on this machine is **permanent** — absolute loadavg cannot gate a measurement

**Evidence**, sampled at a single instant:

    $ cut -d' ' -f1-3 /proc/loadavg     ->  3.34 3.87 4.69
    $ pgrep -c cc1plus                  ->  0
    $ pgrep -c 'cc1 '                   ->  0

`loadavg` sits at 3–9 permanently. During this round's closing campaign it went **from 3.03 to 5.41**, while **no cc1plus / cc1 process existed anywhere in the whole interval**.

In other words the load is **not our compiler**; it is resident external processes on the machine.

**The conclusion, in one sentence**: **on this machine absolute loadavg cannot gate a measurement; the A/A floor has to do the controlling.**

Three operating rules follow:

1. **A quiet window is defined as "no compiler process of ours"**, not by loadavg:
   `pgrep -f 'cc1plus|cc1 |make -C bench/ab|g\+\+ .*flashvec'` is empty for 60 s straight.
   That is what campaign.sh uses. This round needed only 5 s to reach a clean window; the old `loadavg < 1.5` criterion would never have opened.
2. **loadavg is recorded, never rejected on**: note the start and end values, and only when the **rise is > 1.0** warn that "this stretch of measurement may have been bitten".
3. **Control a measurement with a measurement**: the effect of contamination on each (op,size) point is quantified directly by the **same-round interleaved A/A floor**, under the rule "floor > 15% → that point is INCONCLUSIVE".
   This round's campaign saw loadavg rise by 2.38, yet no global "contaminated" label was applied; instead **38 points** were judged INCONCLUSIVE because their own floor exceeded 15%.
   A global label would have thrown out the clean points along with the rest; a per-point floor does not.

## 10. The A/A floor table

`docs/research/raw/ab_aa_floor.tsv` holds the **self-check ratio** for every (phase, op, n) — that is, the resolution this machine can deliver when **both sides run the same code**:

    phase  op  n  ratio  pair_iqr_pct  trials  resolvable  floor_dev  in_band  rfo_flags

* A point with `resolvable=0` cannot even measure its own floor (paired IQR >= --cov); **do not judge with it**;
* `floor_dev = |ratio-1|`: any "the candidate is r times faster" conclusion can only be true if `(r-1) > floor_dev`;
* `in_band`: whether the ratio falls inside [0.85, 1.15] — the self-check band the harness brings with it;
* This table is valid only for a build where **both sides use the same source**. The default `make run` is that case; run-worktree is not, so it does not write this table.

## 10.5 INTRA mode: measuring a new entry point that has no baseline

The `fv_vec4_*_x2` / `fv_vec3_*_soa` entry points **do not exist in git HEAD at all**.

Comparing them against "nothing" measures "from zero to something", not "how much faster". So the registry tags them `AB_OPF_NO_BASELINE`, and:

* Normal A/B mode **refuses** to count them in; naming one with --op switches to INTRA automatically.
* `--intra` mode compares **inside the candidate**: the new entry point vs the reference entry point of the same byte count named in the registry's `AB_INTRA_PAIRS`.
  The pairings are vec4 ↔ vec2_x4, vec3 ↔ vec2_soa, and vec4_soa ↔ vec2_soa.
* **Normalized by byte**: within one tier sweep every op moves the same byte count, element count = n*8/element width. So ratio = t(ref)/t(new) answers "who is faster for the same amount of work".
* The report gives both ns/element and GB/s and **labels them explicitly**: `INTRA (comparison inside the candidate, not an improvement over HEAD)`.

Element width lives in the registry's **ELOP column** with values 8 / 12 / 16, and is not inferred from the output form. The reason: dot/len produce float/NONE output while the element width may still be 12/16 B.

## 10.6 Three verdict levels: PASS / FAIL / INCONCLUSIVE

**FAIL gets misread as "this optimization does not work", when in practice it usually means "this machine cannot measure it".** The two have to be kept apart:

| Level | Trigger condition |
|---|---|
| FAIL | R1 bit-level cross-check fails; R5 sees other sizes regress by more than 5%; the median ratio is clearly below 1 and beyond the floor (a measured slowdown) |
| INCONCLUSIVE | gain ≤ A/A floor; A/A floor > 15% (the point is unresolvable); same direction but short of 1.15; blacklisted tier (n=131072) |
| PASS | R1 passes + all four ratios ≥ 1.15 + smallest gain > A/A p90 floor + floor ≤ 15% + no R5 regression |

## 10.7 Per-size criterion: a floor, not a constant band

In `--selfcheck`, "the median ratio of each tier" originally had to fall inside a hard-coded `[0.90, 1.10]`.

Measured: 2 of 4 runs exited 3 over that, while **the checksum gate passed 182/182 every time**. What failed was the constant band, not the harness. Changed to:

1. The criterion is **that tier's measured A/A floor p90**, produced by `--aa-each-round` in the same batch of runs and therefore from the same distribution.
   That is, `|tier median ratio - 1| <= tier floor p90`. The constant band survives only as the degraded path for when there is no same-round A/A data.
2. **n=131072 is excluded from the per-size criterion** (`--tier-exclude=131072`, configurable). It sits on the L3 boundary and its own A/A floor ranges over 0.50–2.22; judging with it amounts to rolling dice.
   It is still measured, still recorded, and **still part of the gate**.
3. **Backstop**: if a tier's own floor exceeds 20% → that tier is marked INCONCLUSIVE, neither PASS nor FAIL. Because "cannot measure" and "measured a shift" are two different things.

## 10.8 Negative result: `_mask_inplace` has no measurable throughput advantage (5 densities × 4 sizes, all INCONCLUSIVE)

The campaign measured `fv_vec2_add_x4_mask_inplace`: densities 1% / 5% / 20% / 50% / 95%, sizes 8192 / 65536 / 1M / 4M.

**All 20 points are INCONCLUSIVE**, gains **-0.10 … +0.00**, same-round A/A floor 1.06–1.56.

**The reason is in the traffic accounting.** An in-place write is `a[i] += b[i]`; the target array **has already been read** before the write, so the cache lines involved are already in L1/L2.

A normal store **generates no RFO** in that situation. So what in-place saves over a separate output is not 8 B/element — it is 0.

This row's RFO class in the registry is `AB_RFO_HIT`: real traffic 25 B/e = nominal 25 B/e; the separate-output version is `AB_RFO_MISS` at 33 B/e real.

Sparsity only changes how many elements get written, not whether each written element needs an RFO, so a density sweep cannot turn up a difference either.

**The lesson for the harness**, which is the real point: an optimization assumption that looks self-evident — "sparse in-place writes are faster" — simply does not hold up in the traffic accounting.

**Build the harness first, then optimize**, otherwise this entry point gets written into the contract on the intuition that "it should be faster", and the next person cites it as fact.

The value of a harness is not that it proves an optimization works, but that it **picks out the invalid optimization assumptions before the code is written**.

## 11. The five acceptance criteria (verify subcommand)

`make -C bench/ab verify OP=<op> N=<size>` runs everything in one go and prints the verdict table:

| Rule | Criterion | Data source |
|---|---|---|
| R1 | bit-level cross-check (all tail lengths n=N-7…N+3, both builds) | `--gate-only` twice + the gates of four timing runs |
| R2 | target-size median ratio within one process >= 1.15x | ratio(AB, build#1) |
| R3 | **same-direction reproduction**: AB/BA (order swapped) × two independent builds (placement swapped), all four ratios >= 1.15x and in the same direction | 4 runs |
| R4 | the gain (r-1) **exceeds** that point's A/A floor | the **same-round interleaved A/A**, see below |
| R5 | no other size regresses by more than 5% | one multi-size sweep; expected L2/L3 regressions are declared explicitly with `EXPECT_REGRESS=65536,262144` |

R4's floor is measured like this: every A/B round is followed immediately by a pair of A/A, in the same process, on the same data, in the same residency state.

Across 2 builds × 2 orders × ROUNDS rounds that is ≥16 samples; p50/p90/max are taken and p90 is the criterion.

On FAIL it prints **which rule failed and by how much**, e.g. `ratio=0.999 < 1.15 (no measurable gain at this size within one process)`, instead of settling for a red cross.

Knobs: `MIN_GAIN` (default 1.15), `MAX_REGRESS` (default 0.95), `CHECK_SIZES`, `EXPECT_REGRESS`.

**What makes the two builds "independent"**: the same source, with only code placement changed (`-falign-functions=64 -falign-loops=32 -falign-jumps=32`) and zero change in semantics.

This exists specifically to turn off one-off factors like code placement, which is the main source of the A/A floor.

**Why A/A has to be same-round interleaved (the fix behind R4)**: the old version compared a floor from one set of runs against a gain from another set, and the two are not the same distribution.

Now each round inside each process measures A/B first and then a pair of A/A (`--aa-each-round`), so floor and gain share one stretch of frequency drift, one data set, and one residency state.

Measured effect: on one op with no real difference, the old A/A produced 1.208 and "passed" R4; the new version pins that point's floor p90 at 0.0805 with a gain of -0.2393 and judges it INCONCLUSIVE outright, because the gain sits below the floor.

**n=131072 is a blacklisted tier**: its measured A/A ratio ranges over 0.50–2.22 (L3 boundary + 4K aliasing). It may be measured and recorded, but **must not be used as a criterion**.

Both verify and campaign treat it as BLACKLIST.

## 12. RFO accounting: nominal traffic vs real traffic

GB/s used to be computed from **nominal bytes** (rd+wr). But a normal store on a write miss first **reads** the whole cache line in.

That is RFO, read-for-ownership, and the real traffic is "one write plus one read":

    real_wr = wr * 2   normal store (write-allocate miss)
    real_wr = wr        NT store (_mm256_stream_ps: bypasses cache, generates no RFO)
    real_wr = wr       in-place update (target already read before the write, line already in cache)
    real_wr = wr       libc decides (memcpy/memset: labelled libc?, not filled in by force)

Real GB/s = (rd + real_wr) / ns_per_elem, reported alongside nominal GB/s. In the matrix, `gbs` and `real_gbs` are separate metrics.

**The physical meaning**: the real benefit of an NT store is not "each byte is faster" but **25% fewer bytes moved**. The same add operator is nominal 24 B/e either way:

| | nominal B/e | real B/e | real GB/s | nominal/real |
|---|---|---|---|---|
| normal store | 24 | 32 | the same value | 0.75 |
| NT store | 24 | 24 | the same value | 1.00 |

So the **real** GB/s of the two should come out roughly equal, both at DRAM bandwidth, while the **nominal** GB/s differ by 1.33x.

If nominal differs by 1.33 and real also differs by 1.33, the gain does not come from traffic. That is one ruler for spotting a false positive.

Conversely: NT store only delivers that 25% when **the working set is far larger than L3 and the write stream is clean and one-shot**. When the working set fits in L2/L3 it is slower than the normal version.

**Measured 0.13–0.51x, i.e. 2–8 times slower.** At n=8192 all six entry points are judged FAIL (slower):

| Entry point | Ratio |
|---|---|
| add | 0.217 |
| sub | 0.301 |
| scale | 0.130 |
| add_soa | 0.229 |
| scale_soa | 0.138 |
| normalize | 0.506 |

These numbers come from the campaign verdict table `docs/research/raw/ab_campaign_verdict.txt`.

The old 0.40x in the header contract was a single-point measurement of a single entry point and understated how large the penalty is.

Uncertain points are **labelled, not filled in by force**:

| Label | Meaning |
|---|---|
| `mixed?` | pipeline intermediates may be absorbed by L3; 2x is an upper bound |
| `libc?` | memcpy/memset pick NT themselves |
| `miss*` | no NT entry point on this side; degrades to a normal store |

`--list-ops` prints each op's rd/wr/real byte counts and its class.

## 13. Where these criteria still leak (said up front)

1. **The A/A floor measures only the noise on the A side.** The floor comes from an A-vs-A ratio, while the jitter of an A/B ratio also contains the B side — another form, another object-file layout.
   When the two sides' noise differs, the floor underestimates the true A/B jitter, so it is a **lower bound**, not an estimate.
2. **The floor itself can double between rounds.** Three rounds of the 65536 tier measured A/A floors of 0.130 / 0.169 / 0.218, so the same tier can flip between "judgeable" and "INCONCLUSIVE".
   So "this tier is judgeable" is itself an unstable conclusion; only a **per-point PASS with a small floor at that point** counts as a hard result.
3. **R5 looks at only 5 sizes** (by default) and can miss a trap at a size in between.
4. **For libc and mixed ops the real traffic is an interval, not a number**; this table gives only endpoints and marks them with ?. No "exact" conclusion built on real_gbs holds.
5. **The A/A floor depends on build form**: when the baseline is library-defined and the candidate header-inlined, the "same code" on the two sides is not the same object-file layout either, so form differences leak into the floor.
   verify's two builds change code placement only, never form, so this point is **not** covered.
6. **R1 in verify covers only this op**, not the whole library's gate — that belongs to the top-level make test-all. It answers "this change did not break this op's bit-level semantics".
