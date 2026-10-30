# flashvec P0 / P1 Implementation Log and Acceptance

> English | [简体中文](zh-CN/implementation-log.md)

> Corresponds to the rework roadmap in §5 of `docs/optimization-report.md`.
> This log records what was actually done, what the evidence is, and what is still missing. It does not
> repeat the analysis in the report.

## 1. Scope of delivery

| Phase | Item | Status |
|---|---|---|
| P0 | Add `-O3 -fno-math-errno -ffp-contract=off` to the Makefile | done |
| P0 | `min / max / clamp / madd` for `fv_vec2`/`fv_vec3` | done |
| P0 | `fv_aabb2` + `fv_aabb2_merge / overlaps` | done |
| P0 | Fix the `FV_EPS` dimensional bug: add `FV_EPS_SQR`, switch `project` to a squared threshold | done |
| P1 | AoS 4-way interleaving: `_x4` for `add / sub / scale / lerp / madd` | done |
| P1 | Sparse variants: `add_x4_mask`, `add_x4_idx`, `scale_x4_idx` | done |
| P1 | SoA: `_soa` for `add / scale / dot / len / normalize` | done |
| P1 | AVX2 kernels + runtime dispatch; non-AVX2 machines take the scalar fallback | done |
| Refactor | **Batch entry points moved into the header as `static inline`** (default form); same-name `extern` symbols kept for `FV_USE_EXTERN` and old consumers | done |
| P1+ | NT store (`_stream`): 6 entry points; `_mask_inplace`: 1 | done |
| P1+ | `fv_vec4` type (18 scalar functions) + 5 `_x2` AoS batch + 3 `_soa` reduction batch | done |
| P1+ | 5 SoA batch functions for `fv_vec3` (no padding on the three-component form, SoA only) | done |
| P6 | `fv_quat` + `fv_mat3` scalar API (rotation / 3x3) | done |
| P6 | SoA batch: quat multiply / normalize / rotate-vector, mat3 transform | done |
| P6 | General inverses: `fv_quat_inverse` (+checked), `fv_mat3_determinant` / `fv_mat3_inverse` (+checked) | done |
| P6 | Tests: C/C++ x inline/extern/pure-C + FMA consumer, batch bit-identical to per-value | done |
| P6 | Preliminary A/B measurement (`make bench-physics`) | done (preliminary, not five-criterion acceptance) |
| — | Docs: header contract comments, `README.md`, this log, `docs/acceptance-protocol.md` | done |

The `docs/acceptance-protocol.md` in the last row is the document that holds the acceptance criteria.

The changes touch these files:

| File | Notes |
|---|---|
| `src/flashvec_batch.cpp` | new; batch entry points and AVX2 kernels |
| `src/flashvec.cpp` | P0 functions + the `project` fix |
| `include/flashvec.h` | new declarations + contract comments |
| `tests/test_p0.cpp` | |
| `tests/test_batch.cpp` | |
| `tests/test_batch_c.c` | new |
| `Makefile` | compile flags + new targets; existing `all/test/clean` behavior unchanged |
| `verify/` | independent third-party acceptance |

## 2. Acceptance evidence (four independent lines that corroborate each other)

### 2.1 Correctness

| Source | Result |
|---|---|
| Implementer's tests | 0 failures |
| Third-party acceptance (`verify/parity.cpp`) | 12/12 PASS, exit 0 in all three configurations, no ASan/UBSan reports |
| Dispatcher's independent recount, reusing none of the above code | 423/423 pass |
| Negative controls | 15 / 26 FAILs respectively |

What each of the four lines actually tests:

- The implementer's tests run `test_flashvec` / `test_p0` / `test_batch` / `test_batch_c`: 30 / 57 / 1098
  assertions, plus the C batch interface.
- Third-party acceptance: 12 element-wise batch functions, in three configurations: `-O2`,
  `-O3 -march=native`, ASan+UBSan.
- Third-party coverage: every tail length n = 0..17; four mask modes; idx in ascending / descending /
  duplicated / shuffled order; out-of-bounds write canaries.
- Dispatcher's independent recount: 423 checks covering the x4 set of five, mask semantics, reversed idx,
  SoA add/len/normalize, and the `project` threshold.
- Negative controls: the lerp kernel deliberately rewritten to use FMA, and the mask blend deliberately
  inverted; these produced 15 / 26 FAILs respectively.
- Why the negative controls matter: they show the tests really do catch errors, and that the AVX2 kernels
  really run.

### 2.2 Semantics and safety

- **Bit-identical**: the `_x4`/`_soa`/`_mask`/`_idx` variants of `add/sub/scale/lerp/madd`, plus
  `len_soa` and `normalize_soa`, produce bit patterns identical to the scalar API.
- How: these paths use only `vaddps/vsubps/vmulps/vdivps/vsqrtps`, **no FMA and no rsqrt/rcp
  approximation**, and they keep the same reduction order.
- `objdump` survey of the library: `vfmadd` = 0, `rcpps/rsqrtps` = 0, `haddps` = 0.

- **The one exception, `dot_soa`**: 4 accumulators plus a single reduction at the end, so the result
  depends on the reduction order.
- Measured: for n ≤ 1000 the relative error is ≤ 4e-8.
- For n = 8192 a single dot product has absolute error ≤ 4.8e-7 and relative error ≤ 6.2e-8; with
  near-cancellation it can reach ~3e-6.
- The header has been rewritten accordingly: a corpus-dependent value is no longer presented as an upper
  bound.

- **Mask semantics**: where `mask[i] == 0` the original `a[i]` value is explicitly written back; this is
  not "leave it unwritten".
- Acceptance includes a canary check: positions not referenced by `_idx` are rewritten 0 times.

- **Usable from pure C**: `tests/test_batch_c.c` and `tests/test_header_c.c` compile and run under
  `cc -std=c11 -pedantic-errors`.
- **Careful: the accounting changed with the refactor.** The batch entry-point implementation now lives in
  the header itself.
- So a pure C translation unit on x86 with GCC/Clang **will** expand `<immintrin.h>`; measured
  `cc -std=c11 -E | grep -c immintrin` = 52.
- The statement "there are no intrinsics in the header" still holds in only two forms: `FV_USE_EXTERN`
  (declarations only, with an include set byte-for-byte identical to the pre-refactor one) and
  `FV_I_NO_INTRINSICS` (the pure C fallback).
- Sections [1] and [2] of `verify/checks.sh` check exactly those two forms.

- **Existing behavior is untouched**: the original 30 assertions are not edited at all.
- Across four builds (old flags with new sources, new flags, `-march=native`), the bit-pattern hashes of
  the existing scalar API outputs are identical.
- The only change is the `project` threshold fix: |onto| = 1e-4 went from being wrongly judged zero to
  being projected correctly, while 1e-7 is still judged zero.

### 2.3 Performance (this machine: AMD Ryzen 7 5700U, AVX2/FMA, no AVX-512)

N = 8192, 32 KB per array, L2-resident. The table below gives speedups relative to a **scalar loop that
calls the library**.

Data from the third-party `verify/perf.cpp`, median of 5 samples, checksums bit-identical.

| Operation | Scalar API | `_x4` | `_soa` |
|---|---|---|---|
| add | 2.601 ns/elem | 1.240 (2.10x) | 1.369 (1.90x) |
| scale | 2.776 | 1.079 (2.57x) | 1.027 (2.70x) |
| lerp | 3.002 | 1.196 (2.51x) | — |
| dot | 3.087 | — | 0.152 (**20.33x**) |
| len | 3.046 | — | 1.149 (2.65x) |
| normalize | 4.382 | — | 1.306 (3.35x) |
| collision microkernel | 8.481 | — | 1.670 (5.08x) |
| add_mask (50% hit) | 4.648 | 1.223 (3.80x) | — |
| add_idx (50% hit) | 2.248 | 1.670 (1.35x) | — |

At N = 1<<20 the working set is 24 MB, far beyond the 4 MiB of L3 visible to a single thread here.

The chip total is 8 MiB, two CCX of 4 MiB each, so the knee falls at 1.75×10⁵ elements rather than
3×10⁵.

In that range `add` goes 13.977 → `_x4` 9.722 (1.44x) / `_soa` 9.455 (1.48x).

`scale` instead lands at 0.92–0.93x, back at the memory bandwidth ceiling, as predicted in §2 of the
report.

## 3. Three qualifications that must be read together with the numbers

1. **The gain comes mainly from the cross-TU call boundary disappearing, not from intrinsics.**
   - With the same math inlined from the header, `-O3 -march=native` auto-vectorization nearly matches
     hand-written AVX2: on add, `scalar_flat` 2.15x vs `_x4` 2.10x.
   - But **a caller that links a prebuilt static library cannot inline**, and for them the 2.10x of `_x4`
     is real.
   - And `normalize` (3.35x), `dot` (20.33x), `collision` (5.08x) are out of reach even for an inlined
     scalar loop; the scalar path gets only 1.13x / 3.02x / 1.36x respectively.

2. **Size decides everything**: at N = 1<<20 the variants converge to 0.9–1.5x. Quoting a speedup means
   quoting the working-set size as well.

3. **Absolute values are not comparable across fixtures**: this laptop scales frequency (1.64–3.97 GHz)
   and is a shared machine, so the same code can differ by 20–25% between runs.
   - The three fixtures (implementer / third party / dispatcher) give different absolute values; `add_x4`,
     for example, ranges from 0.35 to 1.24 ns/elem.
   - Only **ratios within one run** are meaningful.

## 4. What was deliberately not done (and why)

- **The `_idx` variants have no AVX2 kernel**: they are scalar loops, measured at only 1.08–1.35x. The
  reason is that `vgatherdps` costs about 10 cycles per instruction on Zen 2 (§3.4 of the report).

- **So do not size a broad-phase expectation from the 1.7–3.8x of `_mask`; use 1.1–1.35x.**

- **The tail goes scalar** (`_x4`: 0–3 elements, `_soa`: 0–7), with no masked loads: tail length has no
  performance meaning.

- ~~Batch entry points for the three-component vectors not done~~ → **done this round, via SoA**:
  `fv_vec3_add_soa / scale_soa / dot_soa / len_soa / normalize_soa`.

- `fv_vec3` is 12 B, so packing it into a ymm would idle 1/4 of a lane; hence it is neither packed AoS
  nor is `fv_vec3` given padding (padding would break the ABI).

- For the three-component form the layout is instead "8 elements of the same component in one ymm".

- ~~Four-dimensional vectors not done~~ → **done this round**: a 16 B `fv_vec4` with 18 scalar functions + 5 `_x2` AoS
  batch + 3 `_soa` reduction batch.

- One ymm of `_x2` holds exactly 2 elements; the number means the same as in `_x4`: how many elements fit
  in one 256-bit register.

- **The `_fast` variants (rsqrt + Newton iteration) are not implemented**: they change the last bit and
  must be offered as new functions; and they are measured at only about 10% faster.

- **Runtime dispatch has been verified only on machines with AVX2**: the equivalence test that forces the
  scalar fallback passes, and turning the 11 dispatch sites into `if (0)` leaves 30+57+1098 all green.

- But it has **never run on real non-AVX2 hardware**. That is currently the only path without
  hardware-level evidence.

## 5. Reproduction

*(Comments inside the block are the fixture's own and are reproduced verbatim: 条断言 = assertions,
三档 = three configurations, 纯 C 兼容性 = pure C compatibility, 校验和一致性断言 = checksum-consistency
assertions, 另一个夹具，看趋势 = a different fixture, look at the trend, 指令普查 = instruction survey.)*

~~~bash
make clean && make test-all      # 30 + 57 + 1098 条断言 + 两个 C 兼容性测试
cd verify && make parity         # 第三方位级一致性（O2 / O3-native / ASan 三档）
cd verify && make perf           # 第三方性能对比（含校验和一致性断言）
cd verify && make check          # 结构检查：头文件纯净度、编译选项、派发、指令普查
cd ../bench && make run          # 三档实现 vs 标量基准（另一个夹具，看趋势）
~~~

## 6. Open decisions (mapping to §6 of the report)

| Item | Status |
|---|---|
| `FV_EPS` dimension | Implemented as the split into `FV_EPS` + `FV_EPS_SQR` (the recommended option of decision 1, §6 of the report) |
| Whether to introduce a batch API | Introduced `_x4` + `_soa` (the recommended option of decision 2, §6 of the report) |
| `-march=native` | Not introduced (the recommended option of decision 3); the default build keeps the x86-64 baseline plus function-level `target` attribute dispatch |
| `-fno-math-errno` in the defaults | Added (the recommended option of decision 4) |
| Whether `fv_vec4` / three-component batch are needed | **Decided and implemented**, see below |
| Whether `_fast` variants are needed | Open, pending real performance data (still not done) |
| `_stream` (NT store) | 6 entry points implemented, see below |
| `_mask_inplace` | 1 entry point implemented, see below |

- `fv_vec4` and the three-component batch: the three-component form goes SoA (5 entry points); the
  four-component form uses `_x2` (5) + `_soa` (3) + 18 scalar functions.

- `_stream`: 6 entry points implemented, chosen explicitly by the caller.

- `_stream` readings: measured inside L2 at **0.13–0.51x**, i.e. 2–8× slower; it pays off only when the
  working set far exceeds L3.

- `_stream` ceiling: PASS points beyond L3 at 1.29–1.58x, with ceilings of 1.333 / 1.500 set by the
  entry-point shape (`docs/acceptance-protocol.md` §9).

- `_mask_inplace`: 1 entry point implemented; a sparse in-place add that leaves **values** outside the
  mask unchanged, with the kernel chosen from the mask density.

- `_mask_inplace` verdict: **INCONCLUSIVE in every configuration**; no measured throughput advantage over
  `_mask`.

- Its value is in the form: an in-place write produces no RFO but saves no traffic.

## 7. Post-refactor architecture facts and acceptance criteria (added this round)

1. **There is exactly one copy of the implementation, and it is in the header.**
   - Batch entry points are now `static inline` in `include/flashvec.h`.
   - That header holds the AVX2 kernels, the `target("avx2")` attribute, and the
     `__builtin_cpu_supports("avx2") != 0` dispatch.
   - Measured counts: 32 dispatch sites and 7 `target(` in the header; **0 of each** in `src/*.cpp`.
   - `src/flashvec_batch.cpp` is left with 33 `extern` forwarding symbols, for binary compatibility.
   - Why this change: the cross-TU call boundary is a first-order bottleneck (measured in item 1 of §3),
     and a caller linking a prebuilt static library cannot remove that boundary.

2. **The three forms must give the same result**: default (header inline), `FV_USE_EXTERN` (implemented
   in the library), and `FV_I_NO_INTRINSICS` / pure C fallback.
   - `make test-inline` / `test-vec34-inline` compare the main output streams directly.

3. **The header now conditionally includes `<immintrin.h>`**, guarded by
   `FV_I_INTRINSICS && !FV_USE_EXTERN`.
   - Entry points that could be contracted into FMA carry `#pragma GCC optimize("fp-contract=off")`.
   - `target("avx2")` does not stop a `-mfma` coming from the caller's command line: under
     `-march=native -mfma -ffp-contract=fast`, `vfmadd` still appears.
   - Section [8] of `verify/checks.sh` runs an artifact-level instruction survey under the strictest
     configuration: 7 forbidden mnemonics = 0, and `vmovntps` > 0 as the positive control for NT stores.

4. **The acceptance criteria are now a document of their own**: see `docs/acceptance-protocol.md`.
   - An optimization is accepted only if it satisfies five conditions **simultaneously**:
     - Bit-level comparison.
     - Same-process median ratio at the target size ≥ 1.15x.
     - Same-direction reproduction.
     - Gain above the A/A floor at that point.
     - No other size regresses by more than 5%.
   - "A single-point ratio is not trustworthy" has hard evidence: of 304 large-N measurement points,
     only about 30–35% are resolvable.
   - Across four runs the resolvable counts were 90 / 92 / 93 / 107, always out of 304. The count drifts
     between runs, so cite it together with the round.
   - Writing it as a range does not weaken that evidence; it is where the evidence is strongest: on the
     same machine with the same code, how many points are resolvable itself drifts between 90→107.
   - A single-point ratio for one `(op, n)` is then, of course, even less reproducible.

5. **Reduction-type entry points are not split invariant**; this is a contract fact added this round and
   must be read together with the rest.
   - `fv_vec2_dot_soa` / `fv_vec3_dot_soa` / `fv_vec4_dot_soa`: split the same array into several
     segments, call each one, add the results, and the result is **different** from one call over the
     whole array.
   - Measured: relative difference ≈3.9e-6 at n=1000 and ≈2.4e-5 at n=8192; all 12/12 corpora disagree.
   - In the same experiment `fv_vec2_add_x4` is **bit-identical** in all 12/12 corpora; that is the
     control group showing the experiment works.
   - This is a **design consequence** of the reduction tree changing with the split, not a defect. But it
     is a different statement from "error against the scalar path ≤ 4e-8"; the two must not be mixed.
   - Callers who need split-independent results should use the scalar `fv_vec2_dot` with sequential
     accumulation, or fix the split.
   - Element-wise batch entry points (`add/sub/scale/lerp/madd/len/normalize` and their
     `_x4`/`_soa`/`_mask`/`_idx` variants) are split invariant.

6. **The cache accounting has been corrected.**
   - L3 here is **4 MiB/CCX × 2 CCX = 8 MiB chip total, with 4 MiB visible to a single thread**, read
     from `/sys/.../index3`: `size=4096K shared_cpu_list=0-7`.
   - The old document said "L3 8 MiB shared" and computed the knee as `3×10⁵` elements from that; that
     was wrong.
   - The real knee is about **1.75×10⁵** elements (@24 B/element), and the measured cliff falls between
     2^17 (38–57 GB/s) and 2^18 (5.6–8.8 GB/s).

## 8. P6: quaternion / mat3 physics primitives (added this round)

### 8.1 What landed

| Item | Status |
|---|---|
| `fv_quat` (16 B, x/y/z/w, identity (0,0,0,1)) + `fv_mat3` (36 B, row-major) | done |
| Scalar quat API: `make / add / sub / neg / scale / dot / len_sqr / len / normalize / conjugate / mul / rotate_vec3` | done |
| Scalar mat3 API: `from_rows / transform_vec3 / mul / transpose / from_quat` | done |
| General inverse: `fv_quat_inverse` + `_checked`; `fv_mat3_determinant` / `fv_mat3_inverse` + `_checked` | done |
| SoA batch: `fv_quat_mul_soa`, `fv_quat_normalize_soa`, `fv_quat_rotate_vec3_soa`, `fv_mat3_transform_vec3_soa` | done |
| Tests: C++ scalar + batch parity, C11 inline / `FV_USE_EXTERN` / pure-C, FMA consumer | done |
| Preliminary A/B measurement (`make bench-physics`) | done (preliminary) |

New files: `src/flashvec_quat.cpp`, `src/flashvec_mat3.cpp`, `src/flashvec_physics_batch.cpp`,
`tests/test_quat.cpp`, `tests/test_mat3.cpp`, `tests/test_physics_c.c`, `tests/bench_physics.cpp`.
Edits: `include/flashvec.h` (types, scalar declarations, batch declarations, implementation, public
wrappers), `Makefile` (objects, test/bench targets), `README.md` (this group), this log.

### 8.2 Contract decisions (why these, not other choices)

- **Quaternion normalize maps length <= FV_EPS to the identity (0,0,0,1)**, not to zero: a zero
  quaternion is not a rotation. This is the one deliberate difference from `fv_vec*_normalize`; the
  header states it and `tests/test_quat.cpp` pins it.
- **Unit precondition**: `fv_quat_rotate_vec3` and `fv_mat3_from_quat` assume |q| = 1 and do not
  normalize -- the same rule as a rotation matrix assuming an orthonormal basis. Normalization is a
  separate, explicit call.
- **Composition order**: `fv_quat_mul(a, b)` = a*b applies b first and a second; `fv_mat3_mul(a, b)`
  uses the same order, and `fv_quat_mul(q, conjugate(q))` is the identity for unit q.
- **mat3 layout**: 36 B, row-major, compact 9 floats, no padding. Padding to 48/64 B is a caller-side
  AoS choice; the library's batch form is SoA and does not need a padded struct. The row-major
  convention is the same as `fv_mat4` (and the opposite of GLM).
- **Inverse is the general one** (`conjugate(q)/|q|^2` for a quat, adjugate/det for a 3x3), not just
  the unit/orthonormal shortcut. Degenerate input returns the **zero** value, not the identity: a hidden
  identity would make a degenerate input look valid. The quat threshold is `|q|^2 <= FV_EPS_SQR` (the
  dimensionally correct form of `length <= FV_EPS`); the mat3 threshold is `|det| <= FV_EPS *
  max|entry|^3`, scale-relative, so a uniformly scaled invertible matrix stays invertible. Both
  `_checked` variants return 1/0 and write the zero value on failure; `fv_mat3_inverse_checked` allows
  `out == m`.
- **No slerp / axis-angle / Euler this round**: they need sin/cos and an input-domain contract, which
  is a separate design decision (see section 4).

### 8.3 Evidence

| Source | Result |
|---|---|
| `make test-quat` | 26,321 checks, 0 failures: scalar semantics (including general inverse) + batch bit-parity, n = 0..40 (every tail) + 63 / 1024 / 4097 |
| `make test-mat3` | 14,182 checks, 0 failures: same coverage for the 3x3 API, including determinant/inverse/singular paths |
| `make test-physics-c` | 3,916 checks, 0 failures; the same source also passes under `-DFV_USE_EXTERN` and `-DFV_I_NO_INTRINSICS` |
| `make test-physics-fma` | both C++ tests pass under `-march=native -mfma -ffp-contract=fast`; `vfmadd` = 0 in both consumer objects |
| Vectorization census | quat mul 40, rotate 40, mat3 21 AVX2 instructions including loads/stores, 0 `vpermps`/`vshufps`; normalize 29 (explicit kernel) |
| `make test-all` | green; `test-physics` is appended to the target list |

The AVX2 census is `objdump` of `build/flashvec_physics_batch.o` at the default tier; the FMA census
compiles the same source with `-march=native -mfma -ffp-contract=fast` and counts `vfmadd` in each
kernel (0 in all four). The batch-vs-scalar comparison is `memcmp`/bit equality, not a tolerance.

### 8.4 Preliminary performance

`make bench-physics`: same-process interleaved A/B, 11 paired rounds, `taskset -c 0`, medians of two
runs. The same-round A/A floor was 0.99-1.00 at n=8192 and 0.93-1.07 at n=1M.

| Entry | n = 8192 (L2) | n = 1M (32-60 MiB working set) |
|---|---|---|
| `fv_quat_mul_soa` | 18.0-19.5x | 2.39-2.51x |
| `fv_quat_normalize_soa` | 41.9-42.3x | 4.56-5.13x |
| `fv_quat_rotate_vec3_soa` | 25.8x | 3.39-3.57x |
| `fv_mat3_transform_vec3_soa` | 24.7-25.1x | 3.33-3.80x |

Baseline: per-element library calls (cross-TU, not inlinable). At 8192 the batch runs at
0.4-0.55 ns/element against 8-20 ns for the scalar calls; the removed call boundary is most of the
gain, consistent with the rest of the library. These numbers are **preliminary**: the five-criterion
acceptance of `docs/acceptance-protocol.md` has not been run for P6. Quote them with the harness and
the run.

### 8.5 What remains

- Quaternion interpolation (slerp / shortest-arc lerp), axis-angle and Euler construction (sin/cos +
  domain contract).
- Full `fv_mat4` operations (multiply / inverse); quaternion-from-matrix; batch `fv_mat3_mul` (the
  per-body scalar 3x3 inverse is the current use case, so a batch inverse is not planned yet).
- Absolute speedup for the 3D/4D entries (still unmeasured).
- A hardened P6 A/B fixture under the five acceptance criteria; `bench-physics` is a preliminary
  harness, not the acceptance one.
