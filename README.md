# flashvec

> English | [简体中文](zh-CN/README.md)

A small vector math library for real-time 2D / 3D / 4D computation.

- Every value is passed by value. The library never allocates and holds no global state.
- The header is a pure C interface. It can be included from both C and C++.
- Two APIs: scalar for correctness baselines, batch for hot paths.
- The default path is bit-for-bit reproducible: identical input always yields identical output bits.

## Quick start

```c
#include "flashvec.h"

fv_vec2 a = fv_vec2_make(3, 4);
fv_vec2 n = fv_vec2_normalize(a);   /* (0.6, 0.8) */
fv_real d = fv_vec2_len(a);         /* 5 */

/* Batch: one call processes a whole array */
fv_vec2 pos[N], vel[N];
fv_vec2_add_x4(pos, vel, pos, N);
```

To use it, add `include/` to your header search path. By default the header contains the entire
implementation, so nothing needs to be linked. To use the prebuilt library instead, define
`FV_USE_EXTERN` and link `libflashvec.a`.

## Types

| Type | Size | Notes |
|---|---|---|
| `fv_vec2` | 8 B | |
| `fv_vec3` | 12 B | No padding, so existing layouts stay valid |
| `fv_vec4` | 16 B | Size locked by `static_assert` |
| `fv_aabb2` | 16 B | Axis-aligned bounding box |

## Scalar API

Each vector type has 18 functions (`fv_vec3` has one extra, `cross`):

`make` `add` `sub` `neg` `scale` `lerp` `min` `max` `clamp` `madd`
`dot` `len_sqr` `len` `dist_sqr` `dist` `normalize` `project` `near`

The bounding box adds `fv_aabb2_merge` and `fv_aabb2_overlaps`.

`min`, `max` and `clamp` follow `(a < b) ? a : b` for NaN, not IEEE `fminf`. This is stated in the
header.

## Batch API

33 entry points, grouped by memory layout. The number in a name is how many elements fit in one
256-bit register: four `fv_vec2` (`_x4`) or two `fv_vec4` (`_x2`). Both fill the register
completely, so no lane is wasted.

### AoS (array of structs)

Pass your existing arrays as they are.

| Type | Functions |
|---|---|
| `fv_vec2` | `add_x4` `sub_x4` `scale_x4` `lerp_x4` `madd_x4` |
| `fv_vec4` | `add_x2` `sub_x2` `scale_x2` `lerp_x2` `madd_x2` |

### SoA (component arrays)

Components live in separate arrays. Reductions are only efficient in this layout.

| Type | Functions |
|---|---|
| `fv_vec2` | `add_soa` `scale_soa` `dot_soa` `len_soa` `normalize_soa` |
| `fv_vec3` | `add_soa` `scale_soa` `dot_soa` `len_soa` `normalize_soa` |
| `fv_vec4` | `dot_soa` `len_soa` `normalize_soa` |

3D only has the SoA form. `fv_vec3` is 12 B, so packing two per register leaves 8 B unused; filling
it requires shuffling apart and back together, which gives back the width you saved. Padding the
struct would break the ABI and make it 33% larger. So 3D packs eight elements of the *same
component* into one register, at 100% lane utilization.

### Sparse updates

| Function | Semantics |
|---|---|
| `add_x4_mask` | Add where the mask is non-zero; where it is zero, write back the original `a[i]` |
| `add_x4_mask_inplace` | In-place add. Values of masked-out elements are unchanged; they are not guaranteed to be left physically untouched |
| `add_x4_idx` | Process only the indices in `idx`; duplicates and out-of-order indices are allowed |
| `scale_x4_idx` | Same, for scaling |

The `_idx` variants use a scalar loop rather than `vgatherdps`. On Zen 2, gather is microcoded at
roughly 10 cycles per instruction, which does not pay off.

### Streaming stores

`add_x4_stream` `sub_x4_stream` `scale_x4_stream` `add_soa_stream` `scale_soa_stream`
`normalize_soa_stream`

These write results with non-temporal stores, bypassing the cache and skipping the read-for-
ownership traffic. They pay off only when the result is not read back soon and the working set is
far larger than L3. The library will never switch to them on its own; the choice belongs to the
caller.

## Guarantees

### Bit-for-bit identical

The following entry points are guaranteed to match the corresponding scalar API bit for bit, so
`memcmp` is a valid test:

- the `_x4`, `_x2`, `_soa`, `_mask` and `_idx` variants of `add`, `sub`, `scale`, `lerp` and `madd`
- `len_soa`
- `normalize_soa`
- every `_stream` function against its non-stream counterpart

The implementation uses only `vaddps`, `vsubps`, `vmulps`, `vdivps` and `vsqrtps`. No FMA, no
reciprocal approximations, and reduction order matches the scalar path.

### The one exception: `dot_soa`

`dot_soa` uses four accumulators and reduces once at the end, so its result depends on reduction
order and is not bit-for-bit identical. Below eight elements it degenerates to sequential
accumulation and does match. Use the scalar `dot` when you need exact reproducibility.

The header states two error figures that must not be confused:

| Statement | Meaning |
|---|---|
| relative error ≤ 4e-8 for n ≤ 1000 | difference from the **scalar path** |
| 2.4e-5 between whole and split | difference from a **different split** |

The first is not an upper bound for the second.

### Split invariance

Element-wise entry points (every batch variant of `add`, `sub`, `scale`, `lerp`, `madd`, `len` and
`normalize`) are split invariant: calling them on any number of segments of the same array produces
output bit-identical to a single call.

Reductions are not. `dot_soa`'s reduction tree changes with the split; measured on n = 8192, "one
call" versus "two segments summed" differ by 2.4e-5 in relative terms, and all 12 corpora disagreed.
In the same experiment `add_x4` agreed bit for bit on all 12.

This is inherent to reduction, not a defect. If you need results independent of the split,
accumulate sequentially with the scalar `dot`, or fix the split yourself.

### Lengths and pointers

`n` is arbitrary, including 0 and `NULL`. The main loop handles full blocks and a scalar tail
handles the rest.

`_x4` and `_soa` require the `__restrict` precondition: inputs and outputs must not overlap.

## Performance

Test machine: AMD Ryzen 7 5700U, 8 cores / 16 threads, AVX2 + FMA, no AVX-512.

### Cache boundaries

L3 is **4 MiB per CCX**, 8 MiB across both CCXs. A single thread only reaches its own 4 MiB. This is
visible in `/sys/devices/system/cpu/cpu0/cache/index3`: `size=4096K`, `shared_cpu_list=0-7`.

At 24 B per element for `add`, the knee is at 4 MiB ÷ 24 B ≈ **1.75e5 live elements**. Measurements
agree: 2^17 (3 MB working set) runs at 38–57 GB/s, while 2^18 (6 MB) drops to 5.6–8.8 GB/s.

### Speedup over scalar calls

The baseline is **calling the library's scalar function per element**, which is a cross-TU call that
cannot be inlined. Writing the loop yourself in the same translation unit lets GCC auto-vectorize
it, and the gap shrinks a lot (`dot` comes down to roughly 4x).

N = 8192, 32 KB per array, L2 resident:

| Operation | `_x4` | `_soa` |
|---|---|---|
| add | ~2.1x | ~1.9x |
| scale | ~2.6x | ~2.7x |
| lerp | ~2.5x | — |
| dot | — | ~20x |
| len | — | ~2.7x |
| normalize | — | ~3.4x |
| collision microkernel | — | ~5.1x |
| add_mask (5% / 50% / 95%) | ~1.7x / ~3.8x / ~3.3x | — |
| add_idx (5% / 50%) | ~1.1x / ~1.4x | — |

### Measurements for `_stream`

Data comes from `bench/ab`: same-process interleaved A/B with a same-round A/A floor, AB/BA order
swaps and two independent builds. Ratios are only comparable within one run.

| Scenario | Result |
|---|---|
| 24 MiB working set (n=1M) | add 1.350 / sub 1.293 / scale 1.313 |
| 192 MiB working set (n=8M) | add 1.422 / sub 1.519 / add_soa 1.582 |
| 96 MiB (n=4M) | inconclusive for all; that tier's A/A floor is 1.34–1.61 |
| 3 MiB (L3 boundary) | inconclusive for all; that tier's A/A ratios span 0.50–2.22 |
| Inside L2 (n=8192) | add 0.217 / sub 0.301 / scale 0.130 / add_soa 0.229 / scale_soa 0.138 / normalize_soa 0.506 |

The last row is the counter-example. The data is already cached, so non-temporal stores only pay the
cost of flushing full lines. All six entry points are judged slower in that tier, spanning
0.13–0.51x.

On the ceiling: byte accounting gives expected values of 32/24 = 1.333 for the `add` family and
24/16 = 1.500 for `scale` and `normalize`. How many floats an element moves can only be counted per
entry point. Measured values scatter widely, so a single point exceeding the ceiling by up to 20%
counts as noise.

### Measurements for `_mask_inplace`

Five densities by four sizes, 20 points in total. All are inconclusive, with gains between −0.10 and
+0.00.

The reason is in the traffic accounting. An in-place write touches a target array that has already
been read, so the cache line is already resident and an ordinary store triggers no read-for-
ownership in that situation. In-place saves 0 B per element, not 8. Its value is shape: it removes
an output buffer, not time.

### Reading these numbers

- Always report the working set alongside a speedup. At N = 1<<20 (24 MB) all three implementations
  converge to 1.0–1.5x, limited by this machine's 6–7 GB/s streaming bandwidth.
- Most of the gain comes from removing cross-TU call boundaries, not from intrinsics. Once the same
  math is inlined from the header, `-O3 -march=native` auto-vectorization nearly matches hand-
  written AVX2.
- This machine carries persistent external load, so the same code can vary 20–25% between runs.
  Never compare absolute values across fixtures; only ratios within one run are meaningful.
- The absolute speedup of the 3D and 4D batch entries over scalar is not measured. Do not quote it.
  This round measured only their shape difference: of 65 points, 5 passed, 8 were slower and 52
  were inconclusive, with readings scattered between 0.97 and 1.93.

The acceptance criteria are in `docs/acceptance-protocol.md`. The raw verdict table is
`docs/research/raw/ab_campaign_verdict.txt`, which every `make run` overwrites, so quote it together
with its timestamp.

## Building and testing

```bash
make                         # build/libflashvec.a
make test-all                # 30 + 57 + 1098 assertions + two C compatibility tests

make test-header             # strict compile settings + double inclusion
make test-inline             # default form vs FV_USE_EXTERN form, byte for byte
make test-fma                # still bit-identical for a -march=native -mfma consumer
make test-inline-effects     # survey of what inlines and what does not
make test-stream             # regression for _stream and _mask_inplace
make test-vec4               # 4D scalar and batch
make test-vec34-inline       # 3D/4D batch shape comparison
make test-refactor           # all of the above plus symbol list and disassembly survey
```

### Compiler flags

The default is `-O3 -fno-math-errno -ffp-contract=off`.

`-fno-math-errno` drops the guarantee that math functions set `errno` on failure. Here `sqrt` never
sees a negative argument, so `errno` is never set anyway; adding the flag is what lets loops
containing `sqrtf` vectorize, with no change to a single bit.

`-ffp-contract=off` forbids contracting `a + b*c` into an FMA. Contraction would make the same
expression differ by 1 ulp across translation units or targets. The header adds `#pragma GCC
optimize("fp-contract=off")` as a second line of defence.

`-march=native`, `-mavx2` and `-ffast-math` are not used. The default build stays on the x86-64
baseline. AVX2 is enabled only inside the header, through function-level `target("avx2")` attributes
plus runtime dispatch on `__builtin_cpu_supports`. Machines without AVX2, and every tail, take the
scalar path with identical results.

## Design notes

### C written with a C++ compiler

The header is a pure C interface with `extern "C"` guards. The implementation is compiled as C++ but
only uses the C language subset. Introducing C++ features later requires neither a toolchain change
nor an ABI change.

### Batch entry points live in the header

The cross-TU call boundary is a measured first-order bottleneck, and a caller linking a static
library cannot remove it. So the implementation is written once, in the header, while
`src/flashvec_batch.cpp` keeps same-named `extern` symbols for the `FV_USE_EXTERN` form and for
already-compiled callers.

The cost is that the header conditionally includes `<immintrin.h>` on supported targets. With a pure
C compiler, or with `FV_USE_EXTERN` defined, it remains a declaration-only header.

### The pragma in the header blocks inlining

`target("avx2")` cannot stop the `-mfma` a caller brings on the command line. Measured under
`-march=native -mfma -ffp-contract=fast`, kernels carrying the target attribute still emit `vfmadd`.
Entry points at risk of contraction therefore have to be wrapped in the pragma, and GCC will not
inline a function with different optimization options into its caller.

The consequence is that these entry points (`lerp`, `madd`, `dot`, `len`, `normalize`) cost one
extra translation-unit-local call per batch. The other eight inline completely. `make test-inline-
effects` pins this contract down as a reproducible assertion.

## Not done yet

- `_fast` variants, for example a `rsqrtps` normalize with one Newton iteration, about 2 ulp of
  error and roughly 10% faster. These must be new functions and must not change existing
  semantics.
- Approximate trigonometry. An `atan` polynomial runs about 2.7x faster than glibc, but the input
  domain has to be documented.
- Matrices, quaternions, angles and rotations.
- Absolute speedup measurements for the 3D and 4D batch entries.

Already done, though an earlier version of this file listed them as future work: 3D batch entry
points via SoA, and `fv_vec4` with `_x2` AoS and `_soa` reduction batches.

## Documentation

| Document | Contents |
|---|---|
| [`docs/acceptance-protocol.md`](docs/acceptance-protocol.md) | The five acceptance criteria, the A/A floor, and six directions that measurements have already rejected |
| [`docs/optimization-report.md`](docs/optimization-report.md) | Full optimization analysis with evidence |
| [`docs/implementation-log.md`](docs/implementation-log.md) | What was landed, in what order, with the acceptance evidence for each step |
| [`bench/ab/README.md`](bench/ab/README.md) | The same-process A/B harness: how it works and where it lies |
| [`verify/README.md`](verify/README.md) | Independent verification and benchmarks |
| [`docs/research/`](docs/research/README.md) | Research notes (Chinese, kept for provenance) |

Chinese translations live beside each document under `zh-CN/`. They may lag the English
originals; where they disagree, the English version is authoritative.

## License

Apache License 2.0. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).