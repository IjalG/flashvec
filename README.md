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
| `fv_quat` | 16 B | Rotation; identity (0,0,0,1); same layout as `fv_vec4`, but a distinct type |
| `fv_mat3` | 36 B | Row-major 3x3; the same convention as `fv_mat4`, the opposite of GLM |

## Scalar API

Each vector type has 18 core functions (`fv_vec3` has one extra, `cross`):

`make` `add` `sub` `neg` `scale` `lerp` `min` `max` `clamp` `madd`
`dot` `len_sqr` `len` `dist_sqr` `dist` `normalize` `project` `near`

On top of those, all three types have `addmul` (`x + v*s`), and `fv_vec2` alone has `maddv`
(`a*b + c` with a vector third operand).

The bounding box adds `fv_aabb2_merge` and `fv_aabb2_overlaps`.

The default form defines every scalar body as `static inline` in the header, under the same
`fp-contract=off` pragma as the batch kernels; `FV_USE_EXTERN` switches to the library symbols, and
the two forms are bit-identical (`make test-scalar-parity`). A mixed scalar loop measured with
`make bench-scalar` is **1.70-1.82x** faster in the default form than in the `FV_USE_EXTERN` form,
and a consumer compiled with `-mfma -ffp-contract=fast` still emits **0 `vfmadd`**. The default
form is `static inline`; use `FV_USE_EXTERN` when you need a stable library symbol / function address,
or an ABI boundary. `normalize` keeps its scalar `sqrt`/`div` and zero branch even inline; the
`_fast` family is the opt-in path for that one.

For the exact form, compile with **`-ffp-contract=off`** (the library's documented default). Then
every scalar body inlines: the local 115-item GLM re-run has the scalar median at **1.00x** with
that flag, versus 0.86x (`-O2`) / 0.82x (`-O3n`) under the report's `-ffp-contract=fast` tiers. The difference is the
contract: a body containing `a + b*c` cannot be forced inline under a caller that allows FMA without
letting the caller contract it, so it stays out-of-line; the flag makes the caller's optimization
node match the header's `fp-contract=off` node.

`min`, `max` and `clamp` follow `(a < b) ? a : b` for NaN, not IEEE `fminf`. This is stated in the
header.

## Batch API

52 vector entry points, grouped by memory layout, plus 4 matrix transform entry points. The
number in a name is how many elements fit in one 256-bit register: four `fv_vec2` (`_x4`) or
two `fv_vec4` (`_x2`). Both fill the register completely, so no lane is wasted.

The P6 physics group adds four SoA entries on top of that (quaternion multiply, normalize and
rotate-vector, plus the 3x3 transform), and P7 adds two mat3 batch entries (multiply / inverse) and
four `_fast` normalize entries. They are documented under [Physics primitives](#physics-primitives).

19 of the 52 are `_stream` variants of an operation already in the list, differing only in the
store instruction; see [Streaming stores](#streaming-stores).

The matrix transforms, `fv_mat4_transform_points` and `fv_mat4_transform_dirs`, each with a
`_stream` form, are documented under [Matrix transform](#matrix-transform) below.

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

### Fused, in-place, three operands

A physics step keeps needing the same three expressions, and writing each one as two library
calls pays a second pass over memory.

```c
/* out = normalize(a * s + b) in one pass; bit-identical to scale -> add -> normalize */
void fv_vec2_madd_normalize_soa(const float* ax, const float* ay,
                                const float* bx, const float* by,
                                fv_real s, float* ox, float* oy, size_t n);

/* x += v * s, in place */
void fv_vec2_addmul_x4 (fv_vec2* x, const fv_vec2* v, fv_real s, size_t n);
void fv_vec2_addmul_soa(float* x, float* y, const float* vx, const float* vy,
                        fv_real s, size_t n);
void fv_vec3_addmul_soa(float* x, float* y, float* z,
                        const float* vx, const float* vy, const float* vz,
                        fv_real s, size_t n);
void fv_vec4_addmul_x2 (fv_vec4* x, const fv_vec4* v, fv_real s, size_t n);

/* o = a * b + c; three vector operands, where madd's third operand is a scalar */
void fv_vec2_maddv_x4(const fv_vec2* a, const fv_vec2* b, const fv_vec2* c, fv_vec2* o, size_t n);
```

The per-value forms are `fv_vec2_addmul`, `fv_vec3_addmul`, `fv_vec4_addmul` and
`fv_vec2_maddv`; every batch entry above is bit-identical to its per-value counterpart.

The in-place write produces no RFO, because the store lands on the cache line the read has just
brought in. That is why the in-place family is the one with a ceiling above 2x: two passes cost
56 B/element, in place costs 24 B/element.

`x` is read and written, so `v` must not overlap it. `n == 0` touches no pointer.

The fused pipeline and `maddv` write a separate output, so both have a `_stream` twin. The
in-place `addmul` group has none on purpose: that store already lands on a line the core owns.

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

There are 21 `_stream` variants: 19 of the vector entries, plus both matrix transforms. The
vector ones cover `add` / `sub` / `scale` / `lerp` / `madd` in `_x4` and `_x2` form, the SoA forms
of `add` / `scale` / `normalize` (and of the fused pipeline) for `fv_vec2`, `add` / `scale` /
`normalize` for `fv_vec3`, `normalize_soa` for `fv_vec4`, and `maddv_x4`.

They write results with non-temporal stores, bypassing the cache and skipping the read-for-
ownership traffic. **Use them only when the working set is at least L3 — 4 MiB per CCX on this
machine — and the result is not read back soon.** Below that boundary they are not slightly
slower but **2–8 times slower**: measured inside L2 at n=8192, the ratio against the ordinary
store is **0.13–0.37**. Beyond L3 the same round measures **1.48–1.60x**.

The library will never switch to them on its own; the choice belongs to the caller.

## Physics primitives

The domain types on top of the vector layer: rotation (`fv_quat`), 3x3 linear algebra (`fv_mat3`),
and the general 4x4 operations on the existing `fv_mat4`. The quaternion is a separate type because it
is not a `fv_vec4` -- its product is a non-commutative composition, its zero is not a rotation, and its
batch form needs four component arrays rather than the vector layout.

### Types

| Type | Size | Notes |
|---|---|---|
| `fv_quat` | 16 B | `x, y, z, w`; identity `(0,0,0,1)`; deliberately layout-compatible with `fv_vec4`, but a distinct type |
| `fv_mat3` | 36 B | Row-major `m[row*3+col]`; the same convention as `fv_mat4` (the opposite of GLM/OpenGL); compact 9-float layout, no padding |

A **unit** quaternion `q = (x,y,z,w)` rotates `v` by `v' = q (v,0) q*`. `fv_quat_mul(a,b)` is the
Hamilton product `a*b`, so it applies `b` first and `a` second; `fv_mat3_mul(a,b)` has the same
composition order. `fv_quat_rotate_vec3` and `fv_mat3_from_quat` assume a unit quaternion and do not
normalize. `fv_quat_normalize` maps length <= `FV_EPS` to the identity -- the one deliberate difference
from `fv_vec*_normalize`, which returns the zero vector: a zero quaternion is not a rotation.

`fv_quat_inverse` is the general inverse `conjugate(q) / |q|^2` (for a unit q it equals the conjugate).
A length <= `FV_EPS` returns the **zero** quaternion, not the identity: the inverse does not exist, and
zero is the signal. `fv_mat3_inverse` is the adjugate / determinant; a matrix is treated as singular when
`|det| <= FV_EPS * max|entry|^3` (a scale-relative test, so `1e-3 * I` stays invertible), and the value
form then returns the zero matrix. The `_checked` forms return 1 / 0 and write the zero value on failure;
`fv_mat3_inverse_checked` accepts `out == m`.

`fv_quat_nlerp` / `fv_quat_slerp` use the shortest-arc convention (a negative dot negates the second
quaternion), normalize the result, and do not clamp `t`; slerp falls back to nlerp near dot = 1.
`fv_quat_from_axis_angle` normalizes the axis (zero axis -> identity), and `fv_quat_from_euler_xyz`
applies X first, then Y, then Z (`q = qz*qy*qx`). These use `sin/cos`, so they are accuracy-bounded, not
bit-reproducible across libm implementations.

`fv_quat_from_mat3` is Shepperd's branch-on-largest-diagonal conversion; it assumes a rotation matrix
and returns a normalized quaternion. `fv_mat4_mul` / `fv_mat4_inverse` are the **general** 4x4
operations (not affine-only); the inverse uses the adjugate/determinant closed form and the same
scale-relative singularity rule (`|det| <= FV_EPS * max|entry|^4`), with the value form returning the
zero matrix and `fv_mat4_inverse_checked` returning 1 / 0 (and allowing `out == m`).

### Scalar and batch API

Scalar: `fv_quat_make / add / sub / neg / scale / dot / len_sqr / len / normalize / conjugate /
inverse / mul / rotate_vec3 / nlerp / slerp / from_axis_angle / from_euler_xyz`;
`fv_mat3_from_rows / transform_vec3 / mul / transpose / from_quat / from_mat3 / determinant / inverse /
inverse_checked`; `fv_mat4_mul / fv_mat4_inverse / fv_mat4_inverse_checked`; `fv_atan_approx`.

Batch (SoA, 8 elements per 256-bit register; a quaternion is four arrays, a 3x3 is nine):

- `fv_quat_mul_soa` -- 4 + 4 input arrays, 4 outputs, `n`
- `fv_quat_normalize_soa` -- 4 input arrays, 4 outputs, `n`
- `fv_quat_rotate_vec3_soa` -- 4 quaternion + 3 vector inputs, 3 outputs, `n`
- `fv_mat3_transform_vec3_soa` -- 9 matrix + 3 vector inputs, 3 outputs, `n`
- `fv_mat3_mul_soa` -- 9 + 9 matrix inputs, 9 outputs, `n`
- `fv_mat3_inverse_soa` -- 9 matrix inputs, 9 outputs, `n` (singular lanes write nine zeros)
- `fv_vec2/vec3/vec4_normalize_soa_fast`, `fv_quat_normalize_soa_fast` -- approximate `rsqrtps` family

The first six are **bit-identical** to their per-value counterparts for every `n` and every tail length,
inside the same `fp-contract=off` pragma: `make test-physics-fma` builds the consumer with
`-march=native -mfma -ffp-contract=fast` and the batch output still matches the scalar one
(`vfmadd = 0`). The usual NaN-payload exception applies; see the guarantees section.

Implementation note: quaternion multiply, rotate-vector and mat3 transform are **plain C** inside a
`target("avx2")` function -- GCC vectorizes them at 256 bits with zero shuffles (measured 40 / 40 / 21
AVX2 instructions including loads and stores, 0 `vpermps`/`vshufps`). Quaternion normalize was the one
loop the compiler refused to vectorize (it kept `sqrtss`/`divss` and a branch), so it has an explicit
AVX2 kernel with a blend to the identity. The mat3 batch mul/inverse now have their own AVX2 kernels:
mul keeps the nine B columns in registers and computes three output rows at a time; inverse does one
`vdivps(1/det)` per eight elements (instead of nine scalar divisions) and blends zero for singular
lanes. Both keep the scalar tail and stay bit-identical to the per-value API. Local
`make bench-mat3` (default inline scalar baseline, post-P0): mul 4.9x at n=8192 and 3.1-3.6x at 1M;
inverse 9.7x and 2.7-3.5x. Against the old `FV_USE_EXTERN` scalar form the same kernels measured 13.3x /
3.5x and 14.7x / 3.2x. Against GLM's per-element loop (local 115-item re-run): `m3_mul_soa`
**5.90x** at `-O2` and 4.18x at `-O3n` (8192), 3.88x/1.27x at 2M; `m3_inverse_soa` **4.59x/4.59x**
at 8192, 2.78x/2.62x at 2M. The report's 0.48x / 0.68x negatives are gone.
Hand-written intrinsics are limited to the places where the compiler failed or the layout demands
it: mat4 transform, quaternion normalize, mat3 mul/inverse, the fused normalize kernels, and the
optional AVX-512 variants. The rule remains: hand-write only after the compiler fails or the
gather/scatter shape makes auto-vectorization impossible.

### Approximate entries

The `_fast` normalize family uses `rsqrtps` + one Newton-Raphson step: relative error about 1-2 ulp on
AVX2, exact `1/sqrtf` on the scalar fallback, and the same zero guard as the exact entry (zeros, or the
quaternion identity). It is a separate family and is **not** bit-identical to the exact entries. The
measured gain is about 1.13x for `vec2` and 1.05x for `vec3` at n=8192; at n=1M it is 1.03-1.08x, i.e.
inside the run's A/A spread, so read that tier as inconclusive. `fv_atan_approx` is the 5-term
polynomial on `[-1,1]` (max abs error ~1.2e-5, NaN outside the domain); it is the intended fast path
for angle comparisons, not a libm replacement.

`fv_vec2_madd_normalize_soa_fast` is the fused form of the same approximate path
(`normalize(a*s + b)`); at n=8192 it is 1.07x over the exact fused entry, and at n=1M it is inside the
A/A spread. Against GLM, the exact fused entry is **6.55x** (`-O2`) / 7.50x (`-O3n`) at 8192 and
1.09x/1.16x at 2M; the report's 0.79x cell was its AVX-512 cloud and is not reproduced here. The exact fused entry also has an **optional AVX-512 backend**: 16-lane kernels are
compiled under `target("avx512f")` and selected at run time by `__builtin_cpu_supports("avx512f")`;
the AVX2 path is unchanged and `FV_I_NO_AVX512` removes the backend. This machine has no AVX-512, so
the codegen is checked (`vrsqrt14ps` / `zmm` present) but the runtime path still needs validation on
AVX-512 hardware.

### Preliminary measurements

`make bench-physics`, same-process interleaved A/B, 11 paired rounds, `taskset -c 0`, medians of two
runs. These tables predate P0: their baseline is the **`FV_USE_EXTERN`** scalar form (a cross-TU call).
Since P0 the default scalar form is header-inline, so a batch-vs-default-scalar ratio is smaller;
the table stays valid for the old/extern baseline. The same-round A/A floor was
0.99-1.00 at n=8192 and 0.93-1.07 at n=1M.

| Entry | n = 8192 (L2) | n = 1M (32-60 MiB working set) |
|---|---|---|
| `fv_quat_mul_soa` | 18.0-19.5x | 2.39-2.51x |
| `fv_quat_normalize_soa` | 41.9-42.3x | 4.56-5.13x |
| `fv_quat_rotate_vec3_soa` | 25.8x | 3.39-3.57x |
| `fv_mat3_transform_vec3_soa` | 24.7-25.1x | 3.33-3.80x |

At n=8192 the batch entries run at 0.4-0.55 ns/element against 8-20 ns for the scalar calls; most of
the gain is the removed call boundary, exactly as the library's measured thesis predicts. At n=1M the
working set is DRAM-bound and the ratios collapse toward the memory ceiling. These are preliminary,
not the five-criterion acceptance of `docs/acceptance-protocol.md`; quote them with the harness.

## Guarantees

### Bit-for-bit identical

The following entry points are guaranteed to match the corresponding scalar API bit for bit, so
`memcmp` is a valid test:

- the `_x4`, `_x2`, `_soa`, `_mask` and `_idx` variants of `add`, `sub`, `scale`, `lerp` and `madd`
- `len_soa`
- `normalize_soa`
- the fused pipeline `madd_normalize_soa` against `scale_soa` then `add_soa` then `normalize_soa`
- `addmul` and `maddv` against their per-value counterparts
- `fv_quat_mul_soa` / `fv_quat_normalize_soa` / `fv_quat_rotate_vec3_soa` against the per-value quaternion API
- `fv_mat3_transform_vec3_soa` / `fv_mat3_mul_soa` / `fv_mat3_inverse_soa` against their per-value counterparts
- every `_stream` function against its non-stream counterpart

The implementation uses only `vaddps`, `vsubps`, `vmulps`, `vdivps` and `vsqrtps`. No FMA, no
reciprocal approximations, and reduction order matches the scalar path. The `_fast` entries and
`fv_atan_approx` are a **separate, explicitly approximate family**: they are not covered by the
bit-for-bit list above and must not be mixed with it.

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

### Matrix transform

```c
typedef struct { float m[16]; } fv_mat4;   /* ROW-major: m[row*4 + col] */

fv_mat4 fv_mat4_identity(void);

void fv_mat4_transform_points(const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n);
void fv_mat4_transform_dirs  (const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n);
```

`transform_points` treats `in[i].w` as 1, so translation applies. `transform_dirs` treats it as
0, so translation is ignored and `out[i].w` is 0.

**This is row-major, which is the opposite of GLM and OpenGL.** A GLM matrix must be transposed
before it goes in.

This is the one entry point where hand-written SIMD actually beats the compiler. The compiler
does vectorize a plain AoS loop, but it emits 440 instructions with 14 `vpermps`; the hand-written
kernel is 107 instructions with 3. At L1/L2-resident sizes that is **1.81x over a scalar loop
compiled at this library's default flags** (7/7 runs in 1.81–1.84 at n=8192), and 2.25x if the
scalar baseline is compiled with `-march=native` instead. Quote the flags with the number: the
baseline's compile tier moves the figure by 30%. Above about 131072 elements it converges to 1.0x,
because at 48 B per element it is DRAM-bound like everything else.

There is no SoA variant. Measured, it would be 1.01x: the compiler already handles SoA perfectly
(published as `bench/tfprobe/probe.cpp`, which anyone can re-run).

#### The one place bit-identical does not reach

If one arithmetic instruction gets **two NaNs with different payloads**, which NaN survives depends
on that instruction's operand order. The compiler picks that order, and it picks differently in
different builds: measured, the same entry point returns different NaN payloads when the caller
compiles it with `-march=native` than when it links the library build. So the scalar tail, the
vector kernel, and the same kernel compiled another way can each return a different — all legal —
NaN bit pattern for one element, meaning a result can change with `n` and with your compiler flags.

Every other input is bit-identical, for every `n` and every tail length. If the corpus contains
at most one NaN bit pattern, the guarantee holds for those too.

## Performance

Test machine: AMD Ryzen 7 5700U, 8 cores / 16 threads, AVX2 + FMA, no AVX-512.

### Cache boundaries

L3 is **4 MiB per CCX**, 8 MiB across both CCXs. A single thread only reaches its own 4 MiB. This is
visible in `/sys/devices/system/cpu/cpu0/cache/index3`: `size=4096K`, `shared_cpu_list=0-7`.

At 24 B per element for `add`, the knee is at 4 MiB ÷ 24 B ≈ **1.75e5 live elements**. Measurements
agree: 2^17 (3 MB working set) runs at 38–57 GB/s, while 2^18 (6 MB) drops to 5.6–8.8 GB/s.

### The machine ceiling

Every figure below is a fraction of what this machine can move, so the ceiling comes first. It is
`bench/membw.c`, one thread, taskset-pinned, best of 5:

| Working set | triad | triad16 (with RFO) | read | memcpy |
|---|---|---|---|---|
| 96 MB | 7.88 GB/s | 10.51 GB/s | 13.28 GB/s | 12.93 GB/s |
| 48 MB | 9.16 GB/s | 12.21 GB/s | 14.31 GB/s | 14.31 GB/s |

`triad` counts logical bytes, 12 B/element. `triad16` is the same run with the write's
read-for-ownership allocated in, 16 B/element.

The element-wise baseline: `fv_vec2_add_x4` on a 96 MB working set runs at **3.33 ns/element**,
about **9.6-10.3 GB/s** against a machine triad16 ceiling of 10.51 GB/s at the same 96 MB working set -- **91-98%** of it. There is no 10x left at
that size; what is left is the traffic itself.

### Speedup over scalar calls

The table below uses the historical **`FV_USE_EXTERN`** baseline: calling the library's scalar
function per element, a cross-TU call that cannot be inlined. Since P0 the default scalar form is
header-inline, so a loop written against `fv_vec2_add` etc. now auto-vectorizes too; the table is
kept as the old-baseline record, not as the default-form ratio. `make bench-scalar` measures the two
scalar forms directly (1.70-1.82x for the mixed workload).

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

### Measurements for the fused and in-place entries

The same harness and conditions: single thread, taskset-pinned, the default tier
`-O3 -fno-math-errno -ffp-contract=off` **without `-march`**, alternating A/B inside one
process, the same-round A/A floor beside every row.

| Change | n=8192 | n=1M | n=8M | A/A floor |
|---|---|---|---|---|
| `_stream` vs an ordinary store | — | 1.48–1.60x | 1.54–1.60x | 0.97–1.07 |
| fused `normalize(a*s+b)` vs three calls | 1.73x | 2.43&ndash;2.55x | 2.53&ndash;2.63x | — |
| the same, plus `_stream` | 1.68x | **3.53–3.68x** | **3.72–3.78x** | — |
| in-place `x += v*s` vs `scale` then `add` | 1.06–2.18x | 2.05–2.17x | 1.92–2.03x | 0.97–1.06 |
| `maddv` vs `mul` then `add` | 1.60–1.67x | 1.61–1.79x | 1.59–1.76x | ~1.00 |
| `maddv_stream` vs `mul` then `add` | 1.16–1.37x | 2.16–2.23x | **2.34–2.41x** | ~1.00 |

The `_stream` row is its own round. The campaign in the previous section read 1.29–1.58x on the
same two sizes and 0.13–0.51x inside L2, against 1.48–1.60x and 0.13–0.37 here; that is the same
conclusion within this machine's run-to-run spread.

Three things this table does not say:

- **`_stream` is slower at small sizes, and not slightly: 2–8 times slower** (0.13–0.37 at
  n=8192). An NT store bypasses the cache, so inside L2 it has no RFO to save and loses the cache
  hit as well. Use it only when the working set is at least L3, 4 MiB per CCX.
- **`maddv_stream` at n=8192 is faster relative and slower absolute**: 1.16–1.37x against the
  two-pass workaround, but 0.99 ns/element against 0.58 for the non-stream form. Inside L2 the
  non-stream form is the one to call.
- **The byte account is a ceiling, and 20% above it is noise.** `add` / `sub` / `lerp` /
  `madd` top out at 32/24 = 1.33x, `scale` / `normalize` / the matrix transform at 1.50x,
  in-place `addmul` at 2.33x (56/24), `maddv` at 1.60x and `maddv_stream` at 2.00x. Do not
  claim a measurement is above the theoretical limit on a single point. The `_stream` row reads
  1.48–1.60x, slightly above the 1.33x ceiling for the `add` family; that excess is real and has
  a reason — an NT store also keeps the output from polluting the cache, which the byte account
  does not count.

### Measurements for `_mask_inplace`

Five densities by four sizes, 20 points in total. All are inconclusive, with gains between −0.10 and
+0.00.

The reason is in the traffic accounting. An in-place write touches a target array that has already
been read, so the cache line is already resident and an ordinary store triggers no read-for-
ownership in that situation. In-place saves 0 B per element, not 8. Its value is shape: it removes
an output buffer, not time.

### 3D and 4D absolute speedup

`make bench-vec34`, same-process interleaved A/B, 11 paired rounds, `taskset -c 0`, one run. Baseline
is the per-element library call; the same-round A/A control sat at 0.97-1.00 at n=8192 and 0.98-1.06
at n=1M (the `vec4_add` 8192 A/A spread was wider, 0.15-1.17, but the A/B ratio is 35x, far outside it).

| Entry | n = 8192 | n = 1M |
|---|---|---|
| `fv_vec3_add_soa` | 8.8x | 2.16x |
| `fv_vec3_scale_soa` | 13.9x | 1.43x |
| `fv_vec3_normalize_soa` | 10.1x | 2.11x |
| `fv_vec4_add_x2` | 35.2x | 2.53x |
| `fv_vec4_scale_x2` | 25.3x | 2.05x |
| `fv_vec4_normalize_soa` | 38.8x | 4.34x |

These are absolute speedups over the **pre-P0 `FV_USE_EXTERN` scalar form**, so they answer the
question the earlier document left open; since P0 the default scalar form is inline, so a fresh
comparison against it is narrower. At n=1M the working set is DRAM-bound and the ratios converge
toward the memory ceiling.

### `_fast` normalize

`make bench-fast`, same protocol. A = `_fast`, B = exact SoA; speedup = exact/fast:

| Entry | n = 8192 | n = 1M |
|---|---|---|
| `fv_vec2_normalize_soa_fast` | 1.125x | 1.03x |
| `fv_vec3_normalize_soa_fast` | 1.054x | 1.04x |

At n=8192 the A/A control was 0.995-0.997; at n=1M it was 0.89-1.08, so the 1M gains are inside the
noise band and should be read as "no worse", not as a measured speedup.

### Reading these numbers

- Always report the working set alongside a speedup. At N = 1<<20 (24 MB) the three implementation
  tiers converge to 1.0–1.5x, limited by this machine's 6–7 GB/s streaming bandwidth. Historical
  batch-vs-scalar ratios are larger because that baseline was the `FV_USE_EXTERN` cross-TU call;
  since P0 the default scalar form is header-inline, so quote the baseline with the number.
- The historical gain came mostly from removing the cross-TU call boundary, not from intrinsics.
  Once the same math is inlined from the header, `-O3 -march=native` auto-vectorization nearly matches
  hand-written AVX2. P0 applies that lesson to the scalar API itself.
- This machine carries persistent external load, so the same code can vary 20–25% between runs.
  Never compare absolute values across fixtures; only ratios within one run are meaningful.
- The 3D / 4D absolute speedups are now measured by `make bench-vec34` (table above); like every
  ratio here, quote them with the harness and the working set, not as a single-point multiplier.
- AoS <-> SoA conversion is measured by `make bench-pack`: the round trip costs 0.28-0.32 ns/element
  at n=8192 and 4.2-4.4 ns/element at n=1M. With an AoS input, `pack + add_soa + unpack` is 2.5x
  slower than `fv_vec2_add_x4` at 8192 and 3.1x slower at 1M; for normalize the conversion still pays
  3.9x at 8192 but loses 2.3x at 1M. SoA entries are for data that is already SoA, or
  cache-resident plus expensive enough to amortize the conversion.

The acceptance criteria are in `docs/acceptance-protocol.md`. The raw verdict table is
`docs/research/raw/ab_campaign_verdict.txt`, which every `make run` overwrites, so quote it together
with its timestamp.

### GLM re-run (local)

The external 115-item GLM matrix was re-run on this Ryzen 5700U (gcc 12.4, GLM 0.9.9.8; the report
used GLM 1.1.0 on a cloud Xeon with AVX-512). Key numbers, default report flags unless noted:

| Tier | P median | B median (8192) | B geomean |
|---|---|---|---|
| `-O2` | 0.86x | **2.94x** | 1.58x |
| `-O3 -march=native` | 0.82x | 1.23x | 1.39x |
| `-O2 -ffp-contract=off` | **1.00x** | 2.64x | 1.56x |
| `-O3n -ffp-contract=off` | **1.00x** | 1.11x | 1.39x |

The report's 0.40x scalar median is gone with the documented exact flag; the remaining default-flag
gap is the FMA/inline contract trade described above. The three batch negatives are now wins:
`m3_mul_soa` 5.90x at `-O2` (report 0.48x), `m3_inverse_soa` 4.59x (report 0.68-0.98x), and
`v2_madd_normalize_soa` 6.55x (report 0.79x). `_stream` flips sign on this machine (1.1-1.7x here,
0.47-1.03x on the report's cloud), which is exactly the portability warning. The full table, the
GLM-version caveat and the remaining scalar losses (`mat3/mat4 inverse`, `project`, `slerp`,
`from_mat3`) are in [`docs/glm-rebench.md`](docs/glm-rebench.md).

### Batch sizing and streaming portability

- The batch API targets **n >= 4**. The external GLM report measured n=1 at roughly tie, n=2 at 0.59x,
  and n=4 and up at 2x or better; below 4 elements the call/guard overhead is the whole cost. Use the
  scalar API for one or two elements.
- `_stream` is not portable across machines in sign. It pays only when the working set is beyond L3
  and the machine is RFO-limited; on a bandwidth-rich, read-bound machine it loses (the report's
  machine: 0.47-1.03x). The library never switches automatically; measure with `bench/membw.c` first.
- Compile tier changes GLM's auto-vectorization far more than flashvec's kernels: quote a ratio only
  with `-O2` vs `-O3 -march=native` fixed, and do not compare single points across tiers.

## Building and testing

```bash
make                         # build/libflashvec.a
make test-all                # full suite: 1.97M checks (tests/count_checks.sh reproduces the total)

make test-header             # strict compile settings + double inclusion
make test-inline             # default form vs FV_USE_EXTERN form, byte for byte
make test-scalar-parity      # scalar API: header-inline vs library symbols, byte for byte
make test-scalar-parity-fma  # scalar API: -mfma consumer still bit-identical, 0 vfmadd
make test-fma                # batch entries: still bit-identical for a -march=native -mfma consumer
make test-inline-effects     # survey of what inlines and what does not
make test-stream             # regression for _stream and _mask_inplace
make test-physics            # quat + mat3, all forms (C/C++ inline/extern/pure-C/FMA)
make test-quat               # quaternion scalar semantics + SoA batch, bit-identical
make test-mat3               # 3x3 scalar semantics + SoA batch, bit-identical
make test-physics-c          # the same contract from a plain C11 translation unit
make test-mat4-ops           # general 4x4 multiply / inverse
make test-fast               # _fast normalize accuracy + fv_atan_approx
make bench-physics           # per-element calls vs SoA batch, same-process A/B
make bench-fast              # exact vs _fast normalize
make bench-vec34             # 3D / 4D absolute speedup over scalar calls
make bench-pack              # AoS <-> SoA conversion cost and the SoA decision rule
make bench-mat3              # mat3 AVX2 mul/inverse vs per-element scalar calls
make bench-scalar            # scalar inline vs FV_USE_EXTERN library calls
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
optimize("fp-contract=off")` as a second line of defence; since P0 it covers the inline scalar
definitions as well, which is why a `-mfma -ffp-contract=fast` consumer still emits no `vfmadd`
for them.

`-march=native`, `-mavx2` and `-ffast-math` are not used. The default build stays on the x86-64
baseline. `-ffast-math` is not merely unnecessary: measured, it gives no speedup anywhere, makes
`normalize_soa` 15–17% slower at L1/L2 sizes, and changes 32.5% of output bit patterns. See
`docs/acceptance-protocol.md` §9 item 7. AVX2 is enabled only inside the header, through function-level
`target("avx2")` attributes plus runtime dispatch on `__builtin_cpu_supports`. The optional AVX-512
backend uses `target("avx512f")` and `__builtin_cpu_supports("avx512f")`; `FV_I_NO_AVX512` removes it.
Machines without AVX2, and every tail, take the scalar path with identical results.

## Design notes

### C written with a C++ compiler

The header is a pure C interface with `extern "C"` guards. The implementation is compiled as C++ but
only uses the C language subset. Introducing C++ features later requires neither a toolchain change
nor an ABI change.

### Entry points live in the header

The cross-TU call boundary was a measured first-order bottleneck that a caller linking a static
library could not remove; moving the bodies into the header is the fix. The batch kernels and, since
P0, the scalar bodies are written in the
header as `static inline`; `src/*.cpp` keeps same-named `extern` symbols for the `FV_USE_EXTERN`
form, for `FV_BUILD`, and for already-compiled callers. The scalar and batch forms are compared
bit-for-bit by `tests/test_scalar_parity.c` and `tests/test_parity.c`.

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

## Scope

The library now covers the vector types and their scalar / batch APIs, quaternions, 3x3 and general
4x4 matrices, quaternion interpolation and angle constructors, the `_fast` approximate normalize
family, `fv_atan_approx`, and the optional AVX-512 path for the fused madd-normalize entry. There is
no pending backlog list in this document.

What is deliberately **outside** the scope, not pending: threading / scheduling and other engine
facilities, allocation or arena ownership (the library never allocates), a double-precision batch
family, and a `_fast` variant of every entry (only normalize and the fused madd-normalize have one).
The AVX-512 backend is opt-in/runtime-dispatched; on this AVX2-only machine its codegen is checked but
its runtime path awaits AVX-512 hardware validation. Those scope choices would change the
value-semantics, no-allocation contract rather than complete it.

## Documentation

| Document | Contents |
|---|---|
| [`docs/acceptance-protocol.md`](docs/acceptance-protocol.md) | The five acceptance criteria, the A/A floor, and eight directions that measurements have already settled |
| [`docs/optimization-report.md`](docs/optimization-report.md) | Full optimization analysis with evidence |
| [`docs/implementation-log.md`](docs/implementation-log.md) | What was landed, in what order, with the acceptance evidence for each step |
| [`docs/glm-rebench.md`](docs/glm-rebench.md) | Local re-run of the external 115-item GLM matrix, with the machine/GLM-version caveats |
| [`bench/ab/README.md`](bench/ab/README.md) | The same-process A/B harness: how it works and where it lies |
| [`verify/README.md`](verify/README.md) | Independent verification and benchmarks |
| [`docs/research/`](docs/research/README.md) | Research notes (Chinese, kept for provenance) |

Chinese translations live beside each document under `zh-CN/`. They may lag the English
originals; where they disagree, the English version is authoritative.

## License

Apache License 2.0. See [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE).