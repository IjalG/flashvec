/*
 * flashvec -- a small vector foundation library
 *
 * Design conventions:
 *   - The header is a pure C (C11) interface; both C and C++ can include it directly;
 *   - The implementation files are compiled as C++, but they write only the C language subset, so when C++ features
 *     are wanted later there is no need to change toolchains and no need to change the ABI;
 *   - Every quantity is a small struct passed by value; the library allocates no memory and has no global state.
 *
 * Implementation form of the batch entry points (P1 refactor: from "in-library implementation" to "header inline implementation")
 *   - Default form: the P1 group's 52 batch entry points (44 vector entry points + 2 matrix point-transform entry points + the 6
 *     integrator entries of the in-place / maddv group) are
 *     defined as static inline in this header's "implementation of the batch entry points" section,
 *     so the caller's translation unit gets the implementation directly, with no call across a translation unit;
 *     the P6 physics group adds four more SoA entries (quaternion multiply / normalize / rotate-vector,
 *     mat3 transform), P7 adds the mat3 multiply / inverse pair and four _fast normalize entries, all in
 *     the same section;
 *   - The same-named extern symbols are still exported by the library (src/flashvec_batch.cpp), so old consumers still link;
 *   - Whether intrinsics appear does not depend on the language, only on "compiler + architecture": x86/x86_64 with GCC
 *     or Clang -- **a C translation unit and a C++ translation unit follow the same rule**.
 *     Everything else (MSVC and other compilers without the target attribute, non-x86 architectures) takes the "auto-vectorizable plain C"
 *     path: such a translation unit never expands a single line of immintrin.h / __m256 / _mm_*.
 *   - Having intrinsics does not mean they run: whether an AVX2 kernel runs is decided at run time by
 *     __builtin_cpu_supports("avx2") != 0; a machine without AVX2, and every tail, take the scalar path,
 *     where the scalar expressions are word-for-word identical to the per-element API in src/flashvec.cpp.
 *   - All four statements above can be verified directly (from the repository root):
 *       cc  -std=c11   -Iinclude -E tests/test_header_c.c    | grep -c immintrin   # > 0
 *       c++ -std=c++17 -Iinclude -E tests/test_header_cpp.cpp | grep -c immintrin  # > 0
 *       cc  -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -c tests/test_header_c.c  # 0 warnings
 *       nm -g --defined-only build/libflashvec.a | grep -c '_x4\|_x2\|_soa'        # 33
 *         (the two matrix point-transform symbol names contain no _x4/_x2/_soa, so they are not counted in this 33)
 *       make test-inline           # default form vs FV_USE_EXTERN vs plain C fallback, byte-for-byte comparison
 *       make test-inline-effects   # consumer at the default -O2: first section inlines fully, second calls once per batch
 *       make test-fma              # a consumer with -march=native -mfma -ffp-contract=fast is still bit-identical
 *
 * Responsibility boundary between the floating-point contract and the compiler options
 *   - **The half of that implementation section which carries contraction risk** is covered by #pragma GCC optimize("fp-contract=off") (GCC,
 *     push and pop paired), forbidding the contraction of a + b*s into an FMA. Its scope is "code generation for those
 *     functions inside this header"; it has no effect on the caller's own code or on the caller's compiler options;
 *   - target("avx2") cannot stop FMA coming from the caller's command line: measured with -march=native -mfma
 *     -ffp-contract=fast, vfmadd still appears inside a kernel carrying target("avx2").
 *     So in the consumer scenario **the pragma is the only line of defence**, not the target attribute;
 *   - Therefore: a caller that uses -ffast-math (including the reassociation and no-NaN assumptions of
 *     -fassociative-math / -ffinite-math-only), or that lets FMA contraction take effect on a compiler
 *     without this pragma (-march=native / -mfma with -ffp-contract=fast), **owns the bit-exactness promise**;
 *   - To have the result guaranteed by the library's own compiler options (-O3 -fno-math-errno -ffp-contract=off),
 *     define FV_USE_EXTERN and call the in-library implementation.
 */
#ifndef FLASHVEC_H
#define FLASHVEC_H

#include <stddef.h> /* size_t: used by the batch entry points */
#include <stdint.h> /* uint8_t / uint32_t: used by the mask and index variants */

/* These two headers are needed only when "this header provides the implementation": the include set of the FV_USE_EXTERN form is then exactly the pre-refactor one. */
#if !defined(FV_USE_EXTERN)
#include <math.h>   /* sqrtf: the scalar tail path of the batch entry points (including the plain C path) */
#include <string.h> /* memcpy: the _mask variants load 4 mask bytes */
#endif

/* ---------------- compiler / platform capability probes ----------------
 * FV_I_X86 (0/1)       : whether this is x86 / x86_64;
 * FV_I_INTRINSICS (0/1): whether the batch entry points may use AVX2 intrinsics -- x86 + GCC/Clang,
 *                        independent of the language (C and C++ translation units follow the same rule); when 0 the whole group takes the plain C path.
 * Both macros are implementation details of this header, kept so that they stay verifiable (tests/ asserts their values directly). */
#if defined(__x86_64__) || defined(__i386__)
#define FV_I_X86 1
#else
#define FV_I_X86 0
#endif

/* FV_I_NO_INTRINSICS: a self-test / degradation switch; defining it forces the plain C path.
 * Why it exists: on an x86 + GCC machine the plain C path would otherwise have no test coverage at all,
 * and it is exactly the path MSVC / non-x86 users really take (tests/test_parity.c's c_pure variant uses it for bit-level comparison). */
#if FV_I_X86 && (defined(__GNUC__) || defined(__clang__)) && !defined(FV_I_NO_INTRINSICS)
#define FV_I_INTRINSICS 1
#else
#define FV_I_INTRINSICS 0
#endif

/* FV_I_AVX512: the optional 512-bit backend. It is compiled whenever the compiler can emit
 * target("avx512f") functions (GCC/Clang on x86) and selected at run time by __builtin_cpu_supports;
 * define FV_I_NO_AVX512 to remove it from a build. The AVX2 path is the default and is unaffected. */
#if FV_I_INTRINSICS && !defined(FV_I_NO_AVX512)
#define FV_I_AVX512 1
#else
#define FV_I_AVX512 0
#endif

#if FV_I_INTRINSICS && !defined(FV_USE_EXTERN)
#include <immintrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Scalar precision: switching everything to double needs only this line (the sqrtf calls in the implementation change with it). */
typedef float fv_real;

/* Two zero-test thresholds of different dimension; do not mix them:
 *   FV_EPS     -- threshold on a linear quantity (length), used by normalization and the near comparison;
 *   FV_EPS_SQR -- threshold on a squared quantity = FV_EPS*FV_EPS; use it when you already hold a squared length or a dot product.
 * Comparing a squared length against FV_EPS loosens the threshold 1000x (the dimensional bug project once had; fixed in P0). */
#define FV_EPS ((fv_real)1e-6)
#define FV_EPS_SQR ((fv_real)1e-12)

typedef struct fv_vec2 {
    fv_real x;
    fv_real y;
} fv_vec2;

typedef struct fv_vec3 {
    fv_real x;
    fv_real y;
    fv_real z;
} fv_vec3;

/* 4D: 4 x fv_real = 16 B, **no padding bytes**, components stored contiguously as x, y, z, w.
 * The number 16 B is the precondition of the batch entry points fv_vec4_*_x2 (2 elements = 32 B = one 256-bit register),
 * so a compile-time assertion pins it down: changing this layout means changing the _x2 kernels in the same step. */
typedef struct fv_vec4 {
    fv_real x;
    fv_real y;
    fv_real z;
    fv_real w;
} fv_vec4;

#if defined(__cplusplus)
static_assert(sizeof(fv_vec4) == 16, "fv_vec4 must be exactly 16 B (no padding): _x2 batches rely on it");
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
_Static_assert(sizeof(fv_vec4) == 16, "fv_vec4 must be exactly 16 B (no padding): _x2 batches rely on it");
#endif

/* Axis-aligned bounding box: lo is the per-component lower bound, hi the per-component upper bound (lo <= hi is not required; see the clamp note). */
typedef struct fv_aabb2 {
    fv_vec2 lo;
    fv_vec2 hi;
} fv_aabb2;

/* 4x4 affine matrix, row-major: m[row * 4 + col]. 64 B, always passed by pointer --
 * by value would copy the whole matrix at every call.
 *
 * Affine precondition: the last row m[12..15] must be (0, 0, 0, 1). fv_mat4_transform_points /
 * _dirs read rows 0..2 only and write the w output as a constant (1 for points, 0 for dirs)
 * instead of evaluating 0*x + 0*y + 0*z + 1*w. That is one vaddps less, and -- the reason it is
 * a contract and not an optimization -- it is immune to 0 * inf = NaN, so a non-finite input
 * never contaminates w. A matrix whose last row is not (0, 0, 0, 1) is not rejected: its last
 * row is simply not read, and the output w is still 1 / 0. */
typedef struct fv_mat4 {
    fv_real m[16];
} fv_mat4;

#if defined(__cplusplus)
static_assert(sizeof(fv_mat4) == 64, "fv_mat4 must be exactly 64 B (16 x fv_real, no padding)");
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
_Static_assert(sizeof(fv_mat4) == 64, "fv_mat4 must be exactly 64 B (16 x fv_real, no padding)");
#endif

/* Identity matrix. Constant construction only, no float arithmetic -- so it is header-inline in
 * every form (default / FV_USE_EXTERN / FV_BUILD). There is nothing the two forms could
 * disagree about, hence no extern symbol and no bit-exactness question to answer. */
static inline fv_mat4 fv_mat4_identity(void)
{
    fv_mat4 r = { { (fv_real)1, (fv_real)0, (fv_real)0, (fv_real)0,
                    (fv_real)0, (fv_real)1, (fv_real)0, (fv_real)0,
                    (fv_real)0, (fv_real)0, (fv_real)1, (fv_real)0,
                    (fv_real)0, (fv_real)0, (fv_real)0, (fv_real)1 } };
    return r;
}

/* ---------------- quaternion (physics pose / rotation) ----------------
 *
 * 4 x fv_real = 16 B, no padding bytes, stored x, y, z, w. The identity rotation is (0, 0, 0, 1).
 *
 * Convention: a unit quaternion q = (x, y, z, w) with w = cos(theta/2) and (x, y, z) = axis * sin(theta/2)
 * rotates a vector by v' = q * (v, 0) * conjugate(q). fv_quat_mul(a, b) is the Hamilton product a*b, so it
 * applies b first and a second -- the same composition order as matrices ((A*B)*v == A*(B*v)).
 *
 * The scalar API does not normalize on entry: fv_quat_rotate_vec3 assumes a unit quaternion, exactly as a
 * rotation-matrix function assumes an orthonormal matrix. Call fv_quat_normalize when the input may have
 * drifted; the batch form fv_quat_rotate_vec3_soa carries the same precondition.
 *
 * Why a dedicated type instead of reusing fv_vec4: the field names and the multiplication semantics are the
 * point. fv_vec4's dot is a 4D inner product; a quaternion's product is a non-commutative composition, and
 * its "length" is the normalization denominator. Keeping them separate makes a wrong call a compile error,
 * not a silently different number. The layout is deliberately identical to fv_vec4 (16 B, x/y/z/w), so the
 * two can be reinterpreted in a pinch, but that is not part of the contract.
 */
typedef struct fv_quat {
    fv_real x;
    fv_real y;
    fv_real z;
    fv_real w;
} fv_quat;

#if defined(__cplusplus)
static_assert(sizeof(fv_quat) == 16, "fv_quat must be exactly 16 B (no padding): fv_vec4 layout compatibility");
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
_Static_assert(sizeof(fv_quat) == 16, "fv_quat must be exactly 16 B (no padding): fv_vec4 layout compatibility");
#endif

/* 3x3 matrix, **row-major**: m[row * 3 + col]. 36 B, 9 x fv_real, no padding.
 *
 * Row-major is the same convention as fv_mat4 (and the same warning applies: it is the opposite of GLM and
 * OpenGL, so a GLM matrix must be transposed before it goes in). The compact 9-float layout is deliberate:
 * padding to 48 B (4x3) or 64 B (4x4) is an AoS SIMD choice for the caller, while the library's batch form
 * takes the nine component arrays (SoA) and needs no padded struct. A caller that wants padded AoS can build
 * its own 12/16-float record; no library type forces that cost on everyone.
 */
typedef struct fv_mat3 {
    fv_real m[9];
} fv_mat3;

#if defined(__cplusplus)
static_assert(sizeof(fv_mat3) == 36, "fv_mat3 must be exactly 36 B (9 x fv_real, no padding)");
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
_Static_assert(sizeof(fv_mat3) == 36, "fv_mat3 must be exactly 36 B (9 x fv_real, no padding)");
#endif

/* Constant construction only (no float arithmetic), so these are header-inline in every form -- the same
 * reason as fv_mat4_identity, and the same absence of a bit-exactness question. */
static inline fv_quat fv_quat_identity(void)
{
    fv_quat r = { (fv_real)0, (fv_real)0, (fv_real)0, (fv_real)1 };
    return r;
}

static inline fv_mat3 fv_mat3_identity(void)
{
    fv_mat3 r = { { (fv_real)1, (fv_real)0, (fv_real)0,
                    (fv_real)0, (fv_real)1, (fv_real)0,
                    (fv_real)0, (fv_real)0, (fv_real)1 } };
    return r;
}

#if defined(FV_USE_EXTERN) || defined(FV_BUILD)
/* ---------------- 2D ---------------- */

fv_vec2 fv_vec2_make(fv_real x, fv_real y);
fv_vec2 fv_vec2_add(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_sub(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_neg(fv_vec2 v);
fv_vec2 fv_vec2_scale(fv_vec2 v, fv_real s);
fv_vec2 fv_vec2_lerp(fv_vec2 a, fv_vec2 b, fv_real t);

/* Per-component minimum / maximum. The definition is min(a,b) = (a < b) ? a : b, max(a,b) = (a > b) ? a : b,
 * not IEEE minNum/maxNum: with a NaN on one side the result is decided by the comparison
 * (min(NaN, x) == x, min(x, NaN) == NaN; when both are NaN the second one is returned). */
fv_vec2 fv_vec2_min(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_max(fv_vec2 a, fv_vec2 b);

/* Per-component clamp: x < lo ? lo : (x > hi ? hi : x).
 * When lo > hi the value comes from that same expression (x < lo yields lo, otherwise x > hi yields hi); nothing more is promised;
 * a NaN input lands by the same comparison rule as above (comparing false against both bounds, it is returned unchanged). */
fv_vec2 fv_vec2_clamp(fv_vec2 v, fv_vec2 lo, fv_vec2 hi);

/* a + b*s (per component). The implementation does not contract into an FMA: bit-identical to the scalar expression a.x + b.x*s. */
fv_vec2 fv_vec2_madd(fv_vec2 a, fv_vec2 b, fv_real s);

fv_real fv_vec2_dot(fv_vec2 a, fv_vec2 b);
fv_real fv_vec2_len_sqr(fv_vec2 v);
fv_real fv_vec2_len(fv_vec2 v);
fv_real fv_vec2_dist_sqr(fv_vec2 a, fv_vec2 b);
fv_real fv_vec2_dist(fv_vec2 a, fv_vec2 b);

/* A zero vector (length <= FV_EPS) returns (0, 0) and never amplifies into a NaN. */
fv_vec2 fv_vec2_normalize(fv_vec2 v);

/* The component of v along the onto direction.
 * Returns (0, 0) when onto's squared length <= FV_EPS_SQR (equivalently length <= FV_EPS).
 * P0 fix: this used to compare the squared length against the linear threshold FV_EPS, making the equivalent length
 * threshold 1e-3, 1000x looser than normalize; it now uses the squared threshold FV_EPS_SQR, aligned with normalize. */
fv_vec2 fv_vec2_project(fv_vec2 v, fv_vec2 onto);

/* Merge point p into the box: lo takes the per-component smaller value, hi the per-component larger value. */
fv_aabb2 fv_aabb2_merge(fv_aabb2 box, fv_vec2 p);

/* Overlap test: touching counts as overlapping (a.hi == b.lo and the like), returning 1; otherwise 0. */
int fv_aabb2_overlaps(fv_aabb2 a, fv_aabb2 b);

/* Per-component approximate equality, returning 1 / 0. */
int fv_vec2_near(fv_vec2 a, fv_vec2 b, fv_real eps);

/* ---------------- 3D ---------------- */

fv_vec3 fv_vec3_make(fv_real x, fv_real y, fv_real z);
fv_vec3 fv_vec3_add(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_sub(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_neg(fv_vec3 v);
fv_vec3 fv_vec3_scale(fv_vec3 v, fv_real s);
fv_vec3 fv_vec3_lerp(fv_vec3 a, fv_vec3 b, fv_real t);

/* The same semantics and NaN rule as fv_vec2_min / max / clamp / madd, applied per component. */
fv_vec3 fv_vec3_min(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_max(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_clamp(fv_vec3 v, fv_vec3 lo, fv_vec3 hi);
fv_vec3 fv_vec3_madd(fv_vec3 a, fv_vec3 b, fv_real s);

fv_real fv_vec3_dot(fv_vec3 a, fv_vec3 b);
fv_real fv_vec3_len_sqr(fv_vec3 v);
fv_real fv_vec3_len(fv_vec3 v);
fv_real fv_vec3_dist_sqr(fv_vec3 a, fv_vec3 b);
fv_real fv_vec3_dist(fv_vec3 a, fv_vec3 b);

fv_vec3 fv_vec3_normalize(fv_vec3 v);
/* Same threshold rule as fv_vec2_project (the squared length is compared against FV_EPS_SQR). */
fv_vec3 fv_vec3_project(fv_vec3 v, fv_vec3 onto);
/* Cross product a x b. */
fv_vec3 fv_vec3_cross(fv_vec3 a, fv_vec3 b);

int fv_vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps);

/* ---------------- 4D ----------------
 *
 * Layout: 4 x fv_real = 16 B, no padding bytes (see the typedef and the compile-time assertion above).
 *
 * The semantics line up with fv_vec2 / fv_vec3 one by one -- the same set of names, the same thresholds, the same NaN rule:
 *   min(a,b) = (a < b) ? a : b, max(a,b) = (a > b) ? a : b, **not** IEEE minNum/maxNum:
 *   with a NaN on one side the result is decided by the comparison (min(NaN, x) == x, min(x, NaN) == NaN, and when
 *   both are NaN the second one is returned); in clamp a component that compares false against both bounds is
 *   returned unchanged. That rule is exactly the one in fv_vec2_min/max/clamp, so the two compare directly.
 *   madd(a,b,s) = a + b*s, with no FMA contraction (bit-identical to the scalar expression).
 *   dot sums left to right: ((x*x' + y*y') + z*z') + w*w';
 *   len = sqrtf(dot(v, v)); normalize returns (0, 0, 0, 0) when len <= FV_EPS;
 *   project compares onto's squared length against the squared threshold FV_EPS_SQR (the same rule as fv_vec2_project);
 *   near is the per-component |a-b| <= eps.
 *   No cross: a 4D cross product is not a standard operation and nothing requires one (YAGNI).
 *
 * Implementation location: src/flashvec_vec4.cpp -- like fv_vec2/fv_vec3 these are in-library extern symbols,
 * not header inline. Two reasons:
 *   1) len / normalize need sqrtf; in the header they would force the FV_USE_EXTERN form to include <math.h> as well,
 *     which would break the existing rule "the FV_USE_EXTERN include set matches the pre-refactor one";
 *   2) the scalar API and the scalar tail of the batch entry points share one -ffp-contract=off translation unit,
 *     which makes bit-exactness easier to state and easier to check. */
fv_vec4 fv_vec4_make(fv_real x, fv_real y, fv_real z, fv_real w);
fv_vec4 fv_vec4_add(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_sub(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_neg(fv_vec4 v);
fv_vec4 fv_vec4_scale(fv_vec4 v, fv_real s);
fv_vec4 fv_vec4_lerp(fv_vec4 a, fv_vec4 b, fv_real t);
fv_vec4 fv_vec4_min(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_max(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_clamp(fv_vec4 v, fv_vec4 lo, fv_vec4 hi);
fv_vec4 fv_vec4_madd(fv_vec4 a, fv_vec4 b, fv_real s);

fv_real fv_vec4_dot(fv_vec4 a, fv_vec4 b);
fv_real fv_vec4_len_sqr(fv_vec4 v);
fv_real fv_vec4_len(fv_vec4 v);
fv_real fv_vec4_dist_sqr(fv_vec4 a, fv_vec4 b);
fv_real fv_vec4_dist(fv_vec4 a, fv_vec4 b);

fv_vec4 fv_vec4_normalize(fv_vec4 v);
fv_vec4 fv_vec4_project(fv_vec4 v, fv_vec4 onto);

int fv_vec4_near(fv_vec4 a, fv_vec4 b, fv_real eps);

/* ---------------- integrator primitives (per-value form) ----------------
 * addmul(x, v, s) = x + v*s and maddv(a, b, c) = a*b + c, per component -- the two expressions an explicit integrator puts
 * in its innermost loop (x += v*dt). The evaluation order is fixed at "multiply first, then add" and is never contracted
 * into an FMA, so the result is bit-identical to the scalar expressions x.x + v.x*s and a.x*b.x + c.x (the same promise
 * fv_vec2_madd carries). The batch forms of the same two expressions are fv_vec2_addmul_x4 / fv_vec2_addmul_soa /
 * fv_vec3_addmul_soa / fv_vec4_addmul_x2 / fv_vec2_maddv_x4 below.
 * Definition: src/flashvec_batch.cpp, under the same fp-contract=off pragma as the rest of the scalar API -- a consumer
 * built with -mfma -ffp-contract=fast must still get the unfused bit pattern from these symbols. */
fv_vec2 fv_vec2_addmul(fv_vec2 x, fv_vec2 v, fv_real s);
fv_vec3 fv_vec3_addmul(fv_vec3 x, fv_vec3 v, fv_real s);
fv_vec4 fv_vec4_addmul(fv_vec4 x, fv_vec4 v, fv_real s);
fv_vec2 fv_vec2_maddv(fv_vec2 a, fv_vec2 b, fv_vec2 c);

/* ---------------- quaternion / 3x3 matrix scalar API (physics primitives, P6) ----------------
 *
 * This group is the first half of "the batch layer already exists; what was missing was the domain
 * types". It adds rotation and 3x3 linear algebra as first-class values with the same contract as the
 * rest of the scalar API: by value, no allocation, no global state, no FMA contraction, and every
 * expression's association order is fixed so the batch form below can be compared with memcmp.
 *
 * Evaluation order that the batch entries reproduce word for word:
 *   fv_quat_mul   x = ((w1*x2 + x1*w2) + y1*z2) - z1*y2, and the three analogous lines;
 *   fv_quat_rotate_vec3   t = 2 * cross(q.xyz, v); d = cross(q.xyz, t); out = (v + q.w*t) + d;
 *   fv_mat3_transform_vec3   out.x = (m00*x + m01*y) + m02*z, and the two analogous lines;
 *   fv_mat3_mul   (a*b)[i][j] = (ai0*b0j + ai1*b1j) + ai2*b2j.
 * These are ordinary left-to-right sums; they are written out in the header's contract comment so a
 * future SIMD kernel has a single place to match, rather than an implicit compiler-dependent order.
 *
 * Semantics that are easy to get wrong and are therefore stated here:
 *   - fv_quat_mul(a, b) applies b first, then a (matrix composition order);
 *   - fv_quat_normalize maps a near-zero quaternion (length <= FV_EPS) to the identity (0,0,0,1),
 *     not to (0,0,0,0): a zero quaternion is not a rotation, and returning the identity keeps the
 *     result usable and NaN-free. This is the one deliberate difference from fv_vec*_normalize, which
 *     returns the zero vector; the reason is the rotation contract, not a change of style;
 *   - fv_quat_rotate_vec3 does **not** normalize q. Like a rotation matrix, q is assumed unit;
 *   - fv_mat3_from_quat assumes a unit quaternion and produces the row-major rotation matrix R with
 *     v' = R*v. It is the matrix form of the same rotation, with its own (also fixed) expression order;
 *   - fv_mat3_mul(a, b) applies b first then a, matching fv_quat_mul and fv_mat4's row-major use.
 *
 * Interpolation and angle construction (P7): fv_quat_nlerp / fv_quat_slerp use the shortest-arc
 * convention (dot < 0 negates the second quaternion) and normalize the result; fv_quat_from_axis_angle
 * and fv_quat_from_euler_xyz are the sin/cos constructors. These are the points where the exactness
 * story changes: they involve trigonometry, so their results are accuracy-bounded (not bit-reproducible
 * across libm implementations) and they have no batch form. The exact, bit-reproducible properties
 * remain the arithmetic/transform entries above.
 */
fv_quat fv_quat_make(fv_real x, fv_real y, fv_real z, fv_real w);
fv_quat fv_quat_add(fv_quat a, fv_quat b);
fv_quat fv_quat_sub(fv_quat a, fv_quat b);
fv_quat fv_quat_neg(fv_quat q);
fv_quat fv_quat_scale(fv_quat q, fv_real s);

fv_real fv_quat_dot(fv_quat a, fv_quat b);
fv_real fv_quat_len_sqr(fv_quat q);
fv_real fv_quat_len(fv_quat q);
/* Length <= FV_EPS returns the identity (0,0,0,1), not (0,0,0,0); see the block comment above. */
fv_quat fv_quat_normalize(fv_quat q);

/* Conjugate (-x, -y, -z, w): the inverse of a unit quaternion is its conjugate. */
fv_quat fv_quat_conjugate(fv_quat q);

/* Inverse of a general quaternion: conjugate(q) / |q|^2. This is the inverse for non-unit q as well;
 * for a unit q it equals the conjugate. Length <= FV_EPS (|q|^2 <= FV_EPS_SQR) returns the **zero**
 * quaternion (0,0,0,0) -- the inverse does not exist, and zero is the signal. It is deliberately not
 * the identity: an identity would hide the degenerate input. The _checked form returns 1 on success,
 * 0 on the degenerate case, and writes the same zero value to *out (out must be non-NULL). */
fv_quat fv_quat_inverse(fv_quat q);
int fv_quat_inverse_checked(fv_quat q, fv_quat *out);

/* Hamilton product a*b: applies b first, then a. */
fv_quat fv_quat_mul(fv_quat a, fv_quat b);
/* Rotate v by the **unit** quaternion q: v' = q * (v, 0) * conjugate(q). No normalization is done here. */
fv_vec3 fv_quat_rotate_vec3(fv_quat q, fv_vec3 v);

/* Shortest-arc normalized linear interpolation: if dot(a,b) < 0, b is negated first, then
 * out = normalize(a + (b - a)*t). t is not clamped (extrapolation is allowed), and the result is
 * always normalized. The sign fix means exactly opposite quaternions become equal, so the usual
 * "opposite at t = 0.5" degeneration does not occur; the only degenerate case is a zero input or a
 * zero interpolated vector, where normalize returns the identity. */
fv_quat fv_quat_nlerp(fv_quat a, fv_quat b, fv_real t);
/* Shortest-arc spherical interpolation. It falls back to nlerp when the angle is tiny (dot above
 * 1 - 1e-5) or when sin(theta) underflows, so the endpoints and the near-identical case stay stable.
 * Inputs are assumed unit (the result is normalized); t is not clamped. */
fv_quat fv_quat_slerp(fv_quat a, fv_quat b, fv_real t);

/* Construction from axis-angle: the axis is normalized internally (zero axis -> the identity), the
 * angle is in radians, and the result is the unit quaternion (axis*sin(angle/2), cos(angle/2)). */
fv_quat fv_quat_from_axis_angle(fv_vec3 axis, fv_real angle);
/* Construction from Euler angles in radians, applied X first, then Y, then Z:
 * q = qz * qy * qx (so qx acts on the vector first). The result is unit. */
fv_quat fv_quat_from_euler_xyz(fv_real x, fv_real y, fv_real z);

/* Row-major 3x3 from three row vectors (r0 = m[0..2], r1 = m[3..5], r2 = m[6..8]). */
fv_mat3 fv_mat3_from_rows(fv_vec3 r0, fv_vec3 r1, fv_vec3 r2);
/* out = m * v, with m row-major and v a column vector. */
fv_vec3 fv_mat3_transform_vec3(const fv_mat3 *m, fv_vec3 v);
/* out = a * b (row-major, applies b first then a). */
fv_mat3 fv_mat3_mul(const fv_mat3 *a, const fv_mat3 *b);
/* Transpose: out.m[i*3+j] = m->m[j*3+i]. */
fv_mat3 fv_mat3_transpose(const fv_mat3 *m);
/* Row-major rotation matrix of a **unit** quaternion: v' = R*v has the same rotation as fv_quat_rotate_vec3. */
fv_mat3 fv_mat3_from_quat(fv_quat q);
/* Inverse mapping: rotation matrix -> unit quaternion (Shepperd's branch-on-largest-diagonal method,
 * so the divisor is never the small one). The matrix must be a rotation (orthonormal, det +1); the
 * result is normalized. A non-rotation still yields a unit quaternion, but it is not that matrix's
 * rotation -- there is no error channel in the by-value form. */
fv_quat fv_quat_from_mat3(const fv_mat3 *m);

/* Determinant, row-major, evaluated as
 * m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20). */
fv_real fv_mat3_determinant(const fv_mat3 *m);
/* Inverse = adjugate / determinant. A matrix is treated as singular when
 * |det| <= FV_EPS * max|entry|^3 (a scale-relative test: uniformly scaling an invertible matrix by
 * a positive factor keeps it invertible). The value form returns the **zero matrix** for a singular
 * input; the _checked form returns 1 on success, 0 on a singular input, and writes the zero matrix
 * either way. In the _checked form *out may alias *m (all reads happen before the write). */
fv_mat3 fv_mat3_inverse(const fv_mat3 *m);
int fv_mat3_inverse_checked(const fv_mat3 *m, fv_mat3 *out);

/* ---------------- 4x4 matrix operations (P7) ----------------
 * fv_mat4 is row-major (m[row*4+col], the same convention as fv_mat3 and the opposite of GLM). The
 * transform entry points above keep the affine precondition (last row (0,0,0,1)); the operations here
 * are general 4x4 operations and do not require it.
 *
 *   mul: out = a*b, row-major, applies b first then a; the inner sum is left to right
 *        ((ai0*b0j + ai1*b1j) + ai2*b2j) + ai3*b3j.
 *   inverse: the adjugate / determinant closed form. A matrix is treated as singular when
 *        |det| <= FV_EPS * max|entry|^4; the value form returns the zero matrix and the _checked form
 *        returns 1/0 (writing the zero matrix on failure). *out may alias *m in the _checked form.
 *
 * General 4x4, not affine-only: an affine-specialized inverse would be cheaper, but "full operations"
 * means a caller can multiply or invert the matrices the library already accepts, including a
 * non-affine one.
 */
fv_mat4 fv_mat4_mul(const fv_mat4 *a, const fv_mat4 *b);
fv_mat4 fv_mat4_inverse(const fv_mat4 *m);
int fv_mat4_inverse_checked(const fv_mat4 *m, fv_mat4 *out);

/* ---------------- approximate math (P7) ----------------
 * fv_atan_approx is a polynomial for atan on the documented domain [-1, 1], evaluated with Horner's
 * rule. Max absolute error is about 1.2e-5 (tested against atanf). Outside [-1, 1] and for NaN it
 * returns NaN: clamping would silently return a wrong angle, and the caller can use the identity
 * atan(x) = sign(x)*pi/2 - atan(1/x) for the rest of the range. This is an approximate entry point:
 * it does not promise the libm result, and there is no batch form.
 */
fv_real fv_atan_approx(fv_real x);

#endif /* FV_USE_EXTERN || FV_BUILD: scalar API declarations */

/* ============================ batch entry points (P1) ============================
 *
 * Why this group exists: the per-element API is a call across translation units, which GCC outright refuses to vectorize
 * (missed: statement clobbers memory); an **early probe** measured the boundary cost at roughly 4.2 cycles/element,
 * 3x the packing gain (early microbenchmark; this round's hardened fixtures re-measured the batch entry points themselves, not this boundary cost, so read it only as an order of magnitude).
 * A batch entry point trades "one call per element" for "one call per batch"; the library then decides whether to use SIMD.
 *
 * Implementation forms (three macros are the whole switchboard)
 *   default         this header provides static inline implementations of the P1 group's 52 batch entry points, the P6 group's four physics entries and the P7 mat3/fast entries (see "implementation of the batch entry points" at the end),
 *                   so the caller's translation unit gets the function bodies directly;
 *   FV_USE_EXTERN   the header carries only the extern declarations below and calls go into libflashvec.a -- the pre-refactor behaviour;
 *                   this is also the entry for a C caller that wants the AVX2 kernels (a C translation unit takes no intrinsics; see the file header);
 *   FV_BUILD        used when building the library itself (src/flashvec_batch.cpp defines it); externally it behaves like FV_USE_EXTERN:
 *                   declarations only, with the definitions in the .cpp where each body is one line forwarding to fv_i_*.
 *   Both forms are one implementation: the .cpp copies no kernel, it only wires fv_i_* to the public names.
 *
 * Layout
 *   _x4   : input and output are AoS arrays of fv_vec2. fv_vec2 is 8 B, so 4 of them are exactly 32 B,
 *           one 256-bit register holds 4 elements, **zero repacking**, and the caller's existing data is passed as is.
 *   _soa  : input and output are separate component arrays (x[], y[]). The reduction-style operations (dot/len/normalize) are most natural here.
 *   _mask : whole-batch read-modify-write; a lane outside the mask outputs a[i]'s original value (see the mask semantics below).
 *   _idx  : processes only the indices named in idx[0..m) and really skips the rest (at the cost of random access).
 *   _x2   : fv_vec4 only (16 B/element). 2 elements = 32 B = one ymm, **exactly full, not one lane
 *           wasted**. The name is isomorphic to vec2's _x4: the number says how many elements fit in one 256-bit
 *           register (vec2 is 8 B/element -> 4; vec4 is 16 B/element -> 2), not how many vectors are handled at once.
 *           The packing problem does not arise here: fv_vec4 has no padding, so two adjacent array elements are 32 B of
 *           contiguous memory, hence each input is one 256-bit loadu and the output is one 256-bit storeu -- as with _x4,
 *           zero repacking. A layout that would need two 128-bit loads plus _mm256_insertf128_ps to assemble
 *           (an element stride other than 16 B, such as AoSoA or an array with a stride) is not the AoS used here.
 *   3D has no AoS batch, only _soa:
 *           fv_vec3 stays 12 B, **no padding added** (padding it to 16 B would break the caller's existing 12 B layout;
 *           the documents carry a clear conclusion). 12 B does not divide 32 B evenly: an AoS batch would either idle
 *           1/4 of the lanes in every ymm (3/4 utilization) or unpack with shuffles and pack again, paying the saved width back.
 *           So the only form of a 3D batch is SoA (fv_vec3_*_soa, signatures of 4-9 pointers + n) --
 *           one ymm carries **8 different elements of the same component**, at 100% lane utilization.
 *           fv_vec3_cross has no batch entry point: the 3 yz swaps would need a lot of shuffling for an unclear gain (YAGNI).
 *
 * Where the bit-exactness promise ends (the one thing callers most need to know)
 *   Guaranteed - bit-identical to the corresponding scalar API (same input, same float bit patterns under memcmp):
 *     fv_vec2_add_x4 / sub_x4 / scale_x4 / lerp_x4 / madd_x4;
 *     fv_vec2_add_x4_mask (including the half that writes a[i] back outside the mask);
 *     fv_vec2_add_x4_idx / scale_x4_idx;
 *     fv_vec2_add_soa / scale_soa / len_soa / normalize_soa;
 *     fv_vec2_madd_normalize_soa / madd_normalize_soa_stream (fused pipeline; identical to scale_soa -> add_soa -> normalize_soa);
 *     fv_vec3_add_soa / scale_soa / len_soa / normalize_soa (3D has the SoA form only);
 *     fv_vec4_add_x2 / sub_x2 / scale_x2 / lerp_x2 / madd_x2;
 *     fv_vec4_len_soa / normalize_soa;
 *     fv_mat4_transform_points / transform_dirs (the matrix reads only its first 3 rows; see the AoS part at the end of this section);
 *     fv_vec2_addmul_x4 / addmul_soa / fv_vec3_addmul_soa / fv_vec4_addmul_x2 and fv_vec2_maddv_x4 / maddv_x4_stream
 *     (in place and three-operand madd: one product, then one sum, inside the pragma, so no contraction is possible;
 *     the per-value fv_vec2_addmul / fv_vec3_addmul / fv_vec4_addmul / fv_vec2_maddv are those same expressions).
 *   Not guaranteed - fv_vec2_dot_soa: it uses 4 independent accumulators plus one final reduction to break the dependency chain
 *     (early microbenchmark x2.25, not re-measured this round), and its reduction order differs from the per-element path, so **the result depends on the reduction order**:
 *     same input gives the same output (determinism), but it is not bit-identical.
 *     The error has no closed-form bound; everything below is measured (a different corpus can differ by two orders of magnitude, so quote it with its conditions):
 *       - n <= 1000 random data, compared after normalizing by scale: relative error <= 4e-8 (tests/test_batch.cpp)
 *       - n = 8192 single dot product, against the per-element path: absolute error <= 4.8e-7, relative error <= 6.2e-8
 *       - near cancellation (a tiny denominator for the relative error) it can reach ~3e-6 relative error
 *     Use the scalar fv_vec2_dot when bit-exact reproduction is required.
 *     For n < 8 it degrades to plain sequential accumulation, where it is bit-identical to the per-element path (the +0.0 vs -0.0 boundary is below).
 *   The same exception applies to fv_vec3_dot_soa / fv_vec4_dot_soa: they use the same reduction
 *     scheme as fv_vec2_dot_soa (several independent accumulators plus one final tree reduction), so same input gives the same output but not a bit-identical one;
 *     the tolerance wording is the same too: after normalizing by scale (sum of absolute terms, accumulated in double) the relative error is < 1e-5
 *     (tests/test_batch.cpp measures <= 4e-8; this round's port to the tests/test_parity.c corpus is reported separately).
 *     Use the scalar fv_vec3_dot / fv_vec4_dot when bit-exact reproduction is required.
 *     Implementation detail: when n is under one 8-lane step the whole span takes the tail path, but that tail is "sum accumulated term by term starting from +0.0",
 *     whereas the scalar API does not add its first term to +0.0; so when every term is -0.0 the two differ by one sign bit
 *     (+0.0 vs -0.0). This matches vec2's dot_soa and falls under the exception stated above.
 *   Split invariance (another property, easily missed):
 *     The per-element batch entry points (add/sub/scale/lerp/madd/len/normalize and their _x4/_soa/_mask/_idx variants)
 *     do have it -- split one array into any number of pieces and call them separately, the output is bit-identical to a single call
 *     (measured: add_x4 12/12 bit-identical).
 *     The reduction entry points (dot_soa, including the vec3/vec4 versions) **do not** have it: the reduction tree changes with the split.
 *     Measured n=8192 "one span vs two pieces" relative difference 2.4e-5, n=1000 about 3.9e-6, 12/12 all different.
 *     When the result must not depend on the split, use sequential accumulation in the scalar fv_vec2_dot.
 *   Where the guarantee comes from: these paths use only vaddps / vsubps / vmulps / vdivps / vsqrtps,
 *   no FMA (vfmadd), no reciprocal approximations (rcpps/rsqrtps), and the reduction order (except dot_soa) matches the scalar one.
 *   The implementation section is split in two by "is there a mul+add contraction risk" (see the note in "implementation of the batch entry points" below):
 *   the first section (add/sub/scale and their mask/idx/soa siblings) has no contractible expression and needs no pragma,
 *   so under the consumer's default options it inlines completely (including the n < 4 / n < 8 tail loops);
 *   the second section (lerp/madd/len/normalize/dot) is covered as a whole by #pragma GCC optimize("fp-contract=off")


 *   (push and pop paired, covering only those functions in this header), at the cost of one intra-TU call per batch.
 *   Other compilers and the caller's own compiler options are outside this pragma's reach (see "responsibility boundary" in the file header).
 *   How it is verified: in objdump -d build/libflashvec.a the counts of vfmadd / vrcpps / vrsqrtps / vhaddps are 0;
 *   tests/test_batch.cpp compares every function at every tail length down to memcmp;
 *   make test-inline compares the output of the default form against the FV_USE_EXTERN form byte for byte.
 *
 * Tail length
 *   n may be any value, including 0. A tail shorter than one vector width (_x4: 0-3 elements, _x2: 0-1,
 *   _soa: 0-7)
 *   takes the scalar path, with semantics identical to the main loop. For n == 0 no pointer is read or written (NULL may be passed).
 *
 * Mask semantics (easy to misread, so they must be spelled out)
 *   mask[i] != 0 gives o[i] = the computed result; mask[i] == 0 leaves o[i] **explicitly written back with a[i]'s original value**.
 *   Not "left unwritten" -- a 256-bit register can only be written back whole, skipping the store is impossible, and skipping would leave
 *   half a line in o that disagrees with a. For "keep o's old value" semantics use the _idx variant,
 *   or point a at that old value before the call.
 *
 * Index semantics
 *   For k in [0, m): i = idx[k], then o[i] = the computed result (reading a[i], b[i]).
 *   idx may repeat, may be out of order, need not be monotonic; a repeated index writes the same value, so the result does not depend on the processing order.
 *   Indices outside idx are neither read nor written; the rest of o keeps its value from before the call.
 *
 * Pointers and aliasing
 *   a / b / o carry __restrict and by contract must not overlap (o sharing an address with a or b is UB);
 *   idx / mask must not overlap them either.
 *
 * Other
 *   - n is measured in elements, not bytes.
 *   - The zero-vector guard in normalize_soa uses a mask instead of a branch, with the same semantics as fv_vec2_normalize
 *     (len <= FV_EPS outputs (0,0)); the implementation replaces the length with 1 before dividing, so a zero-vector lane
 *     never triggers a division by zero (nor sets an MXCSR flag).
 *
 * ---- _stream (NT store) variants: use them only when the working set is clearly larger than L3 ----
 * What it does: the main loop writes results with _mm256_stream_ps (a non-temporal store), bypassing the cache, so it **produces no
 * RFO read of "read the output line in first"**. That is its entire physical meaning -- the ceiling is decided purely by byte accounting,
 * and byte accounting must be **computed per entry shape, never merged**: how many floats one element moves is stated only by that entry's signature.
 *
 *   entry                        B/e read  B/e written  nominal  with RFO  ceiling (real ratio)
 *   add_x4 / sub_x4 / add_soa     16      8      24     32    32/24 = 1.333
 *   scale_x4 / scale_soa           8      8      16     24    24/16 = 1.500
 *   normalize_soa                  8      8      16     24    24/16 = 1.500
 *   dot_soa / len_soa (pure read)   16 / 8    0       -      -    no gain (NT affects writes only)
 *   lerp_x4 / madd_x4             16      8      24     32    32/24 = 1.333
 *   vec4 _x2 add/sub/lerp/madd    32     16      48     64    64/48 = 1.333
 *   vec4 _x2 scale / vec4 norm    16     16      32     48    48/32 = 1.500
 *   vec3 SoA add                  24     12      36     48    48/36 = 1.333
 *   vec3 SoA scale / normalize    12     12      24     36    36/24 = 1.500
 *   mat4 transform points / dirs  16     16      32     48    48/32 = 1.500
 *
 * The "read" column follows the signature: _x4 reads two fv_vec2 arrays, a / b (16 B); add_soa reads ax / ay / bx / by,
 * **four** float arrays, which is 16 B as well (not 8 B) -- so it has the same shape as _x4 and the same 1.333 ceiling.
 * The "written" column works the same way: add_soa writes two arrays, ox / oy (8 B); the scale family reads only a (8 B) and writes only o (8 B).
 * An earlier version described _soa as "2 reads, 2 writes, 24/16 = 1.5x", missing that add_soa reads 4 arrays; that was wrong;
 * the ceiling can only be computed per entry, never carried over as one number for all variants.
 *
 * These are the **expected ceiling** given by byte accounting, not a physical guarantee. Measured PASS points land in 1.29-1.58x, while the 4M tier's paired
 * IQR is 66-98% and the 1M/8M tiers are 20-35%, so **any single point within 20% above the ceiling counts as noise** (for example add_soa
 * measured 1.582 against a ceiling of 1.333, 19% over, inside the noise band); never claim from that "the measurement beat the theoretical ceiling".
 * When reporting GB/s give both figures (nominal / with RFO).
 * When **not** to use it (this matters more than the gain): when the result is read right away, or the working set sits in L2/L3 --
 * the penalty there is **2-8x**, not ten or twenty percent.
 * Measured (Ryzen 7 5700U, L3 = 4 MiB/CCX, cores pinned with taskset; same-process interleaved A/B; per point 4 processes x 2 builds
 * x 2 orders, including the same-round A/A floor):
 *   24 MB working set (n=1M): add 1.350 / sub 1.293 / scale 1.313
 *   192 MB working set (n=8M): add 1.422 / sub 1.519 / add_soa 1.582
 *   inside L2 (n=8192): add 0.217 / scale 0.130 / add_soa 0.229 / normalize 0.506 (2-8x slower)
 * Note: the 1M/8M tiers still have a paired IQR of 20-35%, so the ratios above read only as "of order 1.3-1.6x", and a single point only as
 * "cleared or did not clear the same-round A/A floor". This round's verdict (docs/research/raw/ab_campaign_verdict.txt): n=8192
 * all six _stream points FAIL (slower); n=4M all six points INCONCLUSIVE (floor > 15%) -- that tier's same-round A/A floor is
 * 1.34-1.61 with a paired IQR of 66-98%, unmeasurable and not usable as evidence; n=131072 is the blacklisted tier (see the acceptance protocol);
 * PASS occurs only in the 1M and 8M tiers (add 1.350 / sub 1.293 / scale 1.313; 1.422 / 1.519 / add_soa 1.582).
 * Do not claim a ceiling was beaten from a single run, and do not back out a ceiling from a 4M-tier single point.
 * Conclusion: **the caller must choose explicitly; the library never switches automatically**. The criterion is the working-set byte count (elements x 24 B)
 * being clearly larger than L3; with this machine's 4 MB L3, n <= 131072 (3 MB) should use the normal version.
 * Alignment: the caller **need not** guarantee any alignment. The kernel pushes the output to a 32 B boundary itself (at most 7 leading elements go through
 * normal stores), the body uses NT, and the remaining tail keeps using normal stores; fv_vec2 has an 8 B stride, so if the output
 * array is not even 8 B aligned (the 32 B boundary is unreachable), the kernel falls back to normal storeu throughout -- the result is bit-identical,
 * only the NT gain is gone. SoA's two outputs are handled by two independent kernels (each peels and decides for itself whether NT is possible),
 * so they do not drag each other down when their phases differ; normalize, whose two outputs share one computation, requires
 * ox / oy to reach the 32 B boundary **together**, and falls back to normal storeu throughout if they cannot (allocating both outputs at 32 B is recommended).
 * Ordering: every stream kernel executes _mm_sfence() before returning. An NT store is a **weakly ordered** write; the data sits in a write-combining buffer first;
 * returning without an sfence can leave other cores (and a later DMA / device read) unable to see those writes.
 * This is not a theoretical worry: with the sfence removed, a hand-off experiment ("writer thread does NT stores -> release flag -> reader thread acquires
 * and verifies element by element") run 12 times **read stale values 6 times**; with the sfence all 12 were correct.
 * Note also: a single read detects it only about 50% of the time, so no test that "reads once and checks" can claim to catch this
 * bug reliably -- the proper measurement is many repetitions plus a two-thread hand-off (which is how tests/test_stream.c is written).
 * Result: **bit-identical** to the non-stream version (the same set of vaddps / vsubps / vmulps / vdivps / vsqrtps,
 * only the store instruction changed). tests/test_parity.c does a memcmp for every entry, over all tail lengths n = 0..17, a random large n, and every
 * phase shift of the output pointer, and feeds the output into the cross-form comparison stream.
 * Sectioning: every _stream entry followed the same rule -- a kernel whose arithmetic contains mul+add (lerp / madd / normalize /
 * mat4 transform, i.e. a*b feeding straight into a +) sits in the **second section**, inside the pragma, at the cost of one intra-TU
 * call per batch in a consumer; the pure add/sub/scale ones sit in the first section and inline completely under the consumer's
 * default options. Today that is 10 first-section entries (add/sub/scale _x4, add/scale SoA, vec3 add/scale, vec4 add/sub/scale _x2)
 * and 15 second-section ones (lerp/madd _x4, normalize SoA, vec4 lerp/madd _x2, vec3/vec4 normalize SoA, both mat4
 * transforms, and the six integrator entries: addmul _x4/_soa/v3_soa/_x2 plus maddv _x4 and its _stream form).
 *
 * ---- _mask_inplace: in-place sparse add (a semantic entry point; no measured throughput advantage) ----
 * Why a **new function** instead of changing _mask:
 *   the _mask contract is "whole-batch read-modify-write" -- lanes outside the mask are **explicitly written back with a[i]'s original value** (see "mask semantics"),
 *   and it requires the output pointer o not to overlap a/b (o at the same address as a is UB). An in-place caller is exactly the case that cannot do this:
 *   it has no second buffer. The same arguments in _mask cannot mean both "a new output" and "the old array rewritten in place",
 *   so the two are split into two entry points semantically, rather than relaxing the _mask contract.
 *
 * What is promised is **value semantics**, and only three things:
 *   1) elements with mask[i] != 0 are bit-identical to the scalar a[i] + b[i] (and bit-identical to the hit lanes of the non-in-place
 *      _mask -- both take the same vaddps path);
 *   2) elements with mask[i] == 0 keep their **value** (verifiable: compare their bytes against the value before the call, which is what tests/test_stream.c
 *      does);
 *   3) everything else (whether a non-hit element is physically written back, how many times, what state its cache line is in) **carries no promise at all**.
 * In particular: **this function is not a synchronization primitive** and gives no ordering guarantee between concurrent readers and writers; like fv_vec2_add,
 * concurrent access to the same memory is still a data race (UB). A need for "physically not writing" is not a vector library's job.
 * Why "write only the hit lanes" is not done: it needs VMASKMOVPS (a vmaskmovps maskstore). An early probe on this machine
 * (Zen 2) found this family to be microcoded and slow, 1.2-3.8x slower than a whole-block storeu at the same density, and slower than
 * "block prescan + scalar" at low density (about 2.5x at 1% hits and a 64 MB working set) -- these numbers come from probes taken before the hardened
 * fixtures and were **not re-measured** this round, so read them as an order of magnitude. The reason that does not depend on the numbers at all: keeping it would turn
 * this function's "physically writes only the hit lanes" into a property that **varies with mask density** -- no caller can depend on a storage
 * guarantee that varies with the input distribution while paying a performance cost for it. So neither write strategy is used here.
 *
 * Aliasing: a is read-write, so _mask's o/a aliasing problem does not exist; the only constraint is that b must not overlap a
 *   (nor may mask overlap either of them).
 * The implementation is **two-stage density adaptive**; both stages produce bit-identical results and the choice affects performance only:
 *   sampling: inspect only the first FV_I_MIP_SAMPLE (64) mask bytes and count the non-zero ones (a few dozen cycles, negligible against a whole batch),
 *         getting a density estimate of the right order of magnitude; a wrong estimate costs performance, not semantics.
 *   <= NUM/DEN hits (1/8 by default): block prescan + per-element scalar. A block whose four mask bytes are all zero is skipped whole
 *         (not even a's cache line is touched; at 5% hits these blocks are 95%^4 = 81% of all blocks).
 *   >  NUM/DEN hits: whole-block 32 B load a / load b + blend + storeu back into a. No data-dependent branches,
 *         and the non-hit lanes write a's original value back (the value is unchanged, but the write happens).
 *   The cut-off is controlled by three macros, tunable without touching the function body:
 *     FV_I_MIP_DENSE_NUM / FV_I_MIP_DENSE_DEN  -- density threshold (1/8 by default)
 *     FV_I_MIP_SAMPLE                          -- number of sampled bytes (64 by default)
 *   This routing only decides "which loop runs" and **brings no measurable throughput difference** (see the capability boundary below).
 *
 * Capability boundary (measured; read this): this entry point is filed as "tried and ineffective", **not** as "an optimized entry point".
 *   Measured (hardened fixtures, 5 densities x 4 size tiers, including the same-round A/A floor): against _mask (non-in-place) this entry
 *   has **no measurable throughput advantage**; the gain lies between -0.10 and +0.00, all of it below the A/A floor
 *   (A/A floor 1.06-1.56; 1%/5%/20%/50%/95% x 4 size tiers, every verdict INCONCLUSIVE).
 *   Reason: writing in place produces no RFO -- a has to be read before it is written anyway -- so "in place" saves no traffic.
 *   Its value is the **form** (one output buffer less, one allocation less for the caller), not speed. For throughput use
 *   the non-in-place version and let the compiler see a wider loop.
 * Aside: the block prescan test in _mask is "four mask bytes moved into one uint32 give zero for the whole value",
 *   which **must not** be written as (w & 0x01010101u) == 0 -- that tests "the lowest bit of every byte is 0", so a mask byte of
 *   0x02 / 0x80 / 0xFE, all legal hits, would be misjudged as a miss (tests/test_stream.c has a regression case).
 * ---- fv_mat4_transform_points / fv_mat4_transform_dirs: matrix point transform (AoS) ----
 * Semantics: out[i] = M * in[i], where in / out are the caller's existing AoS arrays (4 floats per point).
 *   points treats in[i].w as 1 (translation column included), dirs treats in[i].w as 0 (translation ignored).
 *   The matrix is passed **by pointer**: 64 B by value would copy it on every call.
 *   out.w is written as a constant (point -> 1, direction -> 0) and takes no part in the computation; the precondition and its reason are in the fv_mat4 typedef above.
 *
 * Why **only** AoS and no SoA version (a measured conclusion, not a lack of time):
 *   bench/tfprobe/probe.cpp, Ryzen 7 5700U / Zen 2 / AVX2, reproduced over three or more runs, with bit-exactness compared at every tier:
 *     layout               N=8192   N=32768   N>=131072   A/A floor  bit-exact
 *     AoS (fv_vec4[])      2.25x    2.03x     ~1.0        0.4%-10%   yes
 *     SoA                  1.01x    0.92x     ~1.0        0.4%       n/a
 *   The reason is the **instruction count** (confirmed with objdump): the compiler does vectorize the naive AoS loop, only badly --
 *   440 instructions and 14 vpermps per iteration, against 107 instructions and 3 vpermps in the hand-written kernel. For SoA it is the other way round:
 *   the compiler already needs 256, a hand-written version would need 387, and hand-written would only be slower. So "add a SoA version for symmetry"
 *   would mean shipping an entry point that **is measurably slower**, and it is not added.
 *   Read 2.25x / 2.03x together with the A/A floor: it means "clearly above the same-round noise floor", not a precise multiplier.
 *
 * Bit-exactness (the same promise as the other non-dot entry points): for any n and any tail length the output is
 *   **bit-identical** to the scalar reference implementation (tests/test_mat4.cpp checks with memcmp; the corpus has -0.0 / +/-inf / NaN / extremes).
 *   Only two conditions hold it up, and both must be kept:
 *     1) the evaluation order is fixed at ((m0*x + m1*y) + m2*z) + m3 -- the scalar tail and the AVX2 kernel use the same order;
 *     2) no FMA appears -- so the kernel must sit in the fp-contract=off section above: with a consumer using -mfma
 *        target("avx2") cannot stop vfmadd (both make test-fma and the disassembly verify this).
 *   dirs **drops** the last term instead of adding m3 * 0: x + (+0.0) would flip -0.0 into +0.0,
 *   and the kernel has no such add either, so the two match.
 *   The one exception (the same one the existing lerp/madd entry points carry, with the same note in tests/test_vec4.c): when both
 *   **operands of one arithmetic instruction are NaN with different bit patterns**, which NaN comes out is decided by that instruction's operand order.
 *   The compiler allocates the operands of a commutative operation freely, and it does: measured, GCC emits `vmulps` with the operands reversed in two
 *   kernels of the same binary, and fv_vec2_normalize_soa returns different NaN payloads for one and the same input when the caller compiles it with
 *   -march=native than when it links the library build (a sweep over 25 NaN-pattern pairs x 3 sizes differed on 164320 component comparisons in the
 *   -march=native tier and on 0 in the -O2 and default tiers). So under such input the scalar tail, the AVX2 kernel, and the same kernel compiled
 *   another way may each produce a different NaN bit pattern -- all of them legal NaNs, not a numerical error.
 *   As long as at most one NaN bit pattern occurs in the corpus, bit-exactness still holds: the hardware-generated indefinite is always the same bit
 *   pattern, and a one-sided NaN propagates from that one operand. The tests/test_mat4.cpp and tests/test_madd_norm.cpp corpora are grouped by this
 *   rule, and the exception is counted and printed separately.
 *
 * Aliasing: m / in / out carry __restrict and by contract must not overlap (the same rule as the other batch entry points).
 *   n == 0 reads and writes no pointer (NULL may be passed): the coefficient vectors are loaded only when bulk >= 2.
 * Layout: 2 fv_vec4 = 32 B = one ymm, loadu / storeu with zero repacking; 3 vpermps broadcast x / y / z
 *   into their 4 lanes each. So one iteration has only 2 memory accesses plus a dozen-odd arithmetic instructions, and memory is the remaining bottleneck:
 *   at n = 131072, 131072 * 32 B = 4 MiB = this machine's L3 per CCX (measured: that tier drops back to ~1.0).
 *
 * ====================================================================== */

/* Two forms of the public names: extern declarations (FV_USE_EXTERN / FV_BUILD) or static inline definitions (default).
 * The two cannot appear together (an extern declaration followed by a static definition in one TU is a hard error), hence the conditional compilation here. */
#if defined(FV_USE_EXTERN) || defined(FV_BUILD)

void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n);
void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n);
void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                      fv_vec2 *__restrict o, size_t n);
void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real t, fv_vec2 *__restrict o, size_t n);
void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real s, fv_vec2 *__restrict o, size_t n);

/* Sparse variant: add only when mask[i] != 0, otherwise o[i] = a[i] (a's original value written back). */
void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                         fv_vec2 *__restrict o, const uint8_t *__restrict mask, size_t n);
/* Sparse variant: for k in [0,m) take i = idx[k] and compute o[i] = a[i] + b[i]; idx may repeat and may be out of order. */
void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                        fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m);
void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                          fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m);

/* SoA component-array entry: 6 pointers + n, each __restrict. */
void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict bx, const float *__restrict by,
                     float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay, fv_real s,
                       float *__restrict ox, float *__restrict oy, size_t n);
/* The only entry point not promised to be bit-identical: 4 accumulators + one final reduction, so the result depends on the reduction order (see above). */
float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict bx, const float *__restrict by, size_t n);
void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                     float *__restrict o, size_t n);
/* Same semantics as the per-element fv_vec2_normalize: len <= FV_EPS outputs (0, 0). */
void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           float *__restrict ox, float *__restrict oy, size_t n);

/* Fused pipeline entry: out = normalize(a * s + b) in **one pass** (the three calls are three passes over memory).
 * Byte accounting per element, from the three signatures:
 *   scale_soa (8 read + 8 written + 8 RFO) + add_soa (16 read + 8 written + 8 RFO) + normalize_soa (8 read + 8 written + 8 RFO) = 80 B,
 *   fused: read a and b (16) + write out (8) + RFO (8) = 32 B, so 80/32 = 2.5x is the nominal ceiling.
 * This is a *composition* entry, not a new arithmetic path: the evaluation order is the one the three calls imply
 * (t = a * s, u = t + b, then the normalize steps), and the zero-vector guard is the mask of fv_vec2_normalize_soa word for word,
 * including the propagation of inf / NaN. The result is therefore **memcmp-identical** to that three-call sequence for every n
 * and every tail length, so the guarantee of the three entries backs it and no separate error analysis applies.
 * The only exception is the library-wide one stated in the matrix section above ("when two different NaN bit patterns meet
 * inside one arithmetic instruction, which NaN comes out follows the instruction's operand order"). Fusion neither adds nor
 * widens it: it is a property of the hardware plus the compiler's freedom to allocate commutative operands, and it already
 * separates two builds of fv_vec2_normalize_soa itself (measured while writing this entry: the inline form and the
 * in-library form disagree on that input, as do two instances of the same intrinsic in one binary). tests/test_madd_norm.cpp
 * keeps its strict corpora at one NaN bit pattern and counts the anything-goes corpus separately, so nothing is hidden.
 * Measured here (Ryzen 7 5700U, L3 4 MiB per CCX, taskset-pinned single core, default tier -O3 -fno-math-errno
 * -ffp-contract=off without -march, same-process alternating A/B, 11 paired repetitions x 3 runs, median of the paired
 * ratios, same-round A/A floor 0.97-1.07): n = 8192 -> 1.79x, n = 1048576 (40 MiB working set) -> 2.52x,
 * n = 8388608 (320 MiB) -> 2.51x. That is the order of the 2.5x the byte accounting predicts; read a single point as
 * "clears the same-round floor", not as a precise multiplier. The _stream form reaches 3.53x / 3.93x at the two latter
 * tiers and 0.96-1.04x at n = 8192, i.e. no gain there -- the known behaviour of an NT store inside L2. */
void fv_vec2_madd_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict bx, const float *__restrict by,
                                fv_real s, float *__restrict ox, float *__restrict oy, size_t n);
/* Same computation, output written with NT stores: the same _stream convention as the entries below (explicit caller choice,
 * the kernel peels the output to 32 B itself, sfence before returning, bit-identical result).
 * Three passes become one, but the result is written only once -- so this form saves one more 8 B/element of RFO than the
 * normal fused form does, and it is the only version whose NT store pays off from a smaller working set. */
void fv_vec2_madd_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict bx, const float *__restrict by,
                                       fv_real s, float *__restrict ox, float *__restrict oy, size_t n);

/* ---------------- memory-traffic entry points (P2) ----------------
 * _stream: results are written with NT stores (bypassing the cache, producing no RFO read). **Chosen explicitly by the caller**;
 *   the library never switches automatically -- when the working set sits in L2/L3 it is 2-8x slower than the normal version (0.13-0.51x at L2 sizes on this machine).
 *   The 32 B alignment of the output pointer is handled by an internal peel, so the caller **need not** guarantee any alignment;
 *   the function has already executed an sfence before returning, so later reads always see the results.
 * _mask_inplace: in-place sparse add; the **value** of an element with mask[i] == 0 stays as it was (the same value semantics as
 *   _mask's "write a[i] back", but with no second buffer required). What is promised is value semantics: it is **not promised** that a non-hit
 *   element is left physically unwritten, and it is not a synchronization primitive. **Measured: no throughput advantage at all over _mask** (gain -0.10..+0.00, every
 *   point below the A/A floor) -- writing in place produces no RFO and saves no traffic; its value is the form, not speed.
 * The detailed contract (with the measured numbers and the capability boundary) is in the last two parts of this header's "batch entry points (P1)" section. */
void fv_vec2_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n);
void fv_vec2_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n);
void fv_vec2_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                             fv_vec2 *__restrict o, size_t n);
void fv_vec2_lerp_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                fv_real t, fv_vec2 *__restrict o, size_t n);
void fv_vec2_madd_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                            fv_real s, fv_vec2 *__restrict o, size_t n);
void fv_vec2_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                            const float *__restrict bx, const float *__restrict by,
                            float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_scale_soa_stream(const float *__restrict ax, const float *__restrict ay, fv_real s,
                              float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  float *__restrict ox, float *__restrict oy, size_t n);
/* In-place sparse add: a[i] += b[i] when mask[i] != 0, otherwise a[i] keeps its value. a is read-write,
 * and the only constraint is that b must not overlap a (nor may mask overlap either of them). */
void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                 const uint8_t *__restrict mask, size_t n);

/* ---------------- 3D / 4D batch entry points (P3) ----------------
 * Form selection (see the "layout" part of the contract above):
 *   3D has SoA only (fv_vec3 stays 12 B, no padding added), with signatures of 4-9 pointers + n
 *   (add_soa 9; scale_soa / dot_soa / normalize_soa 6 each; len_soa 4);
 *   4D AoS uses _x2 (2 x 16 B elements = 32 B = one ymm, full with no lane wasted); SoA is offered only for
 *   the three reduction entry points dot / len / normalize (add/sub/scale already repack nothing in AoS and need no SoA version).
 * The bit-level promise is exactly that of the corresponding vec2 entry points: apart from dot_soa (different reduction order, see above),
 * all of them are bit-identical to the corresponding scalar API. Pointers carry __restrict and by contract must not overlap. */

void fv_vec3_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az,
                     const float *__restrict bx, const float *__restrict by,
                     const float *__restrict bz,
                     float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n);
void fv_vec3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                       const float *__restrict az, fv_real s,
                       float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n);
/* _stream (NT store) variants of the entries above: the same bit-identical result, only the store instruction differs.
 * Same convention as the vec2 _stream entries: caller's explicit choice, internal 32 B peel, sfence before returning. */
void fv_vec3_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                            const float *__restrict az,
                            const float *__restrict bx, const float *__restrict by,
                            const float *__restrict bz,
                            float *__restrict ox, float *__restrict oy, float *__restrict oz,
                            size_t n);
void fv_vec3_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                              const float *__restrict az, fv_real s,
                              float *__restrict ox, float *__restrict oy,
                              float *__restrict oz, size_t n);
/* The only vec3 entry point not promised to be bit-identical: multiple accumulators plus one final reduction (the same scheme as fv_vec2_dot_soa). */
float fv_vec3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, size_t n);
void fv_vec3_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, float *__restrict o, size_t n);
/* Same semantics as fv_vec3_normalize: len <= FV_EPS outputs (0, 0, 0). */
void fv_vec3_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az,
                           float *__restrict ox, float *__restrict oy, float *__restrict oz,
                           size_t n);
void fv_vec3_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az, float *__restrict ox,
                                  float *__restrict oy, float *__restrict oz, size_t n);

/* 4D AoS: 2 elements at a time (16 B per element, no padding, each loadu reads 32 B = 2 contiguous elements). */
void fv_vec4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n);
void fv_vec4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n);
void fv_vec4_scale_x2(const fv_vec4 *__restrict a, fv_real s, fv_vec4 *__restrict o, size_t n);
void fv_vec4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real t,
                     fv_vec4 *__restrict o, size_t n);
void fv_vec4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real s,
                     fv_vec4 *__restrict o, size_t n);

/* _stream (NT store) variants of the entries above: the same bit-identical result, only the store instruction differs.
 * Same convention as the vec2 _stream entries: caller's explicit choice, internal 32 B peel, sfence before returning. */
void fv_vec4_add_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                           fv_vec4 *__restrict o, size_t n);
void fv_vec4_sub_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                           fv_vec4 *__restrict o, size_t n);
void fv_vec4_scale_x2_stream(const fv_vec4 *__restrict a, fv_real s,
                             fv_vec4 *__restrict o, size_t n);
void fv_vec4_lerp_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real t,
                            fv_vec4 *__restrict o, size_t n);
void fv_vec4_madd_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real s,
                            fv_vec4 *__restrict o, size_t n);

/* 4D SoA: 4 component arrays (8 for dot) + n. */
float fv_vec4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az, const float *__restrict aw,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, const float *__restrict bw, size_t n);
void fv_vec4_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, const float *__restrict aw,
                     float *__restrict o, size_t n);
/* Same semantics as fv_vec4_normalize: len <= FV_EPS outputs (0, 0, 0, 0). */
void fv_vec4_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az, const float *__restrict aw,
                           float *__restrict ox, float *__restrict oy,
                           float *__restrict oz, float *__restrict ow, size_t n);
void fv_vec4_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az, const float *__restrict aw,
                                  float *__restrict ox, float *__restrict oy,
                                  float *__restrict oz, float *__restrict ow, size_t n);

/* ---------------- matrix point transform (AoS, 2 fv_vec4 at a time) ----------------
 * out[i] = M * in[i], row-major, evaluating only the first 3 rows (row 3 must be (0,0,0,1); see the fv_mat4 typedef).
 * The two entry points differ only in whether the translation column takes part: points transforms with w = 1 (translation included, output w = 1),
 * dirs transforms with w = 0 (translation ignored, output w = 0).
 * Bit-level promise: bit-identical to the scalar expression ((m0*x + m1*y) + m2*z) + m3 (dirs drops the last term),
 * for any n and any tail length. Why there is only AoS and no SoA version, and the instruction-count evidence for the kernel,
 * are in the "matrix point transform (AoS)" part at the end of this section -- that is a measured conclusion, not a to-do. */
void fv_mat4_transform_points(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                              fv_vec4 *__restrict out, size_t n);
void fv_mat4_transform_dirs(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                            fv_vec4 *__restrict out, size_t n);

/* _stream (NT store) variants of the entries above: the same bit-identical result, only the store instruction differs.
 * Same convention as the vec2 _stream entries: caller's explicit choice, internal 32 B peel, sfence before returning. */
void fv_mat4_transform_points_stream(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                                     fv_vec4 *__restrict out, size_t n);
void fv_mat4_transform_dirs_stream(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                                   fv_vec4 *__restrict out, size_t n);

/* ---------------- integrator primitives (batch form) ----------------
 * Why this group exists: "x[i] += v[i] * s" is the innermost loop of every explicit integrator, and the batch API had no
 * single-pass form for it -- fv_vec2_madd's third operand is a **scalar**, so a*b + c with a *vector* c could only be
 * written as fv_vec2_scale_x4 into a temporary followed by fv_vec2_add_x4. The byte accounting decides the ceiling, and
 * it is computed per entry shape (B/element, RFO = the read-for-ownership of a line that is written):
 *   addmul, in place:      read v 8 + read x 8 + write x 8 = 24 (the in-place store hits the line the read just brought
 *                          in, so no RFO)  vs the two-pass workaround scale_x4 (8 + 8 + 8 RFO = 24) + add_x4
 *                          (16 + 8 + 8 RFO = 32) = 56                                -> 56/24 = 2.33x
 *   maddv, fresh output:   read a 8 + b 8 + c 8 + write o 8 + RFO 8 = 40  vs a hand-written multiply pass into a
 *                          temporary plus add_x4: (16 + 8 + 8 RFO) + (16 + 8 + 8 RFO) = 64   -> 64/40 = 1.60x
 *   maddv _stream:         the same 32 B, RFO dropped                                      -> 64/32 = 2.00x
 * The two workarounds do not carry the same baseline: addmul's scale->add is 56 B, while maddv's multiply pass pays its
 * own RFO as well and is 64 B. (Quoting 56 B for both would understate the maddv baseline.)
 * Measured (Ryzen 7 5700U, L3 4 MiB/CCX, taskset-pinned single core, default tier -O3 -fno-math-errno -ffp-contract=off
 * without -march, same-process interleaved A/B, 11 paired rounds, median of the paired ratios, same-round A/A control;
 * three independent runs, quoted as the range of their medians, because a single point is not a precise multiplier):
 *   n = 8192      addmul 1.73-2.18x   maddv 1.60-1.67x   maddv_stream 1.16-1.37x
 *   n = 1048576   addmul 2.05-2.17x   maddv 1.61-1.79x   maddv_stream 2.16-2.23x
 *   n = 8388608   addmul 1.92-2.03x   maddv 1.59-1.76x   maddv_stream 2.34-2.41x
 *   The same-round A/A control sat at 0.97-1.07 in the quiet rounds; the rounds whose A/A minimum fell to 0.17 / 0.44 are
 *   where the wide 8192 addmul range comes from. Quote a tier, not a point.
 * The in-place form therefore lands at 74-95% of its 2.33x ceiling and maddv at 99-112% of its 1.60x (inside the
 * documented "20% above the ceiling counts as noise" band). The stream form reaches its 2.00x ceiling only from n = 1M up
 * (2.16-2.41x); at 8192 it is 1.16-1.37x, and its absolute cost there is about twice the non-stream fused form's -- an NT
 * store inside L2 is simply the wrong tool, so a caller who is not clearly above L3 must use the plain entry. The ratio
 * against the workaround only says "still faster than two passes", not "this is the form to pick".
 * Read a single point as "cleared the same-round A/A floor"; reproduce with make bench-addmul.
 * Hence addmul has **no** _stream variant: an NT store only avoids the RFO of a cache line the core does not already hold,
 * while an in-place read-modify-write reads (and therefore owns) that line -- a stream form would add the alignment peel
 * and the sfence for nothing. maddv writes a separate output buffer, so it does pay the RFO and the _stream form pays off.
 * In-place contract: x is read and written, v must not overlap x (both carry __restrict), and n == 0 touches no pointer
 * (NULL is allowed for every argument). Every entry in this group evaluates "multiply first, then add" inside the
 * fp-contract=off pragma, so the result is bit-identical to the per-value fv_vec2_addmul / fv_vec3_addmul /
 * fv_vec4_addmul / fv_vec2_maddv above, for every n and every tail length. */
void fv_vec2_addmul_x4(fv_vec2 *__restrict x, const fv_vec2 *__restrict v, fv_real s, size_t n);
void fv_vec2_addmul_soa(float *__restrict x, float *__restrict y,
                        const float *__restrict vx, const float *__restrict vy, fv_real s, size_t n);
void fv_vec3_addmul_soa(float *__restrict x, float *__restrict y, float *__restrict z,
                        const float *__restrict vx, const float *__restrict vy,
                        const float *__restrict vz, fv_real s, size_t n);
void fv_vec4_addmul_x2(fv_vec4 *__restrict x, const fv_vec4 *__restrict v, fv_real s, size_t n);

/* Non-in-place: o[i] = a[i]*b[i] + c[i] -- three vector operands, which the library did not offer before (madd takes a
 * scalar third operand). a / b / c / o must not overlap.
 * The _stream form writes the output with NT stores and follows the same three-phase structure as every other entry of
 * that family: head peel with normal stores up to the 32 B boundary, body _mm256_stream_ps, tail with normal stores, and
 * _mm_sfence() before returning. So the caller needs no alignment (the kernel peels; a pointer that can never reach the
 * boundary falls back to storeu over the whole span) and the result is bit-identical to the non-stream form. */
void fv_vec2_maddv_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                      const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n);
void fv_vec2_maddv_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                             const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n);

/* ---------------- quaternion / mat3 SoA batch entry points (P6) ----------------
 *
 * Form: component arrays (SoA), 8 elements per 256-bit register, like fv_vec3_add_soa. A quaternion is
 * four arrays (qx, qy, qz, qw); a 3x3 matrix is nine (m00..m22), so the transform takes the nine matrix
 * arrays plus three input and three output arrays. This is the layout a physics solver already wants:
 * each lane is one independent body, and there is no per-element packing/shuffling.
 *
 * Bit-exactness: every entry is **bit-identical** to its per-value counterpart for every n and every
 * tail length -- the kernels evaluate the same expressions in the same association order as
 * src/flashvec_quat.cpp / src/flashvec_mat3.cpp (listed in the scalar section above), and sit inside the
 * same fp-contract=off pragma. The one documented NaN-payload exception applies unchanged: when one
 * arithmetic instruction receives two NaN operands with different bit patterns, which payload survives
 * follows the instruction's operand order; at most one NaN bit pattern in the corpus keeps the promise
 * strict (the same rule as the vec/mat4 entries).
 *
 * Implementation form (measured decision, this round): quaternion multiply / rotate and mat3 transform
 * are **plain C** inside target("avx2") functions, not hand-written intrinsics. The project's own
 * evidence is that hand-written SIMD only pays where the compiler fails to vectorize an elementwise loop
 * (mat4 transform: 440 instructions / 14 vpermps vs 107 / 3); these three are plain 4- and 9-term
 * elementwise loops, and GCC vectorizes them at 256 bits with zero shuffles (measured 40 / 40 / 21 AVX2
 * instructions, 0 vpermps/vshufps), so `make bench-physics` measures that version first.
 * fv_quat_normalize_soa is the exception: the compiler kept sqrtss/divss and a branch there (0 ymm
 * instructions), so it has an explicit AVX2 kernel with a masked division and a blend to the identity.
 * The rule is the same as the rest of the library: hand-write only after the compiler fails, and keep
 * the census reproducible (`make bench-physics`, `objdump`).
 *
 * Preconditions: every pointer is __restrict and the arrays must not overlap; n == 0 touches no pointer
 * (NULL is allowed for every argument); fv_quat_rotate_vec3_soa assumes unit quaternions (same rule as
 * the per-value function); fv_quat_normalize_soa maps a near-zero quaternion to the identity.
 */
void fv_quat_mul_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, const float *__restrict aw,
                     const float *__restrict bx, const float *__restrict by,
                     const float *__restrict bz, const float *__restrict bw,
                     float *__restrict ox, float *__restrict oy,
                     float *__restrict oz, float *__restrict ow, size_t n);
void fv_quat_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az, const float *__restrict aw,
                           float *__restrict ox, float *__restrict oy,
                           float *__restrict oz, float *__restrict ow, size_t n);
void fv_quat_rotate_vec3_soa(const float *__restrict qx, const float *__restrict qy,
                             const float *__restrict qz, const float *__restrict qw,
                             const float *__restrict vx, const float *__restrict vy,
                             const float *__restrict vz,
                             float *__restrict ox, float *__restrict oy,
                             float *__restrict oz, size_t n);
void fv_mat3_transform_vec3_soa(const float *__restrict m00, const float *__restrict m01,
                                const float *__restrict m02, const float *__restrict m10,
                                const float *__restrict m11, const float *__restrict m12,
                                const float *__restrict m20, const float *__restrict m21,
                                const float *__restrict m22,
                                const float *__restrict vx, const float *__restrict vy,
                                const float *__restrict vz,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, size_t n);

/* 9x9 SoA pair: out = a*b, nine component arrays per operand. The expressions are the scalar
 * fv_mat3_mul's, lane by lane: out[i][j] = ((ai0*b0j + ai1*b1j) + ai2*b2j). Bit-identical to the
 * per-value function for every n and tail. */
void fv_mat3_mul_soa(const float *__restrict a00, const float *__restrict a01,
                     const float *__restrict a02, const float *__restrict a10,
                     const float *__restrict a11, const float *__restrict a12,
                     const float *__restrict a20, const float *__restrict a21,
                     const float *__restrict a22,
                     const float *__restrict b00, const float *__restrict b01,
                     const float *__restrict b02, const float *__restrict b10,
                     const float *__restrict b11, const float *__restrict b12,
                     const float *__restrict b20, const float *__restrict b21,
                     const float *__restrict b22,
                     float *__restrict o00, float *__restrict o01,
                     float *__restrict o02, float *__restrict o10,
                     float *__restrict o11, float *__restrict o12,
                     float *__restrict o20, float *__restrict o21,
                     float *__restrict o22, size_t n);

/* Per-element inverse with the same scale-relative singularity test as fv_mat3_inverse: a singular
 * lane writes nine zeros and is not reported separately (the batch form has no per-lane status).
 * Bit-identical to the per-value function for every n and tail. */
void fv_mat3_inverse_soa(const float *__restrict m00, const float *__restrict m01,
                         const float *__restrict m02, const float *__restrict m10,
                         const float *__restrict m11, const float *__restrict m12,
                         const float *__restrict m20, const float *__restrict m21,
                         const float *__restrict m22,
                         float *__restrict o00, float *__restrict o01,
                         float *__restrict o02, float *__restrict o10,
                         float *__restrict o11, float *__restrict o12,
                         float *__restrict o20, float *__restrict o21,
                         float *__restrict o22, size_t n);

/* ---------------- _fast approximate normalize entries (P7) ----------------
 * rsqrtps + one Newton-Raphson iteration instead of sqrtps + divps. The result is **approximate**
 * (relative error about 1-2 ulp on AVX2; the scalar fallback uses 1/sqrtf and is exact) and is
 * deliberately not bit-identical to the exact normalize entries -- these are new names in a separate
 * family, and the exact entries above keep their promise unchanged.
 *
 * Zero guard: the same values as the exact entries (a vec entry writes zeros, the quat entry writes
 * the identity (0,0,0,1)) when the squared length is <= FV_EPS_SQR. That guard is a blend, so no NaN
 * escapes from the 0 * inf product. Preconditions (__restrict, non-overlap, n == 0 touches no
 * pointer, NULL allowed) are the same as the exact SoA entries. */
void fv_vec2_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec3_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, size_t n);
void fv_vec4_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az, const float *__restrict aw,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, float *__restrict ow, size_t n);
void fv_quat_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az, const float *__restrict aw,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, float *__restrict ow, size_t n);
/* Fused approximate pipeline: out = normalize(a*s + b) with rsqrtps + one Newton step. Same
 * approximate contract as the normalize_fast entries; the zero guard writes (0,0). */
void fv_vec2_madd_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict bx, const float *__restrict by,
                                    fv_real s, float *__restrict ox, float *__restrict oy,
                                    size_t n);

#endif /* FV_USE_EXTERN || FV_BUILD: the public names are declarations only */

/* ============================ implementation of the batch entry points ============================
 * Compiled only when "this header provides the implementation":
 *   FV_USE_EXTERN -- the whole section disappears (the implementation comes from libflashvec.a and the header stays declaration-only);
 *   FV_BUILD      -- compiled (the library needs it), but without the public names below (defined in src/flashvec_batch.cpp);
 *   default       -- compiled, and a final static inline layer forwards to the public names.
 * There is one implementation only: every body lives in fv_i_*, and the two forms differ only in who wires it to the public names.
 *
 * Why two sections -- the single dividing line is "does this function contain an expression where a product feeds straight into an addition":
 *   first section (no contraction risk, outside the pragma)
 *     add / sub / scale / add_mask / add_idx / scale_idx / add_soa / scale_soa:
 *     the expressions are a single addition or a single multiplication, where FMA contraction means nothing. Being outside the pragma, the consumer's
 *     default compiler options fully inline the entry point (including the n < 4 / n < 8 tail loops) into its own code.
 *     Measured (putting the equivalent kernel under -march=native -mfma -ffp-contract=fast): vfmadd = 0.
 *   second section (contraction risk, inside the pragma as a whole)
 *     lerp = a + (b-a)*t, madd = a + b*s, the x*x + y*y of len / normalize,
 *     and the ax*bx + ay*by accumulation of dot: measured under the same options, the equivalent kernel without the pragma shows 1
 *     vfmadd per occurrence, and 0 everywhere with the pragma. The cost is that these 5 entry points do not inline under the consumer's default options
 *     (one intra-TU call per batch) -- deliberately trading "one call" for "the bit-level promise".
 *
 * Why target("avx2") alone is not enough: the target attribute only sets this function's own instruction-set baseline and cannot control the -mfma /
 * -march=native coming from the caller's command line. Measured: under -march=native -mfma -ffp-contract=fast,
 * a kernel carrying target("avx2") still contracts into vfmadd. In the consumer scenario the pragma is the only line of defence.
 */
#if !defined(FV_USE_EXTERN)

#if FV_I_INTRINSICS
#define FV_I_TARGET_AVX2 __attribute__((target("avx2")))
#if FV_I_AVX512
#define FV_I_TARGET_AVX512 __attribute__((target("avx512f")))
#endif
#endif

/* ======== first section: no mul+add contraction risk (outside the pragma -> fully inlined under the consumer's default options) ======== */

/* Scalar expressions (shared by non-AVX2 machines and every tail): the expressions correspond word for word to the per-element API in
 * src/flashvec.cpp; they are a single addition or multiplication, with no contractible mul+add. */
static inline fv_vec2 fv_i_s_add(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    return r;
}


static inline fv_vec2 fv_i_s_sub(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    return r;
}


static inline fv_vec2 fv_i_s_scale(fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = v.x * s;
    r.y = v.y * s;
    return r;
}

static inline fv_vec4 fv_i_s_add4(fv_vec4 a, fv_vec4 b)
{
    fv_vec4 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    r.w = a.w + b.w;
    return r;
}


static inline fv_vec4 fv_i_s_sub4(fv_vec4 a, fv_vec4 b)
{
    fv_vec4 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
    r.w = a.w - b.w;
    return r;
}


static inline fv_vec4 fv_i_s_scale4(fv_vec4 v, fv_real s)
{
    fv_vec4 r;
    r.x = v.x * s;
    r.y = v.y * s;
    r.z = v.z * s;
    r.w = v.w * s;
    return r;
}

#if FV_I_INTRINSICS

/* AVX2 kernels: only the single operations add / sub / mul, with no contractible mul+add.
 * 4 fv_vec2 = 32 B = one ymm: loadu/storeu consume contiguous elements directly, with zero repacking.
 * fv_vec2 is only 4 B aligned, so loadu/storeu are used throughout (_mm256_load_ps may not be used). */

/* 4 fv_vec2 = 32 B = one ymm: loadu/storeu consume contiguous elements directly, with zero repacking.
 * fv_vec2 is only 4 B aligned, so loadu/storeu are used throughout (_mm256_load_ps may not be used).
 * The target names avx2 only, not fma: FMA instructions are simply unavailable, so contraction is excluded at the instruction level. */

FV_I_TARGET_AVX2 static inline void fv_i_add4_avx2(const fv_vec2 *__restrict a,
                                                   const fv_vec2 *__restrict b,
                                                   fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_sub4_avx2(const fv_vec2 *__restrict a,
                                                   const fv_vec2 *__restrict b,
                                                   fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_sub_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_scale4_avx2(const fv_vec2 *__restrict a, fv_real s,
                                                     fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_mul_ps(va, vs));
    }
}


/* Mask variant: 4 elements at a time (one ymm holds 4 fv_vec2 = 8 float lanes).
 * The key point: the mask is given **per element**, and one element spans two float lanes, so each mask byte
 * must be spread over two 32-bit lanes (x and y share one mask bit); it must not be tiled one byte per lane.
 * Then blendvps picks a (mask 0) or a+b (mask non-zero) and the whole vector is written back to o. */
FV_I_TARGET_AVX2 static inline void fv_i_add_mask4_avx2(const fv_vec2 *__restrict a,
                                                        const fv_vec2 *__restrict b,
                                                        fv_vec2 *__restrict o,
                                                        const uint8_t *__restrict mask, size_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        memcpy(&mb, mask + i, sizeof(mb));
        if (mb == 0u) {
            /* Block prescan: four mask bytes all zero -> the block's whole output is a's original value (the write-a[i]-back-outside-the-mask semantics),
             * so this block needs neither a load of b nor a blend.
             * The test must be "the whole u32 is zero" and must not be (mb & 0x01010101u) == 0: that tests only
             * each byte's lowest bit, so a mask byte of 0x02 / 0x80 / 0xFE -- **even, or with only a high bit set** -- is misjudged as
             * "a miss", while the mask semantics are "a non-zero byte is a hit". Testing a u32 for zero does not depend on byte order. */
            _mm256_storeu_ps((float *)(void *)(o + i),
                             _mm256_loadu_ps((const float *)(const void *)(a + i)));
            continue;
        }
        __m128i b4 = _mm_cvtsi32_si128((int)mb);  /* low 4 bytes = the 4 mask bytes */
        __m128i w = _mm_cvtepu8_epi16(b4);        /* 4 bytes -> 4 words (0..255) */
        __m128i wd = _mm_unpacklo_epi16(w, w);    /* duplicate each word */
        __m256i d = _mm256_cvtepi16_epi32(wd);    /* 8 dwords: m0 m0 m1 m1 m2 m2 m3 m3 */
        __m256 m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(d, zero));
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_blendv_ps(va, _mm256_add_ps(va, vb), m));
    }
}


/* SoA: 8 floats at a time. */
FV_I_TARGET_AVX2 static inline void fv_i_add_soa_avx2(const float *__restrict ax,
                                                      const float *__restrict ay,
                                                      const float *__restrict bx,
                                                      const float *__restrict by,
                                                      float *__restrict ox,
                                                      float *__restrict oy, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_add_ps(_mm256_loadu_ps(ax + i),
                                               _mm256_loadu_ps(bx + i)));
        _mm256_storeu_ps(oy + i, _mm256_add_ps(_mm256_loadu_ps(ay + i),
                                               _mm256_loadu_ps(by + i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_scale_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay, fv_real s,
                                                        float *__restrict ox,
                                                        float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_mul_ps(_mm256_loadu_ps(ax + i), vs));
        _mm256_storeu_ps(oy + i, _mm256_mul_ps(_mm256_loadu_ps(ay + i), vs));
    }
}

/* ---------------- streaming (NT store) kernels: first section ----------------
 * Three-phase structure, correct for any pointer alignment:
 *   head peel -- element by element with **normal stores** up to the 32 B boundary (so any pointer alignment can enter the body);
 *   body      -- _mm256_stream_ps (a write-combining store that bypasses the cache and produces no RFO), requiring 32 B alignment;
 *   tail      -- the part under one vector keeps using normal stores.
 * The address of _mm256_stream_ps must be 32 B aligned or #GP; so the alignment test before the body is a hard requirement,
 * not an optimization. The input side still uses loadu: only the **output** needs alignment.
 * Trailing _mm_sfence(): an NT store is weakly ordered (the data in the WC buffer may still be invisible to other cores),
 * and returning without the sfence lets the caller's subsequent reader thread see stale values. */
FV_I_TARGET_AVX2 static inline void fv_i_add4_stream_avx2(const fv_vec2 *__restrict a,
                                                          const fv_vec2 *__restrict b,
                                                          fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    /* Head: with an 8 B stride, pushing o to the 32 B boundary takes at most 3 steps */
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        /* The 32 B boundary is unreachable: only possible for an fv_vec2 array whose o is not even 8 B aligned (a stride of 8 B cannot
         * remove a 4 B phase offset). The whole span goes to the normal storeu kernel -- those head elements are recomputed,
         * the result is bit-identical, only the NT gain is lost. */
        fv_i_add4_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_sub4_stream_avx2(const fv_vec2 *__restrict a,
                                                          const fv_vec2 *__restrict b,
                                                          fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_sub4_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_sub_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_scale4_stream_avx2(const fv_vec2 *__restrict a, fv_real s,
                                                            fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_scale4_avx2(a, s, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_mul_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)), vs));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
    _mm_sfence();
}


/* The SoA streaming kernel is written per output component array: each _soa streaming entry calls it twice (once for x, once for y).
 * That way ox and oy peel independently and each decides on its own whether NT is possible, so they never drag each other
 * down when their phases differ (sharing one i would require both to align **at the same time**). */
FV_I_TARGET_AVX2 static inline void fv_i_add1_stream_avx2(const float *__restrict a,
                                                          const float *__restrict b,
                                                          float *__restrict o, size_t n)
{
    size_t i = 0;
    /* Float stride 4 B: 7 steps cover every phase within 32 B, so it is always reachable for normally aligned pointers */
    for (; i < n && i < 7u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = a[i] + b[i];
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        /* Unreachable: only possible for the pathological input where o is not even 4 B aligned */
        for (; i + 8 <= n; i += 8) {
            _mm256_storeu_ps(o + i, _mm256_add_ps(_mm256_loadu_ps(a + i),
                                                  _mm256_loadu_ps(b + i)));
        }
    } else {
        for (; i + 8 <= n; i += 8) {
            _mm256_stream_ps(o + i, _mm256_add_ps(_mm256_loadu_ps(a + i),
                                                  _mm256_loadu_ps(b + i)));
        }
    }
    for (; i < n; ++i) {
        o[i] = a[i] + b[i];
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_scale1_stream_avx2(const float *__restrict a, fv_real s,
                                                            float *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 7u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = a[i] * s;
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        for (; i + 8 <= n; i += 8) {
            _mm256_storeu_ps(o + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), vs));
        }
    } else {
        for (; i + 8 <= n; i += 8) {
            _mm256_stream_ps(o + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), vs));
        }
    }
    for (; i < n; ++i) {
        o[i] = a[i] * s;
    }
    _mm_sfence();
}


/* The two in-place sparse add kernels -- the results are **bit-identical** (a hit lane is always the scalar a[i] + b[i]);
 * they differ only in "does a non-hit lane perform a memory access". The routing is decided once per batch by fv_i_mask_inplace_use_dense().
 *
 * Sparse branch (sampled density <= FV_I_MIP_DENSE_NUM/DEN): block prescan + per-element scalar.
 *   A block whose four mask bytes are all zero is skipped whole -- not even a's cache line is touched (at 5% hits, 95%^4
 *   = 81% of 4-element blocks are all-zero), and inside a block with hits only the hit elements get an 8 B read-modify-write.
 *   There is no fast path for "four hits in a row, go vector" here: at a density <= 1/8 an all-hit block has a probability of about 0.02%,
 *   and one more branch for it is not worth it.
 *   **Note**: "which branch pays off" has no measured support -- against _mask this entry is **INCONCLUSIVE at every tier**
 *   (gain -0.10..+0.00, all below the same-round A/A floor; see the capability boundary below). So the routing affects only the shape of the code:
 *   the two paths have **no measurable difference**; the sentence above only declines to add a branch for a 0.02% event, it is not a performance claim. */
FV_I_TARGET_AVX2 static inline void fv_i_branch_add4_prescan_avx2(fv_vec2 *__restrict a,
                                                                  const fv_vec2 *__restrict b,
                                                                  const uint8_t *__restrict mask,
                                                                  size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        size_t k;
        memcpy(&mb, mask + i, sizeof(mb));
        if (mb == 0u) {
            continue; /* none of the four elements hits: neither cache line of a / b is touched at all */
        }
        for (k = 0; k < 4u; ++k) {
            if (mask[i + k] != 0u) {
                a[i + k] = fv_i_s_add(a[i + k], b[i + k]);
            }
        }
    }
}


/* Dense branch (sampled density > FV_I_MIP_DENSE_NUM/DEN): whole-block 32 B load a / load b, blend the result,
 * and one storeu writes it back into a. No data-dependent branches and no VMASKMOVPS (that family is microcoded and slow on Zen 2;
 * an **early probe** measured a maskstore as 1.2-3.8x slower than the whole-vector storeu used here at the same density, but the hardened fixtures this round **did not re-measure it**,
 * so read it as an order of magnitude -- that family has been deleted from the library and this sentence only explains why it is not used).
 * **A non-hit lane is written back with a's original value** (the blend picks va): exactly equivalent under value semantics, but the write does happen physically.
 * So this function promises only "the value of a non-hit element is unchanged", not "it is left physically unwritten" -- see the contract. */
FV_I_TARGET_AVX2 static inline void fv_i_blend_add4_inplace_avx2(fv_vec2 *__restrict a,
                                                                 const fv_vec2 *__restrict b,
                                                                 const uint8_t *__restrict mask,
                                                                 size_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        __m128i b4;
        __m128i w;
        __m128i wd;
        __m256i d;
        __m256 m;
        __m256 va;
        __m256 vb;
        memcpy(&mb, mask + i, sizeof(mb));
        b4 = _mm_cvtsi32_si128((int)mb);
        w = _mm_cvtepu8_epi16(b4);
        wd = _mm_unpacklo_epi16(w, w);
        d = _mm256_cvtepi16_epi32(wd);
        m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(d, zero));
        va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(a + i),
                         _mm256_blendv_ps(va, _mm256_add_ps(va, vb), m));
    }
}


/* ---- fv_vec4 AoS kernels (_x2): 2 x 16 B elements = 32 B = one ymm ----
 * There is **no such thing** here as "assembling two 16 B elements into one ymm": fv_vec4 has no padding bytes,
 * two adjacent array elements are already 32 B of contiguous memory, and the 8 lanes one 256-bit loadu brings in are exactly
 * (a[i].x, a[i].y, a[i].z, a[i].w, a[i+1].x, ...) in the naturally correct order; the output is likewise one 256-bit storeu.
 * So _x2 is **zero repacking**, just like vec2's _x4 (no insertf128-style assembly instruction is needed).
 * The two-128-bit-loads + _mm256_insertf128_ps form is needed only in a layout where the element stride is not 16 B
 * (AoSoA, an array with a stride); only there does assembly cost even come up.
 * Element 16 B, ymm 32 B: exactly 2 elements fill one register, **not one lane wasted** --
 * that is what the name _x2 means (isomorphic to vec2's _x4: the number = how many elements one 256-bit register holds). */
FV_I_TARGET_AVX2 static inline void fv_i_v4_add2_avx2(const fv_vec4 *__restrict a,
                                                    const fv_vec4 *__restrict b,
                                                    fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_sub2_avx2(const fv_vec4 *__restrict a,
                                                    const fv_vec4 *__restrict b,
                                                    fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_sub_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_scale2_avx2(const fv_vec4 *__restrict a, fv_real s,
                                                      fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_mul_ps(va, vs));
    }
}


/* ---- 3D SoA kernels: one independent 8-lane loop per component ----
 * "Does running 3 components over 8 lanes waste lanes?" No: a ymm always carries **8 different elements of the same
 * component** (8x ax[i..i+7]), never "the 3 components of one element plus 1 padding lane".
 * Every component runs at 100% lane utilization; the cost falls on the caller only: the data must be stored split by component,
 * and a tail with n < 8 (at most 7 elements) takes the scalar path. The 12 B AoS layout is the opposite: it does not fit a whole number of ymm,
 * so it either idles 1/4 of the lanes or unpacks with shuffles; hence 3D batches offer SoA only. */
FV_I_TARGET_AVX2 static inline void fv_i_v3_add_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       const float *__restrict bx,
                                                       const float *__restrict by,
                                                       const float *__restrict bz,
                                                       float *__restrict ox,
                                                       float *__restrict oy,
                                                       float *__restrict oz, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_add_ps(_mm256_loadu_ps(ax + i),
                                               _mm256_loadu_ps(bx + i)));
        _mm256_storeu_ps(oy + i, _mm256_add_ps(_mm256_loadu_ps(ay + i),
                                               _mm256_loadu_ps(by + i)));
        _mm256_storeu_ps(oz + i, _mm256_add_ps(_mm256_loadu_ps(az + i),
                                               _mm256_loadu_ps(bz + i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v3_scale_soa_avx2(const float *__restrict ax,
                                                         const float *__restrict ay,
                                                         const float *__restrict az, fv_real s,
                                                         float *__restrict ox,
                                                         float *__restrict oy,
                                                         float *__restrict oz, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_mul_ps(_mm256_loadu_ps(ax + i), vs));
        _mm256_storeu_ps(oy + i, _mm256_mul_ps(_mm256_loadu_ps(ay + i), vs));
        _mm256_storeu_ps(oz + i, _mm256_mul_ps(_mm256_loadu_ps(az + i), vs));
    }
}


/* ---- streaming (NT store) kernels of the first section, 4D AoS _x2 shape ----
 * Word for word the three-phase structure of the vec2 _x4 kernels (head peel -> NT body -> tail); only the stride differs:
 * fv_vec4 is 16 B, so the output phase within 32 B takes two values (0 / 16) and **one** peel step always reaches the boundary.
 * An array that is not even 16 B aligned cannot reach it (a 16 B stride cannot remove an 8 B phase offset): the whole span then
 * goes to the normal storeu kernel -- bit-identical result, only the NT gain is lost. Trailing sfence on every path, for the reason
 * stated above the vec2 _x4 kernels. */
FV_I_TARGET_AVX2 static inline void fv_i_v4_add2_stream_avx2(const fv_vec4 *__restrict a,
                                                             const fv_vec4 *__restrict b,
                                                             fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_add4(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_v4_add2_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_add4(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_sub2_stream_avx2(const fv_vec4 *__restrict a,
                                                             const fv_vec4 *__restrict b,
                                                             fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_sub4(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_v4_sub2_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_sub_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub4(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_scale2_stream_avx2(const fv_vec4 *__restrict a,
                                                               fv_real s,
                                                               fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_scale4(a[i], s);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_v4_scale2_avx2(a, s, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_mul_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)), vs));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale4(a[i], s);
    }
    _mm_sfence();
}


#endif /* FV_I_INTRINSICS: kernels of the first section */

/* The 8 entry points of the first section: run-time dispatch + tail.
 * The whole section has no contractible floating-point expression, so it inlines completely into the caller; the two _idx variants are deliberately scalar
 * (with random access, gather is a microcoded loop on this machine's Zen 2 and usually loses to per-element code; and by contract idx may repeat
 * and may be out of order, with repeated indices writing the same value, independent of the processing order).
 * The kernels carry target("avx2") and so do not match the caller's target; the batch part is therefore still one intra-TU call per batch. */
static inline void fv_i_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                               fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                               fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_sub4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                                 fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_scale4_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                    fv_vec2 *__restrict o, const uint8_t *__restrict mask,
                                    size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add_mask4_avx2(a, b, o, mask, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        if (mask[i] != 0) {
            o[i] = fv_i_s_add(a[i], b[i]);
        } else {
            o[i] = a[i]; /* write a's original value back outside the mask (not "skip the write") */
        }
    }
}


static inline void fv_i_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                   size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                                     fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                     size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_soa(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict bx, const float *__restrict by,
                                float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_add_soa_avx2(ax, ay, bx, by, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}


static inline void fv_i_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                  fv_real s, float *__restrict ox, float *__restrict oy,
                                  size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_scale_soa_avx2(ax, ay, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}

/* ---- streaming entries of the first section ----
 * lerp / madd do have _stream variants (see the second section), but their usual pattern is the opposite of add/sub/scale: an
 * interpolation result is often read again right away by the next step, and then keeping a cached copy pays off much better.
 * The gain of an NT store comes from "this large block just written will not be read again", so for lerp/madd the stream form is
 * the exception, not the default choice: use it only when the result really is a write-once batch.
 * The entry shape is word-for-word identical to the corresponding non-stream entry (bulk + scalar tail); only the kernel called differs,
 * so the tail-length / n == 0 semantics match the non-stream version exactly (for n == 0 no pointer is touched and no sfence is needed). */
static inline void fv_i_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add4_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_sub4_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                                        fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_scale4_stream_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict bx, const float *__restrict by,
                                       float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_add1_stream_avx2(ax, bx, ox, bulk);
            fv_i_add1_stream_avx2(ay, by, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}


static inline void fv_i_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                         fv_real s, float *__restrict ox, float *__restrict oy,
                                         size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_scale1_stream_avx2(ax, s, ox, bulk);
            fv_i_scale1_stream_avx2(ay, s, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}


/* Routing for the in-place sparse add: low density takes "block prescan + per-element" (miss blocks untouched),
 * medium and high density take "whole-block load + blend + storeu" (no data-dependent branches).
 * The test samples only the first 64 mask bytes (a few dozen cycles, negligible against a whole batch); an order-of-magnitude estimate suffices,
 * the two paths give bit-identical results, and a wrong estimate costs performance only, never semantics.
 * The cut-off is set by the three macros below (tunable without touching the function body): the sample length and the density threshold. */
#define FV_I_MIP_DENSE_NUM 1u /* the dense path is taken when the sampled density > NUM/DEN */
#define FV_I_MIP_DENSE_DEN 8u
#define FV_I_MIP_SAMPLE 64u

static inline int fv_i_mask_inplace_use_dense(const uint8_t *mask, size_t n)
{
    const size_t samp = (n < FV_I_MIP_SAMPLE) ? n : FV_I_MIP_SAMPLE;
    size_t hits = 0;
    size_t k;
    for (k = 0; k < samp; ++k) {
        hits += (mask[k] != 0u) ? (size_t)1 : (size_t)0;
    }
    return (hits * FV_I_MIP_DENSE_DEN > samp * FV_I_MIP_DENSE_NUM) ? 1 : 0;
}


/* In-place sparse add: a[i] += b[i] when mask[i] != 0, otherwise a[i] keeps its value.
 * There is no contractible mul+add, so it stays in the first section.
 * The implementation picks one of two kernels by sampled density (see fv_i_mask_inplace_use_dense above):
 * the result does not depend on the path chosen -- a hit lane is bit-identical to the scalar a[i] + b[i] either way.
 * Note that the two branches differ in **storage behaviour**: low density does not touch miss blocks, high density writes the non-hit lanes back with their original value.
 * That is a promise about values, not about storage behaviour (see the " _mask_inplace " part of the contract). */
static inline void fv_i_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                            const uint8_t *__restrict mask, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            if (fv_i_mask_inplace_use_dense(mask, bulk) != 0) {
                fv_i_blend_add4_inplace_avx2(a, b, mask, bulk);
            } else {
                fv_i_branch_add4_prescan_avx2(a, b, mask, bulk);
            }
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        if (mask[i] != 0) {
            a[i] = fv_i_s_add(a[i], b[i]);
        }
    }
}

/* ---- first-section entries for 3D / 4D ----
 * Isomorphic to the vec2 entries above: the bulk part goes to a kernel carrying target("avx2") and the remaining tail uses scalar expressions.
 * There is no expression here where a product feeds straight into an addition (only single +, -, x), so FMA contraction means nothing to them,
 * and they are **deliberately left outside the pragma**: with the consumer's default options the whole section (tail loops included) inlines completely.
 * Tail lengths: _x2 has a bulk of n & ~1 (a tail of 0-1) and _soa a bulk of n & ~7 (a tail of 0-7). */
static inline void fv_i_v4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_add2_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add4(a[i], b[i]);
    }
}


static inline void fv_i_v4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_sub2_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub4(a[i], b[i]);
    }
}


static inline void fv_i_v4_scale_x2(const fv_vec4 *__restrict a, fv_real s, fv_vec4 *__restrict o,
                                  size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_scale2_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale4(a[i], s);
    }
}

/* The vec4 _x2 streaming entries of the first section: word for word the shape of the non-stream entries above
 * (bulk = n & ~1 + scalar tail), only the kernel called differs, so n == 0 / the tail lengths keep the same semantics. */
static inline void fv_i_v4_add_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                         fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_add2_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add4(a[i], b[i]);
    }
}


static inline void fv_i_v4_sub_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                         fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_sub2_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub4(a[i], b[i]);
    }
}


static inline void fv_i_v4_scale_x2_stream(const fv_vec4 *__restrict a, fv_real s,
                                           fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_scale2_stream_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale4(a[i], s);
    }
}


static inline void fv_i_v3_add_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az,
                                 const float *__restrict bx, const float *__restrict by,
                                 const float *__restrict bz,
                                 float *__restrict ox, float *__restrict oy,
                                 float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_add_soa_avx2(ax, ay, az, bx, by, bz, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
        oz[i] = az[i] + bz[i];
    }
}


static inline void fv_i_v3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, fv_real s,
                                   float *__restrict ox, float *__restrict oy,
                                   float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_scale_soa_avx2(ax, ay, az, s, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
        oz[i] = az[i] * s;
    }
}


/* The 3D streaming entries of the first section. add / scale have no mul+add, so they belong here (fully inlined under
 * the consumer's default options) and reuse the per-component float kernels of the SoA family: one call per output array,
 * so ox / oy / oz peel independently and each decides on its own whether NT is possible -- the same convention as
 * fv_i_add_soa_stream / fv_i_scale_soa_stream (sharing one i would require all three to align at the same time). */
static inline void fv_i_v3_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                          const float *__restrict az,
                                          const float *__restrict bx, const float *__restrict by,
                                          const float *__restrict bz,
                                          float *__restrict ox, float *__restrict oy,
                                          float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_add1_stream_avx2(ax, bx, ox, bulk);
            fv_i_add1_stream_avx2(ay, by, oy, bulk);
            fv_i_add1_stream_avx2(az, bz, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
        oz[i] = az[i] + bz[i];
    }
}


static inline void fv_i_v3_scale_soa_stream(const float *__restrict ax,
                                            const float *__restrict ay,
                                            const float *__restrict az, fv_real s,
                                            float *__restrict ox, float *__restrict oy,
                                            float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_scale1_stream_avx2(ax, s, ox, bulk);
            fv_i_scale1_stream_avx2(ay, s, oy, bulk);
            fv_i_scale1_stream_avx2(az, s, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
        oz[i] = az[i] * s;
    }
}


/* ======== second section: mul+add contraction risk (inside the pragma -> one intra-TU call per batch) ======== */

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#if !defined(FV_BUILD)
/* ------------------------------------------------------------------                                                              *
 * Scalar API definitions of the default form (P0): one body per public scalar name, static inline, inside the same
 * fp-contract=off pragma as the batch kernels. The .cpp files keep the extern symbols for FV_USE_EXTERN / ABI and
 * compile with FV_BUILD, which turns this block off; the two forms are compared bit-for-bit by tests/test_scalar_parity.
 * ------------------------------------------------------------------                                                              */
static inline fv_real fv_i_sc_sqrt(fv_real v)
{
    return sqrtf(v);
}

static inline fv_real fv_i_sc_abs(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

/* Per-component clamp: tests the lower bound before the upper one, matching the expression written in the header. */
static inline fv_real fv_i_sc_clamp1(fv_real x, fv_real lo, fv_real hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

/* ---------------- 2D ---------------- */

static inline fv_vec2 fv_vec2_make(fv_real x, fv_real y)
{
    fv_vec2 v;
    v.x = x;
    v.y = y;
    return v;
}

static inline fv_vec2 fv_vec2_add(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x + b.x, a.y + b.y);
}

static inline fv_vec2 fv_vec2_sub(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x - b.x, a.y - b.y);
}

static inline fv_vec2 fv_vec2_neg(fv_vec2 v)
{
    return fv_vec2_make(-v.x, -v.y);
}

static inline fv_vec2 fv_vec2_scale(fv_vec2 v, fv_real s)
{
    return fv_vec2_make(v.x * s, v.y * s);
}

static inline fv_vec2 fv_vec2_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    return fv_vec2_make(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
}

static inline fv_vec2 fv_vec2_min(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y);
}

static inline fv_vec2 fv_vec2_max(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y);
}

static inline fv_vec2 fv_vec2_clamp(fv_vec2 v, fv_vec2 lo, fv_vec2 hi)
{
    return fv_vec2_make(fv_i_sc_clamp1(v.x, lo.x, hi.x),
                        fv_i_sc_clamp1(v.y, lo.y, hi.y));
}

static inline fv_vec2 fv_vec2_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    return fv_vec2_make(a.x + b.x * s,
                        a.y + b.y * s);
}

static inline fv_real fv_vec2_dot(fv_vec2 a, fv_vec2 b)
{
    return a.x * b.x + a.y * b.y;
}

static inline fv_real fv_vec2_len_sqr(fv_vec2 v)
{
    return fv_vec2_dot(v, v);
}

static inline fv_real fv_vec2_len(fv_vec2 v)
{
    return fv_i_sc_sqrt(fv_vec2_len_sqr(v));
}

static inline fv_real fv_vec2_dist_sqr(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_len_sqr(fv_vec2_sub(a, b));
}

static inline fv_real fv_vec2_dist(fv_vec2 a, fv_vec2 b)
{
    return fv_i_sc_sqrt(fv_vec2_dist_sqr(a, b));
}

static inline fv_vec2 fv_vec2_normalize(fv_vec2 v)
{
    fv_real len = fv_vec2_len(v);
    if (len <= FV_EPS) {
        return fv_vec2_make((fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec2_scale(v, inv);
    }
}

static inline fv_vec2 fv_vec2_project(fv_vec2 v, fv_vec2 onto)
{
    fv_real denom = fv_vec2_len_sqr(onto);
    /* denom is a squared length, so the squared threshold FV_EPS_SQR must be used (P0 fix: this used to be the linear threshold FV_EPS). */
    if (fv_i_sc_abs(denom) <= FV_EPS_SQR) {
        return fv_vec2_make((fv_real)0, (fv_real)0);
    }
    return fv_vec2_scale(onto, fv_vec2_dot(v, onto) / denom);
}

static inline int fv_vec2_near(fv_vec2 a, fv_vec2 b, fv_real eps)
{
    if (fv_i_sc_abs(a.x - b.x) > eps) {
        return 0;
    }
    if (fv_i_sc_abs(a.y - b.y) > eps) {
        return 0;
    }
    return 1;
}

/* ---------------- 3D ---------------- */

static inline fv_vec3 fv_vec3_make(fv_real x, fv_real y, fv_real z)
{
    fv_vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static inline fv_vec3 fv_vec3_add(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline fv_vec3 fv_vec3_sub(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline fv_vec3 fv_vec3_neg(fv_vec3 v)
{
    return fv_vec3_make(-v.x, -v.y, -v.z);
}

static inline fv_vec3 fv_vec3_scale(fv_vec3 v, fv_real s)
{
    return fv_vec3_make(v.x * s, v.y * s, v.z * s);
}

static inline fv_vec3 fv_vec3_lerp(fv_vec3 a, fv_vec3 b, fv_real t)
{
    return fv_vec3_make(a.x + (b.x - a.x) * t,
                        a.y + (b.y - a.y) * t,
                        a.z + (b.z - a.z) * t);
}

static inline fv_vec3 fv_vec3_min(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y,
                        a.z < b.z ? a.z : b.z);
}

static inline fv_vec3 fv_vec3_max(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y,
                        a.z > b.z ? a.z : b.z);
}

static inline fv_vec3 fv_vec3_clamp(fv_vec3 v, fv_vec3 lo, fv_vec3 hi)
{
    return fv_vec3_make(fv_i_sc_clamp1(v.x, lo.x, hi.x),
                        fv_i_sc_clamp1(v.y, lo.y, hi.y),
                        fv_i_sc_clamp1(v.z, lo.z, hi.z));
}

static inline fv_vec3 fv_vec3_madd(fv_vec3 a, fv_vec3 b, fv_real s)
{
    return fv_vec3_make(a.x + b.x * s,
                        a.y + b.y * s,
                        a.z + b.z * s);
}

static inline fv_real fv_vec3_dot(fv_vec3 a, fv_vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline fv_real fv_vec3_len_sqr(fv_vec3 v)
{
    return fv_vec3_dot(v, v);
}

static inline fv_real fv_vec3_len(fv_vec3 v)
{
    return fv_i_sc_sqrt(fv_vec3_len_sqr(v));
}

static inline fv_real fv_vec3_dist_sqr(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_len_sqr(fv_vec3_sub(a, b));
}

static inline fv_real fv_vec3_dist(fv_vec3 a, fv_vec3 b)
{
    return fv_i_sc_sqrt(fv_vec3_dist_sqr(a, b));
}

static inline fv_vec3 fv_vec3_normalize(fv_vec3 v)
{
    fv_real len = fv_vec3_len(v);
    if (len <= FV_EPS) {
        return fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec3_scale(v, inv);
    }
}

static inline fv_vec3 fv_vec3_project(fv_vec3 v, fv_vec3 onto)
{
    fv_real denom = fv_vec3_len_sqr(onto);
    /* As in fv_vec2_project: a squared quantity takes the squared threshold. */
    if (fv_i_sc_abs(denom) <= FV_EPS_SQR) {
        return fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)0);
    }
    return fv_vec3_scale(onto, fv_vec3_dot(v, onto) / denom);
}

static inline fv_vec3 fv_vec3_cross(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x);
}

static inline int fv_vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps)
{
    if (fv_i_sc_abs(a.x - b.x) > eps) {
        return 0;
    }
    if (fv_i_sc_abs(a.y - b.y) > eps) {
        return 0;
    }
    if (fv_i_sc_abs(a.z - b.z) > eps) {
        return 0;
    }
    return 1;
}

/* ---------------- 2D bounding box (for broad-phase) ---------------- */

static inline fv_aabb2 fv_aabb2_merge(fv_aabb2 box, fv_vec2 p)
{
    fv_aabb2 r;
    r.lo = fv_vec2_min(box.lo, p);
    r.hi = fv_vec2_max(box.hi, p);
    return r;
}

static inline int fv_aabb2_overlaps(fv_aabb2 a, fv_aabb2 b)
{
    if (a.lo.x > b.hi.x || b.lo.x > a.hi.x) {
        return 0;
    }
    if (a.lo.y > b.hi.y || b.lo.y > a.hi.y) {
        return 0;
    }
    return 1;
}


static inline fv_real fv_i_sc_sqrt4(fv_real v)
{
    return sqrtf(v);
}

static inline fv_real fv_i_sc_abs4(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

/* Per-component clamp: tests the lower bound before the upper one (the same expression as fv_i_sc_clamp1 in src/flashvec.cpp,
 * and with lo > hi it places the value by the order of those same two ifs; a NaN compares false against both bounds and is returned unchanged). */
static inline fv_real fv_i_sc_clamp1_4(fv_real x, fv_real lo, fv_real hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

static inline fv_vec4 fv_vec4_make(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_vec4 v;
    v.x = x;
    v.y = y;
    v.z = z;
    v.w = w;
    return v;
}

static inline fv_vec4 fv_vec4_add(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static inline fv_vec4 fv_vec4_sub(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static inline fv_vec4 fv_vec4_neg(fv_vec4 v)
{
    return fv_vec4_make(-v.x, -v.y, -v.z, -v.w);
}

static inline fv_vec4 fv_vec4_scale(fv_vec4 v, fv_real s)
{
    return fv_vec4_make(v.x * s, v.y * s, v.z * s, v.w * s);
}

static inline fv_vec4 fv_vec4_lerp(fv_vec4 a, fv_vec4 b, fv_real t)
{
    return fv_vec4_make(a.x + (b.x - a.x) * t,
                        a.y + (b.y - a.y) * t,
                        a.z + (b.z - a.z) * t,
                        a.w + (b.w - a.w) * t);
}

/* The NaN rule is exactly that of fv_vec2_min / fv_vec3_min: (a < b) ? a : b, not IEEE minNum. */
static inline fv_vec4 fv_vec4_min(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y,
                        a.z < b.z ? a.z : b.z,
                        a.w < b.w ? a.w : b.w);
}

static inline fv_vec4 fv_vec4_max(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y,
                        a.z > b.z ? a.z : b.z,
                        a.w > b.w ? a.w : b.w);
}

static inline fv_vec4 fv_vec4_clamp(fv_vec4 v, fv_vec4 lo, fv_vec4 hi)
{
    return fv_vec4_make(fv_i_sc_clamp1_4(v.x, lo.x, hi.x),
                        fv_i_sc_clamp1_4(v.y, lo.y, hi.y),
                        fv_i_sc_clamp1_4(v.z, lo.z, hi.z),
                        fv_i_sc_clamp1_4(v.w, lo.w, hi.w));
}

static inline fv_vec4 fv_vec4_madd(fv_vec4 a, fv_vec4 b, fv_real s)
{
    return fv_vec4_make(a.x + b.x * s,
                        a.y + b.y * s,
                        a.z + b.z * s,
                        a.w + b.w * s);
}

/* Left to right: ((x*x' + y*y') + z*z') + w*w' -- the SoA kernels of the batch entry points must reproduce this association order,
 * otherwise the bit-level comparison fails (see the bit-exactness promise boundary in include/flashvec.h). */
static inline fv_real fv_vec4_dot(fv_vec4 a, fv_vec4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static inline fv_real fv_vec4_len_sqr(fv_vec4 v)
{
    return fv_vec4_dot(v, v);
}

static inline fv_real fv_vec4_len(fv_vec4 v)
{
    return fv_i_sc_sqrt4(fv_vec4_len_sqr(v));
}

static inline fv_real fv_vec4_dist_sqr(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_len_sqr(fv_vec4_sub(a, b));
}

static inline fv_real fv_vec4_dist(fv_vec4 a, fv_vec4 b)
{
    return fv_i_sc_sqrt4(fv_vec4_dist_sqr(a, b));
}

static inline fv_vec4 fv_vec4_normalize(fv_vec4 v)
{
    fv_real len = fv_vec4_len(v);
    if (len <= FV_EPS) {
        return fv_vec4_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec4_scale(v, inv);
    }
}

static inline fv_vec4 fv_vec4_project(fv_vec4 v, fv_vec4 onto)
{
    fv_real denom = fv_vec4_len_sqr(onto);
    /* denom is a squared quantity, so the squared threshold FV_EPS_SQR must be used (the same fix as in fv_vec2_project). */
    if (fv_i_sc_abs4(denom) <= FV_EPS_SQR) {
        return fv_vec4_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
    }
    return fv_vec4_scale(onto, fv_vec4_dot(v, onto) / denom);
}

static inline int fv_vec4_near(fv_vec4 a, fv_vec4 b, fv_real eps)
{
    if (fv_i_sc_abs4(a.x - b.x) > eps) {
        return 0;
    }
    if (fv_i_sc_abs4(a.y - b.y) > eps) {
        return 0;
    }
    if (fv_i_sc_abs4(a.z - b.z) > eps) {
        return 0;
    }
    if (fv_i_sc_abs4(a.w - b.w) > eps) {
        return 0;
    }
    return 1;
}


static inline fv_quat fv_quat_make(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_quat q;
    q.x = x;
    q.y = y;
    q.z = z;
    q.w = w;
    return q;
}

static inline fv_quat fv_quat_add(fv_quat a, fv_quat b)
{
    return fv_quat_make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static inline fv_quat fv_quat_sub(fv_quat a, fv_quat b)
{
    return fv_quat_make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static inline fv_quat fv_quat_neg(fv_quat q)
{
    return fv_quat_make(-q.x, -q.y, -q.z, -q.w);
}

static inline fv_quat fv_quat_scale(fv_quat q, fv_real s)
{
    return fv_quat_make(q.x * s, q.y * s, q.z * s, q.w * s);
}

/* Left to right: ((x*x' + y*y') + z*z') + w*w'. The same association order the batch normalize uses. */
static inline fv_real fv_quat_dot(fv_quat a, fv_quat b)
{
    return ((a.x * b.x + a.y * b.y) + a.z * b.z) + a.w * b.w;
}

static inline fv_real fv_quat_len_sqr(fv_quat q)
{
    return fv_quat_dot(q, q);
}

static inline fv_real fv_quat_len(fv_quat q)
{
    return sqrtf(fv_quat_len_sqr(q));
}

/* A zero quaternion is not a rotation; the identity keeps the result usable and NaN-free. This is the
 * deliberate difference from fv_vec*_normalize (which returns the zero vector) and is stated in the
 * header contract. */
static inline fv_quat fv_quat_normalize(fv_quat q)
{
    fv_real len = fv_quat_len(q);
    if (len <= FV_EPS) {
        return fv_quat_identity();
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_quat_scale(q, inv);
    }
}

static inline fv_quat fv_quat_conjugate(fv_quat q)
{
    return fv_quat_make(-q.x, -q.y, -q.z, q.w);
}

/* General inverse: conjugate(q) / |q|^2. The threshold is on the squared length, so it is
 * dimensionally consistent with FV_EPS (length <= FV_EPS <=> len_sqr <= FV_EPS_SQR). A degenerate
 * input writes the zero quaternion and returns 0; the value form returns that zero. */
static inline int fv_quat_inverse_checked(fv_quat q, fv_quat *out)
{
    fv_real l2 = fv_quat_len_sqr(q);
    if (!(l2 > FV_EPS_SQR)) {
        *out = fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / l2;
        fv_quat c = fv_quat_conjugate(q);
        *out = fv_quat_scale(c, inv);
        return 1;
    }
}

static inline fv_quat fv_quat_inverse(fv_quat q)
{
    fv_quat r;
    (void)fv_quat_inverse_checked(q, &r);
    return r;
}

/* Hamilton product a*b: applies b first, then a. Association order fixed as
 * ((w1*x2 + x1*w2) + y1*z2) - z1*y2, and the three analogues -- see the note in include/flashvec.h. */
static inline fv_quat fv_quat_mul(fv_quat a, fv_quat b)
{
    fv_quat r;
    r.x = ((a.w * b.x + a.x * b.w) + a.y * b.z) - a.z * b.y;
    r.y = ((a.w * b.y - a.x * b.z) + a.y * b.w) + a.z * b.x;
    r.z = ((a.w * b.z + a.x * b.y) - a.y * b.x) + a.z * b.w;
    r.w = ((a.w * b.w - a.x * b.x) - a.y * b.y) - a.z * b.z;
    return r;
}

/* v' = q*(v,0)*conjugate(q), evaluated as t = 2*cross(q.xyz,v); d = cross(q.xyz,t); out = (v + q.w*t) + d.
 * q must be unit; this function does not normalize (same precondition as a rotation matrix). */
static inline fv_vec3 fv_quat_rotate_vec3(fv_quat q, fv_vec3 v)
{
    fv_real tx = (fv_real)2 * (q.y * v.z - q.z * v.y);
    fv_real ty = (fv_real)2 * (q.z * v.x - q.x * v.z);
    fv_real tz = (fv_real)2 * (q.x * v.y - q.y * v.x);
    fv_real dx = q.y * tz - q.z * ty;
    fv_real dy = q.z * tx - q.x * tz;
    fv_real dz = q.x * ty - q.y * tx;
    fv_vec3 r;
    r.x = (v.x + q.w * tx) + dx;
    r.y = (v.y + q.w * ty) + dy;
    r.z = (v.z + q.w * tz) + dz;
    return r;
}

/* Rotation matrix -> unit quaternion, Shepperd's method: pick the largest of the four diagonal
 * expressions so the divisor is the largest one (the naive trace-only form loses precision when
 * trace is near -1). The matrix is assumed to be a rotation (orthonormal, det +1); the result is
 * normalized. For a non-rotation this still returns a unit quaternion, but not the matrix's rotation,
 * and the by-value API has no error channel to report that. */
static inline fv_quat fv_quat_from_mat3(const fv_mat3 *m)
{
    fv_real trace = m->m[0] + m->m[4] + m->m[8];
    fv_quat q;
    if (trace > (fv_real)0) {
        fv_real s = sqrtf(trace + (fv_real)1) * (fv_real)2;
        q.w = (fv_real)0.25 * s;
        q.x = (m->m[7] - m->m[5]) / s;
        q.y = (m->m[2] - m->m[6]) / s;
        q.z = (m->m[3] - m->m[1]) / s;
    } else if (m->m[0] > m->m[4] && m->m[0] > m->m[8]) {
        fv_real s = sqrtf((fv_real)1 + m->m[0] - m->m[4] - m->m[8]) * (fv_real)2;
        q.w = (m->m[7] - m->m[5]) / s;
        q.x = (fv_real)0.25 * s;
        q.y = (m->m[1] + m->m[3]) / s;
        q.z = (m->m[2] + m->m[6]) / s;
    } else if (m->m[4] > m->m[8]) {
        fv_real s = sqrtf((fv_real)1 + m->m[4] - m->m[0] - m->m[8]) * (fv_real)2;
        q.w = (m->m[2] - m->m[6]) / s;
        q.x = (m->m[1] + m->m[3]) / s;
        q.y = (fv_real)0.25 * s;
        q.z = (m->m[5] + m->m[7]) / s;
    } else {
        fv_real s = sqrtf((fv_real)1 + m->m[8] - m->m[0] - m->m[4]) * (fv_real)2;
        q.w = (m->m[3] - m->m[1]) / s;
        q.x = (m->m[2] + m->m[6]) / s;
        q.y = (m->m[5] + m->m[7]) / s;
        q.z = (fv_real)0.25 * s;
    }
    return fv_quat_normalize(q);
}

/* Shortest-arc normalized lerp: negating the second quaternion when dot < 0 selects the same rotation
 * with the closer sign, then the components are interpolated and normalized. */
static inline fv_quat fv_quat_nlerp(fv_quat a, fv_quat b, fv_real t)
{
    if (fv_quat_dot(a, b) < (fv_real)0) {
        b = fv_quat_neg(b);
    }
    {
        fv_quat r;
        r.x = a.x + (b.x - a.x) * t;
        r.y = a.y + (b.y - a.y) * t;
        r.z = a.z + (b.z - a.z) * t;
        r.w = a.w + (b.w - a.w) * t;
        return fv_quat_normalize(r);
    }
}

/* Shortest-arc slerp. The small-angle threshold (dot > 0.9995) keeps sin(theta) away from underflow and
 * hands the near-identical case to nlerp, which is stable there. */
static inline fv_quat fv_quat_slerp(fv_quat a, fv_quat b, fv_real t)
{
    fv_real d = fv_quat_dot(a, b);
    if (d < (fv_real)0) {
        b = fv_quat_neg(b);
        d = -d;
    }
    if (d > (fv_real)0.9995 || !(d < (fv_real)1)) {
        return fv_quat_nlerp(a, b, t);
    }
    {
        fv_real theta = acosf(d);
        fv_real st = sinf(theta);
        if (st == (fv_real)0) {
            return fv_quat_nlerp(a, b, t);
        }
        {
            fv_real w1 = sinf((fv_real)(1 - t) * theta) / st;
            fv_real w2 = sinf(t * theta) / st;
            fv_quat r = fv_quat_add(fv_quat_scale(a, w1), fv_quat_scale(b, w2));
            return fv_quat_normalize(r);
        }
    }
}

/* Axis-angle: normalize the axis first (a zero axis has no direction and returns the identity), then
 * the half-angle sine/cosine. */
static inline fv_quat fv_quat_from_axis_angle(fv_vec3 axis, fv_real angle)
{
    fv_real len = fv_vec3_len(axis);
    if (len <= FV_EPS) {
        return fv_quat_identity();
    }
    {
        fv_real inv = (fv_real)1 / len;
        fv_real half = angle * (fv_real)0.5;
        fv_real s = sinf(half);
        fv_real c = cosf(half);
        return fv_quat_make(axis.x * inv * s, axis.y * inv * s, axis.z * inv * s, c);
    }
}

/* Euler XYZ, applied X first, then Y, then Z: q = qz * qy * qx. */
static inline fv_quat fv_quat_from_euler_xyz(fv_real x, fv_real y, fv_real z)
{
    fv_quat qx = fv_quat_from_axis_angle(fv_vec3_make((fv_real)1, (fv_real)0, (fv_real)0), x);
    fv_quat qy = fv_quat_from_axis_angle(fv_vec3_make((fv_real)0, (fv_real)1, (fv_real)0), y);
    fv_quat qz = fv_quat_from_axis_angle(fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)1), z);
    return fv_quat_mul(qz, fv_quat_mul(qy, qx));
}


static inline fv_mat3 fv_mat3_from_rows(fv_vec3 r0, fv_vec3 r1, fv_vec3 r2)
{
    fv_mat3 m;
    m.m[0] = r0.x;
    m.m[1] = r0.y;
    m.m[2] = r0.z;
    m.m[3] = r1.x;
    m.m[4] = r1.y;
    m.m[5] = r1.z;
    m.m[6] = r2.x;
    m.m[7] = r2.y;
    m.m[8] = r2.z;
    return m;
}

static inline fv_vec3 fv_mat3_transform_vec3(const fv_mat3 *m, fv_vec3 v)
{
    fv_vec3 r;
    r.x = (m->m[0] * v.x + m->m[1] * v.y) + m->m[2] * v.z;
    r.y = (m->m[3] * v.x + m->m[4] * v.y) + m->m[5] * v.z;
    r.z = (m->m[6] * v.x + m->m[7] * v.y) + m->m[8] * v.z;
    return r;
}

static inline fv_mat3 fv_mat3_mul(const fv_mat3 *a, const fv_mat3 *b)
{
    fv_mat3 r;
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r.m[i * 3 + j] = (a->m[i * 3 + 0] * b->m[0 * 3 + j] +
                              a->m[i * 3 + 1] * b->m[1 * 3 + j]) +
                             a->m[i * 3 + 2] * b->m[2 * 3 + j];
        }
    }
    return r;
}

static inline fv_mat3 fv_mat3_transpose(const fv_mat3 *m)
{
    fv_mat3 r;
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r.m[i * 3 + j] = m->m[j * 3 + i];
        }
    }
    return r;
}

/* Row-major rotation matrix of a unit quaternion:
 *   R = [[1-2(yy+zz), 2(xy-wz),   2(xz+wy)],
 *        [2(xy+wz),   1-2(xx+zz), 2(yz-wx)],
 *        [2(xz-wy),   2(yz+wx),   1-2(xx+yy)]]
 * The product order inside each entry is fixed (2*(...), then the subtraction from 1) so the function
 * has one documented bit pattern; it is the matrix form of the same rotation as fv_quat_rotate_vec3,
 * but the two are different expressions and are compared approximately, not with memcmp. */
static inline fv_mat3 fv_mat3_from_quat(fv_quat q)
{
    fv_real xx = q.x * q.x;
    fv_real yy = q.y * q.y;
    fv_real zz = q.z * q.z;
    fv_real xy = q.x * q.y;
    fv_real xz = q.x * q.z;
    fv_real yz = q.y * q.z;
    fv_real wx = q.w * q.x;
    fv_real wy = q.w * q.y;
    fv_real wz = q.w * q.z;
    fv_mat3 r;
    r.m[0] = (fv_real)1 - (fv_real)2 * (yy + zz);
    r.m[1] = (fv_real)2 * (xy - wz);
    r.m[2] = (fv_real)2 * (xz + wy);
    r.m[3] = (fv_real)2 * (xy + wz);
    r.m[4] = (fv_real)1 - (fv_real)2 * (xx + zz);
    r.m[5] = (fv_real)2 * (yz - wx);
    r.m[6] = (fv_real)2 * (xz - wy);
    r.m[7] = (fv_real)2 * (yz + wx);
    r.m[8] = (fv_real)1 - (fv_real)2 * (xx + yy);
    return r;
}

static inline fv_real fv_i_sc_abs3(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

/* |det| <= FV_EPS * max|entry|^3 is the singularity test: a cubic quantity against a cubic
 * threshold, relative to the matrix scale, so 1e-3 * I stays invertible while a rank-deficient
 * matrix (det == 0) and a uniformly near-singular one both fall out. !(scale > 0) also catches a
 * NaN entry; a non-finite matrix is treated as singular rather than propagating NaNs into the
 * cofactors. */
static inline fv_real fv_mat3_determinant(const fv_mat3 *m)
{
    return (m->m[0] * (m->m[4] * m->m[8] - m->m[5] * m->m[7]) -
            m->m[1] * (m->m[3] * m->m[8] - m->m[5] * m->m[6])) +
           m->m[2] * (m->m[3] * m->m[7] - m->m[4] * m->m[6]);
}

static inline void fv_i_sc_zero3(fv_mat3 *m)
{
    int k;
    for (k = 0; k < 9; ++k) {
        m->m[k] = (fv_real)0;
    }
}

static inline int fv_i_sc_inv3_impl(const fv_mat3 *m, fv_mat3 *out)
{
    fv_real det = fv_mat3_determinant(m);
    fv_real scale = (fv_real)0;
    int k;
    for (k = 0; k < 9; ++k) {
        fv_real a = fv_i_sc_abs3(m->m[k]);
        if (a > scale) {
            scale = a;
        }
    }
    if (!(scale > (fv_real)0) ||
        !(fv_i_sc_abs3(det) > FV_EPS * ((scale * scale) * scale))) {
        fv_i_sc_zero3(out);
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / det;
        fv_mat3 r;
        r.m[0] = (m->m[4] * m->m[8] - m->m[5] * m->m[7]) * inv;
        r.m[1] = (m->m[2] * m->m[7] - m->m[1] * m->m[8]) * inv;
        r.m[2] = (m->m[1] * m->m[5] - m->m[2] * m->m[4]) * inv;
        r.m[3] = (m->m[5] * m->m[6] - m->m[3] * m->m[8]) * inv;
        r.m[4] = (m->m[0] * m->m[8] - m->m[2] * m->m[6]) * inv;
        r.m[5] = (m->m[2] * m->m[3] - m->m[0] * m->m[5]) * inv;
        r.m[6] = (m->m[3] * m->m[7] - m->m[4] * m->m[6]) * inv;
        r.m[7] = (m->m[1] * m->m[6] - m->m[0] * m->m[7]) * inv;
        r.m[8] = (m->m[0] * m->m[4] - m->m[1] * m->m[3]) * inv;
        *out = r;
        return 1;
    }
}

static inline fv_mat3 fv_mat3_inverse(const fv_mat3 *m)
{
    fv_mat3 r;
    (void)fv_i_sc_inv3_impl(m, &r);
    return r;
}

static inline int fv_mat3_inverse_checked(const fv_mat3 *m, fv_mat3 *out)
{
    return fv_i_sc_inv3_impl(m, out);
}


static inline fv_real fv_i_sc_abs4m(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

static inline void fv_i_sc_zero4m(fv_mat4 *m)
{
    int k;
    for (k = 0; k < 16; ++k) {
        m->m[k] = (fv_real)0;
    }
}

static inline fv_mat4 fv_mat4_mul(const fv_mat4 *a, const fv_mat4 *b)
{
    fv_mat4 r;
    int i;
    int j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            r.m[i * 4 + j] = ((a->m[i * 4 + 0] * b->m[0 * 4 + j] +
                               a->m[i * 4 + 1] * b->m[1 * 4 + j]) +
                              a->m[i * 4 + 2] * b->m[2 * 4 + j]) +
                             a->m[i * 4 + 3] * b->m[3 * 4 + j];
        }
    }
    return r;
}

static inline int fv_i_sc_inv4m_impl(const fv_mat4 *m, fv_mat4 *out)
{
    const fv_real *v = m->m;
    fv_real s0 = v[0] * v[5] - v[4] * v[1];
    fv_real s1 = v[0] * v[6] - v[4] * v[2];
    fv_real s2 = v[0] * v[7] - v[4] * v[3];
    fv_real s3 = v[1] * v[6] - v[5] * v[2];
    fv_real s4 = v[1] * v[7] - v[5] * v[3];
    fv_real s5 = v[2] * v[7] - v[6] * v[3];
    fv_real c5 = v[10] * v[15] - v[14] * v[11];
    fv_real c4 = v[9] * v[15] - v[13] * v[11];
    fv_real c3 = v[9] * v[14] - v[13] * v[10];
    fv_real c2 = v[8] * v[15] - v[12] * v[11];
    fv_real c1 = v[8] * v[14] - v[12] * v[10];
    fv_real c0 = v[8] * v[13] - v[12] * v[9];
    fv_real det = (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2) - s4 * c1 + s5 * c0;
    fv_real scale = (fv_real)0;
    fv_real scale2;
    int k;
    for (k = 0; k < 16; ++k) {
        fv_real a = fv_i_sc_abs4m(v[k]);
        if (a > scale) {
            scale = a;
        }
    }
    scale2 = scale * scale;
    if (!(scale > (fv_real)0) ||
        !(fv_i_sc_abs4m(det) > FV_EPS * (scale2 * scale2))) {
        fv_i_sc_zero4m(out);
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / det;
        fv_mat4 r;
        r.m[0] = (v[5] * c5 - v[6] * c4 + v[7] * c3) * inv;
        r.m[1] = (-v[1] * c5 + v[2] * c4 - v[3] * c3) * inv;
        r.m[2] = (v[13] * s5 - v[14] * s4 + v[15] * s3) * inv;
        r.m[3] = (-v[9] * s5 + v[10] * s4 - v[11] * s3) * inv;
        r.m[4] = (-v[4] * c5 + v[6] * c2 - v[7] * c1) * inv;
        r.m[5] = (v[0] * c5 - v[2] * c2 + v[3] * c1) * inv;
        r.m[6] = (-v[12] * s5 + v[14] * s2 - v[15] * s1) * inv;
        r.m[7] = (v[8] * s5 - v[10] * s2 + v[11] * s1) * inv;
        r.m[8] = (v[4] * c4 - v[5] * c2 + v[7] * c0) * inv;
        r.m[9] = (-v[0] * c4 + v[1] * c2 - v[3] * c0) * inv;
        r.m[10] = (v[12] * s4 - v[13] * s2 + v[15] * s0) * inv;
        r.m[11] = (-v[8] * s4 + v[9] * s2 - v[11] * s0) * inv;
        r.m[12] = (-v[4] * c3 + v[5] * c1 - v[6] * c0) * inv;
        r.m[13] = (v[0] * c3 - v[1] * c1 + v[2] * c0) * inv;
        r.m[14] = (-v[12] * s3 + v[13] * s1 - v[14] * s0) * inv;
        r.m[15] = (v[8] * s3 - v[9] * s1 + v[10] * s0) * inv;
        *out = r;
        return 1;
    }
}

static inline fv_mat4 fv_mat4_inverse(const fv_mat4 *m)
{
    fv_mat4 r;
    (void)fv_i_sc_inv4m_impl(m, &r);
    return r;
}

static inline int fv_mat4_inverse_checked(const fv_mat4 *m, fv_mat4 *out)
{
    return fv_i_sc_inv4m_impl(m, out);
}


static inline fv_real fv_atan_approx(fv_real x)
{
    if (!(x >= (fv_real)-1 && x <= (fv_real)1)) {
        return (fv_real)NAN;
    }
    {
        fv_real x2 = x * x;
        fv_real p = (fv_real)0.0208351;
        p = p * x2 + (fv_real)-0.0851330;
        p = p * x2 + (fv_real)0.1801410;
        p = p * x2 + (fv_real)-0.3302995;
        p = p * x2 + (fv_real)0.9998660;
        return x * p;
    }
}


static inline fv_vec2 fv_vec2_addmul(fv_vec2 x, fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    return r;
}

static inline fv_vec3 fv_vec3_addmul(fv_vec3 x, fv_vec3 v, fv_real s)
{
    fv_vec3 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    r.z = x.z + v.z * s;
    return r;
}

static inline fv_vec4 fv_vec4_addmul(fv_vec4 x, fv_vec4 v, fv_real s)
{
    fv_vec4 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    r.z = x.z + v.z * s;
    r.w = x.w + v.w * s;
    return r;
}

/* a*b + c, product first, then the sum -- the third operand is a vector here, which is the whole point (fv_vec2_madd's
 * third operand is a scalar). */
static inline fv_vec2 fv_vec2_maddv(fv_vec2 a, fv_vec2 b, fv_vec2 c)
{
    fv_vec2 r;
    r.x = a.x * b.x + c.x;
    r.y = a.y * b.y + c.y;
    return r;
}

#endif /* !FV_BUILD: scalar API inline definitions */

#endif


static inline fv_vec2 fv_i_s_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    fv_vec2 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    return r;
}


static inline fv_vec2 fv_i_s_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    fv_vec2 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    return r;
}


/* The same expression as fv_vec2_len: sqrtf(x*x + y*y). */
static inline fv_real fv_i_s_len(fv_real x, fv_real y)
{
    return sqrtf(x * x + y * y);
}


/* The same semantics, expressions and branches as fv_vec2_normalize. */
static inline void fv_i_s_normalize(fv_real x, fv_real y, fv_real *ox, fv_real *oy)
{
    fv_real len = fv_i_s_len(x, y);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
    }
}


/* The scalar form of the fused pipeline: t = a * s, u = t + b, then normalize -- exactly the expressions of
 * fv_i_s_scale -> fv_i_s_add -> fv_i_s_normalize in that order, written once. A float store/reload (what the
 * three-call form does with its intermediate buffer) is exact, so keeping u in a register changes no bit.
 * The parentheses are written out to pin the order down in the source; fp-contract=off keeps a + b*c from becoming an FMA. */
static inline void fv_i_s_madd_norm(fv_real ax, fv_real ay, fv_real bx, fv_real by, fv_real s,
                                    fv_real *ox, fv_real *oy)
{
    fv_real ux = (ax * s) + bx;
    fv_real uy = (ay * s) + by;
    fv_i_s_normalize(ux, uy, ox, oy);
}

static inline fv_vec4 fv_i_s_lerp4(fv_vec4 a, fv_vec4 b, fv_real t)
{
    fv_vec4 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    r.z = a.z + (b.z - a.z) * t;
    r.w = a.w + (b.w - a.w) * t;
    return r;
}


static inline fv_vec4 fv_i_s_madd4(fv_vec4 a, fv_vec4 b, fv_real s)
{
    fv_vec4 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    r.z = a.z + b.z * s;
    r.w = a.w + b.w * s;
    return r;
}


/* 3D length: sqrtf((x*x + y*y) + z*z) -- the same left-to-right order as fv_vec3_len_sqr (= fv_vec3_dot(v,v)).
 * The parentheses are written out in full to pin the association order down in the source, not for the compiler's benefit. */
static inline fv_real fv_i_s_len3(fv_real x, fv_real y, fv_real z)
{
    return sqrtf((x * x + y * y) + z * z);
}


static inline void fv_i_s_normalize3(fv_real x, fv_real y, fv_real z,
                                     fv_real *ox, fv_real *oy, fv_real *oz)
{
    fv_real len = fv_i_s_len3(x, y, z);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        *oz = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
        *oz = z * inv;
    }
}


/* 4D length: sqrtf(((x*x + y*y) + z*z) + w*w) -- the same left-to-right order as fv_vec4_len_sqr. */
static inline fv_real fv_i_s_len4(fv_real x, fv_real y, fv_real z, fv_real w)
{
    return sqrtf(((x * x + y * y) + z * z) + w * w);
}


static inline void fv_i_s_normalize4(fv_real x, fv_real y, fv_real z, fv_real w,
                                     fv_real *ox, fv_real *oy, fv_real *oz, fv_real *ow)
{
    fv_real len = fv_i_s_len4(x, y, z, w);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        *oz = (fv_real)0;
        *ow = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
        *oz = z * inv;
        *ow = w * inv;
    }
}


/* Matrix transform, one element: the reference expression the AVX2 kernel must reproduce.
 * ((m0*x + m1*y) + m2*z) + m3 -- C left-associativity already produces this order; the parens
 * are written out so that reordering it is a visible edit instead of an accident.
 * dirs drops the fourth term instead of adding m3 * 0: x + (+0.0) would turn -0.0 into +0.0.
 * w is a constant, not a dot product -- see the fv_mat4 typedef for why. */
static inline fv_vec4 fv_i_s_mat4_point(const fv_mat4 *__restrict m, fv_vec4 p)
{
    const fv_real x = p.x;
    const fv_real y = p.y;
    const fv_real z = p.z;
    fv_vec4 r;
    r.x = ((m->m[0] * x + m->m[1] * y) + m->m[2] * z) + m->m[3];
    r.y = ((m->m[4] * x + m->m[5] * y) + m->m[6] * z) + m->m[7];
    r.z = ((m->m[8] * x + m->m[9] * y) + m->m[10] * z) + m->m[11];
    r.w = (fv_real)1;
    return r;
}


static inline fv_vec4 fv_i_s_mat4_dir(const fv_mat4 *__restrict m, fv_vec4 p)
{
    const fv_real x = p.x;
    const fv_real y = p.y;
    const fv_real z = p.z;
    fv_vec4 r;
    r.x = (m->m[0] * x + m->m[1] * y) + m->m[2] * z;
    r.y = (m->m[4] * x + m->m[5] * y) + m->m[6] * z;
    r.z = (m->m[8] * x + m->m[9] * y) + m->m[10] * z;
    r.w = (fv_real)0;
    return r;
}


#if FV_I_INTRINSICS

/* lerp = a + (b - a) * t, with the scalar multiply-then-add order; fp-contract=off guarantees it never becomes an FMA. */
FV_I_TARGET_AVX2 static inline void fv_i_lerp4_avx2(const fv_vec2 *__restrict a,
                                                    const fv_vec2 *__restrict b, fv_real t,
                                                    fv_vec2 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
}


/* madd = a + b * s, as above. */
FV_I_TARGET_AVX2 static inline void fv_i_madd4_avx2(const fv_vec2 *__restrict a,
                                                    const fv_vec2 *__restrict b, fv_real s,
                                                    fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
}


/* Same two kernels with the output written by NT stores (vec2 _x4 shape, 8 B stride: at most 3 peel steps).
 * These sit **inside the pragma**: lerp is a + (b - a) * t and madd is a + b * s, so the mul feeds straight into the
 * add and the consumer's -ffp-contract=fast would contract them into vfmadd -- which changes the bits, and the whole
 * point of the _stream family is that only the store instruction differs. */
FV_I_TARGET_AVX2 static inline void fv_i_lerp4_stream_avx2(const fv_vec2 *__restrict a,
                                                           const fv_vec2 *__restrict b, fv_real t,
                                                           fv_vec2 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_lerp(a[i], b[i], t);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_lerp4_avx2(a, b, t, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp(a[i], b[i], t);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_madd4_stream_avx2(const fv_vec2 *__restrict a,
                                                           const fv_vec2 *__restrict b, fv_real s,
                                                           fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_madd(a[i], b[i], s);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_madd4_avx2(a, b, s, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd(a[i], b[i], s);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_len_soa_avx2(const float *__restrict ax,
                                                      const float *__restrict ay,
                                                      float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(_mm256_add_ps(
                                    _mm256_mul_ps(x, x), _mm256_mul_ps(y, y))));
    }
}


/* Normalize: the zero-vector guard is a mask rather than a branch.
 * keep = !(len <= FV_EPS), strictly complementary to the scalar if (len <= FV_EPS) (a NaN takes the else branch,
 * which only matches with _CMP_NLE_UQ: unordered also counts as "not satisfying <=").
 * The dividend of the division is blend-ed to 1.0, so a zero-vector lane produces no inf/NaN and sets no divide-by-zero flag;
 * finally and(keep) clears that lane to +0.0, matching the bit pattern the scalar returns for (0,0). */
FV_I_TARGET_AVX2 static inline void fv_i_normalize_soa_avx2(const float *__restrict ax,
                                                            const float *__restrict ay,
                                                            float *__restrict ox,
                                                            float *__restrict oy, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
    }
}


/* dot: 4 independent accumulators break the dependency chain (early microbenchmark: one accumulator 0.3213 -> four accumulators 0.1428 ns/element,
 * not re-measured with this round's hardened fixtures; this only says "the dependency chain is the bottleneck" and is not a gain promise for dot_soa),
 * plus one final reduction. The x and y components accumulate into the same accumulator, so the reduction order differs from the per-element path --
 * the only deviation allowed for fv_vec2_dot_soa (already stated in this header). */
FV_I_TARGET_AVX2 static inline float fv_i_dot_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict bx,
                                                       const float *__restrict by, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    __m256 s2 = _mm256_setzero_ps();
    __m256 s3 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 32 <= n; i += 32) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 0), _mm256_loadu_ps(bx + i + 0)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 0), _mm256_loadu_ps(by + i + 0))));
        s1 = _mm256_add_ps(s1, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 8), _mm256_loadu_ps(bx + i + 8)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 8), _mm256_loadu_ps(by + i + 8))));
        s2 = _mm256_add_ps(s2, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 16), _mm256_loadu_ps(bx + i + 16)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 16), _mm256_loadu_ps(by + i + 16))));
        s3 = _mm256_add_ps(s3, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 24), _mm256_loadu_ps(bx + i + 24)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 24), _mm256_loadu_ps(by + i + 24))));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i), _mm256_loadu_ps(bx + i)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i), _mm256_loadu_ps(by + i))));
    }

    /* Final reduction: first add the 4 ymm accumulators pairwise, then one 4-lane tree.
     * movehdup + shuffle rather than hadd (hadd is about 2x slower on this machine).
     * q = [q0,q1,q2,q3] -> r = [q0+q1, ..., q2+q3, ...] -> take lane0 + lane2. */
    {
        __m256 t = _mm256_add_ps(_mm256_add_ps(s0, s1), _mm256_add_ps(s2, s3));
        __m128 lo = _mm256_castps256_ps128(t);
        __m128 hi = _mm256_extractf128_ps(t, 1);
        __m128 q = _mm_add_ps(lo, hi);
        __m128 r = _mm_add_ps(q, _mm_movehdup_ps(q));
        sum = _mm_cvtss_f32(_mm_add_ss(r, _mm_shuffle_ps(r, r, 0x02)));
    }

    /* The tail accumulates in order (for n < 8 the whole span comes here, giving results bit-identical to the per-element path). */
    for (; i < n; ++i) {
        sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
    }
    return sum;
}
/* The streaming normalize kernel: the arithmetic path is word-for-word that of fv_i_normalize_soa_avx2 (zero-vector mask guard included),
 * with only the two output stores turned into NT stores. ox and oy share one i, so the peel aims to bring both to the 32 B
 * boundary **at the same time**; when that fails (the two outputs' phase difference within 32 B is not a multiple of 8) it falls back to normal storeu with a bit-identical result.
 * Allocating both ox / oy at 32 B is recommended, so no fallback is ever needed. */
FV_I_TARGET_AVX2 static inline void fv_i_normalize_stream_soa_avx2(const float *__restrict ax,
                                                                  const float *__restrict ay,
                                                                  float *__restrict ox,
                                                                  float *__restrict oy, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    int nt;
    size_t i = 0;

    for (; i < n && i < 7u
           && ((((uintptr_t)(const void *)(ox + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oy + i)) & 31u) != 0u); ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
    nt = ((((uintptr_t)(const void *)(ox + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oy + i)) & 31u) == 0u) ? 1 : 0;

    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        __m256 rx = _mm256_and_ps(_mm256_mul_ps(x, inv), keep);
        __m256 ry = _mm256_and_ps(_mm256_mul_ps(y, inv), keep);
        if (nt != 0) {
            _mm256_stream_ps(ox + i, rx);
            _mm256_stream_ps(oy + i, ry);
        } else {
            _mm256_storeu_ps(ox + i, rx);
            _mm256_storeu_ps(oy + i, ry);
        }
    }
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
    _mm_sfence();
}


/* Fused pipeline kernel: u = a * s + b, out = normalize(u), one pass.
 * Two properties make it bit-identical to scale -> add -> normalize, and both are visible here:
 *   1) order: _mm256_mul_ps first, _mm256_add_ps second, exactly like fv_i_scale_soa_avx2 then fv_i_add_soa_avx2.
 *      Consecutive float stores and reloads of the intermediate are exact, so dropping the intermediate array moves no bit.
 *      The pragma around this section is what keeps mul+add from being contracted into vfmadd (which would change the bits
 *      under a consumer's -mfma -ffp-contract=fast).
 *   2) the normalize part is word for word fv_i_normalize_soa_avx2: len = sqrt(x*x + y*y), keep = len _CMP_NLE_UQ FV_EPS,
 *      the divisor blended to 1.0 (so a zero-vector lane neither divides by zero nor sets a flag), then and(keep) which
 *      turns that lane into +0.0 -- the same mask behaviour, including inf / NaN passing through the comparison. */
FV_I_TARGET_AVX2 static inline void fv_i_madd_norm_soa_avx2(const float *__restrict ax,
                                                            const float *__restrict ay,
                                                            const float *__restrict bx,
                                                            const float *__restrict by,
                                                            fv_real s, float *__restrict ox,
                                                            float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ax + i), vs),
                                 _mm256_loadu_ps(bx + i));
        __m256 y = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ay + i), vs),
                                 _mm256_loadu_ps(by + i));
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
    }
}


#if FV_I_AVX512
/* 16-lane version of the same arithmetic. One vdivps per 16 elements; the zero guard is a mask, so no
 * branch. This path is selected only when __builtin_cpu_supports("avx512f") is true. */
FV_I_TARGET_AVX512 static inline void fv_i_madd_norm_soa_avx512(const float *__restrict ax,
                                                                 const float *__restrict ay,
                                                                 const float *__restrict bx,
                                                                 const float *__restrict by,
                                                                 fv_real s,
                                                                 float *__restrict ox,
                                                                 float *__restrict oy, size_t n)
{
    const __m512 vs = _mm512_set1_ps(s);
    const __m512 eps = _mm512_set1_ps((float)FV_EPS);
    const __m512 one = _mm512_set1_ps((float)1);
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512 x = _mm512_add_ps(_mm512_mul_ps(_mm512_loadu_ps(ax + i), vs),
                                 _mm512_loadu_ps(bx + i));
        __m512 y = _mm512_add_ps(_mm512_mul_ps(_mm512_loadu_ps(ay + i), vs),
                                 _mm512_loadu_ps(by + i));
        __m512 len = _mm512_sqrt_ps(_mm512_add_ps(_mm512_mul_ps(x, x),
                                                  _mm512_mul_ps(y, y)));
        __mmask16 keep = _mm512_cmp_ps_mask(len, eps, _CMP_NLE_UQ);
        __m512 safe = _mm512_mask_blend_ps(keep, one, len);
        __m512 inv = _mm512_div_ps(one, safe);
        _mm512_storeu_ps(ox + i, _mm512_maskz_mul_ps(keep, x, inv));
        _mm512_storeu_ps(oy + i, _mm512_maskz_mul_ps(keep, y, inv));
    }
}
#endif /* FV_I_AVX512 */

/* The streaming fused kernel: arithmetic word for word the one above (so the output is bit-identical), only the two output
 * stores become NT stores. Peel / fallback / sfence follow fv_i_normalize_stream_soa_avx2 -- ox and oy share one i here
 * because one output element needs both inputs, so both have to reach the 32 B boundary together; allocating both outputs
 * at 32 B is recommended so the fallback (normal storeu throughout, same bits) is never taken. */
FV_I_TARGET_AVX2 static inline void fv_i_madd_norm_stream_soa_avx2(const float *__restrict ax,
                                                                   const float *__restrict ay,
                                                                   const float *__restrict bx,
                                                                   const float *__restrict by,
                                                                   fv_real s, float *__restrict ox,
                                                                   float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    int nt;
    size_t i = 0;

    for (; i < n && i < 7u
           && ((((uintptr_t)(const void *)(ox + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oy + i)) & 31u) != 0u); ++i) {
        fv_i_s_madd_norm(ax[i], ay[i], bx[i], by[i], s, &ox[i], &oy[i]);
    }
    nt = ((((uintptr_t)(const void *)(ox + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oy + i)) & 31u) == 0u) ? 1 : 0;

    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ax + i), vs),
                                 _mm256_loadu_ps(bx + i));
        __m256 y = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ay + i), vs),
                                 _mm256_loadu_ps(by + i));
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        __m256 rx = _mm256_and_ps(_mm256_mul_ps(x, inv), keep);
        __m256 ry = _mm256_and_ps(_mm256_mul_ps(y, inv), keep);
        if (nt != 0) {
            _mm256_stream_ps(ox + i, rx);
            _mm256_stream_ps(oy + i, ry);
        } else {
            _mm256_storeu_ps(ox + i, rx);
            _mm256_storeu_ps(oy + i, ry);
        }
    }
    for (; i < n; ++i) {
        fv_i_s_madd_norm(ax[i], ay[i], bx[i], by[i], s, &ox[i], &oy[i]);
    }
    _mm_sfence();
}


/* ---- 4D AoS lerp / madd kernels: 2 elements at a time, zero repacking (like _x2's add/sub/scale) ---- */
FV_I_TARGET_AVX2 static inline void fv_i_v4_lerp2_avx2(const fv_vec4 *__restrict a,
                                                     const fv_vec4 *__restrict b, fv_real t,
                                                     fv_vec4 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_madd2_avx2(const fv_vec4 *__restrict a,
                                                     const fv_vec4 *__restrict b, fv_real s,
                                                     fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
}


/* The same two kernels with NT stores (fv_vec4 _x2 shape: 16 B stride, one peel step). Inside the pragma for the
 * same reason as the vec2 versions: mul feeding straight into add must not become an FMA. */
FV_I_TARGET_AVX2 static inline void fv_i_v4_lerp2_stream_avx2(const fv_vec4 *__restrict a,
                                                              const fv_vec4 *__restrict b,
                                                              fv_real t,
                                                              fv_vec4 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_lerp4(a[i], b[i], t);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_v4_lerp2_avx2(a, b, t, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp4(a[i], b[i], t);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_madd2_stream_avx2(const fv_vec4 *__restrict a,
                                                              const fv_vec4 *__restrict b,
                                                              fv_real s,
                                                              fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_madd4(a[i], b[i], s);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_v4_madd2_avx2(a, b, s, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd4(a[i], b[i], s);
    }
    _mm_sfence();
}


/* ---- the three per-lane SoA primitives: the single place where the association order is decided ----
 * fv_i_v3_mul8 / fv_i_v4_mul8 compute "the per-lane product terms of 8 elements in one ymm", in an association order that corresponds
 * word for word to the scalar API:
 *   vec3: (x*x' + y*y') + z*z'
 *   vec4: ((x*x' + y*y') + z*z') + w*w'
 * The dot / len / normalize kernels all take their values through them, so the association order has to be matched only here, once.
 * The parameters **deliberately do not carry** __restrict: len's call site is fv_i_v3_mul8(ax, ay, az, ax, ay, az, i),
 * where the pointers on both sides point into the same memory.
 * These two functions leave no room for contraction between mul and add -- the mul result feeds straight into the add -- so they must sit
 * inside the pragma (fp-contract=off); otherwise the consumer's -mfma -ffp-contract=fast turns them into vfmadd. */
FV_I_TARGET_AVX2 static inline __m256 fv_i_v3_mul8(const float *ax, const float *ay,
                                                 const float *az,
                                                 const float *bx, const float *by,
                                                 const float *bz, size_t off)
{
    __m256 px = _mm256_mul_ps(_mm256_loadu_ps(ax + off), _mm256_loadu_ps(bx + off));
    __m256 py = _mm256_mul_ps(_mm256_loadu_ps(ay + off), _mm256_loadu_ps(by + off));
    __m256 pz = _mm256_mul_ps(_mm256_loadu_ps(az + off), _mm256_loadu_ps(bz + off));
    return _mm256_add_ps(_mm256_add_ps(px, py), pz);
}


FV_I_TARGET_AVX2 static inline __m256 fv_i_v4_mul8(const float *ax, const float *ay,
                                                 const float *az, const float *aw,
                                                 const float *bx, const float *by,
                                                 const float *bz, const float *bw, size_t off)
{
    __m256 px = _mm256_mul_ps(_mm256_loadu_ps(ax + off), _mm256_loadu_ps(bx + off));
    __m256 py = _mm256_mul_ps(_mm256_loadu_ps(ay + off), _mm256_loadu_ps(by + off));
    __m256 pz = _mm256_mul_ps(_mm256_loadu_ps(az + off), _mm256_loadu_ps(bz + off));
    __m256 pw = _mm256_mul_ps(_mm256_loadu_ps(aw + off), _mm256_loadu_ps(bw + off));
    return _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(px, py), pz), pw);
}


/* Two ymm accumulators -> a scalar. Written like the final reduction of fv_i_dot_soa_avx2 (movehdup + shuffle,
 * no hadd; the reason is recorded in that comment). The reduction order differs from the per-element path, so dot_soa promises only determinism
 * and a small error, not bit-exactness (as the contract states). */
FV_I_TARGET_AVX2 static inline float fv_i_v34_reduce2(__m256 s0, __m256 s1)
{
    __m256 t = _mm256_add_ps(s0, s1);
    __m128 lo = _mm256_castps256_ps128(t);
    __m128 hi = _mm256_extractf128_ps(t, 1);
    __m128 q = _mm_add_ps(lo, hi);
    __m128 r = _mm_add_ps(q, _mm_movehdup_ps(q));
    return _mm_cvtss_f32(_mm_add_ss(r, _mm_shuffle_ps(r, r, 0x02)));
}


FV_I_TARGET_AVX2 static inline void fv_i_v3_len_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(fv_i_v3_mul8(ax, ay, az, ax, ay, az, i)));
    }
}


/* Normalize: the zero-vector guard is a mask rather than a branch, with the same semantics as fv_vec3_normalize.
 * keep = !(len <= FV_EPS) (_CMP_NLE_UQ: unordered also counts as "not satisfying <=", so a NaN takes the else branch,
 * strictly complementary to the scalar if (len <= FV_EPS)); the dividend is blend-ed to 1.0 on a zero-vector lane,
 * so no inf/NaN appears and the MXCSR divide-by-zero flag stays clear; finally and(keep) clears that lane to +0.0. */
FV_I_TARGET_AVX2 static inline void fv_i_v3_normalize_soa_avx2(const float *__restrict ax,
                                                             const float *__restrict ay,
                                                             const float *__restrict az,
                                                             float *__restrict ox,
                                                             float *__restrict oy,
                                                             float *__restrict oz, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v3_mul8(ax, ay, az, ax, ay, az, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
        _mm256_storeu_ps(oz + i, _mm256_and_ps(_mm256_mul_ps(z, inv), keep));
    }
}


/* The streaming 3D normalize kernel: the arithmetic path is word for word that of fv_i_v3_normalize_soa_avx2 (same
 * fv_i_v3_mul8 association order, same zero-vector mask guard), with only the three output stores turned into NT stores.
 * ox / oy / oz share one i, so the peel aims to bring all three to the 32 B boundary **at the same time**; when that fails
 * (the phase differences within 32 B are not all multiples of 8) it falls back to normal storeu for the whole span, with a
 * bit-identical result. Allocating the three outputs at 32 B is recommended, so no fallback is ever needed.
 * Inside the pragma: x*x + y*y + z*z is a mul feeding straight into an add. */
FV_I_TARGET_AVX2 static inline void fv_i_v3_normalize_stream_soa_avx2(const float *__restrict ax,
                                                                     const float *__restrict ay,
                                                                     const float *__restrict az,
                                                                     float *__restrict ox,
                                                                     float *__restrict oy,
                                                                     float *__restrict oz, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    int nt;
    size_t i = 0;

    for (; i < n && i < 7u
           && ((((uintptr_t)(const void *)(ox + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oy + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oz + i)) & 31u) != 0u); ++i) {
        fv_i_s_normalize3(ax[i], ay[i], az[i], &ox[i], &oy[i], &oz[i]);
    }
    nt = ((((uintptr_t)(const void *)(ox + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oy + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oz + i)) & 31u) == 0u) ? 1 : 0;

    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v3_mul8(ax, ay, az, ax, ay, az, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        __m256 rx = _mm256_and_ps(_mm256_mul_ps(x, inv), keep);
        __m256 ry = _mm256_and_ps(_mm256_mul_ps(y, inv), keep);
        __m256 rz = _mm256_and_ps(_mm256_mul_ps(z, inv), keep);
        if (nt != 0) {
            _mm256_stream_ps(ox + i, rx);
            _mm256_stream_ps(oy + i, ry);
            _mm256_stream_ps(oz + i, rz);
        } else {
            _mm256_storeu_ps(ox + i, rx);
            _mm256_storeu_ps(oy + i, ry);
            _mm256_storeu_ps(oz + i, rz);
        }
    }
    for (; i < n; ++i) {
        fv_i_s_normalize3(ax[i], ay[i], az[i], &ox[i], &oy[i], &oz[i]);
    }
    _mm_sfence();
}


/* dot: two independent accumulators break the dependency chain (16 elements per iteration), plus one final reduction.
 * Over 3 components each iteration has 3 vmulps + 2 vaddps available in parallel, while the accumulation chain has only 1 vaddps per iteration,
 * so two accumulators already keep that chain off the critical path (no timing this round, hence no measured numbers).
 * The reduction order differs from the per-element path -- the only deviation allowed for fv_vec3_dot_soa (already stated in this header's contract). */
FV_I_TARGET_AVX2 static inline float fv_i_v3_dot_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict bx,
                                                        const float *__restrict by,
                                                        const float *__restrict bz, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        s0 = _mm256_add_ps(s0, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i));
        s1 = _mm256_add_ps(s1, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i + 8));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i));
    }

    sum = fv_i_v34_reduce2(s0, s1);

    /* The tail accumulates in order, with expressions word-for-word the left-to-right order of fv_vec3_dot. */
    for (; i < n; ++i) {
        sum = sum + ((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]);
    }
    return sum;
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_len_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       const float *__restrict aw,
                                                       float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(fv_i_v4_mul8(ax, ay, az, aw, ax, ay, az, aw, i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_normalize_soa_avx2(const float *__restrict ax,
                                                             const float *__restrict ay,
                                                             const float *__restrict az,
                                                             const float *__restrict aw,
                                                             float *__restrict ox,
                                                             float *__restrict oy,
                                                             float *__restrict oz,
                                                             float *__restrict ow, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v4_mul8(ax, ay, az, aw, ax, ay, az, aw, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
        _mm256_storeu_ps(oz + i, _mm256_and_ps(_mm256_mul_ps(z, inv), keep));
        _mm256_storeu_ps(ow + i, _mm256_and_ps(_mm256_mul_ps(w, inv), keep));
    }
}


/* The streaming 4D normalize kernel: the arithmetic path is word for word that of fv_i_v4_normalize_soa_avx2, with only
 * the four output stores turned into NT stores. All four outputs share one i, so the peel must bring all four to the 32 B
 * boundary at the same time; if it cannot, the whole span goes through normal storeu (bit-identical). Inside the pragma:
 * ((x*x' + y*y') + z*z') + w*w' is a mul feeding straight into an add. */
FV_I_TARGET_AVX2 static inline void fv_i_v4_normalize_stream_soa_avx2(const float *__restrict ax,
                                                                     const float *__restrict ay,
                                                                     const float *__restrict az,
                                                                     const float *__restrict aw,
                                                                     float *__restrict ox,
                                                                     float *__restrict oy,
                                                                     float *__restrict oz,
                                                                     float *__restrict ow, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    int nt;
    size_t i = 0;

    for (; i < n && i < 7u
           && ((((uintptr_t)(const void *)(ox + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oy + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oz + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(ow + i)) & 31u) != 0u); ++i) {
        fv_i_s_normalize4(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
    nt = ((((uintptr_t)(const void *)(ox + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oy + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oz + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(ow + i)) & 31u) == 0u) ? 1 : 0;

    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v4_mul8(ax, ay, az, aw, ax, ay, az, aw, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        __m256 rx = _mm256_and_ps(_mm256_mul_ps(x, inv), keep);
        __m256 ry = _mm256_and_ps(_mm256_mul_ps(y, inv), keep);
        __m256 rz = _mm256_and_ps(_mm256_mul_ps(z, inv), keep);
        __m256 rw = _mm256_and_ps(_mm256_mul_ps(w, inv), keep);
        if (nt != 0) {
            _mm256_stream_ps(ox + i, rx);
            _mm256_stream_ps(oy + i, ry);
            _mm256_stream_ps(oz + i, rz);
            _mm256_stream_ps(ow + i, rw);
        } else {
            _mm256_storeu_ps(ox + i, rx);
            _mm256_storeu_ps(oy + i, ry);
            _mm256_storeu_ps(oz + i, rz);
            _mm256_storeu_ps(ow + i, rw);
        }
    }
    for (; i < n; ++i) {
        fv_i_s_normalize4(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline float fv_i_v4_dot_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict aw,
                                                        const float *__restrict bx,
                                                        const float *__restrict by,
                                                        const float *__restrict bz,
                                                        const float *__restrict bw, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        s0 = _mm256_add_ps(s0, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i));
        s1 = _mm256_add_ps(s1, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i + 8));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i));
    }

    sum = fv_i_v34_reduce2(s0, s1);

    for (; i < n; ++i) {
        sum = sum + (((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]) + aw[i] * bw[i]);
    }
    return sum;
}


/* ---- Matrix transform: 2 fv_vec4 per ymm (32 B, zero reshuffle on load and store) ----
 * The coefficient vectors duplicate the 3x4 part per 128-bit half, so lanes 0..3 serve point i
 * and lanes 4..7 serve point i+1: c_k = (m[0*4+k], m[1*4+k], m[2*4+k], -, same again).
 * Three _mm256_permutevar8x32_ps broadcasts turn the two AoS points into x / y / z vectors --
 * that is the whole gain: the naive AoS loop the compiler auto-vectorizes needs 14 vpermps and
 * 440 instructions per iteration, this one 3 and ~107 (objdump on bench/tfprobe).
 * Lanes 3 and 7 (the w output) are overwritten by vblendps with a constant and never computed:
 * 0 * inf is NaN, so computing w would make it NaN for a non-finite input, while the scalar
 * tail writes 1 / 0 unconditionally. Overwriting keeps both paths bit-identical.
 * The evaluation order is exactly ((c0*x + c1*y) + c2*z) + c3 -- one vmulps per term, one
 * vaddps per step, no contraction. This section is under fp-contract=off, which is the only
 * thing that stops -mfma on the consumer command line (target("avx2") does not). */
FV_I_TARGET_AVX2 static inline void fv_i_mat4_points2_avx2(const fv_mat4 *__restrict m,
                                                           const fv_vec4 *__restrict in,
                                                           fv_vec4 *__restrict out, size_t n)
{
    const __m256 c0 = _mm256_setr_ps(m->m[0], m->m[4], m->m[8], (fv_real)0,
                                     m->m[0], m->m[4], m->m[8], (fv_real)0);
    const __m256 c1 = _mm256_setr_ps(m->m[1], m->m[5], m->m[9], (fv_real)0,
                                     m->m[1], m->m[5], m->m[9], (fv_real)0);
    const __m256 c2 = _mm256_setr_ps(m->m[2], m->m[6], m->m[10], (fv_real)0,
                                     m->m[2], m->m[6], m->m[10], (fv_real)0);
    const __m256 c3 = _mm256_setr_ps(m->m[3], m->m[7], m->m[11], (fv_real)0,
                                     m->m[3], m->m[7], m->m[11], (fv_real)0);
    const __m256 wone = _mm256_set1_ps((fv_real)1);
    const __m256i ix = _mm256_setr_epi32(0, 0, 0, 0, 4, 4, 4, 4);
    const __m256i iy = _mm256_setr_epi32(1, 1, 1, 1, 5, 5, 5, 5);
    const __m256i iz = _mm256_setr_epi32(2, 2, 2, 2, 6, 6, 6, 6);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        const __m256 p = _mm256_loadu_ps((const float *)(const void *)(in + i));
        __m256 r = _mm256_mul_ps(c0, _mm256_permutevar8x32_ps(p, ix));
        r = _mm256_add_ps(r, _mm256_mul_ps(c1, _mm256_permutevar8x32_ps(p, iy)));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2, _mm256_permutevar8x32_ps(p, iz)));
        r = _mm256_add_ps(r, c3);
        _mm256_storeu_ps((float *)(void *)(out + i), _mm256_blend_ps(r, wone, 0x88));
    }
}


/* Same kernel without the translation column; the w lanes are blended to +0.0. */
FV_I_TARGET_AVX2 static inline void fv_i_mat4_dirs2_avx2(const fv_mat4 *__restrict m,
                                                         const fv_vec4 *__restrict in,
                                                         fv_vec4 *__restrict out, size_t n)
{
    const __m256 c0 = _mm256_setr_ps(m->m[0], m->m[4], m->m[8], (fv_real)0,
                                     m->m[0], m->m[4], m->m[8], (fv_real)0);
    const __m256 c1 = _mm256_setr_ps(m->m[1], m->m[5], m->m[9], (fv_real)0,
                                     m->m[1], m->m[5], m->m[9], (fv_real)0);
    const __m256 c2 = _mm256_setr_ps(m->m[2], m->m[6], m->m[10], (fv_real)0,
                                     m->m[2], m->m[6], m->m[10], (fv_real)0);
    const __m256 wzero = _mm256_setzero_ps();
    const __m256i ix = _mm256_setr_epi32(0, 0, 0, 0, 4, 4, 4, 4);
    const __m256i iy = _mm256_setr_epi32(1, 1, 1, 1, 5, 5, 5, 5);
    const __m256i iz = _mm256_setr_epi32(2, 2, 2, 2, 6, 6, 6, 6);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        const __m256 p = _mm256_loadu_ps((const float *)(const void *)(in + i));
        __m256 r = _mm256_mul_ps(c0, _mm256_permutevar8x32_ps(p, ix));
        r = _mm256_add_ps(r, _mm256_mul_ps(c1, _mm256_permutevar8x32_ps(p, iy)));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2, _mm256_permutevar8x32_ps(p, iz)));
        _mm256_storeu_ps((float *)(void *)(out + i), _mm256_blend_ps(r, wzero, 0x88));
    }
}


/* ---- the same two kernels with NT stores (fv_vec4 _x2 shape: 16 B stride, one peel step) ----
 * The arithmetic is copied verbatim from fv_i_mat4_points2_avx2 / fv_i_mat4_dirs2_avx2 -- the matrix terms are
 * m*x + ... (a mul feeding straight into an add), so both sit inside the pragma; only the store instruction differs. */
FV_I_TARGET_AVX2 static inline void fv_i_mat4_points2_stream_avx2(const fv_mat4 *__restrict m,
                                                                  const fv_vec4 *__restrict in,
                                                                  fv_vec4 *__restrict out,
                                                                  size_t n)
{
    const __m256 c0 = _mm256_setr_ps(m->m[0], m->m[4], m->m[8], (fv_real)0,
                                     m->m[0], m->m[4], m->m[8], (fv_real)0);
    const __m256 c1 = _mm256_setr_ps(m->m[1], m->m[5], m->m[9], (fv_real)0,
                                     m->m[1], m->m[5], m->m[9], (fv_real)0);
    const __m256 c2 = _mm256_setr_ps(m->m[2], m->m[6], m->m[10], (fv_real)0,
                                     m->m[2], m->m[6], m->m[10], (fv_real)0);
    const __m256 c3 = _mm256_setr_ps(m->m[3], m->m[7], m->m[11], (fv_real)0,
                                     m->m[3], m->m[7], m->m[11], (fv_real)0);
    const __m256 wone = _mm256_set1_ps((fv_real)1);
    const __m256i ix = _mm256_setr_epi32(0, 0, 0, 0, 4, 4, 4, 4);
    const __m256i iy = _mm256_setr_epi32(1, 1, 1, 1, 5, 5, 5, 5);
    const __m256i iz = _mm256_setr_epi32(2, 2, 2, 2, 6, 6, 6, 6);
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(out + i)) & 31u) != 0u; ++i) {
        out[i] = fv_i_s_mat4_point(m, in[i]);
    }
    if ((((uintptr_t)(const void *)(out + i)) & 31u) != 0u) {
        fv_i_mat4_points2_avx2(m, in, out, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        const __m256 p = _mm256_loadu_ps((const float *)(const void *)(in + i));
        __m256 r = _mm256_mul_ps(c0, _mm256_permutevar8x32_ps(p, ix));
        r = _mm256_add_ps(r, _mm256_mul_ps(c1, _mm256_permutevar8x32_ps(p, iy)));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2, _mm256_permutevar8x32_ps(p, iz)));
        r = _mm256_add_ps(r, c3);
        _mm256_stream_ps((float *)(void *)(out + i), _mm256_blend_ps(r, wone, 0x88));
    }
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_point(m, in[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_mat4_dirs2_stream_avx2(const fv_mat4 *__restrict m,
                                                                const fv_vec4 *__restrict in,
                                                                fv_vec4 *__restrict out,
                                                                size_t n)
{
    const __m256 c0 = _mm256_setr_ps(m->m[0], m->m[4], m->m[8], (fv_real)0,
                                     m->m[0], m->m[4], m->m[8], (fv_real)0);
    const __m256 c1 = _mm256_setr_ps(m->m[1], m->m[5], m->m[9], (fv_real)0,
                                     m->m[1], m->m[5], m->m[9], (fv_real)0);
    const __m256 c2 = _mm256_setr_ps(m->m[2], m->m[6], m->m[10], (fv_real)0,
                                     m->m[2], m->m[6], m->m[10], (fv_real)0);
    const __m256 wzero = _mm256_setzero_ps();
    const __m256i ix = _mm256_setr_epi32(0, 0, 0, 0, 4, 4, 4, 4);
    const __m256i iy = _mm256_setr_epi32(1, 1, 1, 1, 5, 5, 5, 5);
    const __m256i iz = _mm256_setr_epi32(2, 2, 2, 2, 6, 6, 6, 6);
    size_t i = 0;
    for (; i < n && i < 1u && (((uintptr_t)(const void *)(out + i)) & 31u) != 0u; ++i) {
        out[i] = fv_i_s_mat4_dir(m, in[i]);
    }
    if ((((uintptr_t)(const void *)(out + i)) & 31u) != 0u) {
        fv_i_mat4_dirs2_avx2(m, in, out, n);
        _mm_sfence();
        return;
    }
    for (; i + 2 <= n; i += 2) {
        const __m256 p = _mm256_loadu_ps((const float *)(const void *)(in + i));
        __m256 r = _mm256_mul_ps(c0, _mm256_permutevar8x32_ps(p, ix));
        r = _mm256_add_ps(r, _mm256_mul_ps(c1, _mm256_permutevar8x32_ps(p, iy)));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2, _mm256_permutevar8x32_ps(p, iz)));
        _mm256_stream_ps((float *)(void *)(out + i), _mm256_blend_ps(r, wzero, 0x88));
    }
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_dir(m, in[i]);
    }
    _mm_sfence();
}


#endif /* FV_I_INTRINSICS: kernels of the second section */


static inline void fv_i_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                fv_real t, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_lerp4_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp(a[i], b[i], t);
    }
}


static inline void fv_i_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                fv_real s, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_madd4_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd(a[i], b[i], s);
    }
}


/* Streaming lerp / madd (second section): the entry shape is word for word that of the non-stream entries above
 * (bulk = n & ~3 + the same scalar tail helper), only the kernel differs. They stay inside the pragma because their
 * bodies call fv_i_s_lerp / fv_i_s_madd, which contain mul+add. */
static inline void fv_i_lerp_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                       fv_real t, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_lerp4_stream_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp(a[i], b[i], t);
    }
}


static inline void fv_i_madd_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                       fv_real s, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_madd4_stream_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd(a[i], b[i], s);
    }
}


static inline float fv_i_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict bx, const float *__restrict by,
                                 size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_dot_soa_avx2(ax, ay, bx, by, n);
    }
#endif
    /* No AVX2: sequential summation with one accumulator, bit-identical to the per-element path. */
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
        }
        return sum;
    }
}


static inline void fv_i_len_soa(const float *__restrict ax, const float *__restrict ay,
                                float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_len_soa_avx2(ax, ay, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len(ax[i], ay[i]);
    }
}


static inline void fv_i_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                      float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_normalize_soa_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
}


/* The streaming entry must stay in the second section: its scalar tail contains x*x + y*y (via fv_i_s_normalize),
 * and outside the pragma the caller's -ffp-contract=fast would contract it into an FMA. */
static inline void fv_i_normalize_soa_stream(const float *__restrict ax,
                                             const float *__restrict ay,
                                             float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_normalize_stream_soa_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
}

/* Fused entry. The bulk split is the one every entry of the three-call chain uses (n & ~7) and the tail is the same
 * scalar helper, so each element is computed in the same order the three calls would compute it:
 * bit-identical by construction, not by tolerance. bulk != 0 also keeps n == 0 from entering a kernel. */
static inline void fv_i_madd_norm_soa(const float *__restrict ax, const float *__restrict ay,
                                      const float *__restrict bx, const float *__restrict by,
                                      fv_real s, float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_AVX512
    if (__builtin_cpu_supports("avx512f") != 0) {
        size_t bulk = n & ~(size_t)15;
        if (bulk != 0) {
            fv_i_madd_norm_soa_avx512(ax, ay, bx, by, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
#if FV_I_INTRINSICS
    if (i == 0 && __builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_madd_norm_soa_avx2(ax, ay, bx, by, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_madd_norm(ax[i], ay[i], bx[i], by[i], s, &ox[i], &oy[i]);
    }
}


/* The streaming fused entry must stay in the second section as well: its scalar tail goes through fv_i_s_madd_norm,
 * which carries both the mul+add of a*s + b and the x*x + y*y of the length. */
static inline void fv_i_madd_norm_soa_stream(const float *__restrict ax,
                                             const float *__restrict ay,
                                             const float *__restrict bx,
                                             const float *__restrict by, fv_real s,
                                             float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_madd_norm_stream_soa_avx2(ax, ay, bx, by, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_madd_norm(ax[i], ay[i], bx[i], by[i], s, &ox[i], &oy[i]);
    }
}

static inline void fv_i_v4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                 fv_real t, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_lerp2_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp4(a[i], b[i], t);
    }
}


static inline void fv_i_v4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                 fv_real s, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_madd2_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd4(a[i], b[i], s);
    }
}


/* Streaming 4D lerp / madd (second section, same reason as the vec2 versions: the tail helper contains mul+add). */
static inline void fv_i_v4_lerp_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                          fv_real t, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_lerp2_stream_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp4(a[i], b[i], t);
    }
}


static inline void fv_i_v4_madd_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                          fv_real s, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_madd2_stream_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd4(a[i], b[i], s);
    }
}


static inline float fv_i_v3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az,
                                  const float *__restrict bx, const float *__restrict by,
                                  const float *__restrict bz, size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_v3_dot_soa_avx2(ax, ay, az, bx, by, bz, n);
    }
#endif
    /* No AVX2: sequential summation with one accumulator, bit-identical to the per-element path (the +0.0/-0.0 boundary is in the contract). */
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + ((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]);
        }
        return sum;
    }
}


static inline void fv_i_v3_len_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az, float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_len_soa_avx2(ax, ay, az, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len3(ax[i], ay[i], az[i]);
    }
}


static inline void fv_i_v3_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict az,
                                       float *__restrict ox, float *__restrict oy,
                                       float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_normalize_soa_avx2(ax, ay, az, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize3(ax[i], ay[i], az[i], &ox[i], &oy[i], &oz[i]);
    }
}


/* Streaming 3D normalize (second section: its scalar tail, and the kernel's x*x + y*y + z*z, are mul+add).
 * Same bulk split and same scalar helper as the non-stream entry, so every element follows the identical order. */
static inline void fv_i_v3_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                const float *__restrict az,
                                                float *__restrict ox, float *__restrict oy,
                                                float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_normalize_stream_soa_avx2(ax, ay, az, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize3(ax[i], ay[i], az[i], &ox[i], &oy[i], &oz[i]);
    }
}


static inline float fv_i_v4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az, const float *__restrict aw,
                                  const float *__restrict bx, const float *__restrict by,
                                  const float *__restrict bz, const float *__restrict bw,
                                  size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_v4_dot_soa_avx2(ax, ay, az, aw, bx, by, bz, bw, n);
    }
#endif
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + (((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]) + aw[i] * bw[i]);
        }
        return sum;
    }
}


static inline void fv_i_v4_len_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az, const float *__restrict aw,
                                 float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v4_len_soa_avx2(ax, ay, az, aw, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len4(ax[i], ay[i], az[i], aw[i]);
    }
}


static inline void fv_i_v4_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict az, const float *__restrict aw,
                                       float *__restrict ox, float *__restrict oy,
                                       float *__restrict oz, float *__restrict ow, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v4_normalize_soa_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize4(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}


/* Streaming 4D normalize (second section, same reason as the 3D version). */
static inline void fv_i_v4_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                const float *__restrict az,
                                                const float *__restrict aw,
                                                float *__restrict ox, float *__restrict oy,
                                                float *__restrict oz, float *__restrict ow,
                                                size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v4_normalize_stream_soa_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize4(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}


/* Matrix transform: runtime dispatch on AVX2, 2 elements per ymm; everything else -- no AVX2,
 * and every tail length -- runs the scalar helper above, which *is* the reference expression.
 * bulk >= 2 also keeps n == 0 / n == 1 from touching the matrix at all: the kernel loads the
 * coefficient vectors before its loop, so it must not be entered for an empty batch. */
static inline void fv_i_mat4_tf_points(const fv_mat4 *__restrict m,
                                        const fv_vec4 *__restrict in,
                                        fv_vec4 *__restrict out, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_mat4_points2_avx2(m, in, out, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_point(m, in[i]);
    }
}


static inline void fv_i_mat4_tf_dirs(const fv_mat4 *__restrict m,
                                      const fv_vec4 *__restrict in,
                                      fv_vec4 *__restrict out, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_mat4_dirs2_avx2(m, in, out, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_dir(m, in[i]);
    }
}


/* Streaming matrix transform (second section: the scalar helpers are m*x + ... , i.e. mul feeding straight into add).
 * bulk >= 2 still keeps n == 0 / n == 1 from entering the kernel, which loads the coefficient vectors before its loop. */
static inline void fv_i_mat4_tf_points_stream(const fv_mat4 *__restrict m,
                                              const fv_vec4 *__restrict in,
                                              fv_vec4 *__restrict out, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_mat4_points2_stream_avx2(m, in, out, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_point(m, in[i]);
    }
}


static inline void fv_i_mat4_tf_dirs_stream(const fv_mat4 *__restrict m,
                                            const fv_vec4 *__restrict in,
                                            fv_vec4 *__restrict out, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_mat4_dirs2_stream_avx2(m, in, out, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        out[i] = fv_i_s_mat4_dir(m, in[i]);
    }
}


/* ======== integrator primitives: in-place x += v*s and the vector madd o = a*b + c ========
 * Every expression here is "a product feeding an addition", so the whole group sits inside the pragma. The byte
 * accounting, the reason addmul needs no _stream form and the in-place contract are stated in the declaration block. */

/* Scalar forms: one product, one sum, in that order -- the same shape as fv_i_s_madd above. */
static inline fv_vec2 fv_i_s_addmul(fv_vec2 x, fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    return r;
}


static inline fv_vec4 fv_i_s_addmul4(fv_vec4 x, fv_vec4 v, fv_real s)
{
    fv_vec4 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    r.z = x.z + v.z * s;
    r.w = x.w + v.w * s;
    return r;
}


static inline fv_vec2 fv_i_s_maddv(fv_vec2 a, fv_vec2 b, fv_vec2 c)
{
    fv_vec2 r;
    r.x = a.x * b.x + c.x;
    r.y = a.y * b.y + c.y;
    return r;
}


#if FV_I_INTRINSICS

/* 4 fv_vec2 (32 B) per step, in place: the load of x and the store back to x address the same cache line, so the line is
 * already owned when it is written -- no RFO, and the in-place shape is what makes the 24 B/element accounting possible. */
FV_I_TARGET_AVX2 static inline void fv_i_addmul4_avx2(fv_vec2 *__restrict x,
                                                      const fv_vec2 *__restrict v,
                                                      fv_real s, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 vx = _mm256_loadu_ps((const float *)(const void *)(x + i));
        __m256 vv = _mm256_loadu_ps((const float *)(const void *)(v + i));
        _mm256_storeu_ps((float *)(void *)(x + i), _mm256_add_ps(vx, _mm256_mul_ps(vv, vs)));
    }
}


/* SoA in place: each component array is its own read-modify-write stream, 8 floats per step. */
FV_I_TARGET_AVX2 static inline void fv_i_addmul_soa_avx2(float *__restrict x, float *__restrict y,
                                                         const float *__restrict vx,
                                                         const float *__restrict vy,
                                                         fv_real s, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(x + i, _mm256_add_ps(_mm256_loadu_ps(x + i),
                                              _mm256_mul_ps(_mm256_loadu_ps(vx + i), vs)));
        _mm256_storeu_ps(y + i, _mm256_add_ps(_mm256_loadu_ps(y + i),
                                              _mm256_mul_ps(_mm256_loadu_ps(vy + i), vs)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_addmul_v3_soa_avx2(float *__restrict x,
                                                            float *__restrict y,
                                                            float *__restrict z,
                                                            const float *__restrict vx,
                                                            const float *__restrict vy,
                                                            const float *__restrict vz,
                                                            fv_real s, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(x + i, _mm256_add_ps(_mm256_loadu_ps(x + i),
                                              _mm256_mul_ps(_mm256_loadu_ps(vx + i), vs)));
        _mm256_storeu_ps(y + i, _mm256_add_ps(_mm256_loadu_ps(y + i),
                                              _mm256_mul_ps(_mm256_loadu_ps(vy + i), vs)));
        _mm256_storeu_ps(z + i, _mm256_add_ps(_mm256_loadu_ps(z + i),
                                              _mm256_mul_ps(_mm256_loadu_ps(vz + i), vs)));
    }
}


/* fv_vec4: 16 B/element, so 2 elements are exactly one ymm -- same shape as the _x2 entries above. */
FV_I_TARGET_AVX2 static inline void fv_i_v4_addmul2_avx2(fv_vec4 *__restrict x,
                                                         const fv_vec4 *__restrict v,
                                                         fv_real s, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 vx = _mm256_loadu_ps((const float *)(const void *)(x + i));
        __m256 vv = _mm256_loadu_ps((const float *)(const void *)(v + i));
        _mm256_storeu_ps((float *)(void *)(x + i), _mm256_add_ps(vx, _mm256_mul_ps(vv, vs)));
    }
}


/* o = a*b + c into a fresh buffer: this one pays the RFO of o (40 B/element), which is what the _stream variant below
 * removes. The product feeds the sum directly, so the store order is "mulps first, then addps" -- never an FMA. */
FV_I_TARGET_AVX2 static inline void fv_i_maddv4_avx2(const fv_vec2 *__restrict a,
                                                     const fv_vec2 *__restrict b,
                                                     const fv_vec2 *__restrict c,
                                                     fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vc = _mm256_loadu_ps((const float *)(const void *)(c + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(_mm256_mul_ps(va, vb), vc));
    }
}


/* The same kernel with the output written by NT stores. Word for word the three-phase structure of the other _stream
 * kernels: peel with normal stores to the 32 B boundary (8 B stride -> at most 3 steps), body _mm256_stream_ps, scalar
 * tail; if the boundary is unreachable (o not even 8 B aligned) the whole span takes the normal storeu kernel, which
 * recomputes the peeled elements -- bit-identical, only the NT gain is lost. _mm_sfence() on every exit. */
FV_I_TARGET_AVX2 static inline void fv_i_maddv4_stream_avx2(const fv_vec2 *__restrict a,
                                                           const fv_vec2 *__restrict b,
                                                           const fv_vec2 *__restrict c,
                                                           fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_maddv(a[i], b[i], c[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_maddv4_avx2(a, b, c, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vc = _mm256_loadu_ps((const float *)(const void *)(c + i));
        _mm256_stream_ps((float *)(void *)(o + i), _mm256_add_ps(_mm256_mul_ps(va, vb), vc));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_maddv(a[i], b[i], c[i]);
    }
    _mm_sfence();
}

#endif /* FV_I_INTRINSICS: the AVX2 kernels of the integrator primitives */


/* ---- dispatch: same shape as every other batch entry (bulk + scalar tail, n == 0 touches nothing) ---- */
static inline void fv_i_addmul_x4(fv_vec2 *__restrict x, const fv_vec2 *__restrict v, fv_real s,
                                  size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_addmul4_avx2(x, v, s, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        x[i] = fv_i_s_addmul(x[i], v[i], s);
    }
}


static inline void fv_i_addmul_soa(float *__restrict x, float *__restrict y,
                                   const float *__restrict vx, const float *__restrict vy,
                                   fv_real s, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_addmul_soa_avx2(x, y, vx, vy, s, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        x[i] = x[i] + vx[i] * s;
        y[i] = y[i] + vy[i] * s;
    }
}


static inline void fv_i_addmul_v3_soa(float *__restrict x, float *__restrict y,
                                      float *__restrict z, const float *__restrict vx,
                                      const float *__restrict vy, const float *__restrict vz,
                                      fv_real s, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_addmul_v3_soa_avx2(x, y, z, vx, vy, vz, s, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        x[i] = x[i] + vx[i] * s;
        y[i] = y[i] + vy[i] * s;
        z[i] = z[i] + vz[i] * s;
    }
}


static inline void fv_i_v4_addmul_x2(fv_vec4 *__restrict x, const fv_vec4 *__restrict v,
                                     fv_real s, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_addmul2_avx2(x, v, s, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        x[i] = fv_i_s_addmul4(x[i], v[i], s);
    }
}


static inline void fv_i_maddv_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                 const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_maddv4_avx2(a, b, c, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_maddv(a[i], b[i], c[i]);
    }
}


static inline void fv_i_maddv_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                        const fv_vec2 *__restrict c, fv_vec2 *__restrict o,
                                        size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_maddv4_stream_avx2(a, b, c, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_maddv(a[i], b[i], c[i]);
    }
}


/* ---------------- quaternion / mat3 SoA batch: scalar helpers, kernels, dispatch (P6) ----------------
 * The expressions below are the single reference for bit-exactness: the per-value API in
 * src/flashvec_quat.cpp / src/flashvec_mat3.cpp writes the same lines, and this section's tail and
 * target("avx2") loop both call these helpers. The multiply / rotate / mat3-transform kernels are
 * deliberately plain C: they are elementwise loops with no reductions, so GCC's vectorizer is the first
 * implementation to measure; a hand-written kernel is added only if the instruction count / vpermps
 * evidence shows it failing (the mat4 rule). Quat normalize is that exception: the compiler kept
 * sqrtss/divss and a branch (0 ymm instructions), so it has an explicit kernel below. The pragma around
 * this section is what makes the result independent of a consumer's -mfma / -ffp-contract=fast, exactly
 * as it does for lerp / madd / normalize.
 */

/* Hamilton product, association order fixed: ((w1*x2 + x1*w2) + y1*z2) - z1*y2, and the three analogues. */
static inline void fv_i_s_quat_mul(fv_real ax, fv_real ay, fv_real az, fv_real aw,
                                   fv_real bx, fv_real by, fv_real bz, fv_real bw,
                                   fv_real *ox, fv_real *oy, fv_real *oz, fv_real *ow)
{
    *ox = ((aw * bx + ax * bw) + ay * bz) - az * by;
    *oy = ((aw * by - ax * bz) + ay * bw) + az * bx;
    *oz = ((aw * bz + ax * by) - ay * bx) + az * bw;
    *ow = ((aw * bw - ax * bx) - ay * by) - az * bz;
}

/* The same semantics as fv_quat_normalize: length <= FV_EPS -> identity (0,0,0,1). */
static inline void fv_i_s_quat_normalize(fv_real x, fv_real y, fv_real z, fv_real w,
                                         fv_real *ox, fv_real *oy, fv_real *oz, fv_real *ow)
{
    fv_real len = sqrtf(((x * x + y * y) + z * z) + w * w);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        *oz = (fv_real)0;
        *ow = (fv_real)1;
    } else {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
        *oz = z * inv;
        *ow = w * inv;
    }
}

/* v' = q*(v,0)*conjugate(q), written as t = 2*cross(q.xyz,v); d = cross(q.xyz,t); out = (v + q.w*t) + d.
 * q is assumed unit; no normalization is performed. */
static inline void fv_i_s_quat_rotate(fv_real qx, fv_real qy, fv_real qz, fv_real qw,
                                      fv_real vx, fv_real vy, fv_real vz,
                                      fv_real *ox, fv_real *oy, fv_real *oz)
{
    fv_real tx = (fv_real)2 * (qy * vz - qz * vy);
    fv_real ty = (fv_real)2 * (qz * vx - qx * vz);
    fv_real tz = (fv_real)2 * (qx * vy - qy * vx);
    fv_real dx = qy * tz - qz * ty;
    fv_real dy = qz * tx - qx * tz;
    fv_real dz = qx * ty - qy * tx;
    *ox = (vx + qw * tx) + dx;
    *oy = (vy + qw * ty) + dy;
    *oz = (vz + qw * tz) + dz;
}

/* out.x = (m00*vx + m01*vy) + m02*vz, and the two analogous rows. */
static inline void fv_i_s_mat3_tf3(fv_real m00, fv_real m01, fv_real m02,
                                   fv_real m10, fv_real m11, fv_real m12,
                                   fv_real m20, fv_real m21, fv_real m22,
                                   fv_real vx, fv_real vy, fv_real vz,
                                   fv_real *ox, fv_real *oy, fv_real *oz)
{
    *ox = (m00 * vx + m01 * vy) + m02 * vz;
    *oy = (m10 * vx + m11 * vy) + m12 * vz;
    *oz = (m20 * vx + m21 * vy) + m22 * vz;
}

#if FV_I_INTRINSICS

/* Plain C inside target("avx2"): GCC vectorizes the loop at 256 bits. Eight independent elements per
 * iteration (one per float lane); the helpers are static inline and the stores are to distinct arrays,
 * so the loop has no loop-carried dependency. */
FV_I_TARGET_AVX2 static inline void fv_i_quat_mul_soa_avx2(const float *__restrict ax,
                                                           const float *__restrict ay,
                                                           const float *__restrict az,
                                                           const float *__restrict aw,
                                                           const float *__restrict bx,
                                                           const float *__restrict by,
                                                           const float *__restrict bz,
                                                           const float *__restrict bw,
                                                           float *__restrict ox, float *__restrict oy,
                                                           float *__restrict oz, float *__restrict ow,
                                                           size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_i_s_quat_mul(ax[i], ay[i], az[i], aw[i], bx[i], by[i], bz[i], bw[i],
                        &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}

/* Hand-written because the compiler did not vectorize the plain-C form (measured: the loop kept
 * sqrtss/divss and a branch, 0 ymm instructions). The arithmetic is still the scalar expression, lane
 * by lane: ((x*x + y*y) + z*z) + w*w, then sqrt, then 1/len, then the four products. The zero-vector
 * guard is a blend with the identity, not a branch: a lane with len <= FV_EPS gets (0,0,0,1), every
 * other lane gets q/len. The division is evaluated unconditionally; for len == 0 it produces inf,
 * whose product is NaN, but that lane is replaced by the identity before the store, so no NaN escapes.
 * (The default fenv has FP exceptions masked; the library does not promise fenv flags.) */
FV_I_TARGET_AVX2 static inline void fv_i_quat_normalize_soa_avx2(const float *__restrict ax,
                                                                 const float *__restrict ay,
                                                                 const float *__restrict az,
                                                                 const float *__restrict aw,
                                                                 float *__restrict ox,
                                                                 float *__restrict oy,
                                                                 float *__restrict oz,
                                                                 float *__restrict ow, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((float)1);
    const __m256 zero = _mm256_setzero_ps();
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len2 = _mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y));
        len2 = _mm256_add_ps(len2, _mm256_mul_ps(z, z));
        len2 = _mm256_add_ps(len2, _mm256_mul_ps(w, w));
        {
            __m256 len = _mm256_sqrt_ps(len2);
            __m256 guard = _mm256_cmp_ps(len, eps, _CMP_LE_OQ);
            __m256 inv = _mm256_div_ps(one, len);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, inv), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, inv), zero, guard));
            _mm256_storeu_ps(oz + i, _mm256_blendv_ps(_mm256_mul_ps(z, inv), zero, guard));
            _mm256_storeu_ps(ow + i, _mm256_blendv_ps(_mm256_mul_ps(w, inv), one, guard));
        }
    }
}

FV_I_TARGET_AVX2 static inline void fv_i_quat_rotate_vec3_soa_avx2(const float *__restrict qx,
                                                                  const float *__restrict qy,
                                                                  const float *__restrict qz,
                                                                  const float *__restrict qw,
                                                                  const float *__restrict vx,
                                                                  const float *__restrict vy,
                                                                  const float *__restrict vz,
                                                                  float *__restrict ox,
                                                                  float *__restrict oy,
                                                                  float *__restrict oz, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_i_s_quat_rotate(qx[i], qy[i], qz[i], qw[i], vx[i], vy[i], vz[i],
                           &ox[i], &oy[i], &oz[i]);
    }
}

FV_I_TARGET_AVX2 static inline void fv_i_mat3_tf_soa_avx2(const float *__restrict m00,
                                                          const float *__restrict m01,
                                                          const float *__restrict m02,
                                                          const float *__restrict m10,
                                                          const float *__restrict m11,
                                                          const float *__restrict m12,
                                                          const float *__restrict m20,
                                                          const float *__restrict m21,
                                                          const float *__restrict m22,
                                                          const float *__restrict vx,
                                                          const float *__restrict vy,
                                                          const float *__restrict vz,
                                                          float *__restrict ox, float *__restrict oy,
                                                          float *__restrict oz, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_i_s_mat3_tf3(m00[i], m01[i], m02[i], m10[i], m11[i], m12[i],
                        m20[i], m21[i], m22[i], vx[i], vy[i], vz[i],
                        &ox[i], &oy[i], &oz[i]);
    }
}

#endif /* FV_I_INTRINSICS: quaternion / mat3 plain-C targets */

static inline void fv_i_quat_mul_soa(const float *__restrict ax, const float *__restrict ay,
                                     const float *__restrict az, const float *__restrict aw,
                                     const float *__restrict bx, const float *__restrict by,
                                     const float *__restrict bz, const float *__restrict bw,
                                     float *__restrict ox, float *__restrict oy,
                                     float *__restrict oz, float *__restrict ow, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_quat_mul_soa_avx2(ax, ay, az, aw, bx, by, bz, bw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_quat_mul(ax[i], ay[i], az[i], aw[i], bx[i], by[i], bz[i], bw[i],
                        &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}

static inline void fv_i_quat_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                           const float *__restrict az, const float *__restrict aw,
                                           float *__restrict ox, float *__restrict oy,
                                           float *__restrict oz, float *__restrict ow, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_quat_normalize_soa_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_quat_normalize(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}

static inline void fv_i_quat_rotate_vec3_soa(const float *__restrict qx,
                                             const float *__restrict qy,
                                             const float *__restrict qz,
                                             const float *__restrict qw,
                                             const float *__restrict vx,
                                             const float *__restrict vy,
                                             const float *__restrict vz,
                                             float *__restrict ox, float *__restrict oy,
                                             float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_quat_rotate_vec3_soa_avx2(qx, qy, qz, qw, vx, vy, vz, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_quat_rotate(qx[i], qy[i], qz[i], qw[i], vx[i], vy[i], vz[i],
                           &ox[i], &oy[i], &oz[i]);
    }
}

static inline void fv_i_mat3_transform_vec3_soa(const float *__restrict m00,
                                                const float *__restrict m01,
                                                const float *__restrict m02,
                                                const float *__restrict m10,
                                                const float *__restrict m11,
                                                const float *__restrict m12,
                                                const float *__restrict m20,
                                                const float *__restrict m21,
                                                const float *__restrict m22,
                                                const float *__restrict vx,
                                                const float *__restrict vy,
                                                const float *__restrict vz,
                                                float *__restrict ox, float *__restrict oy,
                                                float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_mat3_tf_soa_avx2(m00, m01, m02, m10, m11, m12, m20, m21, m22,
                                  vx, vy, vz, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_mat3_tf3(m00[i], m01[i], m02[i], m10[i], m11[i], m12[i],
                        m20[i], m21[i], m22[i], vx[i], vy[i], vz[i],
                        &ox[i], &oy[i], &oz[i]);
    }
}



/* ---------------- mat3 batch mul / inverse (P7) ----------------
 * AVX2 kernels: 8 elements per iteration. mul uses the 9 B columns as live registers and computes
 * the three output rows one at a time (9 B + 3 A + 3 O = 15 ymm, inside the 16-register file);
 * inverse computes the determinant and cofactors once per 8 elements, does a single vdivps(1/det)
 * instead of nine scalar divisions, and blends the zero matrix for singular lanes. Both fall back to
 * the scalar helper for the tail and are bit-identical to the per-value functions (the expressions
 * and association order are the scalar ones, lane by lane).
 */

static inline void fv_i_s_mat3_mul9(const float *__restrict a, const float *__restrict b,
                                    float *__restrict o)
{
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            o[i * 3 + j] = (a[i * 3 + 0] * b[0 * 3 + j] + a[i * 3 + 1] * b[1 * 3 + j]) +
                           a[i * 3 + 2] * b[2 * 3 + j];
        }
    }
}

/* The same determinant / scale threshold / cofactor expressions as fv_mat3_inverse. */
static inline int fv_i_s_mat3_inv9(const float *__restrict m, float *__restrict o)
{
    fv_real det = (m[0] * (m[4] * m[8] - m[5] * m[7]) -
                   m[1] * (m[3] * m[8] - m[5] * m[6])) +
                  m[2] * (m[3] * m[7] - m[4] * m[6]);
    fv_real scale = (fv_real)0;
    int k;
    for (k = 0; k < 9; ++k) {
        fv_real a = m[k] < (fv_real)0 ? -m[k] : m[k];
        if (a > scale) {
            scale = a;
        }
    }
    if (!(scale > (fv_real)0) ||
        !((det < (fv_real)0 ? -det : det) > FV_EPS * ((scale * scale) * scale))) {
        for (k = 0; k < 9; ++k) {
            o[k] = (fv_real)0;
        }
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / det;
        o[0] = (m[4] * m[8] - m[5] * m[7]) * inv;
        o[1] = (m[2] * m[7] - m[1] * m[8]) * inv;
        o[2] = (m[1] * m[5] - m[2] * m[4]) * inv;
        o[3] = (m[5] * m[6] - m[3] * m[8]) * inv;
        o[4] = (m[0] * m[8] - m[2] * m[6]) * inv;
        o[5] = (m[2] * m[3] - m[0] * m[5]) * inv;
        o[6] = (m[3] * m[7] - m[4] * m[6]) * inv;
        o[7] = (m[1] * m[6] - m[0] * m[7]) * inv;
        o[8] = (m[0] * m[4] - m[1] * m[3]) * inv;
        return 1;
    }
}

#if FV_I_INTRINSICS

FV_I_TARGET_AVX2 static inline void fv_i_mat3_mul_soa_avx2(                                                           const float *__restrict a00,
                                                           const float *__restrict a01,
                                                           const float *__restrict a02,
                                                           const float *__restrict a10,
                                                           const float *__restrict a11,
                                                           const float *__restrict a12,
                                                           const float *__restrict a20,
                                                           const float *__restrict a21,
                                                           const float *__restrict a22,
                                                           const float *__restrict b00,
                                                           const float *__restrict b01,
                                                           const float *__restrict b02,
                                                           const float *__restrict b10,
                                                           const float *__restrict b11,
                                                           const float *__restrict b12,
                                                           const float *__restrict b20,
                                                           const float *__restrict b21,
                                                           const float *__restrict b22,
                                                           float *__restrict o00,
                                                           float *__restrict o01,
                                                           float *__restrict o02,
                                                           float *__restrict o10,
                                                           float *__restrict o11,
                                                           float *__restrict o12,
                                                           float *__restrict o20,
                                                           float *__restrict o21,
                                                           float *__restrict o22,
                                                           size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 vb00 = _mm256_loadu_ps(b00 + i);
        __m256 vb01 = _mm256_loadu_ps(b01 + i);
        __m256 vb02 = _mm256_loadu_ps(b02 + i);
        __m256 vb10 = _mm256_loadu_ps(b10 + i);
        __m256 vb11 = _mm256_loadu_ps(b11 + i);
        __m256 vb12 = _mm256_loadu_ps(b12 + i);
        __m256 vb20 = _mm256_loadu_ps(b20 + i);
        __m256 vb21 = _mm256_loadu_ps(b21 + i);
        __m256 vb22 = _mm256_loadu_ps(b22 + i);
        {
            __m256 va00 = _mm256_loadu_ps(a00 + i);
            __m256 va01 = _mm256_loadu_ps(a01 + i);
            __m256 va02 = _mm256_loadu_ps(a02 + i);
            _mm256_storeu_ps(o00 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va00, vb00), _mm256_mul_ps(va01, vb10)), _mm256_mul_ps(va02, vb20)));
            _mm256_storeu_ps(o01 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va00, vb01), _mm256_mul_ps(va01, vb11)), _mm256_mul_ps(va02, vb21)));
            _mm256_storeu_ps(o02 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va00, vb02), _mm256_mul_ps(va01, vb12)), _mm256_mul_ps(va02, vb22)));
        }
        {
            __m256 va10 = _mm256_loadu_ps(a10 + i);
            __m256 va11 = _mm256_loadu_ps(a11 + i);
            __m256 va12 = _mm256_loadu_ps(a12 + i);
            _mm256_storeu_ps(o10 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va10, vb00), _mm256_mul_ps(va11, vb10)), _mm256_mul_ps(va12, vb20)));
            _mm256_storeu_ps(o11 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va10, vb01), _mm256_mul_ps(va11, vb11)), _mm256_mul_ps(va12, vb21)));
            _mm256_storeu_ps(o12 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va10, vb02), _mm256_mul_ps(va11, vb12)), _mm256_mul_ps(va12, vb22)));
        }
        {
            __m256 va20 = _mm256_loadu_ps(a20 + i);
            __m256 va21 = _mm256_loadu_ps(a21 + i);
            __m256 va22 = _mm256_loadu_ps(a22 + i);
            _mm256_storeu_ps(o20 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va20, vb00), _mm256_mul_ps(va21, vb10)), _mm256_mul_ps(va22, vb20)));
            _mm256_storeu_ps(o21 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va20, vb01), _mm256_mul_ps(va21, vb11)), _mm256_mul_ps(va22, vb21)));
            _mm256_storeu_ps(o22 + i, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(va20, vb02), _mm256_mul_ps(va21, vb12)), _mm256_mul_ps(va22, vb22)));
        }
    }
    for (; i < n; ++i) {
        float a[9];
        float b[9];
        float o[9];
        a[0] = a00[i]; a[1] = a01[i]; a[2] = a02[i];
        a[3] = a10[i]; a[4] = a11[i]; a[5] = a12[i];
        a[6] = a20[i]; a[7] = a21[i]; a[8] = a22[i];
        b[0] = b00[i]; b[1] = b01[i]; b[2] = b02[i];
        b[3] = b10[i]; b[4] = b11[i]; b[5] = b12[i];
        b[6] = b20[i]; b[7] = b21[i]; b[8] = b22[i];
        fv_i_s_mat3_mul9(a, b, o);
        o00[i] = o[0]; o01[i] = o[1]; o02[i] = o[2];
        o10[i] = o[3]; o11[i] = o[4]; o12[i] = o[5];
        o20[i] = o[6]; o21[i] = o[7]; o22[i] = o[8];
    }
}

FV_I_TARGET_AVX2 static inline void fv_i_mat3_inv_soa_avx2(const float *__restrict pm00,
                                                          const float *__restrict pm01,
                                                          const float *__restrict pm02,
                                                          const float *__restrict pm10,
                                                          const float *__restrict pm11,
                                                          const float *__restrict pm12,
                                                          const float *__restrict pm20,
                                                          const float *__restrict pm21,
                                                          const float *__restrict pm22,
                                                          float *__restrict o00,
                                                          float *__restrict o01,
                                                          float *__restrict o02,
                                                          float *__restrict o10,
                                                          float *__restrict o11,
                                                          float *__restrict o12,
                                                          float *__restrict o20,
                                                          float *__restrict o21,
                                                          float *__restrict o22,
                                                          size_t n)
{
    const __m256 absmask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    const __m256 zero = _mm256_setzero_ps();
    const __m256 one = _mm256_set1_ps((float)1);
    const __m256 eps = _mm256_set1_ps((float)FV_EPS);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 m00 = _mm256_loadu_ps(pm00 + i);
        __m256 m01 = _mm256_loadu_ps(pm01 + i);
        __m256 m02 = _mm256_loadu_ps(pm02 + i);
        __m256 m10 = _mm256_loadu_ps(pm10 + i);
        __m256 m11 = _mm256_loadu_ps(pm11 + i);
        __m256 m12 = _mm256_loadu_ps(pm12 + i);
        __m256 m20 = _mm256_loadu_ps(pm20 + i);
        __m256 m21 = _mm256_loadu_ps(pm21 + i);
        __m256 m22 = _mm256_loadu_ps(pm22 + i);
        __m256 det = _mm256_sub_ps(_mm256_mul_ps(m00, _mm256_sub_ps(_mm256_mul_ps(m11, m22), _mm256_mul_ps(m12, m21))),
                                   _mm256_mul_ps(m01, _mm256_sub_ps(_mm256_mul_ps(m10, m22), _mm256_mul_ps(m12, m20))));
        det = _mm256_add_ps(det, _mm256_mul_ps(m02, _mm256_sub_ps(_mm256_mul_ps(m10, m21), _mm256_mul_ps(m11, m20))));
        {
            __m256 s = _mm256_max_ps(_mm256_and_ps(m00, absmask), _mm256_and_ps(m01, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m02, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m10, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m11, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m12, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m20, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m21, absmask));
            s = _mm256_max_ps(s, _mm256_and_ps(m22, absmask));
            {
                __m256 scale2 = _mm256_mul_ps(s, s);
                __m256 thresh = _mm256_mul_ps(eps, _mm256_mul_ps(scale2, s));
                __m256 ok = _mm256_and_ps(_mm256_cmp_ps(s, zero, _CMP_GT_OQ),
                                          _mm256_cmp_ps(_mm256_and_ps(det, absmask), thresh, _CMP_GT_OQ));
                __m256 inv = _mm256_div_ps(one, det);
                _mm256_storeu_ps(o00 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m11, m22), _mm256_mul_ps(m12, m21)), inv), ok));
                _mm256_storeu_ps(o01 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m02, m21), _mm256_mul_ps(m01, m22)), inv), ok));
                _mm256_storeu_ps(o02 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m01, m12), _mm256_mul_ps(m02, m11)), inv), ok));
                _mm256_storeu_ps(o10 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m12, m20), _mm256_mul_ps(m10, m22)), inv), ok));
                _mm256_storeu_ps(o11 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m00, m22), _mm256_mul_ps(m02, m20)), inv), ok));
                _mm256_storeu_ps(o12 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m02, m10), _mm256_mul_ps(m00, m12)), inv), ok));
                _mm256_storeu_ps(o20 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m10, m21), _mm256_mul_ps(m11, m20)), inv), ok));
                _mm256_storeu_ps(o21 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m01, m20), _mm256_mul_ps(m00, m21)), inv), ok));
                _mm256_storeu_ps(o22 + i, _mm256_and_ps(_mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(m00, m11), _mm256_mul_ps(m01, m10)), inv), ok));
            }
        }
    }
    for (; i < n; ++i) {
        float m[9];
        float o[9];
        m[0] = pm00[i]; m[1] = pm01[i]; m[2] = pm02[i];
        m[3] = pm10[i]; m[4] = pm11[i]; m[5] = pm12[i];
        m[6] = pm20[i]; m[7] = pm21[i]; m[8] = pm22[i];
        (void)fv_i_s_mat3_inv9(m, o);
        o00[i] = o[0]; o01[i] = o[1]; o02[i] = o[2];
        o10[i] = o[3]; o11[i] = o[4]; o12[i] = o[5];
        o20[i] = o[6]; o21[i] = o[7]; o22[i] = o[8];
    }
}

#endif /* FV_I_INTRINSICS: mat3 batch kernels */

static inline void fv_i_mat3_mul_soa(                     const float *__restrict a00,
                     const float *__restrict a01,
                     const float *__restrict a02,
                     const float *__restrict a10,
                     const float *__restrict a11,
                     const float *__restrict a12,
                     const float *__restrict a20,
                     const float *__restrict a21,
                     const float *__restrict a22,
                     const float *__restrict b00,
                     const float *__restrict b01,
                     const float *__restrict b02,
                     const float *__restrict b10,
                     const float *__restrict b11,
                     const float *__restrict b12,
                     const float *__restrict b20,
                     const float *__restrict b21,
                     const float *__restrict b22,
                     float *__restrict o00,
                     float *__restrict o01,
                     float *__restrict o02,
                     float *__restrict o10,
                     float *__restrict o11,
                     float *__restrict o12,
                     float *__restrict o20,
                     float *__restrict o21,
                     float *__restrict o22,
                     size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_mat3_mul_soa_avx2(a00, a01, a02, a10, a11, a12, a20, a21, a22, b00, b01, b02, b10, b11, b12, b20, b21, b22, o00, o01, o02, o10, o11, o12, o20, o21, o22, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        float a[9];
        float b[9];
        float o[9];
        a[0] = a00[i]; a[1] = a01[i]; a[2] = a02[i];
        a[3] = a10[i]; a[4] = a11[i]; a[5] = a12[i];
        a[6] = a20[i]; a[7] = a21[i]; a[8] = a22[i];
        b[0] = b00[i]; b[1] = b01[i]; b[2] = b02[i];
        b[3] = b10[i]; b[4] = b11[i]; b[5] = b12[i];
        b[6] = b20[i]; b[7] = b21[i]; b[8] = b22[i];
        fv_i_s_mat3_mul9(a, b, o);
        o00[i] = o[0]; o01[i] = o[1]; o02[i] = o[2];
        o10[i] = o[3]; o11[i] = o[4]; o12[i] = o[5];
        o20[i] = o[6]; o21[i] = o[7]; o22[i] = o[8];
    }
}

static inline void fv_i_mat3_inverse_soa(                         const float *__restrict m00,
                         const float *__restrict m01,
                         const float *__restrict m02,
                         const float *__restrict m10,
                         const float *__restrict m11,
                         const float *__restrict m12,
                         const float *__restrict m20,
                         const float *__restrict m21,
                         const float *__restrict m22,
                         float *__restrict o00,
                         float *__restrict o01,
                         float *__restrict o02,
                         float *__restrict o10,
                         float *__restrict o11,
                         float *__restrict o12,
                         float *__restrict o20,
                         float *__restrict o21,
                         float *__restrict o22,
                         size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_mat3_inv_soa_avx2(m00, m01, m02, m10, m11, m12, m20, m21, m22, o00, o01, o02, o10, o11, o12, o20, o21, o22, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        float m[9];
        float o[9];
        m[0] = m00[i]; m[1] = m01[i]; m[2] = m02[i];
        m[3] = m10[i]; m[4] = m11[i]; m[5] = m12[i];
        m[6] = m20[i]; m[7] = m21[i]; m[8] = m22[i];
        (void)fv_i_s_mat3_inv9(m, o);
        o00[i] = o[0]; o01[i] = o[1]; o02[i] = o[2];
        o10[i] = o[3]; o11[i] = o[4]; o12[i] = o[5];
        o20[i] = o[6]; o21[i] = o[7]; o22[i] = o[8];
    }
}

/* ---------------- _fast approximate normalize: rsqrtps + one Newton-Raphson step (P7) ----------------
 * The exact normalize entries above stay as they are; these are separate names in a separate contract.
 * One Newton iteration squares the ~2^-12 rsqrtps error, so the relative error is about 1e-7 (1-2 ulp);
 * the scalar fallback uses 1/sqrtf and is exact. The zero guard is the same value as the exact entry
 * (zeros, or the identity for the quaternion), applied by a blend so no NaN escapes. No FMA is used
 * (the section is inside the fp-contract=off pragma), and the kernels only exist when FV_I_INTRINSICS. */

#if FV_I_INTRINSICS

/* vec2: rsqrtps + one Newton-Raphson step; the zero guard is a blend, not a branch. */
FV_I_TARGET_AVX2 static inline void fv_i_vec2_norm_fast_avx2(                                                        const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        float *__restrict ox,
                                                        float *__restrict oy,
                                                        size_t n)
{
    const __m256 eps2 = _mm256_set1_ps((float)FV_EPS_SQR);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 half = _mm256_set1_ps((float)0.5);
    const __m256 one_point_five = _mm256_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len2 = _mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y));
        __m256 r = _mm256_rsqrt_ps(len2);
        __m256 t = _mm256_mul_ps(_mm256_mul_ps(half, len2), _mm256_mul_ps(r, r));
        r = _mm256_mul_ps(r, _mm256_sub_ps(one_point_five, t));
        {
            __m256 guard = _mm256_cmp_ps(len2, eps2, _CMP_LE_OQ);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, r), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, r), zero, guard));
        }
    }
}

/* Fused madd_normalize fast: (a*s + b), then rsqrt+NR normalize with the same squared guard. The
 * scalar fallback calls the exact fused scalar helper, so accuracy is never worse. */
FV_I_TARGET_AVX2 static inline void fv_i_v2_madd_norm_fast_avx2(const float *__restrict ax,
                                                                const float *__restrict ay,
                                                                const float *__restrict bx,
                                                                const float *__restrict by,
                                                                fv_real s, float *__restrict ox,
                                                                float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    const __m256 eps2 = _mm256_set1_ps((float)FV_EPS_SQR);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 half = _mm256_set1_ps((float)0.5);
    const __m256 one_point_five = _mm256_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ax + i), vs), _mm256_loadu_ps(bx + i));
        __m256 y = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(ay + i), vs), _mm256_loadu_ps(by + i));
        __m256 len2 = _mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y));
        __m256 r = _mm256_rsqrt_ps(len2);
        __m256 t = _mm256_mul_ps(_mm256_mul_ps(half, len2), _mm256_mul_ps(r, r));
        r = _mm256_mul_ps(r, _mm256_sub_ps(one_point_five, t));
        {
            __m256 guard = _mm256_cmp_ps(len2, eps2, _CMP_LE_OQ);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, r), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, r), zero, guard));
        }
    }
}

#if FV_I_AVX512
FV_I_TARGET_AVX512 static inline void fv_i_v2_madd_norm_fast_avx512(const float *__restrict ax,
                                                                    const float *__restrict ay,
                                                                    const float *__restrict bx,
                                                                    const float *__restrict by,
                                                                    fv_real s,
                                                                    float *__restrict ox,
                                                                    float *__restrict oy, size_t n)
{
    const __m512 vs = _mm512_set1_ps(s);
    const __m512 eps2 = _mm512_set1_ps((float)FV_EPS_SQR);
    const __m512 zero = _mm512_setzero_ps();
    const __m512 half = _mm512_set1_ps((float)0.5);
    const __m512 one_point_five = _mm512_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512 x = _mm512_add_ps(_mm512_mul_ps(_mm512_loadu_ps(ax + i), vs), _mm512_loadu_ps(bx + i));
        __m512 y = _mm512_add_ps(_mm512_mul_ps(_mm512_loadu_ps(ay + i), vs), _mm512_loadu_ps(by + i));
        __m512 len2 = _mm512_add_ps(_mm512_mul_ps(x, x), _mm512_mul_ps(y, y));
        __m512 r = _mm512_rsqrt14_ps(len2);
        __m512 t = _mm512_mul_ps(_mm512_mul_ps(half, len2), _mm512_mul_ps(r, r));
        r = _mm512_mul_ps(r, _mm512_sub_ps(one_point_five, t));
        {
            __mmask16 guard = _mm512_cmp_ps_mask(len2, eps2, _CMP_LE_OQ);
            _mm512_storeu_ps(ox + i, _mm512_mask_blend_ps(guard, _mm512_mul_ps(x, r), zero));
            _mm512_storeu_ps(oy + i, _mm512_mask_blend_ps(guard, _mm512_mul_ps(y, r), zero));
        }
    }
}
#endif /* FV_I_AVX512 */

/* vec3: rsqrtps + one Newton-Raphson step; the zero guard is a blend, not a branch. */
FV_I_TARGET_AVX2 static inline void fv_i_vec3_norm_fast_avx2(                                                        const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        float *__restrict ox,
                                                        float *__restrict oy,
                                                        float *__restrict oz,
                                                        size_t n)
{
    const __m256 eps2 = _mm256_set1_ps((float)FV_EPS_SQR);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 half = _mm256_set1_ps((float)0.5);
    const __m256 one_point_five = _mm256_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 len2 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y)), _mm256_mul_ps(z, z));
        __m256 r = _mm256_rsqrt_ps(len2);
        __m256 t = _mm256_mul_ps(_mm256_mul_ps(half, len2), _mm256_mul_ps(r, r));
        r = _mm256_mul_ps(r, _mm256_sub_ps(one_point_five, t));
        {
            __m256 guard = _mm256_cmp_ps(len2, eps2, _CMP_LE_OQ);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, r), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, r), zero, guard));
            _mm256_storeu_ps(oz + i, _mm256_blendv_ps(_mm256_mul_ps(z, r), zero, guard));
        }
    }
}

/* vec4: rsqrtps + one Newton-Raphson step; the zero guard is a blend, not a branch. */
FV_I_TARGET_AVX2 static inline void fv_i_vec4_norm_fast_avx2(                                                        const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict aw,
                                                        float *__restrict ox,
                                                        float *__restrict oy,
                                                        float *__restrict oz,
                                                        float *__restrict ow,
                                                        size_t n)
{
    const __m256 eps2 = _mm256_set1_ps((float)FV_EPS_SQR);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 half = _mm256_set1_ps((float)0.5);
    const __m256 one_point_five = _mm256_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y)), _mm256_mul_ps(z, z)), _mm256_mul_ps(w, w));
        __m256 r = _mm256_rsqrt_ps(len2);
        __m256 t = _mm256_mul_ps(_mm256_mul_ps(half, len2), _mm256_mul_ps(r, r));
        r = _mm256_mul_ps(r, _mm256_sub_ps(one_point_five, t));
        {
            __m256 guard = _mm256_cmp_ps(len2, eps2, _CMP_LE_OQ);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, r), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, r), zero, guard));
            _mm256_storeu_ps(oz + i, _mm256_blendv_ps(_mm256_mul_ps(z, r), zero, guard));
            _mm256_storeu_ps(ow + i, _mm256_blendv_ps(_mm256_mul_ps(w, r), zero, guard));
        }
    }
}

/* quat: rsqrtps + one Newton-Raphson step; the zero guard is a blend, not a branch. */
FV_I_TARGET_AVX2 static inline void fv_i_quat_norm_fast_avx2(                                                        const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict aw,
                                                        float *__restrict ox,
                                                        float *__restrict oy,
                                                        float *__restrict oz,
                                                        float *__restrict ow,
                                                        size_t n)
{
    const __m256 eps2 = _mm256_set1_ps((float)FV_EPS_SQR);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 one = _mm256_set1_ps((float)1);
    const __m256 half = _mm256_set1_ps((float)0.5);
    const __m256 one_point_five = _mm256_set1_ps((float)1.5);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len2 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(x, x), _mm256_mul_ps(y, y)), _mm256_mul_ps(z, z)), _mm256_mul_ps(w, w));
        __m256 r = _mm256_rsqrt_ps(len2);
        __m256 t = _mm256_mul_ps(_mm256_mul_ps(half, len2), _mm256_mul_ps(r, r));
        r = _mm256_mul_ps(r, _mm256_sub_ps(one_point_five, t));
        {
            __m256 guard = _mm256_cmp_ps(len2, eps2, _CMP_LE_OQ);
            _mm256_storeu_ps(ox + i, _mm256_blendv_ps(_mm256_mul_ps(x, r), zero, guard));
            _mm256_storeu_ps(oy + i, _mm256_blendv_ps(_mm256_mul_ps(y, r), zero, guard));
            _mm256_storeu_ps(oz + i, _mm256_blendv_ps(_mm256_mul_ps(z, r), zero, guard));
            _mm256_storeu_ps(ow + i, _mm256_blendv_ps(_mm256_mul_ps(w, r), one, guard));
        }
    }
}

#endif /* FV_I_INTRINSICS: _fast normalize kernels */

static inline void fv_i_vec2_madd_normalize_soa_fast(const float *__restrict ax,
                                                     const float *__restrict ay,
                                                     const float *__restrict bx,
                                                     const float *__restrict by, fv_real s,
                                                     float *__restrict ox, float *__restrict oy,
                                                     size_t n)
{
    size_t i = 0;
#if FV_I_AVX512
    if (__builtin_cpu_supports("avx512f") != 0) {
        size_t bulk = n & ~(size_t)15;
        if (bulk != 0) {
            fv_i_v2_madd_norm_fast_avx512(ax, ay, bx, by, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
#if FV_I_INTRINSICS
    if (i == 0 && __builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v2_madd_norm_fast_avx2(ax, ay, bx, by, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_madd_norm(ax[i], ay[i], bx[i], by[i], s, &ox[i], &oy[i]);
    }
}


static inline void fv_i_vec2_normalize_soa_fast(                                                   const float *__restrict ax,
                                                   const float *__restrict ay,
                                                   float *__restrict ox,
                                                   float *__restrict oy,
                                                   size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_vec2_norm_fast_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_real len2 = ax[i] * ax[i] + ay[i] * ay[i];
        if (len2 <= FV_EPS_SQR) {
            ox[i] = (fv_real)0;
            oy[i] = (fv_real)0;
        } else {
            fv_real inv = (fv_real)1 / sqrtf(len2);
            ox[i] = ax[i] * inv;
            oy[i] = ay[i] * inv;
        }
    }
}

static inline void fv_i_vec3_normalize_soa_fast(                                                   const float *__restrict ax,
                                                   const float *__restrict ay,
                                                   const float *__restrict az,
                                                   float *__restrict ox,
                                                   float *__restrict oy,
                                                   float *__restrict oz,
                                                   size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_vec3_norm_fast_avx2(ax, ay, az, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_real len2 = (ax[i] * ax[i] + ay[i] * ay[i]) + az[i] * az[i];
        if (len2 <= FV_EPS_SQR) {
            ox[i] = (fv_real)0;
            oy[i] = (fv_real)0;
            oz[i] = (fv_real)0;
        } else {
            fv_real inv = (fv_real)1 / sqrtf(len2);
            ox[i] = ax[i] * inv;
            oy[i] = ay[i] * inv;
            oz[i] = az[i] * inv;
        }
    }
}

static inline void fv_i_vec4_normalize_soa_fast(                                                   const float *__restrict ax,
                                                   const float *__restrict ay,
                                                   const float *__restrict az,
                                                   const float *__restrict aw,
                                                   float *__restrict ox,
                                                   float *__restrict oy,
                                                   float *__restrict oz,
                                                   float *__restrict ow,
                                                   size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_vec4_norm_fast_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_real len2 = ((ax[i] * ax[i] + ay[i] * ay[i]) + az[i] * az[i]) + aw[i] * aw[i];
        if (len2 <= FV_EPS_SQR) {
            ox[i] = (fv_real)0;
            oy[i] = (fv_real)0;
            oz[i] = (fv_real)0;
            ow[i] = (fv_real)0;
        } else {
            fv_real inv = (fv_real)1 / sqrtf(len2);
            ox[i] = ax[i] * inv;
            oy[i] = ay[i] * inv;
            oz[i] = az[i] * inv;
            ow[i] = aw[i] * inv;
        }
    }
}

static inline void fv_i_quat_normalize_soa_fast(                                                   const float *__restrict ax,
                                                   const float *__restrict ay,
                                                   const float *__restrict az,
                                                   const float *__restrict aw,
                                                   float *__restrict ox,
                                                   float *__restrict oy,
                                                   float *__restrict oz,
                                                   float *__restrict ow,
                                                   size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_quat_norm_fast_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_real len2 = ((ax[i] * ax[i] + ay[i] * ay[i]) + az[i] * az[i]) + aw[i] * aw[i];
        if (len2 <= FV_EPS_SQR) {
            ox[i] = (fv_real)0;
            oy[i] = (fv_real)0;
            oz[i] = (fv_real)0;
            ow[i] = (fv_real)1;
        } else {
            fv_real inv = (fv_real)1 / sqrtf(len2);
            ox[i] = ax[i] * inv;
            oy[i] = ay[i] * inv;
            oz[i] = az[i] * inv;
            ow[i] = aw[i] * inv;
        }
    }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

/* ======== public names: static inline definitions in the default form ========
 * A body is one forwarding line with no floating-point expression, so all of them stay outside the pragma (the smaller the scope, the better);
 * forwarding carries no semantics -- every number is decided by fv_i_* -- so the two forms (these, and the extern definitions in the .cpp)
 * are two ways of wiring up one and the same implementation. */
#if !defined(FV_BUILD)

static inline void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                  fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4(a, b, o, n);
}


static inline void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                  fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4(a, b, o, n);
}


static inline void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                                    fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4(a, s, o, n);
}


static inline void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_real t, fv_vec2 *__restrict o, size_t n)
{
    fv_i_lerp_x4(a, b, t, o, n);
}


static inline void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_real s, fv_vec2 *__restrict o, size_t n)
{
    fv_i_madd_x4(a, b, s, o, n);
}


static inline void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                       fv_vec2 *__restrict o, const uint8_t *__restrict mask,
                                       size_t n)
{
    fv_i_add_x4_mask(a, b, o, mask, n);
}


static inline void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                      size_t m)
{
    fv_i_add_x4_idx(a, b, o, idx, m);
}


static inline void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                                        fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                        size_t m)
{
    fv_i_scale_x4_idx(a, s, o, idx, m);
}


static inline void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict bx, const float *__restrict by,
                                   float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa(ax, ay, bx, by, ox, oy, n);
}


static inline void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                     fv_real s, float *__restrict ox, float *__restrict oy,
                                     size_t n)
{
    fv_i_scale_soa(ax, ay, s, ox, oy, n);
}


static inline float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict bx, const float *__restrict by,
                                    size_t n)
{
    return fv_i_dot_soa(ax, ay, bx, by, n);
}


static inline void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   float *__restrict o, size_t n)
{
    fv_i_len_soa(ax, ay, o, n);
}


static inline void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                         float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_normalize_soa(ax, ay, ox, oy, n);
}


static inline void fv_vec2_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                         fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4_stream(a, b, o, n);
}


static inline void fv_vec2_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                         fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4_stream(a, b, o, n);
}


static inline void fv_vec2_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                                           fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4_stream(a, s, o, n);
}


static inline void fv_vec2_lerp_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                          fv_real t, fv_vec2 *__restrict o, size_t n)
{
    fv_i_lerp_x4_stream(a, b, t, o, n);
}


static inline void fv_vec2_madd_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                          fv_real s, fv_vec2 *__restrict o, size_t n)
{
    fv_i_madd_x4_stream(a, b, s, o, n);
}


static inline void fv_vec2_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                          const float *__restrict bx, const float *__restrict by,
                                          float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa_stream(ax, ay, bx, by, ox, oy, n);
}


static inline void fv_vec2_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                            fv_real s, float *__restrict ox, float *__restrict oy,
                                            size_t n)
{
    fv_i_scale_soa_stream(ax, ay, s, ox, oy, n);
}


static inline void fv_vec2_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                float *__restrict ox, float *__restrict oy,
                                                size_t n)
{
    fv_i_normalize_soa_stream(ax, ay, ox, oy, n);
}


static inline void fv_vec2_madd_normalize_soa(const float *__restrict ax,
                                         const float *__restrict ay,
                                         const float *__restrict bx,
                                         const float *__restrict by, fv_real s,
                                         float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_madd_norm_soa(ax, ay, bx, by, s, ox, oy, n);
}


static inline void fv_vec2_madd_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                const float *__restrict bx,
                                                const float *__restrict by, fv_real s,
                                                float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_madd_norm_soa_stream(ax, ay, bx, by, s, ox, oy, n);
}


static inline void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                               const uint8_t *__restrict mask, size_t n)
{
    fv_i_add_x4_mask_inplace(a, b, mask, n);
}

/* Public names: static inline forwarding in the default form (isomorphic to the 22 vec2 forwarders above).
 * A body is one forwarding line with no floating-point expression at all, so all of them stay outside the pragma. */
static inline void fv_vec3_add_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, const float *__restrict bx,
                                   const float *__restrict by, const float *__restrict bz,
                                   float *__restrict ox, float *__restrict oy,
                                   float *__restrict oz, size_t n)
{
    fv_i_v3_add_soa(ax, ay, az, bx, by, bz, ox, oy, oz, n);
}


static inline void fv_vec3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                     const float *__restrict az, fv_real s,
                                     float *__restrict ox, float *__restrict oy,
                                     float *__restrict oz, size_t n)
{
    fv_i_v3_scale_soa(ax, ay, az, s, ox, oy, oz, n);
}


static inline void fv_vec3_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                          const float *__restrict az,
                                          const float *__restrict bx, const float *__restrict by,
                                          const float *__restrict bz,
                                          float *__restrict ox, float *__restrict oy,
                                          float *__restrict oz, size_t n)
{
    fv_i_v3_add_soa_stream(ax, ay, az, bx, by, bz, ox, oy, oz, n);
}


static inline void fv_vec3_scale_soa_stream(const float *__restrict ax,
                                            const float *__restrict ay,
                                            const float *__restrict az, fv_real s,
                                            float *__restrict ox, float *__restrict oy,
                                            float *__restrict oz, size_t n)
{
    fv_i_v3_scale_soa_stream(ax, ay, az, s, ox, oy, oz, n);
}


static inline void fv_vec4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                  fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_add_x2(a, b, o, n);
}


static inline void fv_vec4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                  fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_sub_x2(a, b, o, n);
}


static inline void fv_vec4_scale_x2(const fv_vec4 *__restrict a, fv_real s,
                                    fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_scale_x2(a, s, o, n);
}


static inline void fv_vec4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                   fv_real t, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_lerp_x2(a, b, t, o, n);
}


static inline void fv_vec4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                   fv_real s, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_madd_x2(a, b, s, o, n);
}


static inline void fv_vec4_add_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                         fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_add_x2_stream(a, b, o, n);
}


static inline void fv_vec4_sub_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                         fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_sub_x2_stream(a, b, o, n);
}


static inline void fv_vec4_scale_x2_stream(const fv_vec4 *__restrict a, fv_real s,
                                           fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_scale_x2_stream(a, s, o, n);
}


static inline void fv_vec4_lerp_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                          fv_real t, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_lerp_x2_stream(a, b, t, o, n);
}


static inline void fv_vec4_madd_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                          fv_real s, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_madd_x2_stream(a, b, s, o, n);
}


static inline float fv_vec3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict az, const float *__restrict bx,
                                    const float *__restrict by, const float *__restrict bz,
                                    size_t n)
{
    return fv_i_v3_dot_soa(ax, ay, az, bx, by, bz, n);
}


static inline void fv_vec3_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, float *__restrict o, size_t n)
{
    fv_i_v3_len_soa(ax, ay, az, o, n);
}


static inline void fv_vec3_normalize_soa(const float *__restrict ax,
                                         const float *__restrict ay,
                                         const float *__restrict az, float *__restrict ox,
                                         float *__restrict oy, float *__restrict oz, size_t n)
{
    fv_i_v3_normalize_soa(ax, ay, az, ox, oy, oz, n);
}


static inline void fv_vec3_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                const float *__restrict az,
                                                float *__restrict ox, float *__restrict oy,
                                                float *__restrict oz, size_t n)
{
    fv_i_v3_normalize_soa_stream(ax, ay, az, ox, oy, oz, n);
}


static inline float fv_vec4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict az, const float *__restrict aw,
                                    const float *__restrict bx, const float *__restrict by,
                                    const float *__restrict bz, const float *__restrict bw,
                                    size_t n)
{
    return fv_i_v4_dot_soa(ax, ay, az, aw, bx, by, bz, bw, n);
}


static inline void fv_vec4_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, const float *__restrict aw,
                                   float *__restrict o, size_t n)
{
    fv_i_v4_len_soa(ax, ay, az, aw, o, n);
}


static inline void fv_vec4_normalize_soa(const float *__restrict ax,
                                         const float *__restrict ay,
                                         const float *__restrict az,
                                         const float *__restrict aw, float *__restrict ox,
                                         float *__restrict oy, float *__restrict oz,
                                         float *__restrict ow, size_t n)
{
    fv_i_v4_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
}


static inline void fv_vec4_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                const float *__restrict az,
                                                const float *__restrict aw,
                                                float *__restrict ox, float *__restrict oy,
                                                float *__restrict oz, float *__restrict ow,
                                                size_t n)
{
    fv_i_v4_normalize_soa_stream(ax, ay, az, aw, ox, oy, oz, ow, n);
}


static inline void fv_mat4_transform_points(const fv_mat4 *__restrict m,
                                           const fv_vec4 *__restrict in,
                                           fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_points(m, in, out, n);
}


static inline void fv_mat4_transform_dirs(const fv_mat4 *__restrict m,
                                         const fv_vec4 *__restrict in,
                                         fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_dirs(m, in, out, n);
}


static inline void fv_mat4_transform_points_stream(const fv_mat4 *__restrict m,
                                                   const fv_vec4 *__restrict in,
                                                   fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_points_stream(m, in, out, n);
}


static inline void fv_mat4_transform_dirs_stream(const fv_mat4 *__restrict m,
                                                 const fv_vec4 *__restrict in,
                                                 fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_dirs_stream(m, in, out, n);
}


static inline void fv_vec2_addmul_x4(fv_vec2 *__restrict x, const fv_vec2 *__restrict v,
                                     fv_real s, size_t n)
{
    fv_i_addmul_x4(x, v, s, n);
}


static inline void fv_vec2_addmul_soa(float *__restrict x, float *__restrict y,
                                     const float *__restrict vx, const float *__restrict vy,
                                     fv_real s, size_t n)
{
    fv_i_addmul_soa(x, y, vx, vy, s, n);
}


static inline void fv_vec3_addmul_soa(float *__restrict x, float *__restrict y,
                                     float *__restrict z, const float *__restrict vx,
                                     const float *__restrict vy, const float *__restrict vz,
                                     fv_real s, size_t n)
{
    fv_i_addmul_v3_soa(x, y, z, vx, vy, vz, s, n);
}


static inline void fv_vec4_addmul_x2(fv_vec4 *__restrict x, const fv_vec4 *__restrict v,
                                     fv_real s, size_t n)
{
    fv_i_v4_addmul_x2(x, v, s, n);
}


static inline void fv_vec2_maddv_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                    const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n)
{
    fv_i_maddv_x4(a, b, c, o, n);
}


static inline void fv_vec2_maddv_x4_stream(const fv_vec2 *__restrict a,
                                           const fv_vec2 *__restrict b,
                                           const fv_vec2 *__restrict c,
                                           fv_vec2 *__restrict o, size_t n)
{
    fv_i_maddv_x4_stream(a, b, c, o, n);
}


static inline void fv_quat_mul_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, const float *__restrict aw,
                                   const float *__restrict bx, const float *__restrict by,
                                   const float *__restrict bz, const float *__restrict bw,
                                   float *__restrict ox, float *__restrict oy,
                                   float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_quat_mul_soa(ax, ay, az, aw, bx, by, bz, bw, ox, oy, oz, ow, n);
}


static inline void fv_quat_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                         const float *__restrict az, const float *__restrict aw,
                                         float *__restrict ox, float *__restrict oy,
                                         float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_quat_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
}


static inline void fv_quat_rotate_vec3_soa(const float *__restrict qx,
                                           const float *__restrict qy,
                                           const float *__restrict qz,
                                           const float *__restrict qw,
                                           const float *__restrict vx,
                                           const float *__restrict vy,
                                           const float *__restrict vz,
                                           float *__restrict ox, float *__restrict oy,
                                           float *__restrict oz, size_t n)
{
    fv_i_quat_rotate_vec3_soa(qx, qy, qz, qw, vx, vy, vz, ox, oy, oz, n);
}


static inline void fv_mat3_transform_vec3_soa(const float *__restrict m00,
                                              const float *__restrict m01,
                                              const float *__restrict m02,
                                              const float *__restrict m10,
                                              const float *__restrict m11,
                                              const float *__restrict m12,
                                              const float *__restrict m20,
                                              const float *__restrict m21,
                                              const float *__restrict m22,
                                              const float *__restrict vx,
                                              const float *__restrict vy,
                                              const float *__restrict vz,
                                              float *__restrict ox, float *__restrict oy,
                                              float *__restrict oz, size_t n)
{
    fv_i_mat3_transform_vec3_soa(m00, m01, m02, m10, m11, m12, m20, m21, m22,
                                 vx, vy, vz, ox, oy, oz, n);
}

static inline void fv_mat3_mul_soa(                                   const float *__restrict a00,
                                   const float *__restrict a01,
                                   const float *__restrict a02,
                                   const float *__restrict a10,
                                   const float *__restrict a11,
                                   const float *__restrict a12,
                                   const float *__restrict a20,
                                   const float *__restrict a21,
                                   const float *__restrict a22,
                                   const float *__restrict b00,
                                   const float *__restrict b01,
                                   const float *__restrict b02,
                                   const float *__restrict b10,
                                   const float *__restrict b11,
                                   const float *__restrict b12,
                                   const float *__restrict b20,
                                   const float *__restrict b21,
                                   const float *__restrict b22,
                                   float *__restrict o00,
                                   float *__restrict o01,
                                   float *__restrict o02,
                                   float *__restrict o10,
                                   float *__restrict o11,
                                   float *__restrict o12,
                                   float *__restrict o20,
                                   float *__restrict o21,
                                   float *__restrict o22,
                                   size_t n)
{
    fv_i_mat3_mul_soa(a00, a01, a02, a10, a11, a12, a20, a21, a22, b00, b01, b02, b10, b11, b12, b20, b21, b22, o00, o01, o02, o10, o11, o12, o20, o21, o22, n);
}


static inline void fv_mat3_inverse_soa(                                       const float *__restrict m00,
                                       const float *__restrict m01,
                                       const float *__restrict m02,
                                       const float *__restrict m10,
                                       const float *__restrict m11,
                                       const float *__restrict m12,
                                       const float *__restrict m20,
                                       const float *__restrict m21,
                                       const float *__restrict m22,
                                       float *__restrict o00,
                                       float *__restrict o01,
                                       float *__restrict o02,
                                       float *__restrict o10,
                                       float *__restrict o11,
                                       float *__restrict o12,
                                       float *__restrict o20,
                                       float *__restrict o21,
                                       float *__restrict o22,
                                       size_t n)
{
    fv_i_mat3_inverse_soa(m00, m01, m02, m10, m11, m12, m20, m21, m22, o00, o01, o02, o10, o11, o12, o20, o21, o22, n);
}

static inline void fv_vec2_normalize_soa_fast(                                const float *__restrict ax,
                                const float *__restrict ay,
                                float *__restrict ox,
                                float *__restrict oy,
                                size_t n)
{
    fv_i_vec2_normalize_soa_fast(ax, ay, ox, oy, n);
}

static inline void fv_vec3_normalize_soa_fast(                                const float *__restrict ax,
                                const float *__restrict ay,
                                const float *__restrict az,
                                float *__restrict ox,
                                float *__restrict oy,
                                float *__restrict oz,
                                size_t n)
{
    fv_i_vec3_normalize_soa_fast(ax, ay, az, ox, oy, oz, n);
}

static inline void fv_vec4_normalize_soa_fast(                                const float *__restrict ax,
                                const float *__restrict ay,
                                const float *__restrict az,
                                const float *__restrict aw,
                                float *__restrict ox,
                                float *__restrict oy,
                                float *__restrict oz,
                                float *__restrict ow,
                                size_t n)
{
    fv_i_vec4_normalize_soa_fast(ax, ay, az, aw, ox, oy, oz, ow, n);
}

static inline void fv_quat_normalize_soa_fast(                                const float *__restrict ax,
                                const float *__restrict ay,
                                const float *__restrict az,
                                const float *__restrict aw,
                                float *__restrict ox,
                                float *__restrict oy,
                                float *__restrict oz,
                                float *__restrict ow,
                                size_t n)
{
    fv_i_quat_normalize_soa_fast(ax, ay, az, aw, ox, oy, oz, ow, n);
}

static inline void fv_vec2_madd_normalize_soa_fast(const float *__restrict ax,
                                                   const float *__restrict ay,
                                                   const float *__restrict bx,
                                                   const float *__restrict by, fv_real s,
                                                   float *__restrict ox, float *__restrict oy,
                                                   size_t n)
{
    fv_i_vec2_madd_normalize_soa_fast(ax, ay, bx, by, s, ox, oy, n);
}

#endif /* !FV_BUILD: public names of the default form */

#endif /* !FV_USE_EXTERN: implementation section */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FLASHVEC_H */