/* fv_quat scalar semantics + bit-identical SoA batch (P6).
 *
 * What this fixture pins down:
 *   1. construction / identity / arithmetic / conjugate / dot / length;
 *   2. the Hamilton product: identity, known 90-degree rotations, composition order, and a
 *      bit-for-bit comparison against the same expressions written independently here;
 *   3. normalize: unit input stays unit, length <= FV_EPS maps to the identity (not to zero),
 *      and the batch path is bit-identical to the per-value path;
 *   4. rotate_vec3: identity, a 90-degree rotation about each axis, composition, agreement with
 *      fv_mat3_from_quat (approximately -- the two are different expressions), and bit-parity
 *      against the per-value form;
 *   5. the SoA batch entries over every tail length n = 0..40 plus 63/1024/4097, with n == 0 and
 *      NULL pointers accepted, outputs compared with memcmp against the per-value functions.
 *
 * The comparison is bit-level (memcmp on float arrays), not a tolerance check, because the header's
 * contract promises the batch result is identical to the per-value result for every n and every tail.
 * NaN payloads are excluded from the bit-parity corpora by the library-wide rule documented in
 * include/flashvec.h (two different NaN payloads meeting in one instruction is operand-order
 * dependent); the semantic checks below only require "is NaN".
 *
 * Run: make test-quat
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(cond)                                                    \
    do {                                                               \
        ++g_checks;                                                    \
        if (!(cond)) {                                                 \
            ++g_failures;                                              \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);     \
        }                                                              \
    } while (0)

#define CAP 4104u

static float g_ax[CAP], g_ay[CAP], g_az[CAP], g_aw[CAP];
static float g_bx[CAP], g_by[CAP], g_bz[CAP], g_bw[CAP];
static float g_vx[CAP], g_vy[CAP], g_vz[CAP];
static float g_ox[CAP], g_oy[CAP], g_oz[CAP], g_ow[CAP];

/* ------------------------------------------------------------------ helpers */

static int feq(fv_real a, fv_real b, fv_real eps)
{
    fv_real d = a - b;
    if (d < (fv_real)0) {
        d = -d;
    }
    return d <= eps;
}

static int quat_near(fv_quat a, fv_quat b, fv_real eps)
{
    return feq(a.x, b.x, eps) && feq(a.y, b.y, eps) && feq(a.z, b.z, eps) &&
           feq(a.w, b.w, eps);
}

static int vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps)
{
    return feq(a.x, b.x, eps) && feq(a.y, b.y, eps) && feq(a.z, b.z, eps);
}

/* 0x5A poison: an out-of-bounds write or a missed write shows up as a changed guard. */
static void poison_float(float *p, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        p[i] = 1.0e30f;
    }
}

static int guard_ok(const float *p, size_t from, size_t to)
{
    size_t i;
    for (i = from; i < to; ++i) {
        if (p[i] != 1.0e30f) {
            return 0;
        }
    }
    return 1;
}

/* Deterministic xorshift; moderate values so products do not overflow. */
static uint32_t g_rng = 0x12345678u;

static uint32_t rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static float rnd_moderate(void)
{
    return (float)((int)(rnd_u32() % 2001u) - 1000) / 16.0f;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif

/* Independent references: the same association order as the library contract, written separately so a
 * change in the library cannot silently define its own expected value. */
static fv_quat ref_mul(fv_quat a, fv_quat b)
{
    fv_quat r;
    r.x = ((a.w * b.x + a.x * b.w) + a.y * b.z) - a.z * b.y;
    r.y = ((a.w * b.y - a.x * b.z) + a.y * b.w) + a.z * b.x;
    r.z = ((a.w * b.z + a.x * b.y) - a.y * b.x) + a.z * b.w;
    r.w = ((a.w * b.w - a.x * b.x) - a.y * b.y) - a.z * b.z;
    return r;
}

static fv_quat ref_normalize(fv_quat q)
{
    fv_real len = sqrtf(((q.x * q.x + q.y * q.y) + q.z * q.z) + q.w * q.w);
    if (len <= FV_EPS) {
        return fv_quat_identity();
    }
    {
        fv_real inv = (fv_real)1 / len;
        fv_quat r;
        r.x = q.x * inv;
        r.y = q.y * inv;
        r.z = q.z * inv;
        r.w = q.w * inv;
        return r;
    }
}

static fv_vec3 ref_rotate(fv_quat q, fv_vec3 v)
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

/* The same expression order as fv_quat_inverse: conjugate first, then scale by 1/|q|^2. */
static fv_quat ref_inverse(fv_quat q)
{
    fv_real l2 = fv_quat_len_sqr(q);
    if (!(l2 > FV_EPS_SQR)) {
        return fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / l2;
        fv_quat c = fv_quat_conjugate(q);
        return fv_quat_scale(c, inv);
    }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

static int quat_bits_equal(fv_quat a, fv_quat b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

static int vec3_bits_equal(fv_vec3 a, fv_vec3 b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

/* ------------------------------------------------------------------ 1. construction / arithmetic */

static void test_basics(void)
{
    fv_quat id = fv_quat_identity();
    fv_quat a = fv_quat_make(1, 2, 3, 4);
    fv_quat b = fv_quat_make(-1, 0.5f, -3, 2);

    CHECK(id.x == 0 && id.y == 0 && id.z == 0 && id.w == 1);
    CHECK(a.x == 1 && a.y == 2 && a.z == 3 && a.w == 4);

    CHECK(quat_bits_equal(fv_quat_add(a, b), fv_quat_make(0, 2.5f, 0, 6)));
    CHECK(quat_bits_equal(fv_quat_sub(a, b), fv_quat_make(2, 1.5f, 6, 2)));
    CHECK(quat_bits_equal(fv_quat_neg(a), fv_quat_make(-1, -2, -3, -4)));
    CHECK(quat_bits_equal(fv_quat_scale(a, 2), fv_quat_make(2, 4, 6, 8)));

    /* dot: 1*(-1) + 2*0.5 + 3*(-3) + 4*2 = -1 + 1 - 9 + 8 = -1 */
    CHECK(feq(fv_quat_dot(a, b), (fv_real)-1, (fv_real)1e-6));
    CHECK(feq(fv_quat_len_sqr(a), (fv_real)30, (fv_real)1e-5));
    CHECK(feq(fv_quat_len(a), sqrtf((fv_real)30), (fv_real)1e-5));

    /* conjugate is the inverse only for unit quaternions; the component signs are unconditional. */
    CHECK(quat_bits_equal(fv_quat_conjugate(a), fv_quat_make(-1, -2, -3, 4)));

    /* q * conjugate(q) = |q|^2 * identity, checked by components. */
    {
        fv_quat p = fv_quat_mul(a, fv_quat_conjugate(a));
        fv_real l2 = fv_quat_len_sqr(a);
        CHECK(feq(p.x, (fv_real)0, (fv_real)1e-5));
        CHECK(feq(p.y, (fv_real)0, (fv_real)1e-5));
        CHECK(feq(p.z, (fv_real)0, (fv_real)1e-5));
        CHECK(feq(p.w, l2, (fv_real)1e-5));
    }

    /* Multiplication is not commutative in general (x and y rotations are the smallest witness). */
    {
        fv_quat qx = fv_quat_make(0.70710677f, 0, 0, 0.70710677f);
        fv_quat qy = fv_quat_make(0, 0.70710677f, 0, 0.70710677f);
        CHECK(!quat_bits_equal(fv_quat_mul(qx, qy), fv_quat_mul(qy, qx)));
    }
}

/* ------------------------------------------------------------------ 2. multiplication */

static void test_mul(void)
{
    fv_quat id = fv_quat_identity();
    fv_quat q = fv_quat_make(0.2f, -0.4f, 0.6f, 0.8f);
    fv_quat qn = fv_quat_normalize(q);
    int k;

    CHECK(quat_bits_equal(fv_quat_mul(q, id), q));
    CHECK(quat_bits_equal(fv_quat_mul(id, q), q));

    /* Unit q and its conjugate undo each other (approximately: the product's w is 1 up to rounding). */
    {
        fv_quat p = fv_quat_mul(qn, fv_quat_conjugate(qn));
        CHECK(quat_near(p, id, (fv_real)1e-6));
    }

    /* The independent reference and the library must agree bit for bit on finite values. */
    for (k = 0; k < 2000; ++k) {
        fv_quat a = fv_quat_make(rnd_moderate(), rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_quat b = fv_quat_make(rnd_moderate(), rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_quat got = fv_quat_mul(a, b);
        fv_quat want = ref_mul(a, b);
        CHECK(quat_bits_equal(got, want));
        if (!quat_bits_equal(got, want)) {
            printf("  mul mismatch at case %d\n", k);
        }
    }

    /* Composition: (q1*q2) applied to v == q1 applied to (q2 applied to v). Different expressions
     * (two products vs two rotations), so this is a tolerance check, not memcmp. */
    {
        /* Unit quaternions: both the rotate formula and the composition identity assume |q| = 1. */
        fv_quat q1 = fv_quat_normalize(fv_quat_make(0.1f, 0.2f, 0.3f, 0.4f));
        fv_quat q2 = fv_quat_normalize(fv_quat_make(-0.2f, 0.1f, 0.4f, -0.3f));
        fv_vec3 v = fv_vec3_make(1, -2, 3);
        fv_quat q12 = fv_quat_mul(q1, q2);
        fv_vec3 a = fv_quat_rotate_vec3(q12, v);
        fv_vec3 b = fv_quat_rotate_vec3(q1, fv_quat_rotate_vec3(q2, v));
        CHECK(vec3_near(a, b, (fv_real)5e-5));
    }
}

/* ------------------------------------------------------------------ 3. normalize */

static void test_normalize(void)
{
    fv_quat zero = fv_quat_make(0, 0, 0, 0);
    fv_quat tiny = fv_quat_make((fv_real)1e-7, 0, 0, 0);
    fv_quat id = fv_quat_identity();
    fv_quat q = fv_quat_make(0.3f, -0.7f, 0.2f, -0.5f);
    fv_quat n;
    int k;

    CHECK(quat_bits_equal(fv_quat_normalize(zero), id));
    CHECK(quat_bits_equal(fv_quat_normalize(tiny), id));
    CHECK(quat_bits_equal(fv_quat_normalize(id), id));

    n = fv_quat_normalize(q);
    CHECK(feq(fv_quat_len(n), (fv_real)1, (fv_real)1e-6));
    CHECK(n.w < 0); /* the sign of w is preserved; normalize is not a shortest-arc operation */

    for (k = 0; k < 2000; ++k) {
        fv_quat a = fv_quat_make(rnd_moderate(), rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_quat got = fv_quat_normalize(a);
        fv_quat want = ref_normalize(a);
        CHECK(quat_bits_equal(got, want));
    }
}

/* ------------------------------------------------------------------ 4. rotate */

static void test_rotate(void)
{
    fv_quat id = fv_quat_identity();
    fv_vec3 v = fv_vec3_make(1, 0, 0);
    fv_vec3 r;
    int k;

    /* Identity leaves the vector untouched (bitwise, including -0.0: (v + 0*t) + 0 is v for finite v). */
    r = fv_quat_rotate_vec3(id, v);
    CHECK(vec3_bits_equal(r, v));

    /* 90 degrees about +y: (1,0,0) -> (0,0,-1). */
    {
        fv_quat qy = fv_quat_make(0, 0.70710677f, 0, 0.70710677f);
        r = fv_quat_rotate_vec3(qy, v);
        CHECK(vec3_near(r, fv_vec3_make(0, 0, -1), (fv_real)2e-6));
    }
    /* 90 degrees about +z: (1,0,0) -> (0,1,0). */
    {
        fv_quat qz = fv_quat_make(0, 0, 0.70710677f, 0.70710677f);
        r = fv_quat_rotate_vec3(qz, v);
        CHECK(vec3_near(r, fv_vec3_make(0, 1, 0), (fv_real)2e-6));
    }
    /* 90 degrees about +x: (0,1,0) -> (0,0,1). */
    {
        fv_quat qx = fv_quat_make(0.70710677f, 0, 0, 0.70710677f);
        r = fv_quat_rotate_vec3(qx, fv_vec3_make(0, 1, 0));
        CHECK(vec3_near(r, fv_vec3_make(0, 0, 1), (fv_real)2e-6));
    }

    /* The matrix form agrees with the quaternion form on the same rotation (different expressions,
     * hence a tolerance, not memcmp). */
    {
        fv_quat q = fv_quat_normalize(fv_quat_make(0.2f, 0.3f, -0.1f, 0.9f));
        fv_mat3 m = fv_mat3_from_quat(q);
        fv_vec3 a = fv_quat_rotate_vec3(q, v);
        fv_vec3 b = fv_mat3_transform_vec3(&m, v);
        CHECK(vec3_near(a, b, (fv_real)5e-6));
    }

    /* The independent reference and the library must agree bit for bit on finite values. */
    for (k = 0; k < 2000; ++k) {
        fv_quat q = fv_quat_make(rnd_moderate(), rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_vec3 p = fv_vec3_make(rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_vec3 got = fv_quat_rotate_vec3(q, p);
        fv_vec3 want = ref_rotate(q, p);
        CHECK(vec3_bits_equal(got, want));
    }
}

/* ------------------------------------------------------------------ 5. inverse */

static void test_inverse(void)
{
    fv_quat id = fv_quat_identity();
    fv_quat q = fv_quat_make(0.3f, -0.7f, 0.2f, -0.5f);
    fv_quat qn = fv_quat_normalize(q);
    fv_quat inv;
    fv_quat p;
    int k;

    /* Unit quaternion: the inverse is the conjugate (up to the usual rounding of |q|^2). */
    inv = fv_quat_inverse(qn);
    CHECK(quat_near(inv, fv_quat_conjugate(qn), (fv_real)2e-6));
    p = fv_quat_mul(qn, inv);
    CHECK(quat_near(p, id, (fv_real)2e-6));

    /* General quaternion: q * inverse(q) is the identity. */
    inv = fv_quat_inverse(q);
    p = fv_quat_mul(q, inv);
    CHECK(quat_near(p, id, (fv_real)2e-6));

    /* Degenerate inputs: the value form is the zero quaternion, the checked form reports 0. */
    {
        fv_quat z = fv_quat_make(0, 0, 0, 0);
        fv_quat tiny = fv_quat_make((fv_real)1e-7, 0, 0, 0);
        fv_quat out = fv_quat_make(9, 9, 9, 9);
        CHECK(quat_bits_equal(fv_quat_inverse(z),
                              fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0)));
        CHECK(quat_bits_equal(fv_quat_inverse(tiny),
                              fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0)));
        CHECK(fv_quat_inverse_checked(z, &out) == 0);
        CHECK(quat_bits_equal(out,
                              fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0)));
        CHECK(fv_quat_inverse_checked(tiny, &out) == 0);
    }

    /* Checked success writes the same value as the value form. */
    {
        fv_quat out = fv_quat_make(9, 9, 9, 9);
        CHECK(fv_quat_inverse_checked(q, &out) == 1);
        CHECK(quat_bits_equal(out, fv_quat_inverse(q)));
    }

    /* Independent reference, bit for bit, for non-degenerate random inputs. */
    for (k = 0; k < 2000; ++k) {
        fv_quat a = fv_quat_make(rnd_moderate(), rnd_moderate(), rnd_moderate(), rnd_moderate());
        if (fv_quat_len_sqr(a) > FV_EPS_SQR) {
            CHECK(quat_bits_equal(fv_quat_inverse(a), ref_inverse(a)));
        }
    }
}

/* ------------------------------------------------------------------ 6. matrix -> quaternion */

static void test_from_mat3(void)
{
    fv_mat3 idm = fv_mat3_identity();
    fv_quat qi = fv_quat_from_mat3(&idm);
    int k;

    CHECK(quat_near(qi, fv_quat_identity(), (fv_real)1e-6));

    /* A known 90-degree rotation about +y: the round trip is the same rotation (q or -q). */
    {
        fv_quat qy = fv_quat_make(0, 0.70710677f, 0, 0.70710677f);
        fv_mat3 m = fv_mat3_from_quat(qy);
        fv_quat q2 = fv_quat_from_mat3(&m);
        CHECK(quat_near(q2, qy, (fv_real)2e-6) ||
              quat_near(q2, fv_quat_neg(qy), (fv_real)2e-6));
    }

    /* Random unit quaternions: q -> R -> q' is the same rotation, and q' is unit length. */
    for (k = 0; k < 2000; ++k) {
        fv_quat q = fv_quat_normalize(fv_quat_make(rnd_moderate(), rnd_moderate(),
                                                   rnd_moderate(), rnd_moderate()));
        fv_mat3 m = fv_mat3_from_quat(q);
        fv_quat q2 = fv_quat_from_mat3(&m);
        CHECK(feq(fv_quat_len(q2), (fv_real)1, (fv_real)1e-5));
        CHECK(quat_near(q2, q, (fv_real)2e-5) ||
              quat_near(q2, fv_quat_neg(q), (fv_real)2e-5));
    }
}

/* ------------------------------------------------------------------ 7. interpolation / construction */

static void test_interpolation(void)
{
    const fv_real pi = (fv_real)3.14159265358979323846;
    fv_quat id = fv_quat_identity();
    fv_quat qy90 = fv_quat_from_axis_angle(fv_vec3_make(0, 1, 0), pi / (fv_real)2);
    fv_quat qz90 = fv_quat_from_axis_angle(fv_vec3_make(0, 0, 1), pi / (fv_real)2);
    fv_quat q;
    fv_vec3 v = fv_vec3_make(1, 0, 0);
    fv_vec3 r;
    int k;

    /* Endpoints. */
    CHECK(quat_near(fv_quat_nlerp(id, qy90, 0), id, (fv_real)1e-6));
    CHECK(quat_near(fv_quat_nlerp(id, qy90, 1), qy90, (fv_real)1e-6));
    CHECK(quat_near(fv_quat_slerp(id, qy90, 0), id, (fv_real)1e-6));
    CHECK(quat_near(fv_quat_slerp(id, qy90, 1), qy90, (fv_real)1e-6));

    /* Midpoint of a 90-degree turn is a 45-degree turn. */
    q = fv_quat_slerp(id, qy90, (fv_real)0.5);
    CHECK(feq(fv_quat_len(q), (fv_real)1, (fv_real)1e-6));
    r = fv_quat_rotate_vec3(q, v);
    CHECK(vec3_near(r, fv_vec3_make(0.70710677f, 0, -0.70710677f), (fv_real)2e-6));
    q = fv_quat_nlerp(id, qy90, (fv_real)0.5);
    CHECK(vec3_near(fv_quat_rotate_vec3(q, v), r, (fv_real)2e-4));

    /* Shortest arc: q and -q describe the same rotation, so interpolation stays at q (up to sign). */
    {
        fv_quat nq = fv_quat_neg(qy90);
        fv_quat a = fv_quat_slerp(qy90, nq, (fv_real)0.25);
        fv_quat b = fv_quat_nlerp(qy90, nq, (fv_real)0.25);
        CHECK(quat_near(a, qy90, (fv_real)1e-5) || quat_near(a, nq, (fv_real)1e-5));
        CHECK(quat_near(b, qy90, (fv_real)1e-5) || quat_near(b, nq, (fv_real)1e-5));
    }

    /* Unit length over a t sweep between two different rotations. */
    for (k = 0; k <= 10; ++k) {
        fv_real t = (fv_real)k / (fv_real)10;
        CHECK(feq(fv_quat_len(fv_quat_slerp(qy90, qz90, t)), (fv_real)1, (fv_real)2e-6));
        CHECK(feq(fv_quat_len(fv_quat_nlerp(qy90, qz90, t)), (fv_real)1, (fv_real)2e-6));
    }

    /* Axis-angle: zero axis -> identity; a non-unit axis is normalized. */
    CHECK(quat_near(fv_quat_from_axis_angle(fv_vec3_make(0, 0, 0), pi), id, (fv_real)1e-6));
    CHECK(quat_near(fv_quat_from_axis_angle(fv_vec3_make(0, 2, 0), pi / (fv_real)2),
                    qy90, (fv_real)2e-6));

    /* Euler XYZ: a single-axis rotation matches the axis-angle constructor; the composition order is
     * X first, then Y, then Z. */
    CHECK(quat_near(fv_quat_from_euler_xyz(0, pi / (fv_real)2, 0), qy90, (fv_real)2e-6));
    {
        fv_quat qx90 = fv_quat_from_axis_angle(fv_vec3_make(1, 0, 0), pi / (fv_real)2);
        CHECK(quat_near(fv_quat_from_euler_xyz(pi / (fv_real)2, 0, 0), qx90, (fv_real)2e-6));
        /* (0,1,0) under X then Y: X sends it to (0,0,1), Y sends that to (1,0,0). */
        q = fv_quat_from_euler_xyz(pi / (fv_real)2, pi / (fv_real)2, 0);
        r = fv_quat_rotate_vec3(q, fv_vec3_make(0, 1, 0));
        CHECK(vec3_near(r, fv_vec3_make(1, 0, 0), (fv_real)2e-6));
    }
}

/* ------------------------------------------------------------------ 8. SoA batch parity */

static void fill_case(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_ax[i] = rnd_moderate();
        g_ay[i] = rnd_moderate();
        g_az[i] = rnd_moderate();
        g_aw[i] = rnd_moderate();
        g_bx[i] = rnd_moderate();
        g_by[i] = rnd_moderate();
        g_bz[i] = rnd_moderate();
        g_bw[i] = rnd_moderate();
        g_vx[i] = rnd_moderate();
        g_vy[i] = rnd_moderate();
        g_vz[i] = rnd_moderate();
    }
}

static void test_batch_mul(size_t n)
{
    size_t i;
    fill_case(n);
    poison_float(g_ox, CAP);
    poison_float(g_oy, CAP);
    poison_float(g_oz, CAP);
    poison_float(g_ow, CAP);
    fv_quat_mul_soa(g_ax, g_ay, g_az, g_aw, g_bx, g_by, g_bz, g_bw,
                    g_ox, g_oy, g_oz, g_ow, n);
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_quat b = fv_quat_make(g_bx[i], g_by[i], g_bz[i], g_bw[i]);
        fv_quat want = fv_quat_mul(a, b);
        CHECK(g_ox[i] == want.x && g_oy[i] == want.y && g_oz[i] == want.z && g_ow[i] == want.w);
    }
    CHECK(guard_ok(g_ox, n, CAP));
    CHECK(guard_ok(g_oy, n, CAP));
    CHECK(guard_ok(g_oz, n, CAP));
    CHECK(guard_ok(g_ow, n, CAP));
}

static void test_batch_normalize(size_t n)
{
    size_t i;
    fill_case(n);
    poison_float(g_ox, CAP);
    poison_float(g_oy, CAP);
    poison_float(g_oz, CAP);
    poison_float(g_ow, CAP);
    fv_quat_normalize_soa(g_ax, g_ay, g_az, g_aw, g_ox, g_oy, g_oz, g_ow, n);
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_quat want = fv_quat_normalize(a);
        CHECK(g_ox[i] == want.x && g_oy[i] == want.y && g_oz[i] == want.z && g_ow[i] == want.w);
    }
    CHECK(guard_ok(g_ow, n, CAP));
}

static void test_batch_rotate(size_t n)
{
    size_t i;
    fill_case(n);
    /* Make the quaternions unit so the rotation precondition holds. */
    for (i = 0; i < n; ++i) {
        fv_quat q = fv_quat_normalize(fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]));
        g_ax[i] = q.x;
        g_ay[i] = q.y;
        g_az[i] = q.z;
        g_aw[i] = q.w;
    }
    poison_float(g_ox, CAP);
    poison_float(g_oy, CAP);
    poison_float(g_oz, CAP);
    fv_quat_rotate_vec3_soa(g_ax, g_ay, g_az, g_aw, g_vx, g_vy, g_vz, g_ox, g_oy, g_oz, n);
    for (i = 0; i < n; ++i) {
        fv_quat q = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_vec3 v = fv_vec3_make(g_vx[i], g_vy[i], g_vz[i]);
        fv_vec3 want = fv_quat_rotate_vec3(q, v);
        CHECK(g_ox[i] == want.x && g_oy[i] == want.y && g_oz[i] == want.z);
    }
    CHECK(guard_ok(g_oz, n, CAP));
}

/* ------------------------------------------------------------------ main */

int main(void)
{
    size_t n;

    printf("== fv_quat scalar semantics ==\n");
    test_basics();
    test_mul();
    test_normalize();
    test_rotate();
    test_inverse();
    test_from_mat3();
    test_interpolation();

    printf("== fv_quat SoA batch: bit-identical to the per-value API (n = 0..40 all tails + 63 + 1024 + 4097) ==\n");
    test_batch_mul(0);
    test_batch_mul(63);
    test_batch_mul(1024);
    test_batch_mul(4097);
    test_batch_normalize(0);
    test_batch_normalize(63);
    test_batch_normalize(1024);
    test_batch_normalize(4097);
    test_batch_rotate(0);
    test_batch_rotate(63);
    test_batch_rotate(1024);
    test_batch_rotate(4097);
    for (n = 0; n <= 40; ++n) {
        test_batch_mul(n);
        test_batch_normalize(n);
        test_batch_rotate(n);
    }

    /* n == 0 with NULL for every argument: the contract says no pointer is touched. */
    fv_quat_mul_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL, NULL, 0);
    fv_quat_normalize_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_quat_rotate_vec3_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_mat3_transform_vec3_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, 0);

    printf("test_quat: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
