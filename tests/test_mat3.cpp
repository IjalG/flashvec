/* fv_mat3 scalar semantics + bit-identical SoA batch (P6).
 *
 * Coverage:
 *   - identity / from_rows / transform_vec3 / mul / transpose: semantics and simple known values;
 *   - mul is associative in the composition sense and matches an independent triple-loop reference;
 *   - transpose twice is the original (bit for bit -- it only moves components);
 *   - from_quat is the matrix form of the same rotation (tolerance check against fv_quat_rotate_vec3);
 *   - the SoA batch entry is bit-identical to the per-value function for n = 0..40, 63, 1024, 4097,
 *     with n == 0 and NULL pointers accepted and no guard-byte damage.
 *
 * Run: make test-mat3
 */
#include "flashvec.h"

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

static float g_m00[CAP], g_m01[CAP], g_m02[CAP];
static float g_m10[CAP], g_m11[CAP], g_m12[CAP];
static float g_m20[CAP], g_m21[CAP], g_m22[CAP];
static float g_vx[CAP], g_vy[CAP], g_vz[CAP];
static float g_ox[CAP], g_oy[CAP], g_oz[CAP];

/* 9 x CAP component planes for the batch mul / inverse. */
static float g_a[9][CAP];
static float g_b[9][CAP];
static float g_om[9][CAP];
static float g_oi[9][CAP];

static int feq(fv_real a, fv_real b, fv_real eps)
{
    fv_real d = a - b;
    if (d < (fv_real)0) {
        d = -d;
    }
    return d <= eps;
}

static int vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps)
{
    return feq(a.x, b.x, eps) && feq(a.y, b.y, eps) && feq(a.z, b.z, eps);
}

static int vec3_bits_equal(fv_vec3 a, fv_vec3 b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

static int mat3_bits_equal(fv_mat3 a, fv_mat3 b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

static int same_bits_float(fv_real a, fv_real b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

/* Component equality that treats +0 and -0 as equal: inverse involves a multiply by 1/det, whose
 * sign can turn a zero cofactor into -0; that is not a numerical difference. */
static int mat3_component_equal(fv_mat3 a, fv_mat3 b)
{
    int k;
    for (k = 0; k < 9; ++k) {
        if (a.m[k] != b.m[k]) {
            return 0;
        }
    }
    return 1;
}

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

static uint32_t g_rng = 0x9E3779B9u;

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

/* Independent references, same association order as the library contract. */
static fv_vec3 ref_tf(const fv_mat3 *m, fv_vec3 v)
{
    fv_vec3 r;
    r.x = (m->m[0] * v.x + m->m[1] * v.y) + m->m[2] * v.z;
    r.y = (m->m[3] * v.x + m->m[4] * v.y) + m->m[5] * v.z;
    r.z = (m->m[6] * v.x + m->m[7] * v.y) + m->m[8] * v.z;
    return r;
}

static fv_mat3 ref_mul(const fv_mat3 *a, const fv_mat3 *b)
{
    fv_mat3 r;
    int i, j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r.m[i * 3 + j] = (a->m[i * 3 + 0] * b->m[0 * 3 + j] +
                              a->m[i * 3 + 1] * b->m[1 * 3 + j]) +
                             a->m[i * 3 + 2] * b->m[2 * 3 + j];
        }
    }
    return r;
}

static fv_real ref_det(const fv_mat3 *m)
{
    return (m->m[0] * (m->m[4] * m->m[8] - m->m[5] * m->m[7]) -
            m->m[1] * (m->m[3] * m->m[8] - m->m[5] * m->m[6])) +
           m->m[2] * (m->m[3] * m->m[7] - m->m[4] * m->m[6]);
}

/* Same cofactor/determinant expression order as fv_mat3_inverse; no threshold (the caller only
 * compares it on matrices the library accepted). */
static fv_mat3 ref_inverse(const fv_mat3 *m)
{
    fv_real det = ref_det(m);
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
    return r;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

static fv_mat3 make_mat(fv_real m00, fv_real m01, fv_real m02,
                        fv_real m10, fv_real m11, fv_real m12,
                        fv_real m20, fv_real m21, fv_real m22)
{
    fv_mat3 m;
    m.m[0] = m00; m.m[1] = m01; m.m[2] = m02;
    m.m[3] = m10; m.m[4] = m11; m.m[5] = m12;
    m.m[6] = m20; m.m[7] = m21; m.m[8] = m22;
    return m;
}

static void test_basics(void)
{
    fv_mat3 id = fv_mat3_identity();
    fv_vec3 v = fv_vec3_make(1, -2, 3);
    fv_vec3 r;

    CHECK(id.m[0] == 1 && id.m[1] == 0 && id.m[2] == 0);
    CHECK(id.m[3] == 0 && id.m[4] == 1 && id.m[5] == 0);
    CHECK(id.m[6] == 0 && id.m[7] == 0 && id.m[8] == 1);

    r = fv_mat3_transform_vec3(&id, v);
    CHECK(vec3_bits_equal(r, v));

    /* from_rows is exactly the row-major packing. */
    {
        fv_mat3 m = fv_mat3_from_rows(fv_vec3_make(1, 2, 3),
                                      fv_vec3_make(4, 5, 6),
                                      fv_vec3_make(7, 8, 9));
        CHECK(mat3_bits_equal(m, make_mat(1, 2, 3, 4, 5, 6, 7, 8, 9)));
        CHECK(fv_mat3_transform_vec3(&m, fv_vec3_make(1, 0, 0)).x == 1);
        CHECK(fv_mat3_transform_vec3(&m, fv_vec3_make(1, 0, 0)).y == 4);
        CHECK(fv_mat3_transform_vec3(&m, fv_vec3_make(1, 0, 0)).z == 7);
        /* The independent reference agrees bit for bit on the same matrix. */
        CHECK(vec3_bits_equal(fv_mat3_transform_vec3(&m, fv_vec3_make(2, -3, 4)),
                              ref_tf(&m, fv_vec3_make(2, -3, 4))));
    }

    /* Scale and shear, with exact expected values. */
    {
        fv_mat3 s = make_mat(2, 0, 0, 0, 3, 0, 0, 0, 4);
        r = fv_mat3_transform_vec3(&s, fv_vec3_make(1, 2, 3));
        CHECK(r.x == 2 && r.y == 6 && r.z == 12);

        fv_mat3 sh = make_mat(1, 1, 0, 0, 1, 0, 0, 0, 1);
        r = fv_mat3_transform_vec3(&sh, fv_vec3_make(2, 3, 4));
        CHECK(r.x == 5 && r.y == 3 && r.z == 4);
    }

    /* transpose twice is the original, bit for bit (component moves only). */
    {
        fv_mat3 m = make_mat(1, 2, 3, 4, 5, 6, 7, 8, 9);
        fv_mat3 t = fv_mat3_transpose(&m);
        CHECK(t.m[0] == 1 && t.m[1] == 4 && t.m[2] == 7);
        CHECK(t.m[3] == 2 && t.m[4] == 5 && t.m[5] == 8);
        CHECK(t.m[6] == 3 && t.m[7] == 6 && t.m[8] == 9);
        CHECK(mat3_bits_equal(fv_mat3_transpose(&t), m));
    }

    /* Identity is both a left and a right identity for mul. */
    {
        fv_mat3 m = make_mat(1, 2, 3, 4, 5, 6, 7, 8, 9);
        CHECK(mat3_bits_equal(fv_mat3_mul(&id, &m), m));
        CHECK(mat3_bits_equal(fv_mat3_mul(&m, &id), m));
    }
}

static void test_mul(void)
{
    fv_mat3 a = make_mat(1, 2, 3, 4, 5, 6, 7, 8, 9);
    fv_mat3 b = make_mat(-1, 0.5f, 2, 3, -2, 1, 0.25f, 4, -3);
    int k;

    /* Independent reference, bit for bit. */
    CHECK(mat3_bits_equal(fv_mat3_mul(&a, &b), ref_mul(&a, &b)));

    /* Composition order: (A*B)*v == A*(B*v). Different expressions -> tolerance. */
    {
        fv_vec3 v = fv_vec3_make(1, -2, 3);
        fv_vec3 left = fv_mat3_transform_vec3(&a, fv_mat3_transform_vec3(&b, v));
        fv_mat3 ab = fv_mat3_mul(&a, &b);
        fv_vec3 right = fv_mat3_transform_vec3(&ab, v);
        CHECK(vec3_near(left, right, (fv_real)2e-5));
    }

    for (k = 0; k < 2000; ++k) {
        fv_mat3 x = make_mat(rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_mat3 y = make_mat(rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate());
        CHECK(mat3_bits_equal(fv_mat3_mul(&x, &y), ref_mul(&x, &y)));
    }
}

static void test_from_quat(void)
{
    fv_quat id = fv_quat_identity();
    fv_vec3 v = fv_vec3_make(1, 0, 0);

    /* Identity quaternion -> identity matrix, bit for bit. */
    CHECK(mat3_bits_equal(fv_mat3_from_quat(id), fv_mat3_identity()));

    /* 90 degrees about +y: the same rotation as the quaternion path. */
    {
        fv_quat qy = fv_quat_make(0, 0.70710677f, 0, 0.70710677f);
        fv_mat3 m = fv_mat3_from_quat(qy);
        fv_vec3 a = fv_mat3_transform_vec3(&m, v);
        fv_vec3 b = fv_quat_rotate_vec3(qy, v);
        CHECK(vec3_near(a, b, (fv_real)2e-6));
        CHECK(vec3_near(a, fv_vec3_make(0, 0, -1), (fv_real)2e-6));
    }

    /* A non-axis rotation: matrix and quaternion forms agree within rounding. */
    {
        fv_quat q = fv_quat_normalize(fv_quat_make(0.2f, 0.3f, -0.1f, 0.9f));
        fv_mat3 m = fv_mat3_from_quat(q);
        fv_vec3 a = fv_mat3_transform_vec3(&m, v);
        fv_vec3 b = fv_quat_rotate_vec3(q, v);
        CHECK(vec3_near(a, b, (fv_real)5e-6));
    }
}

static void test_inverse(void)
{
    fv_mat3 id = fv_mat3_identity();
    fv_mat3 zero;
    fv_mat3 inv;
    fv_mat3 p;
    int k;

    memset(&zero, 0, sizeof(zero));

    CHECK(fv_mat3_determinant(&id) == 1);
    CHECK(mat3_bits_equal(fv_mat3_inverse(&id), id));

    /* A known invertible matrix: det = 9, and m * inverse(m) is the identity. */
    {
        fv_mat3 m = make_mat(4, 7, 2, 3, 6, 1, 2, 5, 3);
        fv_real det = fv_mat3_determinant(&m);
        CHECK(feq(det, (fv_real)9, (fv_real)1e-5));
        inv = fv_mat3_inverse(&m);
        p = fv_mat3_mul(&m, &inv);
        CHECK(vec3_near(fv_mat3_transform_vec3(&p, fv_vec3_make(1, 0, 0)), fv_vec3_make(1, 0, 0), (fv_real)1e-5));
        CHECK(vec3_near(fv_mat3_transform_vec3(&p, fv_vec3_make(0, 1, 0)), fv_vec3_make(0, 1, 0), (fv_real)1e-5));
        CHECK(vec3_near(fv_mat3_transform_vec3(&p, fv_vec3_make(0, 0, 1)), fv_vec3_make(0, 0, 1), (fv_real)1e-5));
    }

    /* Negative determinant is still invertible. */
    {
        fv_mat3 m = make_mat(1, 0, 0, 0, 1, 0, 0, 0, -1);
        CHECK(fv_mat3_determinant(&m) == -1);
        CHECK(fv_mat3_inverse_checked(&m, &inv) == 1);
        CHECK(mat3_component_equal(inv, m));
    }

    /* Singular inputs: zero matrix, rank-deficient matrix, and a scale small enough to fall under
     * the relative threshold. The value form is the zero matrix, the checked form reports 0. */
    {
        fv_mat3 rank1 = make_mat(1, 2, 3, 2, 4, 6, 0, 1, 0);
        fv_mat3 small = make_mat(1, 0, 0, 0, 1, 0, 0, 0, (fv_real)1e-7);
        fv_mat3 out = make_mat(9, 9, 9, 9, 9, 9, 9, 9, 9);
        CHECK(fv_mat3_determinant(&rank1) == 0);
        CHECK(mat3_bits_equal(fv_mat3_inverse(&zero), zero));
        CHECK(mat3_bits_equal(fv_mat3_inverse(&rank1), zero));
        CHECK(mat3_bits_equal(fv_mat3_inverse(&small), zero));
        CHECK(fv_mat3_inverse_checked(&zero, &out) == 0);
        CHECK(mat3_bits_equal(out, zero));
        CHECK(fv_mat3_inverse_checked(&rank1, &out) == 0);
        CHECK(fv_mat3_inverse_checked(&small, &out) == 0);
    }

    /* Just above the threshold the same shape is accepted (det 1e-5 > 1e-6). */
    {
        fv_mat3 near_singular = make_mat(1, 0, 0, 0, 1, 0, 0, 0, (fv_real)1e-5);
        CHECK(fv_mat3_inverse_checked(&near_singular, &inv) == 1);
    }

    /* out may alias m: the implementation reads everything before writing. */
    {
        fv_mat3 m = make_mat(4, 7, 2, 3, 6, 1, 2, 5, 3);
        fv_mat3 orig = m;
        CHECK(fv_mat3_inverse_checked(&m, &m) == 1);
        p = fv_mat3_mul(&orig, &m);
        CHECK(vec3_near(fv_mat3_transform_vec3(&p, fv_vec3_make(1, 0, 0)), fv_vec3_make(1, 0, 0), (fv_real)1e-5));
    }

    /* Independent reference, bit for bit, and m * inverse(m) = I for well-conditioned inputs. */
    for (k = 0; k < 2000; ++k) {
        fv_mat3 m = make_mat(rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate(),
                             rnd_moderate(), rnd_moderate(), rnd_moderate());
        fv_real det = fv_mat3_determinant(&m);
        CHECK(same_bits_float(det, ref_det(&m)));
        if (fv_mat3_inverse_checked(&m, &inv) == 1) {
            CHECK(mat3_bits_equal(inv, ref_inverse(&m)));
            if (fabsf((float)det) > 0.05f) {
                p = fv_mat3_mul(&m, &inv);
                CHECK(vec3_near(fv_mat3_transform_vec3(&p, fv_vec3_make(1, 0, 0)),
                                fv_vec3_make(1, 0, 0), (fv_real)1e-3));
            }
        }
    }
}

static void poison_planes(float p[9][CAP])
{
    int k;
    size_t i;
    for (k = 0; k < 9; ++k) {
        for (i = 0; i < CAP; ++i) {
            p[k][i] = 1.0e30f;
        }
    }
}

static int planes_guard_ok(float p[9][CAP], size_t from)
{
    int k;
    size_t i;
    for (k = 0; k < 9; ++k) {
        for (i = from; i < CAP; ++i) {
            if (p[k][i] != 1.0e30f) {
                return 0;
            }
        }
    }
    return 1;
}

static void fill_planes(size_t n)
{
    size_t i;
    int k;
    for (k = 0; k < 9; ++k) {
        for (i = 0; i < n; ++i) {
            g_a[k][i] = rnd_moderate();
            g_b[k][i] = rnd_moderate();
        }
    }
}

static void test_batch_mul(size_t n)
{
    size_t i;
    int k;
    fill_planes(n);
    /* mul does not have a singular path; all entries are finite random values. */
    poison_planes(g_om);
    fv_mat3_mul_soa(g_a[0], g_a[1], g_a[2], g_a[3], g_a[4], g_a[5], g_a[6], g_a[7], g_a[8],
                    g_b[0], g_b[1], g_b[2], g_b[3], g_b[4], g_b[5], g_b[6], g_b[7], g_b[8],
                    g_om[0], g_om[1], g_om[2], g_om[3], g_om[4], g_om[5], g_om[6], g_om[7],
                    g_om[8], n);
    for (i = 0; i < n; ++i) {
        fv_mat3 a;
        fv_mat3 b;
        fv_mat3 want;
        for (k = 0; k < 9; ++k) {
            a.m[k] = g_a[k][i];
            b.m[k] = g_b[k][i];
        }
        want = fv_mat3_mul(&a, &b);
        for (k = 0; k < 9; ++k) {
            CHECK(g_om[k][i] == want.m[k]);
        }
    }
    CHECK(planes_guard_ok(g_om, n));
}

static void test_batch_inverse(size_t n)
{
    size_t i;
    int k;
    fill_planes(n);
    /* Every 11th lane is the zero matrix: the batch must take the same singular path as the scalar
     * function and write nine zeros. */
    for (i = 0; i < n; i += 11) {
        for (k = 0; k < 9; ++k) {
            g_a[k][i] = 0.0f;
        }
    }
    poison_planes(g_oi);
    fv_mat3_inverse_soa(g_a[0], g_a[1], g_a[2], g_a[3], g_a[4], g_a[5], g_a[6], g_a[7], g_a[8],
                        g_oi[0], g_oi[1], g_oi[2], g_oi[3], g_oi[4], g_oi[5], g_oi[6], g_oi[7],
                        g_oi[8], n);
    for (i = 0; i < n; ++i) {
        fv_mat3 a;
        fv_mat3 want;
        for (k = 0; k < 9; ++k) {
            a.m[k] = g_a[k][i];
        }
        want = fv_mat3_inverse(&a);
        for (k = 0; k < 9; ++k) {
            CHECK(g_oi[k][i] == want.m[k]);
        }
    }
    CHECK(planes_guard_ok(g_oi, n));
}

static void fill_case(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_m00[i] = rnd_moderate(); g_m01[i] = rnd_moderate(); g_m02[i] = rnd_moderate();
        g_m10[i] = rnd_moderate(); g_m11[i] = rnd_moderate(); g_m12[i] = rnd_moderate();
        g_m20[i] = rnd_moderate(); g_m21[i] = rnd_moderate(); g_m22[i] = rnd_moderate();
        g_vx[i] = rnd_moderate(); g_vy[i] = rnd_moderate(); g_vz[i] = rnd_moderate();
    }
}

static void test_batch(size_t n)
{
    size_t i;
    fill_case(n);
    poison_float(g_ox, CAP);
    poison_float(g_oy, CAP);
    poison_float(g_oz, CAP);
    fv_mat3_transform_vec3_soa(g_m00, g_m01, g_m02, g_m10, g_m11, g_m12,
                               g_m20, g_m21, g_m22, g_vx, g_vy, g_vz,
                               g_ox, g_oy, g_oz, n);
    for (i = 0; i < n; ++i) {
        fv_mat3 m = make_mat(g_m00[i], g_m01[i], g_m02[i],
                             g_m10[i], g_m11[i], g_m12[i],
                             g_m20[i], g_m21[i], g_m22[i]);
        fv_vec3 v = fv_vec3_make(g_vx[i], g_vy[i], g_vz[i]);
        fv_vec3 want = fv_mat3_transform_vec3(&m, v);
        CHECK(g_ox[i] == want.x && g_oy[i] == want.y && g_oz[i] == want.z);
    }
    CHECK(guard_ok(g_ox, n, CAP));
    CHECK(guard_ok(g_oy, n, CAP));
    CHECK(guard_ok(g_oz, n, CAP));
}

int main(void)
{
    size_t n;

    printf("== fv_mat3 scalar semantics ==\n");
    test_basics();
    test_mul();
    test_from_quat();
    test_inverse();

    printf("== fv_mat3 SoA batch: bit-identical to the per-value API (n = 0..40 all tails + 63 + 1024 + 4097) ==\n");
    test_batch(0);
    test_batch(63);
    test_batch(1024);
    test_batch(4097);
    test_batch_mul(0);
    test_batch_mul(63);
    test_batch_mul(1024);
    test_batch_mul(4097);
    test_batch_inverse(0);
    test_batch_inverse(63);
    test_batch_inverse(1024);
    test_batch_inverse(4097);
    for (n = 0; n <= 40; ++n) {
        test_batch(n);
        test_batch_mul(n);
        test_batch_inverse(n);
    }

    fv_mat3_transform_vec3_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_mat3_mul_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_mat3_inverse_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);

    printf("test_mat3: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
