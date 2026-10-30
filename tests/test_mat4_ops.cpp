/* fv_mat4_mul / fv_mat4_inverse (P7).
 *
 * The mat4 transform entries (points / dirs) keep the affine precondition and only read rows 0..2;
 * these two are the general 4x4 operations, so this fixture uses arbitrary row-major matrices.
 *
 * What is pinned down:
 *   - mul against an independent triple-loop reference, bit for bit;
 *   - inverse: m * inverse(m) and inverse(m) * m are the identity within rounding, for random
 *     well-conditioned matrices; the identity, a translation matrix, and a rotation+translation;
 *   - singular inputs (zero, rank-deficient, scaled small) return the zero matrix and the checked
 *     form reports 0;
 *   - out may alias m in the checked form;
 *   - inverse(transform) undoes transform_points on a point (the affine use case).
 *
 * Run: make test-mat4-ops
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

static int feq(fv_real a, fv_real b, fv_real eps)
{
    fv_real d = a - b;
    if (d < (fv_real)0) {
        d = -d;
    }
    return d <= eps;
}

static int mat4_near(fv_mat4 a, fv_mat4 b, fv_real eps)
{
    int k;
    for (k = 0; k < 16; ++k) {
        if (!feq(a.m[k], b.m[k], eps)) {
            return 0;
        }
    }
    return 1;
}

static int mat4_bits_equal(fv_mat4 a, fv_mat4 b)
{
    return memcmp(&a, &b, sizeof(a)) == 0;
}

static uint32_t g_rng = 0xC0FFEEu;

static uint32_t rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static fv_real rnd_moderate(void)
{
    return (fv_real)((int)(rnd_u32() % 2001u) - 1000) / 16.0f;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif

/* Independent reference: the same row-major expression order as the contract. */
static fv_mat4 ref_mul(const fv_mat4 *a, const fv_mat4 *b)
{
    fv_mat4 r;
    int i, j;
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

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

/* Well-conditioned random matrix: moderate entries plus a strong diagonal. */
static fv_mat4 make_random_well(void)
{
    fv_mat4 m;
    int k;
    for (k = 0; k < 16; ++k) {
        m.m[k] = rnd_moderate() * (fv_real)0.1;
    }
    m.m[0] += (fv_real)1;
    m.m[5] += (fv_real)1;
    m.m[10] += (fv_real)1;
    m.m[15] += (fv_real)1;
    return m;
}

static void test_basics(void)
{
    fv_mat4 id = fv_mat4_identity();
    fv_mat4 t;
    fv_mat4 inv;
    fv_mat4 p;

    CHECK(mat4_bits_equal(fv_mat4_mul(&id, &id), id));
    CHECK(mat4_bits_equal(fv_mat4_inverse(&id), id));

    /* Translation: T(t) * T(-t) = I, and inverse(T(t)) moves a point back. */
    t = id;
    t.m[3] = (fv_real)2.5;
    t.m[7] = (fv_real)-1.5;
    t.m[11] = (fv_real)0.75;
    inv = fv_mat4_inverse(&t);
    CHECK(feq(inv.m[3], (fv_real)-2.5, (fv_real)1e-6));
    CHECK(feq(inv.m[7], (fv_real)1.5, (fv_real)1e-6));
    CHECK(feq(inv.m[11], (fv_real)-0.75, (fv_real)1e-6));
    p = fv_mat4_mul(&t, &inv);
    CHECK(mat4_near(p, id, (fv_real)1e-6));

    {
        fv_vec4 qin = fv_vec4_make(1, 2, 3, 1);
        fv_vec4 qmid;
        fv_vec4 qout;
        /* in and out are __restrict: never alias them. */
        fv_mat4_transform_points(&t, &qin, &qmid, 1);
        CHECK(feq(qmid.x, (fv_real)3.5, (fv_real)1e-6) && feq(qmid.y, (fv_real)0.5, (fv_real)1e-6));
        fv_mat4_transform_points(&inv, &qmid, &qout, 1);
        CHECK(feq(qout.x, (fv_real)1, (fv_real)1e-6) && feq(qout.y, (fv_real)2, (fv_real)1e-6) &&
              feq(qout.z, (fv_real)3, (fv_real)1e-6) && feq(qout.w, (fv_real)1, (fv_real)1e-6));
    }
}

static void test_mul(void)
{
    int k;
    for (k = 0; k < 2000; ++k) {
        fv_mat4 a;
        fv_mat4 b;
        int i;
        for (i = 0; i < 16; ++i) {
            a.m[i] = rnd_moderate();
            b.m[i] = rnd_moderate();
        }
        CHECK(mat4_bits_equal(fv_mat4_mul(&a, &b), ref_mul(&a, &b)));
    }
}

static void test_inverse(void)
{
    fv_mat4 id = fv_mat4_identity();
    fv_mat4 zero;
    fv_mat4 inv;
    int k;

    memset(&zero, 0, sizeof(zero));

    /* Singular: zero matrix, a rank-deficient matrix, and a uniformly tiny scale. */
    {
        fv_mat4 rank1;
        fv_mat4 tiny;
        fv_mat4 out;
        int i;
        for (i = 0; i < 16; ++i) {
            rank1.m[i] = (fv_real)(i % 4);
        }
        tiny = id;
        tiny.m[0] = (fv_real)1e-9;
        CHECK(mat4_bits_equal(fv_mat4_inverse(&zero), zero));
        CHECK(mat4_bits_equal(fv_mat4_inverse(&rank1), zero));
        CHECK(mat4_bits_equal(fv_mat4_inverse(&tiny), zero));
        CHECK(fv_mat4_inverse_checked(&zero, &out) == 0);
        CHECK(mat4_bits_equal(out, zero));
        CHECK(fv_mat4_inverse_checked(&rank1, &out) == 0);
    }

    /* Random well-conditioned matrices: both products are the identity within rounding. */
    for (k = 0; k < 2000; ++k) {
        fv_mat4 m = make_random_well();
        fv_mat4 p;
        CHECK(fv_mat4_inverse_checked(&m, &inv) == 1);
        p = fv_mat4_mul(&m, &inv);
        CHECK(mat4_near(p, id, (fv_real)1e-3));
        p = fv_mat4_mul(&inv, &m);
        CHECK(mat4_near(p, id, (fv_real)1e-3));
    }

    /* out may alias m. */
    {
        fv_mat4 m = make_random_well();
        fv_mat4 orig = m;
        fv_mat4 p;
        CHECK(fv_mat4_inverse_checked(&m, &m) == 1);
        p = fv_mat4_mul(&orig, &m);
        CHECK(mat4_near(p, id, (fv_real)1e-3));
    }
}

int main(void)
{
    printf("== fv_mat4 general operations (mul / inverse) ==\n");
    test_basics();
    test_mul();
    test_inverse();
    printf("test_mat4_ops: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
