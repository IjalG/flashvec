/* fv_quat / fv_mat3 from a plain C11 translation unit (P6).
 *
 * What this fixture is for: the header is a pure C interface, and the new types/functions must compile
 * and behave identically under a C compiler -- not only under C++. It checks:
 *   - the scalar API is callable and links from C (extern symbols in libflashvec.a);
 *   - the header's static inline SoA batch entries are usable from C;
 *   - the batch results are bit-identical to the per-value results (memcmp, not tolerance);
 *   - n == 0 with NULL pointers is accepted;
 *   - the identity / rotation / transform semantics are the documented ones.
 *
 * Run: make test-physics-c
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

#define CAP 80u

static float ax[CAP], ay[CAP], az[CAP], aw[CAP];
static float bx[CAP], by[CAP], bz[CAP], bw[CAP];
static float vx[CAP], vy[CAP], vz[CAP];
static float ox[CAP], oy[CAP], oz[CAP], ow[CAP];
static float m00[CAP], m01[CAP], m02[CAP];
static float m10[CAP], m11[CAP], m12[CAP];
static float m20[CAP], m21[CAP], m22[CAP];

static int same_bits(float a, float b)
{
    uint32_t ua;
    uint32_t ub;
    memcpy(&ua, &a, sizeof(ua));
    memcpy(&ub, &b, sizeof(ub));
    return ua == ub;
}

static float seq(int i, int mul)
{
    return (float)((i * mul) % 23) / 8.0f - 1.0f;
}

static void test_quat_scalar(void)
{
    fv_quat id = fv_quat_identity();
    fv_quat a = fv_quat_make(1, 2, 3, 4);
    fv_quat b = fv_quat_make(-1, 0.5f, -3, 2);
    fv_quat p = fv_quat_mul(a, b);
    fv_quat c = fv_quat_conjugate(a);

    CHECK(id.x == 0 && id.y == 0 && id.z == 0 && id.w == 1);
    CHECK(c.x == -1 && c.y == -2 && c.z == -3 && c.w == 4);
    CHECK(same_bits(p.x, ((a.w * b.x + a.x * b.w) + a.y * b.z) - a.z * b.y));
    CHECK(same_bits(p.y, ((a.w * b.y - a.x * b.z) + a.y * b.w) + a.z * b.x));
    CHECK(same_bits(p.z, ((a.w * b.z + a.x * b.y) - a.y * b.x) + a.z * b.w));
    CHECK(same_bits(p.w, ((a.w * b.w - a.x * b.x) - a.y * b.y) - a.z * b.z));

    /* Normalize of zero is the identity, not zero. */
    {
        fv_quat z = fv_quat_make(0, 0, 0, 0);
        fv_quat n = fv_quat_normalize(z);
        CHECK(n.x == 0 && n.y == 0 && n.z == 0 && n.w == 1);
    }

    /* 90 degrees about +y maps (1,0,0) to (0,0,-1). */
    {
        fv_quat qy = fv_quat_make(0, sqrtf(0.5f), 0, sqrtf(0.5f));
        fv_vec3 r = fv_quat_rotate_vec3(qy, fv_vec3_make(1, 0, 0));
        CHECK(fabsf(r.x) < 1e-6f);
        CHECK(fabsf(r.y) < 1e-6f);
        CHECK(fabsf(r.z + 1.0f) < 1e-6f);
    }

    /* q * q^-1 is the identity up to rounding. */
    {
        fv_quat q = fv_quat_normalize(fv_quat_make(0.2f, -0.3f, 0.4f, 0.8f));
        fv_quat u = fv_quat_mul(q, fv_quat_conjugate(q));
        CHECK(fabsf(u.x) < 1e-6f && fabsf(u.y) < 1e-6f && fabsf(u.z) < 1e-6f);
        CHECK(fabsf(u.w - 1.0f) < 1e-6f);
    }

    /* General inverse: q * inverse(q) is the identity; the zero quaternion is reported as invalid. */
    {
        fv_quat q = fv_quat_make(0.2f, -0.3f, 0.4f, 0.8f);
        fv_quat inv = fv_quat_inverse(q);
        fv_quat u = fv_quat_mul(q, inv);
        fv_quat z = fv_quat_make(0, 0, 0, 0);
        fv_quat out = fv_quat_make(9, 9, 9, 9);
        CHECK(fabsf(u.x) < 1e-6f && fabsf(u.y) < 1e-6f && fabsf(u.z) < 1e-6f);
        CHECK(fabsf(u.w - 1.0f) < 1e-6f);
        CHECK(fv_quat_inverse_checked(z, &out) == 0);
        CHECK(out.x == 0 && out.y == 0 && out.z == 0 && out.w == 0);
    }
}

static void test_mat3_scalar(void)
{
    fv_mat3 id = fv_mat3_identity();
    fv_mat3 m = fv_mat3_from_rows(fv_vec3_make(1, 2, 3),
                                  fv_vec3_make(4, 5, 6),
                                  fv_vec3_make(7, 8, 9));
    fv_vec3 r = fv_mat3_transform_vec3(&m, fv_vec3_make(1, 0, 0));

    CHECK(id.m[0] == 1 && id.m[4] == 1 && id.m[8] == 1);
    CHECK(r.x == 1 && r.y == 4 && r.z == 7);

    /* Transpose twice is the original. */
    {
        fv_mat3 t = fv_mat3_transpose(&m);
        fv_mat3 tt = fv_mat3_transpose(&t);
        CHECK(memcmp(&m, &tt, sizeof(m)) == 0);
    }

    /* Identity quaternion -> identity matrix. */
    {
        fv_mat3 qm = fv_mat3_from_quat(fv_quat_identity());
        CHECK(memcmp(&qm, &id, sizeof(id)) == 0);
    }

    /* Inverse: identity, a known matrix (det = 9), and the singular/checked path. */
    {
        fv_mat3 inv;
        fv_mat3 known = fv_mat3_from_rows(fv_vec3_make(4, 7, 2),
                                          fv_vec3_make(3, 6, 1),
                                          fv_vec3_make(2, 5, 3));
        fv_mat3 zero;
        fv_mat3 p;
        memset(&zero, 0, sizeof(zero));
        CHECK(fv_mat3_determinant(&id) == 1);
        CHECK(fv_mat3_inverse_checked(&id, &inv) == 1);
        CHECK(memcmp(&inv, &id, sizeof(id)) == 0);
        CHECK(fabsf(fv_mat3_determinant(&known) - 9.0f) < 1e-5f);
        CHECK(fv_mat3_inverse_checked(&known, &inv) == 1);
        p = fv_mat3_mul(&known, &inv);
        CHECK(fabsf(fv_mat3_transform_vec3(&p, fv_vec3_make(1, 0, 0)).x - 1.0f) < 1e-5f);
        CHECK(fv_mat3_inverse_checked(&zero, &inv) == 0);
        CHECK(memcmp(&inv, &zero, sizeof(zero)) == 0);
    }
}

static void fill(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        ax[i] = seq((int)i, 3);
        ay[i] = seq((int)i, 5);
        az[i] = seq((int)i, 7);
        aw[i] = seq((int)i, 11);
        bx[i] = seq((int)i, 13);
        by[i] = seq((int)i, 17);
        bz[i] = seq((int)i, 19);
        bw[i] = seq((int)i, 23);
        vx[i] = seq((int)i, 2);
        vy[i] = seq((int)i, 3);
        vz[i] = seq((int)i, 5);
        m00[i] = seq((int)i, 3); m01[i] = seq((int)i, 5); m02[i] = seq((int)i, 7);
        m10[i] = seq((int)i, 11); m11[i] = seq((int)i, 13); m12[i] = seq((int)i, 17);
        m20[i] = seq((int)i, 19); m21[i] = seq((int)i, 23); m22[i] = seq((int)i, 29);
    }
}

static void test_batch(size_t n)
{
    size_t i;
    fill(n);

    fv_quat_mul_soa(ax, ay, az, aw, bx, by, bz, bw, ox, oy, oz, ow, n);
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(ax[i], ay[i], az[i], aw[i]);
        fv_quat b = fv_quat_make(bx[i], by[i], bz[i], bw[i]);
        fv_quat want = fv_quat_mul(a, b);
        CHECK(same_bits(ox[i], want.x) && same_bits(oy[i], want.y) &&
              same_bits(oz[i], want.z) && same_bits(ow[i], want.w));
    }

    fv_quat_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(ax[i], ay[i], az[i], aw[i]);
        fv_quat want = fv_quat_normalize(a);
        CHECK(same_bits(ox[i], want.x) && same_bits(oy[i], want.y) &&
              same_bits(oz[i], want.z) && same_bits(ow[i], want.w));
    }

    /* Start from unit quaternions: the rotate entry has the same precondition as the per-value one. */
    for (i = 0; i < n; ++i) {
        fv_quat q = fv_quat_normalize(fv_quat_make(ax[i], ay[i], az[i], aw[i]));
        ax[i] = q.x;
        ay[i] = q.y;
        az[i] = q.z;
        aw[i] = q.w;
    }
    fv_quat_rotate_vec3_soa(ax, ay, az, aw, vx, vy, vz, ox, oy, oz, n);
    for (i = 0; i < n; ++i) {
        fv_quat q = fv_quat_make(ax[i], ay[i], az[i], aw[i]);
        fv_vec3 v = fv_vec3_make(vx[i], vy[i], vz[i]);
        fv_vec3 want = fv_quat_rotate_vec3(q, v);
        CHECK(same_bits(ox[i], want.x) && same_bits(oy[i], want.y) && same_bits(oz[i], want.z));
    }

    fv_mat3_transform_vec3_soa(m00, m01, m02, m10, m11, m12, m20, m21, m22,
                               vx, vy, vz, ox, oy, oz, n);
    for (i = 0; i < n; ++i) {
        fv_mat3 m;
        fv_vec3 v = fv_vec3_make(vx[i], vy[i], vz[i]);
        fv_vec3 want;
        m.m[0] = m00[i]; m.m[1] = m01[i]; m.m[2] = m02[i];
        m.m[3] = m10[i]; m.m[4] = m11[i]; m.m[5] = m12[i];
        m.m[6] = m20[i]; m.m[7] = m21[i]; m.m[8] = m22[i];
        want = fv_mat3_transform_vec3(&m, v);
        CHECK(same_bits(ox[i], want.x) && same_bits(oy[i], want.y) && same_bits(oz[i], want.z));
    }
}

int main(void)
{
    size_t n;

    printf("== fv_quat / fv_mat3 from C11 ==\n");
    test_quat_scalar();
    test_mat3_scalar();

    test_batch(0);
    test_batch(1);
    test_batch(7);
    test_batch(8);
    test_batch(9);
    test_batch(63);
    test_batch(64);
    for (n = 0; n <= 40; ++n) {
        test_batch(n);
    }

    /* n == 0 touches no pointer. */
    fv_quat_mul_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL, NULL, 0);
    fv_quat_normalize_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_quat_rotate_vec3_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_mat3_transform_vec3_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, 0);

    printf("test_physics_c: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
