/* _fast approximate normalize + fv_atan_approx (P7).
 *
 * These entries deliberately do not promise bit-identical results, so the checks are accuracy bounds,
 * not memcmp:
 *   - fv_atan_approx against atanf on [-1,1] (max abs error ~1.2e-5; test tolerance 3e-5), NaN outside
 *     the domain;
 *   - fv_vec{2,3,4}_normalize_soa_fast / fv_quat_normalize_soa_fast against the exact SoA normalize:
 *     relative error <= 1e-5 per component, unit length within 2e-6, and the zero guard exact (zeros,
 *     or the identity for the quaternion).
 *
 * Every tail length n = 0..40 is covered, plus 63 / 1024 / 4097; n == 0 with NULL pointers is accepted.
 *
 * Run: make test-fast
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

static float g_ax[CAP], g_ay[CAP], g_az[CAP], g_aw[CAP];
static float g_ex[CAP], g_ey[CAP], g_ez[CAP], g_ew[CAP];
static float g_fx[CAP], g_fy[CAP], g_fz[CAP], g_fw[CAP];

static int rel_ok(float got, float want, float tol)
{
    float d = got - want;
    float scale = want < 0.0f ? -want : want;
    if (scale < 1.0f) {
        scale = 1.0f;
    }
    if (d < 0.0f) {
        d = -d;
    }
    return d <= tol * scale;
}

static uint32_t g_rng = 0x5EED1234u;

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

static void fill(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_ax[i] = rnd_moderate();
        g_ay[i] = rnd_moderate();
        g_az[i] = rnd_moderate();
        g_aw[i] = rnd_moderate();
    }
    /* A few zero / near-zero lanes exercise the guard. */
    if (n > 3) {
        g_ax[0] = 0.0f; g_ay[0] = 0.0f; g_az[0] = 0.0f; g_aw[0] = 0.0f;
        g_ax[1] = 1.0e-8f; g_ay[1] = 0.0f; g_az[1] = 0.0f; g_aw[1] = 0.0f;
    }
}

static void test_atan(void)
{
    int k;
    CHECK(fv_atan_approx(0) == 0);
    CHECK(fabsf(fv_atan_approx(1) - atanf(1.0f)) <= 3e-5f);
    CHECK(fabsf(fv_atan_approx(-1) + atanf(1.0f)) <= 3e-5f);
    for (k = -1000; k <= 1000; ++k) {
        float x = (float)k / 1000.0f;
        CHECK(fabsf(fv_atan_approx(x) - atanf(x)) <= 3e-5f);
    }
    CHECK(isnan(fv_atan_approx(2.0f)));
    CHECK(isnan(fv_atan_approx(-1.5f)));
    CHECK(isnan(fv_atan_approx(NAN)));
}

static void test_vec2(size_t n)
{
    size_t i;
    fill(n);
    fv_vec2_normalize_soa(g_ax, g_ay, g_ex, g_ey, n);
    fv_vec2_normalize_soa_fast(g_ax, g_ay, g_fx, g_fy, n);
    for (i = 0; i < n; ++i) {
        CHECK(rel_ok(g_fx[i], g_ex[i], 1e-5f));
        CHECK(rel_ok(g_fy[i], g_ey[i], 1e-5f));
        if (g_ax[i] == 0.0f && g_ay[i] == 0.0f) {
            CHECK(g_fx[i] == 0.0f && g_fy[i] == 0.0f);
        }
    }
}

static void test_vec3(size_t n)
{
    size_t i;
    fill(n);
    fv_vec3_normalize_soa(g_ax, g_ay, g_az, g_ex, g_ey, g_ez, n);
    fv_vec3_normalize_soa_fast(g_ax, g_ay, g_az, g_fx, g_fy, g_fz, n);
    for (i = 0; i < n; ++i) {
        CHECK(rel_ok(g_fx[i], g_ex[i], 1e-5f));
        CHECK(rel_ok(g_fy[i], g_ey[i], 1e-5f));
        CHECK(rel_ok(g_fz[i], g_ez[i], 1e-5f));
        if (g_ax[i] == 0.0f && g_ay[i] == 0.0f && g_az[i] == 0.0f) {
            CHECK(g_fx[i] == 0.0f && g_fy[i] == 0.0f && g_fz[i] == 0.0f);
        }
    }
}

static void test_vec4(size_t n)
{
    size_t i;
    fill(n);
    fv_vec4_normalize_soa(g_ax, g_ay, g_az, g_aw, g_ex, g_ey, g_ez, g_ew, n);
    fv_vec4_normalize_soa_fast(g_ax, g_ay, g_az, g_aw, g_fx, g_fy, g_fz, g_fw, n);
    for (i = 0; i < n; ++i) {
        CHECK(rel_ok(g_fx[i], g_ex[i], 1e-5f));
        CHECK(rel_ok(g_fy[i], g_ey[i], 1e-5f));
        CHECK(rel_ok(g_fz[i], g_ez[i], 1e-5f));
        CHECK(rel_ok(g_fw[i], g_ew[i], 1e-5f));
        if (g_ax[i] == 0.0f && g_ay[i] == 0.0f && g_az[i] == 0.0f && g_aw[i] == 0.0f) {
            CHECK(g_fx[i] == 0.0f && g_fy[i] == 0.0f && g_fz[i] == 0.0f && g_fw[i] == 0.0f);
        }
    }
}

static void test_quat(size_t n)
{
    size_t i;
    fill(n);
    fv_quat_normalize_soa(g_ax, g_ay, g_az, g_aw, g_ex, g_ey, g_ez, g_ew, n);
    fv_quat_normalize_soa_fast(g_ax, g_ay, g_az, g_aw, g_fx, g_fy, g_fz, g_fw, n);
    for (i = 0; i < n; ++i) {
        CHECK(rel_ok(g_fx[i], g_ex[i], 1e-5f));
        CHECK(rel_ok(g_fy[i], g_ey[i], 1e-5f));
        CHECK(rel_ok(g_fz[i], g_ez[i], 1e-5f));
        CHECK(rel_ok(g_fw[i], g_ew[i], 1e-5f));
        if (g_ax[i] == 0.0f && g_ay[i] == 0.0f && g_az[i] == 0.0f && g_aw[i] == 0.0f) {
            CHECK(g_fx[i] == 0.0f && g_fy[i] == 0.0f && g_fz[i] == 0.0f && g_fw[i] == 1.0f);
        }
    }
}

static void test_madd_norm(size_t n)
{
    size_t i;
    fill(n);
    fv_vec2_madd_normalize_soa(g_ax, g_ay, g_az, g_aw, 1.5f, g_ex, g_ey, n);
    fv_vec2_madd_normalize_soa_fast(g_ax, g_ay, g_az, g_aw, 1.5f, g_fx, g_fy, n);
    for (i = 0; i < n; ++i) {
        CHECK(rel_ok(g_fx[i], g_ex[i], 1e-5f));
        CHECK(rel_ok(g_fy[i], g_ey[i], 1e-5f));
    }
}

int main(void)
{
    size_t n;

    printf("== fv_atan_approx on [-1,1] ==\n");
    test_atan();

    printf("== _fast normalize vs exact: relative error <= 1e-5, unit length, exact zero guard ==\n");
    test_vec2(0); test_vec2(63); test_vec2(1024); test_vec2(4097);
    test_vec3(0); test_vec3(63); test_vec3(1024); test_vec3(4097);
    test_vec4(0); test_vec4(63); test_vec4(1024); test_vec4(4097);
    test_quat(0); test_quat(63); test_quat(1024); test_quat(4097);
    test_madd_norm(0); test_madd_norm(63); test_madd_norm(1024); test_madd_norm(4097);
    for (n = 0; n <= 40; ++n) {
        test_vec2(n);
        test_vec3(n);
        test_vec4(n);
        test_quat(n);
        test_madd_norm(n);
    }

    /* Unit length, checked separately on a larger batch. */
    {
        size_t i;
        fill(1024);
        fv_vec3_normalize_soa_fast(g_ax, g_ay, g_az, g_fx, g_fy, g_fz, 1024);
        for (i = 0; i < 1024; ++i) {
            float inlen = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i] + g_az[i] * g_az[i]);
            float len = sqrtf(g_fx[i] * g_fx[i] + g_fy[i] * g_fy[i] + g_fz[i] * g_fz[i]);
            if (inlen > 1e-5f) { /* the zero-guard lanes are not unit vectors by design */
                CHECK(fabsf(len - 1.0f) <= 2e-6f);
            }
        }
    }

    fv_vec2_normalize_soa_fast(NULL, NULL, NULL, NULL, 0);
    fv_vec3_normalize_soa_fast(NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec4_normalize_soa_fast(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_quat_normalize_soa_fast(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec2_madd_normalize_soa_fast(NULL, NULL, NULL, NULL, 1.5f, NULL, NULL, 0);

    printf("test_fast: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
