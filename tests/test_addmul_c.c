/* C caller acceptance for the integrator primitives (P4 group).
 * The C translation unit is the one that cannot take intrinsics (see the file header of include/flashvec.h): it sees the
 * header's dispatch and the AVX2 target attribute, and the per-value functions come from libflashvec.a. This fixture
 * proves the group is callable from plain C with the same bit-level result as the hand-written scalar expression, and it
 * covers the two shapes the C++ fixture reaches only indirectly: every tail length n = 0..17 and a misaligned in-place
 * destination.
 *
 * Usage: test_addmul_c        (make test-addmul-c)
 */
#include "flashvec.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CAPL 512u
#define GUARD 8u

static fv_vec2 g_x2[CAPL + GUARD], g_v2[CAPL + GUARD], g_a2[CAPL + GUARD], g_b2[CAPL + GUARD];
static fv_vec2 g_c2[CAPL + GUARD], g_o2[CAPL + GUARD], g_r2[CAPL + GUARD], g_k2[CAPL + GUARD];
static fv_vec4 g_x4[CAPL + GUARD], g_v4[CAPL + GUARD], g_r4[CAPL + GUARD];
static float g_sx[CAPL + GUARD], g_sy[CAPL + GUARD], g_sz[CAPL + GUARD];
static float g_tx[CAPL + GUARD], g_ty[CAPL + GUARD], g_tz[CAPL + GUARD];
static float g_ux[CAPL + GUARD], g_uy[CAPL + GUARD], g_uz[CAPL + GUARD];

static unsigned long g_checks;
static unsigned long g_failures;
static uint32_t g_rng = 0x2545F491u;

static uint32_t rnd(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static uint32_t bits(float f)
{
    uint32_t u;
    memcpy(&u, &f, sizeof u);
    return u;
}

static float rndv(void)
{
    return (float)((double)(rnd() % 2000001u) / 1000000.0 - 1.0) * 3.0f;
}

static void fill_vec2(fv_vec2 *dst, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i].x = rndv();
        dst[i].y = rndv();
        if ((rnd() & 31u) == 0u) {
            dst[i].x = 0.0f; /* +0.0 / -0.0 boundary */
            dst[i].y = -0.0f;
        }
    }
}

static void fill_vec4(fv_vec4 *dst, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i].x = rndv();
        dst[i].y = rndv();
        dst[i].z = rndv();
        dst[i].w = rndv();
        if ((rnd() & 31u) == 0u) {
            dst[i].z = 0.0f;
            dst[i].w = -0.0f;
        }
    }
}

static void fill_f(float *dst, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i] = rndv();
    }
}

static void cmp_f(const float *got, const float *want, size_t n, const char *what)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        ++g_checks;
        if (bits(got[i]) != bits(want[i])) {
            if (g_failures < 8) {
                printf("  FAIL %s [%zu] got %08x want %08x\n", what, i, bits(got[i]), bits(want[i]));
            }
            ++g_failures;
        }
    }
}

int main(void)
{
    size_t n, ph;
    float s;

    printf("test_addmul_c: 1 - every tail length n = 0..17, three scales, in-place and maddv\n");
    for (n = 0; n <= 17; ++n) {
        static const float kS[] = { 1.0f, -0.5f, 1.0e20f };
        size_t k;
        for (k = 0; k < sizeof kS / sizeof kS[0]; ++k) {
            size_t i;
            s = kS[k];
            fill_vec2(g_x2, n);
            fill_vec2(g_v2, n);
            for (i = 0; i < n; ++i) {
                g_r2[i].x = g_x2[i].x + g_v2[i].x * s;
                g_r2[i].y = g_x2[i].y + g_v2[i].y * s;
            }
            fv_vec2_addmul_x4(g_x2, g_v2, s, n);
            for (i = 0; i < n; ++i) {
                ++g_checks;
                if (bits(g_x2[i].x) != bits(g_r2[i].x) || bits(g_x2[i].y) != bits(g_r2[i].y)) {
                    printf("  FAIL addmul_x4 (C) n=%zu [%zu]\n", n, i);
                    ++g_failures;
                }
            }

            fill_vec4(g_x4, n);
            fill_vec4(g_v4, n);
            for (i = 0; i < n; ++i) {
                g_r4[i].x = g_x4[i].x + g_v4[i].x * s;
                g_r4[i].y = g_x4[i].y + g_v4[i].y * s;
                g_r4[i].z = g_x4[i].z + g_v4[i].z * s;
                g_r4[i].w = g_x4[i].w + g_v4[i].w * s;
            }
            fv_vec4_addmul_x2(g_x4, g_v4, s, n);
            for (i = 0; i < n; ++i) {
                ++g_checks;
                if (bits(g_x4[i].x) != bits(g_r4[i].x) || bits(g_x4[i].y) != bits(g_r4[i].y) ||
                    bits(g_x4[i].z) != bits(g_r4[i].z) || bits(g_x4[i].w) != bits(g_r4[i].w)) {
                    printf("  FAIL addmul_x2 (C) n=%zu [%zu]\n", n, i);
                    ++g_failures;
                }
            }

            fill_vec2(g_a2, n);
            fill_vec2(g_b2, n);
            fill_vec2(g_c2, n);
            for (i = 0; i < n; ++i) {
                g_r2[i].x = g_a2[i].x * g_b2[i].x + g_c2[i].x;
                g_r2[i].y = g_a2[i].y * g_b2[i].y + g_c2[i].y;
            }
            memset(g_o2, 0x5A, (n + GUARD) * sizeof *g_o2);
            fv_vec2_maddv_x4(g_a2, g_b2, g_c2, g_o2, n);
            memset(g_k2, 0x5A, (n + GUARD) * sizeof *g_k2);
            fv_vec2_maddv_x4_stream(g_a2, g_b2, g_c2, g_k2, n);
            for (i = 0; i < n; ++i) {
                ++g_checks;
                if (bits(g_o2[i].x) != bits(g_r2[i].x) || bits(g_o2[i].y) != bits(g_r2[i].y)) {
                    printf("  FAIL maddv_x4 (C) n=%zu [%zu]\n", n, i);
                    ++g_failures;
                }
                ++g_checks;
                if (bits(g_k2[i].x) != bits(g_r2[i].x) || bits(g_k2[i].y) != bits(g_r2[i].y)) {
                    printf("  FAIL maddv_x4_stream (C) n=%zu [%zu]\n", n, i);
                    ++g_failures;
                }
            }
        }
    }

    printf("test_addmul_c: 2 - SoA in place (vec2 / vec3) vs the component expressions\n");
    s = 0.375f;
    fill_f(g_sx, 300);
    fill_f(g_sy, 300);
    fill_f(g_sz, 300);
    fill_f(g_tx, 300);
    fill_f(g_ty, 300);
    fill_f(g_tz, 300);
    {
        size_t i;
        for (i = 0; i < 300; ++i) {
            g_ux[i] = g_sx[i] + g_tx[i] * s;
            g_uy[i] = g_sy[i] + g_ty[i] * s;
        }
    }
    fv_vec2_addmul_soa(g_sx, g_sy, g_tx, g_ty, s, 300);
    cmp_f(g_sx, g_ux, 300, "addmul_soa x (C)");
    cmp_f(g_sy, g_uy, 300, "addmul_soa y (C)");
    {
        size_t i;
        for (i = 0; i < 300; ++i) {
            g_ux[i] = g_sx[i] + g_tx[i] * s;
            g_uy[i] = g_sy[i] + g_ty[i] * s;
            g_uz[i] = g_sz[i] + g_tz[i] * s;
        }
    }
    fv_vec3_addmul_soa(g_sx, g_sy, g_sz, g_tx, g_ty, g_tz, s, 300);
    cmp_f(g_sx, g_ux, 300, "addmul_v3_soa x (C)");
    cmp_f(g_sy, g_uy, 300, "addmul_v3_soa y (C)");
    cmp_f(g_sz, g_uz, 300, "addmul_v3_soa z (C)");

    printf("test_addmul_c: 3 - misaligned in-place destination (phases 1..7)\n");
    for (ph = 1; ph <= 7; ++ph) {
        size_t i;
        size_t nn = 129;
        fill_vec2(g_x2 + ph, nn);
        fill_vec2(g_v2 + ph, nn);
        memcpy(g_k2, g_x2 + ph, nn * sizeof *g_k2);
        for (i = 0; i < nn; ++i) {
            g_r2[i].x = g_k2[i].x + g_v2[ph + i].x * 2.0f;
            g_r2[i].y = g_k2[i].y + g_v2[ph + i].y * 2.0f;
        }
        fv_vec2_addmul_x4(g_x2 + ph, g_v2 + ph, 2.0f, nn);
        for (i = 0; i < nn; ++i) {
            ++g_checks;
            if (bits(g_x2[ph + i].x) != bits(g_r2[i].x) || bits(g_x2[ph + i].y) != bits(g_r2[i].y)) {
                printf("  FAIL misaligned addmul_x4 (C) phase=%zu [%zu]\n", ph, i);
                ++g_failures;
            }
        }
    }

    printf("test_addmul_c: 4 - n == 0 with NULL, and the per-value forms\n");
    {
        fv_vec2 *nx2 = (fv_vec2 *)0;
        const fv_vec2 *nv2 = (const fv_vec2 *)0;
        float *fx = (float *)0;
        float *fy = (float *)0;
        const float *cx = (const float *)0;
        const float *cy = (const float *)0;
        fv_vec2_maddv_x4(nv2, nv2, nv2, nx2, 0);
        fv_vec2_maddv_x4_stream(nv2, nv2, nv2, nx2, 0);
        fv_vec2_addmul_x4(nx2, nv2, 1.0f, 0);
        fv_vec2_addmul_soa(fx, fy, cx, cy, 1.0f, 0);
        ++g_checks;
    }
    {
        fv_vec2 x2, v2;
        fv_vec3 x3, v3;
        fv_vec4 x4, v4;
        x2.x = 1.5f; x2.y = -2.25f;
        v2.x = 0.5f; v2.y = 4.0f;
        {
            fv_vec2 r = fv_vec2_addmul(x2, v2, 3.0f);
            ++g_checks;
            if (bits(r.x) != bits(x2.x + v2.x * 3.0f) || bits(r.y) != bits(x2.y + v2.y * 3.0f)) {
                printf("  FAIL fv_vec2_addmul (C)\n");
                ++g_failures;
            }
        }
        x3.x = 1.0f; x3.y = 2.0f; x3.z = -3.0f;
        v3.x = 0.25f; v3.y = -0.5f; v3.z = 8.0f;
        {
            fv_vec3 r = fv_vec3_addmul(x3, v3, 1.5f);
            ++g_checks;
            if (bits(r.x) != bits(x3.x + v3.x * 1.5f) || bits(r.y) != bits(x3.y + v3.y * 1.5f) ||
                bits(r.z) != bits(x3.z + v3.z * 1.5f)) {
                printf("  FAIL fv_vec3_addmul (C)\n");
                ++g_failures;
            }
        }
        x4.x = 1.0f; x4.y = 2.0f; x4.z = 3.0f; x4.w = 4.0f;
        v4.x = -1.0f; v4.y = 0.5f; v4.z = 0.0f; v4.w = -0.0f;
        {
            fv_vec4 r = fv_vec4_addmul(x4, v4, 2.0f);
            ++g_checks;
            if (bits(r.x) != bits(x4.x + v4.x * 2.0f) || bits(r.y) != bits(x4.y + v4.y * 2.0f) ||
                bits(r.z) != bits(x4.z + v4.z * 2.0f) || bits(r.w) != bits(x4.w + v4.w * 2.0f)) {
                printf("  FAIL fv_vec4_addmul (C)\n");
                ++g_failures;
            }
        }
        {
            fv_vec2 r = fv_vec2_maddv(x2, v2, x2);
            ++g_checks;
            if (bits(r.x) != bits(x2.x * v2.x + x2.x) || bits(r.y) != bits(x2.y * v2.y + x2.y)) {
                printf("  FAIL fv_vec2_maddv (C)\n");
                ++g_failures;
            }
        }
    }

    printf("test_addmul_c: C API ok (%lu checks, %lu failures)\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
