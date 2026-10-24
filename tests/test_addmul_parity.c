/* Cross-form bit-level parity stream for the integrator primitives (P4 group):
 *   fv_vec2_addmul_x4 / fv_vec2_addmul_soa / fv_vec3_addmul_soa / fv_vec4_addmul_x2 /
 *   fv_vec2_maddv_x4 / fv_vec2_maddv_x4_stream   plus the per-value fv_vec2_addmul / fv_vec3_addmul /
 *   fv_vec4_addmul / fv_vec2_maddv.
 *
 * Same approach as tests/test_parity.c, tests/test_mat4_parity.c and tests/test_madd_norm_parity.c: one source file is
 * compiled six ways -- C and C++ x inline form and FV_USE_EXTERN, the pure-C fallback (-DFV_I_NO_INTRINSICS), and a
 * hostile consumer (-march=native -mfma -ffp-contract=fast) -- and every build writes its cases as one binary record
 * stream; the streams are then compared byte for byte. Because the source is identical, the inputs, the call order and
 * the arguments are identical too.
 *
 * Why this fixture has to exist: the pure-C fallback runs no AVX2 kernel at all (every element takes the scalar tail of
 * the dispatch), and under -ffp-contract=fast only the header's pragma keeps the consumer's -mfma from turning
 * "x + v*s" and "a*b + c" into vfmadd. If either of those changes a single bit, two streams differ.
 *
 * Record format (native byte order; three u32 header, then the raw payload):
 *   [case_id][n][nbytes] followed by the payload: the destination buffers of every entry at that case, each written with
 *   its GUARD poisoned elements included, so both an out-of-range write and a missed write make two files differ.
 *   The per-value API is recorded once at the end (four fv_vec2 / fv_vec3 / fv_vec4 results, raw bytes).
 *
 * The corpora keep at most one NaN bit pattern (the case the bit-exactness promise covers, see include/flashvec.h).
 *
 * Usage: test_addmul_parity <out.bin>   (called by make test-addmul-inline / test-addmul-fma)
 */
#include "flashvec.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CAPL  4104u
#define GUARD 8u
#define POISON 0x5Au

static fv_vec2 g_x2[CAPL + GUARD], g_v2[CAPL + GUARD], g_a2[CAPL + GUARD], g_b2[CAPL + GUARD];
static fv_vec2 g_c2[CAPL + GUARD], g_o2[CAPL + GUARD], g_k2[CAPL + GUARD];
static fv_vec4 g_x4[CAPL + GUARD], g_v4[CAPL + GUARD];
static float g_sx[CAPL + GUARD], g_sy[CAPL + GUARD], g_sz[CAPL + GUARD];
static float g_tx[CAPL + GUARD], g_ty[CAPL + GUARD], g_tz[CAPL + GUARD];

static unsigned int g_cases;
static unsigned int g_failures;
static uint32_t g_rng = 0x9E3779B9u;

enum { CORP_FINITE, CORP_NAN1 };

static uint32_t rnd(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static float from_bits(uint32_t u)
{
    float f;
    memcpy(&f, &u, sizeof f);
    return f;
}

static float rnd_unit(void)
{
    return (float)((double)(rnd() % 2000001u) / 1000000.0 - 1.0) * 2.0f;
}

static const uint32_t kFinite[] = {
    0x00000000u, 0x80000000u, 0x3f800000u, 0xbf800000u, 0x40600000u, 0xc0100000u,
    0x00800000u, 0x80800000u, 0x7f7fffffu, 0xff7fffffu, 0x7149f2cau, 0xf149f2cau,
    0x00000001u, 0x80000001u
};
#define NFINITE (sizeof kFinite / sizeof kFinite[0])

#define NAN_PATTERN 0x7fc00000u

static float gen(int corpus)
{
    if (corpus == CORP_NAN1) {
        return from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
    }
    return (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
}

static void fill_vec2(fv_vec2 *dst, size_t n, int corpus)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i].x = gen(corpus);
        dst[i].y = gen(corpus);
    }
}

static void fill_vec4(fv_vec4 *dst, size_t n, int corpus)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i].x = gen(corpus);
        dst[i].y = gen(corpus);
        dst[i].z = gen(corpus);
        dst[i].w = gen(corpus);
    }
}

static void fill_f(float *dst, size_t n, int corpus)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        dst[i] = gen(corpus);
    }
}

static FILE *g_out;

static void record(uint32_t id, size_t n, const void *payload, size_t nbytes)
{
    uint32_t hdr[3];
    hdr[0] = id;
    hdr[1] = (uint32_t)n;
    hdr[2] = (uint32_t)nbytes;
    if (g_out != NULL) {
        if (fwrite(hdr, sizeof hdr, 1, g_out) != 1 || fwrite(payload, 1, nbytes, g_out) != nbytes) {
            printf("addmul_parity: write error at case %u\n", id);
            ++g_failures;
        }
    }
}

int main(int argc, char **argv)
{
    static const size_t kSizes[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                     21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                     40, 63, 1023, 4097 };
    static const size_t kPhase[] = { 0, 1, 3, 5 };
    static const float kScale[] = { 1.0f, -2.5f, 0.125f };
    uint32_t id = 1;
    size_t i, j, c;
    float s;

    g_out = NULL;
    if (argc > 1) {
        g_out = fopen(argv[1], "wb");
        if (g_out == NULL) {
            printf("addmul_parity: cannot open %s\n", argv[1]);
            return 2;
        }
    }

    for (i = 0; i < sizeof kSizes / sizeof kSizes[0]; ++i) {
        size_t n = kSizes[i];
        for (j = 0; j < sizeof kPhase / sizeof kPhase[0]; ++j) {
            size_t ph = kPhase[j];
            for (c = 0; c < 2; ++c) {
                s = kScale[(i + j + c) % (sizeof kScale / sizeof kScale[0])];
                /* in-place vec2 */
                fill_vec2(g_x2 + ph, n, (int)c);
                fill_vec2(g_v2 + ph, n, (int)c);
                memset(g_x2 + ph + n, POISON, GUARD * sizeof *g_x2);
                fv_vec2_addmul_x4(g_x2 + ph, g_v2 + ph, s, n);
                record(id++, n, g_x2 + ph, (n + GUARD) * sizeof *g_x2);

                /* in-place SoA vec2, then vec3 stacked on the same arrays */
                fill_f(g_sx + ph, n, (int)c);
                fill_f(g_sy + ph, n, (int)c);
                fill_f(g_tx + ph, n, (int)c);
                fill_f(g_ty + ph, n, (int)c);
                memset(g_sx + ph + n, POISON, GUARD * sizeof *g_sx);
                memset(g_sy + ph + n, POISON, GUARD * sizeof *g_sy);
                fv_vec2_addmul_soa(g_sx + ph, g_sy + ph, g_tx + ph, g_ty + ph, s, n);
                record(id++, n, g_sx + ph, (n + GUARD) * sizeof *g_sx);
                record(id++, n, g_sy + ph, (n + GUARD) * sizeof *g_sy);

                fill_f(g_sz + ph, n, (int)c);
                fill_f(g_tz + ph, n, (int)c);
                memset(g_sz + ph + n, POISON, GUARD * sizeof *g_sz);
                fv_vec3_addmul_soa(g_sx + ph, g_sy + ph, g_sz + ph, g_tx + ph, g_ty + ph, g_tz + ph, s, n);
                record(id++, n, g_sx + ph, (n + GUARD) * sizeof *g_sx);
                record(id++, n, g_sy + ph, (n + GUARD) * sizeof *g_sy);
                record(id++, n, g_sz + ph, (n + GUARD) * sizeof *g_sz);

                /* in-place vec4 */
                fill_vec4(g_x4 + ph, n, (int)c);
                fill_vec4(g_v4 + ph, n, (int)c);
                memset(g_x4 + ph + n, POISON, GUARD * sizeof *g_x4);
                fv_vec4_addmul_x2(g_x4 + ph, g_v4 + ph, s, n);
                record(id++, n, g_x4 + ph, (n + GUARD) * sizeof *g_x4);

                /* maddv, normal and streaming, at the same output phase */
                fill_vec2(g_a2 + ph, n, (int)c);
                fill_vec2(g_b2 + ph, n, (int)c);
                fill_vec2(g_c2 + ph, n, (int)c);
                memset(g_o2 + ph, POISON, (n + GUARD) * sizeof *g_o2);
                fv_vec2_maddv_x4(g_a2 + ph, g_b2 + ph, g_c2 + ph, g_o2 + ph, n);
                record(id++, n, g_o2 + ph, (n + GUARD) * sizeof *g_o2);
                memset(g_k2 + ph, POISON, (n + GUARD) * sizeof *g_k2);
                fv_vec2_maddv_x4_stream(g_a2 + ph, g_b2 + ph, g_c2 + ph, g_k2 + ph, n);
                record(id++, n, g_k2 + ph, (n + GUARD) * sizeof *g_k2);

                ++g_cases;
            }
        }
    }

    /* Per-value API: a fixed set of inputs, recorded raw so every form must produce the same bits. */
    {
        float pv[16];
        memset(pv, 0, sizeof pv); /* every byte of the record must be defined, or two forms differ on stack garbage */
        fv_vec2 a2, b2, c2, r2;
        fv_vec3 a3, b3, r3;
        fv_vec4 a4, b4, r4;
        a2.x = 1.5f; a2.y = -0.0f;
        b2.x = -3.25f; b2.y = 2.0f;
        c2.x = 0.125f; c2.y = 1.0e-30f;
        r2 = fv_vec2_addmul(a2, b2, 0.75f);
        memcpy(pv + 0, &r2, sizeof r2);
        r2 = fv_vec2_maddv(a2, b2, c2);
        memcpy(pv + 2, &r2, sizeof r2);
        a3.x = 2.0f; a3.y = -4.0f; a3.z = 0.0f;
        b3.x = 0.5f; b3.y = 8.0f; b3.z = -1.0f;
        r3 = fv_vec3_addmul(a3, b3, -2.0f);
        memcpy(pv + 4, &r3, sizeof r3);
        a4.x = 1.0f; a4.y = -1.0f; a4.z = 3.0f; a4.w = -0.0f;
        b4.x = 7.0f; b4.y = 0.25f; b4.z = -2.0f; b4.w = 4.0f;
        r4 = fv_vec4_addmul(a4, b4, 3.0f);
        memcpy(pv + 8, &r4, sizeof r4);
        record(id++, 0, pv, sizeof pv);
        ++g_cases;
    }

    if (g_out != NULL) {
        fclose(g_out);
    }
    printf("addmul_parity: cases=%u failures=%u\n", g_cases, g_failures);
    return g_failures == 0 ? 0 : 1;
}
