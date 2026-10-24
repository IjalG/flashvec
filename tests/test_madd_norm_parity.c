/* Cross-form bit-level parity stream for the fused pipeline entry points
 * fv_vec2_madd_normalize_soa / _stream (same approach as tests/test_parity.c and tests/test_mat4_parity.c).
 *
 * One source file is compiled six ways -- C and C++ x inline form and FV_USE_EXTERN, the pure-C fallback
 * (-DFV_I_NO_INTRINSICS), and a hostile-flags consumer (-march=native -mfma -ffp-contract=fast) -- and every
 * build writes its cases as one binary record stream; the streams are then compared byte for byte. Because the
 * source is identical, the inputs, the call order and the arguments are identical too.
 *
 * Why this fixture has to exist next to tests/test_madd_norm.cpp: the pure-C fallback has no AVX2 kernel at all
 * (everything runs the scalar tail), and under -ffp-contract=fast only the header's pragma keeps the consumer's
 * -mfma from turning "a * s + b" into vfmadd. If either of those changes a single bit, two streams differ.
 *
 * Record format (native byte order; three u32 header, then the raw payload):
 *   [case_id][n][nbytes] followed by the payload: (n + GUARD) floats of ox and (n + GUARD) floats of oy, for both
 *   the normal and the streaming entry. Every output buffer is memset to 0x5A before the call, so both an
 *   out-of-range write and a missed write make two files differ.
 *
 * The corpora keep at most one NaN bit pattern (the case the bit-exactness promise covers, see
 * include/flashvec.h); the NaN-payload boundary is quantified in tests/test_madd_norm.cpp instead.
 *
 * Usage: test_madd_norm_parity <out.bin>   (called by make test-maddnorm-inline / test-maddnorm-fma)
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CAPL  4104u
#define GUARD 8u
#define POISON 0x5Au

static float g_ax[CAPL + GUARD], g_ay[CAPL + GUARD];
static float g_bx[CAPL + GUARD], g_by[CAPL + GUARD];
static float g_ox[CAPL + GUARD], g_oy[CAPL + GUARD];
static float g_sx[CAPL + GUARD], g_sy[CAPL + GUARD];

static unsigned int g_cases;
static unsigned int g_failures;
static uint32_t g_rng = 0xC0FFEE11u;

enum { CORP_FINITE, CORP_NAN1, CORP_ZERO_MAKE };

static uint32_t rnd(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static float rnd_unit(void)
{
    return (float)((double)(rnd() % 2000001u) / 1000000.0 - 1.0) * 2.0f;
}

static float from_bits(uint32_t u)
{
    float f;
    memcpy(&f, &u, sizeof f);
    return f;
}

static const uint32_t kFinite[] = {
    0x00000000u, 0x80000000u, 0x3f800000u, 0xbf800000u, 0x40600000u, 0xc0100000u,
    0x00800000u, 0x80800000u, 0x7f7fffffu, 0xff7fffffu, 0x7149f2cau, 0xf149f2cau,
    0x0da24260u, 0x8da24260u, 0x00000001u
};
#define NFINITE (sizeof kFinite / sizeof kFinite[0])

#define NAN_PATTERN 0x7fc00000u

static void fill(size_t n, int corpus, float s)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        switch (corpus) {
        case CORP_FINITE:
            g_ax[i] = (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
            g_ay[i] = (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
            g_bx[i] = (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
            g_by[i] = (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
            break;
        case CORP_NAN1:
            g_ax[i] = from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
            g_ay[i] = from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
            g_bx[i] = from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
            g_by[i] = from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
            break;
        case CORP_ZERO_MAKE:
        default:
            g_ax[i] = rnd_unit();
            g_ay[i] = rnd_unit();
            g_bx[i] = -(g_ax[i] * s);
            g_by[i] = -(g_ay[i] * s);
            break;
        }
    }
}

static void put_u32(FILE *f, uint32_t v)
{
    if (fwrite(&v, sizeof v, 1u, f) != 1u) {
        ++g_failures;
    }
}

static void put_arr(FILE *f, const float *p, size_t n)
{
    if (fwrite((const void *)p, sizeof(float), n, f) != n) {
        ++g_failures;
    }
}

/* One case: both fused forms on the same span, written side by side so a difference in either shows up. */
static void run_case(FILE *f, size_t n, int corpus, float s, size_t in_off, size_t out_off,
                     uint32_t id)
{
    const size_t len = n + GUARD;

    fill(n + in_off, corpus, s);
    memset(g_ox, POISON, sizeof g_ox);
    memset(g_oy, POISON, sizeof g_oy);
    memset(g_sx, POISON, sizeof g_sx);
    memset(g_sy, POISON, sizeof g_sy);

    fv_vec2_madd_normalize_soa(g_ax + in_off, g_ay + in_off, g_bx + in_off, g_by + in_off, s,
                               g_ox + out_off, g_oy + out_off, n);
    fv_vec2_madd_normalize_soa_stream(g_ax + in_off, g_ay + in_off, g_bx + in_off, g_by + in_off,
                                      s, g_sx + out_off, g_sy + out_off, n);

    put_u32(f, id);
    put_u32(f, (uint32_t)n);
    put_u32(f, (uint32_t)(len * sizeof(float)));
    put_arr(f, g_ox, len);
    put_arr(f, g_oy, len);
    put_arr(f, g_sx, len);
    put_arr(f, g_sy, len);
    ++g_cases;
}

/* The same span with the two outputs at different phases inside 32 B: the shared-i peel cannot align both,
 * so the streaming kernel falls back to normal storeu. Still bit-identical, still part of the stream. */
static void run_phase_mismatch(FILE *f, size_t n, float s, uint32_t id)
{
    const size_t len = n + GUARD;
    size_t i;

    fill(n, CORP_NAN1, s);
    memset(g_sx, POISON, sizeof g_sx);
    memset(g_sy, POISON, sizeof g_sy);
    memset(g_ox, POISON, sizeof g_ox);
    memset(g_oy, POISON, sizeof g_oy);

    fv_vec2_madd_normalize_soa(g_ax, g_ay, g_bx, g_by, s, g_ox, g_oy, n);
    fv_vec2_madd_normalize_soa_stream(g_ax, g_ay, g_bx, g_by, s, g_sx, g_sy + 1u, n);

    put_u32(f, id);
    put_u32(f, (uint32_t)n);
    put_u32(f, (uint32_t)(len * sizeof(float)));
    put_arr(f, g_ox, len);
    put_arr(f, g_oy, len);
    /* the streaming y output sits one float further on: shift it back so the record shape stays the same */
    for (i = 0; i < len; ++i) {
        g_sy[i] = g_sy[i + 1u];
    }
    put_arr(f, g_sx, len);
    put_arr(f, g_sy, len);
    ++g_cases;
}

int main(int argc, char **argv)
{
    static const size_t tails[] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u, 12u,
                                    13u, 14u, 15u, 16u, 17u, 18u, 19u, 20u, 21u, 22u, 23u, 24u,
                                    63u, 1023u, 4096u };
    static const size_t phase_tails[] = { 0u, 7u, 8u, 333u };
    static const size_t offs[] = { 0u, 1u, 3u, 7u };
    FILE *f;
    uint32_t id = 0u;
    size_t k;
    size_t c;
    size_t a;
    size_t b;

    if (argc < 2) {
        printf("usage: %s <out.bin>\n", argv[0]);
        return 2;
    }
    f = fopen(argv[1], "wb");
    if (f == NULL) {
        printf("cannot open %s\n", argv[1]);
        return 2;
    }

    for (c = 0; c < 3u; ++c) {
        for (k = 0; k < sizeof tails / sizeof tails[0]; ++k) {
            run_case(f, tails[k], (int)c, 1.5f, 0u, 0u, id++);
        }
    }

    /* s sweep over one span: includes +0 and -0, where the whole output is the zero-vector mask */
    {
        static const float sv[] = { 1.0f, -1.0f, 3.25f, 1e30f, 1e-30f, 0.0f, -0.0f };
        size_t j;
        for (j = 0; j < sizeof sv / sizeof sv[0]; ++j) {
            run_case(f, 41u, CORP_FINITE, sv[j], 0u, 0u, id++);
            run_case(f, 41u, CORP_NAN1, sv[j], 0u, 0u, id++);
            run_case(f, 40u, CORP_ZERO_MAKE, sv[j], 0u, 0u, id++);
        }
    }

    /* pointer phases: the streaming kernel peels its output, and both outputs may sit at different phases */
    for (a = 0; a < sizeof offs / sizeof offs[0]; ++a) {
        for (b = 0; b < sizeof offs / sizeof offs[0]; ++b) {
            for (k = 0; k < sizeof phase_tails / sizeof phase_tails[0]; ++k) {
                run_case(f, phase_tails[k], CORP_NAN1, 1.5f, offs[a], offs[b], id++);
            }
        }
    }
    run_phase_mismatch(f, 333u, 1.5f, id++);
    run_phase_mismatch(f, 8u, -2.0f, id++);

    if (fclose(f) != 0) {
        ++g_failures;
    }
    printf("madd_norm_parity: cases=%u failures=%u\n", g_cases, g_failures);
    return g_failures == 0 ? 0 : 1;
}
