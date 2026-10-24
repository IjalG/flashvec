/* C translation-unit test for the fused pipeline entry points fv_vec2_madd_normalize_soa / _stream.
 *
 * Same criterion as tests/test_madd_norm.cpp, in a plain C11 caller: the fused output is memcmp-compared
 * against the three-call sequence fv_vec2_scale_soa -> fv_vec2_add_soa -> fv_vec2_normalize_soa, and the
 * _stream form against the normal fused form. Not a tolerance check -- one differing bit fails.
 *
 * The corpora stay at one NaN bit pattern (or none), which is what the library's bit-exactness promise
 * covers; the anything-goes corpus lives in the C++ test, where the payload difference is counted. So the
 * two tests divide the work: this one proves the C form of the include is strict, the C++ one also quantifies
 * the documented NaN-payload boundary.
 *
 * Sections: n = 0..40 (every tail length) x four corpora; s sweep; large n = 1023 / 4096 / 8192 / 12345;
 * the 8 output/input phase shifts; the _stream phase-mismatch fallback; n == 0 with NULL.
 *
 * Run: make test-maddnorm-c
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

#define MAXN 12352u
#define ROOM 32u
#define POISON 0x5Au

static float g_ax[MAXN + ROOM], g_ay[MAXN + ROOM];
static float g_bx[MAXN + ROOM], g_by[MAXN + ROOM];
static float g_tx[MAXN + ROOM], g_ty[MAXN + ROOM];
static float g_ux[MAXN + ROOM], g_uy[MAXN + ROOM];
static float g_rx[MAXN + ROOM], g_ry[MAXN + ROOM];
static float g_fx[MAXN + ROOM], g_fy[MAXN + ROOM];
static float g_sx[MAXN + ROOM], g_sy[MAXN + ROOM];

enum { CORP_FINITE, CORP_INF, CORP_NAN1, CORP_ZERO_MAKE, CORP_COUNT };

static uint32_t g_rng = 0x1234ABCDu;

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
    0x00000000u, 0x80000000u, /* +-0.0 */
    0x3f800000u, 0xbf800000u, /* +-1.0 */
    0x40600000u, 0xc0100000u, /* 3.5, -2.25 */
    0x00800000u, 0x80800000u, /* +-FLT_MIN */
    0x7f7fffffu, 0xff7fffffu, /* +-FLT_MAX */
    0x7149f2cau, 0xf149f2cau, /* +-1e30 */
    0x0da24260u, 0x8da24260u, /* +-1e-30 */
    0x00000001u               /* smallest denormal */
};
#define NFINITE (sizeof kFinite / sizeof kFinite[0])

static const uint32_t kInf[] = { 0x7f800000u, 0xff800000u };
#define NINF (sizeof kInf / sizeof kInf[0])

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
        case CORP_INF:
            g_ax[i] = (rnd() & 3u) == 0u ? from_bits(kInf[rnd() % NINF]) : rnd_unit();
            g_ay[i] = (rnd() & 3u) == 0u ? from_bits(kInf[rnd() % NINF]) : rnd_unit();
            g_bx[i] = (rnd() & 3u) == 0u ? from_bits(kInf[rnd() % NINF]) : rnd_unit();
            g_by[i] = (rnd() & 3u) == 0u ? from_bits(kInf[rnd() % NINF]) : rnd_unit();
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

static int same(const float *a, const float *b, size_t n)
{
    return memcmp((const void *)a, (const void *)b, n * sizeof(float)) == 0;
}

static int poison_tail(const float *p, size_t n)
{
    const unsigned char *b = (const unsigned char *)(const void *)p;
    size_t i;
    const size_t to = n + 16u < MAXN + ROOM ? n + 16u : MAXN + ROOM;
    for (i = n * sizeof(float); i < to * sizeof(float); ++i) {
        if (b[i] != POISON) {
            return 0;
        }
    }
    return 1;
}

static void run_case(size_t n, int corpus, float s)
{
    fill(n, corpus, s);
    memset(g_rx, POISON, sizeof g_rx);
    memset(g_ry, POISON, sizeof g_ry);
    memset(g_fx, POISON, sizeof g_fx);
    memset(g_fy, POISON, sizeof g_fy);
    memset(g_sx, POISON, sizeof g_sx);
    memset(g_sy, POISON, sizeof g_sy);

    fv_vec2_scale_soa(g_ax, g_ay, s, g_tx, g_ty, n);
    fv_vec2_add_soa(g_tx, g_ty, g_bx, g_by, g_ux, g_uy, n);
    fv_vec2_normalize_soa(g_ux, g_uy, g_rx, g_ry, n);

    fv_vec2_madd_normalize_soa(g_ax, g_ay, g_bx, g_by, s, g_fx, g_fy, n);
    fv_vec2_madd_normalize_soa_stream(g_ax, g_ay, g_bx, g_by, s, g_sx, g_sy, n);

    CHECK(same(g_fx, g_rx, n));
    CHECK(same(g_fy, g_ry, n));
    CHECK(same(g_sx, g_rx, n));
    CHECK(same(g_sy, g_ry, n));
    CHECK(poison_tail(g_fx, n) && poison_tail(g_fy, n));
    CHECK(poison_tail(g_sx, n) && poison_tail(g_sy, n));
}

int main(void)
{
    static const float sv[] = { 1.5f, -1.0f, 0.25f, 0.0f };
    static const size_t sizes[] = { 1023u, 4096u, 8192u, 12345u };
    size_t corpus;
    size_t k;
    size_t n;
    size_t off_in;
    size_t off_out;

    for (corpus = 0; corpus < (size_t)CORP_COUNT; ++corpus) {
        for (n = 0; n <= 40; ++n) {
            run_case(n, (int)corpus, sv[n % (sizeof sv / sizeof sv[0])]);
        }
        for (k = 0; k < sizeof sizes / sizeof sizes[0]; ++k) {
            run_case(sizes[k], (int)corpus, 1.25f);
        }
    }

    /* output / input phase shifts on one corpus, including two outputs whose 32 B phases differ */
    fill(MAXN, CORP_NAN1, 1.5f);
    for (off_in = 0; off_in < 8u; ++off_in) {
        fv_vec2_scale_soa(g_ax + off_in, g_ay + off_in, 1.5f, g_tx, g_ty, 997u);
        fv_vec2_add_soa(g_tx, g_ty, g_bx + off_in, g_by + off_in, g_ux, g_uy, 997u);
        fv_vec2_normalize_soa(g_ux, g_uy, g_rx, g_ry, 997u);
        for (off_out = 0; off_out < 8u; ++off_out) {
            fv_vec2_madd_normalize_soa(g_ax + off_in, g_ay + off_in, g_bx + off_in, g_by + off_in,
                                       1.5f, g_fx + off_out, g_fy + off_out, 997u);
            CHECK(same(g_fx + off_out, g_rx, 997u));
            CHECK(same(g_fy + off_out, g_ry, 997u));
            fv_vec2_madd_normalize_soa_stream(g_ax + off_in, g_ay + off_in, g_bx + off_in,
                                              g_by + off_in, 1.5f, g_sx + off_out,
                                              g_sy + off_out, 997u);
            CHECK(same(g_sx + off_out, g_rx, 997u));
            CHECK(same(g_sy + off_out, g_ry, 997u));
        }
        fv_vec2_madd_normalize_soa_stream(g_ax + off_in, g_ay + off_in, g_bx + off_in,
                                          g_by + off_in, 1.5f, g_sx, g_sy + 1u, 997u);
        CHECK(same(g_sx, g_rx, 997u));
        CHECK(same(g_sy + 1u, g_ry, 997u));
    }

    /* n == 0 with NULL: no pointer is touched */
    fv_vec2_scale_soa(NULL, NULL, 2.0f, NULL, NULL, 0u);
    fv_vec2_add_soa(NULL, NULL, NULL, NULL, NULL, NULL, 0u);
    fv_vec2_normalize_soa(NULL, NULL, NULL, NULL, 0u);
    fv_vec2_madd_normalize_soa(NULL, NULL, NULL, NULL, 2.0f, NULL, NULL, 0u);
    fv_vec2_madd_normalize_soa_stream(NULL, NULL, NULL, NULL, 2.0f, NULL, NULL, 0u);
    CHECK(1);

    printf("test_madd_norm_c: C API ok: fused = scale->add->normalize (memcmp), _stream = fused  "
           "(%d checks, %d failures)\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
