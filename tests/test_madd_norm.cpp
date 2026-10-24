/* Bit-level test for the fused pipeline entry points fv_vec2_madd_normalize_soa / _stream.
 *
 * The criterion is the one the header contract states: the fused output must be memcmp-identical to
 * the three-call sequence
 *     fv_vec2_scale_soa(a, s -> t);  fv_vec2_add_soa(t, b -> u);  fv_vec2_normalize_soa(u -> out)
 * for every n, every tail length and every input, and the _stream form must equal the normal form.
 * One differing bit fails -- this is not a tolerance check.
 *
 * Why that is enough: those three entries carry the library's bit-exactness promise themselves (they are
 * covered by tests/test_batch.cpp and tests/test_parity.c), so nothing here repeats an error analysis. What
 * this test has to show is only that the composition keeps the evaluation order and the zero-vector mask of
 * those three entries, element by element.
 *
 * The one boundary, which the fused entry inherits rather than adds (header contract, "the one exception"
 * note in the matrix section): when **two different NaN bit patterns meet inside one arithmetic
 * instruction**, which NaN comes out is decided by that instruction's operand order, and the compiler is
 * free to allocate the operands of a commutative vector operation. Measured while writing this test:
 * fv_i_normalize_soa_avx2 and fv_i_madd_norm_soa_avx2 compile the same intrinsic expression
 * _mm256_mul_ps(x, inv) into vmulps with the operands in opposite order (objdump), so for such input the
 * two forms can differ in NaN payload/sign -- both are legal NaNs, not a numerical error. The corpus is
 * therefore grouped exactly like tests/test_mat4.cpp groups it:
 *   - STRICT corpora contain at most one NaN bit pattern (or no NaN at all): every byte must match;
 *   - one ANYTHING-GOES corpus may contain several NaN patterns: a difference is only tolerated when both
 *     sides are NaN in that component, and the tolerated count is printed, never swallowed.
 *
 * Sections:
 *   1 - n = 0..40: every tail length of the 8-elements-per-step kernel, five strict corpora
 *   2 - s sweep on one fixed length: 1.0 / -1.0 / 0.5 / 3.25 / 1e30 / 1e-30 / +0.0 / -0.0
 *   3 - large n: 1023 / 4096 / 4097 / 8192 / 12345 (the last two are not powers of two on purpose)
 *   4 - input and output phase shifts 0..7, plus two outputs whose phases differ (the _stream fallback)
 *   5 - path independence: an element's result must not depend on where it sits relative to the
 *       bulk/tail boundary
 *   6 - the anything-goes corpus (several NaN patterns): counted, NaN-vs-NaN only
 *   7 - n == 0 with NULL arguments touches nothing
 *
 * Run: make test-maddnorm
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0;
static int g_failures = 0;
static int g_tolerated = 0; /* section 6: NaN-vs-NaN differences, counted and printed */
static size_t g_compared = 0; /* component comparisons performed, so the tolerated count has a denominator */

#define CHECK(cond)                                                    \
    do {                                                               \
        ++g_checks;                                                    \
        if (!(cond)) {                                                 \
            ++g_failures;                                              \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);     \
        }                                                              \
    } while (0)

/* Longest span used here (large n) plus room for the pointer phase shifts. */
#define MAXN 20040u
#define ROOM 64u
#define POISON 0x5Au

static float g_ax[MAXN + ROOM], g_ay[MAXN + ROOM];
static float g_bx[MAXN + ROOM], g_by[MAXN + ROOM];
static float g_tx[MAXN + ROOM], g_ty[MAXN + ROOM]; /* scale -> add intermediates */
static float g_ux[MAXN + ROOM], g_uy[MAXN + ROOM];
static float g_rx[MAXN + ROOM], g_ry[MAXN + ROOM]; /* three-call reference */
static float g_fx[MAXN + ROOM], g_fy[MAXN + ROOM]; /* fused */
static float g_sx[MAXN + ROOM], g_sy[MAXN + ROOM]; /* fused, _stream */

enum {
    CORP_FINITE,      /* finite values only: no NaN input, no inf input -> strict */
    CORP_INF,         /* finite + -inf, still no NaN input (so only the hardware indefinite can appear) -> strict */
    CORP_NAN1,        /* finite + exactly one NaN bit pattern, no inf input -> strict */
    CORP_ZERO_MAKE,   /* b = -(a * s), so a * s + b is exactly +0: exercises the zero-vector guard -> strict */
    CORP_NEAR_ZERO,   /* values straddling len <= FV_EPS -> strict */
    CORP_ANY,         /* several NaN bit patterns + infinities: NaN-vs-NaN differences are counted */
    CORP_COUNT
};

static uint32_t g_rng = 0x9E3779B9u;

static uint32_t rnd(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/* [-2, 2), in steps small enough that products and sums stay representable */
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

static uint32_t to_bits(float f)
{
    uint32_t u;
    memcpy(&u, &f, sizeof u);
    return u;
}

static int is_nan_bits(uint32_t u)
{
    return (u & 0x7f800000u) == 0x7f800000u && (u & 0x007fffffu) != 0u;
}

/* Bit patterns, not decimal literals: FLT_MIN / FLT_MAX / the denormal / the NaN payloads must be exact. */
static const uint32_t kFinite[] = {
    0x00000000u, /* +0.0 */
    0x80000000u, /* -0.0 */
    0x3f800000u, /* 1.0 */
    0xbf800000u, /* -1.0 */
    0x40600000u, /* 3.5 */
    0xc0100000u, /* -2.25 */
    0x00800000u, /* FLT_MIN, the smallest normal */
    0x80800000u, /* -FLT_MIN */
    0x7f7fffffu, /* FLT_MAX */
    0xff7fffffu, /* -FLT_MAX */
    0x7149f2cau, /* 1e30 */
    0xf149f2cau, /* -1e30 */
    0x0da24260u, /* 1e-30 */
    0x8da24260u, /* -1e-30 */
    0x00000001u  /* the smallest denormal */
};
#define NFINITE (sizeof kFinite / sizeof kFinite[0])

static const uint32_t kInf[] = {
    0x7f800000u, /* +inf */
    0xff800000u  /* -inf */
};
#define NINF (sizeof kInf / sizeof kInf[0])

/* One quiet NaN pattern only. A second pattern appears only in CORP_ANY. */
#define NAN_PATTERN 0x7fc00000u

static const uint32_t kAny[] = {
    0x7f800000u, /* +inf */
    0xff800000u, /* -inf */
    0x7fc00000u, /* qNaN */
    0x7fc00001u, /* qNaN, another payload */
    0xffc00000u  /* -qNaN */
};
#define NANY (sizeof kAny / sizeof kAny[0])

/* Around FV_EPS (1e-6): some lanes fall below the zero-vector threshold, some just above. */
static const uint32_t kNear[] = {
    0x00000000u, 0x80000000u, /* +-0 */
    0x322bcc77u,              /* ~1e-8 */
    0xb22bcc77u,
    0x33d6bf95u,              /* ~1e-7 */
    0xb3d6bf95u,
    0x34000000u,              /* 1.19e-7 */
    0x358637bdu,              /* 1e-6 */
    0xb58637bdu,
    0x3a83126fu,              /* 1e-3 */
    0x00000001u               /* smallest denormal */
};
#define NNEAR (sizeof kNear / sizeof kNear[0])

static void poison_outputs(void)
{
    memset(g_tx, POISON, sizeof g_tx);
    memset(g_ty, POISON, sizeof g_ty);
    memset(g_ux, POISON, sizeof g_ux);
    memset(g_uy, POISON, sizeof g_uy);
    memset(g_rx, POISON, sizeof g_rx);
    memset(g_ry, POISON, sizeof g_ry);
    memset(g_fx, POISON, sizeof g_fx);
    memset(g_fy, POISON, sizeof g_fy);
    memset(g_sx, POISON, sizeof g_sx);
    memset(g_sy, POISON, sizeof g_sy);
}

static int still_poison(const float *p, size_t from, size_t to)
{
    const unsigned char *b = (const unsigned char *)(const void *)p;
    size_t i;
    for (i = from * sizeof(float); i < to * sizeof(float); ++i) {
        if (b[i] != POISON) {
            return 0;
        }
    }
    return 1;
}

/* No entry may write past n: the tail of every output buffer must still hold the poison pattern. */
static int guards_intact(size_t n)
{
    const size_t to = n + 64u < MAXN + ROOM ? n + 64u : MAXN + ROOM;
    return still_poison(g_rx, n, to) && still_poison(g_ry, n, to) &&
           still_poison(g_fx, n, to) && still_poison(g_fy, n, to) &&
           still_poison(g_sx, n, to) && still_poison(g_sy, n, to);
}

/* strict: every byte must match.
 * tolerant: a difference is accepted only when both sides are NaN in that component (see the file header),
 *           and every accepted difference is counted in g_tolerated. */
static int same_f32(const float *a, const float *b, size_t n, int tolerant)
{
    size_t i;
    g_compared += n;
    for (i = 0; i < n; ++i) {
        const uint32_t ua = to_bits(a[i]);
        const uint32_t ub = to_bits(b[i]);
        if (ua != ub) {
            if (!tolerant || !is_nan_bits(ua) || !is_nan_bits(ub)) {
                return 0;
            }
            ++g_tolerated;
        }
    }
    return 1;
}

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
            /* b = -(a * s) bit for bit, so a * s + b is exactly +0 in both components and the zero-vector
             * guard decides the lane (len == 0 <= FV_EPS). */
            g_ax[i] = rnd_unit();
            g_ay[i] = rnd_unit();
            g_bx[i] = -(g_ax[i] * s);
            g_by[i] = -(g_ay[i] * s);
            break;
        case CORP_NEAR_ZERO:
            g_ax[i] = from_bits(kNear[rnd() % NNEAR]);
            g_ay[i] = from_bits(kNear[rnd() % NNEAR]);
            g_bx[i] = from_bits(kNear[rnd() % NNEAR]);
            g_by[i] = from_bits(kNear[rnd() % NNEAR]);
            break;
        case CORP_ANY:
        default:
            g_ax[i] = from_bits((rnd() & 1u) != 0u ? kAny[rnd() % NANY] : kFinite[rnd() % NFINITE]);
            g_ay[i] = from_bits((rnd() & 1u) != 0u ? kAny[rnd() % NANY] : kFinite[rnd() % NFINITE]);
            g_bx[i] = from_bits((rnd() & 1u) != 0u ? kAny[rnd() % NANY] : kFinite[rnd() % NFINITE]);
            g_by[i] = from_bits((rnd() & 1u) != 0u ? kAny[rnd() % NANY] : kFinite[rnd() % NFINITE]);
            break;
        }
    }
}

/* Everything the contract says, for one (n, corpus, s): three-call reference first, then both fused forms. */
static void run_case(size_t n, int corpus, float s)
{
    const int tolerant = (corpus == CORP_ANY) ? 1 : 0;

    fill(n, corpus, s);
    poison_outputs();

    fv_vec2_scale_soa(g_ax, g_ay, s, g_tx, g_ty, n);
    fv_vec2_add_soa(g_tx, g_ty, g_bx, g_by, g_ux, g_uy, n);
    fv_vec2_normalize_soa(g_ux, g_uy, g_rx, g_ry, n);

    fv_vec2_madd_normalize_soa(g_ax, g_ay, g_bx, g_by, s, g_fx, g_fy, n);
    fv_vec2_madd_normalize_soa_stream(g_ax, g_ay, g_bx, g_by, s, g_sx, g_sy, n);

    CHECK(same_f32(g_fx, g_rx, n, tolerant));
    CHECK(same_f32(g_fy, g_ry, n, tolerant));
    CHECK(same_f32(g_sx, g_rx, n, tolerant));
    CHECK(same_f32(g_sy, g_ry, n, tolerant));
    CHECK(guards_intact(n));
}

static void run_tails(void)
{
    static const struct { int corpus; float s; } cfg[] = {
        { CORP_FINITE, 1.5f },
        { CORP_INF, 1.5f },
        { CORP_NAN1, 1.5f },
        { CORP_ZERO_MAKE, 2.0f },
        { CORP_NEAR_ZERO, 1.0f }
    };
    size_t k;
    size_t n;
    for (k = 0; k < sizeof cfg / sizeof cfg[0]; ++k) {
        for (n = 0; n <= 40; ++n) {
            run_case(n, cfg[k].corpus, cfg[k].s);
        }
    }
}

static void run_s_sweep(void)
{
    static const float sv[] = { 1.0f, -1.0f, 0.5f, 3.25f, 1e30f, 1e-30f, 0.0f, -0.0f };
    size_t k;
    for (k = 0; k < sizeof sv / sizeof sv[0]; ++k) {
        run_case(37, CORP_FINITE, sv[k]);
        run_case(37, CORP_INF, sv[k]);
        run_case(37, CORP_NAN1, sv[k]);
        run_case(40, CORP_ZERO_MAKE, sv[k]);
        run_case(37, CORP_ANY, sv[k]);
    }
}

static void run_large(void)
{
    static const size_t sizes[] = { 1023u, 4096u, 4097u, 8192u, 12345u };
    size_t k;
    int corpus;
    for (corpus = 0; corpus < CORP_COUNT; ++corpus) {
        for (k = 0; k < sizeof sizes / sizeof sizes[0]; ++k) {
            run_case(sizes[k], corpus, 1.25f);
        }
    }
}

/* Phase shifts: the input pointers may be misaligned, and the _stream kernel peels its output to 32 B.
 * The reference for a shifted input is simply the reference of the unshifted run read at the same index:
 * element i of the shifted arrays is element i + off_in of the original corpus. */
static void run_phase_shifts(void)
{
    const size_t n = 999u;
    size_t off_in;
    size_t off_out;
    int corpus;

    for (corpus = 0; corpus < 2; ++corpus) {
        const int c = corpus == 0 ? CORP_INF : CORP_ANY;
        const int tolerant = corpus == 0 ? 0 : 1;

        fill(MAXN, c, 1.5f);
        for (off_in = 0; off_in < 8u; ++off_in) {
            poison_outputs();
            fv_vec2_scale_soa(g_ax + off_in, g_ay + off_in, 1.5f, g_tx, g_ty, n);
            fv_vec2_add_soa(g_tx, g_ty, g_bx + off_in, g_by + off_in, g_ux, g_uy, n);
            fv_vec2_normalize_soa(g_ux, g_uy, g_rx, g_ry, n);

            for (off_out = 0; off_out < 8u; ++off_out) {
                fv_vec2_madd_normalize_soa(g_ax + off_in, g_ay + off_in, g_bx + off_in,
                                           g_by + off_in, 1.5f, g_fx + off_out, g_fy + off_out, n);
                CHECK(same_f32(g_fx + off_out, g_rx, n, tolerant));
                CHECK(same_f32(g_fy + off_out, g_ry, n, tolerant));

                fv_vec2_madd_normalize_soa_stream(g_ax + off_in, g_ay + off_in, g_bx + off_in,
                                                  g_by + off_in, 1.5f, g_sx + off_out,
                                                  g_sy + off_out, n);
                CHECK(same_f32(g_sx + off_out, g_rx, n, tolerant));
                CHECK(same_f32(g_sy + off_out, g_ry, n, tolerant));
            }

            /* Two outputs whose phase within 32 B differs: the shared-i peel cannot align both, so the
             * kernel must fall back to normal storeu. The result is bit-identical either way. */
            fv_vec2_madd_normalize_soa_stream(g_ax + off_in, g_ay + off_in, g_bx + off_in,
                                              g_by + off_in, 1.5f, g_sx, g_sy + 1u, n);
            CHECK(same_f32(g_sx, g_rx, n, tolerant));
            CHECK(same_f32(g_sy + 1u, g_ry, n, tolerant));
        }
    }
}

/* The same element must produce the same bits whether the 8-wide body or the scalar tail handles it. */
static void run_path_independence(void)
{
    static const size_t sizes[] = { 7u, 8u, 9u, 15u, 16u, 17u };
    size_t k;
    size_t idx;

    for (k = 0; k < sizeof sizes / sizeof sizes[0]; ++k) {
        const size_t len = sizes[k];
        fill(len, CORP_NAN1, 1.5f);
        poison_outputs();
        fv_vec2_madd_normalize_soa(g_ax, g_ay, g_bx, g_by, 1.5f, g_fx, g_fy, len);
        for (idx = 0; idx < len; ++idx) {
            const float one_x = g_fx[idx];
            const float one_y = g_fy[idx];
            fv_vec2_madd_normalize_soa(g_ax + idx, g_ay + idx, g_bx + idx, g_by + idx, 1.5f,
                                       g_sx, g_sy, 1u);
            CHECK(memcmp(&g_sx[0], &one_x, sizeof(float)) == 0);
            CHECK(memcmp(&g_sy[0], &one_y, sizeof(float)) == 0);
        }
    }
}

/* Several NaN bit patterns per element: only NaN-vs-NaN payload differences may appear, and they are printed. */
static void run_anything_goes(void)
{
    static const size_t sizes[] = { 8u, 9u, 17u, 64u, 1000u, 4097u };
    size_t k;
    for (k = 0; k < sizeof sizes / sizeof sizes[0]; ++k) {
        run_case(sizes[k], CORP_ANY, 1.5f);
    }
}

static void run_zero_length(void)
{
    /* n == 0: bulk == 0, so no kernel is entered and the tail loop does not run -- NULL is legal. */
    fv_vec2_scale_soa(NULL, NULL, 2.0f, NULL, NULL, 0u);
    fv_vec2_add_soa(NULL, NULL, NULL, NULL, NULL, NULL, 0u);
    fv_vec2_normalize_soa(NULL, NULL, NULL, NULL, 0u);
    fv_vec2_madd_normalize_soa(NULL, NULL, NULL, NULL, 2.0f, NULL, NULL, 0u);
    fv_vec2_madd_normalize_soa_stream(NULL, NULL, NULL, NULL, 2.0f, NULL, NULL, 0u);
    CHECK(1); /* reaching here means neither call dereferenced anything */
}

int main(void)
{
    run_tails();
    run_s_sweep();
    run_large();
    run_phase_shifts();
    run_path_independence();
    run_anything_goes();
    run_zero_length();

    printf("test_madd_norm: 1 - bit-level vs scale->add->normalize, n = 0..40 (every tail length), 5 strict corpora\n");
    printf("test_madd_norm: 2 - s sweep (1.0 / -1.0 / 0.5 / 3.25 / 1e30 / 1e-30 / +0.0 / -0.0)\n");
    printf("test_madd_norm: 3 - bit-level at n = 1023 / 4096 / 4097 / 8192 / 12345\n");
    printf("test_madd_norm: 4 - input/output phase shifts 0..7 and the _stream phase-mismatch fallback\n");
    printf("test_madd_norm: 5 - path independence (element bits do not depend on the bulk/tail split)\n");
    printf("test_madd_norm: 6 - anything-goes corpus: NaN-vs-NaN payload diffs tolerated = %d "
           "(out of %lu component comparisons)\n", g_tolerated, (unsigned long)g_compared);
    printf("test_madd_norm: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
