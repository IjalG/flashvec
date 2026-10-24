/* Bit-level acceptance fixture for the integrator primitives (P4 group):
 *   in-place batch: fv_vec2_addmul_x4 / fv_vec2_addmul_soa / fv_vec3_addmul_soa / fv_vec4_addmul_x2
 *   vector madd   : fv_vec2_maddv_x4 / fv_vec2_maddv_x4_stream
 *   per-value     : fv_vec2_addmul / fv_vec3_addmul / fv_vec4_addmul / fv_vec2_maddv
 *
 * Every comparison here is a bit comparison (memcmp of the float payload), not an epsilon test: the promise of this group
 * is "multiply first, then add, never contracted into an FMA", so the expected value is the hand-written scalar
 * expression, and one differing bit is a failure. The reference expressions below are compiled with the Makefile's
 * -ffp-contract=off (CXXFLAGS), i.e. the same rule the header's pragma enforces for the library side.
 *
 * What is covered:
 *   1..5   every entry, every tail length n = 0..40, three corpora (finite incl. -0.0 / subnormals / FLT_MAX, one NaN
 *          bit pattern, +-inf) x the s sweep; the reference is computed before an in-place call writes x
 *   6      large n and sizes that are not powers of two
 *   7      n == 0 with NULL for every argument (no pointer may be touched)
 *   8      in-place alias discipline: x and v are different arrays (adjacent and separated), v is left untouched
 *   9      misaligned in-place destination (4 B phase: the 32 B boundary is unreachable, which only matters for _stream)
 *   10     maddv _stream output phases 0/4/8/16/32/64 and the odd byte offsets 1/3/7 (no #GP, same bits, guards intact)
 *   11     path independence: an element's bits do not depend on whether it was handled by the AVX2 body or the scalar tail
 *   12     split invariance: one call over [0,n) equals two calls over [0,k) and [k,n)
 *   13     the per-value API equals the written-out scalar expression
 *   14     maddv writes nothing outside o[0,n) (poisoned guard elements stay poisoned)
 *
 * The NaN rule (the library-wide one, see include/flashvec.h): when two *different* NaN bit patterns meet inside one
 * arithmetic instruction, which payload comes out is decided by the instruction and by the compiler's freedom to swap
 * commutative operands. The strict corpora therefore carry at most one NaN pattern per comparison; the anything-goes
 * corpus (14) only requires "both sides NaN".
 *
 * Usage: test_addmul        (make test-addmul)
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum { CAP = 20000u, GUARD = 8u };

static fv_vec2 g_x2[CAP + GUARD], g_v2[CAP + GUARD], g_a2[CAP + GUARD], g_b2[CAP + GUARD];
static fv_vec2 g_c2[CAP + GUARD], g_o2[CAP + GUARD], g_r2[CAP + GUARD], g_s2[CAP + GUARD];
static fv_vec4 g_x4[CAP + GUARD], g_v4[CAP + GUARD], g_r4[CAP + GUARD];

static float g_sx[CAP + GUARD], g_sy[CAP + GUARD], g_sz[CAP + GUARD];
static float g_tx[CAP + GUARD], g_ty[CAP + GUARD], g_tz[CAP + GUARD];
static float g_ux[CAP + GUARD], g_uy[CAP + GUARD], g_uz[CAP + GUARD];

static unsigned long g_checks;
static unsigned long g_failures;
static uint32_t g_rng = 0x51ED270Bu;

enum { CORP_FINITE, CORP_NAN1, CORP_INF };
enum { GUARD_BYTE = 0x5A };

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

/* Finite corpus: random values with a sprinkle of the bit patterns that break naive kernels. */
static const uint32_t kFinite[] = {
    0x00000000u, 0x80000000u, 0x3f800000u, 0xbf800000u, 0x40600000u, 0xc0100000u,
    0x00800000u, 0x80800000u, 0x7f7fffffu, 0xff7fffffu, 0x00000001u, 0x80000001u,
    0x0da24260u, 0x8da24260u, 0x7149f2cau, 0xf149f2cau
};
#define NFINITE (sizeof kFinite / sizeof kFinite[0])

#define NAN_PATTERN 0x7fc00001u
#define INF_PATTERN 0x7f800000u

static float gen(int corpus)
{
    switch (corpus) {
    case CORP_NAN1:
        return from_bits((rnd() & 1u) != 0u ? NAN_PATTERN : kFinite[rnd() % NFINITE]);
    case CORP_INF:
        return from_bits((rnd() & 3u) == 0u ? ((rnd() & 1u) != 0u ? INF_PATTERN : (INF_PATTERN | 0x80000000u))
                                             : kFinite[rnd() % NFINITE]);
    case CORP_FINITE:
    default:
        return (rnd() & 7u) == 0u ? from_bits(kFinite[rnd() % NFINITE]) : rnd_unit();
    }
}

static const float kSweep[] = { 1.0f, -1.0f, 0.5f, 3.25f, 1.0e30f, 1.0e-30f, 0.0f, -0.0f };
#define NSWEEP (sizeof kSweep / sizeof kSweep[0])

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

static void poison_vec2(fv_vec2 *dst, size_t n)
{
    memset(dst, GUARD_BYTE, n * sizeof *dst);
}

/* ---- comparison helpers: one element = one check ---- */
static void expect_vec2(const fv_vec2 *got, const fv_vec2 *want, size_t n, const char *what, size_t off)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        ++g_checks;
        if (bits(got[i].x) != bits(want[i].x) || bits(got[i].y) != bits(want[i].y)) {
            if (g_failures < 8) {
                printf("  FAIL %s [%zu] got (%08x %08x) want (%08x %08x)\n", what, off + i,
                       bits(got[i].x), bits(got[i].y), bits(want[i].x), bits(want[i].y));
            }
            ++g_failures;
        }
    }
}

static void expect_vec4(const fv_vec4 *got, const fv_vec4 *want, size_t n, const char *what, size_t off)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        ++g_checks;
        if (bits(got[i].x) != bits(want[i].x) || bits(got[i].y) != bits(want[i].y) ||
            bits(got[i].z) != bits(want[i].z) || bits(got[i].w) != bits(want[i].w)) {
            if (g_failures < 8) {
                printf("  FAIL %s [%zu]\n", what, off + i);
            }
            ++g_failures;
        }
    }
}

static void expect_f(const float *got, const float *want, size_t n, const char *what, size_t off)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        ++g_checks;
        if (bits(got[i]) != bits(want[i])) {
            if (g_failures < 8) {
                printf("  FAIL %s [%zu] got %08x want %08x\n", what, off + i, bits(got[i]), bits(want[i]));
            }
            ++g_failures;
        }
    }
}

/* ---- references (the expressions the promise names) ---- */
static void ref_addmul_vec2(fv_vec2 *r, const fv_vec2 *x, const fv_vec2 *v, fv_real s, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        r[i].x = x[i].x + v[i].x * s;
        r[i].y = x[i].y + v[i].y * s;
    }
}

static void ref_addmul_vec4(fv_vec4 *r, const fv_vec4 *x, const fv_vec4 *v, fv_real s, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        r[i].x = x[i].x + v[i].x * s;
        r[i].y = x[i].y + v[i].y * s;
        r[i].z = x[i].z + v[i].z * s;
        r[i].w = x[i].w + v[i].w * s;
    }
}

static void ref_maddv(fv_vec2 *r, const fv_vec2 *a, const fv_vec2 *b, const fv_vec2 *c, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        r[i].x = a[i].x * b[i].x + c[i].x;
        r[i].y = a[i].y * b[i].y + c[i].y;
    }
}

/* One pass of every entry over a fixed n, corpus and scale. Returns nothing; adds n checks per entry. */
static void run_all_entries(size_t n, int corpus, float s, const char *tag)
{
    char what[96];

    /* ---- in-place AoS vec2 ---- */
    fill_vec2(g_x2, n, corpus);
    fill_vec2(g_v2, n, corpus);
    ref_addmul_vec2(g_r2, g_x2, g_v2, s, n);
    poison_vec2(g_x2 + n, GUARD);
    fv_vec2_addmul_x4(g_x2, g_v2, s, n);
    snprintf(what, sizeof what, "%s addmul_x4", tag);
    expect_vec2(g_x2, g_r2, n, what, 0);

    /* ---- in-place SoA vec2 ---- */
    fill_f(g_sx, n, corpus);
    fill_f(g_sy, n, corpus);
    fill_f(g_tx, n, corpus);
    fill_f(g_ty, n, corpus);
    {
        size_t i;
        for (i = 0; i < n; ++i) {
            g_ux[i] = g_sx[i] + g_tx[i] * s;
            g_uy[i] = g_sy[i] + g_ty[i] * s;
        }
    }
    fv_vec2_addmul_soa(g_sx, g_sy, g_tx, g_ty, s, n);
    snprintf(what, sizeof what, "%s addmul_soa(vec2)", tag);
    expect_f(g_sx, g_ux, n, what, 0);
    expect_f(g_sy, g_uy, n, what, 0);

    /* ---- in-place SoA vec3 ---- */
    fill_f(g_sz, n, corpus);
    fill_f(g_tz, n, corpus);
    {
        size_t i;
        for (i = 0; i < n; ++i) {
            g_ux[i] = g_sx[i] + g_tx[i] * s;
            g_uy[i] = g_sy[i] + g_ty[i] * s;
            g_uz[i] = g_sz[i] + g_tz[i] * s;
        }
    }
    fv_vec3_addmul_soa(g_sx, g_sy, g_sz, g_tx, g_ty, g_tz, s, n);
    snprintf(what, sizeof what, "%s addmul_v3_soa", tag);
    expect_f(g_sx, g_ux, n, what, 0);
    expect_f(g_sy, g_uy, n, what, 0);
    expect_f(g_sz, g_uz, n, what, 0);

    /* ---- in-place AoS vec4 ---- */
    fill_vec4(g_x4, n, corpus);
    fill_vec4(g_v4, n, corpus);
    ref_addmul_vec4(g_r4, g_x4, g_v4, s, n);
    fv_vec4_addmul_x2(g_x4, g_v4, s, n);
    snprintf(what, sizeof what, "%s addmul_x2", tag);
    expect_vec4(g_x4, g_r4, n, what, 0);

    /* ---- maddv, normal and streaming ---- */
    fill_vec2(g_a2, n, corpus);
    fill_vec2(g_b2, n, corpus);
    fill_vec2(g_c2, n, corpus);
    ref_maddv(g_r2, g_a2, g_b2, g_c2, n);
    poison_vec2(g_o2, n + GUARD);
    fv_vec2_maddv_x4(g_a2, g_b2, g_c2, g_o2, n);
    snprintf(what, sizeof what, "%s maddv_x4", tag);
    expect_vec2(g_o2, g_r2, n, what, 0);

    fill_vec2(g_s2, n, corpus); /* scratch: the stream form must agree with the plain one bit for bit */
    memcpy(g_s2, g_o2, n * sizeof *g_s2);
    poison_vec2(g_o2, n + GUARD);
    fv_vec2_maddv_x4_stream(g_a2, g_b2, g_c2, g_o2, n);
    snprintf(what, sizeof what, "%s maddv_x4_stream", tag);
    expect_vec2(g_o2, g_r2, n, what, 0);
    snprintf(what, sizeof what, "%s maddv_x4_stream vs non-stream", tag);
    expect_vec2(g_o2, g_s2, n, what, 0);
}

int main(void)
{
    size_t n;
    size_t si;
    int corpus;
    static const char *kCorpusName[] = { "finite", "nan1", "inf" };
    char tag[64];

    printf("test_addmul: 1 - addmul_x4 / addmul_soa / addmul_v3_soa / addmul_x2 in place, n = 0..40, 3 corpora x 8 scales\n");
    for (n = 0; n <= 40; ++n) {
        for (corpus = CORP_FINITE; corpus <= CORP_INF; ++corpus) {
            for (si = 0; si < NSWEEP; ++si) {
                snprintf(tag, sizeof tag, "n=%zu %s s=%g", n, kCorpusName[corpus], (double)kSweep[si]);
                run_all_entries(n, corpus, kSweep[si], tag);
            }
        }
    }

    printf("test_addmul: 2 - maddv_x4 / _stream, n = 0..40, 3 corpora\n");
    for (n = 0; n <= 40; ++n) {
        for (corpus = CORP_FINITE; corpus <= CORP_INF; ++corpus) {
            fill_vec2(g_a2, n, corpus);
            fill_vec2(g_b2, n, corpus);
            fill_vec2(g_c2, n, corpus);
            ref_maddv(g_r2, g_a2, g_b2, g_c2, n);
            poison_vec2(g_o2, n + GUARD);
            fv_vec2_maddv_x4(g_a2, g_b2, g_c2, g_o2, n);
            snprintf(tag, sizeof tag, "n=%zu %s maddv", n, kCorpusName[corpus]);
            expect_vec2(g_o2, g_r2, n, tag, 0);
            snprintf(tag, sizeof tag, "n=%zu %s maddv guard", n, kCorpusName[corpus]);
            {
                size_t i;
                for (i = n; i < n + GUARD; ++i) {
                    ++g_checks;
                    if (bits(g_o2[i].x) != 0x5a5a5a5au || bits(g_o2[i].y) != 0x5a5a5a5au) {
                        if (g_failures < 8) {
                            printf("  FAIL %s [%zu] written outside o[0,n)\n", tag, i);
                        }
                        ++g_failures;
                    }
                }
            }
            poison_vec2(g_o2, n + GUARD);
            fv_vec2_maddv_x4_stream(g_a2, g_b2, g_c2, g_o2, n);
            snprintf(tag, sizeof tag, "n=%zu %s maddv_stream", n, kCorpusName[corpus]);
            expect_vec2(g_o2, g_r2, n, tag, 0);
            {
                size_t i;
                for (i = n; i < n + GUARD; ++i) {
                    ++g_checks;
                    if (bits(g_o2[i].x) != 0x5a5a5a5au || bits(g_o2[i].y) != 0x5a5a5a5au) {
                        if (g_failures < 8) {
                            printf("  FAIL %s [%zu] written outside o[0,n) (stream)\n", tag, i);
                        }
                        ++g_failures;
                    }
                }
            }
        }
    }

    printf("test_addmul: 3 - large n and sizes that are not powers of two\n");
    {
        static const size_t kBig[] = { 63, 100, 255, 1000, 1023, 1024, 4095, 4096, 4097, 8192, 12345, 16001 };
        size_t k;
        for (k = 0; k < sizeof kBig / sizeof kBig[0]; ++k) {
            snprintf(tag, sizeof tag, "big n=%zu", kBig[k]);
            run_all_entries(kBig[k], CORP_FINITE, 0.0625f, tag);
            run_all_entries(kBig[k], CORP_FINITE, -3.0e-3f, tag);
        }
        /* an odd n that is not in the table and not near a vector multiple */
        run_all_entries(7999, CORP_FINITE, 1.0f, "odd n=7999");
    }

    printf("test_addmul: 4 - n == 0 with NULL pointers (no pointer may be touched)\n");
    {
        /* distinct null variables: passing one NULL expression to two __restrict parameters is what -Wrestrict reports,
         * and the point here is the n == 0 contract, not aliasing */
        fv_vec2 *nx2 = (fv_vec2 *)0;
        const fv_vec2 *na2 = (const fv_vec2 *)0;
        const fv_vec2 *nb2 = (const fv_vec2 *)0;
        const fv_vec2 *nc2 = (const fv_vec2 *)0;
        fv_vec4 *nx4 = (fv_vec4 *)0;
        const fv_vec4 *nv4 = (const fv_vec4 *)0;
        float *fx = (float *)0;
        float *fy = (float *)0;
        float *fz = (float *)0;
        const float *cvx = (const float *)0;
        const float *cvy = (const float *)0;
        const float *cvz = (const float *)0;
        fv_vec2_maddv_x4(na2, nb2, nc2, nx2, 0);
        fv_vec2_maddv_x4_stream(na2, nb2, nc2, nx2, 0);
        fv_vec2_addmul_x4(nx2, na2, 1.0f, 0);
        fv_vec4_addmul_x2(nx4, nv4, 1.0f, 0);
        fv_vec2_addmul_soa(fx, fy, cvx, cvy, 1.0f, 0);
        fv_vec3_addmul_soa(fx, fy, fz, cvx, cvy, cvz, 1.0f, 0);
        ++g_checks;
        printf("  6 entries accepted n == 0 with NULL (no crash)\n");
    }

    printf("test_addmul: 5 - in-place alias discipline: x and v are different arrays, v is untouched\n");
    {
        /* Adjacent: v starts exactly one element past the end of x. Separated: a gap of 4 elements. */
        static fv_vec2 bufx[64 + 64 + 8];
        static fv_vec2 vcopy[64];
        fv_vec2 ref[64], xsave[64];
        size_t i;
        int mode;
        for (mode = 0; mode < 2; ++mode) {
            size_t gap = (mode == 0) ? 0u : 4u;
            fv_vec2 *x = bufx;
            fv_vec2 *v = bufx + 64 + gap;
            fill_vec2(x, 64, CORP_FINITE);
            fill_vec2(v, 64, CORP_FINITE);
            memcpy(xsave, x, sizeof xsave);
            memcpy(vcopy, v, sizeof(fv_vec2) * 64);
            ref_addmul_vec2(ref, xsave, vcopy, 0.75f, 64);
            fv_vec2_addmul_x4(x, v, 0.75f, 64);
            expect_vec2(x, ref, 64, mode == 0 ? "adjacent x//v addmul" : "separated x..v addmul", 0);
            for (i = 0; i < 64; ++i) {
                ++g_checks;
                if (bits(v[i].x) != bits(vcopy[i].x) || bits(v[i].y) != bits(vcopy[i].y)) {
                    if (g_failures < 8) {
                        printf("  FAIL the in-place form modified v[%zu]\n", i);
                    }
                    ++g_failures;
                }
            }
        }
    }

    printf("test_addmul: 6 - misaligned in-place destination (4 B phase) and SoA output phases 1..7\n");
    {
        size_t ph;
        for (ph = 1; ph <= 7; ++ph) {
            fv_vec2 *x = g_x2 + ph;
            fv_vec2 *v = g_v2 + ph;
            size_t nn = 257;
            fill_vec2(x, nn, CORP_FINITE);
            fill_vec2(v, nn, CORP_FINITE);
            memcpy(g_s2, x, nn * sizeof *x);
            {
                size_t i;
                for (i = 0; i < nn; ++i) {
                    g_r2[i].x = g_s2[i].x + v[i].x * 1.5f;
                    g_r2[i].y = g_s2[i].y + v[i].y * 1.5f;
                }
            }
            fv_vec2_addmul_x4(x, v, 1.5f, nn);
            expect_vec2(x, g_r2, nn, "misaligned addmul_x4", 0);

            fill_f(g_sx + ph, nn, CORP_FINITE);
            fill_f(g_tx + ph, nn, CORP_FINITE);
            fill_f(g_sy + ph, nn, CORP_FINITE);
            fill_f(g_ty + ph, nn, CORP_FINITE);
            {
                size_t i;
                for (i = 0; i < nn; ++i) {
                    g_ux[i] = g_sx[ph + i] + g_tx[ph + i] * 0.5f;
                }
            }
            fv_vec2_addmul_soa(g_sx + ph, g_sy + ph, g_tx + ph, g_ty + ph, 0.5f, nn);
            expect_f(g_sx + ph, g_ux, nn, "misaligned addmul_soa", 0);
        }
    }

    printf("test_addmul: 7 - maddv _stream output pointer phases 0/4/8/16/32/64 and odd 1/3/7\n");
    {
        static unsigned char raw[4096 + 128];
        static const size_t kPhase[] = { 0, 4, 8, 16, 32, 64, 1, 3, 7 };
        size_t k;
        size_t nn = 251; /* odd: forces body + scalar tail after the peel */
        fill_vec2(g_a2, nn, CORP_FINITE);
        fill_vec2(g_b2, nn, CORP_FINITE);
        fill_vec2(g_c2, nn, CORP_FINITE);
        ref_maddv(g_r2, g_a2, g_b2, g_c2, nn);
        for (k = 0; k < sizeof kPhase / sizeof kPhase[0]; ++k) {
            size_t off = kPhase[k];
            fv_vec2 *o;
            unsigned char before[128];
            memset(raw, GUARD_BYTE, sizeof raw);
            o = (fv_vec2 *)(void *)(raw + off);
            memcpy(before, raw, sizeof raw < sizeof before ? sizeof raw : sizeof before);
            fv_vec2_maddv_x4_stream(g_a2, g_b2, g_c2, o, nn);
            {
                const unsigned char *got = raw + off;
                size_t i;
                for (i = 0; i < nn; ++i) {
                    ++g_checks;
                    if (memcmp(got + i * sizeof(fv_vec2), &g_r2[i], sizeof(fv_vec2)) != 0) {
                        if (g_failures < 8) {
                            printf("  FAIL maddv_stream phase %zu element %zu differs\n", off, i);
                        }
                        ++g_failures;
                    }
                }
                /* the bytes before the output and after its end stay poisoned */
                ++g_checks;
                if (off != 0 && memcmp(raw, before, off) != 0) {
                    printf("  FAIL maddv_stream phase %zu wrote before o\n", off);
                    ++g_failures;
                }
                for (i = off + nn * sizeof(fv_vec2); i + 8 <= sizeof raw; i += 8) {
                    ++g_checks;
                    if (memcmp(raw + i, "\x5a\x5a\x5a\x5a\x5a\x5a\x5a\x5a", 8) != 0) {
                        if (g_failures < 8) {
                            printf("  FAIL maddv_stream phase %zu wrote past o[0,n)\n", off);
                        }
                        ++g_failures;
                        break;
                    }
                }
            }
            /* and the non-stream form must produce the same bytes at the same phase */
            memset(raw, GUARD_BYTE, sizeof raw);
            o = (fv_vec2 *)(void *)(raw + off);
            fv_vec2_maddv_x4(g_a2, g_b2, g_c2, o, nn);
            ++g_checks;
            if (memcmp(raw + off, g_r2, nn * sizeof(fv_vec2)) != 0) {
                printf("  FAIL maddv (non-stream) phase %zu differs\n", off);
                ++g_failures;
            }
        }
    }

    printf("test_addmul: 8 - path independence: AVX2 lane vs scalar tail give the same bits\n");
    {
        size_t nn = 1003; /* 250 x 4 + 3: three elements fall to the scalar tail */
        fill_vec2(g_x2, nn, CORP_FINITE);
        fill_vec2(g_v2, nn, CORP_FINITE);
        memcpy(g_s2, g_x2, nn * sizeof *g_x2); /* the pre-call values, for the reference */
        fv_vec2_addmul_x4(g_x2, g_v2, 2.5f, nn);
        {
            size_t i;
            for (i = 0; i < nn; ++i) {
                fv_vec2 one_x, one_v, ref;
                one_x = g_s2[i];
                one_v = g_v2[i];
                ref.x = one_x.x + one_v.x * 2.5f;
                ref.y = one_x.y + one_v.y * 2.5f;
                ++g_checks;
                if (bits(g_x2[i].x) != bits(ref.x) || bits(g_x2[i].y) != bits(ref.y)) {
                    if (g_failures < 8) {
                        printf("  FAIL element %zu differs between the bulk path and the scalar path\n", i);
                    }
                    ++g_failures;
                }
            }
        }
    }

    printf("test_addmul: 9 - split invariance (one call vs two calls)\n");
    {
        size_t nn = 1000, k;
        for (k = 1; k < nn; k += 137) {
            fv_vec2 whole[1000], split[1000];
            fill_vec2(whole, nn, CORP_FINITE);
            fill_vec2(g_v2, nn, CORP_FINITE);
            memcpy(split, whole, sizeof whole);
            fv_vec2_addmul_x4(whole, g_v2, 0.5f, nn);
            fv_vec2_addmul_x4(split, g_v2, 0.5f, k);
            fv_vec2_addmul_x4(split + k, g_v2 + k, 0.5f, nn - k);
            expect_vec2(split, whole, nn, "split invariance addmul", 0);
        }
    }

    printf("test_addmul: 10 - per-value API equals the written-out expression\n");
    {
        size_t i;
        for (i = 0; i < 200; ++i) {
            fv_vec2 x2, v2, c2, r2, e2;
            fv_vec3 x3, v3, r3, e3;
            fv_vec4 x4, v4, r4, e4;
            float s = gen(CORP_FINITE);
            x2.x = gen(CORP_FINITE); x2.y = gen(CORP_FINITE);
            v2.x = gen(CORP_FINITE); v2.y = gen(CORP_FINITE);
            c2.x = gen(CORP_FINITE); c2.y = gen(CORP_FINITE);
            r2 = fv_vec2_addmul(x2, v2, s);
            e2.x = x2.x + v2.x * s;
            e2.y = x2.y + v2.y * s;
            ++g_checks;
            if (bits(r2.x) != bits(e2.x) || bits(r2.y) != bits(e2.y)) {
                printf("  FAIL fv_vec2_addmul\n");
                ++g_failures;
            }
            r2 = fv_vec2_maddv(x2, v2, c2);
            e2.x = x2.x * v2.x + c2.x;
            e2.y = x2.y * v2.y + c2.y;
            ++g_checks;
            if (bits(r2.x) != bits(e2.x) || bits(r2.y) != bits(e2.y)) {
                printf("  FAIL fv_vec2_maddv\n");
                ++g_failures;
            }
            x3.x = gen(CORP_FINITE); x3.y = gen(CORP_FINITE); x3.z = gen(CORP_FINITE);
            v3.x = gen(CORP_FINITE); v3.y = gen(CORP_FINITE); v3.z = gen(CORP_FINITE);
            r3 = fv_vec3_addmul(x3, v3, s);
            e3.x = x3.x + v3.x * s;
            e3.y = x3.y + v3.y * s;
            e3.z = x3.z + v3.z * s;
            ++g_checks;
            if (bits(r3.x) != bits(e3.x) || bits(r3.y) != bits(e3.y) ||
                bits(r3.z) != bits(e3.z)) {
                printf("  FAIL fv_vec3_addmul\n");
                ++g_failures;
            }
            x4.x = gen(CORP_FINITE); x4.y = gen(CORP_FINITE);
            x4.z = gen(CORP_FINITE); x4.w = gen(CORP_FINITE);
            v4.x = gen(CORP_FINITE); v4.y = gen(CORP_FINITE);
            v4.z = gen(CORP_FINITE); v4.w = gen(CORP_FINITE);
            r4 = fv_vec4_addmul(x4, v4, s);
            e4.x = x4.x + v4.x * s;
            e4.y = x4.y + v4.y * s;
            e4.z = x4.z + v4.z * s;
            e4.w = x4.w + v4.w * s;
            ++g_checks;
            if (bits(r4.x) != bits(e4.x) || bits(r4.y) != bits(e4.y) ||
                bits(r4.z) != bits(e4.z) || bits(r4.w) != bits(e4.w)) {
                printf("  FAIL fv_vec4_addmul\n");
                ++g_failures;
            }
        }
        /* the per-value form and the batch form of the same expression must agree bit for bit */
        fill_vec2(g_x2, 64, CORP_FINITE);
        fill_vec2(g_v2, 64, CORP_FINITE);
        ref_addmul_vec2(g_r2, g_x2, g_v2, 3.25f, 64);
        fv_vec2_addmul_x4(g_x2, g_v2, 3.25f, 64);
        expect_vec2(g_x2, g_r2, 64, "per-value reference vs addmul_x4", 0);
        {
            size_t i;
            for (i = 0; i < 64; ++i) {
                fv_vec2 one = fv_vec2_addmul(g_r2[i], g_v2[i], 3.25f);
                ++g_checks;
                if (bits(one.x) != bits(g_r2[i].x + g_v2[i].x * 3.25f) ||
                    bits(one.y) != bits(g_r2[i].y + g_v2[i].y * 3.25f)) {
                    printf("  FAIL per-value vs reference at %zu\n", i);
                    ++g_failures;
                }
            }
        }
    }

    printf("test_addmul: 11 - anything-goes corpus: both sides NaN counts as equal (payload diffs tolerated)\n");
    {
        size_t i, nn = 4096, tolerated = 0, bad = 0;
        fill_vec2(g_a2, nn, CORP_NAN1);
        fill_vec2(g_b2, nn, CORP_NAN1);
        fill_vec2(g_c2, nn, CORP_NAN1);
        ref_maddv(g_r2, g_a2, g_b2, g_c2, nn);
        fv_vec2_maddv_x4(g_a2, g_b2, g_c2, g_o2, nn);
        for (i = 0; i < nn; ++i) {
            int k;
            for (k = 0; k < 2; ++k) {
                float g = k == 0 ? g_o2[i].x : g_o2[i].y;
                float w = k == 0 ? g_r2[i].x : g_r2[i].y;
                ++g_checks;
                if (bits(g) == bits(w)) {
                    continue;
                }
                if (isnan(g) && isnan(w)) {
                    ++tolerated; /* different payload of a NaN: allowed by the library-wide rule */
                } else {
                    if (bad < 8) {
                        printf("  FAIL NaN corpus element %zu component %d got %08x want %08x\n", i, k,
                               bits(g), bits(w));
                    }
                    ++g_failures;
                    ++bad;
                }
            }
        }
        printf("test_addmul: NaN-payload differences tolerated = %zu (both sides NaN)\n", tolerated);
    }

    printf("test_addmul: %lu checks, %lu failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
