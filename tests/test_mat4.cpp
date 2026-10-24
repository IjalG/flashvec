/* Bit-level test for fv_mat4_transform_points / fv_mat4_transform_dirs.
 *
 * "Bit-level" means the whole output buffer is compared with memcmp against an independent
 * scalar reference that uses no intrinsics. Not a tolerance check -- a 1-ulp difference fails.
 *
 * The reference evaluation order is ((m0*x + m1*y) + m2*z) + m3 (dirs: no fourth term), which is
 * what the AVX2 kernel does lane by lane. That, plus no FMA on either side, is why they agree
 * bit-for-bit. The reference is wrapped in the same fp-contract=off guard as the header's
 * implementation section: without it, -ffp-contract=fast (-mfma) would fuse the reference's own
 * mul+add pairs into vfmadd and the comparison would measure compiler flags, not the library.
 *
 * The one documented boundary (section 3, and the header contract): when one addition gets
 * **two NaN addends with different bit patterns**, the resulting NaN payload/sign depends on the
 * instruction operand order, which C does not pin down and intrinsics do. Both paths then produce
 * a legal NaN, possibly different ones. Everywhere else -- including a NaN operand meeting an
 * infinity, or infinities only, or a single NaN pattern anywhere -- the comparison is strict.
 * Sections 1, 2 and 6 therefore run corpora that contain at most one NaN bit pattern, and assert
 * zero payload differences; section 3 runs the anything-goes corpus and only requires "NaN vs NaN"
 * while counting and reporting the payload differences instead of hiding them.
 *
 * Coverage:
 *   - n = 0..40: every tail length of the 2-elements-per-ymm loop, three strict corpora
 *   - n = 1024 / 4097 / 8192: the bulk loop under a real length (4097 is deliberately odd)
 *   - values: -0.0, +0.0, FLT_MAX/FLT_MIN, smallest denormal, +-inf, NaN, random
 *   - matrices that themselves contain inf / NaN / denormals
 *   - the input w component is deliberate garbage: both entries must ignore it
 *   - section 6 compares the AVX2 path against the scalar tail for the *same element*, i.e. the
 *     result must not depend on where the element sits relative to the bulk/tail boundary
 *   - section 6 exercises that boundary on purpose: it asserts "NaN on both paths" and counts the
 *     payload difference; the tolerated count is printed at the end, never swallowed
 *   - semantics: identity, translation (points move, dirs do not), scaling, w output, the last
 *     matrix row is never read, n == 0 touches no pointer
 *
 * Run: make test-mat4
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0;
static int g_failures = 0;
static int g_nan_payload = 0; /* section 3 only: counted and printed, never hidden */

#define CHECK(cond)                                                    \
    do {                                                               \
        ++g_checks;                                                    \
        if (!(cond)) {                                                 \
            ++g_failures;                                              \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);     \
        }                                                              \
    } while (0)

#define CAP   8200u
#define GUARD 4u

static fv_vec4 g_in[CAP + GUARD];
static fv_vec4 g_got[CAP + GUARD];
static fv_vec4 g_ref[CAP + GUARD];

/* ------------------------------------------------------------------ helpers */

static uint32_t bits_u32(fv_real f)
{
    uint32_t b;
    memcpy(&b, &f, sizeof(b));
    return b;
}

static int bits_eq(fv_real a, fv_real b)
{
    return bits_u32(a) == bits_u32(b);
}

static int is_nan_bits(uint32_t b)
{
    return (b & 0x7F800000u) == 0x7F800000u && (b & 0x007FFFFFu) != 0u;
}

static int bytes_eq(const void *p, const void *q, size_t n)
{
    return memcmp(p, q, n) == 0;
}

/* 0x5A poison: "wrote what it should not have" shows up as a guard byte that moved. */
static void poison_all(void)
{
    memset(g_got, 0x5A, sizeof(g_got));
}

static int guard_intact(size_t n)
{
    const unsigned char *p = (const unsigned char *)(const void *)(g_got + n);
    size_t i;
    for (i = 0; i < GUARD * sizeof(fv_vec4); ++i) {
        if (p[i] != 0x5A) {
            return 0;
        }
    }
    return 1;
}

/* Elements that differ only because both sides are NaN with different bits. Returns the number
 * of components in that class; any other difference is reported and left to the caller. */
static size_t nan_only_diff(const fv_vec4 *got, const fv_vec4 *want, size_t n)
{
    size_t i;
    size_t diffs = 0;
    for (i = 0; i < n; ++i) {
        const fv_real *a = &got[i].x;
        const fv_real *b = &want[i].x;
        int c;
        for (c = 0; c < 4; ++c) {
            const uint32_t ba = bits_u32(a[c]);
            const uint32_t bb = bits_u32(b[c]);
            if (ba == bb) {
                continue;
            }
            if (is_nan_bits(ba) && is_nan_bits(bb)) {
                ++diffs;
                continue;
            }
            printf("  diff at element %u component %d: got 0x%08lx want 0x%08lx\n",
                   (unsigned)i, c, (unsigned long)ba, (unsigned long)bb);
            return (size_t)-1;
        }
    }
    return diffs;
}

static unsigned int g_rng = 0x9E3779B9u;

static unsigned int rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/* Random finite value: exponent 1..254, so it is never NaN and never inf. */
static fv_real rnd_finite(void)
{
    const unsigned int b = rnd_u32();
    const unsigned int bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    fv_real f;
    memcpy(&f, &bits, sizeof(f));
    return f;
}

static fv_real rnd_moderate(void)
{
    return (fv_real)((int)(rnd_u32() % 2001u) - 1000) / 16.0f;
}

/* ------------------------------------------------------------------ reference */

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif

static void ref_points(const fv_mat4 *m, const fv_vec4 *in, fv_vec4 *out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        const fv_real x = in[i].x;
        const fv_real y = in[i].y;
        const fv_real z = in[i].z;
        fv_vec4 r;
        r.x = ((m->m[0] * x + m->m[1] * y) + m->m[2] * z) + m->m[3];
        r.y = ((m->m[4] * x + m->m[5] * y) + m->m[6] * z) + m->m[7];
        r.z = ((m->m[8] * x + m->m[9] * y) + m->m[10] * z) + m->m[11];
        r.w = (fv_real)1;
        out[i] = r;
    }
}

static void ref_dirs(const fv_mat4 *m, const fv_vec4 *in, fv_vec4 *out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        const fv_real x = in[i].x;
        const fv_real y = in[i].y;
        const fv_real z = in[i].z;
        fv_vec4 r;
        r.x = (m->m[0] * x + m->m[1] * y) + m->m[2] * z;
        r.y = (m->m[4] * x + m->m[5] * y) + m->m[6] * z;
        r.z = (m->m[8] * x + m->m[9] * y) + m->m[10] * z;
        r.w = (fv_real)0;
        out[i] = r;
    }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

/* ------------------------------------------------------------------ corpus
 * The three strict corpora each contain at most ONE NaN bit pattern (or none at all), which is
 * exactly the condition under which the scalar tail and the AVX2 kernel must agree bit-for-bit.
 * The anything-goes corpus is used only where NaN payload differences are tolerated and counted. */

static const fv_real kExtreme[] = {
    0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f, 123456.789f,
    FLT_MAX, -FLT_MAX, FLT_MIN, -FLT_MIN,
    0x1p-149f, -0x1p-149f,                        /* smallest positive / negative denormal */
    (fv_real)INFINITY, (fv_real)-INFINITY
};
#define NEXTREME ((unsigned)(sizeof(kExtreme) / sizeof(kExtreme[0])))

static const fv_real kNanMix[] = {
    (fv_real)NAN, 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, 1.0e30f, -1.0e30f
};
#define NNANMIX ((unsigned)(sizeof(kNanMix) / sizeof(kNanMix[0])))

enum { CORP_FINITE = 0, CORP_EXTREME = 1, CORP_NAN = 2, CORP_MIXED = 3 };

static void fill_corpus(fv_vec4 *v, size_t n, int corpus, unsigned int salt)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        /* the w input is ignored by both entries; fill it with junk on purpose */
        v[i].w = ((i & 1u) != 0u) ? (fv_real)-7.5 : (fv_real)3.25;
        switch (corpus) {
        case CORP_FINITE:
            v[i].x = rnd_moderate();
            v[i].y = rnd_moderate();
            v[i].z = rnd_moderate();
            break;
        case CORP_EXTREME:
            v[i].x = kExtreme[(i + salt) % NEXTREME];
            v[i].y = kExtreme[(i * 3u + salt + 1u) % NEXTREME];
            v[i].z = kExtreme[(i * 7u + salt + 2u) % NEXTREME];
            break;
        case CORP_NAN:
            v[i].x = kNanMix[(i + salt) % NNANMIX];
            v[i].y = kNanMix[(i * 3u + salt + 1u) % NNANMIX];
            v[i].z = kNanMix[(i * 7u + salt + 2u) % NNANMIX];
            break;
        default:
            v[i].x = kExtreme[(i + salt) % NEXTREME];
            v[i].y = kNanMix[(i * 5u + salt + 3u) % NNANMIX];
            v[i].z = rnd_finite();
            break;
        }
    }
}

static void make_matrix(fv_mat4 *m)
{
    int r;
    for (r = 0; r < 3; ++r) {
        int c;
        for (c = 0; c < 4; ++c) {
            m->m[r * 4 + c] = rnd_moderate();
        }
    }
    m->m[12] = (fv_real)0;
    m->m[13] = (fv_real)0;
    m->m[14] = (fv_real)0;
    m->m[15] = (fv_real)1;
}

/* Coefficients are infinite / denormal / signed zero, but never NaN: an infinite coefficient
 * multiplies a NaN input only in the anything-goes corpus, where payloads are not compared. */
static void make_wild_matrix(fv_mat4 *m)
{
    make_matrix(m);
    m->m[1] = (fv_real)INFINITY;
    m->m[6] = -0.0f;
    m->m[10] = 0x1p-149f;
    m->m[3] = -0.0f;
    m->m[7] = FLT_MAX;
}

/* ------------------------------------------------------------------ driver */

/* strict == 1: any difference fails. strict == 0: differing NaN payloads are counted. */
static void check_case(const fv_mat4 *m, size_t n, int strict, const char *label)
{
    size_t d;

    poison_all();
    fv_mat4_transform_points(m, g_in, g_got, n);
    ref_points(m, g_in, g_ref, n);
    d = bytes_eq(g_got, g_ref, n * sizeof(fv_vec4)) ? 0u : nan_only_diff(g_got, g_ref, n);
    if (d == (size_t)-1) {
        ++g_checks;
        ++g_failures;
        printf("FAIL points  %s n=%u: value differs\n", label, (unsigned)n);
    } else if (strict == 0) {
        ++g_checks;
        g_nan_payload += (int)d;
    } else {
        ++g_checks;
        if (d != 0u) {
            ++g_failures;
            printf("FAIL points  %s n=%u: %u NaN payload diffs (corpus must allow only one NaN pattern)\n",
                   label, (unsigned)n, (unsigned)d);
        }
    }
    CHECK(guard_intact(n));

    poison_all();
    fv_mat4_transform_dirs(m, g_in, g_got, n);
    ref_dirs(m, g_in, g_ref, n);
    d = bytes_eq(g_got, g_ref, n * sizeof(fv_vec4)) ? 0u : nan_only_diff(g_got, g_ref, n);
    if (d == (size_t)-1) {
        ++g_checks;
        ++g_failures;
        printf("FAIL dirs    %s n=%u: value differs\n", label, (unsigned)n);
    } else if (strict == 0) {
        ++g_checks;
        g_nan_payload += (int)d;
    } else {
        ++g_checks;
        if (d != 0u) {
            ++g_failures;
            printf("FAIL dirs    %s n=%u: %u NaN payload diffs (corpus must allow only one NaN pattern)\n",
                   label, (unsigned)n, (unsigned)d);
        }
    }
    CHECK(guard_intact(n));
}

static void run_tails(void)
{
    static const char *kName[] = { "finite", "extreme", "nan" };
    static const int kCorpus[] = { CORP_FINITE, CORP_EXTREME, CORP_NAN };
    fv_mat4 m;
    size_t n;
    int c;

    printf("test_mat4: 1 - bit-level vs scalar reference, n = 0..40 (every tail length)\n");
    for (c = 0; c < 3; ++c) {
        make_matrix(&m);
        for (n = 0; n <= 40; ++n) {
            fill_corpus(g_in, n + GUARD, kCorpus[c], (unsigned)n);
            check_case(&m, n, 1, kName[c]);
        }
    }

    printf("test_mat4: 2 - same, matrix with inf / denormal / signed-zero coefficients\n");
    {
        fv_mat4 wild;
        make_wild_matrix(&wild);
        for (n = 0; n <= 40; ++n) {
            fill_corpus(g_in, n + GUARD, CORP_EXTREME, (unsigned)(n + 3u));
            check_case(&wild, n, 1, "wild-matrix");
        }
    }

    printf("test_mat4: 3 - anything-goes corpus: NaN payload differences are counted, not hidden\n");
    {
        fv_mat4 wild;
        make_wild_matrix(&wild);
        for (n = 0; n <= 40; ++n) {
            fill_corpus(g_in, n + GUARD, CORP_MIXED, (unsigned)(n + 5u));
            check_case(&wild, n, 0, "mixed");
        }
    }
}

static void run_large(void)
{
    static const size_t kBig[] = { 1024u, 4097u, 8192u };
    fv_mat4 m;
    size_t k;

    printf("test_mat4: 4 - bit-level at n = 1024 / 4097 / 8192\n");
    for (k = 0; k < sizeof(kBig) / sizeof(kBig[0]); ++k) {
        const size_t n = kBig[k];
        make_matrix(&m);
        fill_corpus(g_in, n + GUARD, CORP_FINITE, (unsigned)n);
        check_case(&m, n, 1, "big-finite");
        fill_corpus(g_in, n + GUARD, CORP_EXTREME, (unsigned)n);
        check_case(&m, n, 1, "big-extreme");
        fill_corpus(g_in, n + GUARD, CORP_NAN, (unsigned)n);
        check_case(&m, n, 1, "big-nan");
    }
}

/* The property that actually matters to a caller: the value of an element must not depend on
 * whether it happened to land in the bulk region or in the scalar tail. */
static void run_path_agreement(void)
{
    static const int kCorpus[] = { CORP_FINITE, CORP_EXTREME, CORP_NAN };
    fv_vec4 pair[2], single[1], out_pair[2], out_single[1];
    fv_mat4 m;
    size_t k;
    int c;

    printf("test_mat4: 5 - AVX2 lane vs scalar tail for the same element (path independence)\n");
    for (c = 0; c < 3; ++c) {
        make_matrix(&m);
        for (k = 0; k < 16; ++k) {
            fill_corpus(pair, 2, kCorpus[c], (unsigned)k);
            single[0] = pair[0];
            pair[1] = pair[0];
            fv_mat4_transform_points(&m, pair, out_pair, 2);     /* element 0 through AVX2 */
            fv_mat4_transform_points(&m, single, out_single, 1); /* the same through the tail */
            CHECK(bytes_eq(out_pair, out_single, sizeof(fv_vec4)));
            fv_mat4_transform_dirs(&m, pair, out_pair, 2);
            fv_mat4_transform_dirs(&m, single, out_single, 1);
            CHECK(bytes_eq(out_pair, out_single, sizeof(fv_vec4)));
        }
    }
}

/* The documented boundary, exercised on purpose: one addition whose two addends are NaNs with
 * different bit patterns. Both paths return a NaN -- that part is asserted, and it is the part a
 * caller can rely on; *which* NaN depends on the instruction operand order, which C does not pin
 * down, so a payload difference is counted and printed instead of asserted. Everything else in
 * this file is byte-strict. */
static void run_nan_boundary(void)
{
    fv_mat4 m = fv_mat4_identity();
    fv_vec4 in2[2];
    fv_vec4 via_avx2[2];
    fv_vec4 via_tail[1];

    printf("test_mat4: 6 - documented NaN boundary: two NaN addends, different bit patterns\n");
    m.m[8] = 1.0f;  /* z row = (1*NaN) + (1*(-NaN)) -> NaN + NaN */
    m.m[9] = 1.0f;
    m.m[10] = 1.0f;
    m.m[11] = 0.0f;
    in2[0].x = (fv_real)NAN;
    in2[0].y = -(fv_real)NAN;
    in2[0].z = 0.0f;
    in2[0].w = 0.0f;
    in2[1] = in2[0];
    fv_mat4_transform_points(&m, in2, via_avx2, 2); /* element 0 through the AVX2 kernel */
    fv_mat4_transform_points(&m, in2, via_tail, 1); /* the same element through the scalar tail */
    CHECK(is_nan_bits(bits_u32(via_avx2[0].z)));
    CHECK(is_nan_bits(bits_u32(via_tail[0].z)));
    if (!bytes_eq(via_avx2, via_tail, sizeof(fv_vec4))) {
        ++g_nan_payload;
    }
}

static fv_vec4 mk(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_vec4 r;
    r.x = x;
    r.y = y;
    r.z = z;
    r.w = w;
    return r;
}

static void run_semantics(void)
{
    fv_mat4 id = fv_mat4_identity();
    fv_mat4 t;
    fv_mat4 s;
    fv_vec4 in[4];
    fv_vec4 o[4];
    fv_vec4 o2[4];
    size_t i;

    printf("test_mat4: 7 - semantics: identity / translation / scale / w convention\n");

    for (i = 0; i < 16; ++i) {
        CHECK(bits_eq(id.m[i], fv_mat4_identity().m[i]));
    }
    CHECK(bits_eq(id.m[0], (fv_real)1) && bits_eq(id.m[5], (fv_real)1) &&
          bits_eq(id.m[10], (fv_real)1) && bits_eq(id.m[15], (fv_real)1));

    /* --- identity: out == in bit for bit for ordinary values; w becomes 1 / 0. --- */
    in[0] = mk(1.5f, -2.25f, 0.125f, 9.0f);
    in[1] = mk(1e20f, -1e-20f, 7.0f, -3.0f);
    in[2] = mk(0.0f, -1.0f, 2.0f, 0.0f);
    in[3] = mk(FLT_MAX, FLT_MIN, 0.0f, 1.0f);
    fv_mat4_transform_points(&id, in, o, 4);
    for (i = 0; i < 4; ++i) {
        CHECK(bits_eq(o[i].x, in[i].x));
        CHECK(bits_eq(o[i].y, in[i].y));
        CHECK(bits_eq(o[i].z, in[i].z));
        CHECK(bits_eq(o[i].w, (fv_real)1));
    }
    fv_mat4_transform_dirs(&id, in, o, 4);
    for (i = 0; i < 4; ++i) {
        CHECK(bits_eq(o[i].x, in[i].x));
        CHECK(bits_eq(o[i].y, in[i].y));
        CHECK(bits_eq(o[i].z, in[i].z));
        CHECK(bits_eq(o[i].w, (fv_real)0));
    }

    /* -0.0 differs between the two entries, and that is the reference expression, not a bug:
     *   points: (((-0.0) + (-0.0)) + (-0.0)) + 0.0 = +0.0   (the translation term)
     *   dirs  : (((-0.0) + (-0.0)) + (-0.0))       = -0.0
     * Pinned so that a future "optimization" that hides the difference gets caught. */
    in[0] = mk(-0.0f, -0.0f, -0.0f, 0.0f);
    fv_mat4_transform_points(&id, in, o, 1);
    CHECK(bits_eq(o[0].x, 0.0f));
    CHECK(bits_eq(o[0].x, -0.0f) == 0);
    fv_mat4_transform_dirs(&id, in, o, 1);
    CHECK(bits_eq(o[0].x, -0.0f));

    /* --- pure translation: points move by t, dirs do not move at all. --- */
    t = fv_mat4_identity();
    t.m[3] = 10.0f;
    t.m[7] = -20.0f;
    t.m[11] = 0.5f;
    in[0] = mk(1.0f, 2.0f, 3.0f, 0.0f);
    in[1] = mk(-4.5f, 0.25f, 8.0f, 1.0f);
    fv_mat4_transform_points(&t, in, o, 2);
    CHECK(bits_eq(o[0].x, 11.0f) && bits_eq(o[0].y, -18.0f) && bits_eq(o[0].z, 3.5f));
    CHECK(bits_eq(o[1].x, 5.5f) && bits_eq(o[1].y, -19.75f) && bits_eq(o[1].z, 8.5f));
    CHECK(bits_eq(o[0].w, (fv_real)1) && bits_eq(o[1].w, (fv_real)1));
    fv_mat4_transform_dirs(&t, in, o, 2);
    CHECK(bits_eq(o[0].x, 1.0f) && bits_eq(o[0].y, 2.0f) && bits_eq(o[0].z, 3.0f));
    CHECK(bits_eq(o[1].x, -4.5f) && bits_eq(o[1].y, 0.25f) && bits_eq(o[1].z, 8.0f));
    CHECK(bits_eq(o[0].w, (fv_real)0) && bits_eq(o[1].w, (fv_real)0));

    /* --- pure scale by powers of two: exact, so bit equality is the right assertion. --- */
    s = fv_mat4_identity();
    s.m[0] = 2.0f;
    s.m[5] = 4.0f;
    s.m[10] = 8.0f;
    in[0] = mk(1.5f, -3.0f, 0.25f, 0.0f);
    in[1] = mk(0.0f, 1.0f, -2.0f, 0.0f);
    fv_mat4_transform_points(&s, in, o, 2);
    CHECK(bits_eq(o[0].x, 3.0f) && bits_eq(o[0].y, -12.0f) && bits_eq(o[0].z, 2.0f));
    CHECK(bits_eq(o[1].x, 0.0f) && bits_eq(o[1].y, 4.0f) && bits_eq(o[1].z, -16.0f));
    fv_mat4_transform_dirs(&s, in, o, 2);
    CHECK(bits_eq(o[0].x, 3.0f) && bits_eq(o[0].y, -12.0f) && bits_eq(o[0].z, 2.0f));
    CHECK(bits_eq(o[1].x, 0.0f) && bits_eq(o[1].y, 4.0f) && bits_eq(o[1].z, -16.0f));

    /* --- w is written, never computed: a dot product would give NaN here (0 * inf). --- */
    in[0] = mk((fv_real)INFINITY, (fv_real)NAN, -0x1p-149f, 7.0f);
    in[1] = mk((fv_real)-INFINITY, -0.0f, FLT_MAX, -7.0f);
    fv_mat4_transform_points(&t, in, o, 2);
    CHECK(bits_eq(o[0].w, (fv_real)1));
    CHECK(bits_eq(o[1].w, (fv_real)1));
    fv_mat4_transform_dirs(&t, in, o, 2);
    CHECK(bits_eq(o[0].w, (fv_real)0));
    CHECK(bits_eq(o[1].w, (fv_real)0));

    /* --- a non-affine last row must not change anything: it is never read. --- */
    {
        fv_mat4 a = fv_mat4_identity();
        fv_mat4 b = fv_mat4_identity();
        fv_vec4 src[4], dst_a[4], dst_b[4];
        a.m[0] = 1.25f;
        a.m[4] = -0.5f;
        a.m[8] = 2.0f;
        b.m[0] = 1.25f;
        b.m[4] = -0.5f;
        b.m[8] = 2.0f;
        b.m[12] = (fv_real)INFINITY;
        b.m[13] = (fv_real)NAN;
        b.m[14] = -1.0f;
        b.m[15] = 100.0f;
        (void)o2;
        fill_corpus(src, 4, CORP_EXTREME, 3u);
        fv_mat4_transform_points(&a, src, dst_a, 4);
        fv_mat4_transform_points(&b, src, dst_b, 4);
        CHECK(bytes_eq(dst_a, dst_b, 4 * sizeof(fv_vec4)));
        fv_mat4_transform_dirs(&a, src, dst_a, 4);
        fv_mat4_transform_dirs(&b, src, dst_b, 4);
        CHECK(bytes_eq(dst_a, dst_b, 4 * sizeof(fv_vec4)));
    }

    /* --- n == 0 touches no pointer at all: NULL is legal. --- */
    fv_mat4_transform_points(NULL, NULL, NULL, 0);
    fv_mat4_transform_dirs(NULL, NULL, NULL, 0);
    CHECK(1); /* reaching here means neither call dereferenced anything */
}

int main(void)
{
    poison_all();
    fill_corpus(g_in, CAP + GUARD, CORP_FINITE, 0u);

    run_tails();
    run_large();
    run_path_agreement();
    run_nan_boundary();
    run_semantics();

    printf("test_mat4: NaN payload diffs tolerated (both sides NaN, sections 3 and 6) = %d\n",
           g_nan_payload);
    printf("test_mat4: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
