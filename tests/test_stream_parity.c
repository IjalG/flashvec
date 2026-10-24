/* Bit-level cross-check of the 13 _stream (NT store) entry points added for vec2 / vec3 / vec4 / mat4.
 *
 * Two layers, both in this one file (so every build form runs the identical corpus in the identical order):
 *   1) in-process: for every case the stream output is compared with the non-stream output over n + GUARD elements (both
 *      buffers poisoned first), so the NT body, the head peel, the tail and any out-of-range store are all caught.
 *      This runs for BOTH corpora below.
 *   2) cross-form: the stream outputs of the recorded corpus are appended to a binary record stream, and the Makefile
 *      cmp-s the streams of C inline / C++ inline / C extern (FV_USE_EXTERN) / C++ extern / C pure fallback
 *      (-DFV_I_NO_INTRINSICS) and of the FMA consumer (-march=native -mfma -ffp-contract=fast). Layer 2 is what proves the
 *      pragma still blocks the FMA a consumer's command line brings in: lerp / madd / normalize / mat4 are mul+add, so a
 *      missing pragma shows up as a different record stream, never as a tolerance.
 *
 * Two corpora, and why they are split:
 *   - RECORDED (finite values, the same generator tests/test_parity.c uses): n = 0..40 (every tail length of _x4 / _x2 /
 *     _soa) + 63 + 100 + 1023 + 4097; output pointer shifted by 0, 1, 3, 7, 8, 16 elements (0..7 walk through every 32 B
 *     phase of the stride). Every form must produce the same bytes.
 *   - EDGE VALUES (in-process only): -0.0f, +-inf, NaN, FLT_MIN / FLT_MAX, denormals, exact zero vectors and near-zero
 *     vectors, over n = 0..40 + the same large lengths and all 12 phases. These stay out of the cross-form stream on
 *     purpose: with inf / NaN operands the sign bit of a produced NaN depends on the operand order of the generated
 *     mul/add, and the pure-C fallback lets GCC auto-vectorize the scalar loop (mulps / addps), which can pick a different
 *     operand order than the hand-written kernel. Both results are NaN and neither value is part of any promise, but a
 *     cross-form cmp would report it. The stream-vs-non-stream comparison inside one build is unaffected (identical code
 *     on both sides), so the edge values are still verified where the contract actually lives.
 *
 * Record format (native endianness): [case_id][n][off][nbytes] + payload, the same shape tests/test_parity.c uses.
 *
 * Usage: test_stream_parity <main.bin>
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>   /* INFINITY / NAN: the FV_USE_EXTERN form of the header pulls in no math.h by design */
#include <stdio.h>
#include <string.h>

#define CAP 8300u
#define GUARD 8u

static fv_vec2 g_v2a[CAP + GUARD], g_v2b[CAP + GUARD];
static fv_vec2 g_v2r[CAP + GUARD], g_v2s[CAP + GUARD];
static fv_vec4 g_v4a[CAP + GUARD], g_v4b[CAP + GUARD];
static fv_vec4 g_v4r[CAP + GUARD], g_v4s[CAP + GUARD];
static float g_ax[CAP + GUARD], g_ay[CAP + GUARD], g_az[CAP + GUARD], g_aw[CAP + GUARD];
static float g_bx[CAP + GUARD], g_by[CAP + GUARD], g_bz[CAP + GUARD], g_bw[CAP + GUARD];
static float g_r1[CAP + GUARD], g_r2[CAP + GUARD], g_r3[CAP + GUARD], g_r4[CAP + GUARD];
static float g_s1[CAP + GUARD], g_s2[CAP + GUARD], g_s3[CAP + GUARD], g_s4[CAP + GUARD];
static fv_mat4 g_m;

static FILE *g_f;
static unsigned int g_cases;
static unsigned int g_checks;
static unsigned int g_fail;
static unsigned int g_rng = 0x9E3779B9u;

static void chk(int ok, const char *what, size_t n, size_t off)
{
    ++g_checks;
    if (ok == 0) {
        ++g_fail;
        printf("FAIL %-40s n=%u off=%u\n", what, (unsigned int)n, (unsigned int)off);
    }
}

static void rec(unsigned int id, size_t n, size_t off, const void *p, size_t bytes)
{
    unsigned int hdr[4];
    hdr[0] = id;
    hdr[1] = (unsigned int)n;
    hdr[2] = (unsigned int)off;
    hdr[3] = (unsigned int)bytes;
    fwrite(hdr, sizeof(hdr), 1, g_f);
    fwrite(p, 1, bytes, g_f);
    ++g_cases;
}

static void poison(void *p, size_t bytes)
{
    memset(p, 0x5A, bytes);
}

static unsigned int rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/* A finite float with a random sign and a random exponent in 1..254 (never 0 / inf / NaN): the generator
 * tests/test_parity.c uses for its cross-form corpus. */
static float rnd_float(void)
{
    const unsigned int b = rnd_u32();
    const unsigned int bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    float f;
    memcpy(&f, &bits, sizeof(f));
    return f;
}

/* Values a vector kernel can mishandle: signed zero, infinities, NaN, the float range extremes, denormals. */
static float edge_value(size_t i)
{
    static const float kEdge[16] = {
        0.0f, -0.0f, 1.0f, -1.0f,
        (float)INFINITY, -(float)INFINITY, (float)NAN, 3.25f,
        FLT_MIN, -FLT_MIN, FLT_MAX, -FLT_MAX,
        1.0e-30f, -1.0e-30f, 0.5f, -2.5f
    };
    return kEdge[i % 16u];
}

static void fill_inputs(size_t n, int edge)
{
    size_t i;
    for (i = 0; i < n + GUARD; ++i) {
        const float ax = (edge != 0) ? edge_value(i) : rnd_float();
        const float ay = (edge != 0) ? ((i % 5u == 0u) ? 0.0f : edge_value(i * 3u + 1u)) : rnd_float();
        const float az = (edge != 0) ? ((i % 7u == 0u) ? 1.0e-30f : edge_value(i * 5u + 2u)) : rnd_float();
        const float aw = (edge != 0) ? edge_value(i * 7u + 3u) : rnd_float();
        const float bx = (edge != 0) ? edge_value(i * 11u + 4u) : rnd_float();
        const float by = (edge != 0) ? edge_value(i * 13u + 5u) : rnd_float();
        const float bz = (edge != 0) ? edge_value(i * 17u + 6u) : rnd_float();
        const float bw = (edge != 0) ? edge_value(i * 19u + 7u) : rnd_float();
        g_v2a[i].x = ax;
        g_v2a[i].y = ay;
        g_v2b[i].x = bx;
        g_v2b[i].y = by;
        g_v4a[i].x = ax;
        g_v4a[i].y = ay;
        g_v4a[i].z = az;
        g_v4a[i].w = aw;
        g_v4b[i].x = bx;
        g_v4b[i].y = by;
        g_v4b[i].z = bz;
        g_v4b[i].w = bw;
        g_ax[i] = ax;
        g_ay[i] = ay;
        g_az[i] = az;
        g_aw[i] = aw;
        g_bx[i] = bx;
        g_by[i] = by;
        g_bz[i] = bz;
        g_bw[i] = bw;
    }
}

static void set_matrix(void)
{
    static const float kM[12] = { 1.5f, -0.5f, 0.25f, 7.0f,
                                  -2.0f, 3.25f, 0.0f, -1.0f,
                                  1.0e-30f, -0.0f, 0.5f, 2.0f };
    size_t i;
    for (i = 0; i < 12u; ++i) {
        g_m.m[i] = kM[i];
    }
    g_m.m[12] = 0.0f;
    g_m.m[13] = 0.0f;
    g_m.m[14] = 0.0f;
    g_m.m[15] = 1.0f;
}

/* One (n, phase) case: all 13 new entries, stream vs non-stream, plus the recorded stream when record != 0. */
static void run_case(size_t n, size_t off, const float s, const float t, int record)
{
    const size_t b2 = (n + 4u) * sizeof(fv_vec2);
    const size_t b4 = (n + 4u) * sizeof(fv_vec4);
    const size_t bf = (n + 4u) * sizeof(float);
    fv_vec2 *const v2r = g_v2r + off;
    fv_vec2 *const v2s = g_v2s + off;
    fv_vec4 *const v4r = g_v4r + off;
    fv_vec4 *const v4s = g_v4s + off;
    float *const f1r = g_r1 + off;
    float *const f2r = g_r2 + off;
    float *const f3r = g_r3 + off;
    float *const f4r = g_r4 + off;
    float *const f1s = g_s1 + off;
    float *const f2s = g_s2 + off;
    float *const f3s = g_s3 + off;
    float *const f4s = g_s4 + off;
    const unsigned int base = (unsigned int)(0x7000u + (off * 61u + n % 97u) * 16u);

    /* vec2 _x4 */
    poison(g_v2r, sizeof(g_v2r));
    poison(g_v2s, sizeof(g_v2s));
    fv_vec2_lerp_x4(g_v2a, g_v2b, t, v2r, n);
    fv_vec2_lerp_x4_stream(g_v2a, g_v2b, t, v2s, n);
    chk(memcmp(v2r, v2s, b2) == 0, "lerp_x4_stream == lerp_x4", n, off);
    if (record != 0) {
        rec(base + 1u, n, off, v2s, b2);
    }

    poison(g_v2r, sizeof(g_v2r));
    poison(g_v2s, sizeof(g_v2s));
    fv_vec2_madd_x4(g_v2a, g_v2b, s, v2r, n);
    fv_vec2_madd_x4_stream(g_v2a, g_v2b, s, v2s, n);
    chk(memcmp(v2r, v2s, b2) == 0, "madd_x4_stream == madd_x4", n, off);
    if (record != 0) {
        rec(base + 2u, n, off, v2s, b2);
    }

    /* vec3 SoA */
    poison(g_r1, sizeof(g_r1));
    poison(g_r2, sizeof(g_r2));
    poison(g_r3, sizeof(g_r3));
    poison(g_s1, sizeof(g_s1));
    poison(g_s2, sizeof(g_s2));
    poison(g_s3, sizeof(g_s3));
    fv_vec3_add_soa(g_ax, g_ay, g_az, g_bx, g_by, g_bz, f1r, f2r, f3r, n);
    fv_vec3_add_soa_stream(g_ax, g_ay, g_az, g_bx, g_by, g_bz, f1s, f2s, f3s, n);
    chk(memcmp(f1r, f1s, bf) == 0 && memcmp(f2r, f2s, bf) == 0 && memcmp(f3r, f3s, bf) == 0,
        "add_soa_stream (vec3) == add_soa", n, off);
    if (record != 0) {
        rec(base + 3u, n, off, f1s, bf);
        rec(base + 3u, n, off, f2s, bf);
        rec(base + 3u, n, off, f3s, bf);
    }

    poison(g_r1, sizeof(g_r1));
    poison(g_r2, sizeof(g_r2));
    poison(g_r3, sizeof(g_r3));
    poison(g_s1, sizeof(g_s1));
    poison(g_s2, sizeof(g_s2));
    poison(g_s3, sizeof(g_s3));
    fv_vec3_scale_soa(g_ax, g_ay, g_az, s, f1r, f2r, f3r, n);
    fv_vec3_scale_soa_stream(g_ax, g_ay, g_az, s, f1s, f2s, f3s, n);
    chk(memcmp(f1r, f1s, bf) == 0 && memcmp(f2r, f2s, bf) == 0 && memcmp(f3r, f3s, bf) == 0,
        "scale_soa_stream (vec3) == scale_soa", n, off);
    if (record != 0) {
        rec(base + 4u, n, off, f1s, bf);
        rec(base + 4u, n, off, f2s, bf);
        rec(base + 4u, n, off, f3s, bf);
    }

    poison(g_r1, sizeof(g_r1));
    poison(g_r2, sizeof(g_r2));
    poison(g_r3, sizeof(g_r3));
    poison(g_s1, sizeof(g_s1));
    poison(g_s2, sizeof(g_s2));
    poison(g_s3, sizeof(g_s3));
    fv_vec3_normalize_soa(g_ax, g_ay, g_az, f1r, f2r, f3r, n);
    fv_vec3_normalize_soa_stream(g_ax, g_ay, g_az, f1s, f2s, f3s, n);
    chk(memcmp(f1r, f1s, bf) == 0 && memcmp(f2r, f2s, bf) == 0 && memcmp(f3r, f3s, bf) == 0,
        "normalize_soa_stream (vec3) == normalize_soa", n, off);
    if (record != 0) {
        rec(base + 5u, n, off, f1s, bf);
        rec(base + 5u, n, off, f2s, bf);
        rec(base + 5u, n, off, f3s, bf);
    }

    /* vec4 _x2 */
    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_vec4_add_x2(g_v4a, g_v4b, v4r, n);
    fv_vec4_add_x2_stream(g_v4a, g_v4b, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "add_x2_stream == add_x2", n, off);
    if (record != 0) {
        rec(base + 6u, n, off, v4s, b4);
    }

    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_vec4_sub_x2(g_v4a, g_v4b, v4r, n);
    fv_vec4_sub_x2_stream(g_v4a, g_v4b, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "sub_x2_stream == sub_x2", n, off);
    if (record != 0) {
        rec(base + 7u, n, off, v4s, b4);
    }

    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_vec4_scale_x2(g_v4a, s, v4r, n);
    fv_vec4_scale_x2_stream(g_v4a, s, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "scale_x2_stream == scale_x2", n, off);
    if (record != 0) {
        rec(base + 8u, n, off, v4s, b4);
    }

    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_vec4_lerp_x2(g_v4a, g_v4b, t, v4r, n);
    fv_vec4_lerp_x2_stream(g_v4a, g_v4b, t, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "lerp_x2_stream == lerp_x2", n, off);
    if (record != 0) {
        rec(base + 9u, n, off, v4s, b4);
    }

    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_vec4_madd_x2(g_v4a, g_v4b, s, v4r, n);
    fv_vec4_madd_x2_stream(g_v4a, g_v4b, s, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "madd_x2_stream == madd_x2", n, off);
    if (record != 0) {
        rec(base + 10u, n, off, v4s, b4);
    }

    /* vec4 SoA normalize */
    poison(g_r1, sizeof(g_r1));
    poison(g_r2, sizeof(g_r2));
    poison(g_r3, sizeof(g_r3));
    poison(g_r4, sizeof(g_r4));
    poison(g_s1, sizeof(g_s1));
    poison(g_s2, sizeof(g_s2));
    poison(g_s3, sizeof(g_s3));
    poison(g_s4, sizeof(g_s4));
    fv_vec4_normalize_soa(g_ax, g_ay, g_az, g_aw, f1r, f2r, f3r, f4r, n);
    fv_vec4_normalize_soa_stream(g_ax, g_ay, g_az, g_aw, f1s, f2s, f3s, f4s, n);
    chk(memcmp(f1r, f1s, bf) == 0 && memcmp(f2r, f2s, bf) == 0 && memcmp(f3r, f3s, bf) == 0
            && memcmp(f4r, f4s, bf) == 0,
        "normalize_soa_stream (vec4) == normalize_soa", n, off);
    if (record != 0) {
        rec(base + 11u, n, off, f1s, bf);
        rec(base + 11u, n, off, f2s, bf);
        rec(base + 11u, n, off, f3s, bf);
        rec(base + 11u, n, off, f4s, bf);
    }

    /* mat4 */
    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_mat4_transform_points(&g_m, g_v4a, v4r, n);
    fv_mat4_transform_points_stream(&g_m, g_v4a, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "mat4_transform_points_stream", n, off);
    if (record != 0) {
        rec(base + 12u, n, off, v4s, b4);
    }

    poison(g_v4r, sizeof(g_v4r));
    poison(g_v4s, sizeof(g_v4s));
    fv_mat4_transform_dirs(&g_m, g_v4a, v4r, n);
    fv_mat4_transform_dirs_stream(&g_m, g_v4a, v4s, n);
    chk(memcmp(v4r, v4s, b4) == 0, "mat4_transform_dirs_stream", n, off);
    if (record != 0) {
        rec(base + 13u, n, off, v4s, b4);
    }
}

int main(int argc, char **argv)
{
    static const size_t kBigLens[] = { 63u, 100u, 1023u, 4097u };
    static const size_t kRecOffs[] = { 0u, 1u, 3u, 7u, 8u, 16u };
    static const size_t kEdgeOffs[] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 16u, 32u, 64u };
    static const float kS[6] = { 0.0f, -0.0f, 1.0f, -3.25f, 1.0e-30f, 0.5f };
    static const float kT[8] = { 0.0f, -0.0f, 1.0f, 0.5f, -1.0f, 1.0e-30f,
                                 (float)INFINITY, (float)NAN };
    size_t lens[64];
    int nlens = 0;
    int li;
    int oi;
    size_t i;

    if (argc < 2) {
        printf("usage: %s <main.bin>\n", argv[0]);
        return 2;
    }
    g_f = fopen(argv[1], "wb");
    if (g_f == NULL) {
        printf("cannot open %s\n", argv[1]);
        return 2;
    }

    for (i = 0; i <= 40u; ++i) {
        lens[nlens++] = i;
    }
    for (i = 0; i < sizeof(kBigLens) / sizeof(kBigLens[0]); ++i) {
        lens[nlens++] = kBigLens[i];
    }
    set_matrix();

    /* Corpus 1: finite values, recorded -> the five / six build forms must agree byte for byte. */
    for (li = 0; li < nlens; ++li) {
        const size_t n = lens[li];
        fill_inputs(n, 0);
        for (oi = 0; oi < (int)(sizeof(kRecOffs) / sizeof(kRecOffs[0])); ++oi) {
            run_case(n, kRecOffs[oi], kS[(size_t)(li + oi) % 6u], kT[(size_t)(li + oi) % 8u], 1);
        }
    }

    /* Corpus 2: edge values, in-process only (see the note at the top of this file). */
    {
        unsigned int edge_checks = 0;
        for (li = 0; li < nlens; ++li) {
            const size_t n = lens[li];
            fill_inputs(n, 1);
            for (oi = 0; oi < (int)(sizeof(kEdgeOffs) / sizeof(kEdgeOffs[0])); ++oi) {
                const unsigned int before = g_checks;
                run_case(n, kEdgeOffs[oi], kS[(size_t)(li + oi) % 6u],
                         kT[(size_t)(li + oi) % 8u], 0);
                edge_checks += g_checks - before;
            }
        }
        printf("stream_parity: edge-value corpus %u checks (in-process, not recorded)\n", edge_checks);
    }

    /* The 32 B boundary is unreachable for these phases: the kernels must fall back to normal stores, bit-identically. */
    {
        const size_t n = 200u;
        const size_t bytes = (n + 1u) * sizeof(float);
        float *const px = g_r1;
        float *const py = g_r2 + 4;
        float *const pz = g_r3 + 8;
        float *const pw = g_r4 + 12;
        float *const qx = g_s1;
        float *const qy = g_s2 + 4;
        float *const qz = g_s3 + 8;
        float *const qw = g_s4 + 12;

        fill_inputs(n, 1);
        poison(g_r1, sizeof(g_r1));
        poison(g_r2, sizeof(g_r2));
        poison(g_r3, sizeof(g_r3));
        poison(g_s1, sizeof(g_s1));
        poison(g_s2, sizeof(g_s2));
        poison(g_s3, sizeof(g_s3));
        fv_vec3_normalize_soa(g_ax, g_ay, g_az, px, py, pz, n);
        fv_vec3_normalize_soa_stream(g_ax, g_ay, g_az, qx, qy, qz, n);
        chk(memcmp(px, qx, bytes) == 0 && memcmp(py, qy, bytes) == 0 && memcmp(pz, qz, bytes) == 0,
            "normalize_soa_stream (vec3, mismatched phases)", n, 4u);

        poison(g_r1, sizeof(g_r1));
        poison(g_r2, sizeof(g_r2));
        poison(g_r3, sizeof(g_r3));
        poison(g_r4, sizeof(g_r4));
        poison(g_s1, sizeof(g_s1));
        poison(g_s2, sizeof(g_s2));
        poison(g_s3, sizeof(g_s3));
        poison(g_s4, sizeof(g_s4));
        fv_vec4_normalize_soa(g_ax, g_ay, g_az, g_aw, px, py, pz, pw, n);
        fv_vec4_normalize_soa_stream(g_ax, g_ay, g_az, g_aw, qx, qy, qz, qw, n);
        chk(memcmp(px, qx, bytes) == 0 && memcmp(py, qy, bytes) == 0 && memcmp(pz, qz, bytes) == 0
                && memcmp(pw, qw, bytes) == 0,
            "normalize_soa_stream (vec4, mismatched phases)", n, 4u);
    }

    /* n == 0 must touch nothing: NULL is a legal argument for every one of the 13 entries. */
    fv_vec2_lerp_x4_stream(NULL, NULL, 1.0f, NULL, 0);
    fv_vec2_madd_x4_stream(NULL, NULL, 1.0f, NULL, 0);
    fv_vec3_add_soa_stream(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec3_scale_soa_stream(NULL, NULL, NULL, 1.0f, NULL, NULL, NULL, 0);
    fv_vec3_normalize_soa_stream(NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec4_add_x2_stream(NULL, NULL, NULL, 0);
    fv_vec4_sub_x2_stream(NULL, NULL, NULL, 0);
    fv_vec4_scale_x2_stream(NULL, 1.0f, NULL, 0);
    fv_vec4_lerp_x2_stream(NULL, NULL, 1.0f, NULL, 0);
    fv_vec4_madd_x2_stream(NULL, NULL, 1.0f, NULL, 0);
    fv_vec4_normalize_soa_stream(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_mat4_transform_points_stream(NULL, NULL, NULL, 0);
    fv_mat4_transform_dirs_stream(NULL, NULL, NULL, 0);
    ++g_checks;

    fclose(g_f);
    printf("stream_parity: cases=%u  checks=%u  failures=%u\n", g_cases, g_checks, g_fail);
    return g_fail == 0u ? 0 : 1;
}
