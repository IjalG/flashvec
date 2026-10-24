/* 矩阵点变换在纯 C 编译单元里的验收（与 tests/test_batch_c.c 同一分工：C 只验「声明可用、
 * 语义正确、与标量参考逐位一致」这几件事，深度的形态对拍在 tests/test_mat4.cpp 与
 * tests/test_mat4_parity.c 里）。
 *
 * 为什么 C 单元要单独测：fv_vec4 按值返回、fv_mat4 按指针传、__restrict 参数 —— 这些在 C 与
 * C++ 里的 ABI 细节不同（结构体返回的寄存器约定、extern "C" 包装），库的头文件是 C11 接口，
 * 所以必须有一个只用 C 编译器的调用方来证明它真的能被 C 消费。
 *
 * 参考实现与头文件里的标量尾巴用同一个求值顺序 ((m0*x + m1*y) + m2*z) + m3（dirs 少最后一项），
 * 并且同样关掉 fp-contract（否则 -ffp-contract=fast 会把参考自己的 mul+add 融成 FMA，
 * 比的就不是库了）。语料按「最多一种 NaN 位模式」分组，理由见 tests/test_mat4.cpp 文件头。
 *
 * 运行：make test-mat4-c
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(cond)                                                \
    do {                                                           \
        ++g_checks;                                                \
        if (!(cond)) {                                             \
            ++g_failures;                                          \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
        }                                                          \
    } while (0)

#define CAP   2048
#define GUARD 2

static fv_vec4 g_in[CAP + GUARD];
static fv_vec4 g_got[CAP + GUARD];
static fv_vec4 g_ref[CAP + GUARD];

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

static const fv_real kExtreme[] = {
    0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f,
    FLT_MAX, -FLT_MAX, FLT_MIN, -FLT_MIN, 0x1p-149f, -0x1p-149f,
    (fv_real)INFINITY, (fv_real)-INFINITY
};
#define NEXTREME ((unsigned)(sizeof(kExtreme) / sizeof(kExtreme[0])))

static const fv_real kNanMix[] = { (fv_real)NAN, 0.0f, 1.0f, -1.0f, 0.5f, 1.0e30f };
#define NNANMIX ((unsigned)(sizeof(kNanMix) / sizeof(kNanMix[0])))

static unsigned int g_rng = 0x1234567u;

static unsigned int rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static fv_real rnd_moderate(void)
{
    return (fv_real)((int)(rnd_u32() % 2001u) - 1000) / 16.0f;
}

/* 0 = finite, 1 = extreme (no NaN), 2 = single NaN pattern */
static void fill_corpus(fv_vec4 *v, size_t n, int corpus, unsigned int salt)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        v[i].w = (fv_real)5.5; /* ignored by both entries */
        if (corpus == 0) {
            v[i].x = rnd_moderate();
            v[i].y = rnd_moderate();
            v[i].z = rnd_moderate();
        } else if (corpus == 1) {
            v[i].x = kExtreme[(i + salt) % NEXTREME];
            v[i].y = kExtreme[(i * 3u + salt + 1u) % NEXTREME];
            v[i].z = kExtreme[(i * 7u + salt + 2u) % NEXTREME];
        } else {
            v[i].x = kNanMix[(i + salt) % NNANMIX];
            v[i].y = kNanMix[(i * 3u + salt + 1u) % NNANMIX];
            v[i].z = kNanMix[(i * 7u + salt + 2u) % NNANMIX];
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

static fv_vec4 mk(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_vec4 r;
    r.x = x;
    r.y = y;
    r.z = z;
    r.w = w;
    return r;
}

int main(void)
{
    static const size_t kBig[] = { 1024u, 2047u };
    fv_mat4 m;
    fv_mat4 id;
    fv_mat4 t;
    fv_mat4 s;
    fv_vec4 in[4], o[4];
    size_t n;
    size_t k;
    int corpus;

    /* ---- 1：n = 0..40 全尾长，三种语料，逐位对拍 ---- */
    for (corpus = 0; corpus < 3; ++corpus) {
        make_matrix(&m);
        for (n = 0; n <= 40; ++n) {
            fill_corpus(g_in, n + GUARD, corpus, (unsigned)n);
            memset(g_got, 0x5A, sizeof(g_got));
            fv_mat4_transform_points(&m, g_in, g_got, n);
            ref_points(&m, g_in, g_ref, n);
            CHECK(memcmp(g_got, g_ref, n * sizeof(fv_vec4)) == 0);
            memset(g_got, 0x5A, sizeof(g_got));
            fv_mat4_transform_dirs(&m, g_in, g_got, n);
            ref_dirs(&m, g_in, g_ref, n);
            CHECK(memcmp(g_got, g_ref, n * sizeof(fv_vec4)) == 0);
        }
    }

    /* ---- 2：大 n（跨主循环 + 尾巴） ---- */
    make_matrix(&m);
    for (k = 0; k < sizeof(kBig) / sizeof(kBig[0]); ++k) {
        const size_t big = kBig[k];
        fill_corpus(g_in, big + GUARD, 0, 0u);
        fv_mat4_transform_points(&m, g_in, g_got, big);
        ref_points(&m, g_in, g_ref, big);
        CHECK(memcmp(g_got, g_ref, big * sizeof(fv_vec4)) == 0);
        fill_corpus(g_in, big + GUARD, 1, (unsigned)big);
        fv_mat4_transform_dirs(&m, g_in, g_got, big);
        ref_dirs(&m, g_in, g_ref, big);
        CHECK(memcmp(g_got, g_ref, big * sizeof(fv_vec4)) == 0);
    }

    /* ---- 3：语义 ---- */
    id = fv_mat4_identity();
    in[0] = mk(1.5f, -2.25f, 0.125f, 9.0f);
    in[1] = mk(-1e20f, 3.0f, -0.5f, -1.0f);
    fv_mat4_transform_points(&id, in, o, 2);
    CHECK(bits_eq(o[0].x, 1.5f) && bits_eq(o[0].y, -2.25f) && bits_eq(o[0].z, 0.125f));
    CHECK(bits_eq(o[1].x, -1e20f) && bits_eq(o[1].y, 3.0f) && bits_eq(o[1].z, -0.5f));
    CHECK(bits_eq(o[0].w, (fv_real)1) && bits_eq(o[1].w, (fv_real)1));
    fv_mat4_transform_dirs(&id, in, o, 2);
    CHECK(bits_eq(o[0].x, 1.5f) && bits_eq(o[0].y, -2.25f) && bits_eq(o[0].z, 0.125f));
    CHECK(bits_eq(o[0].w, (fv_real)0) && bits_eq(o[1].w, (fv_real)0));

    t = fv_mat4_identity();
    t.m[3] = 10.0f;
    t.m[7] = -20.0f;
    t.m[11] = 0.5f;
    in[0] = mk(1.0f, 2.0f, 3.0f, 0.0f);
    fv_mat4_transform_points(&t, in, o, 1);
    CHECK(bits_eq(o[0].x, 11.0f) && bits_eq(o[0].y, -18.0f) && bits_eq(o[0].z, 3.5f));
    fv_mat4_transform_dirs(&t, in, o, 1);
    CHECK(bits_eq(o[0].x, 1.0f) && bits_eq(o[0].y, 2.0f) && bits_eq(o[0].z, 3.0f));

    s = fv_mat4_identity();
    s.m[0] = 2.0f;
    s.m[5] = 4.0f;
    s.m[10] = 8.0f;
    in[0] = mk(1.5f, -3.0f, 0.25f, 0.0f);
    fv_mat4_transform_points(&s, in, o, 1);
    CHECK(bits_eq(o[0].x, 3.0f) && bits_eq(o[0].y, -12.0f) && bits_eq(o[0].z, 2.0f));

    /* n == 0：不读不写任何指针 */
    fv_mat4_transform_points(NULL, NULL, NULL, 0);
    fv_mat4_transform_dirs(NULL, NULL, NULL, 0);
    CHECK(1);

    printf("C mat4 API ok: points dirs identity translate scale w  (%d checks, %d failures)\n",
           g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
