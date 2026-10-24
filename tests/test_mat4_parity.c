/* 矩阵点变换的跨形态位级对拍（与 tests/test_parity.c、tests/test_vec34_parity.c 同一套做法）。
 *
 * 做法：同一个源文件编 5 份（C / C++ × 默认 inline / FV_USE_EXTERN × 纯 C 回退）+ 1 份敌意
 * 选项的消费者（-march=native -mfma -ffp-contract=fast），每份把所有用例的输出按固定格式写成
 * 二进制记录流，然后逐字节 cmp。同一份源文件保证输入、调用顺序、参数完全一致。
 *
 * 为什么这份夹具必须存在（而不是只靠 tests/test_mat4.cpp）：纯 C 回退形态（-DFV_I_NO_INTRINSICS）
 * 与 FMA 消费者是两条**不同的代码生成**：前者没有 AVX2 内核、全部走标量尾巴；后者命令行里带着
 * -mfma，能穿透 target("avx2") 的只有头文件里那两段 pragma。这两份流只要有一个字节不同，
 * 就说明矩阵内核的位级承诺被编译选项或形态分支破坏了。
 *
 * 记录格式（本机字节序，三个 u32 头 + 原始字节）：
 *   [case_id][n][nbytes] 后跟 payload。payload 是 (n + GUARD) 个 fv_vec4，
 *   调用前 memset 成 0x5A，所以「越界写」与「该写没写」都会让两份文件不同。
 *
 * 语料分三组，每组最多一种 NaN 位模式 —— 两侧同时为 NaN 且位模式不同时，payload 由指令操作数
 * 顺序决定，契约不承诺（见 include/flashvec.h 的 AoS 一节与 tests/test_mat4.cpp 的文件头），
 * 所以那种组合不进这条流，进位级严格相等的那一档。
 *
 * 用法：test_mat4_parity <out.bin>（make test-mat4-inline / test-mat4-fma 调用）
 */
#include "flashvec.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CAPL  4104u
#define GUARD 8u

static fv_vec4 g_in[CAPL + GUARD];
static fv_vec4 g_out[CAPL + GUARD];
static fv_vec4 g_ref[CAPL + GUARD];

static unsigned int g_cases;
static unsigned int g_failures;
static unsigned int g_rng = 0xDEADBEEFu;

static uint32_t rnd_u32(void)
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

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif

/* 参考实现：与头文件里标量尾巴逐字同序 ((m0*x + m1*y) + m2*z) + m3；dirs 没有第四项。 */
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

/* 只有一种 NaN 位模式（+NaN），且不含 inf，不会额外生成 indefinite（0xffc00000） */
static const fv_real kNanMix[] = { (fv_real)NAN, 0.0f, 1.0f, -1.0f, 0.5f, 1.0e30f, -1.0e30f };
#define NNANMIX ((unsigned)(sizeof(kNanMix) / sizeof(kNanMix[0])))

enum { CORP_FINITE = 0, CORP_EXTREME = 1, CORP_NAN = 2 };

static void fill_corpus(fv_vec4 *v, size_t n, int corpus, unsigned int salt)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        v[i].w = (fv_real)((i & 1u) != 0u ? -4.5 : 6.25); /* 两条入口都忽略输入 w */
        if (corpus == CORP_FINITE) {
            v[i].x = rnd_moderate();
            v[i].y = rnd_moderate();
            v[i].z = rnd_moderate();
        } else if (corpus == CORP_EXTREME) {
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

static void make_wild_matrix(fv_mat4 *m)
{
    make_matrix(m);
    m->m[1] = (fv_real)INFINITY;
    m->m[6] = -0.0f;
    m->m[10] = 0x1p-149f;
    m->m[3] = -0.0f;
    m->m[7] = FLT_MAX;
}

static void emit(FILE *f, uint32_t id, size_t n)
{
    const uint32_t hdr[3] = { id, (uint32_t)n, (uint32_t)((n + GUARD) * sizeof(fv_vec4)) };
    if (fwrite(hdr, sizeof(hdr), 1, f) != 1u) {
        ++g_failures;
    }
    if (fwrite(g_out, sizeof(fv_vec4), n + GUARD, f) != n + GUARD) {
        ++g_failures;
    }
}

/* 两块输出都在调用前毒化成 0x5A：比较 (n + GUARD) 个元素，所以「库越界写」与
 * 「参考实现越界写」都会被抓到，两边都必须只碰前 n 个。 */
static void poison_both(void)
{
    memset(g_out, 0x5A, sizeof(g_out));
    memset(g_ref, 0x5A, sizeof(g_ref));
}

static void run_case(FILE *f, const fv_mat4 *m, size_t n, int corpus, uint32_t id)
{
    fill_corpus(g_in, n + GUARD, corpus, (unsigned)n);

    poison_both();
    fv_mat4_transform_points(m, g_in, g_out, n);
    ref_points(m, g_in, g_ref, n);
    if (memcmp(g_out, g_ref, (n + GUARD) * sizeof(fv_vec4)) != 0) {
        ++g_failures;
        printf("FAIL points n=%u corpus=%d id=%u\n", (unsigned)n, corpus, (unsigned)id);
    }
    emit(f, id, n);

    poison_both();
    fv_mat4_transform_dirs(m, g_in, g_out, n);
    ref_dirs(m, g_in, g_ref, n);
    if (memcmp(g_out, g_ref, (n + GUARD) * sizeof(fv_vec4)) != 0) {
        ++g_failures;
        printf("FAIL dirs   n=%u corpus=%d id=%u\n", (unsigned)n, corpus, (unsigned)id);
    }
    emit(f, id + 1u, n);

    g_cases += 2u;
}

int main(int argc, char **argv)
{
    static const size_t kTail[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                                    32, 33, 34, 35, 36, 37, 38, 39, 40 };
    static const size_t kBig[] = { 63u, 64u, 65u, 1024u, 4097u };
    fv_mat4 normal;
    fv_mat4 wild;
    FILE *f;
    size_t k;
    int corpus;
    uint32_t id = 1u;

    if (argc < 2) {
        printf("usage: %s <out.bin>\n", argv[0]);
        return 2;
    }
    f = fopen(argv[1], "wb");
    if (f == NULL) {
        printf("cannot open %s\n", argv[1]);
        return 2;
    }

    make_matrix(&normal);
    make_wild_matrix(&wild);

    for (corpus = 0; corpus < 3; ++corpus) {
        /* Corpus 2 carries a NaN operand, so it must not meet an infinite coefficient either:
         * inf * 0 yields the hardware indefinite (0xffc00000) which, added to a propagated +NaN,
         * is the "two different NaN patterns in one add" case whose payload C does not pin down.
         * Corpora 0/1 contain no NaN operand, so every NaN they can produce is that same
         * indefinite and the strict byte comparison stays meaningful. */
        const fv_mat4 *big_m = (corpus == CORP_NAN) ? &normal : &wild;
        for (k = 0; k < sizeof(kTail) / sizeof(kTail[0]); ++k) {
            run_case(f, &normal, kTail[k], corpus, id);
            id += 2u;
        }
        for (k = 0; k < sizeof(kBig) / sizeof(kBig[0]); ++k) {
            run_case(f, big_m, kBig[k], corpus, id);
            id += 2u;
        }
    }

    if (fclose(f) != 0) {
        ++g_failures;
    }
    printf("mat4_parity: cases=%u failures=%u\n", g_cases, g_failures);
    return g_failures == 0 ? 0 : 1;
}
