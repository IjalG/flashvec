/* P2（内存流量入口）的专项验收：tests/test_parity.c 负责「与旧入口逐位相同」的主对拍，
 * 这个文件补三类主对拍不容易命中的事：
 *
 *   1) 对齐 peel：把输出指针逐个元素平移（fv_vec2 侧 0..3 个元素、SoA 侧 0..7 个 float），
 *      再加一个「4 B 对齐、8 B 不齐」的 fv_vec2 数组 —— 让
 *      头部普通 store / NT 主体 / 尾巴 / storeu 退化路径 四条都被真的走到。
 *      任何一条写错位置或写错值，memcmp 都会连同 GUARD 一起报出来。
 *   2) _mask_inplace 的「不写」语义：未命中元素必须**逐字节**保持原值（不是写回 a[i] 之外
 *      的什么东西），命中元素必须与标量加法逐位相同，GUARD 区一个字节都不能动。
 *   3) sfence：NT store 是弱序写，写完必须 sfence 才对其它核可见。这里开两个线程做
 *      写-读交接，读线程读到的每个元素都必须已经是新值。
 *
 * 编译（严格档）：cc -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -O2 -Iinclude -pthread
 */
#include "flashvec.h"

#include <pthread.h>
#include <stdio.h>
#include <string.h>

#define CAP 4104u
#define GUARD 8u

static fv_vec2 g_a[CAP + GUARD], g_b[CAP + GUARD];
static fv_vec2 g_ref[CAP + GUARD], g_out[CAP + GUARD];
static float g_ax[CAP + GUARD], g_ay[CAP + GUARD], g_bx[CAP + GUARD], g_by[CAP + GUARD];
static float g_rx[CAP + GUARD], g_ry[CAP + GUARD], g_sx[CAP + GUARD], g_sy[CAP + GUARD];
static unsigned char g_mask[CAP + GUARD];

static unsigned int g_rng = 0x9E3779B9u;
static int g_checks;
static int g_failures;

static unsigned int rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static float rnd_float(void)
{
    const unsigned int b = rnd_u32();
    const unsigned int bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    float f;
    memcpy(&f, &bits, 4);
    return f;
}

#define CHECK(cond, what, n, off)                                              \
    do {                                                                       \
        ++g_checks;                                                            \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            printf("FAIL %s  n=%u off=%u  (%s:%d)\n", (what), (unsigned)(n),   \
                   (unsigned)(off), __FILE__, __LINE__);                       \
        }                                                                      \
    } while (0)

static void fill_inputs(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_a[i].x = (i % 7u == 0u) ? 0.0f : rnd_float();
        g_a[i].y = (i % 11u == 0u) ? -0.0f : rnd_float();
        g_b[i].x = rnd_float();
        g_b[i].y = rnd_float();
        g_ax[i] = g_a[i].x;
        g_ay[i] = g_a[i].y;
        g_bx[i] = g_b[i].x;
        g_by[i] = g_b[i].y;
    }
}

static void poison(void *p, size_t bytes)
{
    memset(p, 0x5A, bytes);
}

/* 1) peel / 尾巴 / 退化路径：输出指针逐个相位平移后，结果必须与非 stream 版逐字节相同 */
static void test_stream_peel(void)
{
    static const size_t kLens[] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 15u, 16u, 17u, 63u,
                                    100u, 1027u, 4096u, 4097u };
    const int nlens = (int)(sizeof(kLens) / sizeof(kLens[0]));
    int li;
    unsigned int off;

    for (li = 0; li < nlens; ++li) {
        const size_t n = kLens[li];
        fill_inputs(n);
        for (off = 0; off < 4u; ++off) {   /* fv_vec2 步长 8 B：0..3 个元素覆盖 32 B 内全部相位 */
            {
                fv_vec2 *r = g_ref + off;
                fv_vec2 *o = g_out + off;
                const size_t bytes = (n + GUARD - off) * sizeof(fv_vec2);

                poison(g_ref, sizeof(g_ref));
                poison(g_out, sizeof(g_out));
                fv_vec2_add_x4(g_a, g_b, r, n);
                fv_vec2_add_x4_stream(g_a, g_b, o, n);
                CHECK(memcmp(r, o, bytes) == 0, "add_x4_stream", n, off);

                poison(g_ref, sizeof(g_ref));
                poison(g_out, sizeof(g_out));
                fv_vec2_sub_x4(g_a, g_b, r, n);
                fv_vec2_sub_x4_stream(g_a, g_b, o, n);
                CHECK(memcmp(r, o, bytes) == 0, "sub_x4_stream", n, off);

                poison(g_ref, sizeof(g_ref));
                poison(g_out, sizeof(g_out));
                fv_vec2_scale_x4(g_a, -3.25f, r, n);
                fv_vec2_scale_x4_stream(g_a, -3.25f, o, n);
                CHECK(memcmp(r, o, bytes) == 0, "scale_x4_stream", n, off);
            }
        }
        for (off = 0; off < 8u; ++off) {   /* float 步长 4 B：0..7 个 float 覆盖 32 B 内全部相位 */
            float *rx = g_rx + off;
            float *ry = g_ry + off;
            float *sx = g_sx + off;
            float *sy = g_sy + off;
            const size_t bytes = (n + GUARD - off) * sizeof(float);

            poison(g_rx, sizeof(g_rx));
            poison(g_ry, sizeof(g_ry));
            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            fv_vec2_add_soa(g_ax, g_ay, g_bx, g_by, rx, ry, n);
            fv_vec2_add_soa_stream(g_ax, g_ay, g_bx, g_by, sx, sy, n);
            CHECK(memcmp(rx, sx, bytes) == 0 && memcmp(ry, sy, bytes) == 0, "add_soa_stream", n,
                  off);

            poison(g_rx, sizeof(g_rx));
            poison(g_ry, sizeof(g_ry));
            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            fv_vec2_scale_soa(g_ax, g_ay, 1.0e-6f, rx, ry, n);
            fv_vec2_scale_soa_stream(g_ax, g_ay, 1.0e-6f, sx, sy, n);
            CHECK(memcmp(rx, sx, bytes) == 0 && memcmp(ry, sy, bytes) == 0, "scale_soa_stream", n,
                  off);

            poison(g_rx, sizeof(g_rx));
            poison(g_ry, sizeof(g_ry));
            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            fv_vec2_normalize_soa(g_ax, g_ay, rx, ry, n);
            fv_vec2_normalize_soa_stream(g_ax, g_ay, sx, sy, n);
            CHECK(memcmp(rx, sx, bytes) == 0 && memcmp(ry, sy, bytes) == 0,
                  "normalize_soa_stream", n, off);
        }
    }

    /* normalize 的 ox / oy 相位不匹配（一个能对齐、另一个永远差 8 B）：
     * 上面那种 peel 推不到「两者同时对齐」，内核必须退回普通 storeu —— 结果仍要逐位相同。
     * 相位组合：ox + 0 个 float、oy + 2 个 float（差 8 B，无论 i 怎么走都相差 8 B）。 */
    {
        static const size_t kLens2[] = { 4u, 8u, 9u, 100u, 4097u };
        int li2;
        for (li2 = 0; li2 < (int)(sizeof(kLens2) / sizeof(kLens2[0])); ++li2) {
            const size_t n2 = kLens2[li2];
            float *rx = g_rx;
            float *ry = g_ry + 2;
            float *sx = g_sx;
            float *sy = g_sy + 2;
            const size_t bytes = (n2 + 2u) * sizeof(float);
            fill_inputs(n2);
            poison(g_rx, sizeof(g_rx));
            poison(g_ry, sizeof(g_ry));
            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            fv_vec2_normalize_soa(g_ax, g_ay, rx, ry, n2);
            fv_vec2_normalize_soa_stream(g_ax, g_ay, sx, sy, n2);
            CHECK(memcmp(rx, sx, bytes) == 0 && memcmp(ry, sy, bytes) == 0,
                  "normalize_soa_stream (mismatched phases)", n2, 2u);
        }
    }

    /* 4 B 对齐、8 B 不齐的 fv_vec2 数组：32 B 边界不可达 → storeu 退化路径。
     * 用 float 缓冲区 + 4 B 偏移造出来（对 fv_vec2 的 8 B 步长来说这就是相位 4）。 */
    {
        static float fa[2u * (CAP + GUARD)], fb[2u * (CAP + GUARD)];
        static float fr[2u * (CAP + GUARD)], fo[2u * (CAP + GUARD)];
        const size_t n = 64u;
        size_t i;
        for (i = 0; i < 2u * n; ++i) {
            fa[i] = rnd_float();
            fb[i] = rnd_float();
        }
        {
            const fv_vec2 *a = (const fv_vec2 *)(const void *)(fa + 1);
            const fv_vec2 *b = (const fv_vec2 *)(const void *)(fb + 1);
            fv_vec2 *r = (fv_vec2 *)(void *)(fr + 1);
            fv_vec2 *o = (fv_vec2 *)(void *)(fo + 1);
            poison(fr, sizeof(fr));
            poison(fo, sizeof(fo));
            fv_vec2_add_x4(a, b, r, n);
            fv_vec2_add_x4_stream(a, b, o, n);
            CHECK(memcmp(fr, fo, (n + 1u) * sizeof(fv_vec2)) == 0,
                  "add_x4_stream (4 B-aligned fallback)", n, 1u);
        }
    }

    /* n == 0：契约说「不读不写任何指针」，所以连 NULL 都可以传 */
    fv_vec2_add_x4_stream(NULL, NULL, NULL, 0);
    fv_vec2_sub_x4_stream(NULL, NULL, NULL, 0);
    fv_vec2_scale_x4_stream(NULL, 1.0f, NULL, 0);
    fv_vec2_add_soa_stream(NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec2_scale_soa_stream(NULL, NULL, 1.0f, NULL, NULL, 0);
    fv_vec2_normalize_soa_stream(NULL, NULL, NULL, NULL, 0);
    fv_vec2_add_x4_mask_inplace(NULL, NULL, NULL, 0);
    ++g_checks;
}

/* 2) _mask_inplace：命中才写，未命中的元素逐字节不变 */
static void test_mask_inplace(void)
{
    static const size_t kLens[] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 63u, 100u, 4096u, 4097u };
    const int nlens = (int)(sizeof(kLens) / sizeof(kLens[0]));
    static const unsigned int kSteps[] = { 1u, 2u, 20u, 0u }; /* 100% / 50% / 5% / 0% */
    int li;
    int di;

    for (li = 0; li < nlens; ++li) {
        const size_t n = kLens[li];
        fill_inputs(n);
        for (di = 0; di < 4; ++di) {
            const unsigned int step = kSteps[di];
            size_t i;
            for (i = 0; i < n; ++i) {
                g_mask[i] = (step == 0u) ? 0u
                                         : (unsigned char)((i % (size_t)step == 0u) ? 1u : 0u);
            }
            /* 参考值按标量语义构造 */
            for (i = 0; i < CAP + GUARD; ++i) {
                g_ref[i] = g_a[i];
                g_out[i] = g_a[i];
            }
            for (i = 0; i < n; ++i) {
                if (g_mask[i] != 0u) {
                    g_ref[i] = fv_vec2_add(g_a[i], g_b[i]);
                }
            }
            fv_vec2_add_x4_mask_inplace(g_out, g_b, g_mask, n);
            CHECK(memcmp(g_ref, g_out, (CAP + GUARD) * sizeof(fv_vec2)) == 0,
                  "add_x4_mask_inplace", n, step);
            /* 未命中的元素另外再单查一次「逐字节没动过」：语义要求的是「不写」，
             * 而不是「写回一个恰好相同的值」。 */
            {
                int untouched = 1;
                for (i = 0; i < n; ++i) {
                    if (g_mask[i] == 0u
                        && memcmp(&g_out[i], &g_a[i], sizeof(fv_vec2)) != 0) {
                        untouched = 0;
                    }
                }
                CHECK(untouched, "add_x4_mask_inplace keeps unset lanes", n, step);
            }
        }
    }

    /* mask 全部命中时，结果必须与 fv_vec2_add_x4 逐位相同（即与批量加法同一条路径） */
    {
        const size_t n = 4097u;
        size_t i;
        fill_inputs(n);
        for (i = 0; i < n; ++i) {
            g_mask[i] = 0xFFu;
        }
        for (i = 0; i < CAP + GUARD; ++i) {
            g_out[i] = g_a[i];
        }
        fv_vec2_add_x4_mask_inplace(g_out, g_b, g_mask, n);
        fv_vec2_add_x4(g_a, g_b, g_ref, n);
        CHECK(memcmp(g_ref, g_out, n * sizeof(fv_vec2)) == 0,
              "add_x4_mask_inplace == add_x4 (all hit)", n, 0u);
    }
}

/* 2b) 掩码字节只取「非零但不是 1」的值（0x02 / 0x80 / 0xFE）：
 *       块级预扫如果用「按位与 0x01010101」判零，0x02/0x80/0xFE 都会被误判成没命中，
 *       而契约是「字节非零就算命中」。两个入口都要用标量参考值盯住。 */
static void test_mask_nonzero_bytes(void)
{
    static const size_t kLens[] = { 1u, 3u, 4u, 5u, 8u, 17u, 63u, 4097u };
    static const unsigned char kVals[4] = { 0x02u, 0x80u, 0xFEu, 0x01u };
    const int nlens = (int)(sizeof(kLens) / sizeof(kLens[0]));
    int li;
    int phase;

    for (li = 0; li < nlens; ++li) {
        const size_t n = kLens[li];
        size_t i;
        fill_inputs(n);
        for (phase = 0; phase < 2; ++phase) {
            for (i = 0; i < n; ++i) {
                g_mask[i] = (i % 5u == 4u) ? 0u : kVals[(i + (size_t)phase) % 4u];
            }
            /* _mask：掩码外写回 a[i] */
            poison(g_ref, sizeof(g_ref));
            poison(g_out, sizeof(g_out));
            fv_vec2_add_x4_mask(g_a, g_b, g_out, g_mask, n);
            for (i = 0; i < n; ++i) {
                g_ref[i] = (g_mask[i] != 0u) ? fv_vec2_add(g_a[i], g_b[i]) : g_a[i];
            }
            CHECK(memcmp(g_ref, g_out, n * sizeof(fv_vec2)) == 0,
                  "add_x4_mask (non-1 mask bytes)", n, (unsigned)phase);

            /* _mask_inplace：命中才写 */
            for (i = 0; i < CAP + GUARD; ++i) {
                g_ref[i] = g_a[i];
                g_out[i] = g_a[i];
            }
            for (i = 0; i < n; ++i) {
                if (g_mask[i] != 0u) {
                    g_ref[i] = fv_vec2_add(g_a[i], g_b[i]);
                }
            }
            fv_vec2_add_x4_mask_inplace(g_out, g_b, g_mask, n);
            CHECK(memcmp(g_ref, g_out, (CAP + GUARD) * sizeof(fv_vec2)) == 0,
                  "add_x4_mask_inplace (non-1 mask bytes)", n, (unsigned)phase);
        }
    }
}

/* 3) sfence：写线程用 stream 变体写，读线程在握手之后读，必须全部是新值 */
#define SF_N 8192u
static fv_vec2 g_sf_a[SF_N], g_sf_b[SF_N], g_sf_o[SF_N], g_sf_exp[SF_N];
static volatile int g_sf_seq;
static volatile int g_sf_bad;

static void *sf_writer(void *arg)
{
    int iters = *(const int *)arg;
    int it;
    for (it = 1; it <= iters; ++it) {
        fv_vec2_add_x4_stream(g_sf_a, g_sf_b, g_sf_o, SF_N);
        /* release：读线程看到 g_sf_seq 之后去读 g_sf_o，必须已经看到 sfence 之后的数据 */
        __atomic_store_n(&g_sf_seq, it, __ATOMIC_RELEASE);
        while (__atomic_load_n(&g_sf_seq, __ATOMIC_ACQUIRE) != 0) {
            /* 等读线程确认 */
        }
    }
    return NULL;
}

static void *sf_reader(void *arg)
{
    int iters = *(const int *)arg;
    int it;
    for (it = 1; it <= iters; ++it) {
        size_t i;
        while (__atomic_load_n(&g_sf_seq, __ATOMIC_ACQUIRE) != it) {
            /* 等写线程 */
        }
        for (i = 0; i < SF_N; ++i) {
            if (g_sf_o[i].x != g_sf_exp[i].x || g_sf_o[i].y != g_sf_exp[i].y) {
                g_sf_bad = 1;
                break;
            }
        }
        __atomic_store_n(&g_sf_seq, 0, __ATOMIC_RELEASE);
    }
    return NULL;
}

static void test_sfence_visibility(void)
{
    const int iters = 64;
    pthread_t tw;
    pthread_t tr;
    size_t i;

    for (i = 0; i < SF_N; ++i) {
        g_sf_a[i].x = rnd_float();
        g_sf_a[i].y = rnd_float();
        g_sf_b[i].x = rnd_float();
        g_sf_b[i].y = rnd_float();
        g_sf_o[i].x = -1.0f;
        g_sf_o[i].y = -1.0f;
        g_sf_exp[i].x = g_sf_a[i].x + g_sf_b[i].x;
        g_sf_exp[i].y = g_sf_a[i].y + g_sf_b[i].y;
    }
    g_sf_seq = 0;
    g_sf_bad = 0;
    if (pthread_create(&tw, NULL, sf_writer, (void *)&iters) != 0) {
        printf("SKIP sfence test: pthread_create failed\n");
        return;
    }
    if (pthread_create(&tr, NULL, sf_reader, (void *)&iters) != 0) {
        printf("SKIP sfence test: pthread_create failed\n");
        pthread_join(tw, NULL);
        return;
    }
    pthread_join(tw, NULL);
    pthread_join(tr, NULL);
    ++g_checks;
    if (g_sf_bad != 0) {
        ++g_failures;
        printf("FAIL sfence: 读线程看到了旧值（NT store 之后没有正确 fence）\n");
    }
}

int main(void)
{
    test_stream_peel();
    test_mask_inplace();
    test_mask_nonzero_bytes();
    test_sfence_visibility();
    printf("test_stream: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}