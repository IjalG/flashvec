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

#include <float.h>
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


/* ============ the 13 new _stream entries: alignment x tail length x special values ============
 * Each case compares "stream output" against "non-stream output" over n + 4 elements starting at the shifted output
 * pointer (both buffers are poisoned first), so a wrong NT body, a missed peel element, an out-of-range tail store or an
 * overwritten GUARD all show up as a mismatch.
 * Coverage: n = 0..40 (every tail length of _x4 / _x2 / _soa) + 63/64/65/100/1023/4096/4097/8191;
 * the output pointer shifts by 0..8, 16, 32 and 64 elements -- 0..7 walk through every 32 B phase of the element stride,
 * the larger shifts move the allocation itself. Inputs cycle through -0.0f, +-inf, NaN, FLT_MIN/FLT_MAX, denormals,
 * exact zero vectors (the normalize zero-vector guard) and near-zero vectors. */
#define N_CAP 8300u
#define N_GUARD 8u

static fv_vec2 n_v2a[N_CAP + N_GUARD], n_v2b[N_CAP + N_GUARD];
static fv_vec2 n_v2r[N_CAP + N_GUARD], n_v2s[N_CAP + N_GUARD];
static fv_vec4 n_v4a[N_CAP + N_GUARD], n_v4b[N_CAP + N_GUARD];
static fv_vec4 n_v4r[N_CAP + N_GUARD], n_v4s[N_CAP + N_GUARD];
static float n_ax[N_CAP + N_GUARD], n_ay[N_CAP + N_GUARD], n_az[N_CAP + N_GUARD];
static float n_aw[N_CAP + N_GUARD], n_bx[N_CAP + N_GUARD], n_by[N_CAP + N_GUARD];
static float n_bz[N_CAP + N_GUARD], n_bw[N_CAP + N_GUARD];
static float n_r1[N_CAP + N_GUARD], n_r2[N_CAP + N_GUARD], n_r3[N_CAP + N_GUARD];
static float n_r4[N_CAP + N_GUARD];
static float n_s1[N_CAP + N_GUARD], n_s2[N_CAP + N_GUARD], n_s3[N_CAP + N_GUARD];
static float n_s4[N_CAP + N_GUARD];
static fv_mat4 n_m;

/* One cycle over the values a vector kernel can mishandle: signed zero, infinities, NaN, the float range extremes,
 * denormals, exact zeros and near-zeros. */
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

static void fill_edge(size_t n)
{
    size_t i;
    for (i = 0; i < n + N_GUARD; ++i) {
        const float ax = edge_value(i);
        const float ay = (i % 5u == 0u) ? 0.0f : edge_value(i * 3u + 1u);      /* zero vector every 5 */
        const float az = (i % 7u == 0u) ? 1.0e-30f : edge_value(i * 5u + 2u); /* near-zero every 7 */
        const float aw = edge_value(i * 7u + 3u);
        n_v2a[i].x = ax;
        n_v2a[i].y = ay;
        n_v2b[i].x = edge_value(i * 11u + 4u);
        n_v2b[i].y = edge_value(i * 13u + 5u);
        n_v4a[i].x = ax;
        n_v4a[i].y = ay;
        n_v4a[i].z = az;
        n_v4a[i].w = aw;
        n_v4b[i].x = edge_value(i * 11u + 4u);
        n_v4b[i].y = edge_value(i * 13u + 5u);
        n_v4b[i].z = edge_value(i * 17u + 6u);
        n_v4b[i].w = edge_value(i * 19u + 7u);
        n_ax[i] = ax;
        n_ay[i] = ay;
        n_az[i] = az;
        n_aw[i] = aw;
        n_bx[i] = n_v4b[i].x;
        n_by[i] = n_v4b[i].y;
        n_bz[i] = n_v4b[i].z;
        n_bw[i] = n_v4b[i].w;
    }
}

static void set_matrix(void)
{
    /* Finite coefficients on purpose: a NaN in the matrix would turn every output into NaN and hide every other difference. */
    static const float kM[12] = { 1.5f, -0.5f, 0.25f, 7.0f,
                                  -2.0f, 3.25f, 0.0f, -1.0f,
                                  1.0e-30f, -0.0f, 0.5f, 2.0f };
    size_t i;
    for (i = 0; i < 12u; ++i) {
        n_m.m[i] = kM[i];
    }
    n_m.m[12] = 0.0f;
    n_m.m[13] = 0.0f;
    n_m.m[14] = 0.0f;
    n_m.m[15] = 1.0f;
}

static void test_new_stream_entries(void)
{
    static const size_t kBigLens[] = { 63u, 64u, 65u, 100u, 1023u, 4096u, 4097u, 8191u };
    static const size_t kOffs[] = { 0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 16u, 32u, 64u };
    static const float kS[6] = { 0.0f, -0.0f, 1.0f, -3.25f, 1.0e-30f, 0.5f };
    static const float kT[8] = { 0.0f, -0.0f, 1.0f, 0.5f, -1.0f, 1.0e-30f,
                                 (float)INFINITY, (float)NAN };
    size_t lens[64];
    int nlens = 0;
    int li;
    int oi;
    size_t i;

    for (i = 0; i <= 40u; ++i) {
        lens[nlens++] = i;
    }
    for (i = 0; i < sizeof(kBigLens) / sizeof(kBigLens[0]); ++i) {
        lens[nlens++] = kBigLens[i];
    }
    set_matrix();

    for (li = 0; li < nlens; ++li) {
        const size_t n = lens[li];
        fill_edge(n);
        for (oi = 0; oi < (int)(sizeof(kOffs) / sizeof(kOffs[0])); ++oi) {
            const size_t off = kOffs[oi];
            const size_t nbytes2 = (n + 4u) * sizeof(fv_vec2);
            const size_t nbytes4 = (n + 4u) * sizeof(fv_vec4);
            const size_t nbytesf = (n + 4u) * sizeof(float);
            fv_vec2 *const v2r = n_v2r + off;
            fv_vec2 *const v2s = n_v2s + off;
            fv_vec4 *const v4r = n_v4r + off;
            fv_vec4 *const v4s = n_v4s + off;
            float *const f1r = n_r1 + off;
            float *const f2r = n_r2 + off;
            float *const f3r = n_r3 + off;
            float *const f4r = n_r4 + off;
            float *const f1s = n_s1 + off;
            float *const f2s = n_s2 + off;
            float *const f3s = n_s3 + off;
            float *const f4s = n_s4 + off;
            const float s = kS[(size_t)(li + oi) % 6u];
            const float t = kT[(size_t)(li + oi) % 8u];

            /* ---- vec2 _x4: lerp / madd ---- */
            poison(n_v2r, sizeof(n_v2r));
            poison(n_v2s, sizeof(n_v2s));
            fv_vec2_lerp_x4(n_v2a, n_v2b, t, v2r, n);
            fv_vec2_lerp_x4_stream(n_v2a, n_v2b, t, v2s, n);
            CHECK(memcmp(v2r, v2s, nbytes2) == 0, "lerp_x4_stream", n, off);

            poison(n_v2r, sizeof(n_v2r));
            poison(n_v2s, sizeof(n_v2s));
            fv_vec2_madd_x4(n_v2a, n_v2b, s, v2r, n);
            fv_vec2_madd_x4_stream(n_v2a, n_v2b, s, v2s, n);
            CHECK(memcmp(v2r, v2s, nbytes2) == 0, "madd_x4_stream", n, off);

            /* ---- vec3 SoA: add / scale / normalize ---- */
            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            fv_vec3_add_soa(n_ax, n_ay, n_az, n_bx, n_by, n_bz, f1r, f2r, f3r, n);
            fv_vec3_add_soa_stream(n_ax, n_ay, n_az, n_bx, n_by, n_bz, f1s, f2s, f3s, n);
            CHECK(memcmp(f1r, f1s, nbytesf) == 0 && memcmp(f2r, f2s, nbytesf) == 0
                      && memcmp(f3r, f3s, nbytesf) == 0,
                  "add_soa_stream (vec3)", n, off);

            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            fv_vec3_scale_soa(n_ax, n_ay, n_az, s, f1r, f2r, f3r, n);
            fv_vec3_scale_soa_stream(n_ax, n_ay, n_az, s, f1s, f2s, f3s, n);
            CHECK(memcmp(f1r, f1s, nbytesf) == 0 && memcmp(f2r, f2s, nbytesf) == 0
                      && memcmp(f3r, f3s, nbytesf) == 0,
                  "scale_soa_stream (vec3)", n, off);

            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            fv_vec3_normalize_soa(n_ax, n_ay, n_az, f1r, f2r, f3r, n);
            fv_vec3_normalize_soa_stream(n_ax, n_ay, n_az, f1s, f2s, f3s, n);
            CHECK(memcmp(f1r, f1s, nbytesf) == 0 && memcmp(f2r, f2s, nbytesf) == 0
                      && memcmp(f3r, f3s, nbytesf) == 0,
                  "normalize_soa_stream (vec3)", n, off);

            /* ---- vec4 _x2: add / sub / scale / lerp / madd ---- */
            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_vec4_add_x2(n_v4a, n_v4b, v4r, n);
            fv_vec4_add_x2_stream(n_v4a, n_v4b, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "add_x2_stream", n, off);

            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_vec4_sub_x2(n_v4a, n_v4b, v4r, n);
            fv_vec4_sub_x2_stream(n_v4a, n_v4b, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "sub_x2_stream", n, off);

            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_vec4_scale_x2(n_v4a, s, v4r, n);
            fv_vec4_scale_x2_stream(n_v4a, s, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "scale_x2_stream", n, off);

            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_vec4_lerp_x2(n_v4a, n_v4b, t, v4r, n);
            fv_vec4_lerp_x2_stream(n_v4a, n_v4b, t, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "lerp_x2_stream", n, off);

            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_vec4_madd_x2(n_v4a, n_v4b, s, v4r, n);
            fv_vec4_madd_x2_stream(n_v4a, n_v4b, s, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "madd_x2_stream", n, off);

            /* ---- vec4 SoA: normalize ---- */
            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_r4, sizeof(n_r4));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            poison(n_s4, sizeof(n_s4));
            fv_vec4_normalize_soa(n_ax, n_ay, n_az, n_aw, f1r, f2r, f3r, f4r, n);
            fv_vec4_normalize_soa_stream(n_ax, n_ay, n_az, n_aw, f1s, f2s, f3s, f4s, n);
            CHECK(memcmp(f1r, f1s, nbytesf) == 0 && memcmp(f2r, f2s, nbytesf) == 0
                      && memcmp(f3r, f3s, nbytesf) == 0 && memcmp(f4r, f4s, nbytesf) == 0,
                  "normalize_soa_stream (vec4)", n, off);

            /* ---- mat4 transform: points / dirs ---- */
            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_mat4_transform_points(&n_m, n_v4a, v4r, n);
            fv_mat4_transform_points_stream(&n_m, n_v4a, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "mat4_transform_points_stream", n, off);

            poison(n_v4r, sizeof(n_v4r));
            poison(n_v4s, sizeof(n_v4s));
            fv_mat4_transform_dirs(&n_m, n_v4a, v4r, n);
            fv_mat4_transform_dirs_stream(&n_m, n_v4a, v4s, n);
            CHECK(memcmp(v4r, v4s, nbytes4) == 0, "mat4_transform_dirs_stream", n, off);
        }
    }
}

/* Phase combinations the peel cannot fix: the 32 B boundary stays unreachable, so the kernel must fall back to normal
 * stores for the whole span (bit-identical result, no #GP). Three shapes:
 *   - fv_vec4 at a 4 B offset (not even 16 B aligned: a 16 B stride cannot remove a 4 B phase);
 *   - fv_vec2 at a 4 B offset (the same shape tests/test_stream.c already covers for add/sub/scale);
 *   - three / four SoA outputs whose phases differ by 8 B, so they can never all reach the boundary. */
static void test_new_stream_fallback_phases(void)
{
    static float fa[2u * 520u], fb[2u * 520u], fr[2u * 520u], fs[2u * 520u];
    static const size_t kFO[3] = { 1u, 2u, 3u };  /* fv_vec4 at a 4 / 8 / 12 B offset: not 16 B aligned */
    const size_t n = 200u;
    size_t i;
    int fi;

    for (i = 0; i < 2u * n; ++i) {
        fa[i] = edge_value(i);
        fb[i] = edge_value(i + 3u);
    }
    set_matrix();

    for (fi = 0; fi < 3; ++fi) {
        const size_t fo = kFO[fi];
        const fv_vec4 *a = (const fv_vec4 *)(const void *)(fa + fo);
        const fv_vec4 *b = (const fv_vec4 *)(const void *)(fb + fo);
        fv_vec4 *r = (fv_vec4 *)(void *)(fr + fo);
        fv_vec4 *o = (fv_vec4 *)(void *)(fs + fo);
        const size_t bytes = (n + fo) * sizeof(fv_vec4);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec4_add_x2(a, b, r, n);
        fv_vec4_add_x2_stream(a, b, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "add_x2_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec4_sub_x2(a, b, r, n);
        fv_vec4_sub_x2_stream(a, b, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "sub_x2_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec4_scale_x2(a, -3.25f, r, n);
        fv_vec4_scale_x2_stream(a, -3.25f, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "scale_x2_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec4_lerp_x2(a, b, 0.5f, r, n);
        fv_vec4_lerp_x2_stream(a, b, 0.5f, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "lerp_x2_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec4_madd_x2(a, b, 1.0e-30f, r, n);
        fv_vec4_madd_x2_stream(a, b, 1.0e-30f, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "madd_x2_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_mat4_transform_points(&n_m, a, r, n);
        fv_mat4_transform_points_stream(&n_m, a, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "mat4_transform_points_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_mat4_transform_dirs(&n_m, a, r, n);
        fv_mat4_transform_dirs_stream(&n_m, a, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "mat4_transform_dirs_stream (misaligned fallback)", n, fo);
    }
    {
        /* fv_vec2 at a 4 B offset (not 8 B aligned): the 32 B boundary is unreachable for an 8 B stride. */
        const size_t fo = 1u;
        const fv_vec2 *a = (const fv_vec2 *)(const void *)(fa + fo);
        const fv_vec2 *b = (const fv_vec2 *)(const void *)(fb + fo);
        fv_vec2 *r = (fv_vec2 *)(void *)(fr + fo);
        fv_vec2 *o = (fv_vec2 *)(void *)(fs + fo);
        const size_t bytes = (n + fo) * sizeof(fv_vec2);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec2_lerp_x4(a, b, 0.25f, r, n);
        fv_vec2_lerp_x4_stream(a, b, 0.25f, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "lerp_x4_stream (misaligned fallback)", n, fo);

        poison(fr, sizeof(fr));
        poison(fs, sizeof(fs));
        fv_vec2_madd_x4(a, b, 2.0f, r, n);
        fv_vec2_madd_x4_stream(a, b, 2.0f, o, n);
        CHECK(memcmp(fr, fs, bytes) == 0, "madd_x4_stream (misaligned fallback)", n, fo);
    }
    {
        /* ox aligned, oy + 8 B, oz + 16 B (vec3) and additionally ow + 24 B (vec4): no i makes them align together. */
        static const size_t kLens[] = { 1u, 7u, 8u, 9u, 100u, 4097u };
        int li;
        for (li = 0; li < (int)(sizeof(kLens) / sizeof(kLens[0])); ++li) {
            const size_t n2 = kLens[li];
            const size_t bytes = (n2 + 1u) * sizeof(float);
            fill_edge(n2);

            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            fv_vec3_normalize_soa(n_ax, n_ay, n_az, n_r1, n_r2 + 4, n_r3 + 8, n2);
            fv_vec3_normalize_soa_stream(n_ax, n_ay, n_az, n_s1, n_s2 + 4, n_s3 + 8, n2);
            CHECK(memcmp(n_r1, n_s1, bytes) == 0 && memcmp(n_r2, n_s2, bytes) == 0
                      && memcmp(n_r3, n_s3, bytes) == 0,
                  "normalize_soa_stream (vec3, mismatched phases)", n2, 4u);

            poison(n_r1, sizeof(n_r1));
            poison(n_r2, sizeof(n_r2));
            poison(n_r3, sizeof(n_r3));
            poison(n_r4, sizeof(n_r4));
            poison(n_s1, sizeof(n_s1));
            poison(n_s2, sizeof(n_s2));
            poison(n_s3, sizeof(n_s3));
            poison(n_s4, sizeof(n_s4));
            fv_vec4_normalize_soa(n_ax, n_ay, n_az, n_aw, n_r1, n_r2 + 4, n_r3 + 8, n_r4 + 12, n2);
            fv_vec4_normalize_soa_stream(n_ax, n_ay, n_az, n_aw, n_s1, n_s2 + 4, n_s3 + 8,
                                         n_s4 + 12, n2);
            CHECK(memcmp(n_r1, n_s1, bytes) == 0 && memcmp(n_r2, n_s2, bytes) == 0
                      && memcmp(n_r3, n_s3, bytes) == 0 && memcmp(n_r4, n_s4, bytes) == 0,
                  "normalize_soa_stream (vec4, mismatched phases)", n2, 4u);
        }
    }

    /* n == 0: the contract says no pointer is read or written, so NULL is legal for every one of the 13 new entries. */
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
}

int main(void)
{
    test_stream_peel();
    test_mask_inplace();
    test_mask_nonzero_bytes();
    test_sfence_visibility();
    test_new_stream_entries();
    test_new_stream_fallback_phases();
    printf("test_stream: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}