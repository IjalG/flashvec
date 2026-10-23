/* P1 批量入口的测试：位级一致性（memcmp）、掩码语义、索引语义、尾长、n == 0。
 *
 * 「位级一致」怎么验：参考值由**对应的标量 API** 逐个元素算出
 * （fv_vec2_add / sub / scale / lerp / madd / len / normalize），
 * 再用 memcmp 比较整个输出缓冲区的字节——任何 1 ulp 的差异都会失败。
 * 输出缓冲区在每次调用前被 memset 成 0x5A，所以「该写没写」同样会失败；
 * 另外检查紧接着 n 之后那个元素是否仍是 0x5A（越界写会失败）。
 *
 * 运行：make test-batch
 */
#include "flashvec.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);       \
        }                                                               \
    } while (0)

#define CAP 1024
#define POISON_BITS 0x5A5A5A5Au

/* 尾长覆盖：0-9 连续，再加 4 / 8 / 16 的整倍数与 ±1 */
static const size_t kLens[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33 };
#define NLENS ((int)(sizeof(kLens) / sizeof(kLens[0])))

/* 标量参数（lerp 的 t、scale/madd 的 s）：覆盖 0、-0.0、±1、小数、极小值 */
static const float kParams[] = { 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, 1.0e-6f, -3.25f };
#define NPARAMS ((int)(sizeof(kParams) / sizeof(kParams[0])))

static fv_vec2 g_a[CAP], g_b[CAP], g_o[CAP], g_ref[CAP];
static float g_ax[CAP], g_ay[CAP], g_bx[CAP], g_by[CAP], g_ox[CAP], g_oy[CAP];

static int bits_eq(float a, float b)
{
    uint32_t x, y;
    memcpy(&x, &a, 4);
    memcpy(&y, &b, 4);
    return x == y;
}

static int bytes_eq(const void *p, const void *q, size_t n)
{
    return memcmp(p, q, n) == 0;
}

static void poison(void *p, size_t n)
{
    memset(p, 0x5A, n);
}

static int still_poison(float v)
{
    uint32_t b;
    memcpy(&b, &v, 4);
    return b == POISON_BITS;
}

static uint32_t g_rng = 0x9E3779B9u;

static uint32_t rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/* 定点值：零向量、-0.0、负数、极小值（含次正规）、极大值、1e-6 附近 */
static const float kInteresting[] = {
    0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f, 1.0e-6f, 1.0e-12f,
    1.0e-30f, -1.0e-30f, 4.9e-40f, -4.9e-40f, 3.0e38f, -3.0e38f,
    2.0f, -0.125f, 1024.0f, -7.75f, 1.1754944e-38f, -2.5e-7f
};
#define NINTEREST ((int)(sizeof(kInteresting) / sizeof(kInteresting[0])))

/* 有限随机值：指数限制在 [1,254] 避开 NaN/Inf，符号与尾数随机 */
static float rnd_float(void)
{
    const uint32_t b = rnd_u32();
    const uint32_t bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    float f;
    memcpy(&f, &bits, 4);
    return f;
}

/* 中等量级的随机值（给 dot 用，避免累加溢出）：[-10, 10]，步长 1e-3 */
static float rnd_moderate(void)
{
    return (float)((int32_t)(rnd_u32() % 20001u) - 10000) / 1000.0f;
}

static void fill_vec2(fv_vec2 *v, size_t n, uint32_t seed)
{
    size_t i;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        v[i].x = (i % 3u == 0u) ? kInteresting[i % (size_t)NINTEREST] : rnd_float();
        v[i].y = (i % 3u == 1u) ? kInteresting[(i + 7u) % (size_t)NINTEREST] : rnd_float();
    }
}

static void to_soa(const fv_vec2 *v, float *x, float *y, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        x[i] = v[i].x;
        y[i] = v[i].y;
    }
}

static int soa_bits_eq(const float *ox, const float *oy, const fv_vec2 *ref, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        if (!bits_eq(ox[i], ref[i].x) || !bits_eq(oy[i], ref[i].y)) {
            return 0;
        }
    }
    return 1;
}

/* len_soa 只有一个输出数组，只比 x 一路。 */
static int soa_x_bits_eq(const float *ox, const fv_vec2 *ref, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        if (!bits_eq(ox[i], ref[i].x)) {
            return 0;
        }
    }
    return 1;
}

static int have_avx2(void)
{
#if defined(__x86_64__) || defined(__i386__)
    return __builtin_cpu_supports("avx2") != 0;
#else
    return 0;
#endif
}

int main(void)
{
    int li, pi;
    size_t i;

    /* ---------------- 0 · 比较工具本身够不够敏感（对照） ----------------
     * 如果比较器发现不了 1 ulp，下面所有「逐位一致」都是空话。 */
    {
        float one = 1.0f;
        float next = nextafterf(1.0f, 2.0f);
        CHECK(!bits_eq(one, next));                 /* 相邻可表示数：必须报不同 */
        CHECK(!bits_eq(0.0f, -0.0f));               /* 符号位也是位模式的一部分 */
        CHECK(!bytes_eq(&one, &next, 4));
        CHECK(bits_eq(one, one));
    }

    /* ---------------- 1 · _x4：五个逐元素运算的逐位一致 ---------------- */
    for (li = 0; li < NLENS; ++li) {
        const size_t n = kLens[li];
        fill_vec2(g_a, n, 0x100000u + (uint32_t)li);
        fill_vec2(g_b, n, 0x200000u + (uint32_t)li);

        for (pi = 0; pi < NPARAMS; ++pi) {
            const float s = kParams[pi];

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_add(g_a[i], g_b[i]);
            poison(g_o, sizeof(g_o));
            fv_vec2_add_x4(g_a, g_b, g_o, n);
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_sub(g_a[i], g_b[i]);
            poison(g_o, sizeof(g_o));
            fv_vec2_sub_x4(g_a, g_b, g_o, n);
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_scale(g_a[i], s);
            poison(g_o, sizeof(g_o));
            fv_vec2_scale_x4(g_a, s, g_o, n);
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_lerp(g_a[i], g_b[i], s);
            poison(g_o, sizeof(g_o));
            fv_vec2_lerp_x4(g_a, g_b, s, g_o, n);
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_madd(g_a[i], g_b[i], s);
            poison(g_o, sizeof(g_o));
            fv_vec2_madd_x4(g_a, g_b, s, g_o, n);
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));
        }

        /* 未对齐输入输出（偏移 1 个元素）：证明用的是 loadu/storeu，
         * 同时证明 n 之前和 n 之后的元素都没被动过。 */
        {
            fill_vec2(g_a, n + 1, 0x300000u + (uint32_t)li);
            fill_vec2(g_b, n + 1, 0x400000u + (uint32_t)li);
            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_add(g_a[i + 1], g_b[i + 1]);
            poison(g_o, sizeof(g_o));
            fv_vec2_add_x4(g_a + 1, g_b + 1, g_o + 1, n);
            CHECK(bytes_eq(g_o + 1, g_ref, n * sizeof(fv_vec2)));
            CHECK(still_poison(g_o[0].x) && still_poison(g_o[0].y));
            if (n + 1 < CAP) {
                CHECK(still_poison(g_o[n + 1].x) && still_poison(g_o[n + 1].y));
            }
        }
    }

    /* ---------------- 2 · _x4：越界写检查（大数组） ---------------- */
    {
        const size_t n = 1000;
        fill_vec2(g_a, n, 0x555501u);
        fill_vec2(g_b, n, 0x555502u);
        for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_add(g_a[i], g_b[i]);
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4(g_a, g_b, g_o, n);
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));
        CHECK(still_poison(g_o[n].x));
        CHECK(still_poison(g_o[n].y));
        CHECK(still_poison(g_o[n + 1].x));
    }

    /* ---------------- 3 · 掩码变体：语义 + 位级一致 ---------------- */
    for (li = 0; li < NLENS; ++li) {
        const size_t n = kLens[li];
        uint8_t mask[CAP];
        int masked_lanes_ok = 1;

        fill_vec2(g_a, n, 0x600000u + (uint32_t)li);
        fill_vec2(g_b, n, 0x700000u + (uint32_t)li);
        /* 0 / 非 0 交替（跨 8 元素主循环边界），保证掩码里两种值都有 */
        for (i = 0; i < n; ++i) {
            mask[i] = (uint8_t)((i % 3u == 0u) ? 0u : (uint32_t)(i % 7u + 1u));
        }
        for (i = 0; i < n; ++i) {
            g_ref[i] = (mask[i] != 0u) ? fv_vec2_add(g_a[i], g_b[i]) : g_a[i];
        }
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4_mask(g_a, g_b, g_o, mask, n);
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

        /* 掩码外必须逐位等于 a[i] 的原值（不是「不写」，也不是写 o 的旧值） */
        for (i = 0; i < n; ++i) {
            if (mask[i] == 0u &&
                (!bits_eq(g_o[i].x, g_a[i].x) || !bits_eq(g_o[i].y, g_a[i].y))) {
                masked_lanes_ok = 0;
            }
        }
        CHECK(masked_lanes_ok);
        CHECK(still_poison(g_o[n].x));
    }

    /* 全 0 掩码：输出必须逐位等于 a；全 1 掩码：必须等于逐元素加法 */
    {
        const size_t n = 17;
        uint8_t mask[CAP];

        fill_vec2(g_a, n, 0x800001u);
        fill_vec2(g_b, n, 0x800002u);

        memset(mask, 0, sizeof(mask));
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4_mask(g_a, g_b, g_o, mask, n);
        CHECK(bytes_eq(g_o, g_a, n * sizeof(fv_vec2)));

        memset(mask, 1, sizeof(mask));
        for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_add(g_a[i], g_b[i]);
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4_mask(g_a, g_b, g_o, mask, n);
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));
    }

    /* ---------------- 4 · 索引变体：重复 / 乱序 / m == 0 ---------------- */
    {
        static const uint32_t dup[] = { 5u, 1u, 1u, 0u, 15u, 3u, 3u, 3u, 8u, 2u, 7u, 1u };
        static const uint32_t perm[] = { 3u, 0u, 7u, 2u, 6u, 1u, 5u, 4u };
        static const uint32_t one[] = { 6u };
        const size_t n = 16;
        size_t k;

        fill_vec2(g_a, n, 0x900001u);
        fill_vec2(g_b, n, 0x900002u);

        /* 重复 + 乱序 */
        memcpy(g_o, g_a, n * sizeof(fv_vec2));
        memcpy(g_ref, g_a, n * sizeof(fv_vec2));
        for (k = 0; k < sizeof(dup) / sizeof(dup[0]); ++k) {
            const uint32_t id = dup[k];
            g_ref[id] = fv_vec2_add(g_a[id], g_b[id]);
        }
        fv_vec2_add_x4_idx(g_a, g_b, g_o, dup, sizeof(dup) / sizeof(dup[0]));
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

        /* 乱序（双重写入同一个下标也要得到同一个值） */
        memcpy(g_o, g_a, n * sizeof(fv_vec2));
        memcpy(g_ref, g_a, n * sizeof(fv_vec2));
        for (k = 0; k < sizeof(perm) / sizeof(perm[0]); ++k) {
            const uint32_t id = perm[k];
            g_ref[id] = fv_vec2_add(g_a[id], g_b[id]);
        }
        fv_vec2_add_x4_idx(g_a, g_b, g_o, perm, sizeof(perm) / sizeof(perm[0]));
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

        /* m == 0：什么都不写 */
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4_idx(g_a, g_b, g_o, dup, 0);
        CHECK(still_poison(g_o[0].x) && still_poison(g_o[0].y) && still_poison(g_o[15].y));

        /* 单个索引：只有那一个元素变 */
        memcpy(g_o, g_a, n * sizeof(fv_vec2));
        memcpy(g_ref, g_a, n * sizeof(fv_vec2));
        g_ref[6] = fv_vec2_add(g_a[6], g_b[6]);
        fv_vec2_add_x4_idx(g_a, g_b, g_o, one, 1);
        CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));

        /* 重复索引写到同一个下标：结果与顺序无关（值相同） */
        {
            static const uint32_t same[] = { 4u, 4u, 4u, 4u };
            memcpy(g_o, g_a, n * sizeof(fv_vec2));
            memcpy(g_ref, g_a, n * sizeof(fv_vec2));
            g_ref[4] = fv_vec2_add(g_a[4], g_b[4]);
            fv_vec2_add_x4_idx(g_a, g_b, g_o, same, sizeof(same) / sizeof(same[0]));
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));
        }

        /* scale_x4_idx：多个 s，含尾部下标 */
        for (pi = 0; pi < NPARAMS; ++pi) {
            const float s = kParams[pi];
            memcpy(g_o, g_a, n * sizeof(fv_vec2));
            memcpy(g_ref, g_a, n * sizeof(fv_vec2));
            for (k = 0; k < sizeof(dup) / sizeof(dup[0]); ++k) {
                const uint32_t id = dup[k];
                g_ref[id] = fv_vec2_scale(g_a[id], s);
            }
            fv_vec2_scale_x4_idx(g_a, s, g_o, dup, sizeof(dup) / sizeof(dup[0]));
            CHECK(bytes_eq(g_o, g_ref, n * sizeof(fv_vec2)));
        }
    }

    /* ---------------- 5 · SoA：add / scale / len / normalize ---------------- */
    for (li = 0; li < NLENS; ++li) {
        const size_t n = kLens[li];
        fill_vec2(g_a, n, 0xA00000u + (uint32_t)li);
        fill_vec2(g_b, n, 0xB00000u + (uint32_t)li);
        to_soa(g_a, g_ax, g_ay, n);
        to_soa(g_b, g_bx, g_by, n);

        for (pi = 0; pi < NPARAMS; ++pi) {
            const float s = kParams[pi];

            for (i = 0; i < n; ++i) {
                g_ref[i] = fv_vec2_add(fv_vec2_make(g_ax[i], g_ay[i]),
                                       fv_vec2_make(g_bx[i], g_by[i]));
            }
            poison(g_ox, sizeof(g_ox));
            poison(g_oy, sizeof(g_oy));
            fv_vec2_add_soa(g_ax, g_ay, g_bx, g_by, g_ox, g_oy, n);
            CHECK(soa_bits_eq(g_ox, g_oy, g_ref, n));
            CHECK(still_poison(g_ox[n]));

            for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_scale(fv_vec2_make(g_ax[i], g_ay[i]), s);
            poison(g_ox, sizeof(g_ox));
            poison(g_oy, sizeof(g_oy));
            fv_vec2_scale_soa(g_ax, g_ay, s, g_ox, g_oy, n);
            CHECK(soa_bits_eq(g_ox, g_oy, g_ref, n));
        }

        /* len_soa：逐位等于 fv_vec2_len（只有一个输出数组） */
        for (i = 0; i < n; ++i) g_ref[i].x = fv_vec2_len(fv_vec2_make(g_ax[i], g_ay[i]));
        poison(g_ox, sizeof(g_ox));
        fv_vec2_len_soa(g_ax, g_ay, g_ox, n);
        CHECK(soa_x_bits_eq(g_ox, g_ref, n));
        CHECK(still_poison(g_ox[n]));

        /* normalize_soa：逐位等于 fv_vec2_normalize（含零向量 -> (0,0)） */
        for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_normalize(fv_vec2_make(g_ax[i], g_ay[i]));
        poison(g_ox, sizeof(g_ox));
        poison(g_oy, sizeof(g_oy));
        fv_vec2_normalize_soa(g_ax, g_ay, g_ox, g_oy, n);
        CHECK(soa_bits_eq(g_ox, g_oy, g_ref, n));
        CHECK(still_poison(g_ox[n]) && still_poison(g_oy[n]));
    }

    /* 零向量集合：len <= FV_EPS 必须输出 +0.0，且不产生 NaN */
    {
        const size_t n = 9;
        /* 最后一列是次正规值（4.9e-40 < FLT_MIN），用来确认零向量保护不依赖指数范围 */
        const float zx[9] = { 0.0f, -0.0f, 1.0e-7f, -1.0e-7f, 1.0e-8f, 0.0f, -1.0e-30f, 4.9e-40f, -4.9e-40f };
        const float zy[9] = { 0.0f, 0.0f, 0.0f, 0.0f, -1.0e-8f, -0.0f, 1.0e-30f, -4.9e-40f, 4.9e-40f };
        int no_nan = 1;

        for (i = 0; i < n; ++i) g_ref[i] = fv_vec2_normalize(fv_vec2_make(zx[i], zy[i]));
        poison(g_ox, sizeof(g_ox));
        poison(g_oy, sizeof(g_oy));
        fv_vec2_normalize_soa(zx, zy, g_ox, g_oy, n);
        CHECK(soa_bits_eq(g_ox, g_oy, g_ref, n));
        for (i = 0; i < n; ++i) {
            if (isnan(g_ox[i]) || isnan(g_oy[i])) no_nan = 0;
            if (!bits_eq(g_ox[i], 0.0f) || !bits_eq(g_oy[i], 0.0f)) no_nan = 0;
        }
        CHECK(no_nan);
    }

    /* ---------------- 6 · dot_soa：n < 8 逐位一致，n >= 8 用相对误差 ---------------- */
    for (li = 0; li < 8; ++li) {
        const size_t n = kLens[li];
        float ref = 0.0f;
        fill_vec2(g_a, n, 0xC00000u + (uint32_t)li);
        fill_vec2(g_b, n, 0xD00000u + (uint32_t)li);
        to_soa(g_a, g_ax, g_ay, n);
        to_soa(g_b, g_bx, g_by, n);
        /* 与实现里 n < 8 的尾巴用同一个顺序表达式 */
        for (i = 0; i < n; ++i) ref = ref + (g_ax[i] * g_bx[i] + g_ay[i] * g_by[i]);
        CHECK(bits_eq(fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, n), ref));
    }
    {
        static const size_t sizes[] = { 8, 9, 16, 17, 33, 100, 1000 };
        const size_t nsizes = sizeof(sizes) / sizeof(sizes[0]);
        double max_rel = 0.0;
        size_t si;
        for (si = 0; si < nsizes; ++si) {
            const size_t n = sizes[si];
            double dref = 0.0;
            double scale = 0.0; /* 各项绝对值之和：抵消时也不会把相对误差放大成假失败 */
            double got;
            double rel;
            g_rng = 0xE0000000u + (uint32_t)si;
            for (i = 0; i < n; ++i) {
                g_ax[i] = rnd_moderate();
                g_ay[i] = rnd_moderate();
                g_bx[i] = rnd_moderate();
                g_by[i] = rnd_moderate();
            }
            for (i = 0; i < n; ++i) {
                dref += (double)g_ax[i] * (double)g_bx[i] + (double)g_ay[i] * (double)g_by[i];
                scale += fabs((double)g_ax[i] * (double)g_bx[i]) + fabs((double)g_ay[i] * (double)g_by[i]);
            }
            got = (double)fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, n);
            rel = fabs(got - dref) / (scale > 0.0 ? scale : 1.0);
            if (rel > max_rel) max_rel = rel;
            CHECK(rel < 1.0e-5);
        }
        printf("test_batch: dot_soa max_rel=%.3e (scale-normalized, n<=1000)\n", max_rel);
    }
    /* 全零输入 -> 精确 +0.0（n = 1000 走 32 元素主循环 + 尾巴） */
    {
        const size_t n = 1000;
        for (i = 0; i <= n; ++i) { g_ax[i] = 0.0f; g_ay[i] = 0.0f; g_bx[i] = 0.0f; g_by[i] = 0.0f; }
        CHECK(bits_eq(fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, n), 0.0f));
        CHECK(bits_eq(fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, 0), 0.0f));
    }

    /* ---------------- 7 · n == 0：不读不写任何指针（可以传 NULL） ---------------- */
    {
        fv_vec2_add_x4(NULL, NULL, NULL, 0);
        fv_vec2_sub_x4(NULL, NULL, NULL, 0);
        fv_vec2_scale_x4(NULL, 0.5f, NULL, 0);
        fv_vec2_lerp_x4(NULL, NULL, 0.5f, NULL, 0);
        fv_vec2_madd_x4(NULL, NULL, 0.5f, NULL, 0);
        fv_vec2_add_x4_mask(NULL, NULL, NULL, NULL, 0);
        fv_vec2_add_x4_idx(NULL, NULL, NULL, NULL, 0);
        fv_vec2_scale_x4_idx(NULL, 0.5f, NULL, NULL, 0);
        fv_vec2_add_soa(NULL, NULL, NULL, NULL, NULL, NULL, 0);
        fv_vec2_scale_soa(NULL, NULL, 0.5f, NULL, NULL, 0);
        fv_vec2_len_soa(NULL, NULL, NULL, 0);
        fv_vec2_normalize_soa(NULL, NULL, NULL, NULL, 0);
        /* 能执行到这里就说明上面 12 个调用都没有解引用空指针 */
        CHECK(bits_eq(fv_vec2_dot_soa(NULL, NULL, NULL, NULL, 0), 0.0f));
    }

    printf("test_batch: avx2=%d\n", have_avx2());
    printf("test_batch: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
