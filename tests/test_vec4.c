/* fv_vec4 标量 API 的位级参考测试（独立可执行，不依赖批量入口 —— 批量入口的位级对拍
 * 在 test_vec34_parity.c 里，那是等头文件插入完成之后才有的东西）。
 *
 * 三层断言：
 *   1) 结构：sizeof(fv_vec4) == 16、四个分量的偏移 0/4/8/12、数组里相邻两个元素相距 16 B
 *      —— 这三条是批量入口 _x2 的物理前提（2 个 16 B 元素 = 32 B = 一个 ymm）；
 *   2) 表达式：每个函数与 include/flashvec.h 契约里写的表达式逐位相等
 *      （本文件用 Makefile 的 -ffp-contract=off 编译，所以本地参考表达式与库内实现同源无收缩）；
 *   3) NaN 落位：min / max / clamp 与 fv_vec2 / fv_vec3 的同类函数**逐 lane 位级等价** ——
 *      这是契约里「NaN 语义必须与 fv_vec2_min/max/clamp 完全一致」的直接证据，
 *      而不是靠人读代码去比 (a<b)?a:b。
 *
 * 语料：14 个位模式（±0、±1、±inf、两个不同 payload 的 qNaN、最小次正规、最小正规、
 * FLT_MAX、0.5、5、-10）逐对组合，再乘 4 个 lane 位置。
 */
#include "flashvec.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

static int g_checks;
static int g_fail;

static unsigned fbits(float f)
{
    unsigned b;
    memcpy(&b, &f, sizeof(b));
    return b;
}

static float ffrom(unsigned b)
{
    float f;
    memcpy(&f, &b, sizeof(f));
    return f;
}

static int same(float a, float b)
{
    return fbits(a) == fbits(b);
}

/* 对 sqrtf 的结果用：NaN 的 payload 由 libm 实现与硬件 sqrtss 的差异决定，不属于本库契约
 * （批量入口用 sqrtps，与库内 sqrtss 同一条硬件路径，所以批量对拍仍是逐位严格比较）。
 * 这里只要求「同为 NaN」或「位级相同」。 */
static int same_or_both_nan(float a, float b)
{
    if (isnan(a) != 0 && isnan(b) != 0) {
        return 1;
    }
    return same(a, b);
}

static int veq(fv_vec4 a, fv_vec4 b)
{
    return same(a.x, b.x) && same(a.y, b.y) && same(a.z, b.z) && same(a.w, b.w);
}

/* 算术（add/sub/scale/lerp/madd/dot）专用：允许「两侧同为 NaN」。
 * 理由：x86 的两个 NaN 相加/相乘返回**第一个操作数**的 payload，而 add/mul 是可交换指令，
 * 编译器可以任意挑一个操作数在前；所以 NaN ⊕ NaN 的 payload 不是任何实现能承诺的东西。
 * min/max/clamp 不受这条影响 —— 它们按比较返回输入的位模式（payload 原样带出），仍用 veq 严格比较。
 * 有限值、±inf、单侧 NaN 的输入全部走严格比较（下面第 3 节的有限语料 + 第 4 节的单侧 NaN）。 */
static int veq_arith(fv_vec4 a, fv_vec4 b)
{
    return (same(a.x, b.x) || (isnan(a.x) != 0 && isnan(b.x) != 0))
        && (same(a.y, b.y) || (isnan(a.y) != 0 && isnan(b.y) != 0))
        && (same(a.z, b.z) || (isnan(a.z) != 0 && isnan(b.z) != 0))
        && (same(a.w, b.w) || (isnan(a.w) != 0 && isnan(b.w) != 0));
}

#define CHECK(cond)                                                              \
    do {                                                                         \
        ++g_checks;                                                               \
        if (!(cond)) {                                                            \
            ++g_fail;                                                             \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);                \
        }                                                                         \
    } while (0)

/* ---------------- lane 访问（用来把 NaN 逐个放进 4 个位置） ---------------- */

static void set_lane(fv_vec4 *v, int lane, float x)
{
    switch (lane) {
    case 0: v->x = x; break;
    case 1: v->y = x; break;
    case 2: v->z = x; break;
    default: v->w = x; break;
    }
}

static float get_lane(fv_vec4 v, int lane)
{
    switch (lane) {
    case 0: return v.x;
    case 1: return v.y;
    case 2: return v.z;
    default: return v.w;
    }
}

/* ---------------- 二维 / 三维的同名语义（对齐证据的参照物） ---------------- */

static float v2_min1(float a, float b)
{
    return fv_vec2_min(fv_vec2_make(a, a), fv_vec2_make(b, b)).x;
}

static float v2_max1(float a, float b)
{
    return fv_vec2_max(fv_vec2_make(a, a), fv_vec2_make(b, b)).x;
}

static float v2_clamp1(float x, float lo, float hi)
{
    return fv_vec2_clamp(fv_vec2_make(x, x), fv_vec2_make(lo, lo), fv_vec2_make(hi, hi)).x;
}

static float v3_min1(float a, float b)
{
    return fv_vec3_min(fv_vec3_make(a, a, a), fv_vec3_make(b, b, b)).x;
}

static float v3_max1(float a, float b)
{
    return fv_vec3_max(fv_vec3_make(a, a, a), fv_vec3_make(b, b, b)).x;
}

static float v3_clamp1(float x, float lo, float hi)
{
    return fv_vec3_clamp(fv_vec3_make(x, x, x), fv_vec3_make(lo, lo, lo),
                         fv_vec3_make(hi, hi, hi)).x;
}

/* ---------------- 参考表达式（与契约为同一份写法） ---------------- */

static fv_vec4 ref_add(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static fv_vec4 ref_sub(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static fv_vec4 ref_neg(fv_vec4 v)
{
    return fv_vec4_make(-v.x, -v.y, -v.z, -v.w);
}

static fv_vec4 ref_scale(fv_vec4 v, float s)
{
    return fv_vec4_make(v.x * s, v.y * s, v.z * s, v.w * s);
}

static fv_vec4 ref_lerp(fv_vec4 a, fv_vec4 b, float t)
{
    return fv_vec4_make(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                        a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t);
}

static fv_vec4 ref_madd(fv_vec4 a, fv_vec4 b, float s)
{
    return fv_vec4_make(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s, a.w + b.w * s);
}

static fv_vec4 ref_min(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y,
                        a.z < b.z ? a.z : b.z, a.w < b.w ? a.w : b.w);
}

static fv_vec4 ref_max(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y,
                        a.z > b.z ? a.z : b.z, a.w > b.w ? a.w : b.w);
}

static float ref_clamp1(float x, float lo, float hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

static fv_vec4 ref_clamp(fv_vec4 v, fv_vec4 lo, fv_vec4 hi)
{
    return fv_vec4_make(ref_clamp1(v.x, lo.x, hi.x), ref_clamp1(v.y, lo.y, hi.y),
                        ref_clamp1(v.z, lo.z, hi.z), ref_clamp1(v.w, lo.w, hi.w));
}

/* dot 的结合顺序写全：((x*x' + y*y') + z*z') + w*w'。 */
static float ref_dot(fv_vec4 a, fv_vec4 b)
{
    return ((a.x * b.x + a.y * b.y) + a.z * b.z) + a.w * b.w;
}

static float ref_len(fv_vec4 v)
{
    return sqrtf(ref_dot(v, v));
}

static float ref_dist(fv_vec4 a, fv_vec4 b)
{
    fv_vec4 d = ref_sub(a, b);
    return sqrtf(ref_dot(d, d));
}

/* ---------------- 语料 ---------------- */

static const unsigned kBits[] = {
    0x00000000u, /* +0.0 */
    0x80000000u, /* -0.0 */
    0x3F800000u, /* 1.0 */
    0xBF800000u, /* -1.0 */
    0x7F800000u, /* +inf */
    0xFF800000u, /* -inf */
    0x7FC12345u, /* qNaN（正，带 payload） */
    0xFFC0BEEFu, /* qNaN（负，另一个 payload） */
    0x00000001u, /* 最小次正规 */
    0x00800000u, /* 最小正规 */
    0x7F7FFFFFu, /* FLT_MAX */
    0x3F000000u, /* 0.5 */
    0x40A00000u, /* 5.0 */
    0xC1200000u  /* -10.0 */
};
#define NBITS ((int)(sizeof(kBits) / sizeof(kBits[0])))

/* 有限值 + ±inf 的语料：这个集合里没有 NaN，所以第 3 节的算术比较是**严格位级**的，
 * 没有任何「同为 NaN 就放过」的余地。 */
static const unsigned kFiniteBits[] = {
    0x00000000u, /* +0.0 */
    0x80000000u, /* -0.0 */
    0x3F800000u, /* 1.0 */
    0xBF800000u, /* -1.0 */
    0x7F800000u, /* +inf */
    0xFF800000u, /* -inf */
    0x00000001u, /* 最小次正规 */
    0x00800000u, /* 最小正规 */
    0x7F7FFFFFu, /* FLT_MAX */
    0x3F000000u, /* 0.5 */
    0x40A00000u, /* 5.0 */
    0xC1200000u  /* -10.0 */
};
#define NFINITE ((int)(sizeof(kFiniteBits) / sizeof(kFiniteBits[0])))

/* clamp 的三元组合用较小的语料（14^3 太大，6^3 × 4 lane 已经覆盖了每种落位）。 */
static const unsigned kClampBits[] = {
    0x7FC12345u, /* qNaN */
    0x00000000u, /* +0.0 */
    0x3F800000u, /* 1.0 */
    0xBF800000u, /* -1.0 */
    0x7F800000u, /* +inf */
    0xFF800000u  /* -inf */
};
#define NCLAMP ((int)(sizeof(kClampBits) / sizeof(kClampBits[0])))

int main(void)
{
    fv_vec4 a;
    fv_vec4 b;
    fv_vec4 r;
    int i;
    int j;
    int lane;

    /* ---------- 1 · 结构：_x2 的物理前提 ---------- */
    {
        fv_vec4 arr[2];
        CHECK(sizeof(fv_vec4) == 16u);
        CHECK(offsetof(fv_vec4, x) == 0u);
        CHECK(offsetof(fv_vec4, y) == 4u);
        CHECK(offsetof(fv_vec4, z) == 8u);
        CHECK(offsetof(fv_vec4, w) == 12u);
        /* 相邻两个元素 = 32 B 连续：_x2 内核靠这一条用一条 256-bit loadu 取两个元素。 */
        CHECK((size_t)((const char *)(const void *)&arr[1] - (const char *)(const void *)&arr[0])
              == 16u);
    }

    /* ---------- 2 · 硬编码期望（防止参考表达式本身写错） ---------- */
    {
        fv_vec4 A = fv_vec4_make(1.0f, 2.0f, 3.0f, 4.0f);
        fv_vec4 B = fv_vec4_make(10.0f, 20.0f, 30.0f, 40.0f);

        CHECK(veq(fv_vec4_add(A, B), fv_vec4_make(11.0f, 22.0f, 33.0f, 44.0f)));
        CHECK(veq(fv_vec4_sub(B, A), fv_vec4_make(9.0f, 18.0f, 27.0f, 36.0f)));
        CHECK(veq(fv_vec4_neg(A), fv_vec4_make(-1.0f, -2.0f, -3.0f, -4.0f)));
        CHECK(veq(fv_vec4_scale(A, 2.0f), fv_vec4_make(2.0f, 4.0f, 6.0f, 8.0f)));
        CHECK(veq(fv_vec4_lerp(A, B, 0.0f), A));
        CHECK(veq(fv_vec4_lerp(A, B, 1.0f), B));
        CHECK(veq(fv_vec4_lerp(A, B, 0.5f), fv_vec4_make(5.5f, 11.0f, 16.5f, 22.0f)));
        CHECK(veq(fv_vec4_madd(A, B, 1.0f), fv_vec4_make(11.0f, 22.0f, 33.0f, 44.0f)));
        CHECK(same(fv_vec4_dot(A, B), 1.0f * 10.0f + 2.0f * 20.0f + 3.0f * 30.0f + 4.0f * 40.0f));
        CHECK(same(fv_vec4_len_sqr(A), 1.0f + 4.0f + 9.0f + 16.0f));
        CHECK(fv_vec4_near(A, A, 0.0f) == 1);
        CHECK(fv_vec4_near(A, B, 0.5f) == 0);
        CHECK(fv_vec4_near(A, B, 36.0f) == 1); /* 最大分量差正好 36：<= eps 为真 */
        CHECK(fv_vec4_near(A, B, 35.9f) == 0); /* 小一点点就失败 —— 逐分量比较，不是整体范数 */
    }

    /* ---------- 3a · 逐分量算术：有限 + ±inf 语料，严格位级（12 × 12 = 144 对） ---------- */
    for (i = 0; i < NFINITE; ++i) {
        for (j = 0; j < NFINITE; ++j) {
            const float xi = ffrom(kFiniteBits[i]);
            const float xj = ffrom(kFiniteBits[j]);
            a = fv_vec4_make(xi, xj, xi, xj);
            b = fv_vec4_make(xj, xi, xj, xi);

            CHECK(veq(fv_vec4_add(a, b), ref_add(a, b)));
            CHECK(veq(fv_vec4_sub(a, b), ref_sub(a, b)));
            CHECK(veq(fv_vec4_neg(a), ref_neg(a)));
            CHECK(veq(fv_vec4_scale(a, 0.5f), ref_scale(a, 0.5f)));
            CHECK(veq(fv_vec4_scale(a, -3.25f), ref_scale(a, -3.25f)));
            CHECK(veq(fv_vec4_lerp(a, b, 0.25f), ref_lerp(a, b, 0.25f)));
            CHECK(veq(fv_vec4_lerp(a, b, -0.0f), ref_lerp(a, b, -0.0f)));
            CHECK(veq(fv_vec4_madd(a, b, 1.0e-6f), ref_madd(a, b, 1.0e-6f)));
            CHECK(veq(fv_vec4_min(a, b), ref_min(a, b)));
            CHECK(veq(fv_vec4_max(a, b), ref_max(a, b)));
            CHECK(same(fv_vec4_dot(a, b), ref_dot(a, b)));
            CHECK(same(fv_vec4_len_sqr(a), ref_dot(a, a)));
            CHECK(same_or_both_nan(fv_vec4_len(a), ref_len(a)));
        }
    }

    /* ---------- 3b · 含 NaN 的 14 × 14 语料：算术用「同为 NaN 或位级相同」 ---------- */
    for (i = 0; i < NBITS; ++i) {
        for (j = 0; j < NBITS; ++j) {
            const float xi = ffrom(kBits[i]);
            const float xj = ffrom(kBits[j]);
            a = fv_vec4_make(xi, xj, xi, xj);
            b = fv_vec4_make(xj, xi, xj, xi);

            CHECK(veq_arith(fv_vec4_add(a, b), ref_add(a, b)));
            CHECK(veq_arith(fv_vec4_sub(a, b), ref_sub(a, b)));
            CHECK(veq_arith(fv_vec4_neg(a), ref_neg(a)));
            CHECK(veq_arith(fv_vec4_scale(a, 0.5f), ref_scale(a, 0.5f)));
            CHECK(veq_arith(fv_vec4_scale(a, -3.25f), ref_scale(a, -3.25f)));
            CHECK(veq_arith(fv_vec4_lerp(a, b, 0.25f), ref_lerp(a, b, 0.25f)));
            CHECK(veq_arith(fv_vec4_lerp(a, b, -0.0f), ref_lerp(a, b, -0.0f)));
            CHECK(veq_arith(fv_vec4_madd(a, b, 1.0e-6f), ref_madd(a, b, 1.0e-6f)));
            /* min/max 按比较返回输入的位模式（payload 原样带出），所以严格比较 */
            CHECK(veq(fv_vec4_min(a, b), ref_min(a, b)));
            CHECK(veq(fv_vec4_max(a, b), ref_max(a, b)));
        }
    }

    /* ---------- 4 · NaN 落位：与 fv_vec2 / fv_vec3 逐 lane 位级等价 ---------- */
    for (lane = 0; lane < 4; ++lane) {
        for (i = 0; i < NBITS; ++i) {
            for (j = 0; j < NBITS; ++j) {
                const float xi = ffrom(kBits[i]);
                const float xj = ffrom(kBits[j]);
                fv_vec4 base_a = fv_vec4_make(1.0f, -2.0f, 3.0f, -4.0f);
                fv_vec4 base_b = fv_vec4_make(-1.0f, 2.0f, -3.0f, 4.0f);

                set_lane(&base_a, lane, xi);
                set_lane(&base_b, lane, xj);

                r = fv_vec4_min(base_a, base_b);
                CHECK(same(get_lane(r, lane), v2_min1(xi, xj)));
                CHECK(same(get_lane(r, lane), v3_min1(xi, xj)));

                r = fv_vec4_max(base_a, base_b);
                CHECK(same(get_lane(r, lane), v2_max1(xi, xj)));
                CHECK(same(get_lane(r, lane), v3_max1(xi, xj)));
            }
        }
    }

    /* ---------- 5 · clamp：corners 三元组合 + 与 vec2/vec3 逐 lane 位级等价 ---------- */
    for (lane = 0; lane < 4; ++lane) {
        for (i = 0; i < NCLAMP; ++i) {
            for (j = 0; j < NCLAMP; ++j) {
                int k;
                for (k = 0; k < NCLAMP; ++k) {
                    const float xv = ffrom(kClampBits[i]);
                    const float lo = ffrom(kClampBits[j]);
                    const float hi = ffrom(kClampBits[k]);
                    fv_vec4 v = fv_vec4_make(1.0f, -2.0f, 3.0f, -4.0f);
                    fv_vec4 lov = fv_vec4_make(-1.0f, -2.0f, -3.0f, -4.0f);
                    fv_vec4 hiv = fv_vec4_make(1.0f, 2.0f, 3.0f, 4.0f);

                    set_lane(&v, lane, xv);
                    set_lane(&lov, lane, lo);
                    set_lane(&hiv, lane, hi);

                    r = fv_vec4_clamp(v, lov, hiv);
                    CHECK(same(get_lane(r, lane), ref_clamp1(xv, lo, hi)));
                    CHECK(same(get_lane(r, lane), v2_clamp1(xv, lo, hi)));
                    CHECK(same(get_lane(r, lane), v3_clamp1(xv, lo, hi)));
                    CHECK(veq(r, ref_clamp(v, lov, hiv)));
                }
            }
        }
    }
    /* lo > hi 时按同一表达式落位（先比下界）：x < lo 给 lo，否则 x > hi 给 hi。 */
    {
        fv_vec4 v = fv_vec4_make(5.0f, 5.0f, 5.0f, 5.0f);
        fv_vec4 lo = fv_vec4_make(10.0f, 10.0f, 10.0f, 10.0f);
        fv_vec4 hi = fv_vec4_make(0.0f, 0.0f, 0.0f, 0.0f);
        CHECK(veq(fv_vec4_clamp(v, lo, hi), fv_vec4_make(10.0f, 10.0f, 10.0f, 10.0f)));
    }

    /* ---------- 6 · dot / len / dist（含 NaN、inf、FLT_MAX 的溢出组合） ---------- */
    for (i = 0; i < NBITS; ++i) {
        for (j = 0; j < NBITS; ++j) {
            const float xi = ffrom(kBits[i]);
            const float xj = ffrom(kBits[j]);
            a = fv_vec4_make(xi, xj, xi, xj);
            b = fv_vec4_make(xj, xi, xj, xi);

            CHECK(same(fv_vec4_dot(a, b), ref_dot(a, b)));
            CHECK(same(fv_vec4_len_sqr(a), ref_dot(a, a)));
            CHECK(same_or_both_nan(fv_vec4_len(a), ref_len(a)));
            CHECK(same(fv_vec4_dist_sqr(a, b), ref_dot(ref_sub(a, b), ref_sub(a, b))));
            CHECK(same_or_both_nan(fv_vec4_dist(a, b), ref_dist(a, b)));
        }
    }
    /* 零向量的 len_sqr / len 正好是 +0.0 */
    {
        fv_vec4 zero = fv_vec4_make(0.0f, 0.0f, 0.0f, 0.0f);
        CHECK(same(fv_vec4_len_sqr(zero), 0.0f));
        CHECK(same(fv_vec4_len(zero), 0.0f));
        CHECK(same(fv_vec4_dot(zero, zero), 0.0f));
    }

    /* ---------- 7 · normalize：零向量保护与「不会放大出 NaN」 ---------- */
    {
        fv_vec4 zero = fv_vec4_make(0.0f, 0.0f, 0.0f, 0.0f);
        fv_vec4 tiny = fv_vec4_make(5.0e-7f, -5.0e-7f, 0.0f, 0.0f); /* len < FV_EPS */
        fv_vec4 unit = fv_vec4_make(0.0f, 0.0f, 0.0f, 2.0f);
        fv_vec4 q = fv_vec4_make(1.0f, 1.0f, 1.0f, 1.0f);
        fv_vec4 n;

        CHECK(veq(fv_vec4_normalize(zero), zero));
        CHECK(veq(fv_vec4_normalize(tiny), zero)); /* len <= FV_EPS -> (0,0,0,0) */
        CHECK(veq(fv_vec4_normalize(unit), fv_vec4_make(0.0f, 0.0f, 0.0f, 1.0f)));
        CHECK(same(fv_vec4_len(fv_vec4_normalize(q)), 1.0f)); /* 2.0 / 2.0 */

        /* NaN 输入 -> 走 else 分支（len <= FV_EPS 为假），inv = 1/NaN = NaN，
         * 于是每个分量都是 NaN（连 0 * NaN 也是 NaN）—— 与 fv_vec3_normalize 同一条 IEEE 路径。 */
        n = fv_vec4_normalize(fv_vec4_make(ffrom(0x7FC12345u), 0.0f, 0.0f, 0.0f));
        CHECK(isnan(n.x) != 0);
        CHECK(isnan(n.y) != 0);
        /* inf 输入：len = inf，inv = 1/inf = 0，inf * 0 = NaN（不是零向量、也不是 inf） */
        n = fv_vec4_normalize(fv_vec4_make(ffrom(0x7F800000u), 0.0f, 0.0f, 0.0f));
        CHECK(isnan(n.x) != 0);
    }

    /* ---------- 8 · project：平方量阈值 FV_EPS_SQR 与投影表达式 ---------- */
    {
        fv_vec4 zero = fv_vec4_make(0.0f, 0.0f, 0.0f, 0.0f);
        fv_vec4 v = fv_vec4_make(1.0f, 0.0f, 0.0f, 0.0f);
        fv_vec4 onto = fv_vec4_make(0.0f, 0.0f, 3.0f, 4.0f); /* len_sqr = 25 */
        fv_vec4 p = fv_vec4_project(v, onto);
        fv_vec4 refp = fv_vec4_scale(onto, fv_vec4_dot(v, onto) / fv_vec4_len_sqr(onto));

        CHECK(veq(p, refp));
        CHECK(veq(fv_vec4_project(v, zero), zero)); /* denom == 0 <= FV_EPS_SQR -> 零向量 */
        /* onto = (1e-6, 0, 0, 0)：len_sqr = 1e-12 <= FV_EPS_SQR -> 零向量（平方量阈值） */
        CHECK(veq(fv_vec4_project(v, fv_vec4_make(1.0e-6f, 0.0f, 0.0f, 0.0f)), zero));
        /* onto = (2e-6, 0, 0, 0)：len_sqr = 4e-12 > FV_EPS_SQR -> 正常投影（不是零向量） */
        {
            fv_vec4 onto2 = fv_vec4_make(2.0e-6f, 0.0f, 0.0f, 0.0f);
            fv_vec4 got = fv_vec4_project(v, onto2);
            fv_vec4 want = fv_vec4_scale(onto2, fv_vec4_dot(v, onto2) / fv_vec4_len_sqr(onto2));

            CHECK(veq(got, want));
            CHECK(got.x > 0.99f && got.x < 1.01f); /* 2e-6 * (2e-6 / 4e-12) ≈ 1 */
        }
    }

    /* ---------- 9 · near：逐分量 |a-b| <= eps ---------- */
    {
        fv_vec4 p = fv_vec4_make(1.0f, 2.0f, 3.0f, 4.0f);
        fv_vec4 q = fv_vec4_make(1.0f, 2.0f, 3.0f, 4.0000005f);
        fv_vec4 q2 = fv_vec4_make(1.0f, 2.0f, 3.0f, 100.0f);

        CHECK(fv_vec4_near(p, p, 0.0f) == 1);
        CHECK(fv_vec4_near(p, q, 1.0e-6f) == 1);
        CHECK(fv_vec4_near(p, q2, 1.0e-6f) == 0);
        CHECK(fv_vec4_near(p, q2, 96.0f) == 1);
        CHECK(fv_vec4_near(p, q2, 95.9f) == 0);
        /* eps < 0：|a-b| >= 0 > eps（除非差值本身为 NaN），因此只有 NaN 差值的分量算「不大于」 */
        CHECK(fv_vec4_near(p, p, -1.0f) == 0);
        /* NaN 分量：|NaN - NaN| = NaN，NaN > eps 为假 -> 该分量不判失败
         * （与 fv_vec2_near 同一条比较路径；两侧都是 NaN 所以其余分量也相等） */
        {
            fv_vec4 nanv = fv_vec4_make(ffrom(0x7FC12345u), 0.0f, 0.0f, 0.0f);
            CHECK(fv_vec4_near(nanv, nanv, 0.0f) == 1);
        }
    }

    printf("test_vec4: %d checks, %d failures\n", g_checks, g_fail);
    return g_fail == 0 ? 0 : 1;
}
