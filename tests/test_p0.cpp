/* P0 的单元测试：min / max / clamp / madd、AABB merge / overlaps、FV_EPS 量纲修正。
 * 风格与 tests/test_flashvec.cpp 一致（纯 C 风格断言，不依赖测试框架）；
 * 既有测试文件的 30 条断言一个字都没有改动。 */
#include "flashvec.h"

#include <math.h>
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

/* 位模式相等：+0.0 与 -0.0 必须区分开，所以不能直接用 ==。 */
static int bits_eq(float a, float b)
{
    uint32_t x, y;
    memcpy(&x, &a, 4);
    memcpy(&y, &b, 4);
    return x == y;
}

int main(void)
{
    const fv_vec2 zero = fv_vec2_make(0.0f, 0.0f);
    const float nan = (float)NAN;

    /* ---------------- vec2 min / max：基本语义 ---------------- */
    {
        fv_vec2 a = fv_vec2_make(1.0f, -2.0f);
        fv_vec2 b = fv_vec2_make(3.0f, -5.0f);
        fv_vec2 lo = fv_vec2_min(a, b);
        fv_vec2 hi = fv_vec2_max(a, b);

        CHECK(bits_eq(lo.x, 1.0f) && bits_eq(lo.y, -5.0f));
        CHECK(bits_eq(hi.x, 3.0f) && bits_eq(hi.y, -2.0f));
        /* 参数相等时必须原样返回（含位模式） */
        CHECK(bits_eq(fv_vec2_min(a, a).x, 1.0f) && bits_eq(fv_vec2_max(a, a).y, -2.0f));
        /* 对称性：min(a,b) == -max(-a,-b) */
        {
            fv_vec2 m = fv_vec2_min(a, b);
            fv_vec2 n = fv_vec2_neg(fv_vec2_max(fv_vec2_neg(a), fv_vec2_neg(b)));
            CHECK(bits_eq(m.x, n.x) && bits_eq(m.y, n.y));
        }
    }

    /* ---------------- min / max：±0.0 与 NaN ---------------- */
    {
        fv_vec2 pz = fv_vec2_make(0.0f, -0.0f);
        fv_vec2 nz = fv_vec2_make(-0.0f, 0.0f);
        /* 定义是 (a < b) ? a : b：+0.0 < -0.0 为假 -> 取 b 的 -0.0 */
        CHECK(bits_eq(fv_vec2_min(pz, nz).x, -0.0f));
        CHECK(bits_eq(fv_vec2_max(pz, nz).x, -0.0f));
        CHECK(bits_eq(fv_vec2_min(pz, nz).y, 0.0f));

        /* NaN 不参与比较：(a < b) 为假 -> 取 b */
        {
            fv_vec2 na = fv_vec2_make(nan, 1.0f);
            fv_vec2 v2 = fv_vec2_make(2.0f, 2.0f);
            CHECK(bits_eq(fv_vec2_min(na, v2).x, 2.0f));   /* 左 NaN -> 取右 */
            CHECK(isnan(fv_vec2_min(v2, na).x));           /* 右 NaN -> 取右（NaN） */
            CHECK(bits_eq(fv_vec2_max(na, v2).x, 2.0f));   /* 左 NaN -> 取右（max 也是 (a>b)?a:b） */
            CHECK(bits_eq(fv_vec2_max(v2, na).y, 2.0f));
        }
    }

    /* ---------------- vec2 clamp ---------------- */
    {
        const fv_vec2 lo = fv_vec2_make(-1.0f, -1.0f);
        const fv_vec2 hi = fv_vec2_make(1.0f, 1.0f);
        fv_vec2 inside = fv_vec2_clamp(fv_vec2_make(0.5f, -0.25f), lo, hi);
        fv_vec2 outside = fv_vec2_clamp(fv_vec2_make(-9.0f, 9.0f), lo, hi);

        CHECK(bits_eq(inside.x, 0.5f) && bits_eq(inside.y, -0.25f));
        CHECK(bits_eq(outside.x, -1.0f) && bits_eq(outside.y, 1.0f));
        /* 恰好在边界上：原样返回（不折半、不漂移） */
        CHECK(bits_eq(fv_vec2_clamp(lo, lo, hi).x, -1.0f));
        CHECK(bits_eq(fv_vec2_clamp(hi, lo, hi).y, 1.0f));
        /* lo == hi：夹成那个值 */
        CHECK(bits_eq(fv_vec2_clamp(fv_vec2_make(5.0f, -5.0f), zero, zero).x, 0.0f));
        CHECK(bits_eq(fv_vec2_clamp(fv_vec2_make(5.0f, -5.0f), zero, zero).y, 0.0f));
        /* lo > hi（调用方错误）：按同一表达式落位——x < lo 给 lo，否则 x > hi 给 hi。
         * 取 v = (2, 0)、边界 lo' = (1,1)、hi' = (-1,-1)：x 走上界分支，y 走下界分支。 */
        {
            fv_vec2 bad = fv_vec2_clamp(fv_vec2_make(2.0f, 0.0f), hi, lo);
            CHECK(bits_eq(bad.x, -1.0f) && bits_eq(bad.y, 1.0f));
        }
        /* NaN 与两个边界比较都为假 -> 原样返回 NaN */
        CHECK(isnan(fv_vec2_clamp(fv_vec2_make(nan, 0.0f), lo, hi).x));
        CHECK(bits_eq(fv_vec2_clamp(fv_vec2_make(nan, 0.0f), lo, hi).y, 0.0f));
        /* 逐位等于 (x < lo ? lo : (x > hi ? hi : x)) */
        {
            const float x = 0.1f;
            const float want = (x < lo.x) ? lo.x : ((x > hi.x) ? hi.x : x);
            CHECK(bits_eq(fv_vec2_clamp(fv_vec2_make(x, x), lo, hi).x, want));
        }
    }

    /* ---------------- vec2 / vec3 madd：逐位等于 a + b*s ---------------- */
    {
        const fv_vec2 a = fv_vec2_make(3.5f, -7.25f);
        const fv_vec2 b = fv_vec2_make(1.5f, 2.0e-3f);
        const fv_vec3 a3 = fv_vec3_make(1.0f / 3.0f, -2.0f, 1.0e-30f);
        const fv_vec3 b3 = fv_vec3_make(7.0f, 0.1f, -3.0f);
        const float s = 0.1f;

        CHECK(bits_eq(fv_vec2_madd(a, b, s).x, (float)(a.x + b.x * s)));
        CHECK(bits_eq(fv_vec2_madd(a, b, s).y, (float)(a.y + b.y * s)));
        CHECK(bits_eq(fv_vec3_madd(a3, b3, s).x, (float)(a3.x + b3.x * s)));
        CHECK(bits_eq(fv_vec3_madd(a3, b3, s).z, (float)(a3.z + b3.z * s)));
        /* s = 0 时逐位等于 a */
        CHECK(bits_eq(fv_vec2_madd(a, b, 0.0f).x, a.x));
        CHECK(bits_eq(fv_vec3_madd(a3, b3, 0.0f).y, a3.y));
        /* NaN 输入按 IEEE 规则传播（不做特殊处理） */
        CHECK(isnan(fv_vec2_madd(fv_vec2_make(nan, 0.0f), b, s).x));
        CHECK(bits_eq(fv_vec2_madd(fv_vec2_make(nan, 0.0f), b, s).y, (float)(0.0f + b.y * s)));
    }

    /* ---------------- vec3 min / max / clamp ---------------- */
    {
        const fv_vec3 a = fv_vec3_make(1.0f, -2.0f, 0.5f);
        const fv_vec3 b = fv_vec3_make(3.0f, -5.0f, -0.5f);
        fv_vec3 lo = fv_vec3_min(a, b);
        fv_vec3 hi = fv_vec3_max(a, b);
        fv_vec3 cl = fv_vec3_clamp(fv_vec3_make(9.0f, -0.25f, -9.0f), lo, hi);

        CHECK(bits_eq(lo.x, 1.0f) && bits_eq(lo.y, -5.0f) && bits_eq(lo.z, -0.5f));
        CHECK(bits_eq(hi.x, 3.0f) && bits_eq(hi.y, -2.0f) && bits_eq(hi.z, 0.5f));
        CHECK(bits_eq(cl.x, 3.0f) && bits_eq(cl.y, -2.0f) && bits_eq(cl.z, -0.5f));
        CHECK(bits_eq(fv_vec3_min(fv_vec3_make(nan, 0.0f, 0.0f), b).x, 3.0f)); /* 左 NaN -> 取右 */
        CHECK(bits_eq(fv_vec3_max(a, a).z, 0.5f));
    }

    /* ---------------- AABB ---------------- */
    {
        const fv_aabb2 box = { { 0.0f, 0.0f }, { 1.0f, 1.0f } };
        fv_aabb2 grown = fv_aabb2_merge(box, fv_vec2_make(-2.0f, 0.5f));
        fv_aabb2 same = fv_aabb2_merge(box, fv_vec2_make(0.5f, 0.5f));
        fv_aabb2 touching = { { 1.0f, 1.0f }, { 2.0f, 2.0f } };
        fv_aabb2 apart = { { 1.5f, 1.0f }, { 2.0f, 2.0f } };
        fv_aabb2 apart_y = { { 0.0f, 1.5f }, { 1.0f, 2.0f } };
        fv_aabb2 inside = { { 0.25f, 0.25f }, { 0.75f, 0.75f } };
        fv_aabb2 degenerate = { { 0.5f, 0.5f }, { 0.5f, 0.5f } };
        fv_aabb2 elsewhere = { { 9.0f, 9.0f }, { 9.0f, 9.0f } };

        /* merge 只扩不放 */
        CHECK(bits_eq(grown.lo.x, -2.0f) && bits_eq(grown.lo.y, 0.0f));
        CHECK(bits_eq(grown.hi.x, 1.0f) && bits_eq(grown.hi.y, 1.0f));
        /* 内部点：盒子四个分量逐位不变 */
        CHECK(memcmp(&same, &box, sizeof(fv_aabb2)) == 0);
        /* 合并一个极端点后两个方向都被撑开 */
        {
            fv_aabb2 wide = fv_aabb2_merge(fv_aabb2_merge(box, fv_vec2_make(-1.0f, 0.0f)),
                                           fv_vec2_make(0.0f, 4.0f));
            CHECK(bits_eq(wide.lo.x, -1.0f) && bits_eq(wide.lo.y, 0.0f));
            CHECK(bits_eq(wide.hi.x, 1.0f) && bits_eq(wide.hi.y, 4.0f));
        }
        /* overlaps：相交 / 相切 / 分离 / 包含 / 退化 */
        CHECK(fv_aabb2_overlaps(box, inside) == 1);
        CHECK(fv_aabb2_overlaps(box, touching) == 1);   /* 相切算相交 */
        CHECK(fv_aabb2_overlaps(box, degenerate) == 1);
        CHECK(fv_aabb2_overlaps(box, box) == 1);
        CHECK(fv_aabb2_overlaps(apart, box) == 0);
        CHECK(fv_aabb2_overlaps(box, apart_y) == 0);
        CHECK(fv_aabb2_overlaps(degenerate, degenerate) == 1);
        CHECK(fv_aabb2_overlaps(degenerate, elsewhere) == 0);
        /* 退化盒子与自身相交、与分离的退化盒子不相交 */
        CHECK(fv_aabb2_overlaps(elsewhere, elsewhere) == 1);
    }

    /* ---------------- FV_EPS / FV_EPS_SQR 量纲 ---------------- */
    {
        /* 两个常量差 6 个数量级；用来比长度平方时混用会把阈值放宽 1000 倍 */
        CHECK(FV_EPS_SQR < FV_EPS);
        CHECK(FV_EPS * FV_EPS <= FV_EPS_SQR * 1.01f);
        CHECK(FV_EPS * FV_EPS >= FV_EPS_SQR * 0.99f);
    }

    /* ---------------- project 的阈值修正（vec2 / vec3） ---------------- */
    {
        /* onto 长度 1e-4 -> 长度平方 1e-8：旧的线性阈值 FV_EPS=1e-6 会误判成零向量，
         * 新的平方阈值 FV_EPS_SQR=1e-12 不会，投影结果回到正确值附近。 */
        fv_vec2 p2 = fv_vec2_project(fv_vec2_make(1.0f, 1.0f), fv_vec2_make(1.0e-4f, 0.0f));
        CHECK(fv_vec2_near(p2, fv_vec2_make(1.0f, 0.0f), 1.0e-3f));

        fv_vec3 p3 = fv_vec3_project(fv_vec3_make(1.0f, 1.0f, 1.0f),
                                     fv_vec3_make(1.0e-4f, 0.0f, 0.0f));
        CHECK(fv_vec3_near(p3, fv_vec3_make(1.0f, 0.0f, 0.0f), 1.0e-3f));

        /* 长度 1e-7 -> 长度平方 1e-14 <= FV_EPS_SQR：仍然判成零向量 */
        CHECK(fv_vec2_near(fv_vec2_project(fv_vec2_make(1.0f, 1.0f),
                                           fv_vec2_make(1.0e-7f, 0.0f)),
                           zero, 0.0f));
        CHECK(fv_vec3_near(fv_vec3_project(fv_vec3_make(1.0f, 1.0f, 1.0f),
                                           fv_vec3_make(1.0e-7f, 0.0f, 0.0f)),
                           fv_vec3_make(0.0f, 0.0f, 0.0f), 0.0f));
        /* 零向量、垂直向量的行为不变 */
        CHECK(fv_vec2_near(fv_vec2_project(fv_vec2_make(3.0f, 4.0f), zero), zero, 0.0f));
        CHECK(fv_vec2_near(fv_vec2_project(fv_vec2_make(0.0f, 2.0f),
                                           fv_vec2_make(3.0f, 0.0f)),
                           zero, 0.0f));
    }

    printf("test_p0: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
