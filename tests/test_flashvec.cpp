/* flashvec 的单元测试：纯 C 风格断言，不依赖任何测试框架。 */
#include "flashvec.h"

#include <stdio.h>

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

#define EPS ((fv_real)1e-5)

int main(void)
{
    /* 加减、数乘、取负 */
    {
        fv_vec2 a = fv_vec2_make(1, 2);
        fv_vec2 b = fv_vec2_make(3, 5);
        CHECK(fv_vec2_near(fv_vec2_add(a, b), fv_vec2_make(4, 7), EPS));
        CHECK(fv_vec2_near(fv_vec2_sub(b, a), fv_vec2_make(2, 3), EPS));
        CHECK(fv_vec2_near(fv_vec2_neg(a), fv_vec2_make(-1, -2), EPS));
        CHECK(fv_vec2_near(fv_vec2_scale(a, 2.5f), fv_vec2_make(2.5f, 5), EPS));
        CHECK(fv_vec2_near(fv_vec2_scale(a, 0), fv_vec2_make(0, 0), EPS));
    }

    /* 点积、长度、距离 */
    CHECK(fv_vec2_dot(fv_vec2_make(1, 0), fv_vec2_make(0, 1)) == 0);
    CHECK(fv_vec2_dot(fv_vec2_make(1, 2), fv_vec2_make(3, 4)) == 11);
    CHECK(fv_vec2_len_sqr(fv_vec2_make(3, 4)) == 25);
    CHECK(fv_vec2_len(fv_vec2_make(3, 4)) > 4.999f &&
          fv_vec2_len(fv_vec2_make(3, 4)) < 5.001f);
    CHECK(fv_vec2_dist(fv_vec2_make(1, 1), fv_vec2_make(4, 5)) > 4.999f &&
          fv_vec2_dist(fv_vec2_make(1, 1), fv_vec2_make(4, 5)) < 5.001f);
    /* 勾股：3-4-5 是精确整数，但 1-1-sqrt2 只能近似 */
    CHECK(fv_vec2_dist(fv_vec2_make(0, 0), fv_vec2_make(1, 0)) > 0.999f);

    /* 归一化：结果单位长、方向不变；零向量不炸 */
    {
        fv_vec2 n = fv_vec2_normalize(fv_vec2_make(3, 4));
        CHECK(fv_vec2_near(n, fv_vec2_make(0.6f, 0.8f), EPS));
        CHECK(fv_vec2_len(n) > 0.999f && fv_vec2_len(n) < 1.001f);
        CHECK(fv_vec2_near(fv_vec2_normalize(fv_vec2_make(0, 0)),
                           fv_vec2_make(0, 0), EPS));
    }

    /* 投影：v 在 x 轴上的投影是 (v.x, 0) */
    {
        fv_vec2 p = fv_vec2_project(fv_vec2_make(3, 4), fv_vec2_make(1, 0));
        CHECK(fv_vec2_near(p, fv_vec2_make(3, 0), EPS));
        CHECK(fv_vec2_near(fv_vec2_project(fv_vec2_make(3, 4),
                                           fv_vec2_make(0, 0)),
                           fv_vec2_make(0, 0), EPS));
        /* 投影 + 垂直分量 = 原向量 */
        {
            fv_vec2 v = fv_vec2_make(3, 4);
            fv_vec2 onto = fv_vec2_make(1, 1);
            fv_vec2 perp = fv_vec2_sub(v, fv_vec2_project(v, onto));
            CHECK(fv_vec2_dot(perp, onto) < EPS &&
                  fv_vec2_dot(perp, onto) > -EPS);
        }
    }

    /* 线性插值 */
    {
        fv_vec2 a = fv_vec2_make(0, 0);
        fv_vec2 b = fv_vec2_make(10, 20);
        CHECK(fv_vec2_near(fv_vec2_lerp(a, b, 0), a, EPS));
        CHECK(fv_vec2_near(fv_vec2_lerp(a, b, 1), b, EPS));
        CHECK(fv_vec2_near(fv_vec2_lerp(a, b, 0.5f), fv_vec2_make(5, 10), EPS));
    }

    /* 三维：点积、叉积、长度 */
    {
        fv_vec3 ex = fv_vec3_make(1, 0, 0);
        fv_vec3 ey = fv_vec3_make(0, 1, 0);
        fv_vec3 ez = fv_vec3_make(0, 0, 1);
        CHECK(fv_vec3_near(fv_vec3_cross(ex, ey), ez, EPS));
        CHECK(fv_vec3_near(fv_vec3_cross(ey, ex), fv_vec3_neg(ez), EPS));
        CHECK(fv_vec3_near(fv_vec3_cross(ex, ex), fv_vec3_make(0, 0, 0), EPS));
        CHECK(fv_vec3_dot(ex, ey) == 0);
        CHECK(fv_vec3_dot(fv_vec3_make(1, 2, 3), fv_vec3_make(4, 5, 6)) == 32);
        CHECK(fv_vec3_len(fv_vec3_make(2, 3, 6)) > 6.999f &&
              fv_vec3_len(fv_vec3_make(2, 3, 6)) < 7.001f);
        CHECK(fv_vec3_near(fv_vec3_normalize(fv_vec3_make(0, 3, 4)),
                           fv_vec3_make(0, 0.6f, 0.8f), EPS));
        CHECK(fv_vec3_near(fv_vec3_normalize(fv_vec3_make(0, 0, 0)),
                           fv_vec3_make(0, 0, 0), EPS));
        CHECK(fv_vec3_near(fv_vec3_lerp(fv_vec3_make(0, 0, 0),
                                        fv_vec3_make(2, 4, 6), 0.5f),
                           fv_vec3_make(1, 2, 3), EPS));
        CHECK(fv_vec3_dist(fv_vec3_make(0, 0, 0), fv_vec3_make(0, 0, 1)) > 0.999f);
    }

    printf("test_flashvec: %d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
