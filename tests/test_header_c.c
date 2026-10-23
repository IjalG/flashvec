/* 头文件可用性（纯 C 编译单元）—— 验收项 3 与 5。
 *
 * 这个文件证明三件事：
 *   1) 头文件在 cc -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic 下可用；
 *   2) 同一个 TU 里两次 include 不报错（include guard 必须让它成为空操作）；
 *   3) C 编译单元在 x86 + GCC/Clang 下**必须**走 intrinsics 路径 —— 用编译期断言锁住，不靠 grep
 *      （语言不再是 gate 条件；只有不支持 target 属性的编译器/非 x86 才退回纯 C 路径）。
 * 顺便在运行期拿默认（inline）形态的结果与标量 API 逐位对拍：
 * 这里故意**不加** -ffp-contract=off，让头文件自己的 #pragma GCC optimize 去挡 FMA 收缩；
 * 如果那个 pragma 失效，madd / lerp / len / normalize 的位模式就会和标量 API 不同。
 *
 * 编译（严格档）：cc -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -O2 -c tests/test_header_c.c
 * 链接：需要 libflashvec.a 提供标量 API（批量入口在本文件内以 inline 形态展开）。
 */
#include "flashvec.h"
#include "flashvec.h" /* 双重包含：guard 命中，第二次必须是空操作 */

/* 契约（见 flashvec.h 文件头）：intrinsics 的出现条件与语言无关，只与编译器/架构有关。
 * x86 + GCC/Clang 下，C 编译单元与 C++ 编译单元一样必须拿到 intrinsics 路径 —— 用编译期
 * 断言锁住，防止它悄悄退化成只有 SSE2 的纯 C 路径（那是实打实的性能回退）。 */
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#if !FV_I_INTRINSICS
#error "C translation unit on x86 GCC/Clang must take the intrinsics path (contract in flashvec.h)"
#endif
#endif

#include <stdio.h>
#include <string.h>

#define CAP 40

static int bits_eq(float a, float b)
{
    unsigned char ba[4], bb[4];
    memcpy(ba, &a, 4);
    memcpy(bb, &b, 4);
    return memcmp(ba, bb, 4) == 0;
}

/* 有限随机值（避开 NaN/Inf），同时掺入零向量、-0.0、次正规、极值。 */
static const float kInteresting[] = { 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f, 1.0e-6f,
                                      1.0e-12f, 4.9e-40f, -1.0e-30f, 3.0e38f, 2.0f };
#define NINTEREST ((int)(sizeof(kInteresting) / sizeof(kInteresting[0])))

static unsigned int g_rng;

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

static void fill_vec2(fv_vec2 *v, size_t n, unsigned int seed)
{
    size_t i;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        v[i].x = (i % 3u == 0u) ? kInteresting[i % (size_t)NINTEREST] : rnd_float();
        v[i].y = (i % 3u == 1u) ? kInteresting[(i + 5u) % (size_t)NINTEREST] : rnd_float();
    }
}

int main(void)
{
    static fv_vec2 a[CAP], b[CAP], o[CAP];
    static float ax[CAP], ay[CAP], bx[CAP], by[CAP], ox[CAP], oy[CAP];
    unsigned char mask[CAP];
    unsigned int idx[CAP];
    const float params[3] = { 0.5f, -3.25f, 1.0e-6f };
    int checks = 0, failures = 0;
    size_t n, i;
    int pi;

#define CHECK(cond)                                                     \
    do {                                                                \
        ++checks;                                                       \
        if (!(cond)) {                                                  \
            ++failures;                                                 \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);       \
        }                                                               \
    } while (0)

    for (n = 0; n <= 17; ++n) {
        fill_vec2(a, n, 0x1234u + (unsigned int)n);
        fill_vec2(b, n, 0x5678u + (unsigned int)n);
        for (i = 0; i < n; ++i) mask[i] = (unsigned char)((i % 3u == 0u) ? 0u : 1u);

        for (pi = 0; pi < 3; ++pi) {
            const float s = params[pi];
            int ok = 1;

            memset(o, 0x5A, sizeof(o));
            fv_vec2_scale_x4(a, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_scale(a[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok);

            ok = 1;
            memset(o, 0x5A, sizeof(o));
            fv_vec2_madd_x4(a, b, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_madd(a[i], b[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok); /* 关键：FMA 收缩会在这里露馅 */

            ok = 1;
            memset(o, 0x5A, sizeof(o));
            fv_vec2_lerp_x4(a, b, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_lerp(a[i], b[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        {
            int ok = 1;
            memset(o, 0x5A, sizeof(o));
            fv_vec2_add_x4(a, b, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_add(a[i], b[i]);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        {
            int ok = 1;
            memset(o, 0x5A, sizeof(o));
            fv_vec2_add_x4_mask(a, b, o, mask, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = (mask[i] != 0u) ? fv_vec2_add(a[i], b[i]) : a[i];
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        /* 索引：逆序 + 重复，m == 0 时什么都不写 */
        {
            int ok = 1;
            const size_t m = (n > 1u) ? n : 0u;
            for (i = 0; i < m; ++i) idx[i] = (unsigned int)(m - 1u - i); /* 逆序 */
            if (m >= 2u) idx[m - 1u] = idx[0];                          /* 再加一处重复 */
            memset(o, 0x5A, sizeof(o));
            memcpy(o, a, n * sizeof(fv_vec2));
            fv_vec2_add_x4_idx(a, b, o, idx, m);
            for (i = 0; i < m; ++i) {
                const unsigned int k = idx[i];
                const fv_vec2 r = fv_vec2_add(a[k], b[k]);
                if (!bits_eq(o[k].x, r.x) || !bits_eq(o[k].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        /* SoA */
        for (i = 0; i < n; ++i) {
            ax[i] = a[i].x;
            ay[i] = a[i].y;
            bx[i] = b[i].x;
            by[i] = b[i].y;
        }
        {
            int ok_add = 1, ok_scale = 1, ok_len = 1, ok_norm = 1;
            memset(ox, 0x5A, sizeof(ox));
            memset(oy, 0x5A, sizeof(oy));
            fv_vec2_add_soa(ax, ay, bx, by, ox, oy, n);
            fv_vec2_scale_soa(ax, ay, 0.5f, ox + 0, oy + 0, 0); /* n == 0：不读不写 */
            memset(ox, 0x5A, sizeof(ox));
            memset(oy, 0x5A, sizeof(oy));
            fv_vec2_add_soa(ax, ay, bx, by, ox, oy, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 ari = fv_vec2_add(fv_vec2_make(ax[i], ay[i]),
                                               fv_vec2_make(bx[i], by[i]));
                if (!bits_eq(ox[i], ari.x) || !bits_eq(oy[i], ari.y)) ok_add = 0;
            }
            CHECK(ok_add);

            memset(ox, 0x5A, sizeof(ox));
            memset(oy, 0x5A, sizeof(oy));
            fv_vec2_scale_soa(ax, ay, 0.5f, ox, oy, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 ari = fv_vec2_scale(fv_vec2_make(ax[i], ay[i]), 0.5f);
                if (!bits_eq(ox[i], ari.x) || !bits_eq(oy[i], ari.y)) ok_scale = 0;
            }
            CHECK(ok_scale);

            memset(ox, 0x5A, sizeof(ox));
            fv_vec2_len_soa(ax, ay, ox, n);
            for (i = 0; i < n; ++i) {
                if (!bits_eq(ox[i], fv_vec2_len(fv_vec2_make(ax[i], ay[i])))) ok_len = 0;
            }
            CHECK(ok_len);

            memset(ox, 0x5A, sizeof(ox));
            memset(oy, 0x5A, sizeof(oy));
            fv_vec2_normalize_soa(ax, ay, ox, oy, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_normalize(fv_vec2_make(ax[i], ay[i]));
                if (!bits_eq(ox[i], r.x) || !bits_eq(oy[i], r.y)) ok_norm = 0;
            }
            CHECK(ok_norm);
        }

        /* dot：n < 8 时纯顺序累加，必须与标量逐位一致 */
        if (n < 8u) {
            float ref = 0.0f;
            for (i = 0; i < n; ++i) ref = ref + (ax[i] * bx[i] + ay[i] * by[i]);
            CHECK(bits_eq(fv_vec2_dot_soa(ax, ay, bx, by, n), ref));
        }
    }

    printf("header(C11): FV_I_INTRINSICS=%d  %d checks, %d failures\n",
           FV_I_INTRINSICS, checks, failures);
    return failures == 0 ? 0 : 1;
}
