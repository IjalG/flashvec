/* 头文件可用性（C++ 编译单元）—— 验收项 4 与 5。
 *
 * 与 tests/test_header_c.c 同构，但用 c++ -std=c++17 编译：
 *   - 两次 include 不报错；
 *   - -Wall -Wextra -Wpedantic 下无警告；
 *   - x86 + GCC/Clang 的 C++ TU 必须走 intrinsics 路径（编译期断言），
 *     但照样要与标量 API 位级一致 —— 包括 madd / lerp 这两个最容易被 FMA 收缩的入口。
 * 这里同样故意不加 -ffp-contract=off。
 */
#include "flashvec.h"
#include "flashvec.h" /* 双重包含 */

#if defined(__x86_64__) || defined(__i386__)
#if (defined(__GNUC__) || defined(__clang__)) && !FV_I_INTRINSICS
#error "C++ TU on x86 GCC/Clang must expose the intrinsics path (see flashvec.h)"
#endif
#endif

#include <cstdio>
#include <cstring>

#define CAP 40

static int bits_eq(float a, float b)
{
    unsigned int x, y;
    std::memcpy(&x, &a, 4);
    std::memcpy(&y, &b, 4);
    return x == y;
}

static unsigned int g_rng;

static unsigned int rnd_u32()
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

static float rnd_float()
{
    const unsigned int b = rnd_u32();
    const unsigned int bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}

static void fill_vec2(fv_vec2 *v, std::size_t n, unsigned int seed)
{
    static const float interesting[] = { 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f, 1.0e-6f,
                                         1.0e-12f, 4.9e-40f, -1.0e-30f, 3.0e38f };
    const std::size_t nint = sizeof(interesting) / sizeof(interesting[0]);
    std::size_t i;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        v[i].x = (i % 3u == 0u) ? interesting[i % nint] : rnd_float();
        v[i].y = (i % 3u == 1u) ? interesting[(i + 5u) % nint] : rnd_float();
    }
}

int main()
{
    static fv_vec2 a[CAP], b[CAP], o[CAP];
    static float ax[CAP], ay[CAP], bx[CAP], by[CAP], ox[CAP], oy[CAP];
    unsigned char mask[CAP];
    unsigned int idx[CAP];
    const float params[3] = { 0.5f, -3.25f, 1.0e-6f };
    int checks = 0, failures = 0;
    std::size_t n, i;
    int pi;

#define CHECK(cond)                                                     \
    do {                                                                \
        ++checks;                                                       \
        if (!(cond)) {                                                  \
            ++failures;                                                 \
            std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

    for (n = 0; n <= 17; ++n) {
        fill_vec2(a, n, 0x9ABCu + static_cast<unsigned int>(n));
        fill_vec2(b, n, 0xDEF0u + static_cast<unsigned int>(n));
        for (i = 0; i < n; ++i) mask[i] = static_cast<unsigned char>((i % 3u == 0u) ? 0u : 2u);

        for (pi = 0; pi < 3; ++pi) {
            const float s = params[pi];
            int ok_madd = 1, ok_lerp = 1, ok_scale = 1;

            std::memset(o, 0x5A, sizeof(o));
            fv_vec2_madd_x4(a, b, s, o, n);
            std::memset(o, 0x5A, sizeof(o));
            fv_vec2_madd_x4(a, b, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_madd(a[i], b[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok_madd = 0;
            }
            CHECK(ok_madd);

            std::memset(o, 0x5A, sizeof(o));
            fv_vec2_lerp_x4(a, b, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_lerp(a[i], b[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok_lerp = 0;
            }
            CHECK(ok_lerp);

            std::memset(o, 0x5A, sizeof(o));
            fv_vec2_scale_x4(a, s, o, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_scale(a[i], s);
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok_scale = 0;
            }
            CHECK(ok_scale);
        }

        {
            int ok = 1;
            std::memset(o, 0x5A, sizeof(o));
            fv_vec2_add_x4_mask(a, b, o, mask, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = (mask[i] != 0u) ? fv_vec2_add(a[i], b[i]) : a[i];
                if (!bits_eq(o[i].x, r.x) || !bits_eq(o[i].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        {
            int ok = 1;
            const std::size_t m = (n > 1u) ? n : 0u;
            unsigned int tmp;
            for (i = 0; i < m; ++i) idx[i] = static_cast<unsigned int>(m - 1u - i);
            for (i = 0; i < m / 2u; ++i) { /* 逆序之后再乱序打散 */
                const std::size_t j = (i * 7u + 3u) % m;
                tmp = idx[i];
                idx[i] = idx[j];
                idx[j] = tmp;
            }
            std::memcpy(o, a, n * sizeof(fv_vec2));
            fv_vec2_add_x4_idx(a, b, o, idx, m);
            for (i = 0; i < m; ++i) {
                const unsigned int k = idx[i];
                const fv_vec2 r = fv_vec2_add(a[k], b[k]);
                if (!bits_eq(o[k].x, r.x) || !bits_eq(o[k].y, r.y)) ok = 0;
            }
            CHECK(ok);
        }

        for (i = 0; i < n; ++i) {
            ax[i] = a[i].x;
            ay[i] = a[i].y;
            bx[i] = b[i].x;
            by[i] = b[i].y;
        }
        {
            int ok_len = 1, ok_norm = 1, ok_add = 1;
            std::memset(ox, 0x5A, sizeof(ox));
            fv_vec2_len_soa(ax, ay, ox, n);
            for (i = 0; i < n; ++i) {
                if (!bits_eq(ox[i], fv_vec2_len(fv_vec2_make(ax[i], ay[i])))) ok_len = 0;
            }
            CHECK(ok_len);

            std::memset(ox, 0x5A, sizeof(ox));
            std::memset(oy, 0x5A, sizeof(oy));
            fv_vec2_normalize_soa(ax, ay, ox, oy, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_normalize(fv_vec2_make(ax[i], ay[i]));
                if (!bits_eq(ox[i], r.x) || !bits_eq(oy[i], r.y)) ok_norm = 0;
            }
            CHECK(ok_norm);

            std::memset(ox, 0x5A, sizeof(ox));
            std::memset(oy, 0x5A, sizeof(oy));
            fv_vec2_add_soa(ax, ay, bx, by, ox, oy, n);
            for (i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_add(fv_vec2_make(ax[i], ay[i]),
                                              fv_vec2_make(bx[i], by[i]));
                if (!bits_eq(ox[i], r.x) || !bits_eq(oy[i], r.y)) ok_add = 0;
            }
            CHECK(ok_add);
        }

        if (n < 8u) {
            float ref = 0.0f;
            for (i = 0; i < n; ++i) ref = ref + (ax[i] * bx[i] + ay[i] * by[i]);
            CHECK(bits_eq(fv_vec2_dot_soa(ax, ay, bx, by, n), ref));
        }
    }

    std::printf("header(C++): FV_I_INTRINSICS=%d  %d checks, %d failures\n",
                FV_I_INTRINSICS, checks, failures);
    return failures == 0 ? 0 : 1;
}
