/* 三维 SoA / 四维 _x2、_soa 批量入口的位级对拍（与 tests/test_parity.c 同一套做法）。
 *
 * 判据：
 *   1) 每个非 dot 入口，在 n = 0..17 全尾长 + 63 + 1024 + 8192 下，输出与「逐元素调用标量 API」
 *      的结果逐字节 memcmp 相等，比较长度是 (n + GUARD) 个元素 —— 缓冲区在调用前被 memset(0x5A)
 *      毒化，所以越界写和「该写没写」都会被抓到；
 *   2) dot_soa（三维/四维）是契约里唯一的例外：多累加器 + 末端归约，只承诺确定性与小误差。
 *      容差口径沿用 tests/test_batch.cpp：以「各项绝对值之和」为 scale 归一化的相对误差 < 1e-5；
 *      另外 n < 8 时整段只走顺序尾巴，额外要求与标量顺序累加**逐位**一致；
 *   3) 结果按固定格式写进两条二进制流（主流程 / dot 流），供跨形态（C/C++ × inline/extern ×
 *      纯 C 回退 × FMA 消费者）逐字节 cmp —— 这一层抓「同形态内部对不上」以外的问题，
 *      比如某形态的 tail 走了不同表达式。
 *
 * 语料：随机有限值（bit 级随机，指数取 1..254，避免 NaN/inf）+ 尾部特殊值；
 *       NaN 单独一组（单侧 NaN 严格比较；两侧都是 NaN 时 payload 由指令操作数顺序决定，
 *       契约不承诺，见 tests/test_vec4.c 里同样的说明）。
 *
 * 用法：test_vec34_parity <main.bin> <dot.bin>（make test-vec34-inline / test-vec34-fma 调用）
 *
 * 本轮实测（库内 -O3 -fno-math-errno -ffp-contract=off；见报告）：
 *   cases=956 dot_cases=84 failures=0；max_rel vec3=5.650e-08、vec4=3.159e-08。
 */
#include "flashvec.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NV4 8200u   /* 8192 + GUARD 余量 */
#define GUARD 8u

static fv_vec4 g_a4[NV4], g_b4[NV4], g_o4[NV4], g_r4[NV4];
static float g3a[3][NV4], g3b[3][NV4], g3o[3][NV4], g3r[3][NV4];
static float g4a[4][NV4], g4b[4][NV4], g4o[4][NV4], g4r[4][NV4];

static unsigned g_rng;
static unsigned g_case;
static unsigned g_dot_cases;
static unsigned g_fail;
static double g_rel3;
static double g_rel4;

static void chk_bitwise(int ok, const char *what, size_t n)
{
    if (ok == 0) {
        ++g_fail;
        printf("FAIL %-28s n=%u\n", what, (unsigned int)n);
    }
}

static unsigned int rnd_u32(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 17;
    g_rng ^= g_rng << 5;
    return g_rng;
}

/* 有限值：指数取 1..254，永不产生 NaN/inf，所以算术对拍是严格位级的。 */
static float rnd_float(void)
{
    const unsigned int b = rnd_u32();
    const unsigned int bits = (b & 0x807FFFFFu) | ((((b >> 23) % 254u) + 1u) << 23);
    float f;
    memcpy(&f, &bits, 4);
    return f;
}

/* 中等量级、会互相抵消的值：dot 的容差判据只在结果有限时有意义（契约里的数字也是这个口径）。 */
static float rnd_moderate(void)
{
    return (float)((int)(rnd_u32() % 20001u) - 10000) / 1000.0f;
}

static float from_bits(unsigned int b)
{
    float f;
    memcpy(&f, &b, 4);
    return f;
}

static unsigned int fbits(float f)
{
    unsigned int b;
    memcpy(&b, &f, 4);
    return b;
}

static void poison(void *p, size_t bytes)
{
    memset(p, 0x5A, bytes);
}

static void rec(FILE *f, size_t n, unsigned int param, const void *p, size_t bytes)
{
    unsigned int hdr[4];
    hdr[0] = g_case++;
    hdr[1] = (unsigned int)n;
    hdr[2] = param;
    hdr[3] = (unsigned int)bytes;
    if (fwrite(hdr, sizeof(hdr[0]), 4, f) != 4 || (bytes != 0 && fwrite(p, 1, bytes, f) != bytes)) {
        perror("fwrite");
        exit(1);
    }
}

/* ---------------- 语料填充 ---------------- */

static void fill_aos4(fv_vec4 *v, size_t n, unsigned int seed)
{
    size_t i;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        v[i].x = (i % 5u == 0u) ? 0.0f : rnd_float();
        v[i].y = (i % 7u == 0u) ? -0.0f : rnd_float();
        v[i].z = (i % 11u == 0u) ? 1.0e-30f : rnd_float();
        v[i].w = (i % 13u == 0u) ? 3.0e38f : rnd_float();
    }
}

static void fill_soa(float (*a)[NV4], size_t n, unsigned int seed, int comps)
{
    size_t i;
    int c;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        for (c = 0; c < comps; ++c) {
            a[c][i] = (i % (size_t)(5 + c) == 0u) ? 0.0f : rnd_float();
        }
    }
}

/* dot 专用语料：中等量级（±10），乘积累加不溢出，容差判据才有意义。 */
static void fill_moderate(float (*a)[NV4], float (*b)[NV4], size_t n, unsigned int seed,
                          unsigned int seedb, int comps)
{
    size_t i;
    int c;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        for (c = 0; c < comps; ++c) {
            a[c][i] = rnd_moderate();
        }
    }
    g_rng = seedb;
    for (i = 0; i < n; ++i) {
        for (c = 0; c < comps; ++c) {
            b[c][i] = rnd_moderate();
        }
    }
}

/* ---------------- 四维 _x2：与标量 API 逐元素对拍 ---------------- */

static void t_v4_add(size_t n, FILE *fm)
{
    size_t i;
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_add_x2(g_a4, g_b4, g_o4, n);
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_add(g_a4[i], g_b4[i]);
    }
    chk_bitwise(memcmp(g_o4, g_r4, (n + GUARD) * sizeof(fv_vec4)) == 0, "vec4_add_x2", n);
    rec(fm, n, 0x4001u, g_o4, (n + GUARD) * sizeof(fv_vec4));
}

static void t_v4_sub(size_t n, FILE *fm)
{
    size_t i;
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_sub_x2(g_a4, g_b4, g_o4, n);
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_sub(g_a4[i], g_b4[i]);
    }
    chk_bitwise(memcmp(g_o4, g_r4, (n + GUARD) * sizeof(fv_vec4)) == 0, "vec4_sub_x2", n);
    rec(fm, n, 0x4002u, g_o4, (n + GUARD) * sizeof(fv_vec4));
}

static void t_v4_scale(size_t n, float s, FILE *fm)
{
    size_t i;
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_scale_x2(g_a4, s, g_o4, n);
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_scale(g_a4[i], s);
    }
    chk_bitwise(memcmp(g_o4, g_r4, (n + GUARD) * sizeof(fv_vec4)) == 0, "vec4_scale_x2", n);
    rec(fm, n, 0x4003u ^ fbits(s), g_o4, (n + GUARD) * sizeof(fv_vec4));
}

static void t_v4_lerp(size_t n, float t, FILE *fm)
{
    size_t i;
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_lerp_x2(g_a4, g_b4, t, g_o4, n);
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_lerp(g_a4[i], g_b4[i], t);
    }
    chk_bitwise(memcmp(g_o4, g_r4, (n + GUARD) * sizeof(fv_vec4)) == 0, "vec4_lerp_x2", n);
    rec(fm, n, 0x4004u ^ fbits(t), g_o4, (n + GUARD) * sizeof(fv_vec4));
}

static void t_v4_madd(size_t n, float s, FILE *fm)
{
    size_t i;
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_madd_x2(g_a4, g_b4, s, g_o4, n);
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_madd(g_a4[i], g_b4[i], s);
    }
    chk_bitwise(memcmp(g_o4, g_r4, (n + GUARD) * sizeof(fv_vec4)) == 0, "vec4_madd_x2", n);
    rec(fm, n, 0x4005u ^ fbits(s), g_o4, (n + GUARD) * sizeof(fv_vec4));
}

/* ---------------- 三维 SoA ---------------- */

static void t_v3_add(size_t n, FILE *fm)
{
    size_t i;
    int c;
    for (c = 0; c < 3; ++c) {
        poison(g3o[c], sizeof(g3o[c]));
        poison(g3r[c], sizeof(g3r[c]));
    }
    fv_vec3_add_soa(g3a[0], g3a[1], g3a[2], g3b[0], g3b[1], g3b[2],
                    g3o[0], g3o[1], g3o[2], n);
    for (i = 0; i < n; ++i) {
        const fv_vec3 r = fv_vec3_add(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]),
                                      fv_vec3_make(g3b[0][i], g3b[1][i], g3b[2][i]));
        g3r[0][i] = r.x;
        g3r[1][i] = r.y;
        g3r[2][i] = r.z;
    }
    chk_bitwise(memcmp(g3o[0], g3r[0], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[1], g3r[1], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[2], g3r[2], (n + GUARD) * sizeof(float)) == 0,
                "vec3_add_soa", n);
    for (c = 0; c < 3; ++c) {
        rec(fm, n, 0x4101u + (unsigned int)c, g3o[c], (n + GUARD) * sizeof(float));
    }
}

static void t_v3_scale(size_t n, float s, FILE *fm)
{
    size_t i;
    int c;
    for (c = 0; c < 3; ++c) {
        poison(g3o[c], sizeof(g3o[c]));
        poison(g3r[c], sizeof(g3r[c]));
    }
    fv_vec3_scale_soa(g3a[0], g3a[1], g3a[2], s, g3o[0], g3o[1], g3o[2], n);
    for (i = 0; i < n; ++i) {
        const fv_vec3 r = fv_vec3_scale(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]), s);
        g3r[0][i] = r.x;
        g3r[1][i] = r.y;
        g3r[2][i] = r.z;
    }
    chk_bitwise(memcmp(g3o[0], g3r[0], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[1], g3r[1], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[2], g3r[2], (n + GUARD) * sizeof(float)) == 0,
                "vec3_scale_soa", n);
    for (c = 0; c < 3; ++c) {
        rec(fm, n, 0x4102u ^ fbits(s) ^ (unsigned int)c, g3o[c], (n + GUARD) * sizeof(float));
    }
}

static void t_v3_len(size_t n, FILE *fm)
{
    size_t i;
    poison(g3o[0], sizeof(g3o[0]));
    poison(g3r[0], sizeof(g3r[0]));
    fv_vec3_len_soa(g3a[0], g3a[1], g3a[2], g3o[0], n);
    for (i = 0; i < n; ++i) {
        g3r[0][i] = fv_vec3_len(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]));
    }
    chk_bitwise(memcmp(g3o[0], g3r[0], (n + GUARD) * sizeof(float)) == 0, "vec3_len_soa", n);
    rec(fm, n, 0x4103u, g3o[0], (n + GUARD) * sizeof(float));
}

static void t_v3_normalize(size_t n, FILE *fm)
{
    size_t i;
    int c;
    for (c = 0; c < 3; ++c) {
        poison(g3o[c], sizeof(g3o[c]));
        poison(g3r[c], sizeof(g3r[c]));
    }
    fv_vec3_normalize_soa(g3a[0], g3a[1], g3a[2], g3o[0], g3o[1], g3o[2], n);
    for (i = 0; i < n; ++i) {
        const fv_vec3 r = fv_vec3_normalize(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]));
        g3r[0][i] = r.x;
        g3r[1][i] = r.y;
        g3r[2][i] = r.z;
    }
    chk_bitwise(memcmp(g3o[0], g3r[0], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[1], g3r[1], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g3o[2], g3r[2], (n + GUARD) * sizeof(float)) == 0,
                "vec3_normalize_soa", n);
    for (c = 0; c < 3; ++c) {
        rec(fm, n, 0x4104u + (unsigned int)c, g3o[c], (n + GUARD) * sizeof(float));
    }
}

/* ---------------- 四维 SoA：len / normalize（dot 单独一组） ---------------- */

static void t_v4_len(size_t n, FILE *fm)
{
    size_t i;
    poison(g4o[0], sizeof(g4o[0]));
    poison(g4r[0], sizeof(g4r[0]));
    fv_vec4_len_soa(g4a[0], g4a[1], g4a[2], g4a[3], g4o[0], n);
    for (i = 0; i < n; ++i) {
        g4r[0][i] = fv_vec4_len(fv_vec4_make(g4a[0][i], g4a[1][i], g4a[2][i], g4a[3][i]));
    }
    chk_bitwise(memcmp(g4o[0], g4r[0], (n + GUARD) * sizeof(float)) == 0, "vec4_len_soa", n);
    rec(fm, n, 0x4201u, g4o[0], (n + GUARD) * sizeof(float));
}

static void t_v4_normalize(size_t n, FILE *fm)
{
    size_t i;
    int c;
    for (c = 0; c < 4; ++c) {
        poison(g4o[c], sizeof(g4o[c]));
        poison(g4r[c], sizeof(g4r[c]));
    }
    fv_vec4_normalize_soa(g4a[0], g4a[1], g4a[2], g4a[3],
                          g4o[0], g4o[1], g4o[2], g4o[3], n);
    for (i = 0; i < n; ++i) {
        const fv_vec4 r = fv_vec4_normalize(
            fv_vec4_make(g4a[0][i], g4a[1][i], g4a[2][i], g4a[3][i]));
        g4r[0][i] = r.x;
        g4r[1][i] = r.y;
        g4r[2][i] = r.z;
        g4r[3][i] = r.w;
    }
    chk_bitwise(memcmp(g4o[0], g4r[0], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g4o[1], g4r[1], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g4o[2], g4r[2], (n + GUARD) * sizeof(float)) == 0
                    && memcmp(g4o[3], g4r[3], (n + GUARD) * sizeof(float)) == 0,
                "vec4_normalize_soa", n);
    for (c = 0; c < 4; ++c) {
        rec(fm, n, 0x4202u + (unsigned int)c, g4o[c], (n + GUARD) * sizeof(float));
    }
}

/* ---------------- dot：容差 + n < 8 逐位 + 结果流 ---------------- */

static void t_v3_dot(size_t n, FILE *fd, int check)
{
    const float d = fv_vec3_dot_soa(g3a[0], g3a[1], g3a[2], g3b[0], g3b[1], g3b[2], n);
    double dref = 0.0;
    double scale = 0.0;
    double rel;
    size_t i;

    /* n < 8：整段只走顺序尾巴，必须与标量顺序累加逐位相同（语料里没有 -0.0 的乘积
     * 造成的 +0.0/-0.0 边界，那条边界在契约里已写明属于 dot 例外）。
     * 只在中等量级语料上判定：宽语料的乘积会溢出成 ±inf，那时容差没有意义，
     * 宽语料那一趟只把结果写进流，交给跨形态 cmp 去覆盖确定性。 */
    if (check != 0 && n < 8u) {
        float seq = 0.0f;
        for (i = 0; i < n; ++i) {
            seq = seq + fv_vec3_dot(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]),
                                    fv_vec3_make(g3b[0][i], g3b[1][i], g3b[2][i]));
        }
        chk_bitwise(fbits(d) == fbits(seq), "vec3_dot_soa n<8 bitwise", n);
    }
    for (i = 0; i < n; ++i) {
        const double tx = (double)g3a[0][i] * (double)g3b[0][i];
        const double ty = (double)g3a[1][i] * (double)g3b[1][i];
        const double tz = (double)g3a[2][i] * (double)g3b[2][i];
        dref += tx + ty + tz;
        scale += fabs(tx) + fabs(ty) + fabs(tz);
    }
    rel = fabs((double)d - dref) / (scale > 0.0 ? scale : 1.0);
    if (check == 0) {
        rec(fd, n, 0x4311u, &d, sizeof(d)); /* 宽语料：只进流 */
        ++g_dot_cases;
        return;
    }
    if (rel > g_rel3) {
        g_rel3 = rel;
    }
    chk_bitwise(rel < 1.0e-5, "vec3_dot_soa relative error", n);
    rec(fd, n, 0x4301u, &d, sizeof(d));
    ++g_dot_cases;
}

static void t_v4_dot(size_t n, FILE *fd, int check)
{
    const float d = fv_vec4_dot_soa(g4a[0], g4a[1], g4a[2], g4a[3],
                                    g4b[0], g4b[1], g4b[2], g4b[3], n);
    double dref = 0.0;
    double scale = 0.0;
    double rel;
    size_t i;

    if (check != 0 && n < 8u) {
        float seq = 0.0f;
        for (i = 0; i < n; ++i) {
            seq = seq + fv_vec4_dot(fv_vec4_make(g4a[0][i], g4a[1][i], g4a[2][i], g4a[3][i]),
                                    fv_vec4_make(g4b[0][i], g4b[1][i], g4b[2][i], g4b[3][i]));
        }
        chk_bitwise(fbits(d) == fbits(seq), "vec4_dot_soa n<8 bitwise", n);
    }
    for (i = 0; i < n; ++i) {
        const double tx = (double)g4a[0][i] * (double)g4b[0][i];
        const double ty = (double)g4a[1][i] * (double)g4b[1][i];
        const double tz = (double)g4a[2][i] * (double)g4b[2][i];
        const double tw = (double)g4a[3][i] * (double)g4b[3][i];
        dref += tx + ty + tz + tw;
        scale += fabs(tx) + fabs(ty) + fabs(tz) + fabs(tw);
    }
    rel = fabs((double)d - dref) / (scale > 0.0 ? scale : 1.0);
    if (check == 0) {
        rec(fd, n, 0x4312u, &d, sizeof(d));
        ++g_dot_cases;
        return;
    }
    if (rel > g_rel4) {
        g_rel4 = rel;
    }
    chk_bitwise(rel < 1.0e-5, "vec4_dot_soa relative error", n);
    rec(fd, n, 0x4302u, &d, sizeof(d));
    ++g_dot_cases;
}

/* ---------------- 特殊值 / 边界 ---------------- */

/* n == 0：不读不写任何指针（可以传 NULL）。ASan 下这一组能证明真的没解引用。 */
static void t_null_zero(FILE *fd)
{
    fv_vec4_add_x2(NULL, NULL, NULL, 0);
    fv_vec4_sub_x2(NULL, NULL, NULL, 0);
    fv_vec4_scale_x2(NULL, 0.5f, NULL, 0);
    fv_vec4_lerp_x2(NULL, NULL, 0.5f, NULL, 0);
    fv_vec4_madd_x2(NULL, NULL, 0.5f, NULL, 0);
    fv_vec3_add_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec3_scale_soa(NULL, NULL, NULL, 0.5f, NULL, NULL, NULL, 0);
    fv_vec3_len_soa(NULL, NULL, NULL, NULL, 0);
    fv_vec3_normalize_soa(NULL, NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec4_len_soa(NULL, NULL, NULL, NULL, NULL, 0);
    fv_vec4_normalize_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0);
    chk_bitwise(fbits(fv_vec3_dot_soa(NULL, NULL, NULL, NULL, NULL, NULL, 0)) == 0u,
                "vec3_dot_soa NULL n=0", 0);
    chk_bitwise(fbits(fv_vec4_dot_soa(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0)) == 0u,
                "vec4_dot_soa NULL n=0", 0);
    (void)fd;
}

/* 特殊值：±0 / ±inf / 最小次正规 / FLT_MAX / 单侧 NaN。
 * 单侧 NaN 是严格位级比较（只有一个 NaN 操作数时 payload 与顺序无关）；
 * 两侧都是 NaN 时只要求「同为 NaN」——add/mul 是可交换指令，GCC 可以任选操作数在前，
 * x86 双 NaN 返回第一个操作数的 payload，这不在契约承诺里（tests/test_vec4.c 有同样说明）。 */
static void t_specials(FILE *fm)
{
    static const unsigned int kA[] = { 0x00000000u, 0x80000000u, 0x7F800000u, 0xFF800000u,
                                       0x00000001u, 0x7F7FFFFFu, 0x7FC12345u };
    static const unsigned int kB[] = { 0x80000000u, 0x3F800000u, 0xFF800000u, 0x7F800000u,
                                       0x00800000u, 0xBF800000u, 0x00000000u };
    const size_t n = (size_t)(sizeof(kA) / sizeof(kA[0]));
    size_t i;
    int both_nan_ok;

    for (i = 0; i < n; ++i) {
        g_a4[i] = fv_vec4_make(from_bits(kA[i]), from_bits(kB[i]), from_bits(kA[i]), 1.0f);
        g_b4[i] = fv_vec4_make(from_bits(kB[i]), from_bits(kA[i]), 1.0f, from_bits(kB[i]));
    }
    poison(g_o4, sizeof(g_o4));
    poison(g_r4, sizeof(g_r4));
    fv_vec4_add_x2(g_a4, g_b4, g_o4, n);
    both_nan_ok = 1;
    for (i = 0; i < n; ++i) {
        g_r4[i] = fv_vec4_add(g_a4[i], g_b4[i]);
    }
    /* 逐 lane：位级相同，或两侧都是 NaN（payload 不承诺） */
    {
        size_t k;
        int lane;
        for (k = 0; k < n && both_nan_ok != 0; ++k) {
            for (lane = 0; lane < 4; ++lane) {
                float got;
                float want;
                switch (lane) {
                case 0: got = g_o4[k].x; want = g_r4[k].x; break;
                case 1: got = g_o4[k].y; want = g_r4[k].y; break;
                case 2: got = g_o4[k].z; want = g_r4[k].z; break;
                default: got = g_o4[k].w; want = g_r4[k].w; break;
                }
                if (fbits(got) != fbits(want)
                    && !(isnan(got) != 0 && isnan(want) != 0)) {
                    both_nan_ok = 0;
                }
            }
        }
    }
    chk_bitwise(both_nan_ok, "vec4_add_x2 specials", n);
    rec(fm, n, 0x4401u, g_o4, (n + GUARD) * sizeof(fv_vec4));

    /* lerp / madd 的 NaN 路径也过一遍（同样只对「同为 NaN」放宽） */
    poison(g_o4, sizeof(g_o4));
    fv_vec4_madd_x2(g_a4, g_b4, 0.5f, g_o4, n);
    {
        size_t k;
        int ok = 1;
        for (k = 0; k < n; ++k) {
            const fv_vec4 want = fv_vec4_madd(g_a4[k], g_b4[k], 0.5f);
            if ((fbits(g_o4[k].x) != fbits(want.x)
                 && !(isnan(g_o4[k].x) != 0 && isnan(want.x) != 0))
                || (fbits(g_o4[k].y) != fbits(want.y)
                    && !(isnan(g_o4[k].y) != 0 && isnan(want.y) != 0))
                || (fbits(g_o4[k].z) != fbits(want.z)
                    && !(isnan(g_o4[k].z) != 0 && isnan(want.z) != 0))
                || (fbits(g_o4[k].w) != fbits(want.w)
                    && !(isnan(g_o4[k].w) != 0 && isnan(want.w) != 0))) {
                ok = 0;
            }
        }
        chk_bitwise(ok, "vec4_madd_x2 specials", n);
    }
    rec(fm, n, 0x4402u, g_o4, (n + GUARD) * sizeof(fv_vec4));
}

/* 归一化的零向量保护：全 0 / 极小 / 极大 / 单侧 NaN 的输入，输出必须与标量一致，
 * 且零向量 lane 输出精确的 +0.0（不是 NaN、不是 inf）。 */
static void t_normalize_edges(FILE *fm)
{
    static const float kMag[] = { 0.0f, 5.0e-7f, -5.0e-7f, 1.0f, 1.0e30f, -1.0e30f };
    const size_t n = (size_t)(sizeof(kMag) / sizeof(kMag[0]));
    size_t i;
    int c;

    for (i = 0; i < n; ++i) {
        for (c = 0; c < 4; ++c) {
            g4a[c][i] = kMag[i] * (float)(c + 1);
            g4b[c][i] = 0.25f;
            g3a[c < 3 ? c : 0][i] = kMag[i];
            g3b[c < 3 ? c : 0][i] = 0.25f;
        }
    }
    for (c = 0; c < 3; ++c) {
        poison(g3o[c], sizeof(g3o[c]));
    }
    fv_vec3_normalize_soa(g3a[0], g3a[1], g3a[2], g3o[0], g3o[1], g3o[2], n);
    {
        int ok = 1;
        for (i = 0; i < n; ++i) {
            const fv_vec3 r = fv_vec3_normalize(fv_vec3_make(g3a[0][i], g3a[1][i], g3a[2][i]));
            if (fbits(g3o[0][i]) != fbits(r.x) || fbits(g3o[1][i]) != fbits(r.y)
                || fbits(g3o[2][i]) != fbits(r.z)) {
                ok = 0;
            }
        }
        chk_bitwise(ok, "vec3_normalize_soa edges", n);
    }
    for (c = 0; c < 4; ++c) {
        poison(g4o[c], sizeof(g4o[c]));
    }
    fv_vec4_normalize_soa(g4a[0], g4a[1], g4a[2], g4a[3],
                          g4o[0], g4o[1], g4o[2], g4o[3], n);
    {
        int ok = 1;
        for (i = 0; i < n; ++i) {
            const fv_vec4 r = fv_vec4_normalize(
                fv_vec4_make(g4a[0][i], g4a[1][i], g4a[2][i], g4a[3][i]));
            if (fbits(g4o[0][i]) != fbits(r.x) || fbits(g4o[1][i]) != fbits(r.y)
                || fbits(g4o[2][i]) != fbits(r.z) || fbits(g4o[3][i]) != fbits(r.w)) {
                ok = 0;
            }
        }
        chk_bitwise(ok, "vec4_normalize_soa edges", n);
        /* 零向量那一行必须是精确的 +0.0 */
        chk_bitwise(fbits(g4o[0][0]) == 0u && fbits(g4o[1][0]) == 0u
                        && fbits(g4o[2][0]) == 0u && fbits(g4o[3][0]) == 0u,
                    "vec4_normalize_soa zero -> +0.0", n);
    }
    for (c = 0; c < 3; ++c) {
        rec(fm, n, 0x4501u + (unsigned int)c, g3o[c], (n + GUARD) * sizeof(float));
    }
    for (c = 0; c < 4; ++c) {
        rec(fm, n, 0x4502u + (unsigned int)c, g4o[c], (n + GUARD) * sizeof(float));
    }
}

int main(int argc, char **argv)
{
    static const size_t kLens[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                    16, 17, 63, 1024, 8192 };
    static const float kParams[] = { 0.0f, -0.0f, 1.0f, -1.0f, 1.0e-6f, -3.25f };
    const int nlens = (int)(sizeof(kLens) / sizeof(kLens[0]));
    const int nparams = (int)(sizeof(kParams) / sizeof(kParams[0]));
    FILE *fm;
    FILE *fd;
    unsigned int cases_main;
    int li;
    int pi;
    size_t i;

    if (argc < 3) {
        fprintf(stderr, "usage: %s <main.bin> <dot.bin>\n", argv[0]);
        return 2;
    }
    fm = fopen(argv[1], "wb");
    if (fm == NULL) {
        perror(argv[1]);
        return 2;
    }
    fd = fopen(argv[2], "wb");
    if (fd == NULL) {
        perror(argv[2]);
        fclose(fm);
        return 2;
    }

    t_null_zero(fd);

    for (li = 0; li < nlens; ++li) {
        const size_t n = kLens[li];

        fill_aos4(g_a4, n, 0x500000u + (unsigned int)li);
        fill_aos4(g_b4, n, 0x600000u + (unsigned int)li);
        fill_soa(g3a, n, 0x700000u + (unsigned int)li, 3);
        fill_soa(g3b, n, 0x710000u + (unsigned int)li, 3);
        fill_soa(g4a, n, 0x720000u + (unsigned int)li, 4);
        fill_soa(g4b, n, 0x730000u + (unsigned int)li, 4);

        t_v4_add(n, fm);
        t_v4_sub(n, fm);
        t_v4_lerp(n, 0.25f, fm);
        for (pi = 0; pi < nparams; ++pi) {
            t_v4_scale(n, kParams[pi], fm);
            t_v4_madd(n, kParams[pi], fm);
        }

        t_v3_add(n, fm);
        t_v3_len(n, fm);
        t_v3_normalize(n, fm);
        for (pi = 0; pi < nparams; ++pi) {
            t_v3_scale(n, kParams[pi], fm);
        }

        t_v4_len(n, fm);
        t_v4_normalize(n, fm);

        /* 宽语料：只把结果写进流（跨形态 cmp 覆盖确定性），不做容差判定 */
        t_v3_dot(n, fd, 0);
        t_v4_dot(n, fd, 0);

        /* 中等语料：容差判定 + n < 8 的逐位判定 */
        fill_moderate(g3a, g3b, n, 0x740000u + (unsigned int)li, 0x750000u + (unsigned int)li, 3);
        fill_moderate(g4a, g4b, n, 0x760000u + (unsigned int)li, 0x770000u + (unsigned int)li, 4);
        t_v3_dot(n, fd, 1);
        t_v4_dot(n, fd, 1);
    }

    /* 尾长覆盖之外的补充：8192 已经在上面的 kLens 里（主循环，不是单独特例），
     * 这里补两组「输出指针相位偏移」的用例 —— _x2 的 storeu 与 _soa 的 8 lane 主体
     * 都不要求对齐，但相位变化会改变尾块的起点，是尾部逻辑最容易错的地方。 */
    {
        size_t k;
        for (k = 1; k <= 2; ++k) {
            fv_vec4 *o = g_o4 + k;
            fv_vec4 *r = g_r4 + k;
            const size_t n = 1023u;
            fill_aos4(g_a4, n, 0x800000u + (unsigned int)k);
            fill_aos4(g_b4, n, 0x810000u + (unsigned int)k);
            poison(g_o4, sizeof(g_o4));
            poison(g_r4, sizeof(g_r4));
            fv_vec4_add_x2(g_a4, g_b4, o, n);
            for (i = 0; i < n; ++i) {
                r[i] = fv_vec4_add(g_a4[i], g_b4[i]);
            }
            chk_bitwise(memcmp(o, r, (n + GUARD) * sizeof(fv_vec4)) == 0,
                        "vec4_add_x2 phase offset", n);
            rec(fm, n, 0x4601u ^ (unsigned int)k, o, (n + GUARD) * sizeof(fv_vec4));
        }
    }

    t_specials(fm);
    t_normalize_edges(fm);

    cases_main = g_case;
    fclose(fm);
    fclose(fd);
    printf("parity34: cases=%u dot_cases=%u failures=%u  max_rel v3=%.3e v4=%.3e\n",
           cases_main - g_dot_cases, g_dot_cases, g_fail, g_rel3, g_rel4);
    return g_fail == 0u ? 0 : 1;
}
