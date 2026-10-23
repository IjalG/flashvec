/* 位级对拍：默认（inline）形态 vs FV_USE_EXTERN（库内实现）形态 —— 验收项 6。
 *
 * 做法：同一个源文件编四份（C / C++ × inline / extern），每份把所有用例的输出按固定格式
 * 写成二进制记录流，然后逐字节 cmp。同一个源文件保证四份的输入、调用顺序、参数完全一致。
 *
 * 记录格式（本机字节序，四个 u32 头 + 原始字节）：
 *   [case_id][n][param_bits][nbytes] 后跟 payload。
 * payload 从输出缓冲区第 0 个元素开始，长度是 (n + GUARD) 个元素 —— 多出来的 GUARD 个元素
 * 调用前被 memset 成 0x5A，所以「越界写」与「该写没写」都会让两份文件不同。
 *
 * 覆盖：n = 0..17 全尾长 + 一个 63 的长批（跨 x4 的 4、soa 的 8、dot 的 32 三种宽度）；
 *       scale/lerp/madd 各 5 个参数（含 0.0、-0.0、1.0、-1.0、1e-6、-3.25）；
 *       mask 四种模式（全 0 / 全非 0 / 交替 / 伪随机）；
 *       idx 四种模式（重复、逆序、乱序、单元素）+ m == 0。
 *
 * fv_vec2_dot_soa 单独写一个流：它是契约里唯一不保证逐位一致的入口（归约顺序不同），
 * 所以只在「同一归约路径」之间比 —— 见 Makefile 里的 cmp 分组。
 *
 * 用法：test_parity <main.bin> <dot.bin>
 */
#include "flashvec.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CAP 1032
#define GUARD 4u

static fv_vec2 g_a[CAP + GUARD], g_b[CAP + GUARD], g_o[CAP + GUARD];
static float g_ax[CAP + GUARD], g_ay[CAP + GUARD], g_bx[CAP + GUARD], g_by[CAP + GUARD];
static float g_ox[CAP + GUARD], g_oy[CAP + GUARD];
static unsigned char g_mask[CAP + GUARD];
static unsigned int g_idx[CAP + GUARD];

/* P2 新增入口用的第二组缓冲区：一块放旧入口的结果、一块放新入口的结果，逐字节比较。 */
static fv_vec2 g_s1[CAP + GUARD], g_s2[CAP + GUARD];
static float g_sx[CAP + GUARD], g_sy[CAP + GUARD];   /* 旧入口的 SoA 输出 */
static float g_sdx[CAP + GUARD], g_sdy[CAP + GUARD]; /* 新入口的 SoA 输出 */

/* 大 n 的抽查（CAP 之外）：只用来验证 peel / NT 主体 / 尾巴三段在真实长度下的接线。 */
#define CAPL 5200u
static fv_vec2 g_la[CAPL + GUARD], g_lb[CAPL + GUARD], g_lo1[CAPL + GUARD], g_lo2[CAPL + GUARD];
static float g_lax[CAPL + GUARD], g_lay[CAPL + GUARD], g_lbx[CAPL + GUARD], g_lby[CAPL + GUARD];
static float g_lox1[CAPL + GUARD], g_loy1[CAPL + GUARD], g_lox2[CAPL + GUARD], g_loy2[CAPL + GUARD];
static unsigned char g_lmask[CAPL + GUARD];

static unsigned int g_fail;

static void chk_bitwise(int ok, const char *what, size_t n)
{
    if (ok == 0) {
        ++g_fail;
        printf("FAIL %-28s n=%u\n", what, (unsigned int)n);
    }
}

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

/* 中等量级、会互相抵消的随机值：对归约顺序最敏感，专门用来暴露 dot 的路径差异。 */
static float rnd_moderate(void)
{
    return (float)((int)(rnd_u32() % 20001u) - 10000) / 1000.0f;
}

/* 浮点参数按位模式记录：这样两种形态记录的 param 字段也完全一致。 */
static unsigned int fbits(float f)
{
    unsigned int b;
    memcpy(&b, &f, 4);
    return b;
}

static const float kInteresting[] = { 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f,
                                      1.0e-6f, 1.0e-12f, 4.9e-40f, -1.0e-30f, 3.0e38f };
#define NINTEREST ((int)(sizeof(kInteresting) / sizeof(kInteresting[0])))

static void fill_vec2(fv_vec2 *v, size_t n, unsigned int seed)
{
    size_t i;
    g_rng = seed;
    for (i = 0; i < n; ++i) {
        v[i].x = (i % 3u == 0u) ? kInteresting[i % (size_t)NINTEREST] : rnd_float();
        v[i].y = (i % 3u == 1u) ? kInteresting[(i + 5u) % (size_t)NINTEREST] : rnd_float();
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

static unsigned int g_case;
static unsigned int g_dot_cases;

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

static void poison(void *p, size_t bytes)
{
    memset(p, 0x5A, bytes);
}

int main(int argc, char **argv)
{
    static const size_t kLens[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 63 };
    const int nlens = (int)(sizeof(kLens) / sizeof(kLens[0]));
    static const float kParams[] = { 0.0f, -0.0f, 1.0f, -1.0f, 1.0e-6f, -3.25f };
    const int nparams = (int)(sizeof(kParams) / sizeof(kParams[0]));
    FILE *fm;
    FILE *fd;
    unsigned int cases_main;
    int li;
    int pi;
    int mode;
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

    for (li = 0; li < nlens; ++li) {
        const size_t n = kLens[li];

        fill_vec2(g_a, n, 0x100000u + (unsigned int)li);
        fill_vec2(g_b, n, 0x200000u + (unsigned int)li);
        to_soa(g_a, g_ax, g_ay, n);
        to_soa(g_b, g_bx, g_by, n);

        /* x4：add / sub */
        poison(g_o, sizeof(g_o));
        fv_vec2_add_x4(g_a, g_b, g_o, n);
        rec(fm, n, 0u, g_o, (n + GUARD) * sizeof(fv_vec2));

        poison(g_o, sizeof(g_o));
        fv_vec2_sub_x4(g_a, g_b, g_o, n);
        rec(fm, n, 0u, g_o, (n + GUARD) * sizeof(fv_vec2));

        /* x4：scale / lerp / madd（参数逐个） */
        for (pi = 0; pi < nparams; ++pi) {
            const float s = kParams[pi];
            poison(g_o, sizeof(g_o));
            fv_vec2_scale_x4(g_a, s, g_o, n);
            rec(fm, n, fbits(s), g_o, (n + GUARD) * sizeof(fv_vec2));

            poison(g_o, sizeof(g_o));
            fv_vec2_lerp_x4(g_a, g_b, s, g_o, n);
            rec(fm, n, fbits(s), g_o, (n + GUARD) * sizeof(fv_vec2));

            poison(g_o, sizeof(g_o));
            fv_vec2_madd_x4(g_a, g_b, s, g_o, n);
            rec(fm, n, fbits(s), g_o, (n + GUARD) * sizeof(fv_vec2));
        }

        /* mask：四种模式 */
        for (mode = 0; mode < 4; ++mode) {
            for (i = 0; i < n; ++i) {
                switch (mode) {
                case 0: g_mask[i] = 0u; break;                          /* 全 0 */
                case 1: g_mask[i] = 0xFFu; break;                       /* 全非 0 */
                case 2: g_mask[i] = (unsigned char)(i % 2u); break;      /* 交替 */
                default: g_mask[i] = (unsigned char)(rnd_u32() & 1u); break; /* 伪随机 */
                }
            }
            poison(g_o, sizeof(g_o));
            fv_vec2_add_x4_mask(g_a, g_b, g_o, g_mask, n);
            rec(fm, n, (unsigned int)mode, g_o, (n + GUARD) * sizeof(fv_vec2));
        }

        /* idx：重复 / 逆序 / 乱序 / 单元素 / m == 0 */
        for (mode = 0; mode < 5; ++mode) {
            size_t m;
            switch (mode) {
            case 0: /* 重复 + 少量下标 */
                m = (n > 4u) ? 5u : n;
                for (i = 0; i < m; ++i) g_idx[i] = (unsigned int)((i * 3u) % (n == 0u ? 1u : n));
                if (m > 1u) g_idx[m - 1u] = g_idx[0];
                break;
            case 1: /* 逆序 */
                m = n;
                for (i = 0; i < m; ++i) g_idx[i] = (unsigned int)(m - 1u - i);
                break;
            case 2: /* 乱序（步长与 n 互质时是排列） */
                m = n;
                for (i = 0; i < m; ++i) g_idx[i] = (unsigned int)((i * 7u + 5u) % (n == 0u ? 1u : n));
                break;
            case 3: /* 单元素 */
                m = (n > 0u) ? 1u : 0u;
                if (m != 0u) g_idx[0] = (unsigned int)(n - 1u);
                break;
            default: /* m == 0：什么都不写 */
                m = 0u;
                break;
            }
            poison(g_o, sizeof(g_o));
            if (n != 0u) {
                memcpy(g_o, g_a, n * sizeof(fv_vec2));
            }
            fv_vec2_add_x4_idx(g_a, g_b, g_o, g_idx, m);
            rec(fm, n, (unsigned int)mode, g_o, (n + GUARD) * sizeof(fv_vec2));

            poison(g_o, sizeof(g_o));
            if (n != 0u) {
                memcpy(g_o, g_a, n * sizeof(fv_vec2));
            }
            fv_vec2_scale_x4_idx(g_a, 0.5f, g_o, g_idx, m);
            rec(fm, n, (unsigned int)mode, g_o, (n + GUARD) * sizeof(fv_vec2));
        }

        /* SoA：add / scale / len / normalize */
        poison(g_ox, sizeof(g_ox));
        poison(g_oy, sizeof(g_oy));
        fv_vec2_add_soa(g_ax, g_ay, g_bx, g_by, g_ox, g_oy, n);
        rec(fm, n, 0u, g_ox, (n + GUARD) * sizeof(float));
        rec(fm, n, 1u, g_oy, (n + GUARD) * sizeof(float));

        for (pi = 0; pi < nparams; ++pi) {
            const float s = kParams[pi];
            poison(g_ox, sizeof(g_ox));
            poison(g_oy, sizeof(g_oy));
            fv_vec2_scale_soa(g_ax, g_ay, s, g_ox, g_oy, n);
            rec(fm, n, fbits(s), g_ox, (n + GUARD) * sizeof(float));
            rec(fm, n, fbits(s) ^ 0xFFFFFFFFu, g_oy, (n + GUARD) * sizeof(float));
        }

        poison(g_ox, sizeof(g_ox));
        fv_vec2_len_soa(g_ax, g_ay, g_ox, n);
        rec(fm, n, 0u, g_ox, (n + GUARD) * sizeof(float));

        poison(g_ox, sizeof(g_ox));
        poison(g_oy, sizeof(g_oy));
        fv_vec2_normalize_soa(g_ax, g_ay, g_ox, g_oy, n);
        rec(fm, n, 0u, g_ox, (n + GUARD) * sizeof(float));
        rec(fm, n, 1u, g_oy, (n + GUARD) * sizeof(float));

        /* ---- P2：_stream 与 _mask_inplace ----
         * 断言分两层：
         *   1) 进程内逐位比较：新入口的输出与对应旧入口（或标量 API）memcmp 相等，
         *      比较长度含 GUARD，所以「越界写」与「该写没写」都会被抓到；
         *   2) 把新入口的输出也记进主流，让跨形态（C/C++ × inline/extern × 纯 C 回退 ×
         *      FMA 消费者）的 cmp 一起覆盖它们。 */
        {
            /* x4：add / sub / scale（参考值当场重算：g_o 此刻装的是上面 idx 变体的输出） */
            poison(g_s1, sizeof(g_s1));
            poison(g_s2, sizeof(g_s2));
            fv_vec2_add_x4(g_a, g_b, g_s2, n);
            fv_vec2_add_x4_stream(g_a, g_b, g_s1, n);
            chk_bitwise(memcmp(g_s2, g_s1, (n + GUARD) * sizeof(fv_vec2)) == 0,
                        "add_x4_stream == add_x4", n);
            rec(fm, n, 0x5301u, g_s1, (n + GUARD) * sizeof(fv_vec2));

            poison(g_s1, sizeof(g_s1));
            poison(g_s2, sizeof(g_s2));
            fv_vec2_sub_x4(g_a, g_b, g_s2, n);
            fv_vec2_sub_x4_stream(g_a, g_b, g_s1, n);
            chk_bitwise(memcmp(g_s2, g_s1, (n + GUARD) * sizeof(fv_vec2)) == 0,
                        "sub_x4_stream == sub_x4", n);
            rec(fm, n, 0x5302u, g_s1, (n + GUARD) * sizeof(fv_vec2));

            for (pi = 0; pi < nparams; ++pi) {
                const float s = kParams[pi];
                poison(g_s1, sizeof(g_s1));
                poison(g_s2, sizeof(g_s2));
                fv_vec2_scale_x4(g_a, s, g_s2, n);
                fv_vec2_scale_x4_stream(g_a, s, g_s1, n);
                chk_bitwise(memcmp(g_s2, g_s1, (n + GUARD) * sizeof(fv_vec2)) == 0,
                            "scale_x4_stream == scale_x4", n);
                rec(fm, n, 0x5303u ^ fbits(s), g_s1, (n + GUARD) * sizeof(fv_vec2));
            }

            /* SoA：add / scale / normalize（新旧入口各写一次，参考值当场重算，
               不能复用上面 g_ox/g_oy —— 它们此刻装的是 normalize 的结果） */
            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            poison(g_sdx, sizeof(g_sdx));
            poison(g_sdy, sizeof(g_sdy));
            fv_vec2_add_soa(g_ax, g_ay, g_bx, g_by, g_sx, g_sy, n);
            fv_vec2_add_soa_stream(g_ax, g_ay, g_bx, g_by, g_sdx, g_sdy, n);
            chk_bitwise(memcmp(g_sx, g_sdx, (n + GUARD) * sizeof(float)) == 0
                            && memcmp(g_sy, g_sdy, (n + GUARD) * sizeof(float)) == 0,
                        "add_soa_stream == add_soa", n);
            rec(fm, n, 0x5304u, g_sdx, (n + GUARD) * sizeof(float));
            rec(fm, n, 0x5304u, g_sdy, (n + GUARD) * sizeof(float));

            for (pi = 0; pi < nparams; ++pi) {
                const float s = kParams[pi];
                poison(g_sx, sizeof(g_sx));
                poison(g_sy, sizeof(g_sy));
                poison(g_sdx, sizeof(g_sdx));
                poison(g_sdy, sizeof(g_sdy));
                fv_vec2_scale_soa(g_ax, g_ay, s, g_sx, g_sy, n);
                fv_vec2_scale_soa_stream(g_ax, g_ay, s, g_sdx, g_sdy, n);
                chk_bitwise(memcmp(g_sx, g_sdx, (n + GUARD) * sizeof(float)) == 0
                                && memcmp(g_sy, g_sdy, (n + GUARD) * sizeof(float)) == 0,
                            "scale_soa_stream == scale_soa", n);
                rec(fm, n, 0x5305u ^ fbits(s), g_sdx, (n + GUARD) * sizeof(float));
                rec(fm, n, 0x5305u ^ fbits(s), g_sdy, (n + GUARD) * sizeof(float));
            }

            poison(g_sx, sizeof(g_sx));
            poison(g_sy, sizeof(g_sy));
            poison(g_sdx, sizeof(g_sdx));
            poison(g_sdy, sizeof(g_sdy));
            fv_vec2_normalize_soa(g_ax, g_ay, g_sx, g_sy, n);
            fv_vec2_normalize_soa_stream(g_ax, g_ay, g_sdx, g_sdy, n);
            chk_bitwise(memcmp(g_sx, g_sdx, (n + GUARD) * sizeof(float)) == 0
                            && memcmp(g_sy, g_sdy, (n + GUARD) * sizeof(float)) == 0,
                        "normalize_soa_stream == normalize_soa", n);
            rec(fm, n, 0x5306u, g_sdx, (n + GUARD) * sizeof(float));
            rec(fm, n, 0x5306u, g_sdy, (n + GUARD) * sizeof(float));

            /* _mask 的块级预扫：掩码字节只取 0x02 / 0x80 / 0xFE 这类「非零但最低位为 0
             * 或只有高位」的值。任何「按位与 0x01010101 判零」的写法都会把它们误判成
             * 「没命中」（0x02 & 0x01 == 0），所以这里直接对标量参考值比较，不能只靠
             * 跨形态 cmp（那只能证明各形态错得一样）。 */
            for (mode = 0; mode < 2; ++mode) {
                for (i = 0; i < n; ++i) {
                    static const unsigned char kVals[4] = { 0x02u, 0x80u, 0xFEu, 0x00u };
                    g_mask[i] = kVals[(i + (size_t)mode) % 4u];
                }
                poison(g_s1, sizeof(g_s1));
                poison(g_s2, sizeof(g_s2));
                fv_vec2_add_x4_mask(g_a, g_b, g_s1, g_mask, n);
                for (i = 0; i < n; ++i) {
                    g_s2[i] = (g_mask[i] != 0u) ? fv_vec2_add(g_a[i], g_b[i]) : g_a[i];
                }
                chk_bitwise(memcmp(g_s1, g_s2, n * sizeof(fv_vec2)) == 0,
                            "add_x4_mask vs scalar (odd bytes)", n);
                rec(fm, n, 0x5309u ^ (unsigned int)mode, g_s1, (n + GUARD) * sizeof(fv_vec2));
            }

            /* _mask_inplace：语义是「命中才就地写」，所以参考值按标量语义构造：
             * 命中 → fv_vec2_add(a[i], b[i])，未命中 → a[i] 原值（不是 0x5A，也不是不初始化）。 */
            for (mode = 0; mode < 4; ++mode) {
                for (i = 0; i < n; ++i) {
                    switch (mode) {
                    case 0: g_mask[i] = 0u; break;
                    case 1: g_mask[i] = 0xFFu; break;
                    case 2: g_mask[i] = (unsigned char)(i % 2u); break;
                    default: g_mask[i] = (unsigned char)(rnd_u32() % 3u); break;
                    }
                }
                poison(g_s1, sizeof(g_s1));
                poison(g_s2, sizeof(g_s2));
                if (n != 0u) {
                    memcpy(g_s1, g_a, n * sizeof(fv_vec2));
                    memcpy(g_s2, g_a, n * sizeof(fv_vec2));
                }
                for (i = 0; i < n; ++i) {
                    if (g_mask[i] != 0u) {
                        g_s1[i] = fv_vec2_add(g_a[i], g_b[i]);
                    }
                }
                fv_vec2_add_x4_mask_inplace(g_s2, g_b, g_mask, n);
                chk_bitwise(memcmp(g_s1, g_s2, (n + GUARD) * sizeof(fv_vec2)) == 0,
                            "add_x4_mask_inplace semantics", n);
                rec(fm, n, 0x5307u ^ (unsigned int)mode, g_s2, (n + GUARD) * sizeof(fv_vec2));
            }
        }

        /* dot 单独一个流（契约里唯一不保证逐位一致的入口） */
        {
            const float d = fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, n);
            rec(fd, n, 0u, &d, sizeof(d));
            ++g_dot_cases;
        }
    }

    cases_main = g_case;

    /* dot 的敏感语料：中等量级、会互相抵消的随机值 + 更大的 n。
     * 归约顺序不同的两条路径在这里通常就会分出 1 ulp —— 这正是契约里写明的那条例外。 */
    {
        static const size_t kDotLens[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 63, 1000 };
        const int ndot = (int)(sizeof(kDotLens) / sizeof(kDotLens[0]));
        int di;
        for (di = 0; di < ndot; ++di) {
            const size_t nd = kDotLens[di];
            float d;
            g_rng = 0xE0000000u + (unsigned int)di;
            for (i = 0; i < nd; ++i) {
                g_ax[i] = rnd_moderate();
                g_ay[i] = rnd_moderate();
                g_bx[i] = rnd_moderate();
                g_by[i] = rnd_moderate();
            }
            d = fv_vec2_dot_soa(g_ax, g_ay, g_bx, g_by, nd);
            rec(fd, nd, 1u, &d, sizeof(d));
            ++g_dot_cases;
        }
    }
    /* ---- P2：大 n 抽查（CAP 之外的真实长度；含非 8 的倍数，逼出尾巴） ----
     * 覆盖 peel 的各种相位：4096（32 B 对齐、整向量）、4103（偏移 7 元素）、5039（素数长度）。 */
    {
        static const size_t kBig[] = { 4096u, 4103u, 5039u };
        const int nbig = (int)(sizeof(kBig) / sizeof(kBig[0]));
        int bi;
        for (bi = 0; bi < nbig; ++bi) {
            const size_t nb = kBig[bi];
            int off;
            fill_vec2(g_la, nb, 0x300000u + (unsigned int)bi);
            fill_vec2(g_lb, nb, 0x400000u + (unsigned int)bi);
            to_soa(g_la, g_lax, g_lay, nb);
            to_soa(g_lb, g_lbx, g_lby, nb);
            /* 输出指针逐个相位平移（fv_vec2 侧 0..3 元素，SoA 侧 0..7 个 float）：
             * 这样 peel 的每一个分支都被真实执行到，而不是只在「恰好对齐」的情况下过一遍。 */
            for (off = 0; off < 4; ++off) {
                fv_vec2 *o1 = g_lo1 + off;
                fv_vec2 *o2 = g_lo2 + off;
                const size_t bytes = (nb + GUARD - (size_t)off) * sizeof(fv_vec2);

                poison(g_lo1, sizeof(g_lo1));
                poison(g_lo2, sizeof(g_lo2));
                fv_vec2_add_x4(g_la, g_lb, o1, nb);
                fv_vec2_add_x4_stream(g_la, g_lb, o2, nb);
                chk_bitwise(memcmp(o1, o2, bytes) == 0, "big add_x4_stream", nb);

                poison(g_lo1, sizeof(g_lo1));
                poison(g_lo2, sizeof(g_lo2));
                fv_vec2_sub_x4(g_la, g_lb, o1, nb);
                fv_vec2_sub_x4_stream(g_la, g_lb, o2, nb);
                chk_bitwise(memcmp(o1, o2, bytes) == 0, "big sub_x4_stream", nb);

                poison(g_lo1, sizeof(g_lo1));
                poison(g_lo2, sizeof(g_lo2));
                fv_vec2_scale_x4(g_la, -3.25f, o1, nb);
                fv_vec2_scale_x4_stream(g_la, -3.25f, o2, nb);
                chk_bitwise(memcmp(o1, o2, bytes) == 0, "big scale_x4_stream", nb);
            }
            for (off = 0; off < 8; ++off) {
                float *x1 = g_lox1 + off;
                float *y1 = g_loy1 + off;
                float *x2 = g_lox2 + off;
                float *y2 = g_loy2 + off;
                const size_t bytes = (nb + GUARD - (size_t)off) * sizeof(float);

                poison(g_lox1, sizeof(g_lox1));
                poison(g_loy1, sizeof(g_loy1));
                poison(g_lox2, sizeof(g_lox2));
                poison(g_loy2, sizeof(g_loy2));
                fv_vec2_add_soa(g_lax, g_lay, g_lbx, g_lby, x1, y1, nb);
                fv_vec2_add_soa_stream(g_lax, g_lay, g_lbx, g_lby, x2, y2, nb);
                chk_bitwise(memcmp(x1, x2, bytes) == 0 && memcmp(y1, y2, bytes) == 0,
                            "big add_soa_stream", nb);

                poison(g_lox1, sizeof(g_lox1));
                poison(g_loy1, sizeof(g_loy1));
                poison(g_lox2, sizeof(g_lox2));
                poison(g_loy2, sizeof(g_loy2));
                fv_vec2_scale_soa(g_lax, g_lay, 1.0e-6f, x1, y1, nb);
                fv_vec2_scale_soa_stream(g_lax, g_lay, 1.0e-6f, x2, y2, nb);
                chk_bitwise(memcmp(x1, x2, bytes) == 0 && memcmp(y1, y2, bytes) == 0,
                            "big scale_soa_stream", nb);

                poison(g_lox1, sizeof(g_lox1));
                poison(g_loy1, sizeof(g_loy1));
                poison(g_lox2, sizeof(g_lox2));
                poison(g_loy2, sizeof(g_loy2));
                fv_vec2_normalize_soa(g_lax, g_lay, x1, y1, nb);
                fv_vec2_normalize_soa_stream(g_lax, g_lay, x2, y2, nb);
                chk_bitwise(memcmp(x1, x2, bytes) == 0 && memcmp(y1, y2, bytes) == 0,
                            "big normalize_soa_stream", nb);
            }
            /* 就地稀疏加：5% 命中（<=1/8，走稀疏支 = 块级预扫 + 标量）
             * 与 50% 命中（>1/8，走密集支 = blend + storeu）—— 两段都被真实执行到 */
            {
                int d;
                for (d = 0; d < 2; ++d) {
                    const unsigned int step = (d == 0) ? 20u : 2u;
                    size_t k;
                    for (k = 0; k < nb; ++k) {
                        g_lmask[k] = (unsigned char)((k % (size_t)step == 0u) ? 1u : 0u);
                    }
                    poison(g_lo1, sizeof(g_lo1));
                    poison(g_lo2, sizeof(g_lo2));
                    memcpy(g_lo1, g_la, nb * sizeof(fv_vec2));
                    memcpy(g_lo2, g_la, nb * sizeof(fv_vec2));
                    for (k = 0; k < nb; ++k) {
                        if (g_lmask[k] != 0u) {
                            g_lo1[k] = fv_vec2_add(g_la[k], g_lb[k]);
                        }
                    }
                    fv_vec2_add_x4_mask_inplace(g_lo2, g_lb, g_lmask, nb);
                    chk_bitwise(memcmp(g_lo1, g_lo2, (nb + GUARD) * sizeof(fv_vec2)) == 0,
                                "big add_x4_mask_inplace", nb);
                }
            }
        }
    }

    fclose(fm);
    fclose(fd);
    printf("parity: main_cases=%u dot_cases=%u  p2_failures=%u\n", cases_main - g_dot_cases,
           g_dot_cases, g_fail);
    return g_fail == 0u ? 0 : 1;
}