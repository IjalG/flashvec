/* verify/parity.cpp —— flashvec 批量 SIMD 路径的位级一致性验证（验收者 F 独立实现）
 *
 * 做三件事：
 *   1) 每个逐元素型批量函数 vs (a) 库自己的标量 API、(b) 我自己的独立标量参考（ref_strict.cpp，
 *      -ffp-contract=off，不会被收缩成 FMA）逐元素做位模式比较；
 *   2) 覆盖全部尾长 n = 0..17、多种数据类别、掩码/索引语义、越界写入 canary；
 *   3) P0 小改动（min/max/clamp/madd/aabb2）与常量 FV_EPS_SQR 的检查。
 *
 * 比较策略：默认按 IEEE 位模式 memcmp 语义；唯一例外是 NaN —— 若两侧都是 NaN，
 * 即使 payload 不同也判为「相等」（记为 nan_pairs 计数并在输出里说明）。
 * ±0.0 的差异按位模式如实计为不一致，但单独分类为 SIGNED_ZERO。
 * 含 NaN/Inf 的输入数据单独统计（edge 计数），不计入主判定，避免把定义域外行为当失败。
 * 独立 TU，可脱离 perf 单独构建运行：./parity
 */
#include "flashvec.h"
#include "ref_strict.h"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static const size_t MAXN = 64;      /* 数组容量 */
static const size_t MAX_TAIL = 17;  /* 尾长覆盖 0..17 */

static_assert(sizeof(fv_vec2) == 2 * sizeof(float), "fv_vec2 must be 2 packed floats");

/* ------------------------------------------------------------------ 基础工具 */

static uint32_t fbits(float f) { uint32_t u; std::memcpy(&u, &f, 4); return u; }
static float ffrom(uint32_t u) { float f; std::memcpy(&f, &u, 4); return f; }
static bool nan_bits(uint32_t u) { return ((u >> 23) & 0xffu) == 0xffu && (u & 0x7fffffu) != 0u; }
static bool zero_bits(uint32_t u) { return (u << 1) == 0u; }

struct Same { bool equal; bool nan_pair; bool szero; };

static Same same_bits(float x, float y)
{
    Same r;
    r.equal = false; r.nan_pair = false; r.szero = false;
    uint32_t ux = fbits(x), uy = fbits(y);
    if (ux == uy) { r.equal = true; return r; }
    if (nan_bits(ux) && nan_bits(uy)) { r.equal = true; r.nan_pair = true; return r; }
    r.szero = zero_bits(ux) && zero_bits(uy);
    return r;
}

static uint32_t hashu(uint32_t x)
{
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return x;
}

/* ------------------------------------------------------------------ 统计 */

struct Stat {
    const char *label;
    long samples = 0;         /* vs 库标量 API 比较过的 float 标量个数 */
    long samples_strict = 0;  /* vs 独立标量参考 比较过的 float 标量个数 */
    long bad_lib = 0, bad_strict = 0;
    long canary = 0;          /* 越过 n 的写入 */
    long untouched = 0;       /* _idx：o 中未被引用的位置被改写 */
    long edge = 0;            /* 含 NaN/Inf 的输入上的分歧 */
    long nan_pairs = 0;       /* 两侧都是 NaN（payload 可能不同）计为相等 */
    long szero = 0;           /* 仅 ±0 符号位不同 */
    int shown_lib = 0, shown_strict = 0, shown_canary = 0, shown_edge = 0, shown_untouched = 0;
};

static void show_diff(Stat& st, const char* who, size_t i, size_t comps, float got, float want, const char* cls)
{
    size_t elem = comps ? i / comps : i;
    size_t comp = comps ? i % comps : 0;
    printf("        %-22s %-11s elem=%zu %s got=0x%08x (%.9g) want=0x%08x (%.9g) class=%s\n",
           st.label, who, elem, comps == 2 ? (comp == 0 ? ".x" : ".y") : "[f]",
           fbits(got), (double)got, fbits(want), (double)want, cls);
}

/* comps = 1（float 数组）或 2（fv_vec2 数组）；which: 0 = vs 库标量, 1 = vs 独立参考 */
static void cmp_floats(Stat& st, const float* got, const float* want, size_t nfloats, int comps, int which, bool edge)
{
    for (size_t i = 0; i < nfloats; ++i) {
        Same s = same_bits(got[i], want[i]);
        if (which == 0) st.samples++; else st.samples_strict++;
        if (s.nan_pair) st.nan_pairs++;
        if (s.equal) continue;
        const char* cls = s.szero ? "SIGNED_ZERO"
                        : (nan_bits(fbits(got[i])) || nan_bits(fbits(want[i])) ? "NAN_ONE_SIDE" : "VALUE");
        if (s.szero) st.szero++;
        if (edge) {
            st.edge++;
            if (!st.shown_edge) { st.shown_edge = 1; show_diff(st, which ? "strict-ref" : "lib-scalar", i, comps, got[i], want[i], cls); }
            continue;
        }
        if (which == 0) {
            st.bad_lib++;
            if (!st.shown_lib) { st.shown_lib = 1; show_diff(st, "lib-scalar", i, comps, got[i], want[i], cls); }
        } else {
            st.bad_strict++;
            if (!st.shown_strict) { st.shown_strict = 1; show_diff(st, "strict-ref", i, comps, got[i], want[i], cls); }
        }
    }
}

static float sentinel_f(void) { return ffrom(0xCDCDCDCDu); }
static void fill_sentinel(std::vector<fv_vec2>& v) { for (size_t i = 0; i < v.size(); ++i) { v[i].x = sentinel_f(); v[i].y = sentinel_f(); } }
static void fill_sentinel_f(std::vector<float>& v) { for (size_t i = 0; i < v.size(); ++i) v[i] = sentinel_f(); }

static void canary(Stat& st, const std::vector<fv_vec2>& before, const std::vector<fv_vec2>& after, size_t n)
{
    if (std::memcmp(&before[n], &after[n], (MAXN - n) * sizeof(fv_vec2)) != 0) {
        st.canary++;
        if (!st.shown_canary) { st.shown_canary = 1; printf("        [CANARY] %s: wrote o[n..] (n=%zu), out-of-bounds write\n", st.label, n); }
    }
}
static void canary_f(Stat& st, const std::vector<float>& before, const std::vector<float>& after, size_t n)
{
    if (std::memcmp(&before[n], &after[n], (MAXN - n) * sizeof(float)) != 0) {
        st.canary++;
        if (!st.shown_canary) { st.shown_canary = 1; printf("        [CANARY] %s: wrote o[n..] (n=%zu), out-of-bounds write\n", st.label, n); }
    }
}

/* ------------------------------------------------------------------ 数据语料 */

static const float kPool[] = {
    0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -0.5f, 2.0f, -2.0f,
    1e-30f, -1e-30f, 1e-45f, 3e-39f, /* 次正规 */ 1e-7f, -1e-7f,
    9.9999994e-7f, 1e-6f, 1.0000001e-6f, 2e-6f,   /* FV_EPS 邻域 */
    1e-3f, 1e3f, 1e6f, 1e30f, -1e30f, 3.4e38f, -3.4e38f, FLT_MIN,
    123.456f, -987.654f, 7.0f, 1.0000001f, 0.333333343f, 3.14159274f
};
static const size_t kPoolN = sizeof(kPool) / sizeof(kPool[0]);

/* NaN/Inf 语料：单独统计（edge），不计入主判定 */
static const float kEdgePool[] = {
    1e30f, -1e30f, INFINITY, -INFINITY, NAN, -NAN, 0.0f, -0.0f, 1.0f, -1.0f, 1e-6f, 3.4e38f
};
static const size_t kEdgePoolN = sizeof(kEdgePool) / sizeof(kEdgePool[0]);

/* 长度落在 FV_EPS 阈值两侧的定向语料（len_soa / normalize_soa 的判定边界） */
static const float kThresh[] = {
    0.0f, -0.0f, 1e-45f, 1e-30f, 5e-7f, 9.9999994e-7f, 1e-6f, 1.0000001e-6f,
    1.2e-6f, 1.5e-6f, 1e-5f, 1.0f, -1.0f, 1e30f, 3.4e38f, FLT_MIN
};
static const size_t kThreshN = sizeof(kThresh) / sizeof(kThresh[0]);

static const float kS[] = { 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, 2.0f, -3.5f, 1e-30f, 1e30f, 3.4e38f, 1e-7f, 123.456f };
static const size_t kSN = sizeof(kS) / sizeof(kS[0]);

static float pick(const float* pool, size_t pn, uint32_t salt, size_t i)
{
    return pool[hashu(salt ^ (uint32_t)(i * 2654435761u)) % pn];
}

static void fill_mix(std::vector<fv_vec2>& v, size_t n, uint32_t salt, bool edge)
{
    for (size_t i = 0; i < n; ++i)
        v[i] = fv_vec2_make(edge ? pick(kEdgePool, kEdgePoolN, salt, 2 * i)
                                 : pick(kPool, kPoolN, salt, 2 * i),
                            edge ? pick(kEdgePool, kEdgePoolN, salt, 2 * i + 1)
                                 : pick(kPool, kPoolN, salt, 2 * i + 1));
}

static void fill_thresh(std::vector<fv_vec2>& v, size_t n, uint32_t salt)
{
    for (size_t i = 0; i < n; ++i) {
        float s = kThresh[(i + salt) % kThreshN];
        switch ((i + salt) % 4) {
        case 0:  v[i] = fv_vec2_make(s, 0.0f); break;
        case 1:  v[i] = fv_vec2_make(s, s); break;
        case 2:  v[i] = fv_vec2_make(s, -s); break;
        default: v[i] = fv_vec2_make(0.0f, s); break;
        }
    }
}

static void to_soa(const std::vector<fv_vec2>& v, float* x, float* y, size_t n)
{
    for (size_t i = 0; i < n; ++i) { x[i] = v[i].x; y[i] = v[i].y; }
}
static void from_soa(const float* x, const float* y, float* oxy, size_t n)
{
    for (size_t i = 0; i < n; ++i) { oxy[2 * i] = x[i]; oxy[2 * i + 1] = y[i]; }
}

/* ------------------------------------------------------------------ 运行器 */

typedef void (*Batch2)(const fv_vec2*, const fv_vec2*, fv_vec2*, size_t);
typedef void (*BatchS)(const fv_vec2*, fv_real, fv_vec2*, size_t);
typedef void (*BatchL)(const fv_vec2*, const fv_vec2*, fv_real, fv_vec2*, size_t);

static const int DATASETS = 3;  /* 主语料数据集个数（每个 n） */

typedef fv_vec2 (*Scal2)(fv_vec2, fv_vec2);
typedef fv_vec2 (*ScalS)(fv_vec2, fv_real);
typedef fv_vec2 (*ScalL)(fv_vec2, fv_vec2, fv_real);

static void run2(Stat& st, Batch2 batch, Scal2 libref, Scal2 strictref)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN), o(MAXN), rl(MAXN), rs(MAXN), before(MAXN);
    for (int edge = 0; edge <= 1; ++edge) {
        int nds = edge ? 1 : DATASETS;
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < nds; ++d) {
                uint32_t salt = (uint32_t)(0x1000u + (uint32_t)d * 977u);
                fill_mix(a, n, salt, edge != 0);
                fill_mix(b, n, salt + 31u, edge != 0);
                fill_sentinel(o);
                before = o;
                batch(a.data(), b.data(), o.data(), n);
                canary(st, before, o, n);
                for (size_t i = 0; i < n; ++i) { rl[i] = libref(a[i], b[i]); rs[i] = strictref(a[i], b[i]); }
                cmp_floats(st, (const float*)o.data(), (const float*)rl.data(), 2 * n, 2, 0, edge != 0);
                cmp_floats(st, (const float*)o.data(), (const float*)rs.data(), 2 * n, 2, 1, edge != 0);
            }
        }
    }
}

static void runS(Stat& st, BatchS batch, ScalS libref, ScalS strictref)
{
    std::vector<fv_vec2> a(MAXN), o(MAXN), rl(MAXN), rs(MAXN), before(MAXN);
    for (int ds = 0; ds < (int)kSN; ++ds) {
        fv_real s = kS[ds];
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < DATASETS; ++d) {
                uint32_t salt = (uint32_t)(0x2000u + (uint32_t)d * 541u);
                fill_mix(a, n, salt, false);
                fill_sentinel(o);
                before = o;
                batch(a.data(), s, o.data(), n);
                canary(st, before, o, n);
                for (size_t i = 0; i < n; ++i) { rl[i] = libref(a[i], s); rs[i] = strictref(a[i], s); }
                cmp_floats(st, (const float*)o.data(), (const float*)rl.data(), 2 * n, 2, 0, false);
                cmp_floats(st, (const float*)o.data(), (const float*)rs.data(), 2 * n, 2, 1, false);
            }
        }
    }
}

static void runL(Stat& st, BatchL batch, ScalL libref, ScalL strictref)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN), o(MAXN), rl(MAXN), rs(MAXN), before(MAXN);
    for (int ds = 0; ds < (int)kSN; ++ds) {
        fv_real t = kS[ds];
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < DATASETS; ++d) {
                uint32_t salt = (uint32_t)(0x3000u + (uint32_t)d * 113u);
                fill_mix(a, n, salt, false);
                fill_mix(b, n, salt + 71u, false);
                fill_sentinel(o);
                before = o;
                batch(a.data(), b.data(), t, o.data(), n);
                canary(st, before, o, n);
                for (size_t i = 0; i < n; ++i) { rl[i] = libref(a[i], b[i], t); rs[i] = strictref(a[i], b[i], t); }
                cmp_floats(st, (const float*)o.data(), (const float*)rl.data(), 2 * n, 2, 0, false);
                cmp_floats(st, (const float*)o.data(), (const float*)rs.data(), 2 * n, 2, 1, false);
            }
        }
    }
}

/* ---- 掩码变体：o[i] = a[i]+b[i] if mask[i]!=0 else a[i]（逐位等于 a[i] 原值） ---- */
static void run_mask(Stat& st, const char* label)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN), o(MAXN), rl(MAXN), rs(MAXN), before(MAXN);
    std::vector<uint8_t> mask(MAXN);
    for (int edge = 0; edge <= 1; ++edge) {
        int nds = edge ? 1 : 4;
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < nds; ++d) {
                uint32_t salt = (uint32_t)(0x4000u + (uint32_t)d * 313u);
                fill_mix(a, n, salt, edge != 0);
                fill_mix(b, n, salt + 17u, edge != 0);
                for (size_t i = 0; i < n; ++i) {
                    switch (d) {
                    case 0:  mask[i] = 0xff; break;                       /* 全开 */
                    case 1:  mask[i] = 0x00; break;                       /* 全关 */
                    case 2:  mask[i] = (i & 1) ? 0x00 : 0x01; break;      /* 交替 */
                    default: mask[i] = (hashu(salt + (uint32_t)i) & 3u) ? (uint8_t)(1u << (hashu((uint32_t)i) & 7u)) : 0x00; break;
                    }
                }
                fill_sentinel(o);
                before = o;
                fv_vec2_add_x4_mask(a.data(), b.data(), o.data(), mask.data(), n);
                canary(st, before, o, n);
                for (size_t i = 0; i < n; ++i) {
                    /* 参考：掩码外必须逐位等于 a[i]，不是「不写」 */
                    rl[i] = mask[i] ? fv_vec2_add(a[i], b[i]) : a[i];
                    rs[i] = mask[i] ? fvs_add(a[i], b[i]) : a[i];
                }
                cmp_floats(st, (const float*)o.data(), (const float*)rl.data(), 2 * n, 2, 0, edge != 0);
                cmp_floats(st, (const float*)o.data(), (const float*)rs.data(), 2 * n, 2, 1, edge != 0);
            }
        }
    }
    (void)label;
}

/* ---- 索引变体：对 k in [0,m) 逐个应用 i=idx[k]; o[i]=a[i]+b[i]（含重复、乱序） ---- */
static void run_idx(Stat& st, const char* label, bool is_scale)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN), o(MAXN), rl(MAXN), rs(MAXN), before(MAXN);
    std::vector<uint32_t> idx(MAXN);
    const fv_real s = (fv_real)-2.5f;
    for (size_t n = 1; n <= MAX_TAIL; ++n) {
        for (int d = 0; d < 4; ++d) {
            uint32_t salt = (uint32_t)(0x5000u + (uint32_t)d * 197u);
            fill_mix(a, n, salt, false);
            fill_mix(b, n, salt + 23u, false);
            size_t m;
            switch (d) {
            case 0:  /* 顺序、无重复：m = n */
                m = n;
                for (size_t k = 0; k < m; ++k) idx[k] = (uint32_t)k;
                break;
            case 1:  /* 逆序 */
                m = n;
                for (size_t k = 0; k < m; ++k) idx[k] = (uint32_t)(n - 1 - k);
                break;
            case 2:  /* 大量重复 */
                m = n;
                for (size_t k = 0; k < m; ++k) idx[k] = (uint32_t)((k / 2) % n);
                break;
            case 3:  /* 随机乱序 + 重复 + m <= n */
                m = 1 + (hashu(salt) % n);
                for (size_t k = 0; k < m; ++k) idx[k] = hashu(salt + (uint32_t)k * 7919u) % (uint32_t)n;
                break;
            default: m = 0; break;
            }
            for (size_t k = 0; k < m; ++k) if (idx[k] >= n) { printf("        [BUG] generator index out of range\n"); std::exit(3); }

            fill_sentinel(o);
            before = o;
            if (is_scale) fv_vec2_scale_x4_idx(a.data(), s, o.data(), idx.data(), m);
            else          fv_vec2_add_x4_idx(a.data(), b.data(), o.data(), idx.data(), m);

            canary(st, before, o, n);
            /* 未被 idx 引用的 o[i]（i<n）不允许被改写 */
            std::vector<char> touched(n, 0);
            for (size_t k = 0; k < m; ++k) touched[idx[k]] = 1;
            for (size_t i = 0; i < n; ++i) {
                if (!touched[i] && (fbits(o[i].x) != fbits(before[i].x) || fbits(o[i].y) != fbits(before[i].y))) {
                    st.untouched++;
                    if (!st.shown_untouched) {
                        st.shown_untouched = 1;
                        printf("        [UNTOUCHED] %s: o[%zu] not referenced by idx was modified (m=%zu n=%zu) got=0x%08x/0x%08x keep=0x%08x/0x%08x\n",
                               st.label, i, m, n, fbits(o[i].x), fbits(o[i].y), fbits(before[i].x), fbits(before[i].y));
                    }
                }
            }
            /* 参考：按顺序逐个应用（重复索引下最后一次写入胜出，值相同） */
            for (size_t i = 0; i < n; ++i) { rl[i] = before[i]; rs[i] = before[i]; }
            for (size_t k = 0; k < m; ++k) {
                uint32_t i = idx[k];
                rl[i] = is_scale ? fv_vec2_scale(a[i], s) : fv_vec2_add(a[i], b[i]);
                rs[i] = is_scale ? fvs_scale(a[i], s) : fvs_add(a[i], b[i]);
            }
            /* 只比较被引用的位置（其余位置由 UNTOUCHED 检查单独负责） */
            for (size_t i = 0; i < n; ++i) {
                if (!touched[i]) continue;
                cmp_floats(st, (const float*)&o[i], (const float*)&rl[i], 2, 2, 0, false);
                cmp_floats(st, (const float*)&o[i], (const float*)&rs[i], 2, 2, 1, false);
            }
        }
    }
    (void)label;
}

/* ---- SoA 变体 ---- */
static void run_soa_add(Stat& st)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN), o(MAXN), rl(MAXN), rs(MAXN);
    std::vector<float> ax(MAXN), ay(MAXN), bx(MAXN), by(MAXN), ox(MAXN), oy(MAXN);
    std::vector<float> oxb(MAXN), oyb(MAXN);
    for (int edge = 0; edge <= 1; ++edge) {
        int nds = edge ? 1 : DATASETS;
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < nds; ++d) {
                uint32_t salt = (uint32_t)(0x6000u + (uint32_t)d * 89u);
                fill_mix(a, n, salt, edge != 0);
                fill_mix(b, n, salt + 5u, edge != 0);
                to_soa(a, ax.data(), ay.data(), n);
                to_soa(b, bx.data(), by.data(), n);
                fill_sentinel_f(ox); fill_sentinel_f(oy);
                oxb = ox; oyb = oy;
                fv_vec2_add_soa(ax.data(), ay.data(), bx.data(), by.data(), ox.data(), oy.data(), n);
                canary_f(st, oxb, ox, n);
                canary_f(st, oyb, oy, n);
                from_soa(ox.data(), oy.data(), (float*)o.data(), n);
                for (size_t i = 0; i < n; ++i) { rl[i] = fv_vec2_add(a[i], b[i]); rs[i] = fvs_add(a[i], b[i]); }
                cmp_floats(st, (const float*)o.data(), (const float*)rl.data(), 2 * n, 2, 0, edge != 0);
                cmp_floats(st, (const float*)o.data(), (const float*)rs.data(), 2 * n, 2, 1, edge != 0);
            }
        }
    }
}

static void run_soa_scale(Stat& st)
{
    std::vector<fv_vec2> a(MAXN), rl(MAXN), rs(MAXN);
    std::vector<float> ax(MAXN), ay(MAXN), ox(MAXN), oy(MAXN), oxb(MAXN), oyb(MAXN);
    std::vector<float> flat(2 * MAXN);
    for (int ds = 0; ds < (int)kSN; ++ds) {
        fv_real s = kS[ds];
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < DATASETS; ++d) {
                uint32_t salt = (uint32_t)(0x7000u + (uint32_t)d * 61u);
                fill_mix(a, n, salt, false);
                to_soa(a, ax.data(), ay.data(), n);
                fill_sentinel_f(ox); fill_sentinel_f(oy);
                oxb = ox; oyb = oy;
                fv_vec2_scale_soa(ax.data(), ay.data(), s, ox.data(), oy.data(), n);
                canary_f(st, oxb, ox, n);
                canary_f(st, oyb, oy, n);
                from_soa(ox.data(), oy.data(), flat.data(), n);
                for (size_t i = 0; i < n; ++i) { rl[i] = fv_vec2_scale(a[i], s); rs[i] = fvs_scale(a[i], s); }
                cmp_floats(st, flat.data(), (const float*)rl.data(), 2 * n, 2, 0, false);
                cmp_floats(st, flat.data(), (const float*)rs.data(), 2 * n, 2, 1, false);
            }
        }
    }
}

static void run_soa_len(Stat& st)
{
    std::vector<fv_vec2> a(MAXN);
    std::vector<float> ax(MAXN), ay(MAXN), o(MAXN), ob(MAXN), rl(MAXN), rs(MAXN);
    for (int edge = 0; edge <= 1; ++edge) {
        int nds = edge ? 1 : DATASETS;
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < nds; ++d) {
                uint32_t salt = (uint32_t)(0x8000u + (uint32_t)d * 43u);
                /* 一半数据集用阈值定向语料 */
                if (!edge && (d & 1)) fill_thresh(a, n, salt); else fill_mix(a, n, salt, edge != 0);
                to_soa(a, ax.data(), ay.data(), n);
                fill_sentinel_f(o);
                ob = o;
                fv_vec2_len_soa(ax.data(), ay.data(), o.data(), n);
                canary_f(st, ob, o, n);
                for (size_t i = 0; i < n; ++i) { rl[i] = fv_vec2_len(a[i]); rs[i] = fvs_len(a[i]); }
                cmp_floats(st, o.data(), rl.data(), n, 1, 0, edge != 0);
                cmp_floats(st, o.data(), rs.data(), n, 1, 1, edge != 0);
            }
        }
    }
}

static void run_soa_normalize(Stat& st)
{
    std::vector<fv_vec2> a(MAXN), rl(MAXN), rs(MAXN);
    std::vector<float> ax(MAXN), ay(MAXN), ox(MAXN), oy(MAXN), oxb(MAXN), oyb(MAXN);
    std::vector<float> flat(2 * MAXN);
    for (int edge = 0; edge <= 1; ++edge) {
        int nds = edge ? 1 : DATASETS;
        for (size_t n = 0; n <= MAX_TAIL; ++n) {
            for (int d = 0; d < nds; ++d) {
                uint32_t salt = (uint32_t)(0x9000u + (uint32_t)d * 37u);
                if (!edge && (d & 1)) fill_thresh(a, n, salt); else fill_mix(a, n, salt, edge != 0);
                to_soa(a, ax.data(), ay.data(), n);
                fill_sentinel_f(ox); fill_sentinel_f(oy);
                oxb = ox; oyb = oy;
                fv_vec2_normalize_soa(ax.data(), ay.data(), ox.data(), oy.data(), n);
                canary_f(st, oxb, ox, n);
                canary_f(st, oyb, oy, n);
                from_soa(ox.data(), oy.data(), flat.data(), n);
                for (size_t i = 0; i < n; ++i) { rl[i] = fv_vec2_normalize(a[i]); rs[i] = fvs_normalize(a[i]); }
                cmp_floats(st, flat.data(), (const float*)rl.data(), 2 * n, 2, 0, edge != 0);
                cmp_floats(st, flat.data(), (const float*)rs.data(), 2 * n, 2, 1, edge != 0);
            }
        }
    }
}

/* ---- dot_soa：规格允许改变归约顺序 → 只做误差统计，不做位级判定 ---- */
struct DotStat {
    long n_cases = 0; double max_abs = 0.0; double max_rel = 0.0; long exact = 0; long total = 0;
    /* 头文件声称：n < 8 时 dot_soa 退化为纯顺序累加，与逐元素路径逐位一致 → 逐 n 核对该断言 */
    long exact_by_n[18]; long total_by_n[18]; int shown_small = 0;
};
static DotStat g_dot = {};

static void run_soa_dot(void)
{
    std::vector<fv_vec2> a(MAXN), b(MAXN);
    std::vector<float> ax(MAXN), ay(MAXN), bx(MAXN), by(MAXN);
    for (size_t n = 0; n <= MAX_TAIL; ++n) {
        for (int d = 0; d < 8; ++d) {
            uint32_t salt = (uint32_t)(0xa000u + (uint32_t)d * 127u);
            fill_mix(a, n, salt, false);
            fill_mix(b, n, salt + 11u, false);
            to_soa(a, ax.data(), ay.data(), n);
            to_soa(b, bx.data(), by.data(), n);
            float got = fv_vec2_dot_soa(ax.data(), ay.data(), bx.data(), by.data(), n);
            double want = 0.0;
            for (size_t i = 0; i < n; ++i) want += (double)fv_vec2_dot(a[i], b[i]);   /* float 逐个累加，顺序与标量一致 */
            float wantf = 0.0f;
            for (size_t i = 0; i < n; ++i) wantf = (float)(wantf + fv_vec2_dot(a[i], b[i]));
            double abserr = std::fabs((double)got - (double)wantf);
            double rel = (std::fabs((double)wantf) > 0.0) ? abserr / std::fabs((double)wantf) : 0.0;
            g_dot.n_cases++;
            g_dot.total++;
            g_dot.total_by_n[n]++;
            if (fbits(got) == fbits(wantf)) { g_dot.exact++; g_dot.exact_by_n[n]++; }
            else if (n < 8 && !g_dot.shown_small) {
                g_dot.shown_small = 1;
                printf("        fv_vec2_dot_soa n<8 assertion failed: n=%zu got=0x%08x (%.9g) scalar sequential sum=0x%08x (%.9g)\n",
                       n, fbits(got), (double)got, fbits(wantf), (double)wantf);
            }
            if (abserr > g_dot.max_abs) g_dot.max_abs = abserr;
            if (rel > g_dot.max_rel) g_dot.max_rel = rel;
        }
    }
}

/* ------------------------------------------------------------------ 矩阵点变换（P4）
 * fv_mat4_transform_points / _dirs 没有对应的标量 API（它本身就是「一次调用处理 N 个点」的
 * 最小粒度），所以这里只与**独立参考**比较；两条入口的差别只在平移列参与不参与。
 * 参考实现与头文件里标量尾巴逐字同序：((m0*x + m1*y) + m2*z) + m3，dirs 没有第四项。
 * 本文件由 verify/Makefile 以 -ffp-contract=off 编译，这两个 mul+add 不会被收缩成 FMA。 */

static void mat4_ref_points(const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        const float x = in[i].x, y = in[i].y, z = in[i].z;
        out[i].x = ((m->m[0] * x + m->m[1] * y) + m->m[2] * z) + m->m[3];
        out[i].y = ((m->m[4] * x + m->m[5] * y) + m->m[6] * z) + m->m[7];
        out[i].z = ((m->m[8] * x + m->m[9] * y) + m->m[10] * z) + m->m[11];
        out[i].w = 1.0f;
    }
}

static void mat4_ref_dirs(const fv_mat4* m, const fv_vec4* in, fv_vec4* out, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        const float x = in[i].x, y = in[i].y, z = in[i].z;
        out[i].x = (m->m[0] * x + m->m[1] * y) + m->m[2] * z;
        out[i].y = (m->m[4] * x + m->m[5] * y) + m->m[6] * z;
        out[i].z = (m->m[8] * x + m->m[9] * y) + m->m[10] * z;
        out[i].w = 0.0f;
    }
}

/* 三组语料，按「一次加法里最多出现一种 NaN 位模式」设计（理由见 include/flashvec.h 的 AoS 一节）：
 *   kPool     —— 有限值，含 -0.0 / 次正规 / FV_EPS 邻域，没有 inf / NaN；
 *   kMatEdge  —— ±inf / FLT_MAX / 次正规 / -0.0，**没有 NaN**：NaN 只由 inf*0 / inf-inf 生成，
 *                硬件给出的 indefinite 恒为同一位模式，所以严格位比较仍然成立；
 *   kMatNan   —— 只有一种 NaN 位模式（+NAN）且不含 inf：NaN 只会作为操作数传播，同样唯一。
 * 混用两种 NaN 位模式会产生「结果取哪一个 NaN 由指令操作数顺序决定」的格子，那属于契约不承诺的
 * 部分，不该混进这条统计（要验它请用 tests/test_mat4.cpp 的第 6 节，那里单独计数）。 */
static const float kMatEdge[] = {
    0.0f, -0.0f, 1.0f, -1.0f, 0.5f, -3.25f, FLT_MAX, -FLT_MAX, FLT_MIN,
    1e-45f, -1e-45f, INFINITY, -INFINITY
};
static const size_t kMatEdgeN = sizeof(kMatEdge) / sizeof(kMatEdge[0]);

static const float kMatNan[] = { NAN, 0.0f, -0.0f, 1.0f, -1.0f, 0.5f, 1e30f, -1e30f };
static const size_t kMatNanN = sizeof(kMatNan) / sizeof(kMatNan[0]);

#define MAT_DS 5 /* 0..2 = kPool，3 = kMatEdge，4 = kMatNan */

static const float* mat_pool(int ds)
{
    return (ds == 3) ? kMatEdge : (ds == 4) ? kMatNan : kPool;
}
static size_t mat_pool_n(int ds)
{
    return (ds == 3) ? kMatEdgeN : (ds == 4) ? kMatNanN : kPoolN;
}

static void fill_v4(std::vector<fv_vec4>& v, size_t n, uint32_t salt, int ds)
{
    const float* pool = mat_pool(ds);
    const size_t pn = mat_pool_n(ds);
    for (size_t i = 0; i < n; ++i) {
        v[i].x = pick(pool, pn, salt, 4 * i + 0);
        v[i].y = pick(pool, pn, salt, 4 * i + 1);
        v[i].z = pick(pool, pn, salt, 4 * i + 2);
        v[i].w = pick(pool, pn, salt, 4 * i + 3); /* 被两条入口忽略，故意给垃圾值 */
    }
}

static fv_mat4 make_mat(uint32_t salt, int ds)
{
    fv_mat4 m;
    /* ds == 4 的矩阵必须是有限值：语料里只有一种 NaN 位模式，不能再引入 inf 生成的 indefinite */
    const float* pool = kPool;
    const size_t pn = kPoolN;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c) {
            m.m[r * 4 + c] = pick(pool, pn, salt, (size_t)(r * 4 + c));
        }
    }
    if (ds == 3) {
        m.m[1] = INFINITY;  /* 无穷系数：与 ±inf 输入相乘仍是 inf，与 0 相乘得到 indefinite */
        m.m[6] = -0.0f;
        m.m[10] = 1e-45f;
    }
    m.m[12] = 0.0f;
    m.m[13] = 0.0f;
    m.m[14] = 0.0f;
    m.m[15] = 1.0f;
    return m;
}

static void canary_v4(Stat& st, const std::vector<fv_vec4>& before, const std::vector<fv_vec4>& after, size_t n)
{
    if (std::memcmp(&before[n], &after[n], (MAXN - n) * sizeof(fv_vec4)) != 0) {
        st.canary++;
        if (!st.shown_canary) {
            st.shown_canary = 1;
            printf("        [CANARY] %s: wrote o[n..] (n=%zu), out-of-bounds write\n", st.label, n);
        }
    }
}

typedef void (*BatchTF)(const fv_mat4*, const fv_vec4*, fv_vec4*, size_t);

static void run_mat4(Stat& st, BatchTF batch, bool dirs)
{
    /* 0..17 覆盖 _x2 内核的两种尾长（偶数走主循环、奇数留 1 个尾巴），后面补奇数与整批 */
    static const size_t kLens[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                    31, 33, 63, 64 };
    const size_t nlens = sizeof(kLens) / sizeof(kLens[0]);
    std::vector<fv_vec4> a(MAXN), o(MAXN), r(MAXN), rt(MAXN), before(MAXN);

    for (int ds = 0; ds < MAT_DS; ++ds) {
        const bool edge = (ds >= 3);
        for (size_t li = 0; li < nlens; ++li) {
            const size_t n = kLens[li];
            const uint32_t salt = 0x4000u + (uint32_t)ds * 733u;
            fv_mat4 m = make_mat(salt + 17u, ds);
            fill_v4(a, n, salt, ds);
            fill_v4(o, MAXN, salt + 3u, ds); /* 哨兵值 */
            before = o;
            batch(&m, a.data(), o.data(), n);
            canary_v4(st, before, o, n);

            /* which = 0：库自己的标量路径。这个入口没有「标量 API」，最贴近的就是 n = 1 的调用
             * —— bulk = 0，整个调用走标量表达式。它同时验证「同一元素的结果不随 n 改变」，
             * 并且让这一行的 samples 落到与其它入口同一张表里可比的量级（否则会是 0，
             * 看起来像没验证过，而这个项目最不能容忍的就是「看起来验证过了」）。 */
            for (size_t i = 0; i < n; ++i) {
                batch(&m, &a[i], &rt[i], 1);
            }
            cmp_floats(st, (const float*)o.data(), (const float*)rt.data(), 4 * n, 4, 0, edge);

            /* which = 1：独立标量参考（本 TU 的 mat4_ref_*） */
            if (dirs) {
                mat4_ref_dirs(&m, a.data(), r.data(), n);
            } else {
                mat4_ref_points(&m, a.data(), r.data(), n);
            }
            cmp_floats(st, (const float*)o.data(), (const float*)r.data(), 4 * n, 4, 1, edge);
        }
    }
}

/* ------------------------------------------------------------------ P0 检查 */

static int g_fail_p0 = 0;
static void p0_line(const char* what, bool ok, const char* detail)
{
    printf("  %-42s %s  %s\n", what, ok ? "PASS" : "FAIL", detail ? detail : "");
    if (!ok) g_fail_p0++;
}

static void run_p0(void)
{
    printf("\n[P0] constants and small helpers\n");
    bool eps_ok = (FV_EPS == (fv_real)1e-6);
    bool sqr_ok = (FV_EPS_SQR == (fv_real)1e-12) && (FV_EPS_SQR == FV_EPS * FV_EPS);
    char buf[160];
    std::snprintf(buf, sizeof buf, "FV_EPS=%.9g FV_EPS_SQR=%.9g (FV_EPS*FV_EPS=%.9g) bits=FV_EPS_SQR:0x%08x",
                  (double)FV_EPS, (double)FV_EPS_SQR, (double)(FV_EPS * FV_EPS), fbits(FV_EPS_SQR));
    p0_line("FV_EPS / FV_EPS_SQR values", eps_ok && sqr_ok, buf);

    /* min/max/madd：纯值语义，避开 ±0 与 NaN 的未定义边界 */
    int bad = 0, cases = 0;
    const float vals[] = { -1e30f, -2.0f, -1.0f, -0.5f, 0.5f, 1.0f, 2.0f, 1e30f, 1e-30f, -1e-30f };
    for (size_t i = 0; i < sizeof(vals) / sizeof(vals[0]); ++i)
        for (size_t j = 0; j < sizeof(vals) / sizeof(vals[0]); ++j) {
            fv_vec2 a = fv_vec2_make(vals[i], vals[j]);
            fv_vec2 b = fv_vec2_make(vals[j], vals[i]);
            fv_vec2 mn = fv_vec2_min(a, b), mx = fv_vec2_max(a, b);
            fv_vec2 md = fv_vec2_madd(a, b, (fv_real)-1.75f);
            /* 期望 = 逐分量三元比较（规格未定义 ±0/NaN 语义，这里只用非零同号数据） */
            fv_vec2 emn = fv_vec2_make(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y);
            fv_vec2 emx = fv_vec2_make(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y);
            fv_vec2 emd = fv_vec2_make(a.x + b.x * (fv_real)-1.75f, a.y + b.y * (fv_real)-1.75f);
            cases++;
            if (!(mn.x == emn.x && mn.y == emn.y)) bad++;
            if (!(mx.x == emx.x && mx.y == emx.y)) bad++;
            if (!(md.x == emd.x && md.y == emd.y)) bad++;
        }
    std::snprintf(buf, sizeof buf, "cases=%d mismatches=%d", cases, bad);
    p0_line("fv_vec2 min/max/madd", bad == 0, buf);

    /* fv_vec2_clamp：单独一次，含越界与界内两种取值 */
    {
        int cbad = 0, ccases = 0;
        for (size_t i = 0; i < sizeof(vals) / sizeof(vals[0]); ++i)
            for (size_t j = 0; j < sizeof(vals) / sizeof(vals[0]); ++j) {
                fv_vec2 lo = fv_vec2_make(vals[i], vals[j]);
                fv_vec2 hi = fv_vec2_make(-vals[i], -vals[j]);
                if (lo.x > hi.x || lo.y > hi.y) continue;   /* 只测 lo <= hi 的合法区间 */
                fv_vec2 vs[3] = { fv_vec2_make(lo.x, lo.y), fv_vec2_make(hi.x, hi.y), fv_vec2_make(lo.x - 1.0f, hi.y + 1.0f) };
                for (int t = 0; t < 3; ++t) {
                    fv_vec2 v = vs[t];
                    fv_vec2 cl = fv_vec2_clamp(v, lo, hi);
                    float ex = v.x < lo.x ? lo.x : (v.x > hi.x ? hi.x : v.x);
                    float ey = v.y < lo.y ? lo.y : (v.y > hi.y ? hi.y : v.y);
                    ccases++;
                    if (cl.x != ex || cl.y != ey) cbad++;
                }
            }
        std::snprintf(buf, sizeof buf, "cases=%d mismatches=%d", ccases, cbad);
        p0_line("fv_vec2_clamp clamping", cbad == 0, buf);
    }

    bad = 0; cases = 0;
    for (size_t i = 0; i < sizeof(vals) / sizeof(vals[0]); ++i)
        for (size_t j = 0; j < sizeof(vals) / sizeof(vals[0]); ++j)
            for (size_t k = 0; k < sizeof(vals) / sizeof(vals[0]); ++k) {
                fv_vec3 a = fv_vec3_make(vals[i], vals[j], vals[k]);
                fv_vec3 b = fv_vec3_make(vals[k], vals[i], vals[j]);
                fv_vec3 mn = fv_vec3_min(a, b), mx = fv_vec3_max(a, b);
                fv_vec3 md = fv_vec3_madd(a, b, (fv_real)0.5f);
                /* clamp 用 lo=mn, hi=mx（保证 lo<=hi），分别夹 a、b 与一个区间外的向量 */
                fv_vec3 out = fv_vec3_make(mn.x - 1.0f, mx.y + 1.0f, mn.z - 1.0f);
                fv_vec3 cls[3] = { fv_vec3_clamp(a, mn, mx), fv_vec3_clamp(b, mn, mx), fv_vec3_clamp(out, mn, mx) };
                fv_vec3 emn = fv_vec3_make(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y, a.z < b.z ? a.z : b.z);
                fv_vec3 emx = fv_vec3_make(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y, a.z > b.z ? a.z : b.z);
                fv_vec3 emd = fv_vec3_make(a.x + b.x * (fv_real)0.5f, a.y + b.y * (fv_real)0.5f, a.z + b.z * (fv_real)0.5f);
                cases++;
                if (!(mn.x == emn.x && mn.y == emn.y && mn.z == emn.z)) bad++;
                if (!(mx.x == emx.x && mx.y == emx.y && mx.z == emx.z)) bad++;
                if (!(md.x == emd.x && md.y == emd.y && md.z == emd.z)) bad++;
                for (int t = 0; t < 3; ++t) {
                    const fv_vec3* v = (t == 0) ? &a : (t == 1 ? &b : &out);
                    float ex = v->x < mn.x ? mn.x : (v->x > mx.x ? mx.x : v->x);
                    float ey = v->y < mn.y ? mn.y : (v->y > mx.y ? mx.y : v->y);
                    float ez = v->z < mn.z ? mn.z : (v->z > mx.z ? mx.z : v->z);
                    if (cls[t].x != ex || cls[t].y != ey || cls[t].z != ez) bad++;
                }
            }
    std::snprintf(buf, sizeof buf, "cases=%d mismatches=%d", cases, bad);
    p0_line("fv_vec3 min/max/clamp/madd", bad == 0, buf);

    /* aabb2：merge / overlaps */
    {
        fv_aabb2 box = { fv_vec2_make(0, 0), fv_vec2_make(1, 1) };
        struct { float x, y; } pts[] = { { -1, 0.5f }, { 0.5f, -1 }, { 2, 0.5f }, { 0.5f, 2 }, { 0.5f, 0.5f }, { -0.0f, 1.0f } };
        fv_aabb2 boxes[6];
        for (int i = 0; i < 6; ++i) boxes[i] = fv_aabb2_merge(box, fv_vec2_make(pts[i].x, pts[i].y));
        int mergebad = 0;
        for (int i = 0; i < 6; ++i) {
            float elx = pts[i].x < 0.0f ? pts[i].x : 0.0f, ely = pts[i].y < 0.0f ? pts[i].y : 0.0f;
            float ehx = pts[i].x > 1.0f ? pts[i].x : 1.0f, ehy = pts[i].y > 1.0f ? pts[i].y : 1.0f;
            if (boxes[i].lo.x != elx || boxes[i].lo.y != ely || boxes[i].hi.x != ehx || boxes[i].hi.y != ehy) mergebad++;
        }
        std::snprintf(buf, sizeof buf, "6 points, mismatches=%d", mergebad);
        p0_line("fv_aabb2_merge(box,p) = [min(lo,p),max(hi,p)]", mergebad == 0, buf);

        fv_aabb2 A = { fv_vec2_make(0, 0), fv_vec2_make(1, 1) };
        fv_aabb2 B = { fv_vec2_make(0.5f, 0.5f), fv_vec2_make(2, 2) };
        fv_aabb2 C = { fv_vec2_make(2, 2), fv_vec2_make(3, 3) };
        fv_aabb2 D = { fv_vec2_make(1.5f, -1), fv_vec2_make(2.5f, 0.5f) };
        int ov = fv_aabb2_overlaps(A, B) + fv_aabb2_overlaps(A, A) + fv_aabb2_overlaps(A, C) + fv_aabb2_overlaps(A, D);
        std::snprintf(buf, sizeof buf, "A∩B=%d A∩A=%d(far)=%d A∩D=%d  (expected 1,1,0,0)", fv_aabb2_overlaps(A, B), fv_aabb2_overlaps(A, A), fv_aabb2_overlaps(A, C), fv_aabb2_overlaps(A, D));
        p0_line("fv_aabb2_overlaps clear overlap/separate", ov == 2, buf);
        /* 只接触边界的语义规格未定义 → 只报告，不判定 */
        printf("  %-42s INFO  A(lo=1,1) vs C(lo=2,2) touch only at a corner: overlaps=%d (spec undefined, recorded only)\n", "fv_aabb2_overlaps edge touch", fv_aabb2_overlaps(A, C));
    }

    /* ±0 / NaN 边界行为：规格未定义，仅记录 */
    {
        fv_vec2 p = fv_vec2_make(-0.0f, 0.0f), q = fv_vec2_make(0.0f, -0.0f);
        fv_vec2 mn = fv_vec2_min(p, q), mx = fv_vec2_max(p, q);
        printf("  %-42s INFO  min((-0,0),(0,-0))=(0x%08x,0x%08x) max=(0x%08x,0x%08x)\n",
               "min/max signed-zero behavior", fbits(mn.x), fbits(mn.y), fbits(mx.x), fbits(mx.y));
    }

    /* 比较器自检：证明 NaN 策略与 ±0 分类按文档工作 */
    {
        Same a = same_bits(ffrom(0x7fc00000u), ffrom(0x7fc0deadu));
        Same b = same_bits(0.0f, -0.0f);
        Same c = same_bits(1.0f, 1.00000012f);
        Same d = same_bits(ffrom(0x7f800000u), ffrom(0x7f800000u));
        bool ok = a.equal && a.nan_pair && !b.equal && b.szero && !c.equal && !c.szero && d.equal;
        std::snprintf(buf, sizeof buf, "NaN(different payload)→equal=%d nan_pair=%d; (+0,-0)→equal=%d szero=%d; (1,1+1ulp)→equal=%d",
                      (int)a.equal, (int)a.nan_pair, (int)b.equal, (int)b.szero, (int)c.equal);
        p0_line("comparator self-check (NaN / +0-0 / 1ulp)", ok, buf);
    }
}

/* ------------------------------------------------------------------ 汇总 */

static std::vector<Stat*> g_stats;
static Stat g_pool[64];      /* 静态池：避免在 ASan/LSan 下产生我自己造成的泄漏噪音 */
static size_t g_pool_n = 0;

static Stat& reg(const char* label)
{
    if (g_pool_n >= 64) { std::fprintf(stderr, "statistics pool overflow\n"); std::exit(3); }
    Stat& s = g_pool[g_pool_n++];
    s.label = label;
    g_stats.push_back(&s);
    return s;
}

static int summarize(void)
{
    int fails = 0;
    printf("\n[element-wise batch path vs scalar path: bit-identical]\n");
    printf("  %-26s %-8s %-8s %-16s %-18s %s\n", "Function", "Verdict", "samples", "vs lib scalar", "vs independent ref", "note");
    for (size_t k = 0; k < g_stats.size(); ++k) {
        Stat& s = *g_stats[k];
        long hard = s.bad_lib + s.bad_strict + s.canary + s.untouched;
        const char* verdict = hard ? "FAIL" : (s.edge ? "PASS*" : "PASS");
        if (hard) fails++;
        char note[200];
        note[0] = 0;
        if (s.canary) std::snprintf(note + std::strlen(note), sizeof note - std::strlen(note), "out-of-bounds write=%ld ", s.canary);
        if (s.untouched) std::snprintf(note + std::strlen(note), sizeof note - std::strlen(note), "unreferenced slot modified=%ld ", s.untouched);
        if (s.edge) std::snprintf(note + std::strlen(note), sizeof note - std::strlen(note), "NaN/Inf input divergence=%ld ", s.edge);
        if (s.szero) std::snprintf(note + std::strlen(note), sizeof note - std::strlen(note), "±0 difference=%ld ", s.szero);
        if (s.nan_pairs) std::snprintf(note + std::strlen(note), sizeof note - std::strlen(note), "NaN pairs=%ld ", s.nan_pairs);
        printf("  %-26s %-8s %-8ld %-16s %-18s %s\n", s.label, verdict, s.samples,
               s.bad_lib ? "FAIL" : "PASS", s.bad_strict ? "FAIL" : "PASS", note);
    }
    if (g_dot.total) {
        printf("\n[dot_soa: spec allows reordering the reduction → no bit-level verdict]\n");
        printf("  %-26s %-8s %-8s %s\n", "Function", "Verdict", "cases", "note");
        printf("  %-26s %-8s %-8ld bit-identical %ld/%ld, max|abs err|=%.3g, max rel err=%.3g (float input, |value| up to 1e30)\n",
               "fv_vec2_dot_soa", "PASS*", g_dot.n_cases, g_dot.exact, g_dot.total, g_dot.max_abs, g_dot.max_rel);
        printf("  %-26s %s\n", "", "* dot_soa may reorder the reduction, so error statistics replace the bit comparison; every other path gets a bit-level verdict.");
        long small_e = 0, small_t = 0;
        for (size_t n = 0; n <= MAX_TAIL; ++n) { if (n < 8) { small_e += g_dot.exact_by_n[n]; small_t += g_dot.total_by_n[n]; } }
        printf("  %-26s %-8s n<8 bit-identical %ld/%ld  %s\n", "dot_soa n<8 header assert", small_e == small_t ? "PASS" : "FAIL", small_e, small_t,
               small_e == small_t ? "(header line 149 assertion holds)" : "(header line 149 assertion violated, see bit patterns above)");
        printf("  %-26s", "dot_soa exact bits per n");
        for (size_t n = 0; n <= MAX_TAIL; ++n) printf(" %zu:%ld/%ld", n, g_dot.exact_by_n[n], g_dot.total_by_n[n]);
        printf("\n");
    }
    return fails;
}

int main(int argc, char** argv)
{
    bool quick = (argc > 1 && std::strcmp(argv[1], "--p0") == 0);
    printf("flashvec parity — bit-identical verification (independent third party)\n");
    printf("build: %s %s\n", __DATE__, __TIME__);

    if (!quick) {
        Stat& s_add_x4   = reg("fv_vec2_add_x4");
        run2(s_add_x4, fv_vec2_add_x4, fv_vec2_add, fvs_add);
        Stat& s_sub_x4   = reg("fv_vec2_sub_x4");
        run2(s_sub_x4, fv_vec2_sub_x4, fv_vec2_sub, fvs_sub);
        Stat& s_scale_x4 = reg("fv_vec2_scale_x4");
        runS(s_scale_x4, fv_vec2_scale_x4, fv_vec2_scale, fvs_scale);
        Stat& s_lerp_x4  = reg("fv_vec2_lerp_x4");
        runL(s_lerp_x4, fv_vec2_lerp_x4, fv_vec2_lerp, fvs_lerp);
        Stat& s_madd_x4  = reg("fv_vec2_madd_x4");
        runL(s_madd_x4, fv_vec2_madd_x4, fv_vec2_madd, fvs_madd);
        Stat& s_mask     = reg("fv_vec2_add_x4_mask");
        run_mask(s_mask, "mask");
        Stat& s_idx_a    = reg("fv_vec2_add_x4_idx");
        run_idx(s_idx_a, "idx-add", false);
        Stat& s_idx_s    = reg("fv_vec2_scale_x4_idx");
        run_idx(s_idx_s, "idx-scale", true);
        Stat& s_soa_add  = reg("fv_vec2_add_soa");
        run_soa_add(s_soa_add);
        Stat& s_soa_scl  = reg("fv_vec2_scale_soa");
        run_soa_scale(s_soa_scl);
        Stat& s_soa_len  = reg("fv_vec2_len_soa");
        run_soa_len(s_soa_len);
        Stat& s_soa_nrm  = reg("fv_vec2_normalize_soa");
        run_soa_normalize(s_soa_nrm);
        run_soa_dot();
        /* 矩阵点变换：没有对应的标量 API，所以汇总表里「vs lib scalar」一列对它没有意义
         * （0 样本、恒 PASS），判定只看「vs independent ref」一列与 canary。 */
        Stat& s_m4p = reg("fv_mat4_transform_points");
        run_mat4(s_m4p, fv_mat4_transform_points, false);
        Stat& s_m4d = reg("fv_mat4_transform_dirs");
        run_mat4(s_m4d, fv_mat4_transform_dirs, true);
    }

    run_p0();

    int fails = summarize();
    printf("\nverdict: %s (FAIL items=%d)\n", fails ? "bit mismatches / memory-semantics problems" : "all element-wise paths bit-identical", fails);
    printf("note: comparison uses the IEEE bit pattern; the only exception is both sides NaN (payload may differ), which counts as equal and is tracked separately.\n");
    return fails ? 1 : 0;
}
