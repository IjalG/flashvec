/* verify/perf.cpp —— flashvec 批量 SIMD 路径基准（验收者 F 独立实现）
 *
 * 计时：clock_gettime(CLOCK_MONOTONIC)；每项 5 次取样取中位数（排序后第 3 个）。
 * 每次取样内层重复 iters 次（校准到 ~6 ms/取样），事先 warm-up，避免把首次触页/冷 cache 计进去。
 *
 * 数据完整性（这条决定了数字是否可比）：
 *   - 输入 a / b 是不可变输入，任何被计时实现都不写它们（全部是 out-of-place 的 o = f(a,b)）；
 *   - 输出缓冲 o 在使用前清零，从不当输入；
 *   - 校验和 = **单次调用**的输出累加（与 iters 无关）；
 *   - 校验和判定按 **op 形态** 分两类（不看实现名字，见 is_reduction_op()）：
 *       · 逐元素 op：除参考档 scalar_api 外，所有实现的校验和必须**位级相等**；
 *       · 归约 op（dot：结果是一个累加标量，规格与编译器都允许重结合）：用随项数增长的
 *         相对容差 8·eps·√n（见 reduction_tol()）—— 固定 1e-6 在 n 大时数学上不可能成立，
 *         n=1048576 实测 soa 相对偏差 2.35e-05。
 *     越界者该 op 的数据判为无效并计入 g_cs_fail（退出码 1）。
 *
 * 四档实现：
 *   scalar_api  —— 逐元素调用库标量 API（跨 TU 调用，不可内联）
 *   scalar_flat —— 同样的数学在基准 TU 内联展开（编译器可自行向量化），代表「免费午餐」上限
 *   x4          —— 库的 AoS 4 路交错批量 API
 *   soa         —— 库的 SoA 批量 API
 */
#include "flashvec.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/* 当前 CPU 频率（MHz）：用于把 ns/elem 换算成 cycle/elem，避免被频率调节误导 */
static double cpu_mhz(void)
{
    FILE* f = std::fopen("/proc/cpuinfo", "r");
    if (!f) return 0.0;
    char line[256];
    double sum = 0.0; int cnt = 0;
    while (std::fgets(line, sizeof line, f)) {
        double v;
        if (std::sscanf(line, "cpu MHz : %lf", &v) == 1) { sum += v; cnt++; }
    }
    std::fclose(f);
    return cnt ? sum / cnt : 0.0;
}

static void* xalloc(size_t bytes)
{
    void* p = std::aligned_alloc(64, (bytes + 63u) & ~(size_t)63u);
    if (!p) { std::fprintf(stderr, "alloc failed (%zu bytes)\n", bytes); std::exit(2); }
    return p;
}

struct Buf {
    size_t n;
    fv_vec2 *a, *b, *o;                      /* a/b 不可变输入，o 是输出（从不当输入） */
    fv_vec2 *nrm, *vel, *vel0;               /* collision 微内核的数据 */
    float *ax, *ay, *bx, *by, *ox, *oy;      /* a/b 的 SoA 视图 */
    float *nx, *ny, *nox, *noy;
    float *vx, *vy, *vx0, *vy0;
    float *lenout;
};

static uint64_t g_rng = 0x9e3779b97f4a7c15ull;
static double rnd01(void)
{
    g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17;
    return (double)(g_rng >> 11) * (1.0 / 9007199254740992.0);
}

static void init_buf(Buf& B, size_t n)
{
    B.n = n;
    B.a = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.b = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.o = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.nrm = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.vel = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.vel0 = (fv_vec2*)xalloc(n * sizeof(fv_vec2));
    B.ax = (float*)xalloc(n * sizeof(float)); B.ay = (float*)xalloc(n * sizeof(float));
    B.bx = (float*)xalloc(n * sizeof(float)); B.by = (float*)xalloc(n * sizeof(float));
    B.ox = (float*)xalloc(n * sizeof(float)); B.oy = (float*)xalloc(n * sizeof(float));
    B.nx = (float*)xalloc(n * sizeof(float)); B.ny = (float*)xalloc(n * sizeof(float));
    B.nox = (float*)xalloc(n * sizeof(float)); B.noy = (float*)xalloc(n * sizeof(float));
    B.vx = (float*)xalloc(n * sizeof(float)); B.vy = (float*)xalloc(n * sizeof(float));
    B.vx0 = (float*)xalloc(n * sizeof(float)); B.vy0 = (float*)xalloc(n * sizeof(float));
    B.lenout = (float*)xalloc(n * sizeof(float));

    for (size_t i = 0; i < n; ++i) {
        float x = (float)(rnd01() * 4.0 - 2.0), y = (float)(rnd01() * 4.0 - 2.0);
        float u = (float)(rnd01() * 4.0 - 2.0), v = (float)(rnd01() * 4.0 - 2.0);
        B.a[i] = fv_vec2_make(x, y);
        B.b[i] = fv_vec2_make(u, v);
        B.ax[i] = x; B.ay[i] = y; B.bx[i] = u; B.by[i] = v;
        B.o[i] = fv_vec2_make(0.0f, 0.0f);
        /* 约 1.6% 的近零法线（走 len <= FV_EPS 分支） */
        float nxv = (float)(rnd01() * 2.0 - 1.0), nyv = (float)(rnd01() * 2.0 - 1.0);
        if ((i & 63u) == 0u) { nxv = 0.0f; nyv = 0.0f; }
        B.nrm[i] = fv_vec2_make(nxv, nyv);
        B.nx[i] = nxv; B.ny[i] = nyv;
        float vvx = (float)(rnd01() * 2.0 - 1.0), vvy = (float)(rnd01() * 2.0 - 1.0);
        B.vel[i] = fv_vec2_make(vvx, vvy);
        B.vx[i] = vvx; B.vy[i] = vvy;
    }
    std::memcpy(B.vel0, B.vel, n * sizeof(fv_vec2));
    std::memcpy(B.vx0, B.vx, n * sizeof(float));
    std::memcpy(B.vy0, B.vy, n * sizeof(float));
}

/* ------------------------------------------------------------------ 计时与校验和 */

static volatile double g_sink = 0.0;
static int g_cs_fail = 0, g_cs_total = 0;
static double g_ref_cs = 0.0, g_ref_tol = -1.0;
static double g_mhz = 0.0;

static uint64_t dbits(double d) { uint64_t u; std::memcpy(&u, &d, 8); return u; }

static void cs_check(const char* op, const char* impl, double cs, double ref, double tol_rel)
{
    g_cs_total++;
    if (dbits(cs) == dbits(ref)) return;
    double rel = std::fabs(cs - ref) / (std::fabs(ref) > 0.0 ? std::fabs(ref) : 1.0);
    if (tol_rel >= 0.0 && rel <= tol_rel) {
        std::printf("    [checksum] %s/%s not bit-identical, relative error %.3g ≤ %.3g (allowed by spec)\n", op, impl, rel, tol_rel);
        return;
    }
    std::printf("    [checksum] ✗ %s/%s cs=%.17g ref=%.17g rel=%.3g → this round's data is invalid\n", op, impl, cs, ref, rel);
    g_cs_fail++;
}

/* 归约 op：结果是一个累加标量 → 和式顺序可变，要求位级相等本身是错判据 */
static bool is_reduction_op(const char* op) { return std::strcmp(op, "dot") == 0; }

/* float 重结合误差量级 ~ eps·√n（随机符号下部分抵消），留 8 倍余量；
   逐元素 op 不用它（row() 对非归约 op 一律传 -1，仍要求位级相等） */
static double reduction_tol(size_t n) { return 8.0 * (double)FLT_EPSILON * std::sqrt((double)n) + 1e-9; }

template <class F>
static double time_op(size_t n, F&& f, double* cs_out)
{
    *cs_out = f();                                   /* 校验和＝单次调用结果 */
    double c0 = now_ns(); g_sink += f(); double c1 = now_ns();  /* warm-up + 校准：结果必须进 volatile sink，否则纯函数调用会被整个删掉，once≈0 → iters 撞上限 */
    double once = c1 - c0;
    if (once < 1.0) once = 1.0;
    double iters_d = 6.0e6 / once;                    /* ~6 ms 一个取样 */
    if (iters_d < 1.0) iters_d = 1.0;
    if (iters_d > 200000.0) iters_d = 200000.0;
    int iters = (int)iters_d;
    double acc = 0.0;
    double t[5];
    for (int r = 0; r < 5; ++r) {
        double a0 = now_ns();
        for (int k = 0; k < iters; ++k) acc += f();
        double a1 = now_ns();
        t[r] = (a1 - a0) / ((double)n * (double)iters);
    }
    std::sort(t, t + 5);
    g_sink += acc;
    return t[2];
}

static double cs_vec2(const fv_vec2* v, size_t n)
{
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += (double)v[i].x + (double)v[i].y;
    return s;
}
static double cs_soa(const float* x, const float* y, size_t n)
{
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += (double)x[i] + (double)y[i];
    return s;
}
static void zero_o(Buf& B) { std::memset(B.o, 0, B.n * sizeof(fv_vec2)); }

/* 一行结果：ns/elem、相对标量档加速比、校验和；并顺手做校验和一致性断言 */
static void row(const char* op, const char* impl, double ns, double cs, double base_ns, bool avail)
{
    if (!avail) { std::printf("  %-12s %-14s %10s %10s %18s\n", op, impl, "-", "-", "-"); return; }
    std::printf("  %-12s %-14s %10.3f %9.2fx %18.8g", op, impl, ns, base_ns / ns, cs);
    if (g_mhz > 0.0) std::printf(" %8.2f", ns * g_mhz / 1000.0);
    std::printf("\n");
    /* 容差只给归约 op；逐元素 op 即使 g_ref_tol 被设过也要求位级相等 */
    if (std::strcmp(impl, "scalar_api") != 0) cs_check(op, impl, cs, g_ref_cs, is_reduction_op(op) ? g_ref_tol : -1.0);
}

int main(int argc, char** argv)
{
    size_t sizes[2] = { 8192u, (size_t)1 << 20 };
    if (argc > 1 && std::strcmp(argv[1], "--small") == 0) sizes[1] = 65536u;
    g_mhz = cpu_mhz();
    std::printf("flashvec perf — independent benchmark (median of 5 × ~6 ms samples, CLOCK_MONOTONIC)\n");
    std::printf("build: %s %s\n", __DATE__, __TIME__);
    std::printf("starting mean CPU frequency: %.1f MHz (the cycle column converts through it; parallel builds and load shift both frequency and bandwidth)\n", g_mhz);

    for (int si = 0; si < 2; ++si) {
        size_t n = sizes[si];
        Buf B;
        init_buf(B, n);
        std::printf("\n=== N=%zu (%s, working set AoS a/b/o=%.1f MB) ===\n", n, si == 0 ? "L2-resident" : "bandwidth-bound", 3.0 * (double)n * 8.0 / 1048576.0);
        std::printf("  %-12s %-14s %10s %10s %18s %8s\n", "op", "impl", "ns/elem", "speedup", "checksum", "cycle/el");

        /* ---- add ---- */
        g_ref_tol = -1.0;
        double base_cs = 0.0, cs = 0.0;
        double base = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) B.o[i] = fv_vec2_add(B.a[i], B.b[i]); return cs_vec2(B.o, n); }, &base_cs);
        g_ref_cs = base_cs; row("add", "scalar_api", base, base_cs, base, true);
        double f1 = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i], y = B.b[i]; fv_vec2 r; r.x = x.x + y.x; r.y = x.y + y.y; B.o[i] = r; } return cs_vec2(B.o, n); }, &cs);
        row("add", "scalar_flat", f1, cs, base, true);
        double x1 = time_op(n, [&]{ fv_vec2_add_x4(B.a, B.b, B.o, n); return cs_vec2(B.o, n); }, &cs);
        row("add", "x4", x1, cs, base, true);
        double s1 = time_op(n, [&]{ fv_vec2_add_soa(B.ax, B.ay, B.bx, B.by, B.ox, B.oy, n); return cs_soa(B.ox, B.oy, n); }, &cs);
        row("add", "soa", s1, cs, base, true);

        /* ---- scale ---- */
        const fv_real sc = (fv_real)1.0009765625f;
        base = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) B.o[i] = fv_vec2_scale(B.a[i], sc); return cs_vec2(B.o, n); }, &base_cs);
        g_ref_cs = base_cs; row("scale", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i]; fv_vec2 r; r.x = x.x * sc; r.y = x.y * sc; B.o[i] = r; } return cs_vec2(B.o, n); }, &cs);
        row("scale", "scalar_flat", f1, cs, base, true);
        x1 = time_op(n, [&]{ fv_vec2_scale_x4(B.a, sc, B.o, n); return cs_vec2(B.o, n); }, &cs);
        row("scale", "x4", x1, cs, base, true);
        s1 = time_op(n, [&]{ fv_vec2_scale_soa(B.ax, B.ay, sc, B.ox, B.oy, n); return cs_soa(B.ox, B.oy, n); }, &cs);
        row("scale", "soa", s1, cs, base, true);

        /* ---- lerp ---- */
        const fv_real t = (fv_real)0.375f;
        base = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) B.o[i] = fv_vec2_lerp(B.a[i], B.b[i], t); return cs_vec2(B.o, n); }, &base_cs);
        g_ref_cs = base_cs; row("lerp", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i], y = B.b[i]; fv_vec2 r; r.x = x.x + (y.x - x.x) * t; r.y = x.y + (y.y - x.y) * t; B.o[i] = r; } return cs_vec2(B.o, n); }, &cs);
        row("lerp", "scalar_flat", f1, cs, base, true);
        x1 = time_op(n, [&]{ fv_vec2_lerp_x4(B.a, B.b, t, B.o, n); return cs_vec2(B.o, n); }, &cs);
        row("lerp", "x4", x1, cs, base, true);
        row("lerp", "soa(n/a)", 0, 0, base, false);

        /* ---- dot（归约 op：允许重排和式顺序 → 相对容差随 √n 增长） ---- */
        g_ref_tol = reduction_tol(n);
        base = time_op(n, [&]{ float s = 0; for (size_t i = 0; i < n; ++i) s += fv_vec2_dot(B.a[i], B.b[i]); return (double)s; }, &base_cs);
        g_ref_cs = base_cs; row("dot", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{ float s = 0; for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i], y = B.b[i]; s += x.x * y.x + x.y * y.y; } return (double)s; }, &cs);
        row("dot", "scalar_flat", f1, cs, base, true);
        row("dot", "x4(n/a)", 0, 0, base, false);
        s1 = time_op(n, [&]{ return (double)fv_vec2_dot_soa(B.ax, B.ay, B.bx, B.by, n); }, &cs);
        row("dot", "soa", s1, cs, base, true);
        g_ref_tol = -1.0;

        /* ---- len ---- */
        base = time_op(n, [&]{ double s = 0; for (size_t i = 0; i < n; ++i) s += (double)fv_vec2_len(B.a[i]); return s; }, &base_cs);
        g_ref_cs = base_cs; row("len", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{ double s = 0; for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i]; s += (double)std::sqrt(x.x * x.x + x.y * x.y); } return s; }, &cs);
        row("len", "scalar_flat", f1, cs, base, true);
        row("len", "x4(n/a)", 0, 0, base, false);
        s1 = time_op(n, [&]{ fv_vec2_len_soa(B.ax, B.ay, B.lenout, n); double s = 0; for (size_t i = 0; i < n; ++i) s += (double)B.lenout[i]; return s; }, &cs);
        row("len", "soa", s1, cs, base, true);

        /* ---- normalize ---- */
        base = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) B.o[i] = fv_vec2_normalize(B.a[i]); return cs_vec2(B.o, n); }, &base_cs);
        g_ref_cs = base_cs; row("normalize", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{ for (size_t i = 0; i < n; ++i) { fv_vec2 x = B.a[i]; float l = std::sqrt(x.x * x.x + x.y * x.y); fv_vec2 r; if (l <= FV_EPS) { r.x = 0; r.y = 0; } else { float inv = 1.0f / l; r.x = x.x * inv; r.y = x.y * inv; } B.o[i] = r; } return cs_vec2(B.o, n); }, &cs);
        row("normalize", "scalar_flat", f1, cs, base, true);
        row("normalize", "x4(n/a)", 0, 0, base, false);
        s1 = time_op(n, [&]{ fv_vec2_normalize_soa(B.ax, B.ay, B.ox, B.oy, n); return cs_soa(B.ox, B.oy, n); }, &cs);
        row("normalize", "soa", s1, cs, base, true);

        /* ---- collision 微内核：法线归一化 + 速度反射 v' = v - 2(v·n)n ---- */
        base = time_op(n, [&]{
            for (size_t i = 0; i < n; ++i) {
                fv_vec2 nn = fv_vec2_normalize(B.nrm[i]);
                fv_vec2 v = B.vel0[i];
                fv_real d = fv_vec2_dot(v, nn);
                B.vel[i] = fv_vec2_sub(v, fv_vec2_scale(nn, (fv_real)2 * d));
            }
            return cs_vec2(B.vel, n);
        }, &base_cs);
        g_ref_cs = base_cs; row("collision", "scalar_api", base, base_cs, base, true);
        f1 = time_op(n, [&]{
            for (size_t i = 0; i < n; ++i) {
                fv_vec2 x = B.nrm[i];
                float l = std::sqrt(x.x * x.x + x.y * x.y);
                float nxv, nyv;
                if (l <= FV_EPS) { nxv = 0; nyv = 0; } else { float inv = 1.0f / l; nxv = x.x * inv; nyv = x.y * inv; }
                float vx = B.vel0[i].x, vy = B.vel0[i].y;
                float d = vx * nxv + vy * nyv;
                B.vel[i] = fv_vec2_make(vx - nxv * (2.0f * d), vy - nyv * (2.0f * d));
            }
            return cs_vec2(B.vel, n);
        }, &cs);
        row("collision", "scalar_flat", f1, cs, base, true);
        /* SoA 混合档：归一化用库的 normalize_soa（能向量化）；反射阶段接口无法表达
           （逐元素点积 / 逐元素缩放没有批量形式），只能标量循环 → 如实标注 */
        s1 = time_op(n, [&]{
            fv_vec2_normalize_soa(B.nx, B.ny, B.nox, B.noy, n);
            std::memcpy(B.vx, B.vx0, n * sizeof(float));
            std::memcpy(B.vy, B.vy0, n * sizeof(float));
            for (size_t i = 0; i < n; ++i) {
                float d = B.vx[i] * B.nox[i] + B.vy[i] * B.noy[i];
                B.vx[i] = B.vx[i] - B.nox[i] * (2.0f * d);
                B.vy[i] = B.vy[i] - B.noy[i] * (2.0f * d);
            }
            double s = 0; for (size_t i = 0; i < n; ++i) s += (double)B.vx[i] + (double)B.vy[i];
            return s;
        }, &cs);
        row("collision", "soa_hybrid", s1, cs, base, true);
        row("collision", "x4(n/a)", 0, 0, base, false);

        /* ---- 稀疏：掩码 / 索引（两种实现的初始 o 都清零，校验和可直接比对） ---- */
        const size_t NS = n;
        uint8_t* mask = (uint8_t*)xalloc(NS);
        uint32_t* idx = (uint32_t*)xalloc(NS * sizeof(uint32_t));
        const size_t ratios[3] = { 5, 50, 95 };
        for (int r = 0; r < 3; ++r) {
            size_t pct = ratios[r];
            for (size_t i = 0; i < NS; ++i) mask[i] = (rnd01() * 100.0 < (double)pct) ? 0xff : 0x00;
            char op[24], impl[24];
            std::snprintf(op, sizeof op, "add_mask");
            zero_o(B);
            base = time_op(NS, [&]{ for (size_t i = 0; i < NS; ++i) B.o[i] = mask[i] ? fv_vec2_add(B.a[i], B.b[i]) : B.a[i]; return cs_vec2(B.o, NS); }, &base_cs);
            std::snprintf(impl, sizeof impl, "scalar(%zu%%)", pct);
            g_ref_cs = base_cs; row(op, impl, base, base_cs, base, true);
            zero_o(B);
            double m1 = time_op(NS, [&]{ fv_vec2_add_x4_mask(B.a, B.b, B.o, mask, NS); return cs_vec2(B.o, NS); }, &cs);
            std::snprintf(impl, sizeof impl, "x4_mask(%zu%%)", pct);
            row(op, impl, m1, cs, base, true);
        }
        for (int r = 0; r < 2; ++r) {
            size_t pct = ratios[r];
            size_t m = NS * pct / 100;
            for (size_t k = 0; k < m; ++k) idx[k] = (uint32_t)(rnd01() * (double)NS);
            char impl[24];
            zero_o(B);
            base = time_op(NS, [&]{ for (size_t k = 0; k < m; ++k) { uint32_t i = idx[k]; B.o[i] = fv_vec2_add(B.a[i], B.b[i]); } return cs_vec2(B.o, NS); }, &base_cs);
            std::snprintf(impl, sizeof impl, "scalar(%zu%%)", pct);
            g_ref_cs = base_cs; row("add_idx", impl, base, base_cs, base, true);
            zero_o(B);
            double g1 = time_op(NS, [&]{ fv_vec2_add_x4_idx(B.a, B.b, B.o, idx, m); return cs_vec2(B.o, NS); }, &cs);
            std::snprintf(impl, sizeof impl, "x4_idx(%zu%%)", pct);
            row("add_idx", impl, g1, cs, base, true);
        }
        std::free(mask); std::free(idx);
        std::printf("  (ns/elem is divided by N; sparse speedup is relative to the scalar loop at the same density)\n");
    }
    std::printf("\nchecksum consistency: %d/%d passed (reduction op dot: relative tolerance 8*eps*sqrt(n); everything else must be bit-identical)\n",
                g_cs_total - g_cs_fail, g_cs_total);
    std::printf("ending mean CPU frequency: %.1f MHz\n", cpu_mhz());
    std::printf("verdict: %s\n", g_cs_fail ? "checksum mismatch → part of the table above is invalid" : "checksums consistent → all variants measure the same work");
    return g_cs_fail ? 1 : 0;
}
