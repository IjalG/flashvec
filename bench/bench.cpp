/*
 * bench/bench.cpp -- flashvec 三档实现基准（档1 标量 AoS / 档2 自动向量化 SoA / 档3 手写 AVX2 SoA）
 *
 * 口径：
 *   - 每个用例固定重复 FV_REPS 次取一次平均（默认 100），重复 FV_TRIALS 次取中位数（默认 5）；
 *   - ns_elem = 中位数单次遍历耗时 / 元素数；
 *   - 每个用例在重复循环之间插入 asm 内存屏障，防止编译器把不变量提到循环外；
 *   - 计时区外单独跑一次校验和（读结果数组）以防死代码消除；
 *   - 精度对照用的参考实现在 ref.cpp（严格 IEEE，独立编译单元）。
 *
 * 编译与运行见 bench/Makefile。
 */
#include "flashvec.h"
#include "ref.h"

#include <immintrin.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef FV_N
#define FV_N (1u << 20)
#endif
#ifndef FV_NSMALL
#define FV_NSMALL 8192u
#endif
#ifndef FV_NMID
#define FV_NMID (1u << 16)
#endif
#ifndef FV_REPS
#define FV_REPS 100u
#endif
#ifndef FV_REPS_SMALL
#define FV_REPS_SMALL 4000u
#endif
#ifndef FV_REPS_MID
#define FV_REPS_MID 400u
#endif
#ifndef FV_TRIALS
#define FV_TRIALS 5
#endif
#ifndef FV_CHECKS
#define FV_CHECKS FV_NSMALL
#endif
#ifndef FV_CFG_NAME
#define FV_CFG_NAME "unknown"
#endif

#define FV_ALIGN alignas(32)
#define FV_NOINLINE __attribute__((noinline))
#define FV_FWD(x) static_cast<decltype(x)>(x)

static volatile float g_sink;
static float g_dot_scalar;
#define CLOBBER() __asm__ __volatile__("" ::: "memory")

/* ---------------- 数据（AoS 与 SoA 由同一批随机数生成） ---------------- */

FV_ALIGN static fv_vec2 g_a[FV_N], g_b[FV_N], g_v[FV_N], g_o[FV_N];
FV_ALIGN static float g_ax[FV_N], g_ay[FV_N], g_bx[FV_N], g_by[FV_N];
FV_ALIGN static float g_vx[FV_N], g_vy[FV_N], g_ox[FV_N], g_oy[FV_N];
FV_ALIGN static float g_rx[FV_N], g_ry[FV_N];

static uint32_t g_rng = 0x12345678u;

static float frand(float lo, float hi)
{
    g_rng = g_rng * 1664525u + 1013904223u;
    const float u = (float)((g_rng >> 8) & 0xFFFFFFu) / (float)0x1000000u;
    return lo + (hi - lo) * u;
}

static void init_data(void)
{
    for (unsigned i = 0; i < FV_N; ++i) {
        const float ax = frand(-1.0f, 1.0f), ay = frand(-1.0f, 1.0f);
        const float bx = frand(-1.0f, 1.0f), by = frand(-1.0f, 1.0f);
        const float vx = frand(-2.0f, 2.0f), vy = frand(-2.0f, 2.0f);
        g_a[i].x = ax; g_a[i].y = ay;
        g_b[i].x = bx; g_b[i].y = by;
        g_v[i].x = vx; g_v[i].y = vy;
        g_ax[i] = ax; g_ay[i] = ay;
        g_bx[i] = bx; g_by[i] = by;
        g_vx[i] = vx; g_vy[i] = vy;
    }
    /* 每 4096 个元素放一个零向量 / 一对重合点，专门踩 normalize 的 FV_EPS 分支。 */
    for (unsigned i = 0; i < FV_N; i += 4096u) {
        g_a[i].x = 0.0f; g_a[i].y = 0.0f;
        g_ax[i] = 0.0f;  g_ay[i] = 0.0f;
        g_b[i].x = 0.0f; g_b[i].y = 0.0f;
        g_bx[i] = 0.0f;  g_by[i] = 0.0f;
    }
}

static void aos_to_soa(const fv_vec2 *v, float *x, float *y, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) { x[i] = v[i].x; y[i] = v[i].y; }
}

static float sum_f(const float *p, unsigned n)
{
    float s = 0.0f;
    for (unsigned i = 0; i < n; ++i) { s += p[i]; }
    return s;
}

static float sum_aos(const fv_vec2 *v, unsigned n)
{
    float s = 0.0f;
    for (unsigned i = 0; i < n; ++i) { s += v[i].x + v[i].y; }
    return s;
}

/* ---------------- 计时 / 环境 ---------------- */

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static double median_of(const double *v, int n)
{
    double t[16];
    for (int i = 0; i < n; ++i) { t[i] = v[i]; }
    for (int i = 1; i < n; ++i) {
        const double k = t[i];
        int j = i - 1;
        while (j >= 0 && t[j] > k) { t[j + 1] = t[j]; --j; }
        t[j + 1] = k;
    }
    return (n & 1) ? t[n / 2] : 0.5 * (t[n / 2 - 1] + t[n / 2]);
}

static long read_cpu_khz(int cpu)
{
    char path[160];
    snprintf(path, sizeof path, "/sys/devices/system/cpu/cpu%d/cpufreq/scaling_cur_freq", cpu);
    FILE *f = fopen(path, "r");
    if (!f) { return -1; }
    long khz = -1;
    if (fscanf(f, "%ld", &khz) != 1) { khz = -1; }
    fclose(f);
    return khz;
}

static void pin_cpu(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (sched_setaffinity(0, sizeof set, &set) == 0) {
        printf("# pinned cpu=%d running_on=%d\n", cpu, sched_getcpu());
    } else {
        printf("# pin failed cpu=%d (continuing unpinned)\n", cpu);
    }
}

static uint64_t rdtsc_now(void)
{
    unsigned aux;
    return __rdtscp(&aux);
}

static double calibrate_tsc_hz(void)
{
    const double t0 = now_s();
    const uint64_t c0 = rdtsc_now();
    double t1 = t0;
    while (t1 - t0 < 0.05) { t1 = now_s(); }
    const uint64_t c1 = rdtsc_now();
    return (double)(c1 - c0) / (t1 - t0);
}

/* ---------------- 档1：标量 AoS，调用库函数（跨编译单元，不会被内联） ---------------- */

FV_NOINLINE static void t1_add(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) { g_o[i] = fv_vec2_add(g_a[i], g_b[i]); }
}

FV_NOINLINE static void t1_scale(unsigned n, float s)
{
    for (unsigned i = 0; i < n; ++i) { g_o[i] = fv_vec2_scale(g_a[i], s); }
}

FV_NOINLINE static void t1_dot(unsigned n)
{
    float s = 0.0f;
    for (unsigned i = 0; i < n; ++i) { s += fv_vec2_dot(g_a[i], g_b[i]); }
    g_dot_scalar = s;
}

/* 逐元素点积：库没有 fv_vec2 逐分量乘，标量档只能用 fv_vec2_dot 逐元素调一次。 */
FV_NOINLINE static void t1_dot_elem(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float d = fv_vec2_dot(g_a[i], g_b[i]);
        g_o[i].x = d;
        g_o[i].y = d;
    }
}

FV_NOINLINE static void t1_len(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        g_o[i].x = fv_vec2_len(g_a[i]);
        g_o[i].y = 0.0f;
    }
}

FV_NOINLINE static void t1_normalize(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) { g_o[i] = fv_vec2_normalize(g_a[i]); }
}

FV_NOINLINE static void t1_collision(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const fv_vec2 d = fv_vec2_sub(g_a[i], g_b[i]);
        const fv_vec2 nv = fv_vec2_normalize(d);
        const float dn = fv_vec2_dot(g_v[i], nv);
        g_o[i] = fv_vec2_sub(g_v[i], fv_vec2_scale(nv, 2.0f * dn));
    }
}

/* ---------------- 档2：SoA + 最朴素 for 循环，靠编译器自动向量化 ---------------- */

FV_NOINLINE static void t2_add(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        g_ox[i] = g_ax[i] + g_bx[i];
        g_oy[i] = g_ay[i] + g_by[i];
    }
}

FV_NOINLINE static void t2_scale(unsigned n, float s)
{
    for (unsigned i = 0; i < n; ++i) {
        g_ox[i] = g_ax[i] * s;
        g_oy[i] = g_ay[i] * s;
    }
}

FV_NOINLINE static void t2_dot(unsigned n)
{
    float s = 0.0f;
    for (unsigned i = 0; i < n; ++i) {
        s += g_ax[i] * g_bx[i] + g_ay[i] * g_by[i];
    }
    g_dot_scalar = s;
}

FV_NOINLINE static void t2_dot_elem(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float d = g_ax[i] * g_bx[i] + g_ay[i] * g_by[i];
        g_ox[i] = d;
        g_oy[i] = d;
    }
}

FV_NOINLINE static void t2_len(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float l = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i]);
        g_ox[i] = (l <= FV_EPS) ? 0.0f : l;
        g_oy[i] = 0.0f;
    }
}

FV_NOINLINE static void t2_normalize(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float l = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i]);
        const float inv = (l <= FV_EPS) ? 0.0f : 1.0f / l;
        g_ox[i] = g_ax[i] * inv;
        g_oy[i] = g_ay[i] * inv;
    }
}

FV_NOINLINE static void t2_collision(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float dx = g_ax[i] - g_bx[i];
        const float dy = g_ay[i] - g_by[i];
        const float l = sqrtf(dx * dx + dy * dy);
        float nx = 0.0f, ny = 0.0f;
        if (l > FV_EPS) {
            const float inv = 1.0f / l;
            nx = dx * inv;
            ny = dy * inv;
        }
        const float dn = g_vx[i] * nx + g_vy[i] * ny;
        g_ox[i] = g_vx[i] - 2.0f * dn * nx;
        g_oy[i] = g_vy[i] - 2.0f * dn * ny;
    }
}

/* ---------------- 档3：手写 AVX2，SoA，8 float/次 ---------------- */

#if defined(__AVX2__)
#define FV_HAVE_AVX2 1

/* 跨 lane 归约：用 extract + shuffle，不用 vhaddps。 */
static inline float hsum_shuffle(__m256 v)
{
    const __m128 lo = _mm256_castps256_ps128(v);
    const __m128 hi = _mm256_extractf128_ps(v, 1);
    __m128 s = _mm_add_ps(lo, hi);
    s = _mm_add_ps(s, _mm_movehl_ps(s, s));
    s = _mm_add_ss(s, _mm_shuffle_ps(s, s, 0x55));
    return _mm_cvtss_f32(s);
}

/* 对照组：vhaddps 归约，用于量化 hadd 的代价。 */
static inline float hsum_hadd(__m256 v)
{
    const __m128 lo = _mm256_castps256_ps128(v);
    const __m128 hi = _mm256_extractf128_ps(v, 1);
    __m128 s = _mm_hadd_ps(lo, hi);  /* [a0+a1, a2+a3, a4+a5, a6+a7] */
    s = _mm_hadd_ps(s, s);           /* [(a01+a23), (a45+a67), ...] */
    s = _mm_hadd_ps(s, s);           /* 4 个 lane 全部求和 */
    return _mm_cvtss_f32(s);
}

FV_NOINLINE static void t3_add(unsigned n)
{
    const unsigned n8 = n & ~7u;
    for (unsigned i = 0; i < n8; i += 8) {
        _mm256_store_ps(g_ox + i, _mm256_add_ps(_mm256_load_ps(g_ax + i), _mm256_load_ps(g_bx + i)));
        _mm256_store_ps(g_oy + i, _mm256_add_ps(_mm256_load_ps(g_ay + i), _mm256_load_ps(g_by + i)));
    }
    for (unsigned i = n8; i < n; ++i) { g_ox[i] = g_ax[i] + g_bx[i]; g_oy[i] = g_ay[i] + g_by[i]; }
}

FV_NOINLINE static void t3_scale(unsigned n, float s)
{
    const unsigned n8 = n & ~7u;
    const __m256 vs = _mm256_set1_ps(s);
    for (unsigned i = 0; i < n8; i += 8) {
        _mm256_store_ps(g_ox + i, _mm256_mul_ps(_mm256_load_ps(g_ax + i), vs));
        _mm256_store_ps(g_oy + i, _mm256_mul_ps(_mm256_load_ps(g_ay + i), vs));
    }
    for (unsigned i = n8; i < n; ++i) { g_ox[i] = g_ax[i] * s; g_oy[i] = g_ay[i] * s; }
}

FV_NOINLINE static void t3_dot(unsigned n)
{
    const unsigned n8 = n & ~7u;
    __m256 acc = _mm256_setzero_ps();
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 bx = _mm256_load_ps(g_bx + i);
        const __m256 by = _mm256_load_ps(g_by + i);
        acc = _mm256_fmadd_ps(ax, bx, acc);
        acc = _mm256_fmadd_ps(ay, by, acc);
    }
    float s = hsum_shuffle(acc);
    for (unsigned i = n8; i < n; ++i) { s += g_ax[i] * g_bx[i] + g_ay[i] * g_by[i]; }
    g_dot_scalar = s;
}

/* 4 个独立累加器版本：每轮 16 个元素，专治单累加器的 FMA 依赖链。 */
FV_NOINLINE static void t3_dot4(unsigned n)
{
    const unsigned n16 = n & ~15u;
    __m256 a0 = _mm256_setzero_ps(), a1 = _mm256_setzero_ps();
    __m256 a2 = _mm256_setzero_ps(), a3 = _mm256_setzero_ps();
    for (unsigned i = 0; i < n16; i += 16) {
        const __m256 ax0 = _mm256_load_ps(g_ax + i);
        const __m256 ay0 = _mm256_load_ps(g_ay + i);
        const __m256 bx0 = _mm256_load_ps(g_bx + i);
        const __m256 by0 = _mm256_load_ps(g_by + i);
        const __m256 ax1 = _mm256_load_ps(g_ax + i + 8);
        const __m256 ay1 = _mm256_load_ps(g_ay + i + 8);
        const __m256 bx1 = _mm256_load_ps(g_bx + i + 8);
        const __m256 by1 = _mm256_load_ps(g_by + i + 8);
        a0 = _mm256_fmadd_ps(ax0, bx0, a0);
        a1 = _mm256_fmadd_ps(ay0, by0, a1);
        a2 = _mm256_fmadd_ps(ax1, bx1, a2);
        a3 = _mm256_fmadd_ps(ay1, by1, a3);
    }
    float s = hsum_shuffle(_mm256_add_ps(_mm256_add_ps(a0, a1), _mm256_add_ps(a2, a3)));
    for (unsigned i = n16; i < n; ++i) { s += g_ax[i] * g_bx[i] + g_ay[i] * g_by[i]; }
    g_dot_scalar = s;
}

/* 对照组：每 8 个元素就 hadd 归约一次（“naive SIMD”写法）。 */
FV_NOINLINE static void t3_dot_hadd(unsigned n)
{
    const unsigned n8 = n & ~7u;
    float s = 0.0f;
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 bx = _mm256_load_ps(g_bx + i);
        const __m256 by = _mm256_load_ps(g_by + i);
        __m256 m = _mm256_mul_ps(ax, bx);
        m = _mm256_fmadd_ps(ay, by, m);
        s += hsum_hadd(m);
    }
    for (unsigned i = n8; i < n; ++i) { s += g_ax[i] * g_bx[i] + g_ay[i] * g_by[i]; }
    g_dot_scalar = s;
}

FV_NOINLINE static void t3_dot_elem(unsigned n)
{
    const unsigned n8 = n & ~7u;
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 m = _mm256_fmadd_ps(ax, _mm256_load_ps(g_bx + i),
                                        _mm256_mul_ps(ay, _mm256_load_ps(g_by + i)));
        _mm256_store_ps(g_ox + i, m);
        _mm256_store_ps(g_oy + i, m);
    }
    for (unsigned i = n8; i < n; ++i) {
        const float d = g_ax[i] * g_bx[i] + g_ay[i] * g_by[i];
        g_ox[i] = d;
        g_oy[i] = d;
    }
}

FV_NOINLINE static void t3_len(unsigned n)
{
    const unsigned n8 = n & ~7u;
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 l2 = _mm256_fmadd_ps(ax, ax, _mm256_mul_ps(ay, ay));
        const __m256 l = _mm256_sqrt_ps(l2);
        const __m256 keep = _mm256_cmp_ps(l, eps, _CMP_GT_OQ);
        _mm256_store_ps(g_ox + i, _mm256_and_ps(l, keep));
        _mm256_store_ps(g_oy + i, _mm256_setzero_ps());
    }
    for (unsigned i = n8; i < n; ++i) {
        const float l = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i]);
        g_ox[i] = (l <= FV_EPS) ? 0.0f : l;
        g_oy[i] = 0.0f;
    }
}

/* sqrt + div 路径：与标量库逐位一致的候选实现。 */
FV_NOINLINE static void t3_normalize(unsigned n)
{
    const unsigned n8 = n & ~7u;
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 l2 = _mm256_fmadd_ps(ax, ax, _mm256_mul_ps(ay, ay));
        const __m256 l = _mm256_sqrt_ps(l2);
        const __m256 keep = _mm256_cmp_ps(l, eps, _CMP_GT_OQ);
        const __m256 inv = _mm256_div_ps(one, l);
        _mm256_store_ps(g_ox + i, _mm256_and_ps(_mm256_mul_ps(ax, inv), keep));
        _mm256_store_ps(g_oy + i, _mm256_and_ps(_mm256_mul_ps(ay, inv), keep));
    }
    for (unsigned i = n8; i < n; ++i) {
        const float l = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i]);
        const float inv = (l <= FV_EPS) ? 0.0f : 1.0f / l;
        g_ox[i] = g_ax[i] * inv;
        g_oy[i] = g_ay[i] * inv;
    }
}

/* rsqrt + 一次牛顿迭代路径：精度换速度的候选实现。 */
FV_NOINLINE static void t3_normalize_rsqrt(unsigned n)
{
    const unsigned n8 = n & ~7u;
    const __m256 half = _mm256_set1_ps(0.5f);
    const __m256 c15 = _mm256_set1_ps(1.5f);
    const __m256 eps2 = _mm256_set1_ps(FV_EPS * FV_EPS);
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 ax = _mm256_load_ps(g_ax + i);
        const __m256 ay = _mm256_load_ps(g_ay + i);
        const __m256 l2 = _mm256_fmadd_ps(ax, ax, _mm256_mul_ps(ay, ay));
        __m256 r = _mm256_rsqrt_ps(l2);
        const __m256 t = _mm256_mul_ps(_mm256_mul_ps(l2, r), r);
        r = _mm256_mul_ps(r, _mm256_fnmadd_ps(t, half, c15));
        const __m256 keep = _mm256_cmp_ps(l2, eps2, _CMP_GT_OQ);
        _mm256_store_ps(g_ox + i, _mm256_and_ps(_mm256_mul_ps(ax, r), keep));
        _mm256_store_ps(g_oy + i, _mm256_and_ps(_mm256_mul_ps(ay, r), keep));
    }
    for (unsigned i = n8; i < n; ++i) {
        const float l = sqrtf(g_ax[i] * g_ax[i] + g_ay[i] * g_ay[i]);
        const float inv = (l <= FV_EPS) ? 0.0f : 1.0f / l;
        g_ox[i] = g_ax[i] * inv;
        g_oy[i] = g_ay[i] * inv;
    }
}

FV_NOINLINE static void t3_collision(unsigned n)
{
    const unsigned n8 = n & ~7u;
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 two = _mm256_set1_ps(2.0f);
    const __m256 one = _mm256_set1_ps(1.0f);
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 dx = _mm256_sub_ps(_mm256_load_ps(g_ax + i), _mm256_load_ps(g_bx + i));
        const __m256 dy = _mm256_sub_ps(_mm256_load_ps(g_ay + i), _mm256_load_ps(g_by + i));
        const __m256 l2 = _mm256_fmadd_ps(dx, dx, _mm256_mul_ps(dy, dy));
        const __m256 l = _mm256_sqrt_ps(l2);
        const __m256 keep = _mm256_cmp_ps(l, eps, _CMP_GT_OQ);
        const __m256 inv = _mm256_and_ps(_mm256_div_ps(one, l), keep);
        const __m256 nx = _mm256_mul_ps(dx, inv);
        const __m256 ny = _mm256_mul_ps(dy, inv);
        const __m256 vx = _mm256_load_ps(g_vx + i);
        const __m256 vy = _mm256_load_ps(g_vy + i);
        const __m256 dn = _mm256_fmadd_ps(vx, nx, _mm256_mul_ps(vy, ny));
        const __m256 k = _mm256_mul_ps(two, dn);
        _mm256_store_ps(g_ox + i, _mm256_fnmadd_ps(k, nx, vx));
        _mm256_store_ps(g_oy + i, _mm256_fnmadd_ps(k, ny, vy));
    }
    for (unsigned i = n8; i < n; ++i) {
        const float dx = g_ax[i] - g_bx[i], dy = g_ay[i] - g_by[i];
        const float l = sqrtf(dx * dx + dy * dy);
        float nx = 0.0f, ny = 0.0f;
        if (l > FV_EPS) { const float inv = 1.0f / l; nx = dx * inv; ny = dy * inv; }
        const float dn = g_vx[i] * nx + g_vy[i] * ny;
        g_ox[i] = g_vx[i] - 2.0f * dn * nx;
        g_oy[i] = g_vy[i] - 2.0f * dn * ny;
    }
}

/* rsqrt 版碰撞微内核：法线用 rsqrt+NR，其余同 AVX2 版。 */
FV_NOINLINE static void t3_collision_rsqrt(unsigned n)
{
    const unsigned n8 = n & ~7u;
    const __m256 half = _mm256_set1_ps(0.5f);
    const __m256 c15 = _mm256_set1_ps(1.5f);
    const __m256 eps2 = _mm256_set1_ps(FV_EPS * FV_EPS);
    const __m256 two = _mm256_set1_ps(2.0f);
    for (unsigned i = 0; i < n8; i += 8) {
        const __m256 dx = _mm256_sub_ps(_mm256_load_ps(g_ax + i), _mm256_load_ps(g_bx + i));
        const __m256 dy = _mm256_sub_ps(_mm256_load_ps(g_ay + i), _mm256_load_ps(g_by + i));
        const __m256 l2 = _mm256_fmadd_ps(dx, dx, _mm256_mul_ps(dy, dy));
        __m256 r = _mm256_rsqrt_ps(l2);
        const __m256 t = _mm256_mul_ps(_mm256_mul_ps(l2, r), r);
        r = _mm256_mul_ps(r, _mm256_fnmadd_ps(t, half, c15));
        const __m256 keep = _mm256_cmp_ps(l2, eps2, _CMP_GT_OQ);
        const __m256 nx = _mm256_and_ps(_mm256_mul_ps(dx, r), keep);
        const __m256 ny = _mm256_and_ps(_mm256_mul_ps(dy, r), keep);
        const __m256 vx = _mm256_load_ps(g_vx + i);
        const __m256 vy = _mm256_load_ps(g_vy + i);
        const __m256 dn = _mm256_fmadd_ps(vx, nx, _mm256_mul_ps(vy, ny));
        const __m256 k = _mm256_mul_ps(two, dn);
        _mm256_store_ps(g_ox + i, _mm256_fnmadd_ps(k, nx, vx));
        _mm256_store_ps(g_oy + i, _mm256_fnmadd_ps(k, ny, vy));
    }
    for (unsigned i = n8; i < n; ++i) {
        const float dx = g_ax[i] - g_bx[i], dy = g_ay[i] - g_by[i];
        const float l = sqrtf(dx * dx + dy * dy);
        float nx = 0.0f, ny = 0.0f;
        if (l > FV_EPS) { const float inv = 1.0f / l; nx = dx * inv; ny = dy * inv; }
        const float dn = g_vx[i] * nx + g_vy[i] * ny;
        g_ox[i] = g_vx[i] - 2.0f * dn * nx;
        g_oy[i] = g_vy[i] - 2.0f * dn * ny;
    }
}
#endif /* __AVX2__ */

/* ---------------- 计时主体 ---------------- */

template <typename K, typename C>
static void run_case(const char *tier, const char *op, unsigned n, unsigned reps,
                     double bytes_per_elem, K &&kernel, C &&checksum)
{
    kernel(n); /* 预热 + 迫使页表/Cache 就位 */
    CLOBBER();

    double trials[FV_TRIALS];
    float chk = 0.0f;
    for (int t = 0; t < FV_TRIALS; ++t) {
        CLOBBER();
        const double t0 = now_s();
        for (unsigned r = 0; r < reps; ++r) {
            kernel(n);
            CLOBBER();
        }
        const double t1 = now_s();
        chk += checksum(n);
        trials[t] = (t1 - t0) / (double)reps * 1e9;
    }
    g_sink = chk;

    const double med = median_of(trials, FV_TRIALS);
    double tmin = trials[0];
    for (int t = 1; t < FV_TRIALS; ++t) { if (trials[t] < tmin) { tmin = trials[t]; } }
    const double ns_elem = med / (double)n;
    const double gbps = (bytes_per_elem * (double)n) / med;
    printf("CASE cfg=%s tier=%s op=%s n=%u reps=%u med_ns=%.1f min_ns=%.1f ns_elem=%.4f "
           "ns_elem_min=%.4f gbps=%.2f chk=%.6g trials_ns=",
           FV_CFG_NAME, tier, op, n, reps, med, tmin, ns_elem, tmin / (double)n, gbps, (double)chk);
    for (int t = 0; t < FV_TRIALS; ++t) { printf(" %.1f", trials[t]); }
    printf("\n");
    fflush(stdout);
}

static unsigned reps_for(unsigned n)
{
    if (n >= 1000000u) { return FV_REPS; }
    if (n >= 65536u) { return FV_REPS_MID; }
    return FV_REPS_SMALL;
}

/* 打印 /proc/loadavg，用来记录测量时的机器负载（本机是共享环境）。 */
static void print_load(const char *tag)
{
    char buf[160] = "unavailable";
    FILE *f = fopen("/proc/loadavg", "r");
    if (f) { if (fgets(buf, (int)sizeof buf, f) == NULL) { } fclose(f); }
    printf("# loadavg_%s=%s", tag, buf);
    fflush(stdout);
}

/* ---------------- 精度对照 ---------------- */

FV_ALIGN static float g_sx[FV_N], g_sy[FV_N];

/* 存下 tier1 的输出，供后续档位做逐位比较。 */
static void save_out(unsigned n)
{
    memcpy(g_sx, g_ox, n * sizeof(float));
    memcpy(g_sy, g_oy, n * sizeof(float));
}

/* 与 tier1 输出逐位比较：bits_equal=1 表示 IEEE-754 位模式完全一致。 */
static void cmp_bits(const char *tier, const char *op, unsigned n)
{
    const int bx = (memcmp(g_sx, g_ox, n * sizeof(float)) == 0);
    const int by = (memcmp(g_sy, g_oy, n * sizeof(float)) == 0);
    unsigned nd = 0;
    for (unsigned i = 0; i < n; ++i) {
        if (g_sx[i] != g_ox[i] || g_sy[i] != g_oy[i]) { ++nd; }
    }
    printf("IDENT cfg=%s tier=%s op=%s n=%u bits_equal=%d n_diff=%u\n", FV_CFG_NAME, tier, op, n,
           (bx && by) ? 1 : 0, nd);
}

static void diff_report(const char *tier, const char *op, const float *gx, const float *gy,
                        const float *rx, const float *ry, unsigned n)
{
    double max_abs = 0.0, max_rel = 0.0;
    unsigned n_diff = 0, n_nan = 0;
    for (unsigned i = 0; i < n; ++i) {
        const double gxa = gx[i], gya = gy[i];
        if (isnan(gxa) || isnan(gya)) { ++n_nan; continue; }
        const double dx = fabs(gxa - (double)rx[i]);
        const double dy = fabs(gya - (double)ry[i]);
        const double a = dx > dy ? dx : dy;
        const double denom = fmax(fabs((double)rx[i]), fabs((double)ry[i]));
        const double rel = (denom > 1e-9) ? a / denom : a;
        if (a > max_abs) { max_abs = a; }
        if (rel > max_rel) { max_rel = rel; }
        if (a != 0.0) { ++n_diff; }
    }
    printf("DIFF cfg=%s tier=%s op=%s n=%u max_abs=%.3e max_rel=%.3e n_diff=%u n_nan=%u\n",
           FV_CFG_NAME, tier, op, n, max_abs, max_rel, n_diff, n_nan);
    fflush(stdout);
}

static void run_checks(unsigned n)
{
    /* add */
    ref_add(g_ax, g_ay, g_bx, g_by, g_rx, g_ry, n);
    t1_add(n); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "add", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_add(n); diff_report("tier2", "add", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "add", n);
#if defined(FV_HAVE_AVX2)
    t3_add(n); diff_report("tier3", "add", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "add", n);
#endif

    /* scale */
    const float s = 0.5f;
    ref_scale(g_ax, g_ay, s, g_rx, g_ry, n);
    t1_scale(n, s); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "scale", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_scale(n, s); diff_report("tier2", "scale", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "scale", n);
#if defined(FV_HAVE_AVX2)
    t3_scale(n, s); diff_report("tier3", "scale", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "scale", n);
#endif

    /* dot（标量归约），参考值单独比较 */
    const float rdot = ref_dot_sum(g_ax, g_ay, g_bx, g_by, n);
    t1_dot(n);
    printf("DIFF cfg=%s tier=tier1 op=dot n=%u scalar=%.9g ref=%.9g rel=%.3e\n", FV_CFG_NAME, n,
           (double)g_dot_scalar, (double)rdot,
           fabs((double)g_dot_scalar - (double)rdot) / fabs((double)rdot));
    t2_dot(n);
    printf("DIFF cfg=%s tier=tier2 op=dot n=%u scalar=%.9g ref=%.9g rel=%.3e\n", FV_CFG_NAME, n,
           (double)g_dot_scalar, (double)rdot,
           fabs((double)g_dot_scalar - (double)rdot) / fabs((double)rdot));
#if defined(FV_HAVE_AVX2)
    t3_dot(n);
    printf("DIFF cfg=%s tier=tier3 op=dot n=%u scalar=%.9g ref=%.9g rel=%.3e\n", FV_CFG_NAME, n,
           (double)g_dot_scalar, (double)rdot,
           fabs((double)g_dot_scalar - (double)rdot) / fabs((double)rdot));
    t3_dot_hadd(n);
    printf("DIFF cfg=%s tier=tier3hadd op=dot n=%u scalar=%.9g ref=%.9g rel=%.3e\n", FV_CFG_NAME, n,
           (double)g_dot_scalar, (double)rdot,
           fabs((double)g_dot_scalar - (double)rdot) / fabs((double)rdot));
    t3_dot4(n);
    printf("DIFF cfg=%s tier=tier3_4acc op=dot n=%u scalar=%.9g ref=%.9g rel=%.3e\n", FV_CFG_NAME, n,
           (double)g_dot_scalar, (double)rdot,
           fabs((double)g_dot_scalar - (double)rdot) / fabs((double)rdot));
#endif

    /* dot_elem */
    ref_dot_elem(g_ax, g_ay, g_bx, g_by, g_rx, g_ry, n);
    t1_dot_elem(n); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "dot_elem", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_dot_elem(n); diff_report("tier2", "dot_elem", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "dot_elem", n);
#if defined(FV_HAVE_AVX2)
    t3_dot_elem(n); diff_report("tier3", "dot_elem", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "dot_elem", n);
#endif

    /* len */
    ref_len(g_ax, g_ay, g_rx, g_ry, n);
    t1_len(n); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "len", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_len(n); diff_report("tier2", "len", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "len", n);
#if defined(FV_HAVE_AVX2)
    t3_len(n); diff_report("tier3", "len", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "len", n);
#endif

    /* normalize */
    ref_normalize(g_ax, g_ay, g_rx, g_ry, n);
    t1_normalize(n); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "normalize", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_normalize(n); diff_report("tier2", "normalize", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "normalize", n);
#if defined(FV_HAVE_AVX2)
    t3_normalize(n); diff_report("tier3", "normalize", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "normalize", n);
    t3_normalize_rsqrt(n); diff_report("tier3_rsqrt", "normalize", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3_rsqrt", "normalize", n);
#endif

    /* collision */
    ref_collision(g_ax, g_ay, g_bx, g_by, g_vx, g_vy, g_rx, g_ry, n);
    t1_collision(n); aos_to_soa(g_o, g_ox, g_oy, n);
    diff_report("tier1", "collision", g_ox, g_oy, g_rx, g_ry, n); save_out(n);
    t2_collision(n); diff_report("tier2", "collision", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier2", "collision", n);
#if defined(FV_HAVE_AVX2)
    t3_collision(n); diff_report("tier3", "collision", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3", "collision", n);
    t3_collision_rsqrt(n); diff_report("tier3_rsqrt", "collision", g_ox, g_oy, g_rx, g_ry, n); cmp_bits("tier3_rsqrt", "collision", n);
#endif
    fflush(stdout);
}

/* ---------------- 各档计时入口 ---------------- */

static void bench_all(const char *tier, unsigned n)
{
    const unsigned reps = reps_for(n);
    const float s = 0.5f;
    if (strcmp(tier, "tier1") == 0 || strcmp(tier, "all") == 0) {
        run_case("tier1", "add", n, reps, 24.0,
                 [](unsigned m) { t1_add(m); }, [](unsigned m) { return sum_aos(g_o, m); });
        run_case("tier1", "scale", n, reps, 16.0,
                 [s](unsigned m) { t1_scale(m, s); }, [](unsigned m) { return sum_aos(g_o, m); });
        run_case("tier1", "dot", n, reps, 16.0,
                 [](unsigned m) { t1_dot(m); }, [](unsigned) { return g_dot_scalar; });
        run_case("tier1", "dot_elem", n, reps, 20.0,
                 [](unsigned m) { t1_dot_elem(m); }, [](unsigned m) { return sum_aos(g_o, m); });
        run_case("tier1", "len", n, reps, 20.0,
                 [](unsigned m) { t1_len(m); }, [](unsigned m) { return sum_aos(g_o, m); });
        run_case("tier1", "normalize", n, reps, 16.0,
                 [](unsigned m) { t1_normalize(m); }, [](unsigned m) { return sum_aos(g_o, m); });
        run_case("tier1", "collision", n, reps, 32.0,
                 [](unsigned m) { t1_collision(m); }, [](unsigned m) { return sum_aos(g_o, m); });
    }
    if (strcmp(tier, "tier2") == 0 || strcmp(tier, "all") == 0) {
        run_case("tier2", "add", n, reps, 24.0,
                 [](unsigned m) { t2_add(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier2", "scale", n, reps, 16.0,
                 [s](unsigned m) { t2_scale(m, s); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier2", "dot", n, reps, 16.0,
                 [](unsigned m) { t2_dot(m); }, [](unsigned) { return g_dot_scalar; });
        run_case("tier2", "dot_elem", n, reps, 20.0,
                 [](unsigned m) { t2_dot_elem(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier2", "len", n, reps, 20.0,
                 [](unsigned m) { t2_len(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier2", "normalize", n, reps, 16.0,
                 [](unsigned m) { t2_normalize(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier2", "collision", n, reps, 32.0,
                 [](unsigned m) { t2_collision(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
    }
#if defined(FV_HAVE_AVX2)
    if (strcmp(tier, "tier3") == 0 || strcmp(tier, "all") == 0) {
        run_case("tier3", "add", n, reps, 24.0,
                 [](unsigned m) { t3_add(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3", "scale", n, reps, 16.0,
                 [s](unsigned m) { t3_scale(m, s); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3", "dot", n, reps, 16.0,
                 [](unsigned m) { t3_dot(m); }, [](unsigned) { return g_dot_scalar; });
        run_case("tier3hadd", "dot", n, reps, 16.0,
                 [](unsigned m) { t3_dot_hadd(m); }, [](unsigned) { return g_dot_scalar; });
        run_case("tier3_4acc", "dot", n, reps, 16.0,
                 [](unsigned m) { t3_dot4(m); }, [](unsigned) { return g_dot_scalar; });
        run_case("tier3", "dot_elem", n, reps, 20.0,
                 [](unsigned m) { t3_dot_elem(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3", "len", n, reps, 20.0,
                 [](unsigned m) { t3_len(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3", "normalize", n, reps, 16.0,
                 [](unsigned m) { t3_normalize(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3_rsqrt", "normalize", n, reps, 16.0,
                 [](unsigned m) { t3_normalize_rsqrt(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3", "collision", n, reps, 32.0,
                 [](unsigned m) { t3_collision(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
        run_case("tier3_rsqrt", "collision", n, reps, 32.0,
                 [](unsigned m) { t3_collision_rsqrt(m); }, [](unsigned m) { return sum_f(g_ox, m) + sum_f(g_oy, m); });
    }
#endif
}

int main(int argc, char **argv)
{
    const char *tier = "all";
    int cpu = 2;
    if (argc > 1) { tier = argv[1]; }
    if (argc > 2) { cpu = atoi(argv[2]); }

    init_data();
    pin_cpu(cpu);

    const long khz0 = read_cpu_khz(cpu);
    const double tsc_hz = calibrate_tsc_hz();

    printf("# flashvec bench: cfg=%s\n", FV_CFG_NAME);
    printf("# build: %s %s\n", __VERSION__, __DATE__);
    printf("# n_full=%u n_mid=%u n_small=%u reps=%u reps_mid=%u reps_small=%u trials=%d checks=%u\n",
           (unsigned)FV_N, (unsigned)FV_NMID, (unsigned)FV_NSMALL, (unsigned)FV_REPS,
           (unsigned)FV_REPS_MID, (unsigned)FV_REPS_SMALL, (int)FV_TRIALS, (unsigned)FV_CHECKS);
#if defined(FV_HAVE_AVX2)
    printf("# avx2=yes fma=%d\n", (int)(__FMA__ ? 1 : 0));
#else
    printf("# avx2=no\n");
#endif
    printf("# compiled_features: avx2=%d fma=%d fast_math=%d\n",
#if defined(__AVX2__)
           1,
#else
           0,
#endif
#if defined(__FMA__)
           1,
#else
           0,
#endif
#if defined(__FAST_MATH__)
           1
#else
           0
#endif
    );
    printf("# tsc_hz=%.0f (%.3f GHz) cpu_khz_before=%ld\n", tsc_hz, tsc_hz / 1e9, khz0);

    const unsigned check_n = (FV_CHECKS < FV_N) ? (unsigned)FV_CHECKS : (unsigned)FV_N;
    run_checks(check_n);

    print_load("before");
    bench_all(tier, (unsigned)FV_NSMALL);
    bench_all(tier, (unsigned)FV_NMID);
    bench_all(tier, (unsigned)FV_N);
    print_load("after");

    printf("# cpu_khz_after=%ld\n", read_cpu_khz(cpu));
    return 0;
}
