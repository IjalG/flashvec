/* Measurement harness for the _fast approximate normalize entries (P7).
 *
 * Same protocol as tests/bench_physics.cpp: same-process interleaved A/B, 11 paired rounds, median of
 * the paired ratios, and a same-round A/A control. A = fast, B = exact, so the printed "speedup" is
 * exact/fast (> 1 means fast wins). The exact entry is the reference implementation, not a scalar loop:
 * both sides are batch SoA, so the measured difference is the rsqrtps + Newton path vs sqrtps + divps.
 *
 * Usage: bench_fast        (make bench-fast; taskset -c 0 recommended)
 */
#include "flashvec.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static double now_ns(void)
{
    using clock = std::chrono::steady_clock;
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now().time_since_epoch()).count();
}

static size_t g_n;
static volatile double g_sink;
static float *g_ax, *g_ay, *g_az, *g_bx, *g_by;
static float *g_ex, *g_ey, *g_ez;
static float *g_fx, *g_fy, *g_fz;

enum { ROUNDS = 11 };

static void vec2_exact(size_t n) { fv_vec2_normalize_soa(g_ax, g_ay, g_ex, g_ey, n); }
static void vec2_fast(size_t n) { fv_vec2_normalize_soa_fast(g_ax, g_ay, g_fx, g_fy, n); }
static void vec3_exact(size_t n) { fv_vec3_normalize_soa(g_ax, g_ay, g_az, g_ex, g_ey, g_ez, n); }
static void vec3_fast(size_t n) { fv_vec3_normalize_soa_fast(g_ax, g_ay, g_az, g_fx, g_fy, g_fz, n); }
static void madd_exact(size_t n) { fv_vec2_madd_normalize_soa(g_ax, g_ay, g_bx, g_by, 1.5f, g_ex, g_ey, n); }
static void madd_fast(size_t n) { fv_vec2_madd_normalize_soa_fast(g_ax, g_ay, g_bx, g_by, 1.5f, g_fx, g_fy, n); }

using Fn = void (*)(size_t);

static void touch(void)
{
    size_t i;
    size_t step = g_n / 64u + 1u;
    double s = 0.0;
    for (i = 0; i < g_n; i += step) {
        s += (double)g_ex[i] + (double)g_ey[i] + (double)g_ez[i];
        s += (double)g_fx[i] + (double)g_fy[i] + (double)g_fz[i];
    }
    g_sink += s;
}

static double measure(Fn fn)
{
    double t0 = now_ns();
    fn(g_n);
    return now_ns() - t0;
}

static int cmp_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa;
    double b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double median(double *v, int n)
{
    double tmp[ROUNDS];
    memcpy(tmp, v, (size_t)n * sizeof tmp[0]);
    qsort(tmp, (size_t)n, sizeof tmp[0], cmp_double);
    return tmp[n / 2];
}

static void run_pair(const char *label, Fn fast, Fn exact)
{
    double tf[ROUNDS];
    double te[ROUNDS];
    double ratio[ROUNDS];
    int r;
    for (r = 0; r < ROUNDS; ++r) {
        touch();
        tf[r] = measure(fast);
        touch();
        te[r] = measure(exact);
        ratio[r] = tf[r] > 0.0 ? te[r] / tf[r] : 0.0;
    }
    {
        double f = median(tf, ROUNDS) / (double)g_n;
        double e = median(te, ROUNDS) / (double)g_n;
        double med = median(ratio, ROUNDS);
        double lo = ratio[0];
        double hi = ratio[0];
        for (r = 1; r < ROUNDS; ++r) {
            if (ratio[r] < lo) lo = ratio[r];
            if (ratio[r] > hi) hi = ratio[r];
        }
        printf("  %-18s fast  %8.3f ns/element\n", label, f);
        printf("  %-18s exact %8.3f ns/element\n", "", e);
        printf("  %-18s speedup %6.3fx (min %.3f / max %.3f over %d paired rounds)\n",
               "", med, lo, hi, (int)ROUNDS);
    }
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_fast: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios, A/A control in the same round\n",
           (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<float> ax(cap), ay(cap), az(cap), bx(cap), by(cap);
        std::vector<float> ex(cap), ey(cap), ez(cap);
        std::vector<float> fx(cap), fy(cap), fz(cap);
        size_t i;

        for (i = 0; i < cap; ++i) {
            bx[i] = (float)((i * 11) % 17) / 8.0f - 1.0f;
            by[i] = (float)((i * 13) % 19) / 8.0f - 1.0f;
            ax[i] = (float)((i * 3) % 23) / 8.0f - 1.0f;
            ay[i] = (float)((i * 5) % 19) / 8.0f - 1.0f;
            az[i] = (float)((i * 7) % 17) / 8.0f - 1.0f;
        }
        g_n = n;
        g_ax = ax.data(); g_ay = ay.data(); g_az = az.data();
        g_bx = bx.data(); g_by = by.data();
        g_ex = ex.data(); g_ey = ey.data(); g_ez = ez.data();
        g_fx = fx.data(); g_fy = fy.data(); g_fz = fz.data();

        printf("\n== n = %zu (one component array: %.1f KiB) ==\n",
               n, (double)cap * sizeof(float) / 1024.0);
        run_pair("vec2 A/A", vec2_fast, vec2_fast);
        run_pair("vec2 A/B", vec2_fast, vec2_exact);
        run_pair("vec3 A/A", vec3_fast, vec3_fast);
        run_pair("vec3 A/B", vec3_fast, vec3_exact);
        run_pair("madd A/A", madd_fast, madd_fast);
        run_pair("madd A/B", madd_fast, madd_exact);
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
