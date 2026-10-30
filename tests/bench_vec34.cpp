/* Absolute speedup of the 3D / 4D batch entries over per-element library calls (P7).
 *
 * The README used to say this was not measured. This harness measures a representative set at the same
 * two tiers as the other benches (L2-resident n=8192 and a ~L3/DRAM n=1M):
 *   3D SoA : fv_vec3_add_soa / scale_soa / normalize_soa
 *   4D     : fv_vec4_add_x2 / scale_x2 (AoS) and fv_vec4_normalize_soa (SoA)
 *
 * Baseline: one library call per element (cross-TU, not inlinable), the same baseline the project uses
 * elsewhere. Protocol: same-process interleaved A/B, 11 paired rounds, median of the paired ratios, and
 * a same-round A/A control. A = batch, B = scalar, so speedup = scalar/batch.
 *
 * Usage: bench_vec34       (make bench-vec34; taskset -c 0 recommended)
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

static float *g_ax, *g_ay, *g_az;
static float *g_bx, *g_by, *g_bz;
static float *g_ox, *g_oy, *g_oz;
static fv_vec4 *g_va, *g_vb, *g_vo;
static float *g_wx, *g_wy, *g_wz, *g_ww;
static float *g_ux, *g_uy, *g_uz, *g_uw;

enum { ROUNDS = 11 };

/* ---------------- 3D batch / scalar ---------------- */

static void v3_add_batch(size_t n) { fv_vec3_add_soa(g_ax, g_ay, g_az, g_bx, g_by, g_bz, g_ox, g_oy, g_oz, n); }
static void v3_scale_batch(size_t n) { fv_vec3_scale_soa(g_ax, g_ay, g_az, 1.5f, g_ox, g_oy, g_oz, n); }
static void v3_norm_batch(size_t n) { fv_vec3_normalize_soa(g_ax, g_ay, g_az, g_ox, g_oy, g_oz, n); }

static void v3_add_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_vec3 a = fv_vec3_make(g_ax[i], g_ay[i], g_az[i]);
        fv_vec3 b = fv_vec3_make(g_bx[i], g_by[i], g_bz[i]);
        fv_vec3 r = fv_vec3_add(a, b);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
    }
}

static void v3_scale_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_vec3 a = fv_vec3_make(g_ax[i], g_ay[i], g_az[i]);
        fv_vec3 r = fv_vec3_scale(a, 1.5f);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
    }
}

static void v3_norm_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_vec3 a = fv_vec3_make(g_ax[i], g_ay[i], g_az[i]);
        fv_vec3 r = fv_vec3_normalize(a);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
    }
}

/* ---------------- 4D batch / scalar ---------------- */

static void v4_add_batch(size_t n) { fv_vec4_add_x2(g_va, g_vb, g_vo, n); }
static void v4_scale_batch(size_t n) { fv_vec4_scale_x2(g_va, 1.5f, g_vo, n); }
static void v4_norm_batch(size_t n) { fv_vec4_normalize_soa(g_wx, g_wy, g_wz, g_ww, g_ux, g_uy, g_uz, g_uw, n); }

static void v4_add_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_vo[i] = fv_vec4_add(g_va[i], g_vb[i]);
    }
}

static void v4_scale_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_vo[i] = fv_vec4_scale(g_va[i], 1.5f);
    }
}

static void v4_norm_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_vec4 a = fv_vec4_make(g_wx[i], g_wy[i], g_wz[i], g_ww[i]);
        fv_vec4 r = fv_vec4_normalize(a);
        g_ux[i] = r.x;
        g_uy[i] = r.y;
        g_uz[i] = r.z;
        g_uw[i] = r.w;
    }
}

/* ---------------- harness ---------------- */

using Fn = void (*)(size_t);

static void touch(void)
{
    size_t i;
    size_t step = g_n / 64u + 1u;
    double s = 0.0;
    for (i = 0; i < g_n; i += step) {
        s += (double)g_ox[i] + (double)g_oy[i] + (double)g_oz[i];
        s += (double)g_vo[i].x + (double)g_vo[i].y;
        s += (double)g_ux[i] + (double)g_uy[i] + (double)g_uz[i] + (double)g_uw[i];
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

static void run_pair(const char *label, Fn batch, Fn scalar)
{
    double tb[ROUNDS];
    double ts[ROUNDS];
    double ratio[ROUNDS];
    int r;
    for (r = 0; r < ROUNDS; ++r) {
        touch();
        tb[r] = measure(batch);
        touch();
        ts[r] = measure(scalar);
        ratio[r] = tb[r] > 0.0 ? ts[r] / tb[r] : 0.0;
    }
    {
        double b = median(tb, ROUNDS) / (double)g_n;
        double s = median(ts, ROUNDS) / (double)g_n;
        double med = median(ratio, ROUNDS);
        double lo = ratio[0];
        double hi = ratio[0];
        for (r = 1; r < ROUNDS; ++r) {
            if (ratio[r] < lo) lo = ratio[r];
            if (ratio[r] > hi) hi = ratio[r];
        }
        printf("  %-16s batch  %8.3f ns/element\n", label, b);
        printf("  %-16s scalar %8.3f ns/element\n", "", s);
        printf("  %-16s speedup %6.3fx (min %.3f / max %.3f over %d paired rounds)\n",
               "", med, lo, hi, (int)ROUNDS);
    }
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_vec34: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios, A/A control in the same round\n",
           (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<float> ax(cap), ay(cap), az(cap), bx(cap), by(cap), bz(cap);
        std::vector<float> ox(cap), oy(cap), oz(cap);
        std::vector<float> wx(cap), wy(cap), wz(cap), ww(cap);
        std::vector<float> ux(cap), uy(cap), uz(cap), uw(cap);
        std::vector<fv_vec4> va(cap), vb(cap), vo(cap);
        size_t i;

        for (i = 0; i < cap; ++i) {
            ax[i] = (float)((i * 3) % 23) / 8.0f - 1.0f;
            ay[i] = (float)((i * 5) % 19) / 8.0f - 1.0f;
            az[i] = (float)((i * 7) % 17) / 8.0f - 1.0f;
            bx[i] = (float)((i * 11) % 29) / 8.0f - 1.0f;
            by[i] = (float)((i * 13) % 31) / 8.0f - 1.0f;
            bz[i] = (float)((i * 17) % 37) / 8.0f - 1.0f;
            wx[i] = (float)((i * 3) % 23) / 8.0f - 1.0f;
            wy[i] = (float)((i * 5) % 19) / 8.0f - 1.0f;
            wz[i] = (float)((i * 7) % 17) / 8.0f - 1.0f;
            ww[i] = (float)((i * 11) % 29) / 8.0f - 1.0f;
            va[i] = fv_vec4_make((float)((i * 3) % 23) / 8.0f,
                                 (float)((i * 5) % 19) / 8.0f,
                                 (float)((i * 7) % 17) / 8.0f,
                                 (float)((i * 11) % 29) / 8.0f);
            vb[i] = fv_vec4_make((float)((i * 13) % 31) / 8.0f,
                                 (float)((i * 17) % 37) / 8.0f,
                                 (float)((i * 19) % 41) / 8.0f,
                                 (float)((i * 23) % 43) / 8.0f);
        }

        g_n = n;
        g_ax = ax.data(); g_ay = ay.data(); g_az = az.data();
        g_bx = bx.data(); g_by = by.data(); g_bz = bz.data();
        g_ox = ox.data(); g_oy = oy.data(); g_oz = oz.data();
        g_wx = wx.data(); g_wy = wy.data(); g_wz = wz.data(); g_ww = ww.data();
        g_ux = ux.data(); g_uy = uy.data(); g_uz = uz.data(); g_uw = uw.data();
        g_va = va.data(); g_vb = vb.data(); g_vo = vo.data();

        printf("\n== n = %zu ==\n", n);
        run_pair("vec3_add A/A", v3_add_batch, v3_add_batch);
        run_pair("vec3_add A/B", v3_add_batch, v3_add_scalar);
        run_pair("vec3_scale A/A", v3_scale_batch, v3_scale_batch);
        run_pair("vec3_scale A/B", v3_scale_batch, v3_scale_scalar);
        run_pair("vec3_norm A/A", v3_norm_batch, v3_norm_batch);
        run_pair("vec3_norm A/B", v3_norm_batch, v3_norm_scalar);
        run_pair("vec4_add A/A", v4_add_batch, v4_add_batch);
        run_pair("vec4_add A/B", v4_add_batch, v4_add_scalar);
        run_pair("vec4_scale A/A", v4_scale_batch, v4_scale_batch);
        run_pair("vec4_scale A/B", v4_scale_batch, v4_scale_scalar);
        run_pair("vec4_norm A/A", v4_norm_batch, v4_norm_batch);
        run_pair("vec4_norm A/B", v4_norm_batch, v4_norm_scalar);
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
