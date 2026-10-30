/* Measurement harness for the P6 physics primitives (fv_quat / fv_mat3 SoA batch).
 *
 * Same protocol as tests/bench_addmul.cpp: same-process interleaved A/B, 11 paired rounds, median of
 * the paired ratios, and a same-round A/A control that is printed with the real numbers (a gain inside
 * the A/A band is not a gain). The baseline is the per-element library call -- a cross-TU call that
 * cannot be inlined -- which is exactly the consumer the batch entry point exists for.
 *
 * Cases:
 *   quat_mul        fv_quat_mul_soa          vs  fv_quat_mul per element
 *   quat_normalize  fv_quat_normalize_soa    vs  fv_quat_normalize per element
 *   quat_rotate     fv_quat_rotate_vec3_soa  vs  fv_quat_rotate_vec3 per element
 *   mat3_transform  fv_mat3_transform_vec3_soa vs fv_mat3_transform_vec3 per element
 *
 * Read a single point as "cleared the same-round A/A floor", not as a precise multiplier. The numbers
 * are preliminary evidence for the P6 implementation decision (plain C inside target("avx2")); they are
 * not part of the five-criterion acceptance protocol until the hardened A/B fixture is run.
 *
 * Usage: bench_physics        (make bench-physics; taskset -c 0 recommended)
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

static float *g_ax, *g_ay, *g_az, *g_aw;
static float *g_bx, *g_by, *g_bz, *g_bw;
static float *g_vx, *g_vy, *g_vz;
static float *g_ox, *g_oy, *g_oz, *g_ow;
static float *g_m00, *g_m01, *g_m02;
static float *g_m10, *g_m11, *g_m12;
static float *g_m20, *g_m21, *g_m22;

enum { ROUNDS = 11 };

/* ---------------- batch ---------------- */

static void quat_mul_batch(size_t n)
{
    fv_quat_mul_soa(g_ax, g_ay, g_az, g_aw, g_bx, g_by, g_bz, g_bw,
                    g_ox, g_oy, g_oz, g_ow, n);
}

static void quat_normalize_batch(size_t n)
{
    fv_quat_normalize_soa(g_ax, g_ay, g_az, g_aw, g_ox, g_oy, g_oz, g_ow, n);
}

static void quat_rotate_batch(size_t n)
{
    fv_quat_rotate_vec3_soa(g_ax, g_ay, g_az, g_aw, g_vx, g_vy, g_vz,
                            g_ox, g_oy, g_oz, n);
}

static void mat3_transform_batch(size_t n)
{
    fv_mat3_transform_vec3_soa(g_m00, g_m01, g_m02, g_m10, g_m11, g_m12,
                               g_m20, g_m21, g_m22, g_vx, g_vy, g_vz,
                               g_ox, g_oy, g_oz, n);
}

/* ---------------- per-element library calls ---------------- */

static void quat_mul_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_quat b = fv_quat_make(g_bx[i], g_by[i], g_bz[i], g_bw[i]);
        fv_quat r = fv_quat_mul(a, b);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
        g_ow[i] = r.w;
    }
}

static void quat_normalize_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_quat a = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_quat r = fv_quat_normalize(a);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
        g_ow[i] = r.w;
    }
}

static void quat_rotate_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_quat q = fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]);
        fv_vec3 v = fv_vec3_make(g_vx[i], g_vy[i], g_vz[i]);
        fv_vec3 r = fv_quat_rotate_vec3(q, v);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
    }
}

static void mat3_transform_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_mat3 m;
        fv_vec3 v = fv_vec3_make(g_vx[i], g_vy[i], g_vz[i]);
        fv_vec3 r;
        m.m[0] = g_m00[i]; m.m[1] = g_m01[i]; m.m[2] = g_m02[i];
        m.m[3] = g_m10[i]; m.m[4] = g_m11[i]; m.m[5] = g_m12[i];
        m.m[6] = g_m20[i]; m.m[7] = g_m21[i]; m.m[8] = g_m22[i];
        r = fv_mat3_transform_vec3(&m, v);
        g_ox[i] = r.x;
        g_oy[i] = r.y;
        g_oz[i] = r.z;
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
        s += (double)g_ox[i] + (double)g_oy[i] + (double)g_oz[i] + (double)g_ow[i];
    }
    g_sink += s;
}

static double measure(Fn fn)
{
    double t0 = now_ns();
    fn(g_n);
    double t1 = now_ns();
    return t1 - t0;
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

/* A = candidate (batch), B = baseline (per-element calls). The printed "speedup" is B/A, so > 1 means
 * the batch entry wins. The A/A control passes the same function twice. */
static void run_pair(const char *label, Fn a, Fn b, const char *aname, const char *bname)
{
    double ta[ROUNDS];
    double tb[ROUNDS];
    double ratio[ROUNDS];
    int r;
    for (r = 0; r < ROUNDS; ++r) {
        touch();
        ta[r] = measure(a);
        touch();
        tb[r] = measure(b);
        ratio[r] = ta[r] > 0.0 ? tb[r] / ta[r] : 0.0;
    }
    {
        double ns_a = median(ta, ROUNDS) / (double)g_n;
        double ns_b = median(tb, ROUNDS) / (double)g_n;
        double med = median(ratio, ROUNDS);
        double lo = ratio[0];
        double hi = ratio[0];
        for (r = 1; r < ROUNDS; ++r) {
            if (ratio[r] < lo) lo = ratio[r];
            if (ratio[r] > hi) hi = ratio[r];
        }
        printf("  %-22s %-26s %8.3f ns/element\n", label, aname, ns_a);
        printf("  %-22s %-26s %8.3f ns/element\n", "", bname, ns_b);
        printf("  %-22s speedup %-19s %8.3fx  (min %.3f / max %.3f over %d paired rounds)\n",
               "", bname, med, lo, hi, (int)ROUNDS);
    }
}

static void fill(size_t cap)
{
    size_t i;
    for (i = 0; i < cap; ++i) {
        g_ax[i] = (float)((i * 3) % 23) / 8.0f - 1.0f;
        g_ay[i] = (float)((i * 5) % 19) / 8.0f - 1.0f;
        g_az[i] = (float)((i * 7) % 17) / 8.0f - 1.0f;
        g_aw[i] = (float)((i * 11) % 13) / 8.0f - 1.0f;
        g_bx[i] = (float)((i * 13) % 29) / 8.0f - 1.0f;
        g_by[i] = (float)((i * 17) % 31) / 8.0f - 1.0f;
        g_bz[i] = (float)((i * 19) % 37) / 8.0f - 1.0f;
        g_bw[i] = (float)((i * 23) % 41) / 8.0f - 1.0f;
        g_vx[i] = (float)((i * 2) % 11) / 4.0f - 1.0f;
        g_vy[i] = (float)((i * 3) % 7) / 4.0f - 1.0f;
        g_vz[i] = (float)((i * 5) % 13) / 4.0f - 1.0f;
        g_m00[i] = (float)((i * 3) % 17) / 8.0f - 1.0f;
        g_m01[i] = (float)((i * 5) % 19) / 8.0f - 1.0f;
        g_m02[i] = (float)((i * 7) % 23) / 8.0f - 1.0f;
        g_m10[i] = (float)((i * 11) % 29) / 8.0f - 1.0f;
        g_m11[i] = (float)((i * 13) % 31) / 8.0f - 1.0f;
        g_m12[i] = (float)((i * 17) % 37) / 8.0f - 1.0f;
        g_m20[i] = (float)((i * 19) % 41) / 8.0f - 1.0f;
        g_m21[i] = (float)((i * 23) % 43) / 8.0f - 1.0f;
        g_m22[i] = (float)((i * 29) % 47) / 8.0f - 1.0f;
        g_ox[i] = 0.0f;
        g_oy[i] = 0.0f;
        g_oz[i] = 0.0f;
        g_ow[i] = 0.0f;
    }

    /* The rotate entry assumes unit quaternions: normalize the input once, outside the timed loop. */
    for (i = 0; i < cap; ++i) {
        fv_quat q = fv_quat_normalize(fv_quat_make(g_ax[i], g_ay[i], g_az[i], g_aw[i]));
        g_ax[i] = q.x;
        g_ay[i] = q.y;
        g_az[i] = q.z;
        g_aw[i] = q.w;
    }
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_physics: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios, A/A control in the same round\n",
           (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<float> ax(cap), ay(cap), az(cap), aw(cap);
        std::vector<float> bx(cap), by(cap), bz(cap), bw(cap);
        std::vector<float> vx(cap), vy(cap), vz(cap);
        std::vector<float> ox(cap), oy(cap), oz(cap), ow(cap);
        std::vector<float> m00(cap), m01(cap), m02(cap);
        std::vector<float> m10(cap), m11(cap), m12(cap);
        std::vector<float> m20(cap), m21(cap), m22(cap);

        g_n = n;
        g_ax = ax.data(); g_ay = ay.data(); g_az = az.data(); g_aw = aw.data();
        g_bx = bx.data(); g_by = by.data(); g_bz = bz.data(); g_bw = bw.data();
        g_vx = vx.data(); g_vy = vy.data(); g_vz = vz.data();
        g_ox = ox.data(); g_oy = oy.data(); g_oz = oz.data(); g_ow = ow.data();
        g_m00 = m00.data(); g_m01 = m01.data(); g_m02 = m02.data();
        g_m10 = m10.data(); g_m11 = m11.data(); g_m12 = m12.data();
        g_m20 = m20.data(); g_m21 = m21.data(); g_m22 = m22.data();

        fill(cap);

        printf("\n== n = %zu (one quat SoA component array: %.1f KiB) ==\n",
               n, (double)cap * sizeof(float) / 1024.0);
        run_pair("quat_mul A/A", quat_mul_batch, quat_mul_batch, "batch", "batch (control)");
        run_pair("quat_mul A/B", quat_mul_batch, quat_mul_scalar, "batch", "scalar calls");
        run_pair("quat_norm A/A", quat_normalize_batch, quat_normalize_batch, "batch", "batch (control)");
        run_pair("quat_norm A/B", quat_normalize_batch, quat_normalize_scalar, "batch", "scalar calls");
        run_pair("quat_rot A/A", quat_rotate_batch, quat_rotate_batch, "batch", "batch (control)");
        run_pair("quat_rot A/B", quat_rotate_batch, quat_rotate_scalar, "batch", "scalar calls");
        run_pair("mat3_tf A/A", mat3_transform_batch, mat3_transform_batch, "batch", "batch (control)");
        run_pair("mat3_tf A/B", mat3_transform_batch, mat3_transform_scalar, "batch", "scalar calls");
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
