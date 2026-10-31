/* fv_mat3_mul_soa / fv_mat3_inverse_soa: AVX2 kernels vs per-element library calls (P1).
 *
 * Same protocol as the other benches: same-process interleaved A/B, 11 paired rounds, median of the
 * paired ratios, same-round A/A control. A = batch, B = scalar, speedup = scalar/batch.
 *
 * Usage: bench_mat3        (make bench-mat3; taskset -c 0 recommended)
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
static float *g_a[9], *g_b[9], *g_o[9];

enum { ROUNDS = 11 };

static void mul_batch(size_t n)
{
    fv_mat3_mul_soa(g_a[0], g_a[1], g_a[2], g_a[3], g_a[4], g_a[5], g_a[6], g_a[7], g_a[8],
                    g_b[0], g_b[1], g_b[2], g_b[3], g_b[4], g_b[5], g_b[6], g_b[7], g_b[8],
                    g_o[0], g_o[1], g_o[2], g_o[3], g_o[4], g_o[5], g_o[6], g_o[7], g_o[8], n);
}

static void inv_batch(size_t n)
{
    fv_mat3_inverse_soa(g_a[0], g_a[1], g_a[2], g_a[3], g_a[4], g_a[5], g_a[6], g_a[7], g_a[8],
                        g_o[0], g_o[1], g_o[2], g_o[3], g_o[4], g_o[5], g_o[6], g_o[7], g_o[8], n);
}

static void mul_scalar(size_t n)
{
    size_t i;
    int k;
    for (i = 0; i < n; ++i) {
        fv_mat3 a;
        fv_mat3 b;
        fv_mat3 r;
        for (k = 0; k < 9; ++k) {
            a.m[k] = g_a[k][i];
            b.m[k] = g_b[k][i];
        }
        r = fv_mat3_mul(&a, &b);
        for (k = 0; k < 9; ++k) {
            g_o[k][i] = r.m[k];
        }
    }
}

static void inv_scalar(size_t n)
{
    size_t i;
    int k;
    for (i = 0; i < n; ++i) {
        fv_mat3 a;
        fv_mat3 r;
        for (k = 0; k < 9; ++k) {
            a.m[k] = g_a[k][i];
        }
        r = fv_mat3_inverse(&a);
        for (k = 0; k < 9; ++k) {
            g_o[k][i] = r.m[k];
        }
    }
}

using Fn = void (*)(size_t);

static void touch(void)
{
    size_t i;
    size_t step = g_n / 64u + 1u;
    double s = 0.0;
    for (i = 0; i < g_n; i += step) {
        s += (double)g_o[0][i] + (double)g_o[4][i] + (double)g_o[8][i];
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
        printf("  %-14s batch  %8.3f ns/element\n", label, b);
        printf("  %-14s scalar %8.3f ns/element\n", "", s);
        printf("  %-14s speedup %6.3fx (min %.3f / max %.3f)\n", "", med, lo, hi);
    }
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_mat3: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios\n", (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<float> a[9], b[9], o[9];
        size_t idx;
        int k;
        for (k = 0; k < 9; ++k) {
            a[k].resize(cap);
            b[k].resize(cap);
            o[k].resize(cap);
            g_a[k] = a[k].data();
            g_b[k] = b[k].data();
            g_o[k] = o[k].data();
        }
        for (idx = 0; idx < cap; ++idx) {
            for (k = 0; k < 9; ++k) {
                a[k][idx] = (float)((idx * (k + 3)) % 23) / 8.0f - 1.0f;
                b[k][idx] = (float)((idx * (k + 5)) % 19) / 8.0f - 1.0f;
            }
            /* Make the diagonal dominant so inverse lanes are mostly non-singular. */
            a[0][idx] += 1.0f;
            a[4][idx] += 1.0f;
            a[8][idx] += 1.0f;
        }
        g_n = n;
        printf("\n== n = %zu ==\n", n);
        run_pair("mul A/A", mul_batch, mul_batch);
        run_pair("mul A/B", mul_batch, mul_scalar);
        run_pair("inv A/A", inv_batch, inv_batch);
        run_pair("inv A/B", inv_batch, inv_scalar);
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
