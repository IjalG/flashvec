/* Scalar API: header-inline (default) vs in-library (FV_USE_EXTERN) — the P0 measurement.
 *
 * This file is compiled twice, once per form, and linked into one process. The workload is the
 * report's P-group hot set: vec2 add/lerp/madd/dot/normalize, vec3 add/dot, and mat3 transform.
 * Same-process interleaved A/B, 11 paired rounds, median ratio.
 *
 * Usage: bench_scalar        (make bench-scalar; taskset -c 0 recommended)
 */
#include "flashvec.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifdef BENCH_EXTERN
#define BENCH_FN bench_extern
#else
#define BENCH_FN bench_inline
#endif

enum { ROUNDS = 11 };

static volatile double g_sink;

double bench_inline(const fv_vec2 *a2, const fv_vec2 *b2, fv_vec2 *o2, size_t n,
                    const fv_vec3 *a3, const fv_vec3 *b3, const fv_mat3 *m3);
double bench_extern(const fv_vec2 *a2, const fv_vec2 *b2, fv_vec2 *o2, size_t n,
                    const fv_vec3 *a3, const fv_vec3 *b3, const fv_mat3 *m3);

double BENCH_FN(const fv_vec2 *a2, const fv_vec2 *b2, fv_vec2 *o2, size_t n,
                const fv_vec3 *a3, const fv_vec3 *b3, const fv_mat3 *m3)
{
    double s = 0.0;
    size_t i;
    for (i = 0; i < n; ++i) {
        fv_vec2 x = fv_vec2_add(a2[i], b2[i]);
        fv_vec2 y = fv_vec2_lerp(a2[i], b2[i], 0.25f);
        fv_vec2 z = fv_vec2_normalize(fv_vec2_madd(x, y, 0.5f));
        s += (double)fv_vec2_dot(z, a2[i]);
        s += (double)fv_vec2_dot(fv_vec2_normalize(b2[i]), a2[i]);
        o2[i] = z;
    }
    for (i = 0; i < n; ++i) {
        s += (double)fv_vec3_dot(fv_vec3_add(a3[i], b3[i]), a3[i]);
        s += (double)fv_mat3_transform_vec3(&m3[i % 64u], a3[i]).x;
    }
    g_sink += s;
    return s;
}

#ifndef BENCH_NO_MAIN

static double now_ns(void)
{
    using clock = std::chrono::steady_clock;
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now().time_since_epoch()).count();
}

static double measure(double (*fn)(const fv_vec2 *, const fv_vec2 *, fv_vec2 *, size_t,
                                   const fv_vec3 *, const fv_vec3 *, const fv_mat3 *),
                      const fv_vec2 *a2, const fv_vec2 *b2, fv_vec2 *o2, size_t n,
                      const fv_vec3 *a3, const fv_vec3 *b3, const fv_mat3 *m3)
{
    double t0 = now_ns();
    (void)fn(a2, b2, o2, n, a3, b3, m3);
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

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_scalar: default (-O3 -fno-math-errno -ffp-contract=off)\n");
    printf("inline = default header form; extern = FV_USE_EXTERN library calls\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median ratio\n", (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<fv_vec2> a2(cap), b2(cap), o2(cap);
        std::vector<fv_vec3> a3(cap), b3(cap);
        std::vector<fv_mat3> m3(64);
        size_t i;
        int k;
        for (i = 0; i < cap; ++i) {
            a2[i] = fv_vec2_make((float)((i * 3) % 23) / 8.0f - 1.0f,
                                 (float)((i * 5) % 19) / 8.0f - 1.0f);
            b2[i] = fv_vec2_make((float)((i * 7) % 17) / 8.0f - 1.0f,
                                 (float)((i * 11) % 29) / 8.0f - 1.0f);
            a3[i] = fv_vec3_make((float)((i * 3) % 13) / 8.0f - 0.5f,
                                 (float)((i * 5) % 17) / 8.0f - 0.5f,
                                 (float)((i * 7) % 19) / 8.0f - 0.5f);
            b3[i] = fv_vec3_make((float)((i * 11) % 23) / 8.0f - 0.5f,
                                 (float)((i * 13) % 29) / 8.0f - 0.5f,
                                 (float)((i * 17) % 31) / 8.0f - 0.5f);
        }
        for (k = 0; k < 64; ++k) {
            m3[k] = fv_mat3_from_rows(fv_vec3_make(1.0f + (float)k / 64.0f, 0.1f, 0.2f),
                                      fv_vec3_make(-0.1f, 1.0f - (float)k / 128.0f, 0.3f),
                                      fv_vec3_make(0.2f, -0.3f, 1.0f));
        }

        {
            double ti[ROUNDS], te[ROUNDS], ratio[ROUNDS];
            int r;
            for (r = 0; r < ROUNDS; ++r) {
                ti[r] = measure(bench_inline, a2.data(), b2.data(), o2.data(), n,
                                a3.data(), b3.data(), m3.data());
                te[r] = measure(bench_extern, a2.data(), b2.data(), o2.data(), n,
                                a3.data(), b3.data(), m3.data());
                ratio[r] = ti[r] > 0.0 ? te[r] / ti[r] : 0.0;
            }
            {
                double mi = median(ti, ROUNDS) / (double)n;
                double me = median(te, ROUNDS) / (double)n;
                double med = median(ratio, ROUNDS);
                double lo = ratio[0], hi = ratio[0];
                for (r = 1; r < ROUNDS; ++r) {
                    if (ratio[r] < lo) lo = ratio[r];
                    if (ratio[r] > hi) hi = ratio[r];
                }
                printf("\n== n = %zu ==\n", n);
                printf("  inline %8.3f ns/element\n", mi);
                printf("  extern %8.3f ns/element\n", me);
                printf("  extern/inline speedup %6.3fx (min %.3f / max %.3f)\n", med, lo, hi);
            }
        }
    }

    printf("\nchecksum sink = %.6f\n", (double)g_sink);
    return 0;
}

#endif /* BENCH_NO_MAIN */
