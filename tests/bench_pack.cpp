/* AoS <-> SoA pack/unpack cost, and what it does to the SoA entry points (P7).
 *
 * The acceptance protocol used to list this as unmeasured. It measures, at the same two tiers:
 *   pack_unpack : fv_vec2 AoS -> (x[], y[]) -> fv_vec2 AoS, no arithmetic;
 *   add         : fv_vec2_add_x4 directly on AoS vs pack + fv_vec2_add_soa + unpack;
 *   normalize   : per-element fv_vec2_normalize vs pack + fv_vec2_normalize_soa + unpack.
 *
 * The point is the decision rule, not a single number: conversion is a fixed per-element cost, so it
 * can eat the gain of a cheap op (add) but is amortized by an expensive one (normalize).
 *
 * Same-process interleaved A/B, 11 paired rounds, median of the paired ratios, A/A control.
 * Usage: bench_pack        (make bench-pack; taskset -c 0 recommended)
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
static fv_vec2 *g_a, *g_b, *g_o;
static float *g_ax, *g_ay, *g_bx, *g_by, *g_ox, *g_oy;

enum { ROUNDS = 11 };

static void pack(const fv_vec2 *a, float *x, float *y, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        x[i] = a[i].x;
        y[i] = a[i].y;
    }
}

static void unpack(const float *x, const float *y, fv_vec2 *o, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        o[i].x = x[i];
        o[i].y = y[i];
    }
}

static void pack_unpack(size_t n)
{
    pack(g_a, g_ax, g_ay, n);
    unpack(g_ax, g_ay, g_o, n);
}

static void add_aos(size_t n)
{
    fv_vec2_add_x4(g_a, g_b, g_o, n);
}

static void add_soa_total(size_t n)
{
    pack(g_a, g_ax, g_ay, n);
    pack(g_b, g_bx, g_by, n);
    fv_vec2_add_soa(g_ax, g_ay, g_bx, g_by, g_ox, g_oy, n);
    unpack(g_ox, g_oy, g_o, n);
}

static void normalize_scalar(size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        g_o[i] = fv_vec2_normalize(g_a[i]);
    }
}

static void normalize_soa_total(size_t n)
{
    pack(g_a, g_ax, g_ay, n);
    fv_vec2_normalize_soa(g_ax, g_ay, g_ox, g_oy, n);
    unpack(g_ox, g_oy, g_o, n);
}

using Fn = void (*)(size_t);

static void touch(void)
{
    size_t i;
    size_t step = g_n / 64u + 1u;
    double s = 0.0;
    for (i = 0; i < g_n; i += step) {
        s += (double)g_o[i].x + (double)g_o[i].y;
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
        double na = median(ta, ROUNDS) / (double)g_n;
        double nb = median(tb, ROUNDS) / (double)g_n;
        double med = median(ratio, ROUNDS);
        double lo = ratio[0];
        double hi = ratio[0];
        for (r = 1; r < ROUNDS; ++r) {
            if (ratio[r] < lo) lo = ratio[r];
            if (ratio[r] > hi) hi = ratio[r];
        }
        printf("  %-16s %-22s %8.3f ns/element\n", label, aname, na);
        printf("  %-16s %-22s %8.3f ns/element\n", "", bname, nb);
        printf("  %-16s ratio %-18s %8.3fx (min %.3f / max %.3f)\n", "", bname, med, lo, hi);
    }
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u };
    size_t si;

    printf("bench_pack: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios\n", (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8u;
        std::vector<fv_vec2> a(cap), b(cap), o(cap);
        std::vector<float> ax(cap), ay(cap), bx(cap), by(cap), ox(cap), oy(cap);
        size_t i;

        for (i = 0; i < cap; ++i) {
            a[i].x = (float)((i * 3) % 23) / 8.0f - 1.0f;
            a[i].y = (float)((i * 5) % 19) / 8.0f - 1.0f;
            b[i].x = (float)((i * 7) % 17) / 8.0f - 1.0f;
            b[i].y = (float)((i * 11) % 29) / 8.0f - 1.0f;
        }
        g_n = n;
        g_a = a.data(); g_b = b.data(); g_o = o.data();
        g_ax = ax.data(); g_ay = ay.data(); g_bx = bx.data(); g_by = by.data();
        g_ox = ox.data(); g_oy = oy.data();

        printf("\n== n = %zu ==\n", n);
        run_pair("pack A/A", pack_unpack, pack_unpack, "pack+unpack", "pack+unpack (control)");
        run_pair("add A/A", add_aos, add_aos, "AoS add", "AoS add (control)");
        run_pair("add A/B", add_aos, add_soa_total, "AoS add", "pack+SoA+unpack");
        run_pair("norm A/A", normalize_scalar, normalize_scalar, "scalar normalize", "scalar (control)");
        run_pair("norm A/B", normalize_scalar, normalize_soa_total, "scalar normalize", "pack+SoA+unpack");
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
