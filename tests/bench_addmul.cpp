/* Measurement harness for the integrator primitives (P4 group), same protocol as the other measured entries in
 * include/flashvec.h: same-process interleaved A/B, 11 paired repetitions, median of the paired ratios, and a same-round
 * A/A control that is run and printed with the real numbers (a gain inside the A/A band is not a gain).
 *
 * Cases (default tier = the Makefile's -O3 -fno-math-errno -ffp-contract=off, no -march, i.e. the published baseline):
 *   addmul   : in-place fv_vec2_addmul_x4          vs  fv_vec2_scale_x4 into a temporary + fv_vec2_add_x4
 *   maddv    : fv_vec2_maddv_x4                    vs  a*b into a temporary (a scalar loop: the library has no
 *                                                     elementwise multiply) + fv_vec2_add_x4
 *   maddv_s  : fv_vec2_maddv_x4_stream             vs  the same two-pass workaround
 *
 * The two-pass workaround writes its final result with fv_vec2_add_x4(x, tmp, x) -- o aliasing a is outside that entry's
 * __restrict contract, but it is exactly what a consumer without a fused entry writes, and the loop GCC emits for it is
 * the same load-load-store. Both passes produce the byte traffic the header's accounting quotes (56 B/element).
 *
 * Usage: bench_addmul        (make bench-addmul)
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
static fv_vec2 *g_x;
static fv_vec2 *g_v;
static fv_vec2 *g_tmp;
static fv_vec2 *g_a;
static fv_vec2 *g_b;
static fv_vec2 *g_c;
static fv_vec2 *g_o;
static fv_vec2 *g_alias; /* == g_x, through a variable the compiler cannot fold: the aliased call below is deliberate */
static volatile double g_sink;

enum { ROUNDS = 11 };

/* ---------------- addmul ---------------- */
static void addmul_fused(size_t n)
{
    fv_vec2_addmul_x4(g_x, g_v, 0.001f, n);
}

static void addmul_twopass(size_t n)
{
    fv_vec2_scale_x4(g_v, 0.001f, g_tmp, n);
    /* o aliases a: outside _x4's __restrict contract, but this is the second pass a consumer without a fused entry
     * writes, and the loop GCC emits is the same load-load-store either way. */
    fv_vec2_add_x4(g_x, g_tmp, g_alias, n);
}

/* ---------------- maddv ---------------- */
static void maddv_fused(size_t n)
{
    fv_vec2_maddv_x4(g_a, g_b, g_c, g_o, n);
}

static void maddv_fused_stream(size_t n)
{
    fv_vec2_maddv_x4_stream(g_a, g_b, g_c, g_o, n);
}

static void maddv_twopass(size_t n)
{
    /* The library has no elementwise vector multiply, so a consumer writes this loop itself. With __restrict on the local
     * pointers -O3 vectorizes it (SSE at the default tier, since the default tier carries no -march); without them GCC has
     * to assume the three arrays may overlap and emits a run-time alias check. Both variants were measured; the restricted
     * one is the fair baseline, because "the compiler cannot vectorize the workaround" is not a property of the fusion. */
    const fv_vec2 *__restrict ap = g_a;
    const fv_vec2 *__restrict bp = g_b;
    fv_vec2 *__restrict tp = g_tmp;
    size_t i;
    for (i = 0; i < n; ++i) {
        tp[i].x = ap[i].x * bp[i].x;
        tp[i].y = ap[i].y * bp[i].y;
    }
    fv_vec2_add_x4(g_tmp, g_c, g_o, n);
}

using Fn = void (*)(size_t);

static void touch(void)
{
    size_t i, step = g_n / 64 + 1;
    double s = 0.0;
    for (i = 0; i < g_n; i += step) {
        s += (double)g_o[i].x + (double)g_x[i].y + (double)g_tmp[i].x;
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
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double median(double *v, int n)
{
    double tmp[ROUNDS];
    memcpy(tmp, v, (size_t)n * sizeof tmp[0]);
    qsort(tmp, (size_t)n, sizeof tmp[0], cmp_double);
    return tmp[n / 2];
}

/* Runs A and B interleaved for ROUNDS rounds and reports median ns/element plus the median paired ratio B/A.
 * The same function is used for the A/A control (A and B pointing at the same entry). */
static void run_pair(const char *label, Fn a, Fn b, const char *aname, const char *bname)
{
    double ta[ROUNDS], tb[ROUNDS], ratio[ROUNDS];
    int r;
    for (r = 0; r < ROUNDS; ++r) {
        touch();
        ta[r] = measure(a);
        touch();
        tb[r] = measure(b);
        ratio[r] = ta[r] > 0.0 ? tb[r] / ta[r] : 0.0;
    }
    double ns_a = median(ta, ROUNDS) / (double)g_n;
    double ns_b = median(tb, ROUNDS) / (double)g_n;
    double med = median(ratio, ROUNDS);
    double lo = ratio[0], hi = ratio[0];
    for (r = 1; r < ROUNDS; ++r) {
        if (ratio[r] < lo) lo = ratio[r];
        if (ratio[r] > hi) hi = ratio[r];
    }
    printf("  %-22s %-26s %8.3f ns/element\n", label, aname, ns_a);
    printf("  %-22s %-26s %8.3f ns/element\n", "", bname, ns_b);
    printf("  %-22s ratio %-20s %8.3fx  (min %.3f / max %.3f over %d paired rounds)\n", "", bname, med,
           lo, hi, (int)ROUNDS);
}

int main(void)
{
    static const size_t kSizes[] = { 8192u, 1048576u, 8388608u };
    size_t si;

    printf("bench_addmul: default tier (-O3 -fno-math-errno -ffp-contract=off, no -march)\n");
    printf("protocol: same-process interleaved A/B, %d paired rounds, median of the paired ratios, A/A control in the same round\n",
           (int)ROUNDS);

    for (si = 0; si < sizeof kSizes / sizeof kSizes[0]; ++si) {
        size_t n = kSizes[si];
        size_t cap = n + 8;
        g_n = n;
        printf("\n== n = %zu (fv_vec2 arrays: %.1f MiB each) ==\n", n, (double)cap * sizeof(fv_vec2) / (1024.0 * 1024.0));

        {
            std::vector<fv_vec2> x(cap), v(cap), tmp(cap);
            size_t i;
            for (i = 0; i < cap; ++i) {
                x[i].x = (float)(i % 97) * 1e-3f;
                x[i].y = (float)(i % 89) * 1e-3f;
                v[i].x = (float)(i % 13) - 6.0f;
                v[i].y = (float)(i % 17) - 8.0f;
                tmp[i].x = 0.0f;
                tmp[i].y = 0.0f;
            }
            g_x = x.data();
            g_v = v.data();
            g_alias = x.data();
            g_tmp = tmp.data();
            g_a = x.data();
            g_b = v.data();
            g_c = v.data();
            g_o = tmp.data();

            addmul_fused(n);
            printf("  -- addmul: in place vs fv_vec2_scale_x4 + fv_vec2_add_x4 (56 B/element vs 24 B/element) --\n");
            run_pair("addmul A/A", addmul_fused, addmul_fused, "in place", "in place (control)");
            run_pair("addmul A/B", addmul_fused, addmul_twopass, "in place", "scale->add");
        }

        {
            std::vector<fv_vec2> a(cap), b(cap), c(cap), o(cap), tmp(cap);
            size_t i;
            for (i = 0; i < cap; ++i) {
                a[i].x = (float)(i % 31) * 1e-2f;
                a[i].y = (float)(i % 29) * 1e-2f;
                b[i].x = (float)(i % 23) * 1e-2f;
                b[i].y = (float)(i % 19) * 1e-2f;
                c[i].x = (float)(i % 11);
                c[i].y = (float)(i % 7);
                o[i].x = 0.0f;
                o[i].y = 0.0f;
                tmp[i].x = 0.0f;
                tmp[i].y = 0.0f;
            }
            g_a = a.data();
            g_b = b.data();
            g_c = c.data();
            g_o = o.data();
            g_tmp = tmp.data();
            g_x = a.data();
            g_v = b.data();

            maddv_fused(n);
            printf("  -- maddv: fused vs a*b into a temporary + fv_vec2_add_x4 (64 B/element vs 40 B/element; each pass of\n"
                   "     the workaround pays its own RFO, which is why its baseline is 64 B and not 56 B) --\n");
            run_pair("maddv A/A", maddv_fused, maddv_fused, "maddv", "maddv (control)");
            run_pair("maddv A/B", maddv_fused, maddv_twopass, "maddv", "mul->add");
            printf("  -- maddv _stream: NT store vs the same two-pass workaround (64 B/element vs 32 B/element) --\n");
            run_pair("maddv_s A/A", maddv_fused_stream, maddv_fused_stream, "maddv_stream", "maddv_stream (control)");
            run_pair("maddv_s A/B", maddv_fused_stream, maddv_twopass, "maddv_stream", "mul->add");
        }
    }

    printf("\nchecksum sink = %.6f (keeps the timed loops alive)\n", (double)g_sink);
    return 0;
}
