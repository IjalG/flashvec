/*
 * bench/membw.c -- 本机内存带宽实测（用来解释为什么 N=1<<20 的三档数字会撞在一起）
 *
 * 三个内核：triad(c=a+b)、纯写、4 路累加读。每种尺寸取 5 次最小值。
 * 报告里的 GB/s 按“逻辑字节”算（triad = 2 读 + 1 写 = 12 B/元素），另给一列把
 * 写分配的隐性读也算进去的 16 B/元素，方便和 memcpy/memset 的读数对照。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN (8u << 20)

static float *a, *b, *c;

static double now_s(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static void init(unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        a[i] = 1.0f + (float)i;
        b[i] = 0.5f * (float)i;
        c[i] = 0.25f * (float)i;
    }
}

static double best_triad(unsigned n, int reps)
{
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        const double t0 = now_s();
        for (unsigned i = 0; i < n; ++i) { c[i] = a[i] + b[i]; }
        const double t1 = now_s();
        if (t1 - t0 < best) { best = t1 - t0; }
    }
    return best;
}

static double best_write(unsigned n, int reps)
{
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        const double t0 = now_s();
        for (unsigned i = 0; i < n; ++i) { c[i] = 0.0f; }
        const double t1 = now_s();
        if (t1 - t0 < best) { best = t1 - t0; }
    }
    return best;
}

static double best_read(unsigned n, int reps)
{
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        const double t0 = now_s();
        float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f, s3 = 0.0f;
        for (unsigned i = 0; i + 4 <= n; i += 4) {
            s0 += a[i]; s1 += a[i + 1]; s2 += a[i + 2]; s3 += a[i + 3];
        }
        const double t1 = now_s();
        if (t1 - t0 < best) { best = t1 - t0; }
        if (s0 + s1 + s2 + s3 == 12345.0f) { printf("!"); }
    }
    return best;
}

static double best_memcpy(unsigned n, int reps)
{
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        const double t0 = now_s();
        memcpy(c, a, (size_t)n * sizeof(float));
        const double t1 = now_s();
        if (t1 - t0 < best) { best = t1 - t0; }
    }
    return best;
}

int main(void)
{
    a = (float *)aligned_alloc(64, (size_t)MAXN * sizeof(float));
    b = (float *)aligned_alloc(64, (size_t)MAXN * sizeof(float));
    c = (float *)aligned_alloc(64, (size_t)MAXN * sizeof(float));
    if (!a || !b || !c) { return 1; }

    printf("# membw: 每项 5 次取最小值，CPU 绑定由外部 taskset 决定\n");
    printf("%-12s %10s %12s %12s %12s %12s %12s\n", "size", "working_MB", "triad_GBs", "triad16_GBs",
           "write_GBs", "read_GBs", "memcpy_GBs");

    const unsigned sizes[3] = { 1u << 20, 1u << 22, 1u << 23 };
    const int reps[3] = { 50, 20, 10 };
    for (int k = 0; k < 3; ++k) {
        const unsigned n = sizes[k];
        init(n);
        const double wmb = (double)n * 4.0 * 3.0 / (1024.0 * 1024.0);
        const double t_triad = best_triad(n, reps[k]);
        const double t_write = best_write(n, reps[k]);
        const double t_read = best_read(n, reps[k]);
        const double t_copy = best_memcpy(n, reps[k]);
        const double bytes = (double)n * 4.0;
        printf("N=%-10u %10.1f %12.2f %12.2f %12.2f %12.2f %12.2f\n", n, wmb,
               3.0 * bytes / t_triad / 1e9,
               4.0 * bytes / t_triad / 1e9,
               bytes / t_write / 1e9,
               bytes / t_read / 1e9,
               2.0 * bytes / t_copy / 1e9);
    }
    free(a); free(b); free(c);
    return 0;
}
