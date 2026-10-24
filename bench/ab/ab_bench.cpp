/* bench/ab/ab_bench.cpp -- flashvec in-process A/B (interleaved) performance / throughput /
 * call-overhead fixture
 *
 * Why this file exists: the Performance section of the README (item 3) says this machine is
 * shared and its clock drifts, so the same code can differ by 20-25% between runs. "Run two
 * binaries and compare" is therefore unusable: clock, neighbour load and page mappings all
 * change between the two runs. What does hold is **one process, one data set, one time
 * window**, with the two implementations interleaved (A/B/A/B...), then a per-tier median
 * ratio.
 *
 * Three regimes
 *   1) large-N throughput sweep (default 2^12 ... 2^24, 8 tiers): ns/element and GB/s against
 *      memcpy/memset/pure-read/triad ceiling references, plus "percent of memcpy throughput at
 *      the same size"; both warm and cold are reported.
 *   2) small-N call-overhead regime (n = 4/8/16/32/64): ns/call and ns/element both, with
 *      scalar_loop as the control (element-wise scalar API = the pre-refactor / static-library
 *      consumer).
 *   3) the ceiling references span both regimes.
 *
 * Structure
 *   baseline  = git HEAD snapshot (bench/ab/snapshot/, produced by git show; the worktree
 *               src/ is never compiled)
 *   candidate = chosen at compile time: the HEAD snapshot itself (== baseline, fixture
 *               self-check) or any source tree
 *   each variant has its own TU (ab_variant.cpp, -include rename macros plus the chosen
 *   header/source), so the timing loop lives inside the variant TU and both forms -- library
 *   definition and header-inline -- are measured as they really are.
 *   Calls go through the single entry point ab_X_loop, so both sides are perfectly symmetric.
 *
 * Three semantic gates (a failing gate suppresses every ratio and every throughput number)
 *   1) the output is compared bit-for-bit against **this file's own independent scalar
 *      reference implementation** (including mask write-back, unreferenced idx positions and
 *      pipeline intermediates). The single exception is fv_vec2_dot_soa, where the header
 *      explicitly makes no bit-identical promise: it is compared with
 *      max(1e-5*|result|, 1e-6*sum|term|), a bound independent of n and of the corpus and
 *      steadier than a bare relative error.
 *   2) A and B must also agree bit-for-bit (same exception).
 *   3) the input-buffer checksum must not change across the timed region.
 *   --selfcheck adds one more: baseline vs baseline ratios must land in [0.85, 1.15];
 *   otherwise the fixture itself is untrustworthy and fails hard (exit code 3).
 *
 * Measurement conventions
 *   - single-threaded; the library has no threads and no global state. Core pinning is the
 *     Makefile's taskset -c (cpu0-7 share one CCX).
 *   - N counts fv_vec2 elements. One element occupies one float in ax and one in ay.
 *   - B/element = the bytes this function **really reads + writes** (checked by hand function
 *     by function, see the OPS table and the report).
 *     GB/s = (rd+wr) / (ns/element), i.e. bytes actually moved per second (read and written).
 *   - warm = the resident buffers are run repeatedly; cold = one call after
 *     madvise(DONTNEED) on the arrays under test (includes minor page faults; the kernel
 *     zeroes those pages, so cold only serves throughput comparison and stays out of the
 *     bit-level gates). The idx array is never dropped -- without it every scattered write
 *     lands on index 0 and random writes turn into hot writes.
 *   - cache hierarchy as measured here: L1d 32 KB, L2 512 KB, L3 **4 MB/CCX** (index3=4096K,
 *     shared_cpu_list=0-7; the "8 MB L3" line in the docs is wrong, a thread only ever sees
 *     its own 4 MB).
 *   - at least 5 trials per point, median taken, until the coefficient of variation < 5% or
 *     the cap is hit; the actual count goes into the TSV.
 *
 * Output: console tables + TSV (--tsv-out long table, --matrix-out matrix, --small-out small-N
 * table)
 */
#include "snapshot/include/flashvec.h"

#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#include "ab_abi.h"

/* The driver only includes the HEAD snapshot header, and fv_vec4 (the type added with
 * vec3/vec4 in P1) does not exist there yet. The reference implementations get a local
 * definition with the same layout; layout consistency is guaranteed by the variant TU's ABI
 * fingerprint and the 16 B static_assert in the header (and pinned once more here). Once HEAD
 * has fv_vec4, the Makefile auto-detection defines AB_SNAP_HAS_VEC4=1 and this block
 * disappears by itself. */
#ifndef AB_SNAP_HAS_VEC4
#define AB_SNAP_HAS_VEC4 0
#endif
#if !AB_SNAP_HAS_VEC4
struct fv_vec4 {
    fv_real x, y, z, w;
};
static_assert(sizeof(fv_vec4) == 16, "fv_vec4 must be 16 B (the header contract)");
#endif

#ifndef AB_LABEL
#define AB_LABEL "unlabeled"
#endif
#ifndef AB_CAND_DESC
#define AB_CAND_DESC "unknown"
#endif
#ifndef AB_GIT_REV
#define AB_GIT_REV "unknown"
#endif

#define AB_L1_BYTES ((size_t)32 * 1024)
#define AB_L2_BYTES ((size_t)512 * 1024)
#define AB_L3_BYTES ((size_t)4 * 1024 * 1024)
#define AB_BAND_LO 0.85
#define AB_BAND_HI 1.15
/* Lower bound on the A/A floor used by the **tier-level** health check. A tier floor is a
 * finite-sample estimate (15 ops x a few A/A pairs) and occasionally comes out absurdly small:
 * 0.3% was measured on n=16777216. The fixture's own run-to-run reproducibility on the large
 * tiers is 1-2%, so below that "stay within the floor" stops being a meaningful demand and
 * flags perfectly normal runs. The clamp applies to this tier-level check only -- R2/R3/R4 and
 * the all-point median criterion are untouched. See README section 8 (known limitations). */
#define AB_MIN_TIER_FLOOR 0.02
#define AB_ABS_FRAC 1e-6 /* dot_soa absolute bound: 1e-6 * sum|term| */
#define AB_REL_FRAC 1e-5 /* dot_soa relative bound: 1e-5 * |result| */
#define AB_T 0.375f
#define AB_S (-0.25f)

/* =====================================================================
 * 1. Environment probes
 * ===================================================================== */
/* The clock. Default is **this thread's CPU time**, not the wall clock: foreign load keeps
 * 30-50% of this machine busy, and wall time charges the code under test for every second it
 * was descheduled -- more than 20% at the DRAM tiers. CPU time counts only the time this
 * thread really spent on a core; it still includes stalls from neighbours stealing L3 (a real
 * effect) but not the time we were preempted. --clock=mono switches back to the wall clock so
 * the difference can be checked by hand. */
static clockid_t g_clock = CLOCK_THREAD_CPUTIME_ID;
static const char *g_clock_name = "cpu(CLOCK_THREAD_CPUTIME_ID)";

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(g_clock, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

static void cpu_mhz(double *mn, double *mx, double *av)
{
    FILE *f = std::fopen("/proc/cpuinfo", "r");
    *mn = *mx = *av = 0.0;
    if (!f) return;
    char line[256];
    double sum = 0.0, lo = 1e9, hi = 0.0;
    int cnt = 0;
    while (std::fgets(line, sizeof line, f)) {
        double v;
        if (std::sscanf(line, "cpu MHz : %lf", &v) == 1) {
            sum += v;
            lo = std::min(lo, v);
            hi = std::max(hi, v);
            ++cnt;
        }
    }
    std::fclose(f);
    if (cnt) {
        *mn = lo;
        *mx = hi;
        *av = sum / (double)cnt;
    }
}

static double loadavg1(void)
{
    FILE *f = std::fopen("/proc/loadavg", "r");
    if (!f) return -1.0;
    double a = -1.0, b = 0.0, c = 0.0;
    if (std::fscanf(f, "%lf %lf %lf", &a, &b, &c) != 3) a = -1.0;
    std::fclose(f);
    return a;
}

static long proc_kb(const char *key)
{
    FILE *f = std::fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256];
    long v = -1;
    const size_t kl = std::strlen(key);
    while (std::fgets(line, sizeof line, f)) {
        if (std::strncmp(line, key, kl) == 0) {
            long x = 0;
            if (std::sscanf(line + kl, ": %ld", &x) == 1) v = x;
            break;
        }
    }
    std::fclose(f);
    return v;
}

static long mem_available_mib(void)
{
    FILE *f = std::fopen("/proc/meminfo", "r");
    if (!f) return -1;
    char line[256];
    long kb = -1;
    while (std::fgets(line, sizeof line, f)) {
        if (std::strncmp(line, "MemAvailable:", 13) == 0) {
            long x = 0;
            if (std::sscanf(line + 12, ": %ld", &x) == 1) kb = x;
            break;
        }
    }
    std::fclose(f);
    return kb < 0 ? -1 : kb / 1024;
}

static long self_statm(size_t field)
{
    FILE *f = std::fopen("/proc/self/statm", "r");
    if (!f) return -1;
    long a = -1, b = -1, c = -1;
    if (std::fscanf(f, "%ld %ld %ld", &a, &b, &c) != 3) b = -1;
    std::fclose(f);
    return (field == 1) ? a : b;
}

static long self_minflt(void)
{
    FILE *f = std::fopen("/proc/self/stat", "r");
    if (!f) return -1;
    char buf[4096];
    const size_t got = std::fread(buf, 1, sizeof buf - 1, f);
    std::fclose(f);
    buf[got] = '\0';
    char *p = std::strrchr(buf, ')');
    if (!p) return -1;
    int k = 0;
    long minflt = -1;
    for (char *tok = std::strtok(p + 1, " "); tok; tok = std::strtok(nullptr, " ")) {
        if (++k == 8) {
            minflt = std::atol(tok);
            break;
        }
    }
    return minflt;
}

static void *xalloc(size_t bytes)
{
    void *p = nullptr;
    const size_t b = (bytes + 63u) & ~(size_t)63u;
    if (posix_memalign(&p, 64, b ? b : 64) != 0 || !p) {
        std::fprintf(stderr, "allocation failed: %zu bytes\n", bytes);
        std::exit(2);
    }
    return p;
}

/* =====================================================================
 * 2. Buffers and data
 * ===================================================================== */
/* Output form. KIND also fixes the **element width** (see elem_bytes): VEC2/SOA2/SOA4/FLOAT/
 * NONE elements are fv_vec2 (8 B), SOA3 is 12 B (fv_vec3) and VEC4/SOA4V is 16 B (fv_vec4).
 * Within one tier every op moves **the same number of bytes**: element count = sweep size
 * n (in 8 B units) * 8 / elem_bytes, so a vec4 vs vec2 ratio is a byte-normalized comparison
 * by construction. */
/* VEC4M is VEC4 plus a matrix operand: the driver then allocates a 16 B AoS input array and
 * the 16-float matrix, and keeps az/bz/aw/bw (SoA-only) out of the allocation. */
enum OutKind { OUT_VEC2, OUT_SOA2, OUT_SOA3, OUT_SOA4, OUT_SOA4V, OUT_VEC4, OUT_VEC4M,
               OUT_FLOAT, OUT_NONE };

/* Element width is no longer derived from KIND: a reduction op (dot/len) writes a float/NONE
 * yet its elements can still be 12/16 B. Width is part of the op's identity and lives in the
 * ELOP column of every registry row. */
enum Policy { POL_EXACT, POL_REL, POL_NOCHECK };
enum Side { SIDE_A = 0, SIDE_B = 1, SIDE_S = 2 };

struct OutSet {
    void *base;  /* 8n: shared by o / ox,oy / lenv (only one form is live at a time) */
    void *base2; /* 8n: collision (vox,voy) | vec3 oz | vec4 oz,ow */
    fv_vec2 *o;
    float *ox, *oy, *lenv;
    float *vx, *vy;    /* the two halves of base2, used by collision */
    float *oz, *ow;    /* the two halves of base2 for vec3/vec4 SoA (alias vx/vy) */
    float dot;
    double abs_sum; /* sum|term| from the dot_soa reference, feeds the absolute tolerance */
};

struct Buf {
    size_t n;
    size_t alloc_bytes;
    fv_vec2 *a, *b;
    uint8_t *mask;
    uint32_t *idx;
    float *ax, *ay, *bx, *by;
    float *px, *py, *qx, *qy;
    /* Component arrays for wide elements (vec3 12 B / vec4 16 B): allocated only when a wide
     * op is selected, so the default sweep (all vec2) costs not a single extra byte. */
    float *az, *bz, *aw, *bw;
    /* Matrix transforms: 16 B AoS input (one fv_vec4 = 4 floats) plus the 16-float matrix,
     * row-major with its last row (0,0,0,1). Allocated only when such an op is selected -- an
     * AoS element is 16 B, so the array holds nel() = n/2 elements = 2n floats = the same 8n
     * bytes as a/b. The output reuses OutSet::base (also 8n bytes, 16 B capable). */
    float *vin;
    float mat[16];
    OutSet out[2];
    double in_cs;
};

/* mask density (--mask-density): 0.5 = the old behaviour */
static double g_mask_density = 0.5;

static uint64_t g_rng = 0x243f6a8885a308d3ull;
static uint64_t rnd64(void)
{
    g_rng += 0x9e3779b97f4a7c15ull;
    uint64_t z = g_rng;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}
static float rnd_pm1(void)
{
    return (float)((double)(rnd64() >> 11) * (1.0 / 4503599627370496.0) - 1.0);
}

static void prefill(OutSet &O, size_t n);

/* True when a vec3/vec4 op is selected: allocate az/bz/aw/bw and base2 (vec4 SoA needs
 * oz/ow). The default sweep (all vec2) allocates nothing and keeps the old footprint. */
static bool g_wide = false;
/* True when a matrix-transform op is selected: allocate vin + mat (see Buf). Kept separate
 * from g_wide so a 16 B AoS op does not pay for the SoA component arrays (4 x 4n floats plus
 * 16n bytes of base2 at the largest tier). */
static bool g_mat4 = false;

static void alloc_buf(Buf &B, size_t n, bool with_pipe)
{
    std::memset(&B, 0, sizeof B);
    B.n = n;
    B.a = (fv_vec2 *)xalloc(n * sizeof(fv_vec2));
    B.b = (fv_vec2 *)xalloc(n * sizeof(fv_vec2));
    B.mask = (uint8_t *)xalloc(n);
    B.idx = (uint32_t *)xalloc(n * sizeof(uint32_t));
    B.ax = (float *)xalloc(n * sizeof(float));
    B.ay = (float *)xalloc(n * sizeof(float));
    B.bx = (float *)xalloc(n * sizeof(float));
    B.by = (float *)xalloc(n * sizeof(float));
    if (with_pipe) {
        B.px = (float *)xalloc(n * sizeof(float));
        B.py = (float *)xalloc(n * sizeof(float));
        B.qx = (float *)xalloc(n * sizeof(float));
        B.qy = (float *)xalloc(n * sizeof(float));
    }
    if (g_wide) {
        B.az = (float *)xalloc(n * sizeof(float));
        B.bz = (float *)xalloc(n * sizeof(float));
        B.aw = (float *)xalloc(n * sizeof(float));
        B.bw = (float *)xalloc(n * sizeof(float));
    }
    if (g_mat4) B.vin = (float *)xalloc(2u * n * sizeof(float)); /* nel() * 4 floats */
    const bool need2 = with_pipe || g_wide;
    /* The two output sets must be allocated contiguously and exactly n*8 bytes apart. n*8 is
     * a multiple of 4096, so out[0] and out[1] share the same offset relative to a/b modulo
     * 4096 -- otherwise A and B fall into different 4K store/load aliasing relations and can
     * differ by 20% out of nowhere at the L1/L2 tiers (measured, been there). */
    void *blk = xalloc(2u * n * sizeof(fv_vec2));
    void *blk2 = need2 ? xalloc(2u * n * sizeof(fv_vec2)) : nullptr;
    for (int k = 0; k < 2; ++k) {
        B.out[k].base = (char *)blk + (size_t)k * n * sizeof(fv_vec2);
        B.out[k].o = (fv_vec2 *)B.out[k].base;
        B.out[k].ox = (float *)B.out[k].base;
        B.out[k].oy = (float *)B.out[k].base + n;
        B.out[k].lenv = (float *)B.out[k].base;
        B.out[k].dot = 0.0f;
        B.out[k].abs_sum = 0.0;
        if (need2) {
            B.out[k].base2 = (char *)blk2 + (size_t)k * n * sizeof(fv_vec2);
            B.out[k].vx = (float *)B.out[k].base2;
            B.out[k].vy = (float *)B.out[k].base2 + n;
            B.out[k].oz = (float *)B.out[k].base2;
            B.out[k].ow = (float *)B.out[k].base2 + n;
        }
    }
    B.alloc_bytes = n * (2 * sizeof(fv_vec2) + 1 + sizeof(uint32_t) + 4 * sizeof(float));
    B.alloc_bytes += 2u * n * sizeof(fv_vec2);
    if (with_pipe) B.alloc_bytes += n * (4 * sizeof(float) + 2 * sizeof(fv_vec2));
    if (g_mat4) B.alloc_bytes += 2u * n * sizeof(float) + sizeof B.mat;
}

static void free_buf(Buf &B)
{
    std::free(B.a);
    std::free(B.b);
    std::free(B.mask);
    std::free(B.idx);
    std::free(B.ax);
    std::free(B.ay);
    std::free(B.bx);
    std::free(B.by);
    std::free(B.px);
    std::free(B.py);
    std::free(B.qx);
    std::free(B.qy);
    std::free(B.az);
    std::free(B.bz);
    std::free(B.aw);
    std::free(B.bw);
    std::free(B.vin);
    std::free(B.out[0].base); /* both output sets are the two halves of one block */
    std::free(B.out[0].base2);
    std::memset(&B, 0, sizeof B);
}

/* Deterministic refill of the pipeline scratch: collision uses it as normal/velocity input,
 * chain3 as t/u. Because they overwrite each other, px..qy never enter the input checksum. */
static void init_pipe(Buf &B)
{
    if (!B.px) return;
    for (size_t i = 0; i < B.n; ++i) {
        const uint32_t h = (uint32_t)(i * 2654435761u);
        B.px[i] = B.bx[i];
        B.py[i] = B.by[i];
        B.qx[i] = (float)((double)(h & 0xFFFFu) * (1.0 / 65536.0)) - 0.5f;
        B.qy[i] = (float)((double)((h >> 16) & 0xFFFFu) * (1.0 / 65536.0)) - 0.5f;
    }
}

/* The data depends only on n and a fixed seed. One exact zero vector every 512 elements
 * covers normalize's zero-vector guard -- without it an error such as an inverted keep mask
 * would be completely invisible in this fixture. */
static void init_data(Buf &B, size_t n)
{
    g_rng = 0x243f6a8885a308d3ull;
    B.n = n;
    for (size_t i = 0; i < n; ++i) {
        fv_vec2 va, vb;
        if ((i % 512) == 0) {
            va.x = 0.0f;
            va.y = 0.0f;
        } else {
            va.x = rnd_pm1();
            va.y = rnd_pm1();
        }
        vb.x = rnd_pm1();
        vb.y = rnd_pm1();
        B.a[i] = va;
        B.b[i] = vb;
        B.ax[i] = va.x;
        B.ay[i] = va.y;
        B.bx[i] = vb.x;
        B.by[i] = vb.y;
        {
            uint32_t h = (uint32_t)(i * 2654435761u);
            h ^= h >> 15;
            h *= 2246822519u;
            h ^= h >> 13;
            /* Mask density is adjustable (--mask-density) so the sparsity of _mask /
             * _mask_inplace can be swept. The 24-bit mantissa comparison has the same
             * statistics as the old ((h>>11)&1) at density 0.5. */
            B.mask[i] = (uint8_t)(((double)((h >> 8) & 0xFFFFFFu) / 16777216.0 < g_mask_density)
                                      ? 1u
                                      : 0u);
        }
        B.idx[i] = (uint32_t)i;
    }
    for (size_t i = n; i > 1; --i) {
        size_t j = (size_t)(rnd64() % i);
        uint32_t t = B.idx[i - 1];
        B.idx[i - 1] = B.idx[j];
        B.idx[j] = t;
    }
    if (B.az) {
        /* component arrays of wide elements: same deterministic sequence, semantically
         * "component 3/4" */
        for (size_t i = 0; i < n; ++i) {
            const uint32_t h = (uint32_t)(i * 2654435761u);
            B.az[i] = (float)((double)(h & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.bz[i] = (float)((double)((h >> 16) & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.aw[i] = (float)((double)((h ^ 0x5bd1e995u) & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.bw[i] = (float)((double)(((h >> 13) ^ 0x1b873593u) & 0xFFFFu) * (1.0 / 32768.0)) -
                      1.0f;
        }
        /* one exact zero vector every 512 elements (covers normalize's zero guard; vec3 and
         * vec4 get the same treatment) */
        for (size_t i = 0; i < n; i += 512) {
            B.az[i] = 0.0f;
            B.aw[i] = 0.0f;
        }
    }
    init_pipe(B);
    if (B.vin) {
        /* 16 B AoS corpus for the matrix transforms. Drawn last so the vec2/vec3/vec4 corpora
         * above stay bit-identical to what they were before these ops existed.
         * Every component is random, including w: the contract says in.w is *not read* (it is
         * the constant 1 for points, 0 for dirs), so a kernel that accidentally used it would
         * disagree with the reference instead of passing by luck. */
        for (size_t i = 0; i < 2 * n; ++i) B.vin[i] = rnd_pm1();
        for (int k = 0; k < 12; ++k) B.mat[k] = rnd_pm1();
        B.mat[12] = 0.0f; /* affine precondition: last row is (0, 0, 0, 1) */
        B.mat[13] = 0.0f;
        B.mat[14] = 0.0f;
        B.mat[15] = 1.0f;
    }
    for (int k = 0; k < 2; ++k) prefill(B.out[k], n); /* pre-touch: RSS check and warm runs */
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        s += (double)B.a[i].x + (double)B.a[i].y + (double)B.b[i].x + (double)B.b[i].y;
        s += (double)B.mask[i] + (double)B.idx[i];
    }
    if (B.vin) {
        for (size_t i = 0; i < 2 * n; ++i) s += (double)B.vin[i];
        for (int k = 0; k < 16; ++k) s += (double)B.mat[k];
    }
    B.in_cs = s;
}

static void prefill(OutSet &O, size_t n)
{
    std::memset(O.base, 0xDE, n * sizeof(fv_vec2));
    if (O.base2) std::memset(O.base2, 0xDE, n * sizeof(fv_vec2));
    O.dot = -12345.0f;
    O.abs_sum = 0.0;
}

static void drop_pages(Buf &B)
{
    const size_t v2 = B.n * sizeof(fv_vec2);
    madvise(B.a, v2, MADV_DONTNEED);
    madvise(B.b, v2, MADV_DONTNEED);
    madvise(B.mask, B.n, MADV_DONTNEED);
    madvise(B.ax, B.n * sizeof(float), MADV_DONTNEED);
    madvise(B.ay, B.n * sizeof(float), MADV_DONTNEED);
    madvise(B.bx, B.n * sizeof(float), MADV_DONTNEED);
    madvise(B.by, B.n * sizeof(float), MADV_DONTNEED);
    for (int k = 0; k < 2; ++k) {
        madvise(B.out[k].base, v2, MADV_DONTNEED);
        if (B.out[k].base2) madvise(B.out[k].base2, v2, MADV_DONTNEED);
    }
    if (B.px) {
        madvise(B.px, B.n * sizeof(float), MADV_DONTNEED);
        madvise(B.py, B.n * sizeof(float), MADV_DONTNEED);
        madvise(B.qx, B.n * sizeof(float), MADV_DONTNEED);
        madvise(B.qy, B.n * sizeof(float), MADV_DONTNEED);
    }
}

/* =====================================================================
 * 3. Independent reference implementations (pure scalar, sequential, no library batch call)
 * ===================================================================== */
static std::string fbits_str(float f)
{
    unsigned u;
    std::memcpy(&u, &f, sizeof u);
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.9g[0x%08x]", (double)f, u);
    return std::string(buf);
}

static std::string cmp_vec2(const char *what, const fv_vec2 *X, const fv_vec2 *R, size_t n)
{
    size_t bad = 0, first = 0;
    for (size_t i = 0; i < n; ++i) {
        unsigned ux, uy, rx, ry;
        std::memcpy(&ux, &X[i].x, 4);
        std::memcpy(&uy, &X[i].y, 4);
        std::memcpy(&rx, &R[i].x, 4);
        std::memcpy(&ry, &R[i].y, 4);
        if (ux != rx || uy != ry) {
            if (bad == 0) first = i;
            ++bad;
        }
    }
    if (!bad) return std::string();
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "%s mismatch: %zu/%zu element bit patterns differ; first i=%zu got(x,y)=%s,%s want=%s,%s",
                  what, bad, n, first, fbits_str(X[first].x).c_str(),
                  fbits_str(X[first].y).c_str(), fbits_str(R[first].x).c_str(),
                  fbits_str(R[first].y).c_str());
    return std::string(buf);
}

static std::string cmp_floats(const char *what, const float *X, const float *R, size_t n)
{
    size_t bad = 0, first = 0;
    for (size_t i = 0; i < n; ++i) {
        unsigned ux, rx;
        std::memcpy(&ux, &X[i], 4);
        std::memcpy(&rx, &R[i], 4);
        if (ux != rx) {
            if (bad == 0) first = i;
            ++bad;
        }
    }
    if (!bad) return std::string();
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "%s mismatch: %zu/%zu value bit patterns differ; first i=%zu got=%s want=%s",
                  what, bad, n,
                  first, fbits_str(X[first]).c_str(), fbits_str(R[first]).c_str());
    return std::string(buf);
}

static inline void scalar_normalize(float x, float y, float *ox, float *oy)
{
    const float len = std::sqrt(x * x + y * y);
    if (len <= FV_EPS) {
        *ox = 0.0f;
        *oy = 0.0f;
    } else {
        const float inv = 1.0f / len;
        *ox = x * inv;
        *oy = y * inv;
    }
}
static inline void scalar_reflect(float nx, float ny, float vx, float vy, float *ox, float *oy)
{
    const float d = (vx * nx + vy * ny) * 2.0f;
    *ox = vx - d * nx;
    *oy = vy - d * ny;
}


/* ================= vec3 / vec4 references (pure scalar, sequential, no library) ========= */
#if AB_OPS_EXTRA /* these entries are off the registry by default; the references go too */
/* vec3: 12 B elements, three component arrays. The zero guard matches
 * fv_vec2_normalize (len<=FV_EPS -> 0). */
static void ref_v3_add(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.ox[i] = B.ax[i] + B.bx[i];
        O.oy[i] = B.ay[i] + B.by[i];
        O.oz[i] = B.az[i] + B.bz[i];
    }
}
static void ref_v3_scale(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.ox[i] = B.ax[i] * AB_S;
        O.oy[i] = B.ay[i] * AB_S;
        O.oz[i] = B.az[i] * AB_S;
    }
}
static void ref_v3_dot(Buf &B, size_t n, OutSet &O)
{
    double s = 0.0, a = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double t = (double)B.ax[i] * (double)B.bx[i] + (double)B.ay[i] * (double)B.by[i] +
                         (double)B.az[i] * (double)B.bz[i];
        s += t;
        a += std::fabs(t);
    }
    O.dot = (float)s;
    O.abs_sum = a;
}
static void ref_v3_len(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i)
        O.lenv[i] = std::sqrt(B.ax[i] * B.ax[i] + B.ay[i] * B.ay[i] + B.az[i] * B.az[i]);
}
static void ref_v3_norm(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        const float x = B.ax[i], y = B.ay[i], z = B.az[i];
        const float len = std::sqrt(x * x + y * y + z * z);
        if (len <= FV_EPS) {
            O.ox[i] = 0.0f;
            O.oy[i] = 0.0f;
            O.oz[i] = 0.0f;
        } else {
            const float inv = 1.0f / len;
            O.ox[i] = x * inv;
            O.oy[i] = y * inv;
            O.oz[i] = z * inv;
        }
    }
}
/* vec4 AoS (16 B/element): two elements fill one ymm, hence _x2. */
static void ref_v4_add_x2(Buf &B, size_t n, OutSet &O)
{
    const fv_vec4 *a = (const fv_vec4 *)B.a, *b = (const fv_vec4 *)B.b;
    fv_vec4 *o = (fv_vec4 *)O.o;
    for (size_t i = 0; i < n; ++i) {
        o[i].x = a[i].x + b[i].x;
        o[i].y = a[i].y + b[i].y;
        o[i].z = a[i].z + b[i].z;
        o[i].w = a[i].w + b[i].w;
    }
}
static void ref_v4_sub_x2(Buf &B, size_t n, OutSet &O)
{
    const fv_vec4 *a = (const fv_vec4 *)B.a, *b = (const fv_vec4 *)B.b;
    fv_vec4 *o = (fv_vec4 *)O.o;
    for (size_t i = 0; i < n; ++i) {
        o[i].x = a[i].x - b[i].x;
        o[i].y = a[i].y - b[i].y;
        o[i].z = a[i].z - b[i].z;
        o[i].w = a[i].w - b[i].w;
    }
}
static void ref_v4_scale_x2(Buf &B, size_t n, OutSet &O)
{
    const fv_vec4 *a = (const fv_vec4 *)B.a;
    fv_vec4 *o = (fv_vec4 *)O.o;
    for (size_t i = 0; i < n; ++i) {
        o[i].x = a[i].x * AB_S;
        o[i].y = a[i].y * AB_S;
        o[i].z = a[i].z * AB_S;
        o[i].w = a[i].w * AB_S;
    }
}
static void ref_v4_lerp_x2(Buf &B, size_t n, OutSet &O)
{
    const fv_vec4 *a = (const fv_vec4 *)B.a, *b = (const fv_vec4 *)B.b;
    fv_vec4 *o = (fv_vec4 *)O.o;
    for (size_t i = 0; i < n; ++i) {
        o[i].x = a[i].x + (b[i].x - a[i].x) * AB_T;
        o[i].y = a[i].y + (b[i].y - a[i].y) * AB_T;
        o[i].z = a[i].z + (b[i].z - a[i].z) * AB_T;
        o[i].w = a[i].w + (b[i].w - a[i].w) * AB_T;
    }
}
static void ref_v4_madd_x2(Buf &B, size_t n, OutSet &O)
{
    const fv_vec4 *a = (const fv_vec4 *)B.a, *b = (const fv_vec4 *)B.b;
    fv_vec4 *o = (fv_vec4 *)O.o;
    for (size_t i = 0; i < n; ++i) {
        o[i].x = a[i].x + b[i].x * AB_S;
        o[i].y = a[i].y + b[i].y * AB_S;
        o[i].z = a[i].z + b[i].z * AB_S;
        o[i].w = a[i].w + b[i].w * AB_S;
    }
}
static void ref_v4_dot(Buf &B, size_t n, OutSet &O)
{
    double s = 0.0, a = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double t = (double)B.ax[i] * (double)B.bx[i] + (double)B.ay[i] * (double)B.by[i] +
                         (double)B.az[i] * (double)B.bz[i] + (double)B.aw[i] * (double)B.bw[i];
        s += t;
        a += std::fabs(t);
    }
    O.dot = (float)s;
    O.abs_sum = a;
}
static void ref_v4_len(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i)
        O.lenv[i] = std::sqrt(B.ax[i] * B.ax[i] + B.ay[i] * B.ay[i] + B.az[i] * B.az[i] +
                              B.aw[i] * B.aw[i]);
}
static void ref_v4_norm(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        const float x = B.ax[i], y = B.ay[i], z = B.az[i], w = B.aw[i];
        const float len = std::sqrt(x * x + y * y + z * z + w * w);
        if (len <= FV_EPS) {
            O.ox[i] = 0.0f;
            O.oy[i] = 0.0f;
            O.oz[i] = 0.0f;
            O.ow[i] = 0.0f;
        } else {
            const float inv = 1.0f / len;
            O.ox[i] = x * inv;
            O.oy[i] = y * inv;
            O.oz[i] = z * inv;
            O.ow[i] = w * inv;
        }
    }
}

#endif /* AB_OPS_EXTRA */

/* ---- per-op reference implementation (one signature, stored in the OPS table) ---- */
/* Independent scalar references for the matrix transforms. The header promises bit-identical
 * results under exactly two conditions, and both are reproduced here or EXACT could not hold:
 *   1) the evaluation order is ((m0*x + m1*y) + m2*z) + m3;
 *   2) no FMA (the build uses -ffp-contract=off and the kernels live in a matching pragma);
 * dirs simply drops the translation term instead of adding m3 * 0 -- x + (+0.0) turns -0.0
 * into +0.0, and the kernel does not perform that add either. out.w is written as the
 * documented constant rather than computed. */
static void ref_mat4_points(Buf &B, size_t n, OutSet &O)
{
    const float *m = B.mat;
    const float *in = B.vin;
    float *out = (float *)(void *)O.o;
    for (size_t i = 0; i < n; ++i) {
        const float x = in[4 * i + 0];
        const float y = in[4 * i + 1];
        const float z = in[4 * i + 2];
        out[4 * i + 0] = ((m[0] * x + m[1] * y) + m[2] * z) + m[3];
        out[4 * i + 1] = ((m[4] * x + m[5] * y) + m[6] * z) + m[7];
        out[4 * i + 2] = ((m[8] * x + m[9] * y) + m[10] * z) + m[11];
        out[4 * i + 3] = 1.0f;
    }
}
static void ref_mat4_dirs(Buf &B, size_t n, OutSet &O)
{
    const float *m = B.mat;
    const float *in = B.vin;
    float *out = (float *)(void *)O.o;
    for (size_t i = 0; i < n; ++i) {
        const float x = in[4 * i + 0];
        const float y = in[4 * i + 1];
        const float z = in[4 * i + 2];
        out[4 * i + 0] = (m[0] * x + m[1] * y) + m[2] * z;
        out[4 * i + 1] = (m[4] * x + m[5] * y) + m[6] * z;
        out[4 * i + 2] = (m[8] * x + m[9] * y) + m[10] * z;
        out[4 * i + 3] = 0.0f;
    }
}

typedef void (*RefFn)(Buf &B, size_t n, OutSet &O);

static void ref_add_x4(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x + B.b[i].x;
        O.o[i].y = B.a[i].y + B.b[i].y;
    }
}
static void ref_sub_x4(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x - B.b[i].x;
        O.o[i].y = B.a[i].y - B.b[i].y;
    }
}
static void ref_scale_x4(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x * AB_S;
        O.o[i].y = B.a[i].y * AB_S;
    }
}
static void ref_lerp_x4(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x + (B.b[i].x - B.a[i].x) * AB_T;
        O.o[i].y = B.a[i].y + (B.b[i].y - B.a[i].y) * AB_T;
    }
}
static void ref_madd_x4(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x + B.b[i].x * AB_S;
        O.o[i].y = B.a[i].y + B.b[i].y * AB_S;
    }
}
static void ref_add_mask(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        if (B.mask[i] != 0) {
            O.o[i].x = B.a[i].x + B.b[i].x;
            O.o[i].y = B.a[i].y + B.b[i].y;
        } else {
            O.o[i] = B.a[i];
        }
    }
}
/* In-place sparse add: o[i] = a[i] + b[i] where mask[i] != 0, else o[i] = a[i] (unchanged).
 * The driver passes O.o as the writable array (see the AB_OPF_INPLACE branch of exec).
 * For now it only serves the experimental _mask_inplace row; drop the #if below when that row
 * moves into the default table. */
#if AB_OPS_EXTRA
static void ref_add_mask_inplace(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        if (B.mask[i] != 0) {
            O.o[i].x = B.a[i].x + B.b[i].x;
            O.o[i].y = B.a[i].y + B.b[i].y;
        } else {
            O.o[i] = B.a[i];
        }
    }
}
#endif /* AB_OPS_EXTRA */
static void ref_add_idx(Buf &B, size_t n, OutSet &O)
{
    for (size_t k = 0; k < n; ++k) {
        const uint32_t i = B.idx[k];
        O.o[i].x = B.a[i].x + B.b[i].x;
        O.o[i].y = B.a[i].y + B.b[i].y;
    }
}
static void ref_scale_idx(Buf &B, size_t n, OutSet &O)
{
    for (size_t k = 0; k < n; ++k) {
        const uint32_t i = B.idx[k];
        O.o[i].x = B.a[i].x * AB_S;
        O.o[i].y = B.a[i].y * AB_S;
    }
}
static void ref_add_soa(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.ox[i] = B.ax[i] + B.bx[i];
        O.oy[i] = B.ay[i] + B.by[i];
    }
}
static void ref_scale_soa(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.ox[i] = B.ax[i] * AB_S;
        O.oy[i] = B.ay[i] * AB_S;
    }
}
static void ref_dot_soa(Buf &B, size_t n, OutSet &O)
{
    /* Accumulate in double while also recording sum|term|: the tolerance has to follow the
     * magnitude of the terms, not of the result (when the result nearly cancels out a relative
     * error explodes -- measured 7e-5 already at n=2^18). */
    double s = 0.0, a = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double t = (double)B.ax[i] * (double)B.bx[i] + (double)B.ay[i] * (double)B.by[i];
        s += t;
        a += std::fabs(t);
    }
    O.dot = (float)s;
    O.abs_sum = a;
}
static void ref_len_soa(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) O.lenv[i] = std::sqrt(B.ax[i] * B.ax[i] + B.ay[i] * B.ay[i]);
}
static void ref_norm_soa(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) scalar_normalize(B.ax[i], B.ay[i], &O.ox[i], &O.oy[i]);
}
static void ref_collision(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        float x, y;
        scalar_normalize(B.px[i], B.py[i], &x, &y);
        O.ox[i] = x;
        O.oy[i] = y;
        scalar_reflect(x, y, B.qx[i], B.qy[i], &O.vx[i], &O.vy[i]);
    }
}
static void ref_chain3(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        const float tx = B.ax[i] * AB_S, ty = B.ay[i] * AB_S;
        const float ux = tx + B.bx[i], uy = ty + B.by[i];
        scalar_normalize(ux, uy, &O.ox[i], &O.oy[i]);
    }
}
/* ---- ceiling references (bypass the variant TU; A/B runs identical code, so the ratio is a
 * pure noise probe) ---- */
static void mem_memcpy(Buf &B, size_t n, OutSet &O)
{
    std::memcpy(O.o, B.a, n * sizeof(fv_vec2));
}
static void mem_set(Buf &B, size_t n, OutSet &O)
{
    (void)B;
    std::memset(O.o, 0, n * sizeof(fv_vec2));
}
static void mem_read(Buf &B, size_t n, OutSet &O)
{
    const uint32_t *__restrict x = (const uint32_t *)(const void *)B.a;
    uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0;
    const size_t m = 2 * n, lim = m & ~(size_t)7;
    size_t i = 0;
    for (; i < lim; i += 8) {
        s0 += x[i + 0];
        s1 += x[i + 1];
        s2 += x[i + 2];
        s3 += x[i + 3];
        s4 += x[i + 4];
        s5 += x[i + 5];
        s6 += x[i + 6];
        s7 += x[i + 7];
    }
    for (; i < m; ++i) s0 += x[i];
    O.dot = (float)((double)(((s0 + s1) + (s2 + s3)) + ((s4 + s5) + (s6 + s7))) * 1e-12);
}
static void mem_triad(Buf &B, size_t n, OutSet &O)
{
    for (size_t i = 0; i < n; ++i) {
        O.o[i].x = B.a[i].x + B.b[i].x;
        O.o[i].y = B.a[i].y + B.b[i].y;
    }
}

/* =====================================================================
 * 4. The op table
 * ===================================================================== */
struct Op {
    const char *name;
    size_t rd, wr; /* nominal read/write bytes per element (real traffic: real_wr_bytes) */
    size_t elop;   /* element width in bytes: vec2 8 / vec3 12 / vec4 16 */
    OutKind kind;
    Policy pol;
    RefFn ref;       /* reference implementation; null for POL_NOCHECK */
    RefFn local_run; /* implementation that bypasses the variant TU (ceiling); null = ab_X_loop */
    int abi_op;      /* AB_OP_*; -1 means local_run */
    unsigned flags;  /* AB_OPF_* / AB_MF_* */
    int rfo;         /* AB_RFO_*: write-allocate behaviour of the plain store (real traffic) */
    const char *note;/* human-readable basis of the RFO accounting */
};

/* The table is expanded from the ab_ops.h registry: the driver holds no per-op code, so
 * adding an op never touches this file. */
static const Op OPS[] = {
#define AB_LIB_ROW(ID, NAME, KIND, ELOP, POL, REF, FLAGS, RD, WR, RFO, NOTE)                    \
    {NAME, RD, WR, ELOP, OUT_##KIND, POL_##POL, REF, nullptr, AB_OP_##ID, FLAGS, RFO, NOTE},
#define AB_MEM_ROW(NAME, KIND, LOCAL, FLAGS, RD, WR, RFO, NOTE)                                  \
    {NAME, RD, WR, 8, OUT_##KIND, POL_NOCHECK, nullptr, LOCAL, -1, FLAGS, RFO, NOTE},
    AB_LIB_OPS(AB_LIB_ROW)
    AB_LIB_OPS_EXTRA(AB_LIB_ROW)
    AB_MEM_OPS(AB_MEM_ROW)
#undef AB_LIB_ROW
#undef AB_MEM_ROW
};
static const size_t OP_COUNT = sizeof(OPS) / sizeof(OPS[0]);

/* Follows registry additions/removals instead of hard-coding numbers (19 used to be 15 batch
 * entries + 4 ceiling references) */
#define AB_COUNT_ITEM(...) +1
enum {
    N_LIB_OPS = 0 AB_LIB_OPS(AB_COUNT_ITEM) AB_LIB_OPS_EXTRA(AB_COUNT_ITEM),
    N_MEM_OPS = 0 AB_MEM_OPS(AB_COUNT_ITEM)
};
#undef AB_COUNT_ITEM
static_assert(sizeof(OPS) / sizeof(OPS[0]) == N_LIB_OPS + N_MEM_OPS,
              "the OPS table and the registry disagree");

static bool has_flag(const Op &op, unsigned f)
{
    return (op.flags & f) != 0;
}

/* Tiers to exclude from the per-size criterion (comma separated; default 131072 = the L3
 * boundary, where the A/A floor itself is unresolvable) */
static bool tier_excluded(size_t n, const std::string &list)
{
    size_t p = 0;
    while (p <= list.size()) {
        const size_t c = list.find(',', p);
        const std::string t = list.substr(p, c == std::string::npos ? std::string::npos : c - p);
        if (!t.empty() && (size_t)std::strtoul(t.c_str(), nullptr, 0) == n) return true;
        if (c == std::string::npos) break;
        p = c + 1;
    }
    return false;
}

/* How many **of its own elements** this op processes at sweep size n (8 B units).
 * n*8 is the byte volume every op shares, so ops of different element width are
 * byte-normalized. */
static size_t nel(const Op &op, size_t n)
{
    return op.elop ? (n * 8u) / op.elop : n;
}

/* Ceiling references are picked by **role**, not by index (inserting a row cannot shift them) */
static size_t memref_index(unsigned role)
{
    for (size_t k = 0; k < OP_COUNT; ++k)
        if (OPS[k].abi_op < 0 && (OPS[k].flags & role)) return k;
    return OP_COUNT;
}

/* INTRA pairs (AB_INTRA_PAIRS in ab_ops.h): new entry ID -> reference entry ID */
struct IntraPair {
    int a, b;
};
static const IntraPair kIntra[] = {
#define AB_INTRA_ITEM(A, B) {AB_OP_##A, AB_OP_##B},
    AB_INTRA_PAIRS(AB_INTRA_ITEM)
#undef AB_INTRA_ITEM
    {-1, -1} /* sentinel: the table must stay legal even when AB_OPS_EXTRA=0 empties it */
};

static size_t op_index_of_id(int id)
{
    for (size_t k = 0; k < OP_COUNT; ++k)
        if (OPS[k].abi_op == id) return k;
    return OP_COUNT;
}

static size_t intra_ref_of(size_t opi)
{
    for (const IntraPair &p : kIntra)
        if (p.a >= 0 && p.a == OPS[opi].abi_op) return op_index_of_id(p.b);
    return OP_COUNT;
}

static unsigned side_opflags(int side, int abi_op)
{
    if (abi_op < 0) return 0u;
    if (side == SIDE_A) return ab_a_opflags(abi_op);
    if (side == SIDE_B) return ab_b_opflags(abi_op);
    return ab_s_opflags(abi_op);
}

/* The **real** write traffic of the output array (B/element):
 *   NT store                -> wr      (bypasses the cache, no RFO)
 *   target was already read -> wr      (in-place update; the line is already cached)
 *   libc picks for itself   -> wr      (uncertain, flagged libc?, never judged)
 *   plain store             -> wr * 2  (one write + one RFO pulling the line in)
 * Input (read) traffic is always counted 1x. */
static size_t real_wr_bytes(const Op &op, int side)
{
    if (op.wr == 0) return 0;
    if (side_opflags(side, op.abi_op) & AB_SIDE_NT) return op.wr;
    if (op.rfo == AB_RFO_HIT || op.rfo == AB_RFO_LIBC) return op.wr;
    return op.wr * 2;
}

static const char *rfo_side_name(const Op &op, int side)
{
    if (op.wr == 0) return "-";
    if (side_opflags(side, op.abi_op) & AB_SIDE_NT) return "NT";
    switch (op.rfo) {
    case AB_RFO_HIT:
        return "hit";
    case AB_RFO_MIXED:
        return "mixed?";
    case AB_RFO_LIBC:
        return "libc?";
    case AB_RFO_NT:
        return "miss*"; /* this side has no NT entry, so it degrades to a plain store */
    default:
        return "miss";
    }
}

static void run_side(int side, int abi_op, size_t n, long iters, const AbArgs *p, float *dot)
{
    if (side == SIDE_A)
        ab_a_loop(abi_op, n, iters, p, dot);
    else if (side == SIDE_B)
        ab_b_loop(abi_op, n, iters, p, dot);
    else
        ab_s_loop(abi_op, n, iters, p, dot);
}

static AbArgs mk_args(Buf &B, OutSet &O)
{
    AbArgs p{};
    p.a = B.a;
    p.b = B.b;
    p.o = O.o;
    p.mask = B.mask;
    p.idx = B.idx;
    p.ax = B.ax;
    p.ay = B.ay;
    p.bx = B.bx;
    p.by = B.by;
    p.ox = O.ox;
    p.oy = O.oy;
    p.lenv = O.lenv;
    p.ox2 = O.vx;
    p.oy2 = O.vy;
    p.px = B.px;
    p.py = B.py;
    p.qx = B.qx;
    p.qy = B.qy;
    p.az = B.az;
    p.bz = B.bz;
    p.aw = B.aw;
    p.bw = B.bw;
    p.oz = O.oz;
    p.ow = O.ow;
    p.mat = B.mat;
    p.vin = B.vin;
    p.vout = (void *)O.o; /* AoS output, same 8n bytes as the vec2 output buffer */
    return p;
}

/* In-place ops: the output array is first filled with a's contents and then **is** the
 * writable array itself. That way B.a is never written and the "input checksum unchanged"
 * gate holds for in-place ops too. */
static void seed_inplace(const Op &op, Buf &B, OutSet &O)
{
    if (!has_flag(op, AB_OPF_INPLACE)) return;
    std::memcpy(O.base, B.a, B.n * sizeof(fv_vec2));
}

static void exec(const Op &op, int side, Buf &B, OutSet &O, long iters)
{
    if (op.kind == OUT_NONE) O.dot = 0.0f; /* reduction ops use dot as the accumulator */
    if (op.abi_op >= 0) {
        AbArgs p = mk_args(B, O);
        if (has_flag(op, AB_OPF_INPLACE)) p.a = O.o;
        run_side(side, op.abi_op, nel(op, B.n), iters, &p, &O.dot);
        return;
    }
    const size_t m = nel(op, B.n);
    for (long k = 0; k < iters; ++k) op.local_run(B, m, O);
}

/* =====================================================================
 * 5. Checksums and gates
 * ===================================================================== */
static double out_sum(const Op &op, const OutSet &O, size_t n)
{
    double s = 0.0;
    switch (op.kind) {
    case OUT_VEC2:
        for (size_t i = 0; i < n; ++i) s += (double)O.o[i].x + (double)O.o[i].y;
        break;
    case OUT_SOA2:
        for (size_t i = 0; i < n; ++i) s += (double)O.ox[i] + (double)O.oy[i];
        break;
    case OUT_SOA3:
        for (size_t i = 0; i < n; ++i)
            s += (double)O.ox[i] + (double)O.oy[i] + (double)O.oz[i];
        break;
    case OUT_SOA4:
        for (size_t i = 0; i < n; ++i)
            s += (double)O.ox[i] + (double)O.oy[i] + (double)O.vx[i] + (double)O.vy[i];
        break;
    case OUT_SOA4V:
        for (size_t i = 0; i < n; ++i)
            s += (double)O.ox[i] + (double)O.oy[i] + (double)O.oz[i] + (double)O.ow[i];
        break;
    case OUT_VEC4:
    case OUT_VEC4M: {
        const fv_vec4 *v = (const fv_vec4 *)O.o;
        for (size_t i = 0; i < n; ++i)
            s += (double)v[i].x + (double)v[i].y + (double)v[i].z + (double)v[i].w;
        break;
    }
    case OUT_FLOAT:
        for (size_t i = 0; i < n; ++i) s += (double)O.lenv[i];
        break;
    case OUT_NONE:
        s = (double)O.dot;
        break;
    }
    return s;
}

static std::string cmp_out(const Op &op, const OutSet &X, const OutSet &R, size_t n)
{
    switch (op.kind) {
    case OUT_VEC2:
        return cmp_vec2("o[]", X.o, R.o, n);
    case OUT_SOA2: {
        std::string m1 = cmp_floats("ox[]", X.ox, R.ox, n);
        if (!m1.empty()) return m1;
        return cmp_floats("oy[]", X.oy, R.oy, n);
    }
    case OUT_SOA4: {
        std::string m1 = cmp_floats("ox[]", X.ox, R.ox, n);
        if (!m1.empty()) return m1;
        std::string m2 = cmp_floats("oy[]", X.oy, R.oy, n);
        if (!m2.empty()) return m2;
        std::string m3 = cmp_floats("vox[]", X.vx, R.vx, n);
        if (!m3.empty()) return m3;
        return cmp_floats("voy[]", X.vy, R.vy, n);
    }
    case OUT_SOA3: {
        std::string m1 = cmp_floats("ox[]", X.ox, R.ox, n);
        if (!m1.empty()) return m1;
        std::string m2 = cmp_floats("oy[]", X.oy, R.oy, n);
        if (!m2.empty()) return m2;
        return cmp_floats("oz[]", X.oz, R.oz, n);
    }
    case OUT_SOA4V: {
        std::string m1 = cmp_floats("ox[]", X.ox, R.ox, n);
        if (!m1.empty()) return m1;
        std::string m2 = cmp_floats("oy[]", X.oy, R.oy, n);
        if (!m2.empty()) return m2;
        std::string m3 = cmp_floats("oz[]", X.oz, R.oz, n);
        if (!m3.empty()) return m3;
        return cmp_floats("ow[]", X.ow, R.ow, n);
    }
    case OUT_VEC4:
    case OUT_VEC4M:
        /* fv_vec4 is 4 floats, so a bitwise compare is a compare of 4n floats */
        return cmp_floats("o[].xyzw", (const float *)(const void *)X.o,
                          (const float *)(const void *)R.o, 4 * n);
    case OUT_FLOAT:
        return cmp_floats("o[]", X.lenv, R.lenv, n);
    case OUT_NONE: {
        const double tol =
            std::max(AB_REL_FRAC * std::fabs((double)R.dot), AB_ABS_FRAC * R.abs_sum);
        const double d = std::fabs((double)X.dot - (double)R.dot);
        if (d <= tol) return std::string();
        char buf[512];
        std::snprintf(buf, sizeof buf,
                      "out of tolerance: got=%s want=%s |diff|=%.4g > tol=%.4g (=max(%.0e*|r|, "
                      "%.0e*sum|term|))",
                      fbits_str(X.dot).c_str(), fbits_str(R.dot).c_str(), d, tol,
                      (double)AB_REL_FRAC, (double)AB_ABS_FRAC);
        return std::string(buf);
    }
    }
    return std::string();
}

struct Gate {
    bool ok;
    double cs_a, cs_b, cs_ref;
    std::string why;
};

static Gate check_op(const Op &op, Buf &B, size_t n, bool with_scalar)
{
    Gate g;
    g.ok = true;
    g.cs_a = g.cs_b = g.cs_ref = 0.0;
    OutSet &P = B.out[0];
    OutSet &R = B.out[1];

    if (op.pol == POL_NOCHECK) { /* ceiling reference: A/B run the same code, check bits only */
        prefill(P, n);
        exec(op, SIDE_A, B, P, 1);
        prefill(R, n);
        exec(op, SIDE_B, B, R, 1);
        g.cs_a = out_sum(op, P, nel(op, n));
        g.cs_b = out_sum(op, R, nel(op, n));
        std::string m = cmp_out(op, P, R, nel(op, n));
        if (!m.empty()) {
            g.ok = false;
            g.why = "  A and B outputs differ (a ceiling reference should be identical): " + m +
                    "\n";
        }
        g.cs_ref = g.cs_a;
        return g;
    }

    const size_t ne = nel(op, n); /* element count of this op (byte volume is constant) */
    prefill(R, n);
    op.ref(B, ne, R);

    if (has_flag(op, AB_OPF_NO_BASELINE)) {
        /* The baseline (git HEAD) has no such entry at all: side A does not take part, only
         * "candidate vs independent scalar reference" is checked. The most this path can say
         * is "the implementation matches the scalar semantics", never "it is faster than the
         * baseline". */
        const double abs0 = R.abs_sum;
        prefill(P, n);
        exec(op, SIDE_B, B, P, 1);
        g.cs_b = g.cs_a = out_sum(op, P, ne);
        g.cs_ref = out_sum(op, R, ne);
        std::string mb = cmp_out(op, P, R, ne);
        if (!mb.empty()) {
            g.ok = false;
            g.why += "  candidate vs reference: " + mb + "\n";
        }
        R.abs_sum = abs0;
        return g;
    }
    g.cs_ref = out_sum(op, R, ne);
    const double abs_sum = R.abs_sum;

    prefill(P, n);
    seed_inplace(op, B, P);
    exec(op, SIDE_A, B, P, 1);
    g.cs_a = out_sum(op, P, ne);
    std::string m_ar = cmp_out(op, P, R, ne);

    prefill(P, n);
    seed_inplace(op, B, P);
    exec(op, SIDE_B, B, P, 1);
    g.cs_b = out_sum(op, P, ne);
    std::string m_br = cmp_out(op, P, R, ne);

    if (!m_ar.empty() || !m_br.empty()) {
        g.ok = false;
        if (!m_ar.empty()) g.why += "  baseline vs reference: " + m_ar + "\n";
        if (!m_br.empty()) g.why += "  candidate vs reference: " + m_br + "\n";
    }
    if (g.ok && with_scalar && op.kind != OUT_SOA4) { /* the scalar anchor goes through the gate too */
        OutSet S = R;
        S.dot = R.dot;
        prefill(S, n);
        exec(op, SIDE_S, B, S, 1);
        std::string ms = cmp_out(op, S, R, ne);
        if (!ms.empty()) {
            g.ok = false;
            g.why += "  scalar_loop anchor vs reference: " + ms + "\n";
        }
    }
    R.abs_sum = abs_sum;
    return g;
}

/* =====================================================================
 * 6. Timing and statistics
 * ===================================================================== */
static double median(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t m = v.size();
    if (m & 1u) return v[m / 2];
    return 0.5 * (v[m / 2 - 1] + v[m / 2]);
}

/* Interquartile range of the paired ratios (relative to the median). Far more robust than
 * the sample CoV: an occasional multi-x stall on this machine blows up std, so a CoV test
 * declares every point unresolvable while the ratio distribution is actually tight.
 * Resolvability is therefore decided by this function; CoV is still reported (the metric the
 * parent task asked for by name). */
static double iqr_pct(std::vector<double> v)
{
    if (v.size() < 4) return 0.0;
    std::sort(v.begin(), v.end());
    const double q25 = v[v.size() / 4];
    const double q75 = v[(v.size() * 3) / 4];
    const double m = median(v);
    return (m > 0.0) ? (q75 - q25) / m * 100.0 : 0.0;
}

/* Quantile (linear interpolation). The floor needs p50/p90/max -- max alone is not enough. */
static double quantile(std::vector<double> v, double q)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    if (v.size() == 1) return v[0];
    const double pos = q * (double)(v.size() - 1);
    const size_t lo = (size_t)pos;
    const size_t hi = (lo + 1 < v.size()) ? lo + 1 : lo;
    const double f = pos - (double)lo;
    return v[lo] * (1.0 - f) + v[hi] * f;
}

static double cov_pct(const std::vector<double> &v)
{
    if (v.size() < 2) return 0.0;
    double m = 0.0;
    for (double x : v) m += x;
    m /= (double)v.size();
    if (m <= 0.0) return 0.0;
    double s = 0.0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (double)(v.size() - 1)) / m * 100.0;
}

/* Returns ns/element (or ns/call: the caller multiplies/divides as needed) */
static double measure(const Op &op, int side, Buf &B, OutSet &O, long iters)
{
    const double t0 = now_ns();
    exec(op, side, B, O, iters);
    const double t1 = now_ns();
    const size_t m = nel(op, B.n); /* element count from this op width (byte volume constant) */
    return (t1 - t0) / (double)iters / (double)(m ? m : 1);
}

#define AB_F_COV 0x1u
#define AB_F_LOAD 0x2u
#define AB_F_RSS 0x4u

struct Pt {
    const char *op;
    size_t opi;
    size_t n;
    bool cold;
    double a_med, b_med, ratio, ratio_m, pr_min, pr_max, a_cov, b_cov, riqr;
    double a_min, a_max, b_min, b_max;
    double gbs_a, gbs_b, pct_memcpy_a, pct_memcpy_b;
    /* RFO accounting: real traffic = rd + wr*(RFO ? 2 : 1). Computed per side, because one
     * side may use an NT store while the other degrades to a plain store (see *_opflags in
     * ab_variant.cpp). */
    double real_gbs_a, real_gbs_b, real_pct_memcpy_a, real_pct_memcpy_b;
    size_t real_rd, real_wr_a, real_wr_b;
    int rfo; /* AB_RFO_* declaration of the row (whether NT really happened: rfo_side_name) */
    const char *mode;   /* "ab" (A/B) or "intra" (control inside the candidate) */
    const char *ref_op; /* reference op name in intra mode; "-" in ab mode */
    /* A/A pairs interleaved with A/B inside the same round (R4 denominator: same
     * distribution, same data, same time window) */
    int aa_n;
    double aa_p50, aa_p90, aa_max, aa_dev_max;
    int trials;
    long iters;
    double cs_a, cs_b, cs_ref;
    double load0, load1;
    unsigned flags;
    size_t rd, wr;
};

/* Ops covered by the small-N regime: those flagged AB_OPF_SMALL in the registry (the scalar
 * anchor ab_s_loop implements them). A flag rather than an index, so inserting rows or adding
 * ops can never silently change the coverage. */
static bool small_supported(int opi)
{
    return (size_t)opi < OP_COUNT && (OPS[opi].flags & AB_OPF_SMALL) != 0;
}

static std::string checksum_text(const Op &op, double a, double b, double r)
{
    char buf[256];
    if (op.pol == POL_EXACT) {
        std::snprintf(buf, sizeof buf, "cs=%.7g A=B=ref", r);
    } else if (op.pol == POL_REL) {
        std::snprintf(buf, sizeof buf, "cs=%.7g/%.7g/%.7g", a, b, r);
    } else {
        std::snprintf(buf, sizeof buf, "cs=%.7g (A=B)", a);
    }
    return std::string(buf);
}

static const char *tier_of(size_t ws)
{
    if (ws <= AB_L1_BYTES) return "L1d";
    if (ws <= AB_L2_BYTES) return "L2";
    if (ws <= AB_L3_BYTES) return "L3";
    return "DRAM";
}

struct Args {
    std::vector<size_t> sizes;
    std::vector<size_t> small_ns;
    std::string ops_filter;
    std::string op_exact; /* --op=<exact name>: used by verify to avoid substring hits */
    std::string label;
    std::string tsv_out, matrix_out, small_out, floor_out;
    int trials_min, trials_max;
    double cov_target, target_ms, small_target_ms;
    bool selfcheck, do_cold, do_pipe, do_small;
    bool gate_only; /* bit-level gates only, no timing output at all (R1 of verify) */
    bool list_ops;  /* print the registry and exit (verify uses it to test membership) */
    /* one A/A pair right after every A/B round: the R4 denominator must share the A/B
     * distribution (same process, same data, same time window) */
    bool aa_each_round;
    bool intra;     /* INTRA mode: new entry vs a same-byte-volume control op in the candidate */
    double mask_density; /* mask density of _mask / _mask_inplace */
    std::string tier_exclude; /* tiers excluded from the per-size criterion (default 131072) */
    int order;    /* 0=interleaved 1=always A first (default) 2=always B first */
    bool swap_out; /* swap the output buffers used by A/B (diagnose buffer asymmetry) */
};

static bool parse_args(int argc, char **argv, Args &A)
{
    A.sizes = {4096, 16384, 65536, 131072, 262144, 1048576, 4194304, 16777216};
    A.small_ns = {4, 8, 16, 32, 64};
    A.trials_min = 5;
    A.trials_max = 15;
    A.cov_target = 5.0;
    A.target_ms = 8.0;
    A.small_target_ms = 20.0;
    A.selfcheck = false;
    A.do_cold = true;
    A.do_pipe = true;
    A.do_small = true;
    A.gate_only = false;
    A.list_ops = false;
    A.aa_each_round = false;
    A.intra = false;
    A.mask_density = 0.5;
    A.tier_exclude = "131072";
    A.label = AB_LABEL;
    /* Default is a fixed "A then B". Why (measured): at 2^20 order=ab gives ratio 0.995 /
     * CoV 9.7%, order=ba gives 1.030 / 34.5% and interleaving gives 1.305 / 24.3%. A fixed
     * order removes the confounder "each pair enters in one of two different states", so the
     * paired-ratio distribution is no longer a mixture of two populations; and since ab and ba
     * both land near 1.0, order itself carries no systematic bias. Use --order=ba to
     * cross-check. */
    A.order = 1;
    A.swap_out = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto val = [&](const char *pfx) -> std::string {
            const size_t l = std::strlen(pfx);
            return (a.size() > l && a.compare(0, l, pfx) == 0) ? a.substr(l) : std::string();
        };
        auto numlist = [&](const std::string &v, std::vector<size_t> &out) {
            out.clear();
            size_t p = 0;
            while (p <= v.size()) {
                const size_t c = v.find(',', p);
                const std::string t =
                    v.substr(p, c == std::string::npos ? std::string::npos : c - p);
                if (!t.empty()) out.push_back((size_t)std::strtoul(t.c_str(), nullptr, 0));
                if (c == std::string::npos) break;
                p = c + 1;
            }
        };
        std::string v;
        if (a == "--help" || a == "-h") {
            std::printf(
                "usage: %s [options]\n"
                "  --sizes=4096,...      large-N throughput sweep (default 2^12..2^24, 8 tiers)\n"
                "  --small=4,8,16,32,64  small-N call-overhead regime (the n values)\n"
                "  --ops=SUBSTR          only entries whose name matches the substring\n"
                "  --op=EXACT            only this single entry (used by verify)\n"
                "  --gate-only           bit-level gates only (no timing/throughput); exit 4 on\n"
                "                        a failing gate\n"
                "  --list-ops            print the op registry and exit\n"
                "  --floor-out=FILE      write the A/A floor table (meaningful for SELFCHECK)\n"
                "  --trials-min=5 --trials-max=12 --cov=5.0  adaptive sampling\n"
                "  --target-ms=8         per-trial target for large N; --small-target-ms=20\n"
                "  --selfcheck           baseline vs baseline self-check (the per-size\n"
                "                        criterion uses that tier's A/A floor p90)\n"
                "  --aa-each-round       run one A/A pair after every A/B round (denominator of\n"
                "                        R4 and of the per-size criterion)\n"
                "  --tier-exclude=131072 tier left out of the per-size criterion (recorded only)\n"
                "  --no-cold / --no-pipe / --no-small\n"
                "  --label=.. --tsv-out=.. --matrix-out=.. --small-out=..\n",
                argv[0]);
            return false;
        } else if (!(v = val("--sizes=")).empty()) {
            numlist(v, A.sizes);
        } else if (!(v = val("--small=")).empty()) {
            numlist(v, A.small_ns);
        } else if (!(v = val("--ops=")).empty()) {
            A.ops_filter = v;
        } else if (!(v = val("--op=")).empty()) {
            A.op_exact = v;
        } else if (!(v = val("--floor-out=")).empty()) {
            A.floor_out = v;
        } else if (a == "--aa-each-round") {
            A.aa_each_round = true;
        } else if (a == "--intra") {
            A.intra = true;
        } else if (!(v = val("--tier-exclude=")).empty()) {
            A.tier_exclude = v;
        } else if (!(v = val("--mask-density=")).empty()) {
            A.mask_density = std::atof(v.c_str());
            if (A.mask_density < 0.0) A.mask_density = 0.0;
            if (A.mask_density > 1.0) A.mask_density = 1.0;
        } else if (a == "--gate-only") {
            A.gate_only = true;
        } else if (a == "--list-ops") {
            A.list_ops = true;
        } else if (!(v = val("--trials-min=")).empty()) {
            A.trials_min = std::atoi(v.c_str());
        } else if (!(v = val("--trials-max=")).empty()) {
            A.trials_max = std::atoi(v.c_str());
        } else if (!(v = val("--cov=")).empty()) {
            A.cov_target = std::atof(v.c_str());
        } else if (!(v = val("--target-ms=")).empty()) {
            A.target_ms = std::atof(v.c_str());
        } else if (!(v = val("--small-target-ms=")).empty()) {
            A.small_target_ms = std::atof(v.c_str());
        } else if (!(v = val("--label=")).empty()) {
            A.label = v;
        } else if (!(v = val("--tsv-out=")).empty()) {
            A.tsv_out = v;
        } else if (!(v = val("--matrix-out=")).empty()) {
            A.matrix_out = v;
        } else if (!(v = val("--small-out=")).empty()) {
            A.small_out = v;
        } else if (a == "--selfcheck") {
            A.selfcheck = true;
        } else if ((v = val("--clock=")) == "mono") {
            g_clock = CLOCK_MONOTONIC;
            g_clock_name = "wall(CLOCK_MONOTONIC)";
        } else if (v == "cpu") {
            g_clock = CLOCK_THREAD_CPUTIME_ID;
            g_clock_name = "cpu(CLOCK_THREAD_CPUTIME_ID)";
        } else if ((v = val("--order=")) == "ab") {
            A.order = 1;
        } else if (v == "ba") {
            A.order = 2;
        } else if (v == "alt") {
            A.order = 0;
        } else if (a == "--swap-out") {
            A.swap_out = true;
        } else if (a == "--no-cold") {
            A.do_cold = false;
        } else if (a == "--no-pipe") {
            A.do_pipe = false;
        } else if (a == "--no-small") {
            A.do_small = false;
        } else {
            std::fprintf(stderr, "unknown argument: %s (see --help)\n", a.c_str());
            return false;
        }
    }
    if (A.trials_min < 5)
        std::fprintf(stderr, "warning: --trials-min=%d < 5, fewer trials than the required 5\n",
                     A.trials_min);
    return true;
}

/* =====================================================================
 * 7. Main flow
 * ===================================================================== */
int main(int argc, char **argv)
{
    Args args;
    if (!parse_args(argc, argv, args)) return 2;

    /* ABI fingerprint: the fv_vec2 / fv_real layout of all three TUs (A, B, scalar anchor)
     * must match the driver's */
    {
        const unsigned long long want =
            ((unsigned long long)sizeof(fv_vec2) << 32) | (unsigned)sizeof(fv_real);
        const unsigned long long ga = ab_a_abi(), gb = ab_b_abi(), gs = ab_s_abi();
        if (ga != want || gb != want || gs != want) {
            std::printf("#####################################################################\n");
            std::printf("## hard failure: ABI fingerprint mismatch. driver %llu, A %llu, B %llu, "
                        "scalar anchor %llu\n",
                        want, ga, gb, gs);
            std::printf("## some variant TU's fv_vec2 / fv_real differs from the driver's (header "
                        "layout changed?).\n");
            std::printf("## every pointer would then be read as the wrong struct, so nothing can "
                        "be measured; exiting.\n");
            std::printf("#####################################################################\n");
            return 5;
        }
    }

    /* --list-ops: the registry is the only source of truth for "which ops exist"; verify.sh
     * and the docs both use it, so the list is never restated here (a restated list always
     * drifts away from the table). */
    if (args.list_ops) {
        std::printf("# op\tabi_op\trd\twr\treal_rd\treal_wr_a\treal_wr_b\trfo\tflags\telop"
                    "\tintra_ref\tnote\n");
        for (size_t k = 0; k < OP_COUNT; ++k) {
            const Op &op = OPS[k];
            const size_t ir = intra_ref_of(k);
            std::printf("%s\t%d\t%zu\t%zu\t%zu\t%zu\t%zu\t%s\t0x%x\t%zu\t%s\t%s\n",
                        op.name, op.abi_op, op.rd, op.wr, op.rd, real_wr_bytes(op, SIDE_A),
                        real_wr_bytes(op, SIDE_B), rfo_side_name(op, SIDE_A), op.flags,
                        op.elop, ir < OP_COUNT ? OPS[ir].name : "-", op.note);
        }
        return 0;
    }

    /* Naming an entry that does not exist in the baseline switches to INTRA automatically:
     * A/B is meaningless for it (side A has no such function at all). */
    if (!args.op_exact.empty() && !args.intra) {
        for (size_t k = 0; k < OP_COUNT; ++k)
            if (args.op_exact == OPS[k].name && has_flag(OPS[k], AB_OPF_NO_BASELINE)) {
                args.intra = true;
                break;
            }
    }

    std::vector<size_t> sel;
    for (size_t k = 0; k < OP_COUNT; ++k) {
        if (has_flag(OPS[k], AB_OPF_PIPE) && !args.do_pipe) continue;
        /* the new API (vec3/vec4) is absent from the baseline: A/B excludes it, INTRA takes
         * nothing else */
        if (args.intra != has_flag(OPS[k], AB_OPF_NO_BASELINE)) continue;
        if (!args.op_exact.empty()) {
            if (args.op_exact != OPS[k].name) continue;
        } else if (!args.ops_filter.empty() &&
                   !std::strstr(OPS[k].name, args.ops_filter.c_str()))
            continue;
        sel.push_back(k);
    }
    /* only wide SoA elements (vec3 12 B / vec4 16 B) need az/bz/aw/bw and base2; the matrix
     * transforms are 16 B AoS and need only vin + mat, so they must not drag the SoA arrays in */
    for (size_t k : sel) {
        if (OPS[k].kind == OUT_VEC4M)
            g_mat4 = true;
        else if (OPS[k].elop != 8)
            g_wide = true;
    }
    if (sel.empty()) {
        std::fprintf(stderr, "no entry matches --ops=%s%s%s (see the registry with --list-ops)\n",
                     args.ops_filter.c_str(), args.op_exact.empty() ? "" : " --op=",
                     args.op_exact.c_str());
        return 2;
    }
    if (args.intra)
        std::printf("# mode: INTRA (a control inside the candidate; **not** an improvement "
                    "over git HEAD)\n");

    /* --gate-only: bit-level gates only. verify uses it for the all-tail-lengths comparison
     * (n = N-7 .. N+3 covers the 0-3 tail of _x4 and the 0-7 tail of _soa) and it produces no
     * timing numbers at all. */
    if (args.gate_only) {
        std::printf("# gate-only: %zu ops x %zu sizes (bit-level, A/B/independent reference agree)\n",
                    sel.size(), args.sizes.size());
        size_t nfail = 0, ntot = 0;
        for (size_t n : args.sizes) {
            Buf B;
            alloc_buf(B, n, args.do_pipe); /* pipeline ops need the px/py/qx/qy scratch */
            init_data(B, n);
            for (size_t k : sel) {
                const Op &op = OPS[k];
                if (has_flag(op, AB_OPF_PIPE)) init_pipe(B); /* the previous op may have written it */
                const double in_cs = B.in_cs;
                Gate g = check_op(op, B, n, small_supported((int)k));
                bool ok = g.ok;
                if (B.in_cs != in_cs) {
                    ok = false;
                    g.why += "  the implementation wrote the input buffers (a/b/mask/idx)\n";
                }
                ++ntot;
                if (!ok) {
                    ++nfail;
                    std::printf("[GATE FAIL] %s @ n=%zu\n%s", op.name, n, g.why.c_str());
                }
            }
            free_buf(B);
        }
        std::printf("# gate-only result: %zu/%zu passed\n", ntot - nfail, ntot);
        return nfail ? 4 : 0;
    }

    double mn0, mx0, av0;
    cpu_mhz(&mn0, &mx0, &av0);
    const double load_start = loadavg1();
    const long rss0 = self_statm(2);

    std::printf("=====================================================================\n");
    std::printf("flashvec in-process A/B interleaving fixture   label=%s\n", args.label.c_str());
    std::printf("  baseline  : git HEAD snapshot bench/ab/snapshot/  rev=%s\n", AB_GIT_REV);
    std::printf("  candidate : %s\n", AB_CAND_DESC);
    std::printf("  build     : same as the top-level Makefile defaults (no -march; AVX2 by "
                "per-function dispatch)\n");
    std::printf("  entries   : %zu / %zu (registry: AB_LIB_OPS %d%s + AB_MEM_OPS %d\n"
                "              ceiling references)\n",
                sel.size(), OP_COUNT, (int)N_LIB_OPS,
                args.do_pipe ? " (with pipeline)" : " (no pipeline)", (int)N_MEM_OPS);
    std::printf("  large N   :");
    for (size_t n : args.sizes) std::printf(" 2^%.0f", std::log2((double)n));
    std::printf("\n  small N   :");
    for (size_t n : args.small_ns) std::printf(" %zu", n);
    if (!args.do_small) std::printf(" (off)");
    std::printf("\n");
    std::printf("  sampling  : at least %d and at most %d trials per point, until the\n"
                "              [paired-ratio IQR] < %.1f%%; large-N target %.1f ms/trial,\n"
                "              small-N target %.1f ms/trial\n",
                args.trials_min, args.trials_max, args.cov_target, args.target_ms,
                args.small_target_ms);
    std::printf("  cache     : L1d 32 KiB / L2 512 KiB / L3 4 MiB (single CCX, cpu0-7)\n");
    std::printf("  threads   : single-threaded; the library has no threads and no global\n"
                "              state; pinning via taskset in the Makefile\n");
    std::printf("  frequency : governor=performance, observed %.0f-%.0f MHz (mean %.0f) --\n"
                "              not pinned, so absolute values do not carry across sessions\n",
                mn0, mx0, av0);
    std::printf("  load      : loadavg(1m)=%.2f; MemAvailable=%ld MiB; VmSwap=%ld MiB; initial "
                "RSS=%ld MiB\n",
                load_start, mem_available_mib(), proc_kb("VmSwap") / 1024, rss0 * 4 / 1024);
    std::printf("  clock     : %s\n", g_clock_name);
    std::printf("  accounting: GB/s=(bytes read+written)/time; %%memcpy=this point's GB/s /\n"
                "              memcpy at the same size; ratio=median(per-round A/B)\n");
    std::printf("=====================================================================\n");

    {
        const double t0 = now_ns();
        volatile double acc = 0.0;
        while (now_ns() - t0 < 120e6) acc += 1.0;
        (void)acc;
    }

    int verdict = 0; /* 0=pass 3=self-check failed (still written to disk for review) */
    std::vector<Pt> pts;
    bool any_gate_fail = false;
    size_t n_gate = 0, n_gate_fail = 0;

    /* ================= small-N regime: call overhead ================= */
    struct SmallPt {
        const char *op;
        size_t n;
        double a_ns_call, b_ns_call, s_ns_call;
        double ratio, scalar_x;
        double a_cov, b_cov;
        int trials;
        long calls;
        double cs;
    };
    std::vector<SmallPt> small;
    if (args.do_small) {
        Buf SB;
        alloc_buf(SB, args.small_ns.back(), false);
        std::printf("\n=== small-N regime (call overhead; buffers stay in L1, each call touches "
                    "only n elements) ===\n");
        std::printf("  ns/call is the bottleneck metric of this regime; ns/elem dilutes the call "
                    "overhead.\n");
        std::printf("%-4s %-24s %10s %11s %11s %6s %10s %10s %11s %8s %s\n", "n", "op", "calls",
                    "ns/call(A)", "ns/call(B)", "ratio", "ns/el(A)", "ns/el(B)",
                    "scalar/call", "A/scalar", "checksum");
        for (size_t n : args.small_ns) {
            init_data(SB, n);
            for (size_t k : sel) {
                const Op &op = OPS[k];
                if (op.pol == POL_NOCHECK || has_flag(op, AB_OPF_PIPE)) continue;
                if (!small_supported((int)k)) continue;
                Gate g = check_op(op, SB, n, true);
                ++n_gate;
                if (!g.ok) {
                    ++n_gate_fail;
                    any_gate_fail = true;
                    std::printf("  [GATE FAIL] %s @ n=%zu\n%s", op.name, n, g.why.c_str());
                    continue;
                }
                std::vector<double> ta, tb, ts, pr;
                long calls = 1;
                {
                    (void)measure(op, SIDE_A, SB, SB.out[0], 1000); /* warm i-cache/predictors */
                    const double one = measure(op, SIDE_A, SB, SB.out[0], 1) * (double)n;
                    calls = (long)(args.small_target_ms * 1e6 / std::max(one, 1.0));
                    if (calls < 100) calls = 100;
                    if (calls > 50000000L) calls = 50000000L;
                    (void)measure(op, SIDE_A, SB, SB.out[0], calls);
                    (void)measure(op, SIDE_B, SB, SB.out[1], calls);
                    (void)measure(op, SIDE_S, SB, SB.out[0], calls);
                }
                for (int t = 0; t < args.trials_max; ++t) {
                    if (t >= args.trials_min && iqr_pct(pr) < args.cov_target) break;
                    double x, y, z;
                    if (t & 1) {
                        y = measure(op, SIDE_B, SB, SB.out[1], calls);
                        x = measure(op, SIDE_A, SB, SB.out[0], calls);
                        z = measure(op, SIDE_S, SB, SB.out[0], calls);
                    } else {
                        x = measure(op, SIDE_A, SB, SB.out[0], calls);
                        y = measure(op, SIDE_B, SB, SB.out[1], calls);
                        z = measure(op, SIDE_S, SB, SB.out[0], calls);
                    }
                    ta.push_back(x);
                    tb.push_back(y);
                    ts.push_back(z);
                    pr.push_back(x / y);
                }
                SmallPt sp;
                sp.op = op.name;
                sp.n = n;
                sp.a_ns_call = median(ta) * (double)n;
                sp.b_ns_call = median(tb) * (double)n;
                sp.s_ns_call = median(ts) * (double)n;
                sp.ratio = median(pr);
                sp.scalar_x = sp.s_ns_call / sp.a_ns_call;
                sp.a_cov = iqr_pct(pr); /* in small N this field is the paired IQR% */
                sp.b_cov = cov_pct(ta);
                sp.trials = (int)ta.size();
                sp.calls = calls;
                sp.cs = g.cs_ref;
                small.push_back(sp);
                std::printf("%-4zu %-24s %10ld %11.3f %11.3f %6.3f %10.4f %10.4f %11.3f %7.2fx %s\n",
                            n, op.name, calls, sp.a_ns_call, sp.b_ns_call, sp.ratio,
                            sp.a_ns_call / (double)n, sp.b_ns_call / (double)n, sp.s_ns_call,
                            sp.scalar_x, sp.a_cov > args.cov_target ? "CoV!" : "");
                std::fflush(stdout);
            }
        }
        free_buf(SB);
        std::printf("\n");
    }

    /* ================= large-N throughput sweep ================= */
    for (size_t n : args.sizes) {
        const long rss_before = self_statm(2);
        const long mf_before = self_minflt();
        Buf B;
        alloc_buf(B, n, args.do_pipe);
        init_data(B, n);
        const long rss_after = self_statm(2);
        const size_t rss_grew =
            (rss_after > rss_before) ? (size_t)(rss_after - rss_before) * 4096 : 0;
        const bool rss_bad = rss_grew * 10 < B.alloc_bytes * 9;

        std::printf("\n=== size n=%zu  allocated %.1f MiB  RSS +%.1f MiB%s  loadavg=%.2f ===\n", n,
                    (double)B.alloc_bytes / 1048576.0, (double)rss_grew / 1048576.0,
                    rss_bad ? "  [!! RSS below expectation: pages may be swapped, "
                              "throughput unreliable]" : "",
                    loadavg1());

        std::vector<Gate> gates(OP_COUNT);
        for (size_t k : sel) {
            const double in_cs = B.in_cs;
            gates[k] = check_op(OPS[k], B, n, false);
            ++n_gate;
            if (B.in_cs != in_cs) {
                gates[k].ok = false;
                gates[k].why += "  the implementation wrote the input buffers (a/b/mask/idx)\n";
            }
            if (!gates[k].ok) {
                ++n_gate_fail;
                any_gate_fail = true;
                std::printf("  [GATE FAIL] %s @ n=%zu\n%s", OPS[k].name, n, gates[k].why.c_str());
            }
        }
        if (any_gate_fail) {
            std::printf("\n");
            std::printf("#####################################################################\n");
            std::printf("## hard failure: the checksum gate failed (%zu/%zu items).\n", n_gate_fail,
                        n_gate);
            std::printf("## every item above names the function: baseline / candidate / independent\n");
            std::printf("## scalar reference do not compute the same thing there. Any ratio and any\n");
            std::printf("## throughput would be meaningless, so this run prints no timing or\n## throughput results.\n");
            std::printf("#####################################################################\n");
            free_buf(B);
            return 4;
        }

        /* one measurement point: adaptive sampling + interleaving */
        auto run_point = [&](size_t opi, bool cold) -> Pt {
            const Op &op = OPS[opi];
            Pt p{};
            p.op = op.name;
            p.opi = opi;
            p.n = n;
            p.cold = cold;
            p.mode = "ab";
            p.ref_op = "-";
            p.rd = op.rd;
            p.wr = op.wr;
            p.cs_a = gates[opi].cs_a;
            p.cs_b = gates[opi].cs_b;
            p.cs_ref = gates[opi].cs_ref;
            const double l0 = loadavg1();
            std::vector<double> ta, tb, pr, aa;
            long iters = 1;
            int tmax = args.trials_max;
            if (!cold) {
                /* Warm up before calibrating: the first call runs with a cold i-cache and cold
                 * predictors, and using it as the calibration base underestimates the inner
                 * repetition count by an order of magnitude (been there). */
                (void)measure(op, SIDE_A, B, B.out[0], 4);
                /* measure returns ns/element of this op; multiply by the element count to get
                 * the duration of one call */
                const double one = measure(op, SIDE_A, B, B.out[0], 1) * (double)nel(op, n);
                iters = (long)(args.target_ms * 1e6 / std::max(one, 1.0));
                if (iters < 1) iters = 1;
                if (iters > 200000) iters = 200000;
                /* points whose single call already exceeds 30 ms (the DRAM tiers) do not chase
                 * 12 trials: that is pure burn */
                if (one > 30e6) tmax = std::min(args.trials_max, 7);
                (void)measure(op, SIDE_A, B, B.out[0], iters);
                (void)measure(op, SIDE_B, B, B.out[1], iters);
            }
            for (int t = 0; t < tmax; ++t) {
                /* The stop condition looks at the interquartile range of the paired ratios, not
                 * at the sample CoV: the ratio is what gets judged, and a sample std blows up on
                 * one stall -- which would keep burning time long after the ratio settled. */
                if (t >= args.trials_min && iqr_pct(pr) < args.cov_target) break;
                double x, y;
                OutSet &OA = args.swap_out ? B.out[1] : B.out[0];
                OutSet &OB = args.swap_out ? B.out[0] : B.out[1];
                const bool b_first = (args.order == 2) || (args.order == 0 && (t & 1));
                if (b_first) {
                    if (cold) drop_pages(B);
                    y = measure(op, SIDE_B, B, OB, cold ? 1 : iters);
                    if (cold) drop_pages(B);
                    x = measure(op, SIDE_A, B, OA, cold ? 1 : iters);
                } else {
                    if (cold) drop_pages(B);
                    x = measure(op, SIDE_A, B, OA, cold ? 1 : iters);
                    if (cold) drop_pages(B);
                    y = measure(op, SIDE_B, B, OB, cold ? 1 : iters);
                }
                ta.push_back(x);
                tb.push_back(y);
                pr.push_back(x / y);
                if (args.aa_each_round) {
                    /* One more A/A pair inside the same round, the same code and the same data:
                     * this is the R4 denominator. Only while it shares the A/B immersion state
                     * and clock-drift segment is it a same-distribution floor. */
                    const double u = measure(op, SIDE_A, B, OA, cold ? 1 : iters);
                    const double v = measure(op, SIDE_A, B, OB, cold ? 1 : iters);
                    aa.push_back(u / v);
                }
            }
            p.a_med = median(ta);
            p.b_med = median(tb);
            p.ratio = median(pr); /* main criterion: pairing cancels the drift */
            p.ratio_m = p.a_med / p.b_med;
            p.a_min = *std::min_element(ta.begin(), ta.end());
            p.a_max = *std::max_element(ta.begin(), ta.end());
            p.b_min = *std::min_element(tb.begin(), tb.end());
            p.b_max = *std::max_element(tb.begin(), tb.end());
            p.pr_min = *std::min_element(pr.begin(), pr.end());
            p.pr_max = *std::max_element(pr.begin(), pr.end());
            p.a_cov = cov_pct(ta);
            p.b_cov = cov_pct(tb);
            p.trials = (int)ta.size();
            p.iters = iters;
            p.gbs_a = (double)(op.rd + op.wr) / p.a_med;
            p.gbs_b = (double)(op.rd + op.wr) / p.b_med;
            p.real_rd = op.rd;
            p.real_wr_a = real_wr_bytes(op, SIDE_A);
            p.real_wr_b = real_wr_bytes(op, SIDE_B);
            p.rfo = op.rfo;
            p.real_gbs_a = (double)(op.rd + p.real_wr_a) / p.a_med;
            p.real_gbs_b = (double)(op.rd + p.real_wr_b) / p.b_med;
            p.aa_n = (int)aa.size();
            p.aa_p50 = aa.empty() ? 1.0 : median(aa);
            p.aa_p90 = aa.empty() ? 1.0 : quantile(aa, 0.90);
            p.aa_max = aa.empty() ? 1.0 : quantile(aa, 1.0);
            p.aa_dev_max = std::max(std::fabs(p.aa_max - 1.0), std::fabs(quantile(aa, 0.0) - 1.0));
            p.load0 = l0;
            p.load1 = loadavg1();
            p.riqr = iqr_pct(pr);
            p.flags = (p.riqr > args.cov_target) ? AB_F_COV : 0u;
            if (p.load0 > 6.0 || p.load1 > 6.0) p.flags |= AB_F_LOAD;
            if (rss_bad) p.flags |= AB_F_RSS;
            return p;
        };


        /* ---- INTRA: a control inside the candidate (new entry vs an existing entry that
         * moves the same bytes). Both ops run on **side B** (the candidate).
         * ratio = t_ref / t_new, > 1 means the new entry is faster. This is not "relative to
         * the baseline", because the baseline has none of these entries. */
        auto run_pair = [&](size_t inew, size_t iref) -> Pt {
            const Op &nw = OPS[inew], &rf = OPS[iref];
            Pt p{};
            p.op = nw.name;
            p.opi = inew;
            p.n = n;
            p.cold = false;
            p.mode = "intra";
            p.ref_op = rf.name;
            p.rd = nw.rd;
            p.wr = nw.wr;
            p.real_rd = nw.rd;
            p.real_wr_a = real_wr_bytes(nw, SIDE_B);
            p.real_wr_b = real_wr_bytes(rf, SIDE_B);
            p.rfo = nw.rfo;
            p.cs_a = gates[iref].cs_b;
            p.cs_b = gates[inew].cs_b;
            p.cs_ref = gates[inew].cs_ref;
            const double l0 = loadavg1();
            std::vector<double> tn, tr, pr, aa;
            int tmax = args.trials_max;
            (void)measure(nw, SIDE_B, B, B.out[1], 4);
            const double one = measure(nw, SIDE_B, B, B.out[1], 1) * (double)nel(nw, n);
            long iters = (long)(args.target_ms * 1e6 / std::max(one, 1.0));
            if (iters < 1) iters = 1;
            if (iters > 200000) iters = 200000;
            if (one > 30e6) tmax = std::min(args.trials_max, 7);
            (void)measure(nw, SIDE_B, B, B.out[1], iters);
            (void)measure(rf, SIDE_B, B, B.out[0], iters);
            for (int t = 0; t < tmax; ++t) {
                if (t >= args.trials_min && iqr_pct(pr) < args.cov_target) break;
                double x, y;
                if (t & 1) {
                    y = measure(nw, SIDE_B, B, B.out[1], iters);
                    x = measure(rf, SIDE_B, B, B.out[0], iters);
                } else {
                    x = measure(rf, SIDE_B, B, B.out[0], iters);
                    y = measure(nw, SIDE_B, B, B.out[1], iters);
                }
                tr.push_back(x);
                tn.push_back(y);
                /* Normalize by bytes: the two ops have different element widths, so comparing
                 * x/y directly is per element rather than per byte. Total time per call =
                 * ns/element * element count, and both move the same byte volume, so this is
                 * the honest "who is faster for the same work". */
                pr.push_back((x * (double)nel(rf, n)) / (y * (double)nel(nw, n)));
                if (args.aa_each_round) {
                    /* The A/A floor is taken with the **reference op**: it represents this
                     * tier's "same code twice" noise */
                    const double u = measure(rf, SIDE_B, B, B.out[0], iters);
                    const double v = measure(rf, SIDE_B, B, B.out[1], iters);
                    aa.push_back(u / v);
                }
            }
            p.a_med = median(tr);
            p.b_med = median(tn);
            p.ratio = median(pr);
            p.ratio_m = p.a_med / p.b_med;
            p.a_cov = cov_pct(tr);
            p.b_cov = cov_pct(tn);
            p.riqr = iqr_pct(pr);
            p.trials = (int)tr.size();
            p.iters = iters;
            p.gbs_a = (double)(rf.rd + rf.wr) / p.a_med;
            p.gbs_b = (double)(nw.rd + nw.wr) / p.b_med;
            p.real_gbs_a = (double)(rf.rd + p.real_wr_b) / p.a_med;
            p.real_gbs_b = (double)(nw.rd + p.real_wr_a) / p.b_med;
            p.pct_memcpy_a = p.pct_memcpy_b = 0.0;
            p.real_pct_memcpy_a = p.real_pct_memcpy_b = 0.0;
            p.load0 = l0;
            p.load1 = loadavg1();
            p.flags = (p.riqr > args.cov_target) ? AB_F_COV : 0u;
            if (p.load0 > 6.0 || p.load1 > 6.0) p.flags |= AB_F_LOAD;
            p.aa_n = (int)aa.size();
            p.aa_p50 = aa.empty() ? 1.0 : median(aa);
            p.aa_p90 = aa.empty() ? 1.0 : quantile(aa, 0.90);
            p.aa_max = aa.empty() ? 1.0 : quantile(aa, 1.0);
            p.aa_dev_max = std::max(std::fabs(p.aa_max - 1.0), std::fabs(quantile(aa, 0.0) - 1.0));
            return p;
        };

        /* ---- INTRA report (explicitly labelled: a control inside the candidate, not an
         * improvement over HEAD) ---- */
        if (args.intra) {
            std::printf("\n--- INTRA control: inside the candidate (new entry vs an existing entry "
                        "moving the same bytes) ---\n");
            std::printf("--- **this is not an improvement over git HEAD**: these entries do not exist "
                        "in HEAD ---\n");
            std::printf("--- both ops move the same bytes in one tier; ratio = t(ref)/t(new), >1 means "
                        "the new entry is faster ---\n");
            std::printf("%-24s %-22s %7s %7s %10s %10s %6s %9s %9s %-14s %s\n", "op(new)", "ref",
                        "n(new)", "n(ref)", "ns/e(new)", "ns/e(ref)", "ratio", "GB/s(new)",
                        "GB/s(ref)", "AA p50/p90/max", "flags");
            for (size_t k : sel) {
                const Op &nw = OPS[k];
                if (!has_flag(nw, AB_OPF_NO_BASELINE)) continue;
                const size_t iref = intra_ref_of(k);
                if (iref >= OP_COUNT) {
                    std::printf("%-24s %-22s (no INTRA pair in the registry, skipped)\n", nw.name,
                                "-");
                    continue;
                }
                Pt pp = run_pair(k, iref);
                pts.push_back(pp);
                char aas[40];
                std::snprintf(aas, sizeof aas, "%.3f/%.3f/%.3f", pp.aa_p50, pp.aa_p90, pp.aa_max);
                const char *fl = (pp.flags == 0) ? "" : ((pp.flags & AB_F_LOAD) ? "LOAD" : "CoV!");
                std::printf("%-24s %-22s %7zu %7zu %10.4f %10.4f %6.3f %9.2f %9.2f %-14s %s\n",
                            nw.name, OPS[iref].name, nel(nw, n), nel(OPS[iref], n), pp.b_med,
                            pp.a_med, pp.ratio, pp.gbs_b, pp.gbs_a, aas, fl);
                std::fflush(stdout);
            }
            free_buf(B);
            continue;
        }

        /* Ceiling references (the %memcpy denominators; memset/pure-read/triad use the same
         * formula, so two denominators suffice). The real_* figures divide by memcpy's **real**
         * traffic -- at scale memcpy itself is probably an NT store (libc?), so the absolute
         * value of that column is not fully trustworthy, while comparisons within a tier are. */
        double memcpy_warm = 0.0, memcpy_cold = 0.0;
        double memcpy_real_warm = 0.0, memcpy_real_cold = 0.0;
        auto apply_pct = [](Pt &p, double mem, double mem_real) {
            p.pct_memcpy_a = mem > 0.0 ? p.gbs_a / mem * 100.0 : 0.0;
            p.pct_memcpy_b = mem > 0.0 ? p.gbs_b / mem * 100.0 : 0.0;
            p.real_pct_memcpy_a = mem_real > 0.0 ? p.real_gbs_a / mem_real * 100.0 : 0.0;
            p.real_pct_memcpy_b = mem_real > 0.0 ? p.real_gbs_b / mem_real * 100.0 : 0.0;
        };
        {
            Pt p = run_point(memref_index(AB_MF_COPY), false);
            memcpy_warm = p.gbs_a;
            memcpy_real_warm = p.real_gbs_a;
            p.pct_memcpy_a = p.pct_memcpy_b = 100.0;
            p.real_pct_memcpy_a = p.real_pct_memcpy_b = 100.0;
            pts.push_back(p);
        }
        if (args.do_cold) {
            Pt p = run_point(memref_index(AB_MF_COPY), true);
            memcpy_cold = p.gbs_a;
            memcpy_real_cold = p.real_gbs_a;
            p.pct_memcpy_a = p.pct_memcpy_b = 100.0;
            p.real_pct_memcpy_a = p.real_pct_memcpy_b = 100.0;
            pts.push_back(p);
        }
        for (const size_t opi : {memref_index(AB_MF_SET), memref_index(AB_MF_READ),
                                memref_index(AB_MF_TRIAD)}) {
            if (opi >= OP_COUNT) continue;
            if (!args.ops_filter.empty() && !std::strstr(OPS[opi].name, args.ops_filter.c_str()))
                continue;
            Pt w = run_point(opi, false);
            apply_pct(w, memcpy_warm, memcpy_real_warm);
            pts.push_back(w);
        }

        /* Run every warm point first (the data is still intact), then cold (cold empties the
         * pages, so it must not pollute a later warm run) */
        std::printf("%-24s %5s %7s %5s %8s %9s %6s %10s %9s %9s %6s %-15s %3s %-9s %-14s "
                    "%-22s %s\n",
                    "op", "B/e", "realB/e", "tier", "GB/s", "realGB/s", "%mem", "CoV/IQR%",
                    "A_med", "B_med", "ratio", "pair[min,max]", "tr", "rfo",
                    "AA p50/p90/max", "checksum", "flags");
        for (size_t k : sel) {
            const Op &op = OPS[k];
            if (op.pol == POL_NOCHECK) continue;
            if (has_flag(op, AB_OPF_PIPE)) init_pipe(B);
            const size_t ws = (op.rd + op.wr) * nel(op, n);
            Pt pw = run_point(k, false);
            apply_pct(pw, memcpy_warm, memcpy_real_warm);
            /* After many timed calls an in-place op's output state has run away (it keeps
             * accumulating into the same buffer), so A/B is not compared after timing; the
             * bit-level verdict comes from check_op before timing. */
            if (!has_flag(op, AB_OPF_INPLACE)) {
                std::string m =
                    /* Element count, not the sweep size: the output array holds nel() elements,
                     * and for a 16 B element kind 4 * n floats would read twice the buffer and
                     * report every op as broken. */
                    cmp_out(op, B.out[0], B.out[1], nel(op, n));
                if (!m.empty()) {
                    std::printf("  [GATE FAIL] %s @ n=%zu: A/B outputs differ after timing: %s\n",
                                op.name, n, m.c_str());
                    any_gate_fail = true;
                }
            }
            pts.push_back(pw);
            const char *fl = (pw.flags == 0)        ? ""
                             : (pw.flags & AB_F_RSS) ? "RSS?"
                             : (pw.flags & AB_F_LOAD) ? "LOAD"
                                                      : "CoV!";
            char pair[40], rfo[24], aas[40];
            std::snprintf(pair, sizeof pair, "[%.3f,%.3f]", pw.pr_min, pw.pr_max);
            if (pw.aa_n)
                std::snprintf(aas, sizeof aas, "%.3f/%.3f/%.3f", pw.aa_p50, pw.aa_p90,
                              pw.aa_max);
            else
                std::snprintf(aas, sizeof aas, "-");
            const char *ra = rfo_side_name(op, SIDE_A), *rb = rfo_side_name(op, SIDE_B);
            if (std::strcmp(ra, rb) == 0)
                std::snprintf(rfo, sizeof rfo, "%s", ra);
            else
                std::snprintf(rfo, sizeof rfo, "%s/%s", ra, rb);
            std::printf("%-24s %5zu %7zu %5s %8.2f %9.2f %6.1f %5.1f/%-4.1f "
                        "%9.4f %9.4f %6.3f %-15s %3d %-9s %-14s %-22s %s\n",
                        op.name, op.rd + op.wr, pw.real_rd + pw.real_wr_a, tier_of(ws),
                        pw.gbs_a, pw.real_gbs_a, pw.pct_memcpy_a, pw.a_cov, pw.riqr, pw.a_med,
                        pw.b_med, pw.ratio, pair, pw.trials, rfo, aas,
                        checksum_text(op, pw.cs_a, pw.cs_b, pw.cs_ref).c_str(), fl);
            std::fflush(stdout);
        }
        std::printf("  -- ceiling references (same tier, same clock; A/B run identical code, so "
                    "the ratio is noise) --\n");
        std::printf("  %-16s %5s %7s %5s %8s %9s  %s\n", "op", "B/e", "realB/e", "tier", "GB/s",
                    "realGB/s", "ratio  CoV");
        for (const Pt &p : pts) {
            if (p.n != n || p.cold) continue;
            if (OPS[p.opi].pol != POL_NOCHECK) continue;
            std::printf("  %-16s %5zu %7zu %5s %8.2f %9.2f  ratio=%.3f  CoV=%.1f%%  rfo=%s\n",
                        p.op, p.rd + p.wr, p.real_rd + p.real_wr_a,
                        tier_of((p.rd + p.wr) * n), p.gbs_a, p.real_gbs_a, p.ratio, p.a_cov,
                        rfo_side_name(OPS[p.opi], SIDE_A));
        }
        std::printf("  -- cold (one call after madvise(DONTNEED), page faults included) --\n");
        if (args.do_cold) {
            std::printf("  %-24s %8s %9s %6s %-15s\n", "op", "coldGB/s", "coldRGB/s",
                        "%memC", "pair[min,max]");
            for (size_t k : sel) {
                const Op &op = OPS[k];
                if (op.pol == POL_NOCHECK) continue;
                init_pipe(B);
                Pt pc = run_point(k, true);
                apply_pct(pc, memcpy_cold, memcpy_real_cold);
                pts.push_back(pc);
                std::printf("  %-24s %8.2f %9.2f %6.1f [%.3f,%.3f]\n", op.name, pc.gbs_a,
                            pc.real_gbs_a, pc.pct_memcpy_a, pc.pr_min, pc.pr_max);
            }
            for (const size_t opi : {memref_index(AB_MF_SET), memref_index(AB_MF_READ),
                                    memref_index(AB_MF_TRIAD)}) {
                if (opi >= OP_COUNT) continue;
                if (!args.ops_filter.empty() &&
                    !std::strstr(OPS[opi].name, args.ops_filter.c_str()))
                    continue;
                Pt c = run_point(opi, true);
                apply_pct(c, memcpy_cold, memcpy_real_cold);
                pts.push_back(c);
            }
        }
        const long mf = self_minflt();
        std::printf("  (this tier: %ld minor faults; freed %.1f MiB, VmSwap=%ld MiB)\n",
                    (mf > mf_before) ? (mf - mf_before) : -1,
                    (double)B.alloc_bytes / 1048576.0, proc_kb("VmSwap") / 1024);
        free_buf(B);
        if (any_gate_fail) return 4;
    }

    /* ================= summary ================= */
    std::printf("\n=== summary ===\n");
    std::printf("checksum gates: %zu/%zu items passed (A, B and the independent reference\n"
                "agree; dot_soa uses max(%.0e*|r|, %.0e*sum|term|))\n",
                n_gate - n_gate_fail, n_gate, (double)AB_REL_FRAC, (double)AB_ABS_FRAC);
    std::printf("large-N points: %zu (warm+cold)  small-N points: %zu  load: loadavg(1m) %.2f -> "
                "%.2f\n",
                pts.size(), small.size(), load_start, loadavg1());

    auto find_pt = [&](size_t opi, size_t n, bool cold) -> const Pt * {
        for (const Pt &p : pts)
            if (p.op == OPS[opi].name && p.n == n && p.cold == cold) return &p;
        return nullptr;
    };

    for (int cold = 0; cold <= (args.do_cold ? 1 : 0); ++cold) {
        std::printf("\n--- throughput matrix: %s, GB/s (single-threaded; rd+wr both counted) ---\n",
                    cold ? "cold" : "warm");
        std::printf("%-24s", "op");
        for (size_t n : args.sizes) std::printf(" %9.0f", (double)n);
        std::printf("\n");
        for (const size_t opi : sel) {
            if (!find_pt(opi, args.sizes[0], cold != 0)) continue;
            std::printf("%-24s", OPS[opi].name);
            for (size_t n : args.sizes) {
                const Pt *p = find_pt(opi, n, cold != 0);
                if (p)
                    std::printf(" %9.2f", p->gbs_a);
                else
                    std::printf(" %9s", "-");
            }
            std::printf("\n");
        }
        std::printf("\n--- percent of same-size memcpy throughput: %s ---\n", cold ? "cold" : "warm");
        std::printf("%-24s", "op");
        for (size_t n : args.sizes) std::printf(" %9.0f", (double)n);
        std::printf("\n");
        for (const size_t opi : sel) {
            if (!find_pt(opi, args.sizes[0], cold != 0)) continue;
            std::printf("%-24s", OPS[opi].name);
            for (size_t n : args.sizes) {
                const Pt *p = find_pt(opi, n, cold != 0);
                if (p)
                    std::printf(" %8.1f%%", p->pct_memcpy_a);
                else
                    std::printf(" %9s", "-");
            }
            std::printf("\n");
        }
    }

    std::printf("\n--- candidate/baseline ratio matrix (warm; >1 = candidate is faster) ---\n");
    std::printf("%-24s", "op");
    for (size_t n : args.sizes) std::printf(" %9.0f", (double)n);
    std::printf("\n");
    for (const size_t opi : sel) {
        if (!find_pt(opi, args.sizes[0], false)) continue;
        std::printf("%-24s", OPS[opi].name);
        for (size_t n : args.sizes) {
            const Pt *p = find_pt(opi, n, false);
            if (p)
                std::printf(" %9.3f", p->ratio);
            else
                std::printf(" %9s", "-");
        }
        std::printf("\n");
    }

    std::printf("\n--- real memory-traffic matrix (warm, GB/s; after RFO accounting: plain store "
                "writes 2x, NT store 1x) ---\n");
    std::printf("%-24s", "op");
    for (size_t n : args.sizes) std::printf(" %9.0f", (double)n);
    std::printf("\n");
    for (const size_t opi : sel) {
        if (!find_pt(opi, args.sizes[0], false)) continue;
        std::printf("%-24s", OPS[opi].name);
        for (size_t n : args.sizes) {
            const Pt *p = find_pt(opi, n, false);
            if (p)
                std::printf(" %9.2f", p->real_gbs_a);
            else
                std::printf(" %9s", "-");
        }
        std::printf("\n");
    }

    std::printf("\n--- RFO accounting: nominal vs real memory traffic (bytes per element; "
                "real/nom = real/nominal) ---\n");
    std::printf("%-24s %4s %4s %8s %8s %9s  %-9s %s\n", "op", "rd", "wr", "realA", "realB",
                "real/nom", "rfo", "basis");
    for (const size_t opi : sel) {
        const Op &op = OPS[opi];
        const size_t wa = real_wr_bytes(op, SIDE_A), wb = real_wr_bytes(op, SIDE_B);
        const double nom = (double)(op.rd + op.wr);
        const double fac = nom > 0.0 ? (double)(op.rd + wa) / nom : 0.0;
        char rf[24];
        const char *ra = rfo_side_name(op, SIDE_A), *rb = rfo_side_name(op, SIDE_B);
        if (std::strcmp(ra, rb) == 0)
            std::snprintf(rf, sizeof rf, "%s", ra);
        else
            std::snprintf(rf, sizeof rf, "%s/%s", ra, rb);
        std::printf("%-24s %4zu %4zu %8zu %8zu %8.2fx  %-9s %s\n", op.name, op.rd, op.wr,
                    op.rd + wa, op.rd + wb, fac, rf, op.note);
    }
    std::printf("  (real/nom quantifies the real payoff of an NT store: the same add op with a "
                "nominal 24 B/e\n");
    std::printf("   a plain store moves 32 B/e real, an NT store 24 B/e real -- the difference\n"
                "   is traffic, not speed;\n");
    std::printf("   different A/B columns = one side degraded to a plain store, i.e. the baseline "
                "header has no _stream yet.)\n");

    std::printf("\n--- throughput knee (warm): the size where throughput drops to 50%% of its own "
                "peak ---\n");
    for (const size_t opi : sel) {
        double peak = 0.0;
        size_t peak_n = 0;
        for (size_t n : args.sizes) {
            const Pt *p = find_pt(opi, n, false);
            if (p && p->gbs_a > peak) {
                peak = p->gbs_a;
                peak_n = n;
            }
        }
        if (peak <= 0.0) continue;
        std::printf("  %-24s peak %6.2f GB/s @ n=%-8zu", OPS[opi].name, peak, peak_n);
        for (size_t n : args.sizes) {
            const Pt *p = find_pt(opi, n, false);
            if (p && n > peak_n && p->gbs_a < peak * 0.5) {
                std::printf("  halves at n=%zu (%.2f GB/s)", n, p->gbs_a);
                break;
            }
        }
        std::printf("\n");
    }

    if (!small.empty()) {
        std::printf("\n--- small-N regime: ns/call (batch A) and the factor against the scalar "
                    "anchor ---\n");
        std::printf("%-24s", "op");
        for (size_t n : args.small_ns) std::printf(" %10zu", n);
        std::printf("\n");
        for (size_t k : sel) {
            if (!small_supported((int)k)) continue;
            bool any = false;
            for (const SmallPt &s : small)
                if (s.op == OPS[k].name) any = true;
            if (!any) continue;
            std::printf("%-24s", OPS[k].name);
            for (size_t n : args.small_ns) {
                const SmallPt *f = nullptr;
                for (const SmallPt &s : small)
                    if (s.op == OPS[k].name && s.n == n) f = &s;
                if (f)
                    std::printf(" %10.2f", f->a_ns_call);
                else
                    std::printf(" %10s", "-");
            }
            std::printf("\n");
        }
        std::printf("\n--- small-N regime: speedup A / scalar_loop (>1 = the batch entry wins) ---\n");
        std::printf("%-24s", "op");
        for (size_t n : args.small_ns) std::printf(" %10zu", n);
        std::printf("\n");
        for (size_t k : sel) {
            if (!small_supported((int)k)) continue;
            std::printf("%-24s", OPS[k].name);
            for (size_t n : args.small_ns) {
                const SmallPt *f = nullptr;
                for (const SmallPt &s : small)
                    if (s.op == OPS[k].name && s.n == n) f = &s;
                if (f)
                    std::printf(" %9.2fx", f->scalar_x);
                else
                    std::printf(" %10s", "-");
            }
            std::printf("\n");
        }
    }

    {
        size_t tot = 0, ok = 0, cold_bad = 0;
        std::vector<int> tr;
        for (const Pt &p : pts) {
            ++tot;
            if (!(p.flags & AB_F_COV)) ++ok;
            if (p.cold && (p.flags & AB_F_COV)) ++cold_bad;
            tr.push_back(p.trials);
        }
        std::sort(tr.begin(), tr.end());
        std::printf("\nsampling stability: large-N points with paired IQR < %.1f%%\n"
                    "(= resolvable): %zu/%zu (%zu cold below target); trials: median %d max %d\n",
                    args.cov_target, ok, tot, cold_bad, tr[tr.size() / 2], tr.back());
    }

    if (args.selfcheck) {
        /* Two classes of points:
         *   resolvable   = both CoV <= threshold -> the ratio must land in [0.85,1.15], else the
         *                  fixture fails;
         *   unresolvable = either CoV > threshold (external load burned a hole through this
         *                  tier) -> numbers are still reported but never judged.
         * One more guard: resolvable points must be the majority, otherwise the whole machine is
         * unfit for measurement and the fixture fails as well. */
        std::vector<double> wr, cr;
        size_t bad_w = 0, bad_c = 0, unres_w = 0, unres_c = 0, resolvable_w = 0;
        for (const Pt &p : pts) {
            if (OPS[p.opi].pol == POL_NOCHECK) continue;
            const bool unres = (p.flags & AB_F_COV) != 0;
            if (p.cold) {
                cr.push_back(p.ratio);
                if (unres) ++unres_c;
                else if (p.ratio < AB_BAND_LO || p.ratio > AB_BAND_HI) ++bad_c;
            } else {
                wr.push_back(p.ratio);
                if (unres) ++unres_w;
                else {
                    ++resolvable_w;
                    if (p.ratio < AB_BAND_LO || p.ratio > AB_BAND_HI) ++bad_w;
                }
            }
        }
        std::printf("\nself-check (baseline vs baseline; ratios must fall in [%.2f, %.2f]):\n",
                    AB_BAND_LO,
                    AB_BAND_HI);
        if (!wr.empty()) {
            std::sort(wr.begin(), wr.end());
            std::printf("  warm: %zu pts  min=%.4f  p50=%.4f  max=%.4f  resolvable %zu  out of "
                        "band %zu  unresolvable %zu\n",
                        wr.size(), wr.front(), wr[wr.size() / 2], wr.back(), resolvable_w, bad_w,
                        unres_w);
        }
        if (!cr.empty()) {
            std::sort(cr.begin(), cr.end());
            std::printf("  cold: %zu pts  min=%.4f  p50=%.4f  max=%.4f  unresolvable %zu  out of "
                        "band %zu\n",
                        cr.size(), cr.front(), cr[cr.size() / 2], cr.back(), unres_c, bad_c);
        }
        std::printf("  (unresolvable = the sample CoV on either side > %.1f%%, external load burned "
                    "through\n   that tier; these ratios are reported but **never judged**, their "
                    "GB/s magnitude still is)\n", args.cov_target);
        std::vector<double> all = wr;
        all.insert(all.end(), cr.begin(), cr.end());
        std::sort(all.begin(), all.end());
        const double edges[] = {0.90, 0.95, 0.98, 1.02, 1.05, 1.10};
        const char *lbl[] = {"<0.90     ", "[0.90,0.95)", "[0.95,0.98)", "[0.98,1.02)",
                             "[1.02,1.05)", "[1.05,1.10)", ">=1.10    "};
        size_t cnt[7] = {0, 0, 0, 0, 0, 0, 0};
        for (double v : all) {
            int b = 0;
            while (b < 6 && v >= edges[b]) ++b;
            ++cnt[b];
        }
        std::printf("  ratio distribution (large N, warm+cold, %zu pts):\n", all.size());
        for (int b = 0; b < 7; ++b) {
            std::printf("    %s %4zu  ", lbl[b], cnt[b]);
            for (size_t j = 0; j < cnt[b]; ++j) std::printf("#");
            std::printf("\n");
        }
        size_t bad_s = 0, unres_s = 0;
        for (const SmallPt &s : small) {
            if (s.a_cov > args.cov_target || s.b_cov > args.cov_target) {
                ++unres_s;
            } else if (s.ratio < AB_BAND_LO || s.ratio > AB_BAND_HI) {
                ++bad_s;
            }
        }
        /* All-point unbiasedness: over 240 paired measurements the median has to hug 1.0.
         * A few points burned through by noise do not move it (the median barely looks at the
         * tails), which makes this harder evidence of unbiasedness than "every point in band". */
        std::vector<double> allr = wr;
        allr.insert(allr.end(), cr.begin(), cr.end());
        for (const SmallPt &s : small) allr.push_back(s.ratio);
        const double gm = allr.empty() ? 1.0 : median(allr);
        std::printf("  all-point median ratio: %.4f (%zu paired measurements, criterion "
                    "[0.95, 1.05])%s\n",
                    gm, allr.size(), (gm >= 0.95 && gm <= 1.05) ? "  OK" : "  <== outside");
        /* Environment floor: if too few points are resolvable to mean anything, the run's
         * conclusion does not hold. 5% is the "at least a few dozen points can still be
         * judged" floor; with 30-50% external load this machine usually resolves only 5-15%. */
        /* Bias per tier: the median ratio of the 15 ops in one tier is a robust estimate of
         * "are the two sides equivalent in this tier" (median of 15 samples, standard error
         * about 1/4 of the single-point noise). Per-point extremes are still reported, but a
         * fixture-level fault such as "one tier drifted as a whole" is always caught here. */
        /* The per-size criterion is **floor-based**: a tier's median ratio must sit inside the
         * p90 of the A/A floor measured for that tier itself, not inside a hard-coded
         * [0.90, 1.10]. A floor measured in the same runs beats a constant -- but it is still a
         * finite-sample estimate and can come out far too small (0.3% was measured), so the
         * floor actually applied is max(measured, AB_MIN_TIER_FLOOR).
         * Fallback: floor > 20% means the tier cannot be judged and is tagged INCONCLUSIVE
         * (neither PASS nor FAIL). Blacklisted tiers (default 131072, the L3 boundary where A/A
         * runs 0.50-2.22) are recorded and never judged. */
        std::printf("  per-size median ratio (warm / cold) -- criterion = that tier's\n"
                    "  measured A/A floor p90 (no longer a hard-coded [0.90,1.10]):\n");
        bool size_bad = false;
        size_t tier_inconc = 0, tier_judged = 0;
        std::vector<std::string> tier_fail; /* one line per tier that missed its floor */
        for (size_t n : args.sizes) {
            std::vector<double> vw, vc, fw, fc;
            for (const Pt &p : pts) {
                if (p.n != n || OPS[p.opi].pol == POL_NOCHECK) continue;
                (p.cold ? vc : vw).push_back(p.ratio);
                if (p.aa_n) (p.cold ? fc : fw).push_back(p.aa_p90);
            }
            const double mw = vw.empty() ? 1.0 : median(vw);
            const double mc = vc.empty() ? 1.0 : median(vc);
            const bool have_floor = !fw.empty() || !fc.empty();
            const double flw = fw.empty() ? 0.0 : std::fabs(median(fw) - 1.0);
            const double flc = fc.empty() ? 0.0 : std::fabs(median(fc) - 1.0);
            const bool excluded = tier_excluded(n, args.tier_exclude);
            const char *tag = "";
            char tagbuf[96];
            if (excluded) {
                tag = "  [blacklisted tier: recorded, not judged]";
            } else if (!have_floor) {
                /* no same-round A/A data (--aa-each-round was off): fall back to the old
                 * constant band and say so */
                const bool bad = !(mw >= 0.90 && mw <= 1.10) || !(mc >= 0.90 && mc <= 1.10);
                if (bad) {
                    size_bad = true;
                    tag = "  <== tier drifted overall (no same-round A/A, constant band used)";
                }
            } else if (flw > 0.20 || flc > 0.20) {
                ++tier_inconc;
                tag = "  [floor >20% -> INCONCLUSIVE]";
            } else {
                ++tier_judged;
                /* Clamped floor: below AB_MIN_TIER_FLOOR the fixture's own reproducibility
                 * dominates the ratio, so judging against the raw estimate only produces false
                 * alarms (README section 8). The measured floor is still what gets printed. */
                const double ew = std::max(flw, AB_MIN_TIER_FLOOR);
                const double ec = std::max(flc, AB_MIN_TIER_FLOOR);
                const bool bw = std::fabs(mw - 1.0) > ew;
                const bool bc = std::fabs(mc - 1.0) > ec;
                if (bw || bc) {
                    size_bad = true;
                    char d[224];
                    std::snprintf(d, sizeof d,
                                  "n=%zu: warm %+.2f%% (floor %.2f%% -> effective %.2f%%), "
                                  "cold %+.2f%% (floor %.2f%% -> effective %.2f%%)",
                                  n, 100.0 * (mw - 1.0), 100.0 * flw, 100.0 * ew,
                                  100.0 * (mc - 1.0), 100.0 * flc, 100.0 * ec);
                    tier_fail.push_back(d);
                    if (ew > flw || ec > flc)
                        std::snprintf(tagbuf, sizeof tagbuf,
                                      "  <== outside the criterion (effective floor %.1f%%)",
                                      100.0 * std::max(ew, ec));
                    else
                        std::snprintf(tagbuf, sizeof tagbuf, "  <== outside the criterion");
                    tag = tagbuf;
                }
            }
            std::printf("    n=%-9zu warm %6.4f (floor %5.3f, %2zu pts)   cold %6.4f "
                        "(floor %5.3f, %2zu pts)%s\n",
                        n, mw, flw, vw.size(), mc, flc, vc.size(), tag);
        }
        std::printf("    -> judged %zu tiers, INCONCLUSIVE %zu tiers, blacklist %s\n", tier_judged,
                    tier_inconc, args.tier_exclude.c_str());
        const size_t total_res = resolvable_w + bad_w + bad_s;
        const size_t total_all = resolvable_w + bad_w + unres_w + bad_s + unres_s;
        const bool env_bad = total_res * 20 < total_all;
        const bool bias_bad = !(gm >= 0.95 && gm <= 1.05);
        if (bad_w || bad_c || bad_s > 0) {
            std::printf("  (points out of band: warm %zu / cold %zu / small-N %zu. A single point can\n"
                        "   be moved by one-off effects such as code placement or page mapping; "
                        "that\n   is the known resolution limit of this A/A fixture. The "
                        "criterion rests on\n   the per-size and all-point median ratios, not on "
                        "single points.)\n",
                        bad_w, bad_c, bad_s);
        }
        if (env_bad || bias_bad || size_bad) {
            std::printf("\n");
            std::printf("#####################################################################\n");
            if (env_bad || bias_bad) {
                /* These two are genuinely run-level properties: without enough resolvable
                 * points, or with a systematic bias, no ratio of this run can be judged. */
                std::printf("## self-check failed: %s%s\n",
                            env_bad ? "too few resolvable points (<5%) " : "",
                            bias_bad ? "all-point median ratio is more than 5% off 1.0 "
                                       "(systematic bias)" : "");
                std::printf("## resolvable %zu/%zu, loadavg=%.2f. A run in this state cannot be\n"
                            "## judged: re-run before using any ratio from it.\n",
                            total_res, total_all, loadavg1());
            } else {
                /* A tier-level miss is local evidence, not proof that the run is unusable. */
                std::printf("## self-check failed: tier-level deviation from the tier criterion\n");
                for (const std::string &t : tier_fail)
                    std::printf("##   %s\n", t.c_str());
                std::printf("## deviation = tier median ratio - 1; floor = that tier's measured A/A\n"
                            "## p90; effective floor = max(measured, %.0f%%) (see section 8 of\n"
                            "## bench/ab/README.md).\n",
                            100.0 * AB_MIN_TIER_FLOOR);
                std::printf("## The tier(s) listed above are not judged. The other tiers of this\n"
                            "## run stand and their numbers stay usable.\n");
            }
            std::printf("#####################################################################\n");
            verdict = 3;
        } else {
            std::printf("  verdict: fixture self-check passed -- the per-size and the all-point\n"
                        "           median ratios are both in band, so identical code shows no\n"
                        "           systematic difference between the two sides.\n");
        }
    } else {
        std::printf("\ncriterion reminder: this project's speedup threshold is >= 1.15x and must "
                    "always be reported together with the working-set size.\n");
    }

    /* ================= writing the artifacts ================= */
    if (!args.tsv_out.empty()) {
        FILE *f = std::fopen(args.tsv_out.c_str(), "w");
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", args.tsv_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec in-process A/B interleaving fixture  label=%s\n",
                     args.label.c_str());
        std::fprintf(f, "# git_rev=%s  baseline=git HEAD snapshot  candidate=%s\n", AB_GIT_REV,
                     AB_CAND_DESC);
        std::fprintf(f, "# machine=AMD Ryzen 7 5700U 8C/16T AVX2 no AVX-512; L1d 32KiB L2 512KiB "
                        "L3 4MiB/CCX; governor=performance %.0f-%.0f MHz not pinned; "
                        "single-threaded, taskset pinned\n",
                     mn0, mx0);
        std::fprintf(f, "# accounting: gbs=(rd+wr)/ns_per_elem; pct_memcpy=this point's GB/s / "
                        "same-size same-phase memcpy GB/s\n");
        std::fprintf(f, "# RFO accounting: real_gbs=(rd + real_wr)/ns_per_elem where real_wr=wr*"
                        "(plain store 2x / NT 1x / in-place 1x); rfo per op: --list-ops; "
                        "mixed?/libc? = uncertain, flagged as such\n");
        std::fprintf(f, "# ratio=median(per-round A/B) (main criterion); ratio_medAB=median(A)/"
                        "median(B). >1 means the candidate is faster\n");
        std::fprintf(f, "# warm=repeated runs on resident buffers; cold=one call after "
                        "madvise(DONTNEED), minor page faults included\n");
        std::fprintf(f, "# clock=%s\n", g_clock_name);
        std::fprintf(f, "# trials per point >= %d, until CoV<%.1f%% or the cap %d\n",
                     args.trials_min,
                     args.cov_target, args.trials_max);
        std::fprintf(f, "phase\top\tn\tworkingset_bytes\trd_bytes\twr_bytes\ttier\tgbs_a\tgbs_b\t"
                        "pct_memcpy_a\tpct_memcpy_b\ta_median_ns_per_elem\tb_median_ns_per_elem\t"
                        "ratio\tratio_medAB\ta_min\ta_max\tb_min\tb_max\tpair_ratio_min\t"
                        "pair_ratio_max\ta_cov_pct\tb_cov_pct\tpair_iqr_pct\ttrials\titers\t"
                        "cs_a\tcs_b\tcs_ref\tload0\tload1\tflags\t"
                        "real_rd_bytes\treal_wr_a_bytes\treal_wr_b_bytes\treal_gbs_a\t"
                        "real_gbs_b\treal_pct_memcpy_a\treal_pct_memcpy_b\trfo\t"
                        "aa_n\taa_p50\taa_p90\taa_max\taa_dev_max\tmode\tref_op\n");
        for (const Pt &p : pts) {
            /* Byte volume = per-element bytes x element count. A wide element (12 B / 16 B) has
             * nel() elements per sweep slot, not n: using n here overstates ws by 12/8 or 16/8
             * and can push the tier label across a boundary (the console path already uses
             * nel(); this is the same rule). */
            const size_t ws = (p.rd + p.wr) * nel(OPS[p.opi], p.n);
            std::fprintf(f,
                         "%s\t%s\t%zu\t%zu\t%zu\t%zu\t%s\t%.4f\t%.4f\t%.1f\t%.1f\t%.6f\t%.6f\t"
                         "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.2f\t%d\t"
                         "%ld\t%.9g\t%.9g\t%.9g\t%.2f\t%.2f\t%u\t"
                         "%zu\t%zu\t%zu\t%.4f\t%.4f\t%.1f\t%.1f\t%s\t"
                         "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%s\n",
                         p.cold ? "cold" : "warm", p.op, p.n, ws, p.rd, p.wr, tier_of(ws),
                         p.gbs_a, p.gbs_b, p.pct_memcpy_a, p.pct_memcpy_b, p.a_med, p.b_med,
                         p.ratio, p.ratio_m, p.a_min, p.a_max, p.b_min, p.b_max, p.pr_min,
                         p.pr_max, p.a_cov, p.b_cov, p.riqr, p.trials, p.iters, p.cs_a,
                         p.cs_b, p.cs_ref, p.load0, p.load1, p.flags, p.real_rd, p.real_wr_a,
                         p.real_wr_b, p.real_gbs_a, p.real_gbs_b, p.real_pct_memcpy_a,
                         p.real_pct_memcpy_b, rfo_side_name(OPS[p.opi], SIDE_A), p.aa_n,
                         p.aa_p50, p.aa_p90, p.aa_max, p.aa_dev_max,
                         p.mode ? p.mode : "ab", p.ref_op ? p.ref_op : "-");
        }
        std::fclose(f);
        std::printf("\nTSV -> %s\n", args.tsv_out.c_str());
    }

    if (!args.small_out.empty() && !small.empty()) {
        FILE *f = std::fopen(args.small_out.c_str(), "w");
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", args.small_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec small-N regime (call overhead) label=%s rev=%s "
                        "candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# ns/call is the bottleneck metric; scalar_loop = element-wise scalar API "
                        "(the pre-refactor anchor, never inlined across a TU boundary)\n");
        std::fprintf(f, "# A_scalar_x = ns/call(scalar) / ns/call(A); >1 means the batch entry is "
                        "faster\n");
        std::fprintf(f, "n\top\tcalls\tns_per_call_a\tns_per_call_b\tns_per_call_scalar\t"
                        "ns_per_elem_a\tns_per_elem_b\tratio_A_over_B\tA_over_scalar\ta_cov_pct\t"
                        "b_cov_pct\ttrials\tchecksum\n");
        for (const SmallPt &s : small) {
            std::fprintf(f, "%zu\t%s\t%ld\t%.4f\t%.4f\t%.4f\t%.5f\t%.5f\t%.4f\t%.3f\t%.2f\t%.2f\t"
                            "%d\t%.9g\n",
                         s.n, s.op, s.calls, s.a_ns_call, s.b_ns_call, s.s_ns_call,
                         s.a_ns_call / (double)s.n, s.b_ns_call / (double)s.n, s.ratio,
                         s.scalar_x, s.a_cov, s.b_cov, s.trials, s.cs);
        }
        std::fclose(f);
        std::printf("small-N TSV -> %s\n", args.small_out.c_str());
    }

    if (!args.matrix_out.empty()) {
        FILE *f = std::fopen(args.matrix_out.c_str(), "w");
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", args.matrix_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec throughput matrix, baseline snapshot  label=%s rev=%s "
                        "candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# each row is (phase, metric, op), the columns are the sizes; metric: "
                        "gbs=nominal GB/s, real_gbs=real GB/s after RFO accounting, "
                        "pct_memcpy=percent of same-size memcpy, ratio=median(per-round A/B)\n");
        std::fprintf(f, "phase\tmetric\top");
        for (size_t n : args.sizes) std::fprintf(f, "\tn%zu", n);
        std::fprintf(f, "\n");
        const char *metrics[] = {"gbs", "real_gbs", "pct_memcpy", "ratio"};
        for (int cold = 0; cold <= (args.do_cold ? 1 : 0); ++cold) {
            for (int mi = 0; mi < 4; ++mi) {
                for (const size_t opi : sel) {
                    std::fprintf(f, "%s\t%s\t%s", cold ? "cold" : "warm", metrics[mi],
                                 OPS[opi].name);
                    for (size_t n : args.sizes) {
                        const Pt *p = find_pt(opi, n, cold != 0);
                        if (!p)
                            std::fprintf(f, "\t-");
                        else if (mi == 0)
                            std::fprintf(f, "\t%.3f", p->gbs_a);
                        else if (mi == 1)
                            std::fprintf(f, "\t%.3f", p->real_gbs_a);
                        else if (mi == 2)
                            std::fprintf(f, "\t%.1f", p->pct_memcpy_a);
                        else
                            std::fprintf(f, "\t%.4f", p->ratio);
                    }
                    std::fprintf(f, "\n");
                }
            }
        }
        std::fclose(f);
        std::printf("matrix TSV -> %s\n", args.matrix_out.c_str());
    }

    if (!args.floor_out.empty()) {
        FILE *f = std::fopen(args.floor_out.c_str(), "w");
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", args.floor_out.c_str());
            return 2;
        }
        std::fprintf(f, "# A/A floor table (fixture resolution floor)  label=%s rev=%s "
                        "candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# valid only for builds where both sides are **the same** source "
                        "(SELFCHECK / candidate=git HEAD snapshot);\n");
        std::fprintf(f, "# it is the A/A self-check of that build: ratios of the same code, two "
                        "TUs, interleaved inside one process.\n");
        std::fprintf(f, "# usage: a claim that the candidate is r times faster can only be true "
                        "when (r-1) > floor_dev;\n");
        std::fprintf(f, "#        a resolvable=0 point (paired IQR >= %.1f%%) cannot even measure the "
                        "floor, so never judge with it.\n",
                     args.cov_target);
        std::fprintf(f, "# floor_dev=|ratio-1|; in_band=the ratio falls in [%.2f,%.2f]; rfo: "
                        "--list-ops\n",
                     AB_BAND_LO, AB_BAND_HI);
        std::fprintf(f, "phase\top\tn\tratio\tpair_iqr_pct\ttrials\tresolvable\tfloor_dev\t"
                        "in_band\trfo_flags\n");
        for (const Pt &p : pts) {
            const bool resolvable = (p.flags & AB_F_COV) == 0;
            const double dev = std::fabs(p.ratio - 1.0);
            const bool in_band = p.ratio >= AB_BAND_LO && p.ratio <= AB_BAND_HI;
            std::fprintf(f, "%s\t%s\t%zu\t%.6f\t%.2f\t%d\t%d\t%.6f\t%d\t0x%x\n",
                         p.cold ? "cold" : "warm", p.op, p.n, p.ratio, p.riqr, p.trials,
                         resolvable ? 1 : 0, dev, in_band ? 1 : 0, p.flags);
        }
        std::fclose(f);
        std::printf("A/A floor TSV -> %s\n", args.floor_out.c_str());
    }

    return verdict;
}
