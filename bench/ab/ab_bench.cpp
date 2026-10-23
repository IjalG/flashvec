/* bench/ab/ab_bench.cpp —— flashvec 同进程 A/B（交错）性能 / 吞吐 / 调用开销夹具
 *
 * 为什么要有这个文件：README「性能」第 3 条限定说这台机器是共享机、频率会漂移，
 * 同一份代码跨运行可差 20-25%。"跑两个二进制再比"因此不可用：两次运行之间频率、
 * 邻居负载、页映射都变过。站得住的比法是**同一个进程、同一份数据、同一段时间内**
 * 交替跑两个实现（A/B/A/B...），再按档取中位数相除。
 *
 * 三个区制
 *   1) 大 N 吞吐扫描（默认 2^12 … 2^24，8 档）：ns/元素 与 GB/s，配 memcpy/memset/
 *      纯读/triad 上限参考，给出"占同规模 memcpy 吞吐的百分比"；warm 与 cold 都报。
 *   2) 小 N 调用开销区制（n = 4/8/16/32/64）：ns/调用 与 ns/元素 都给，对照组含
 *      scalar_loop（逐元素调标量 API，代表"重构前/链接静态库的调用方"）。
 *   3) 上限参考贯穿两个区制。
 *
 * 结构
 *   baseline  = git HEAD 快照（bench/ab/snapshot/，git show 生成，绝不编工作树 src/）
 *   candidate = 编译期可选：HEAD 快照本身（==baseline，夹具自检）或任意源码树
 *   每个变体各有自己的 TU（ab_variant.cpp，-include 改名宏 + 指定头文件/源文件），
 *   计时循环在变体 TU 内部，因此"库内定义"与"头文件内联"两种形态都能如实测到。
 *   调用只经 ab_X_loop 一个入口，两侧完全对等。
 *
 * 三道语义门禁（任何一道不过都不打印比值/吞吐）
 *   1) 输出与**本文件自带的独立标量参考实现**逐位比较（含掩码写回、idx 未引用位置、
 *      管线中间量）。唯一例外 fv_vec2_dot_soa：头文件明确不承诺逐位一致，按
 *      max(1e-5*|结果|, 1e-6*Σ|逐项积|) 判等——这个界与 n 和语料无关，比裸相对误差稳。
 *   2) A 与 B 的输出也必须逐位相等（例外同上）。
 *   3) 输入缓冲校验和在计时前后不变。
 *   --selfcheck 再加一道：baseline vs baseline 的比值必须落在 [0.85, 1.15]，不落就是
 *   夹具自己不可信，硬失败（退出码 3）。
 *
 * 测量口径
 *   - 单线程；库无多线程无全局状态。绑核由 Makefile 的 taskset -c 负责（cpu0-7 同一 CCX）。
 *   - N 的量纲是 fv_vec2 元素个数。一个元素在 ax/ay 中各占一个 float。
 *   - B/元素 = 该函数**实际读到 + 写到**的字节数（逐函数人工核对，见 OPS 表与报告）。
 *     GB/s = (rd+wr) / (ns/元素)，即每秒真实移动的字节数（读写都算）。
 *   - warm = 已驻留缓冲区反复跑；cold = 对被测数组 madvise(DONTNEED) 后单次调用
 *     （含 minor page fault；页被内核清零，故 cold 只用于吞吐对照，不参与位级门禁）。
 *     idx 数组永不丢弃——丢了之后所有倒排写都落到下标 0，随机写变热点写。
 *   - 缓存层级按本机实测：L1d 32 KB、L2 512 KB、L3 **4 MB/CCX**（index3=4096K，
 *     shared_cpu_list=0-7；文档里那句"8 MB L3"是错的，一个线程只吃得到自己那 4 MB）。
 *   - 每点最少 5 次取样取中位数，到变异系数 < 5% 或到上限为止，实际次数进 TSV。
 *
 * 出口：控制台表格 + TSV（--tsv-out 长表、--matrix-out 矩阵、--small-out 小 N 表）
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

/* 驱动只 include HEAD 快照的头文件，而 fv_vec4（vec3/vec4 的 P1 新增类型）在 HEAD 里
 * 还不存在。为参考实现补一个同布局的本地定义；布局一致性由变体 TU 的 ABI 指纹与
 * 头文件里的 16 B static_assert 保证（这里再钉一次）。等 HEAD 有了 fv_vec4，
 * Makefile 的自动探测会定义 AB_SNAP_HAS_VEC4=1，这段自动消失。 */
#ifndef AB_SNAP_HAS_VEC4
#define AB_SNAP_HAS_VEC4 0
#endif
#if !AB_SNAP_HAS_VEC4
struct fv_vec4 {
    fv_real x, y, z, w;
};
static_assert(sizeof(fv_vec4) == 16, "fv_vec4 必须是 16 B（与头文件契约一致）");
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
#define AB_ABS_FRAC 1e-6 /* dot_soa 的绝对界：1e-6 * Σ|逐项积| */
#define AB_REL_FRAC 1e-5 /* dot_soa 的相对界：1e-5 * |结果| */
#define AB_T 0.375f
#define AB_S (-0.25f)

/* =====================================================================
 * 一、环境探针
 * ===================================================================== */
/* 计时基准。默认用**本线程 CPU 时间**而不是墙钟：这台机器长期被别的负载占 30-50%，
 * 墙钟会把"被调度走了多久"算进被测代码，在 DRAM 档上能占到 20% 以上。
 * CPU 时间只统计本线程真正在核上跑的时长——它仍然包含被邻居抢 L3 造成的停顿（那是真实
 * 影响），但不包含被抢占的时间。--clock=mono 可以切回墙钟自己核对差值。 */
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
        std::fprintf(stderr, "分配失败: %zu 字节\n", bytes);
        std::exit(2);
    }
    return p;
}

/* =====================================================================
 * 二、缓冲区与数据
 * ===================================================================== */
/* 输出形态。KIND 同时决定**元素宽度**（见 elem_bytes）：VEC2/SOA2/SOA4/FLOAT/NONE
 * 的元素是 fv_vec2（8 B），SOA3 是 12 B（fv_vec3），VEC4/SOA4V 是 16 B（fv_vec4）。
 * 同一档扫描里所有 op 处理**相同的字节量**：元素数 = 扫描规模 n(8B 单位) * 8 / elem_bytes，
 * 于是 vec4 与 vec2 的比值天然是按字节归一的对照。 */
enum OutKind { OUT_VEC2, OUT_SOA2, OUT_SOA3, OUT_SOA4, OUT_SOA4V, OUT_VEC4, OUT_FLOAT, OUT_NONE };

/* 元素宽度不再从 KIND 推导：归约型 op（dot/len）的输出是 float/NONE，但元素宽度
 * 仍可能是 12/16 B。宽度是 op 的身份，写在注册表每一行的 ELOP 列里。 */
enum Policy { POL_EXACT, POL_REL, POL_NOCHECK };
enum Side { SIDE_A = 0, SIDE_B = 1, SIDE_S = 2 };

struct OutSet {
    void *base;  /* 8n：o / ox,oy / lenv 共用（同一时刻只用一种形状） */
    void *base2; /* 8n：collision 的 (vox,voy) | vec3 的 oz | vec4 的 oz,ow */
    fv_vec2 *o;
    float *ox, *oy, *lenv;
    float *vx, *vy;    /* = base2 的前后两半，供 collision 用 */
    float *oz, *ow;    /* = base2 的前后两半，供 vec3/vec4 SoA 用（与 vx/vy 同址） */
    float dot;
    double abs_sum; /* dot_soa 参考实现的 Σ|逐项积|，用于绝对容差 */
};

struct Buf {
    size_t n;
    size_t alloc_bytes;
    fv_vec2 *a, *b;
    uint8_t *mask;
    uint32_t *idx;
    float *ax, *ay, *bx, *by;
    float *px, *py, *qx, *qy;
    /* 宽元素（vec3 12 B / vec4 16 B）专用的分量数组：只有选中宽 op 时才分配，
     * 默认扫描（全是 vec2）一个字节都不多花。 */
    float *az, *bz, *aw, *bw;
    OutSet out[2];
    double in_cs;
};

/* 掩码密度（--mask-density）：0.5 = 旧行为 */
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

/* 选中了 vec3/vec4 op 时为真：分配 az/bz/aw/bw 与 base2（vec4 SoA 的 oz/ow 要用）。
 * 默认扫描（全是 vec2）不分配，内存足迹与以前完全一致。 */
static bool g_wide = false;

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
    const bool need2 = with_pipe || g_wide;
    /* 两套输出必须连续分配、相距恰好 n*8 字节。n*8 是 4096 的整数倍，于是 out[0] 与
     * out[1] 相对 a/b 的 mod-4096 偏移完全相同——否则 A 与 B 会落进不同的
     * store/load 4K 别名关系，在 L1/L2 驻留档能凭空差出 20%（实测踩过）。 */
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
    std::free(B.out[0].base); /* 两套输出是同一块内存的前后两半 */
    std::free(B.out[0].base2);
    std::memset(&B, 0, sizeof B);
}

/* 管线 scratch 的确定性重灌：collision 拿它当法线/速度输入，chain3 拿它当 t/u。
 * 因为会互相覆盖，px..qy 不计入输入校验和。 */
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

/* 数据只取决于 n 与固定种子。每 512 个元素放一个精确零向量，专门覆盖 normalize 的
 * 零向量保护分支——否则"keep 掩码取反"这类错误在本夹具里完全不可见。 */
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
            /* 掩码密度可调（--mask-density）：_mask / _mask_inplace 的"稀疏度"要能扫。
             * 用 24 位尾数比较，密度 0.5 时与旧版 ((h>>11)&1) 的统计性质相同。 */
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
        /* 宽元素的分量数组：用同一条确定性序列，语义上是"第 3/4 个分量" */
        for (size_t i = 0; i < n; ++i) {
            const uint32_t h = (uint32_t)(i * 2654435761u);
            B.az[i] = (float)((double)(h & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.bz[i] = (float)((double)((h >> 16) & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.aw[i] = (float)((double)((h ^ 0x5bd1e995u) & 0xFFFFu) * (1.0 / 32768.0)) - 1.0f;
            B.bw[i] = (float)((double)(((h >> 13) ^ 0x1b873593u) & 0xFFFFu) * (1.0 / 32768.0)) -
                      1.0f;
        }
        /* 每 512 个元素放一个精确零向量（覆盖 normalize 的零保护，vec3/vec4 同样适用） */
        for (size_t i = 0; i < n; i += 512) {
            B.az[i] = 0.0f;
            B.aw[i] = 0.0f;
        }
    }
    init_pipe(B);
    for (int k = 0; k < 2; ++k) prefill(B.out[k], n); /* 预触碰：RSS 核对与 warm 口径都靠它 */
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        s += (double)B.a[i].x + (double)B.a[i].y + (double)B.b[i].x + (double)B.b[i].y;
        s += (double)B.mask[i] + (double)B.idx[i];
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
 * 三、独立参考实现（纯标量、顺序、不调用库的批量入口）
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
                  "%s 不一致: 共 %zu/%zu 个元素位模式不同；首个 i=%zu got(x,y)=%s,%s want=%s,%s",
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
                  "%s 不一致: 共 %zu/%zu 个值位模式不同；首个 i=%zu got=%s want=%s", what, bad, n,
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


/* ================= vec3 / vec4 参考实现（纯标量、顺序、不调库） ================= */
#if AB_OPS_EXTRA /* 这些入口默认不在注册表里，参考实现跟着一起编掉 */
/* vec3：元素 12 B，三个分量数组。零向量保护与 fv_vec2_normalize 同语义（len<=FV_EPS → 0）。 */
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
/* vec4 AoS（16 B/元素）：2 个元素装满一个 ymm，所以是 _x2。 */
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

/* ---- 每个 op 的参考实现（签名统一，放进 OPS 表） ---- */
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
/* 就地稀疏加：mask[i] != 0 时 o[i] = a[i] + b[i]，否则 o[i] = a[i]（保持原值）。
 * 驱动把 O.o 当作那个"可写数组"传进去（见 exec 的 AB_OPF_INPLACE 分支）。
 * 暂时只服务实验区的 _mask_inplace 行；把它挪进默认表时去掉下面这个 #if。 */
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
    /* 用 double 累加并同时记下 Σ|逐项积|：容差要按"项的量级"定，不能按结果的相对值定
     * （结果接近抵消时相对误差会爆，n 越大越明显，实测 n=2^18 就能到 7e-5）。 */
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
/* ---- 上限参考（不经过变体 TU；A/B 跑同一段代码，比值即纯噪声探针） ---- */
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
 * 四、op 表
 * ===================================================================== */
struct Op {
    const char *name;
    size_t rd, wr; /* 实际每元素读/写字节（名义；真实流量另算，见 real_wr_bytes） */
    size_t elop;   /* 元素宽度（字节）：vec2 8 / vec3 12 / vec4 16 —— 由 KIND 决定 */
    OutKind kind;
    Policy pol;
    RefFn ref;       /* 参考实现；POL_NOCHECK 时为空 */
    RefFn local_run; /* 不经过变体 TU 的实现（上限参考）；为空则走 ab_X_loop */
    int abi_op;      /* AB_OP_*；-1 表示 local_run */
    unsigned flags;  /* AB_OPF_* / AB_MF_* */
    int rfo;         /* AB_RFO_*：普通 store 的写分配行为（真实流量记账用） */
    const char *note;/* RFO 记账的人读依据 */
};

/* 表由 ab_ops.h 的注册表展开：驱动里没有任何 per-op 代码，加 op 不改这里。 */
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

/* 注册表增删自动跟随，不写死数字（19 是 15 批量入口 + 4 上限参考的老账） */
#define AB_COUNT_ITEM(...) +1
enum {
    N_LIB_OPS = 0 AB_LIB_OPS(AB_COUNT_ITEM) AB_LIB_OPS_EXTRA(AB_COUNT_ITEM),
    N_MEM_OPS = 0 AB_MEM_OPS(AB_COUNT_ITEM)
};
#undef AB_COUNT_ITEM
static_assert(sizeof(OPS) / sizeof(OPS[0]) == N_LIB_OPS + N_MEM_OPS, "OPS 表与注册表不一致");

static bool has_flag(const Op &op, unsigned f)
{
    return (op.flags & f) != 0;
}

/* 逐规模判据里要剔除的档（逗号分隔；默认 131072 = L3 边界，A/A 地板本身不可分辨） */
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

/* 该 op 在扫描规模 n（8 B 单位）下处理多少个**自己的元素**。
 * n*8 是所有 op 共同的字节量，所以不同元素宽度的 op 之间是字节归一的。 */
static size_t nel(const Op &op, size_t n)
{
    return op.elop ? (n * 8u) / op.elop : n;
}

/* 上限参考按**角色**取用，不按下标（往表里插行不会错位） */
static size_t memref_index(unsigned role)
{
    for (size_t k = 0; k < OP_COUNT; ++k)
        if (OPS[k].abi_op < 0 && (OPS[k].flags & role)) return k;
    return OP_COUNT;
}

/* INTRA 对（ab_ops.h 的 AB_INTRA_PAIRS）：新入口 ID -> 参照入口 ID */
struct IntraPair {
    int a, b;
};
static const IntraPair kIntra[] = {
#define AB_INTRA_ITEM(A, B) {AB_OP_##A, AB_OP_##B},
    AB_INTRA_PAIRS(AB_INTRA_ITEM)
#undef AB_INTRA_ITEM
    {-1, -1} /* 哨兵：AB_OPS_EXTRA=0 时表为空也必须合法 */
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

/* 输出数组的**真实**写流量（B/元素）：
 *   NT store        -> wr            （不经过 cache、不产生 RFO）
 *   写目标已被读过   -> wr            （就地更新；行本来就在 cache 里）
 *   libc 自己选      -> wr            （不确定，标注 libc?，不参与判定）
 *   普通 store       -> wr * 2        （写一次 + RFO 把该 line 读进来一次）
 * 输入（读）永远按 1x。 */
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
        return "miss*"; /* 本侧没有 NT 入口，退化成普通 store */
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
    return p;
}

/* 就地算子：输出数组先被灌成 a 的内容，然后**它自己**就是那个可写数组。
 * 这样 B.a 永远不被写，"输入校验和不变"这条门禁对就地算子同样成立。 */
static void seed_inplace(const Op &op, Buf &B, OutSet &O)
{
    if (!has_flag(op, AB_OPF_INPLACE)) return;
    std::memcpy(O.base, B.a, B.n * sizeof(fv_vec2));
}

static void exec(const Op &op, int side, Buf &B, OutSet &O, long iters)
{
    if (op.kind == OUT_NONE) O.dot = 0.0f; /* 归约型算子把 dot 当累加器用 */
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
 * 五、校验和 / 门禁
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
    case OUT_VEC4: {
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
        /* fv_vec4 是 4 个 float，逐位比较等价于把它当 4n 个 float 比 */
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
                      "结果超出容差: got=%s want=%s |diff|=%.4g > tol=%.4g (=max(%.0e*|r|, "
                      "%.0e*Σ|项|))",
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

    if (op.pol == POL_NOCHECK) { /* 上限参考：A/B 同一段代码，只验位级相等 */
        prefill(P, n);
        exec(op, SIDE_A, B, P, 1);
        prefill(R, n);
        exec(op, SIDE_B, B, R, 1);
        g.cs_a = out_sum(op, P, nel(op, n));
        g.cs_b = out_sum(op, R, nel(op, n));
        std::string m = cmp_out(op, P, R, nel(op, n));
        if (!m.empty()) {
            g.ok = false;
            g.why = "  A 与 B 输出不一致（上限参考本应完全相同）: " + m + "\n";
        }
        g.cs_ref = g.cs_a;
        return g;
    }

    const size_t ne = nel(op, n); /* 本 op 的元素数（字节量恒定，宽度不同则元素数不同） */
    prefill(R, n);
    op.ref(B, ne, R);

    if (has_flag(op, AB_OPF_NO_BASELINE)) {
        /* 基线（git HEAD）里根本没有这个入口：A 侧不参与，只验 candidate vs 独立标量参考。
         * 这条路径的结论只能说"实现与标量语义一致"，不能说"比基线快"。 */
        const double abs0 = R.abs_sum;
        prefill(P, n);
        exec(op, SIDE_B, B, P, 1);
        g.cs_b = g.cs_a = out_sum(op, P, ne);
        g.cs_ref = out_sum(op, R, ne);
        std::string mb = cmp_out(op, P, R, ne);
        if (!mb.empty()) {
            g.ok = false;
            g.why += "  candidate vs 参考: " + mb + "\n";
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
        if (!m_ar.empty()) g.why += "  baseline vs 参考: " + m_ar + "\n";
        if (!m_br.empty()) g.why += "  candidate vs 参考: " + m_br + "\n";
    }
    if (g.ok && with_scalar && op.kind != OUT_SOA4) { /* 标量锚点也过一遍门禁 */
        OutSet S = R;
        S.dot = R.dot;
        prefill(S, n);
        exec(op, SIDE_S, B, S, 1);
        std::string ms = cmp_out(op, S, R, ne);
        if (!ms.empty()) {
            g.ok = false;
            g.why += "  scalar_loop 锚点 vs 参考: " + ms + "\n";
        }
    }
    R.abs_sum = abs_sum;
    return g;
}

/* =====================================================================
 * 六、计时与统计
 * ===================================================================== */
static double median(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t m = v.size();
    if (m & 1u) return v[m / 2];
    return 0.5 * (v[m / 2 - 1] + v[m / 2]);
}

/* 配对差比的四分位距（相对中位数）。比"样本 CoV"稳健得多：这台机器上偶发的
 * 多倍停顿会把 std 拉爆，于是 CoV 判所有点都不可信，而实际比值分布很紧。
 * 可分辨性判据因此用本函数，CoV 仍然照报（父任务点名要的指标）。 */
static double iqr_pct(std::vector<double> v)
{
    if (v.size() < 4) return 0.0;
    std::sort(v.begin(), v.end());
    const double q25 = v[v.size() / 4];
    const double q75 = v[(v.size() * 3) / 4];
    const double m = median(v);
    return (m > 0.0) ? (q75 - q25) / m * 100.0 : 0.0;
}

/* 分位数（线性插值）。地板要给 p50/p90/max 三个数，不能只给 max。 */
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

/* 返回 ns/元素（或 ns/调用，调用方自己乘除） */
static double measure(const Op &op, int side, Buf &B, OutSet &O, long iters)
{
    const double t0 = now_ns();
    exec(op, side, B, O, iters);
    const double t1 = now_ns();
    const size_t m = nel(op, B.n); /* 元素数按本 op 的元素宽度算（字节量恒定） */
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
    /* RFO 记账：真实流量 = rd + wr*(RFO?2:1)。两侧分别算，因为一侧可能是 NT store
     * 而另一侧退化成普通 store（见 ab_variant.cpp 的 *_opflags）。 */
    double real_gbs_a, real_gbs_b, real_pct_memcpy_a, real_pct_memcpy_b;
    size_t real_rd, real_wr_a, real_wr_b;
    int rfo; /* AB_RFO_*（行上的声明；实际是否 NT 由 rfo_side_name 报） */
    const char *mode;   /* "ab"（A/B）或 "intra"（candidate 内部对照） */
    const char *ref_op; /* intra 模式下的参照 op 名；ab 模式为 "-" */
    /* 同一轮里与 A/B 交错跑的 A/A 配对（R4 的分母：同分布、同数据、同一段时间） */
    int aa_n;
    double aa_p50, aa_p90, aa_max, aa_dev_max;
    int trials;
    long iters;
    double cs_a, cs_b, cs_ref;
    double load0, load1;
    unsigned flags;
    size_t rd, wr;
};

/* 小 N 区制覆盖的算子：注册表里带 AB_OPF_SMALL 的（标量锚点 ab_s_loop 已实现）。
 * 用标志而不是下标，往表里插行/加 op 都不会让覆盖面悄悄变掉。 */
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
    std::string op_exact; /* --op=<精确名>：verify 用，避免 substring 误匹配 */
    std::string label;
    std::string tsv_out, matrix_out, small_out, floor_out;
    int trials_min, trials_max;
    double cov_target, target_ms, small_target_ms;
    bool selfcheck, do_cold, do_pipe, do_small;
    bool gate_only; /* 只跑位级门禁，不出任何计时结果（verify 的 R1） */
    bool list_ops;  /* 打印注册表后退出（verify 用它判断 op 在不在） */
    /* 每一轮 A/B 之后紧接一对 A/A：R4 的分母必须与 A/B 同分布（同进程/同数据/同一段时间） */
    bool aa_each_round;
    bool intra;     /* INTRA 模式：新入口 vs candidate 内部的同字节量对照 op */
    double mask_density; /* _mask / _mask_inplace 的掩码密度 */
    std::string tier_exclude; /* 逐规模判据里剔除的档（默认 131072） */
    int order;    /* 0=交错(默认) 1=总是 A 先 2=总是 B 先 */
    bool swap_out; /* 交换 A/B 使用的输出缓冲（诊断缓冲不对称） */
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
    /* 默认固定"先 A 后 B"。理由（实测）：在 2^20 上 order=ab 得 ratio 0.995 / CoV 9.7%，
     * order=ba 得 1.030 / 34.5%，交错得 1.305 / 24.3%——固定顺序把"每个配对处在两种不同
     * 进入状态"这个混杂因素去掉了，配对比值分布不再是两个总体的混合；而 ab 与 ba 都
     * 落在 1.0 附近，说明顺序本身没有系统性偏置。要交叉核对就用 --order=ba。 */
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
                "用法: %s [选项]\n"
                "  --sizes=4096,...       大 N 吞吐扫描（默认 2^12..2^24 共 8 档）\n"
                "  --small=4,8,16,32,64   小 N 调用开销区制（n 的取值）\n"
                "  --ops=子串              只看匹配的入口\n"
                "  --op=精确名              只看这一个入口（verify 用）\n"
                "  --gate-only            只跑位级门禁（不出计时/吞吐），不通过则退出码 4\n"
                "  --list-ops             打印 op 注册表后退出\n"
                "  --floor-out=FILE       落盘 A/A 地板表（只对 SELFCHECK 构建有意义）\n"
                "  --trials-min=5 --trials-max=12 --cov=5.0  自适应取样\n"
                "  --target-ms=8          大 N 每取样目标时长；--small-target-ms=20\n"
                "  --selfcheck            baseline vs baseline 自检（逐规模判据用该档 A/A 地板 p90）\n"
                "  --aa-each-round        每轮 A/B 后紧跟一对 A/A（R4 与逐规模判据的分母）\n"
                "  --tier-exclude=131072  逐规模判据剔除的档（只记录不判）\n"
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
            std::fprintf(stderr, "未知参数: %s（--help 看用法）\n", a.c_str());
            return false;
        }
    }
    if (A.trials_min < 5)
        std::fprintf(stderr, "警告: --trials-min=%d < 5，违反「每档至少 5 次取样」的要求\n",
                     A.trials_min);
    return true;
}

/* =====================================================================
 * 七、主流程
 * ===================================================================== */
int main(int argc, char **argv)
{
    Args args;
    if (!parse_args(argc, argv, args)) return 2;

    /* ABI 指纹：A/B/标量锚点三份 TU 的 fv_vec2 / fv_real 布局必须与驱动一致 */
    {
        const unsigned long long want =
            ((unsigned long long)sizeof(fv_vec2) << 32) | (unsigned)sizeof(fv_real);
        const unsigned long long ga = ab_a_abi(), gb = ab_b_abi(), gs = ab_s_abi();
        if (ga != want || gb != want || gs != want) {
            std::printf("#####################################################################\n");
            std::printf("## 硬失败：ABI 指纹不一致。驱动 %llu，A %llu，B %llu，标量锚点 %llu\n",
                        want, ga, gb, gs);
            std::printf("## 意思是某个变体 TU 的 fv_vec2 / fv_real 与驱动不同（头文件被改过布局？）。\n");
            std::printf("## 这种情况下所有指针都会被解释成错误的结构，不能测，直接退出。\n");
            std::printf("#####################################################################\n");
            return 5;
        }
    }

    /* --list-ops：注册表就是"有哪些 op 可测"的唯一真源；verify.sh 与文档都用它，
     * 不在这里口述列表（口述的列表一定会和表漂移）。 */
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

    /* 显式点名一个"基线里不存在"的入口时，自动切到 INTRA 模式：
     * A/B 对它们无意义（A 侧根本没有这个函数）。 */
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
        /* 新 API（vec3/vec4）在基线里不存在：A/B 模式把它们排除，INTRA 模式只要它们 */
        if (args.intra != has_flag(OPS[k], AB_OPF_NO_BASELINE)) continue;
        if (!args.op_exact.empty()) {
            if (args.op_exact != OPS[k].name) continue;
        } else if (!args.ops_filter.empty() &&
                   !std::strstr(OPS[k].name, args.ops_filter.c_str()))
            continue;
        sel.push_back(k);
    }
    /* 宽元素（12 B / 16 B）才需要 az/bz/aw/bw 与 base2 */
    for (size_t k : sel)
        if (OPS[k].elop != 8) g_wide = true;
    if (sel.empty()) {
        std::fprintf(stderr, "没有入口匹配 --ops=%s%s%s（--list-ops 看注册表）\n",
                     args.ops_filter.c_str(), args.op_exact.empty() ? "" : " --op=",
                     args.op_exact.c_str());
        return 2;
    }
    if (args.intra)
        std::printf("# 模式: INTRA（candidate 内部对照；**不是**相对 git HEAD 的改进）\n");

    /* --gate-only：只跑位级门禁。verify 用它做"全尾长"对拍（n = N-7 … N+3 覆盖
     * _x4 的 0-3、_soa 的 0-7 尾巴），不产生任何计时数字。 */
    if (args.gate_only) {
        std::printf("# gate-only：%zu 个 op × %zu 个规模（位级对拍，A/B/独立参考三者一致）\n",
                    sel.size(), args.sizes.size());
        size_t nfail = 0, ntot = 0;
        for (size_t n : args.sizes) {
            Buf B;
            alloc_buf(B, n, args.do_pipe); /* 管线算子要用 px/py/qx/qy scratch */
            init_data(B, n);
            for (size_t k : sel) {
                const Op &op = OPS[k];
                if (has_flag(op, AB_OPF_PIPE)) init_pipe(B); /* 上一个 op 可能写过 scratch */
                const double in_cs = B.in_cs;
                Gate g = check_op(op, B, n, small_supported((int)k));
                bool ok = g.ok;
                if (B.in_cs != in_cs) {
                    ok = false;
                    g.why += "  实现改写了输入缓冲（a/b/mask/idx）\n";
                }
                ++ntot;
                if (!ok) {
                    ++nfail;
                    std::printf("[GATE FAIL] %s @ n=%zu\n%s", op.name, n, g.why.c_str());
                }
            }
            free_buf(B);
        }
        std::printf("# gate-only 结果: %zu/%zu 通过\n", ntot - nfail, ntot);
        return nfail ? 4 : 0;
    }

    double mn0, mx0, av0;
    cpu_mhz(&mn0, &mx0, &av0);
    const double load_start = loadavg1();
    const long rss0 = self_statm(2);

    std::printf("=====================================================================\n");
    std::printf("flashvec 同进程 A/B 交错夹具   label=%s\n", args.label.c_str());
    std::printf("  baseline  : git HEAD 快照 bench/ab/snapshot/  rev=%s\n", AB_GIT_REV);
    std::printf("  candidate : %s\n", AB_CAND_DESC);
    std::printf("  编译选项  : 与顶层 Makefile 默认一致（无 -march；AVX2 靠函数级派发）\n");
    std::printf("  入口      : %zu / %zu（注册表：AB_LIB_OPS %d 个%s + AB_MEM_OPS %d 个上限参考）\n",
                sel.size(), OP_COUNT, (int)N_LIB_OPS, args.do_pipe ? "（含管线）" : "（无管线）",
                (int)N_MEM_OPS);
    std::printf("  大 N 扫描 :");
    for (size_t n : args.sizes) std::printf(" 2^%.0f", std::log2((double)n));
    std::printf("\n  小 N 区制 :");
    for (size_t n : args.small_ns) std::printf(" %zu", n);
    if (!args.do_small) std::printf(" (关闭)");
    std::printf("\n");
    std::printf("  取样      : 每点最少 %d 次、最多 %d 次，到【配对差比 IQR】< %.1f%% 为止；大 N 目标 %.1f "
                "ms/次、小 N 目标 %.1f ms/次\n",
                args.trials_min, args.trials_max, args.cov_target, args.target_ms,
                args.small_target_ms);
    std::printf("  缓存层级  : L1d 32 KiB / L2 512 KiB / L3 4 MiB（单 CCX，cpu0-7）\n");
    std::printf("  线程      : 单线程；库无多线程无全局状态；绑核见 Makefile 的 taskset\n");
    std::printf("  频率      : governor=performance，观测 %.0f-%.0f MHz（均 %.0f）-- 未固定，"
                "绝对值不可跨会话比\n",
                mn0, mx0, av0);
    std::printf("  负载      : loadavg(1m)=%.2f；MemAvailable=%ld MiB；VmSwap=%ld MiB；起始 RSS=%ld MiB\n",
                load_start, mem_available_mib(), proc_kb("VmSwap") / 1024, rss0 * 4 / 1024);
    std::printf("  计时基准  : %s\n", g_clock_name);
    std::printf("  口径      : GB/s=(读字节+写字节)/时间；%%memcpy=该点 GB/s / 同规模 memcpy；"
                "ratio=median(每轮 A/B)\n");
    std::printf("=====================================================================\n");

    {
        const double t0 = now_ns();
        volatile double acc = 0.0;
        while (now_ns() - t0 < 120e6) acc += 1.0;
        (void)acc;
    }

    int verdict = 0; /* 0=通过 3=自检失败（仍然落盘，方便复核） */
    std::vector<Pt> pts;
    bool any_gate_fail = false;
    size_t n_gate = 0, n_gate_fail = 0;

    /* ================= 小 N 区制：调用开销 ================= */
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
        std::printf("\n=== 小 N 区制（调用开销；缓冲区常驻 L1，每次调用只碰 n 个元素） ===\n");
        std::printf("  ns/call 才是这一区制的瓶颈指标；ns/elem 会把调用开销摊薄。\n");
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
                    (void)measure(op, SIDE_A, SB, SB.out[0], 1000); /* 先热 i-cache/分支预测 */
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
                sp.a_cov = iqr_pct(pr); /* 小 N 里这个字段就是配对 IQR% */
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

    /* ================= 大 N 吞吐扫描 ================= */
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

        std::printf("\n=== 规模 n=%zu  分配 %.1f MiB  RSS 增 %.1f MiB%s  loadavg=%.2f ===\n", n,
                    (double)B.alloc_bytes / 1048576.0, (double)rss_grew / 1048576.0,
                    rss_bad ? "  [!! RSS 明显小于预期，可能被换出，本档吞吐不可信]" : "",
                    loadavg1());

        std::vector<Gate> gates(OP_COUNT);
        for (size_t k : sel) {
            const double in_cs = B.in_cs;
            gates[k] = check_op(OPS[k], B, n, false);
            ++n_gate;
            if (B.in_cs != in_cs) {
                gates[k].ok = false;
                gates[k].why += "  实现改写了输入缓冲（a/b/mask/idx）\n";
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
            std::printf("## 硬失败：校验和门禁不通过（%zu/%zu 项）。\n", n_gate_fail, n_gate);
            std::printf("## 上面每一项都指名了函数：baseline / candidate / 独立标量参考实现在该\n");
            std::printf("## 函数上算的不是同一个东西。任何比值与吞吐都没有意义——\n");
            std::printf("## 本次不输出任何计时/吞吐结果。\n");
            std::printf("#####################################################################\n");
            free_buf(B);
            return 4;
        }

        /* 一个测量点：自适应取样 + 交错 */
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
                /* 先预热再校准：第一次调用的 i-cache / 分支预测是冷的，拿它当校准基准
                 * 会把内层重复次数估小一个数量级（踩过）。 */
                (void)measure(op, SIDE_A, B, B.out[0], 4);
                /* measure 返回 ns/本 op 的元素；乘元素数才是"一次调用"的时长 */
                const double one = measure(op, SIDE_A, B, B.out[0], 1) * (double)nel(op, n);
                iters = (long)(args.target_ms * 1e6 / std::max(one, 1.0));
                if (iters < 1) iters = 1;
                if (iters > 200000) iters = 200000;
                /* 单次调用就已经超过 30 ms 的点（DRAM 档）不追 12 次取样：那是纯烧时间 */
                if (one > 30e6) tmax = std::min(args.trials_max, 7);
                (void)measure(op, SIDE_A, B, B.out[0], iters);
                (void)measure(op, SIDE_B, B, B.out[1], iters);
            }
            for (int t = 0; t < tmax; ++t) {
                /* 停止条件看"配对差比的四分位距"，不是样本 CoV：我们要判定的是比值，
                 * 样本 std 被偶发停顿一拉就爆，会在比值已经很稳的时候还继续烧时间。 */
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
                    /* 同一轮、同一段代码、同一份数据里再跑一对 A/A：这就是 R4 的分母。
                     * 它和 A/B 处在同一浸入状态与同一频率漂移段里，才是同分布的地板。 */
                    const double u = measure(op, SIDE_A, B, OA, cold ? 1 : iters);
                    const double v = measure(op, SIDE_A, B, OB, cold ? 1 : iters);
                    aa.push_back(u / v);
                }
            }
            p.a_med = median(ta);
            p.b_med = median(tb);
            p.ratio = median(pr); /* 主判据：配对消漂移 */
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


        /* ---- INTRA：candidate 内部对照（新入口 vs 同字节量的既有入口） ----
         * 两个 op 都走 **B 侧**（candidate）。ratio = t_ref / t_new，>1 表示新入口更快。
         * 这不是"相对基线"，因为基线里没有这些入口。 */
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
                /* 按字节归一：两个 op 的元素宽度不同，x/y 直接比是"每元素"而不是"每字节"。
                 * 一次调用的总时长 = ns/元素 * 元素数，两者处理相同字节量，所以这样才是
                 * "同一份工作量谁更快"。 */
                pr.push_back((x * (double)nel(rf, n)) / (y * (double)nel(nw, n)));
                if (args.aa_each_round) {
                    /* A/A 地板用**参照 op** 打：它代表这一档的"同代码两次"噪声 */
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

        /* ---- INTRA 报表（显式标注：candidate 内部对照，不是相对 HEAD 的改进） ---- */
        if (args.intra) {
            std::printf("\n--- INTRA 对照：candidate 内部（新入口 vs 同字节量的既有入口）---\n");
            std::printf("--- **这不是相对 git HEAD 基线的改进**：这些入口在 HEAD 里不存在。   ---\n");
            std::printf("--- 两个 op 在同一档处理相同字节量；ratio = t(ref)/t(new)，>1 表示新入口更快 ---\n");
            std::printf("%-24s %-22s %7s %7s %10s %10s %6s %9s %9s %-14s %s\n", "op(new)", "ref",
                        "n(new)", "n(ref)", "ns/e(new)", "ns/e(ref)", "ratio", "GB/s(new)",
                        "GB/s(ref)", "AA p50/p90/max", "flags");
            for (size_t k : sel) {
                const Op &nw = OPS[k];
                if (!has_flag(nw, AB_OPF_NO_BASELINE)) continue;
                const size_t iref = intra_ref_of(k);
                if (iref >= OP_COUNT) {
                    std::printf("%-24s %-22s （注册表里没有 INTRA 对，跳过）\n", nw.name, "-");
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

        /* 上限参考（%memcpy 的分母；memset/纯读/triad 用同一条公式，故只需两个分母）。
         * real_* 用 memcpy 的**真实**流量作分母——memcpy 自己在大规模上很可能就是 NT store
         * （libc?），所以这一列的绝对值不可全信，但同一档内比较可用。 */
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

        /* warm 全跑完（数据还完好），再跑 cold（cold 会把页清空，避免污染后面的 warm） */
        std::printf("%-24s %5s %7s %5s %8s %9s %6s %10s %9s %9s %6s %-15s %3s %-9s %-14s %s\n",
                    "op", "B/e", "realB/e", "tier", "GB/s", "realGB/s", "%mem", "CoV/IQR%",
                    "A_med", "B_med", "ratio", "pair[min,max]", "tr", "rfo",
                    "AA p50/p90/max", "flags");
        for (size_t k : sel) {
            const Op &op = OPS[k];
            if (op.pol == POL_NOCHECK) continue;
            if (has_flag(op, AB_OPF_PIPE)) init_pipe(B);
            const size_t ws = (op.rd + op.wr) * nel(op, n);
            Pt pw = run_point(k, false);
            apply_pct(pw, memcpy_warm, memcpy_real_warm);
            /* 就地算子的输出状态在多次计时调用后已经"跑飞"（每次都在同一个缓冲上累加），
             * 所以计时后不再比 A/B；位级结论由计时前的 check_op 给。 */
            if (!has_flag(op, AB_OPF_INPLACE)) {
                std::string m =
                    cmp_out(op, B.out[0], B.out[1], n); /* 计时后重比，顺带读遍输出 */
                if (!m.empty()) {
                    std::printf("  [GATE FAIL] %s @ n=%zu: 计时后 A/B 输出不一致: %s\n",
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
        std::printf("  -- 上限参考（同一档、同一计时器；A/B 同一段代码，比值即噪声） --\n");
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
        std::printf("  -- cold（madvise(DONTNEED) 后单次调用，含 page fault） --\n");
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
        std::printf("  (本档 minor fault %ld 次；释放 %.1f MiB，VmSwap=%ld MiB)\n",
                    (mf > mf_before) ? (mf - mf_before) : -1,
                    (double)B.alloc_bytes / 1048576.0, proc_kb("VmSwap") / 1024);
        free_buf(B);
        if (any_gate_fail) return 4;
    }

    /* ================= 汇总 ================= */
    std::printf("\n=== 汇总 ===\n");
    std::printf("校验和门禁: %zu/%zu 项通过（A、B、独立参考实现三者一致；dot_soa 走 "
                "max(%.0e*|r|, %.0e*Σ|项|)）\n",
                n_gate - n_gate_fail, n_gate, (double)AB_REL_FRAC, (double)AB_ABS_FRAC);
    std::printf("大 N 测量点: %zu（warm+cold）  小 N 测量点: %zu  负载: loadavg(1m) %.2f -> %.2f\n",
                pts.size(), small.size(), load_start, loadavg1());

    auto find_pt = [&](size_t opi, size_t n, bool cold) -> const Pt * {
        for (const Pt &p : pts)
            if (p.op == OPS[opi].name && p.n == n && p.cold == cold) return &p;
        return nullptr;
    };

    for (int cold = 0; cold <= (args.do_cold ? 1 : 0); ++cold) {
        std::printf("\n--- 吞吐矩阵：%s，GB/s（单线程；rd+wr 都算） ---\n", cold ? "cold" : "warm");
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
        std::printf("\n--- 占同规模 memcpy 吞吐的百分比：%s ---\n", cold ? "cold" : "warm");
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

    std::printf("\n--- candidate/baseline 比值矩阵（warm；>1 = candidate 更快） ---\n");
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

    std::printf("\n--- 真实内存流量矩阵（warm，GB/s；RFO 记账后：普通 store 写按 2x、NT store 按 1x） ---\n");
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

    std::printf("\n--- RFO 记账：名义 vs 真实内存流量（每元素字节；real/nom = 真实/名义） ---\n");
    std::printf("%-24s %4s %4s %8s %8s %9s  %-9s %s\n", "op", "rd", "wr", "realA", "realB",
                "real/nom", "rfo", "依据");
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
    std::printf("  （real/nom 就是「NT store 的真实收益」的量化：同一份名义 24 B/e 的 add 算子，\n");
    std::printf("    普通 store 真实 32 B/e、NT store 真实 24 B/e —— 差的是流量不是速度；\n");
    std::printf("    A/B 两列不同 = 一侧退化成普通 store（baseline 头文件还没有 _stream）。）\n");

    std::printf("\n--- 吞吐拐点（warm）：相对自身峰值掉到 50%% 的规模 ---\n");
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
        std::printf("  %-24s 峰值 %6.2f GB/s @ n=%-8zu", OPS[opi].name, peak, peak_n);
        for (size_t n : args.sizes) {
            const Pt *p = find_pt(opi, n, false);
            if (p && n > peak_n && p->gbs_a < peak * 0.5) {
                std::printf("  掉一半 @ n=%zu (%.2f GB/s)", n, p->gbs_a);
                break;
            }
        }
        std::printf("\n");
    }

    if (!small.empty()) {
        std::printf("\n--- 小 N 区制：ns/调用（batch A）与相对标量锚点的倍数 ---\n");
        std::printf("%-6s", "op");
        for (size_t n : args.small_ns) std::printf(" %10zu", n);
        std::printf("\n");
        for (size_t k : sel) {
            if (!small_supported((int)k)) continue;
            bool any = false;
            for (const SmallPt &s : small)
                if (s.op == OPS[k].name) any = true;
            if (!any) continue;
            std::printf("%-6s", OPS[k].name);
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
        std::printf("\n--- 小 N 区制：A / scalar_loop 的加速比（>1 = 批量入口赢） ---\n");
        std::printf("%-6s", "op");
        for (size_t n : args.small_ns) std::printf(" %10zu", n);
        std::printf("\n");
        for (size_t k : sel) {
            if (!small_supported((int)k)) continue;
            std::printf("%-6s", OPS[k].name);
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
        std::printf("\n取样稳定性: 大 N 里配对 IQR < %.1f%%（=可分辨）的点 %zu/%zu（cold 不达标 %zu）；"
                    "所需取样数 中位 %d 最大 %d\n",
                    args.cov_target, ok, tot, cold_bad, tr[tr.size() / 2], tr.back());
    }

    if (args.selfcheck) {
        /* 判定规则分两类：
         *   可分辨   = 两侧 CoV 都 <= 阈值 -> 比值必须落在 [0.85,1.15]，否则夹具判失败；
         *   不可分辨 = 任一侧 CoV > 阈值（外部负载把这一档咬穿了）-> 照报数字，但不进判据。
         * 另加一条护栏：可分辨点必须占多数，否则说明整台机器都不适合测量，一样判失败。 */
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
        std::printf("\n自检（baseline vs baseline，比值必须落在 [%.2f, %.2f]）:\n", AB_BAND_LO,
                    AB_BAND_HI);
        if (!wr.empty()) {
            std::sort(wr.begin(), wr.end());
            std::printf("  warm: %zu 点  min=%.4f  p50=%.4f  max=%.4f  可分辨 %zu 超带 %zu  "
                        "不可分辨 %zu\n",
                        wr.size(), wr.front(), wr[wr.size() / 2], wr.back(), resolvable_w, bad_w,
                        unres_w);
        }
        if (!cr.empty()) {
            std::sort(cr.begin(), cr.end());
            std::printf("  cold: %zu 点  min=%.4f  p50=%.4f  max=%.4f  不可分辨 %zu 超带 %zu\n",
                        cr.size(), cr.front(), cr[cr.size() / 2], cr.back(), unres_c, bad_c);
        }
        std::printf("  （不可分辨 = 两侧取样变异系数 > %.1f%%，外部负载把该档咬穿；这些点的\n"
                    "    比值照报但**不作为判据**，它们的 GB/s 量级仍然可用）\n", args.cov_target);
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
        std::printf("  比值分布（大 N，warm+cold，共 %zu 点）:\n", all.size());
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
        /* 全点无偏性：240 个配对测下来，中位数必须贴着 1.0。这一条不受"个别点被噪声
         * 咬穿"影响（中位数对尾部极不敏感），是比"每点都在带内"更硬的无偏性证据。 */
        std::vector<double> allr = wr;
        allr.insert(allr.end(), cr.begin(), cr.end());
        for (const SmallPt &s : small) allr.push_back(s.ratio);
        const double gm = allr.empty() ? 1.0 : median(allr);
        std::printf("  全点中位比: %.4f（%zu 个配对测量，判据 [0.95, 1.05]）%s\n", gm,
                    allr.size(), (gm >= 0.95 && gm <= 1.05) ? "  OK" : "  <== 超出");
        /* 环境可用性下限：可分辨点若少到没有统计意义，本次结论不成立。5% 是"至少还有
         * 几十个点能判"的下限；本机被外部负载占 30-50% 时通常只有 5-15% 的点可分辨。 */
        /* 按规模看偏置：同一档里 15 个算子的比值中位数，是"这一档两侧是否等价"的
         * 稳健估计（15 个样本的中位数，标准误约等于单点噪声的 1/4）。逐点的极端值
         * 仍然照报，但"某一档整体偏了"这种夹具级故障一定会被这一行抓住。 */
        /* 逐规模判据改成**地板制**：一档的比值中位数必须落在该档自己测出来的 A/A 地板
         * p90 以内，而不是写死的 [0.90, 1.10]。地板是同一批运行里测的，比常数准。
         * 兜底：该档地板本身 > 20% → 这一档不可判，标 INCONCLUSIVE（不 PASS 也不 FAIL）。
         * 黑名单档（默认 131072，L3 边界，A/A 在 0.50–2.22）只记录、不参与判定。 */
        std::printf("  逐规模中位比（warm / cold）—— 判据 = 该档实测 A/A 地板 p90（不再是写死的 "
                    "[0.90,1.10]）:\n");
        bool size_bad = false;
        size_t tier_inconc = 0, tier_judged = 0;
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
            if (excluded) {
                tag = "  [黑名单档：只记录 不判]";
            } else if (!have_floor) {
                /* 没有同轮 A/A 数据（例如没开 --aa-each-round）：退回旧的常数带，并注明 */
                const bool bad = !(mw >= 0.90 && mw <= 1.10) || !(mc >= 0.90 && mc <= 1.10);
                if (bad) {
                    size_bad = true;
                    tag = "  <== 该档整体偏离（无同轮 A/A，退回常数带）";
                }
            } else if (flw > 0.20 || flc > 0.20) {
                ++tier_inconc;
                tag = "  [地板 >20% → INCONCLUSIVE]";
            } else {
                ++tier_judged;
                const bool bw = std::fabs(mw - 1.0) > flw;
                const bool bc = std::fabs(mc - 1.0) > flc;
                if (bw || bc) {
                    size_bad = true;
                    tag = "  <== 该档中位比超出该档 A/A 地板";
                }
            }
            std::printf("    n=%-9zu warm %6.4f (地板 %5.3f, %2zu 点)   cold %6.4f (地板 %5.3f, "
                        "%2zu 点)%s\n",
                        n, mw, flw, vw.size(), mc, flc, vc.size(), tag);
        }
        std::printf("    -> 可判 %zu 档，INCONCLUSIVE %zu 档，黑名单 %s\n", tier_judged,
                    tier_inconc, args.tier_exclude.c_str());
        const size_t total_res = resolvable_w + bad_w + bad_s;
        const size_t total_all = resolvable_w + bad_w + unres_w + bad_s + unres_s;
        const bool env_bad = total_res * 20 < total_all;
        const bool bias_bad = !(gm >= 0.95 && gm <= 1.05);
        if (bad_w || bad_c || bad_s > 0) {
            std::printf("  （逐点超带：warm %zu / cold %zu / 小N %zu。单点被『代码放置/页映射』这类\n"
                        "    一次性因素影响，是这套 A/A 夹具已知的分辨率下限；判据建立在\n"
                        "    『逐规模中位比』与『全点中位比』上，不在单点上。）\n",
                        bad_w, bad_c, bad_s);
        }
        if (env_bad || bias_bad || size_bad) {
            std::printf("\n");
            std::printf("#####################################################################\n");
            std::printf("## 自检失败：%s%s%s\n",
                        env_bad ? "可分辨点太少（<5%）" : "",
                        bias_bad ? "全点中位比偏离 1.0 超 5%（系统性偏置）" : "",
                        size_bad ? "存在超出该档 A/A 地板的规模档" : "");
            std::printf("## 可分辨 %zu/%zu，loadavg=%.2f。夹具或环境不可用，本次所有比值作废。\n",
                        total_res, total_all, loadavg1());
            std::printf("#####################################################################\n");
            verdict = 3;
        } else {
            std::printf("  结论: 夹具自检通过 —— 逐规模中位比与全点中位比都在带内，"
                        "同代码两侧无系统性差异。\n");
        }
    } else {
        std::printf("\n判据提醒：本项目的加速比门槛 >= 1.15x，且必须与工作集大小一起报。\n");
    }

    /* ================= 落盘 ================= */
    if (!args.tsv_out.empty()) {
        FILE *f = std::fopen(args.tsv_out.c_str(), "w");
        if (!f) {
            std::fprintf(stderr, "无法写入 %s\n", args.tsv_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec 同进程 A/B 交错夹具  label=%s\n", args.label.c_str());
        std::fprintf(f, "# git_rev=%s  baseline=git HEAD 快照  candidate=%s\n", AB_GIT_REV,
                     AB_CAND_DESC);
        std::fprintf(f, "# 机器=AMD Ryzen 7 5700U 8C/16T AVX2 无AVX-512; L1d 32KiB L2 512KiB "
                        "L3 4MiB/CCX; governor=performance %.0f-%.0f MHz 未固定; 单线程 taskset 绑核\n",
                     mn0, mx0);
        std::fprintf(f, "# 口径: gbs=(rd+wr)/ns_per_elem; pct_memcpy=该点 GB/s / 同规模同 phase 的 memcpy GB/s\n");
        std::fprintf(f, "# RFO 记账: real_gbs=(rd + real_wr)/ns_per_elem，real_wr=wr*(普通 store 2x / NT 1x / "
                        "就地 1x)；rfo 见 --list-ops，mixed?/libc? = 不确定 已标注\n");
        std::fprintf(f, "# ratio=median(每轮 A/B)（主判据）；ratio_medAB=median(A)/median(B)。>1 表示 candidate 更快\n");
        std::fprintf(f, "# warm=已驻留缓冲区重复跑；cold=madvise(DONTNEED) 后单次调用（含 minor page fault）\n");
        std::fprintf(f, "# clock=%s\n", g_clock_name);
        std::fprintf(f, "# 每点取样 >= %d 次，到 CoV<%.1f%% 或上限 %d 次\n", args.trials_min,
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
            const size_t ws = (p.rd + p.wr) * p.n;
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
            std::fprintf(stderr, "无法写入 %s\n", args.small_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec 小 N 区制（调用开销）label=%s rev=%s candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# ns/call 是瓶颈指标；scalar_loop = 逐元素调标量 API（重构前锚点，跨 TU 不可内联）\n");
        std::fprintf(f, "# A_scalar_x = ns/call(scalar) / ns/call(A)，>1 表示批量入口更快\n");
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
            std::fprintf(stderr, "无法写入 %s\n", args.matrix_out.c_str());
            return 2;
        }
        std::fprintf(f, "# flashvec 吞吐矩阵 baseline 快照  label=%s rev=%s candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# 每行 (phase, metric, op)，列是规模；metric: gbs=名义 GB/s, "
                        "real_gbs=RFO 记账后的真实 GB/s, pct_memcpy=占同规模 memcpy 百分比, "
                        "ratio=median(每轮 A/B)\n");
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
            std::fprintf(stderr, "无法写入 %s\n", args.floor_out.c_str());
            return 2;
        }
        std::fprintf(f, "# A/A 地板表（fixture resolution floor）  label=%s rev=%s candidate=%s\n",
                     args.label.c_str(), AB_GIT_REV, AB_CAND_DESC);
        std::fprintf(f, "# 只有两侧源码**相同**的构建（SELFCHECK / candidate=git HEAD 快照）本表才成立；\n");
        std::fprintf(f, "# 本表是对应构建的 A/A 自检结果：同一份代码、两个 TU、同一进程交替跑出来的比值。\n");
        std::fprintf(f, "# 用法：任何「候选比基线快 r 倍」的结论，只有 (r-1) > floor_dev 才可能为真；\n");
        std::fprintf(f, "#       resolvable=0 的点（配对 IQR >= %.1f%%）连地板都测不准，不要拿它判。\n",
                     args.cov_target);
        std::fprintf(f, "# floor_dev=|ratio-1|；in_band=比值落在 [%.2f,%.2f]；rfo 见 --list-ops\n",
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
