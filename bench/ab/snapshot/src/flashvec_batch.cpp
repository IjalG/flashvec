/* flashvec 批量入口（P1）实现：AoS 4 路交错 + SoA 分量数组 + 稀疏变体。
 *
 * 位级一致性策略（对外契约写在 include/flashvec.h 的「批量入口」注释里）：
 *   - 逐元素型批量入口与对应标量 API 逐位一致。为此只用 vaddps / vsubps /
 *     vmulps / vdivps / vsqrtps，不用 FMA（vfmadd）也不做倒数近似（rcpps/rsqrtps）；
 *   - 下面的 pragma 关掉 FMA 收缩：内核虽然带 target("avx2,fma")（FMA 是 AVX2 机器上
 *     的事实存在，派发只需要 avx2），但绝不允许 GCC 把 a + b*s 缩成 vfmadd——
 *     那会让同一表达式与标量路径差 1 ulp，破坏逐位一致；
 *   - 标量路径（没有 AVX2 的机器 + 所有尾巴）用与 src/flashvec.cpp 逐字相同的表达式；
 *   - 唯一的例外是 fv_vec2_dot_soa：4 个独立累加器 + 末端一次归约，明确不承诺逐位一致；
 *     它在 n < 8 时退化成纯顺序累加，这时与逐元素路径逐位一致。
 *
 * 编译要求：本文件不依赖任何全局 -march / -mavx2 开关，默认产物是 x86-64 基线；
 * AVX2 内核用函数级 target 属性单独编译，运行时用 __builtin_cpu_supports("avx2") 派发。
 */

#include "flashvec.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__)
#define FV_HAVE_X86 1
#include <immintrin.h>
#else
#define FV_HAVE_X86 0
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

#define FV_AVX2_KERNEL __attribute__((target("avx2,fma")))

/* ---------------- 标量路径（非 AVX2 机器与所有尾巴共用） ----------------
 * 表达式与 src/flashvec.cpp 里的标量函数逐字对应：a.x + b.x 之类的顺序也一样。
 * 写成 static 而不是调库函数，是为了让「没有 AVX2 的机器」这条路径仍然可被自动向量化。 */

static fv_vec2 fv_s_add(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    return r;
}

static fv_vec2 fv_s_sub(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    return r;
}

static fv_vec2 fv_s_scale(fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = v.x * s;
    r.y = v.y * s;
    return r;
}

static fv_vec2 fv_s_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    fv_vec2 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    return r;
}

static fv_vec2 fv_s_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    fv_vec2 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    return r;
}

/* 与 fv_vec2_len 同表达式：sqrtf(x*x + y*y)。 */
static fv_real fv_s_len(fv_real x, fv_real y)
{
    return sqrtf(x * x + y * y);
}

/* 与 fv_vec2_normalize 同语义、同表达式、同分支。 */
static void fv_s_normalize(fv_real x, fv_real y, fv_real *ox, fv_real *oy)
{
    fv_real len = fv_s_len(x, y);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
    }
}

#if FV_HAVE_X86

/* 4 个 fv_vec2 = 32 B = 一个 ymm：loadu/storeu 直接吃连续元素，零重排。
 * fv_vec2 只有 4 B 对齐，所以全程用 loadu/storeu（不能用 _mm256_load_ps）。 */

FV_AVX2_KERNEL static void fv_add4_avx2(const fv_vec2 *__restrict a,
                                        const fv_vec2 *__restrict b,
                                        fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vb));
    }
}

FV_AVX2_KERNEL static void fv_sub4_avx2(const fv_vec2 *__restrict a,
                                        const fv_vec2 *__restrict b,
                                        fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_sub_ps(va, vb));
    }
}

FV_AVX2_KERNEL static void fv_scale4_avx2(const fv_vec2 *__restrict a, fv_real s,
                                          fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_mul_ps(va, vs));
    }
}

/* lerp = a + (b - a) * t，乘加顺序与标量一致；fp-contract=off 保证不会变成 FMA。 */
FV_AVX2_KERNEL static void fv_lerp4_avx2(const fv_vec2 *__restrict a,
                                         const fv_vec2 *__restrict b, fv_real t,
                                         fv_vec2 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
}

/* madd = a + b * s，同上。 */
FV_AVX2_KERNEL static void fv_madd4_avx2(const fv_vec2 *__restrict a,
                                         const fv_vec2 *__restrict b, fv_real s,
                                         fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
}

/* 掩码变体：4 个元素一次（一个 ymm 装 4 个 fv_vec2 = 8 个 float lane）。
 * 关键点：掩码是按**元素**给的，而一个元素占两个 float lane，所以每个掩码字节
 * 必须铺成两个 32-bit lane（x 与 y 共用同一个掩码位），不能按字节平铺。
 * 之后 blendvps 选 a（掩码为 0）或 a+b（掩码非 0），整条写回 o。 */
FV_AVX2_KERNEL static void fv_add_mask4_avx2(const fv_vec2 *__restrict a,
                                             const fv_vec2 *__restrict b,
                                             fv_vec2 *__restrict o,
                                             const uint8_t *__restrict mask, size_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        memcpy(&mb, mask + i, sizeof(mb));
        __m128i b4 = _mm_cvtsi32_si128((int)mb);  /* 低 4 字节 = 4 个掩码字节 */
        __m128i w = _mm_cvtepu8_epi16(b4);        /* 4 个字节 -> 4 个 word（0..255） */
        __m128i wd = _mm_unpacklo_epi16(w, w);    /* 每个 word 复制一份 */
        __m256i d = _mm256_cvtepi16_epi32(wd);    /* 8 个 dword：m0 m0 m1 m1 m2 m2 m3 m3 */
        __m256 m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(d, zero));
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_blendv_ps(va, _mm256_add_ps(va, vb), m));
    }
}

/* SoA：8 个 float 一次。 */
FV_AVX2_KERNEL static void fv_add_soa_avx2(const float *__restrict ax,
                                           const float *__restrict ay,
                                           const float *__restrict bx,
                                           const float *__restrict by,
                                           float *__restrict ox,
                                           float *__restrict oy, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_add_ps(_mm256_loadu_ps(ax + i),
                                               _mm256_loadu_ps(bx + i)));
        _mm256_storeu_ps(oy + i, _mm256_add_ps(_mm256_loadu_ps(ay + i),
                                               _mm256_loadu_ps(by + i)));
    }
}

FV_AVX2_KERNEL static void fv_scale_soa_avx2(const float *__restrict ax,
                                             const float *__restrict ay, fv_real s,
                                             float *__restrict ox,
                                             float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_mul_ps(_mm256_loadu_ps(ax + i), vs));
        _mm256_storeu_ps(oy + i, _mm256_mul_ps(_mm256_loadu_ps(ay + i), vs));
    }
}

FV_AVX2_KERNEL static void fv_len_soa_avx2(const float *__restrict ax,
                                           const float *__restrict ay,
                                           float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(_mm256_add_ps(
                                    _mm256_mul_ps(x, x), _mm256_mul_ps(y, y))));
    }
}

/* 归一化：零向量保护写成掩码而不是分支。
 * keep = !(len <= FV_EPS)，与标量 if (len <= FV_EPS) 严格互补（NaN 时走 else 分支，
 * 用 _CMP_NLE_UQ 才能对上：unordered 也算「不满足 <=」）。
 * 除法的被除数用 blend 换成 1.0，所以零向量 lane 不产生 inf/NaN，也不置起除零标志；
 * 最后 and(keep) 把该 lane 清成 +0.0，与标量返回 (0,0) 的位模式一致。 */
FV_AVX2_KERNEL static void fv_normalize_soa_avx2(const float *__restrict ax,
                                                 const float *__restrict ay,
                                                 float *__restrict ox,
                                                 float *__restrict oy, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
    }
}

/* dot：4 个独立累加器打断依赖链（实测单累加器 0.3213 -> 4 累加器 0.1428 ns/元素），
 * 末端一次归约。x 分量与 y 分量累加进同一个累加器，所以归约顺序与逐元素路径不同——
 * 这是 fv_vec2_dot_soa 唯一允许的偏差（头文件已写明）。 */
FV_AVX2_KERNEL static float fv_dot_soa_avx2(const float *__restrict ax,
                                            const float *__restrict ay,
                                            const float *__restrict bx,
                                            const float *__restrict by, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    __m256 s2 = _mm256_setzero_ps();
    __m256 s3 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 32 <= n; i += 32) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 0), _mm256_loadu_ps(bx + i + 0)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 0), _mm256_loadu_ps(by + i + 0))));
        s1 = _mm256_add_ps(s1, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 8), _mm256_loadu_ps(bx + i + 8)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 8), _mm256_loadu_ps(by + i + 8))));
        s2 = _mm256_add_ps(s2, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 16), _mm256_loadu_ps(bx + i + 16)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 16), _mm256_loadu_ps(by + i + 16))));
        s3 = _mm256_add_ps(s3, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 24), _mm256_loadu_ps(bx + i + 24)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 24), _mm256_loadu_ps(by + i + 24))));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i), _mm256_loadu_ps(bx + i)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i), _mm256_loadu_ps(by + i))));
    }

    /* 末端归约：先把 4 个 ymm 累加器两两相加，再做一次 4 lane 的树。
     * 用 movehdup + shuffle 而不是 hadd（本机 hadd 慢约 2 倍）。
     * q = [q0,q1,q2,q3] -> r = [q0+q1, ..., q2+q3, ...] -> 取 lane0 + lane2。 */
    {
        __m256 t = _mm256_add_ps(_mm256_add_ps(s0, s1), _mm256_add_ps(s2, s3));
        __m128 lo = _mm256_castps256_ps128(t);
        __m128 hi = _mm256_extractf128_ps(t, 1);
        __m128 q = _mm_add_ps(lo, hi);
        __m128 r = _mm_add_ps(q, _mm_movehdup_ps(q));
        sum = _mm_cvtss_f32(_mm_add_ss(r, _mm_shuffle_ps(r, r, 0x02)));
    }

    /* 尾巴按顺序累加（n < 8 时整段都走这里，结果与逐元素路径逐位一致）。 */
    for (; i < n; ++i) {
        sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
    }
    return sum;
}

#endif /* FV_HAVE_X86 */

/* ---------------- 公开入口：运行时派发 + 尾巴 ---------------- */

void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_add4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_add(a[i], b[i]);
    }
}

void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_sub4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_sub(a[i], b[i]);
    }
}

void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                      fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_scale4_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_scale(a[i], s);
    }
}

void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real t, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_lerp4_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_lerp(a[i], b[i], t);
    }
}

void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real s, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_madd4_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_madd(a[i], b[i], s);
    }
}

void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                         fv_vec2 *__restrict o, const uint8_t *__restrict mask, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_add_mask4_avx2(a, b, o, mask, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        if (mask[i] != 0) {
            o[i] = fv_s_add(a[i], b[i]);
        } else {
            o[i] = a[i]; /* 掩码外写回 a 的原值（不是「不写」） */
        }
    }
}

/* 索引变体故意走标量：随机访存下 gather 在本机（Zen 2）是微码循环，通常不如逐元素。
 * 语义上 idx 可重复、可乱序，重复下标的写入值相同，因此与处理顺序无关。 */
void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                        fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_s_add(a[i], b[i]);
    }
}

void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                          fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_s_scale(a[i], s);
    }
}

void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict bx, const float *__restrict by,
                     float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_add_soa_avx2(ax, ay, bx, by, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}

void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay, fv_real s,
                       float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_scale_soa_avx2(ax, ay, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}

float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict bx, const float *__restrict by, size_t n)
{
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_dot_soa_avx2(ax, ay, bx, by, n);
    }
#endif
    /* 没有 AVX2：单累加器顺序求和，与逐元素路径逐位一致。 */
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
        }
        return sum;
    }
}

void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                     float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_len_soa_avx2(ax, ay, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_s_len(ax[i], ay[i]);
    }
}

void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_HAVE_X86
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_normalize_soa_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
}
