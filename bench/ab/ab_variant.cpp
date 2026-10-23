/* bench/ab/ab_variant.cpp —— 变体 TU 模板：同一份源码编两次（A / B）
 *
 * 每个变体 TU 自己包含某一版 include/flashvec.h（需要时再包含 src/flashvec_batch.cpp），
 * 并导出**整段计时循环**（ab_X_loop）。这样：
 *   库内定义形态  -> 循环里是跨函数 call（= 链接静态库的调用方拿到的形态）
 *   头文件内联形态 -> 循环里被内联（= 重构想要的形态）
 * 驱动（ab_bench.cpp）只认 ab_X_loop / ab_X_abi / ab_X_opflags 三个符号，完全不知道
 * 批量入口叫什么、有几个。
 *
 * 本文件里**直接写库里的名字**（fv_vec2_add_x4 等）：Makefile 用 -include 强制包含一份
 * 自动生成的改名表（mkrename.sh 从本侧头文件抓出所有 fv_xxx( 名字），于是同一份源码在
 * A 侧落到 ab_a_*、B 侧落到 ab_b_*，不需要任何中性别名，也不需要手工维护改名表。
 *
 * 编译参数（全部由 Makefile 传）：
 *   -include <build>/rename_a.h|<build>/rename_b.h   ← 必须先于头文件生效（宏改名）
 *   -DAB_HDR='"…/flashvec.h"'              ← 本变体的头文件
 *   -DAB_WITH_SRC=1 -DAB_SRC='"…/flashvec_batch.cpp"'  ← 库内定义形态才需要
 *   -DAB_HAVE_STREAM=0|1                   ← 本侧头文件是否已有 _stream 入口（自动探测）
 *   -DAB_LOOP_NAME=ab_a_loop -DAB_ABI_NAME=ab_a_abi -DAB_FLAGS_NAME=ab_a_opflags
 *
 * 加一个 op：在 ab_ops.h 的注册表里加一行 + 在下面加一个 abt_<ID> thunk。就这样。
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include AB_HDR /* 本变体的公开头文件 */

#if defined(AB_WITH_SRC) && AB_WITH_SRC
#include AB_SRC /* 库内定义的形态：批量入口的实现 */
#endif

#include "ab_abi.h"

#if !defined(AB_LOOP_NAME) || !defined(AB_ABI_NAME) || !defined(AB_FLAGS_NAME)
#error "需要 -DAB_LOOP_NAME / -DAB_ABI_NAME / -DAB_FLAGS_NAME"
#endif

/* 本侧头文件有没有新入口（Makefile 自动 grep 得到；没有就退化，不报错） */
#ifndef AB_HAVE_STREAM
#define AB_HAVE_STREAM 0
#endif
#ifndef AB_HAVE_INPLACE
#define AB_HAVE_INPLACE 0
#endif
/* vec3/vec4 入口（基线 HEAD 里没有；Makefile 自动 grep 本侧头文件） */
#ifndef AB_HAVE_WIDE
#define AB_HAVE_WIDE 0
#endif

#define AB_T_ ((fv_real)0.375)
#define AB_S_ ((fv_real)-0.25)

/* 管线里属于"调用方"的那段标量代码：与 ab_bench.cpp 的参考实现逐字相同，
 * 这样 A/B 差异只可能来自库调用本身。 */
static inline void abv_normalize(fv_real x, fv_real y, fv_real *ox, fv_real *oy)
{
    const fv_real len = sqrtf(x * x + y * y);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
    } else {
        const fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
    }
}
static inline void abv_reflect(fv_real nx, fv_real ny, fv_real vx, fv_real vy, fv_real *ox,
                               fv_real *oy)
{
    const fv_real d = (vx * nx + vy * ny) * (fv_real)2;
    *ox = vx - d * nx;
    *oy = vy - d * ny;
}

/* =====================================================================
 * thunk：一次调用。名字必须是 abt_<注册表 ID>（注册表用宏拼出来）。
 * 类型一律 (const AbArgs *, size_t n, float *dot)：dot 只给归约型算子用。
 * ===================================================================== */
typedef void (*AbThunk)(const AbArgs *p, size_t n, float *dot);

static void abt_ADD_X4(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_add_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
}
static void abt_SUB_X4(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_sub_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
}
static void abt_SCALE_X4(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_scale_x4((const fv_vec2 *)p->a, AB_S_, (fv_vec2 *)p->o, n);
}
static void abt_LERP_X4(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_lerp_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, AB_T_, (fv_vec2 *)p->o, n);
}
static void abt_MADD_X4(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_madd_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, AB_S_, (fv_vec2 *)p->o, n);
}
static void abt_ADD_MASK(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_add_x4_mask((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, p->mask,
                        n);
}
static void abt_ADD_IDX(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_add_x4_idx((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, p->idx, n);
}
static void abt_SCALE_IDX(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_scale_x4_idx((const fv_vec2 *)p->a, AB_S_, (fv_vec2 *)p->o, p->idx, n);
}
static void abt_ADD_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_add_soa(p->ax, p->ay, p->bx, p->by, p->ox, p->oy, n);
}
static void abt_SCALE_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_scale_soa(p->ax, p->ay, AB_S_, p->ox, p->oy, n);
}
static void abt_DOT_SOA(const AbArgs *p, size_t n, float *dot)
{
    /* 累加到 *dot：内联形态下编译器才有机会把"纯函数 + 循环不变量"提到循环外
     * （-O3 的 LICM 会这么干），那时 measure 到的就不是 n 次调用了。
     * 驱动在每次 run 前把 dot 清零（见 ab_bench.cpp 的 exec）。 */
    *dot = *dot + fv_vec2_dot_soa(p->ax, p->ay, p->bx, p->by, n);
}
static void abt_LEN_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_len_soa(p->ax, p->ay, p->lenv, n);
}
static void abt_NORM_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_normalize_soa(p->ax, p->ay, p->ox, p->oy, n);
}
static void abt_COLLISION(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec2_normalize_soa(p->px, p->py, p->ox, p->oy, n);
    for (size_t i = 0; i < n; ++i)
        abv_reflect(p->ox[i], p->oy[i], p->qx[i], p->qy[i], &p->ox2[i], &p->oy2[i]);
}
static void abt_CHAIN3(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    /* t = scale(a); u = add(t,b); o = normalize(u)，中间数组用 px/py 与 qx/qy */
    fv_vec2_scale_soa(p->ax, p->ay, AB_S_, p->px, p->py, n);
    fv_vec2_add_soa(p->px, p->py, p->bx, p->by, p->qx, p->qy, n);
    fv_vec2_normalize_soa(p->qx, p->qy, p->ox, p->oy, n);
}

/* ---- 实验区（AB_OPS_EXTRA=1 才编） ---- */
#if AB_OPS_EXTRA
static void abt_ADD_X4_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_add_x4_stream((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#else
    /* 本侧头文件还没有 _stream：退化成普通 store。这一行于是就是
     * "普通 store vs NT store" 的对照（真实流量 32 vs 24 B/e）。 */
    fv_vec2_add_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#endif
}
static void abt_SUB_X4_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_sub_x4_stream((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#else
    fv_vec2_sub_x4((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->o, n);
#endif
}
static void abt_SCALE_X4_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_scale_x4_stream((const fv_vec2 *)p->a, AB_S_, (fv_vec2 *)p->o, n);
#else
    fv_vec2_scale_x4((const fv_vec2 *)p->a, AB_S_, (fv_vec2 *)p->o, n);
#endif
}
static void abt_ADD_SOA_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_add_soa_stream(p->ax, p->ay, p->bx, p->by, p->ox, p->oy, n);
#else
    fv_vec2_add_soa(p->ax, p->ay, p->bx, p->by, p->ox, p->oy, n);
#endif
}
static void abt_SCALE_SOA_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_scale_soa_stream(p->ax, p->ay, AB_S_, p->ox, p->oy, n);
#else
    fv_vec2_scale_soa(p->ax, p->ay, AB_S_, p->ox, p->oy, n);
#endif
}
static void abt_NORM_SOA_STREAM(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_STREAM
    fv_vec2_normalize_soa_stream(p->ax, p->ay, p->ox, p->oy, n);
#else
    fv_vec2_normalize_soa(p->ax, p->ay, p->ox, p->oy, n);
#endif
}
static void abt_ADD_MASK_INPLACE(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
#if AB_HAVE_INPLACE
    fv_vec2_add_x4_mask_inplace((fv_vec2 *)p->a, (const fv_vec2 *)p->b, p->mask, n);
#else
    /* 退化：_mask 的"掩码外写回 a[i]"在 o == a 时与就地语义逐位相同。
     * （形式上 o 与 a 重叠违反 __restrict 契约，但该内核逐 4 元素读-改-写、无跨元素
     *   依赖，实测逐位一致；门禁会兜住任何不一致。） */
    fv_vec2_add_x4_mask((const fv_vec2 *)p->a, (const fv_vec2 *)p->b, (fv_vec2 *)p->a, p->mask,
                        n);
#endif
}

/* ---- vec3 SoA（元素 12 B）---- */
#if AB_HAVE_WIDE
static void abt_V3_ADD_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec3_add_soa(p->ax, p->ay, p->az, p->bx, p->by, p->bz, p->ox, p->oy, p->oz, n);
}
static void abt_V3_SCALE_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec3_scale_soa(p->ax, p->ay, p->az, AB_S_, p->ox, p->oy, p->oz, n);
}
static void abt_V3_DOT_SOA(const AbArgs *p, size_t n, float *dot)
{
    *dot = *dot + fv_vec3_dot_soa(p->ax, p->ay, p->az, p->bx, p->by, p->bz, n);
}
static void abt_V3_LEN_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec3_len_soa(p->ax, p->ay, p->az, p->lenv, n);
}
static void abt_V3_NORM_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec3_normalize_soa(p->ax, p->ay, p->az, p->ox, p->oy, p->oz, n);
}
/* ---- vec4 AoS _x2（元素 16 B）---- */
static void abt_V4_ADD_X2(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_add_x2((const fv_vec4 *)p->a, (const fv_vec4 *)p->b, (fv_vec4 *)p->o, n);
}
static void abt_V4_SUB_X2(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_sub_x2((const fv_vec4 *)p->a, (const fv_vec4 *)p->b, (fv_vec4 *)p->o, n);
}
static void abt_V4_SCALE_X2(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_scale_x2((const fv_vec4 *)p->a, AB_S_, (fv_vec4 *)p->o, n);
}
static void abt_V4_LERP_X2(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_lerp_x2((const fv_vec4 *)p->a, (const fv_vec4 *)p->b, AB_T_, (fv_vec4 *)p->o, n);
}
static void abt_V4_MADD_X2(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_madd_x2((const fv_vec4 *)p->a, (const fv_vec4 *)p->b, AB_S_, (fv_vec4 *)p->o, n);
}
/* ---- vec4 SoA ---- */
static void abt_V4_DOT_SOA(const AbArgs *p, size_t n, float *dot)
{
    *dot = *dot + fv_vec4_dot_soa(p->ax, p->ay, p->az, p->aw, p->bx, p->by, p->bz, p->bw, n);
}
static void abt_V4_LEN_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_len_soa(p->ax, p->ay, p->az, p->aw, p->lenv, n);
}
static void abt_V4_NORM_SOA(const AbArgs *p, size_t n, float *dot)
{
    (void)dot;
    fv_vec4_normalize_soa(p->ax, p->ay, p->az, p->aw, p->ox, p->oy, p->oz, p->ow, n);
}
#else
/* baseline 头文件里没有这些入口：A 侧永远不会被调用（AB_OPF_NO_BASELINE），
 * 这里给空 thunk 只是为了让同一份 ab_variant.cpp 在两侧都能编过。 */
#define AB_WIDE_STUB(ID)                                                                       \
    static void abt_##ID(const AbArgs *p, size_t n, float *dot)                                \
    {                                                                                          \
        (void)p;                                                                               \
        (void)n;                                                                               \
        (void)dot;                                                                             \
    }
AB_WIDE_STUB(V3_ADD_SOA)
AB_WIDE_STUB(V3_SCALE_SOA)
AB_WIDE_STUB(V3_DOT_SOA)
AB_WIDE_STUB(V3_LEN_SOA)
AB_WIDE_STUB(V3_NORM_SOA)
AB_WIDE_STUB(V4_ADD_X2)
AB_WIDE_STUB(V4_SUB_X2)
AB_WIDE_STUB(V4_SCALE_X2)
AB_WIDE_STUB(V4_LERP_X2)
AB_WIDE_STUB(V4_MADD_X2)
AB_WIDE_STUB(V4_DOT_SOA)
AB_WIDE_STUB(V4_LEN_SOA)
AB_WIDE_STUB(V4_NORM_SOA)
#endif /* AB_HAVE_WIDE */
#endif /* AB_OPS_EXTRA */


/* =====================================================================
 * 注册表 -> 本 TU 的分派表（由 ab_ops.h 的同一份表生成，新增 op 只多一个 thunk）
 * ===================================================================== */
struct AbRow {
    int op;
    AbThunk fn;
};

static const AbRow kRows[] = {
#define AB_ROW_ITEM(ID, NAME, KIND, ELOP, POL, REF, FLAGS, RD, WR, RFO, NOTE)                   \
    {AB_OP_##ID, abt_##ID},
    AB_LIB_OPS(AB_ROW_ITEM)
    AB_LIB_OPS_EXTRA(AB_ROW_ITEM)
#undef AB_ROW_ITEM
};

static AbThunk find_thunk(int op)
{
    for (size_t i = 0; i < sizeof(kRows) / sizeof(kRows[0]); ++i)
        if (kRows[i].op == op) return kRows[i].fn;
    return nullptr;
}

extern "C" unsigned long long AB_ABI_NAME(void)
{
    return ((unsigned long long)sizeof(fv_vec2) << 32) | (unsigned)sizeof(fv_real);
}

extern "C" unsigned AB_FLAGS_NAME(int op)
{
    switch (op) {
#define AB_FLAGS_ITEM(ID, NAME, KIND, ELOP, POL, REF, FLAGS, RD, WR, RFO, NOTE)                 \
    case AB_OP_##ID:                                                                            \
        return ((RFO) == AB_RFO_NT && AB_HAVE_STREAM) ? AB_SIDE_NT : 0u;
        AB_LIB_OPS(AB_FLAGS_ITEM)
        AB_LIB_OPS_EXTRA(AB_FLAGS_ITEM)
#undef AB_FLAGS_ITEM
    default:
        return 0u;
    }
}

extern "C" void AB_LOOP_NAME(int op, size_t n, long iters, const AbArgs *p, float *dot)
{
    const AbThunk fn = find_thunk(op);
    if (!fn) return;
    for (long k = 0; k < iters; ++k) fn(p, n, dot);
}
