/* bench/ab/ab_scalar_loop.cpp —— "重构前"锚点：逐元素调标量 API。
 *
 * 本 TU 只包含头文件，标量函数在别处定义 -> 每次调用都是跨 TU 的普通 call。
 * 这既是 ab_scalar_lib.cpp 存在的原因，也是本夹具最需要守住的一条：
 * 锚点必须比"批量入口"笨，否则量出来的是编译器内联的功劳。
 *
 * 覆盖范围只做小 N 区制要用的 6 个算子（add/scale/lerp/dot/len/normalize）：
 * 大 N 区制比的是两个批量实现（A/B），标量锚点由 verify/perf.cpp 负责。
 */
#include "snapshot/include/flashvec.h"

#include "ab_abi.h"

#define AB_S_ ((fv_real)-0.25)
#define AB_T_ ((fv_real)0.375)

extern "C" unsigned long long ab_s_abi(void)
{
    return ((unsigned long long)sizeof(fv_vec2) << 32) | (unsigned)sizeof(fv_real);
}

/* 标量锚点没有 NT store：所有写都是普通 store。 */
extern "C" unsigned ab_s_opflags(int op)
{
    (void)op;
    return 0u;
}

extern "C" void ab_s_loop(int op, size_t n, long iters, const AbArgs *p, float *dot)
{
    const fv_vec2 *a = (const fv_vec2 *)p->a;
    const fv_vec2 *b = (const fv_vec2 *)p->b;
    fv_vec2 *o = (fv_vec2 *)p->o;

    switch (op) {
    case AB_OP_ADD_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_add(a[i], b[i]);
        break;
    case AB_OP_SCALE_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_scale(a[i], AB_S_);
        break;
    case AB_OP_LERP_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_lerp(a[i], b[i], AB_T_);
        break;
    case AB_OP_DOT_SOA: {
        fv_real acc = (fv_real)0;
        for (long k = 0; k < iters; ++k) {
            fv_real s = (fv_real)0;
            for (size_t i = 0; i < n; ++i)
                s = s + fv_vec2_dot(fv_vec2_make(p->ax[i], p->ay[i]),
                                    fv_vec2_make(p->bx[i], p->by[i]));
            acc += s;
        }
        *dot = (float)acc;
        break;
    }
    case AB_OP_LEN_SOA:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) p->lenv[i] = fv_vec2_len(fv_vec2_make(p->ax[i], p->ay[i]));
        break;
    case AB_OP_NORM_SOA:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_normalize(fv_vec2_make(p->ax[i], p->ay[i]));
                p->ox[i] = r.x;
                p->oy[i] = r.y;
            }
        break;
    default:
        break;
    }
}
