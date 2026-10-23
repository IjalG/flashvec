/* fv_vec4 的标量实现 —— 与 src/flashvec.cpp 同一套写法（C++ 编译、只用 C 语言子集、
 * 无全局状态、不分配内存），只是把分量从 2/3 个扩到 4 个。
 *
 * 为什么单独一个编译单元，而不是并进 src/flashvec.cpp：
 *   本轮的工作边界明确不改 src/flashvec.cpp（二维/三维是稳定基线，工作树要留给人工审查）。
 * 为什么不是头文件 inline：见 include/flashvec.h「四维」契约里的两条理由
 *   （sqrtf 的包含集、标量实现与批量尾巴共用同一个 -ffp-contract=off 编译单元）。
 *
 * 位级约定（这是与 include/flashvec.h 契约、以及批量入口逐位对拍的前提）：
 *   - 本单元带 #pragma GCC optimize("fp-contract=off")：madd 的 a + b*s、lerp 的
 *     a + (b-a)*t、dot/len 的 x*x + y*y + z*z + w*w 都不会被收缩成 FMA；
 *   - 逐分量表达式与 src/flashvec.cpp 里二维/三维的写法逐字对应，所以
 *     同一份输入下 fv_vec4_min/max/clamp 的 NaN 落位与 fv_vec2_min/max/clamp 完全一致
 *     （(a<b)?a:b 这条比较规则，不是 IEEE minNum）；
 *   - dot 的求和顺序是左结合的 ((x*x' + y*y') + z*z') + w*w'，len 就是 sqrtf(dot(v,v))。
 */
#include "flashvec.h"

#include <math.h>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

static fv_real fv4_sqrt(fv_real v)
{
    return sqrtf(v);
}

static fv_real fv4_abs(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

/* 逐分量夹取：先比下界再比上界（与 src/flashvec.cpp 的 fv_clamp1 同一表达式，
 * lo > hi 时同样按这两个 if 的顺序落位；NaN 不与边界比较为真，原样返回）。 */
static fv_real fv4_clamp1(fv_real x, fv_real lo, fv_real hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

fv_vec4 fv_vec4_make(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_vec4 v;
    v.x = x;
    v.y = y;
    v.z = z;
    v.w = w;
    return v;
}

fv_vec4 fv_vec4_add(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

fv_vec4 fv_vec4_sub(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

fv_vec4 fv_vec4_neg(fv_vec4 v)
{
    return fv_vec4_make(-v.x, -v.y, -v.z, -v.w);
}

fv_vec4 fv_vec4_scale(fv_vec4 v, fv_real s)
{
    return fv_vec4_make(v.x * s, v.y * s, v.z * s, v.w * s);
}

fv_vec4 fv_vec4_lerp(fv_vec4 a, fv_vec4 b, fv_real t)
{
    return fv_vec4_make(a.x + (b.x - a.x) * t,
                        a.y + (b.y - a.y) * t,
                        a.z + (b.z - a.z) * t,
                        a.w + (b.w - a.w) * t);
}

/* NaN 规则与 fv_vec2_min / fv_vec3_min 完全一致：(a < b) ? a : b，不是 IEEE minNum。 */
fv_vec4 fv_vec4_min(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y,
                        a.z < b.z ? a.z : b.z,
                        a.w < b.w ? a.w : b.w);
}

fv_vec4 fv_vec4_max(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y,
                        a.z > b.z ? a.z : b.z,
                        a.w > b.w ? a.w : b.w);
}

fv_vec4 fv_vec4_clamp(fv_vec4 v, fv_vec4 lo, fv_vec4 hi)
{
    return fv_vec4_make(fv4_clamp1(v.x, lo.x, hi.x),
                        fv4_clamp1(v.y, lo.y, hi.y),
                        fv4_clamp1(v.z, lo.z, hi.z),
                        fv4_clamp1(v.w, lo.w, hi.w));
}

fv_vec4 fv_vec4_madd(fv_vec4 a, fv_vec4 b, fv_real s)
{
    return fv_vec4_make(a.x + b.x * s,
                        a.y + b.y * s,
                        a.z + b.z * s,
                        a.w + b.w * s);
}

/* 左结合：((x*x' + y*y') + z*z') + w*w' —— 批量入口的 SoA 内核必须复现这个结合顺序，
 * 否则位级对拍会挂（见 include/flashvec.h 的 commit 边界）。 */
fv_real fv_vec4_dot(fv_vec4 a, fv_vec4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

fv_real fv_vec4_len_sqr(fv_vec4 v)
{
    return fv_vec4_dot(v, v);
}

fv_real fv_vec4_len(fv_vec4 v)
{
    return fv4_sqrt(fv_vec4_len_sqr(v));
}

fv_real fv_vec4_dist_sqr(fv_vec4 a, fv_vec4 b)
{
    return fv_vec4_len_sqr(fv_vec4_sub(a, b));
}

fv_real fv_vec4_dist(fv_vec4 a, fv_vec4 b)
{
    return fv4_sqrt(fv_vec4_dist_sqr(a, b));
}

fv_vec4 fv_vec4_normalize(fv_vec4 v)
{
    fv_real len = fv_vec4_len(v);
    if (len <= FV_EPS) {
        return fv_vec4_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec4_scale(v, inv);
    }
}

fv_vec4 fv_vec4_project(fv_vec4 v, fv_vec4 onto)
{
    fv_real denom = fv_vec4_len_sqr(onto);
    /* denom 是平方量，必须用平方量阈值 FV_EPS_SQR（与 fv_vec2_project 同一条修正）。 */
    if (fv4_abs(denom) <= FV_EPS_SQR) {
        return fv_vec4_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
    }
    return fv_vec4_scale(onto, fv_vec4_dot(v, onto) / denom);
}

int fv_vec4_near(fv_vec4 a, fv_vec4 b, fv_real eps)
{
    if (fv4_abs(a.x - b.x) > eps) {
        return 0;
    }
    if (fv4_abs(a.y - b.y) > eps) {
        return 0;
    }
    if (fv4_abs(a.z - b.z) > eps) {
        return 0;
    }
    if (fv4_abs(a.w - b.w) > eps) {
        return 0;
    }
    return 1;
}
