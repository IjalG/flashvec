/* The scalar implementation of fv_vec4 -- written the same way as src/flashvec.cpp (compiled as C++, C language subset only,
 * no global state, no allocation), with the component count extended from 2/3 to 4.
 *
 * Why a separate translation unit instead of merging into src/flashvec.cpp:
 *   this round's work boundary explicitly leaves src/flashvec.cpp untouched (2D/3D are the stable baseline and the working tree must stay reviewable by hand).
 * Why not header inline: see the two reasons in the "4D" contract of include/flashvec.h
 *   (the sqrtf include set, and the scalar implementation sharing one -ffp-contract=off translation unit with the batch tail).
 *
 * Bit-level conventions (the premise for the include/flashvec.h contract and for the bit-for-bit comparison against the batch entry points):
 *   - This unit carries #pragma GCC optimize("fp-contract=off"): madd's a + b*s, lerp's
 *     a + (b-a)*t and dot/len's x*x + y*y + z*z + w*w are never contracted into an FMA;
 *   - the per-component expressions correspond word for word to the 2D/3D ones in src/flashvec.cpp, so
 *     on the same input the NaN placement of fv_vec4_min/max/clamp is exactly that of fv_vec2_min/max/clamp
 *     (the (a<b)?a:b comparison rule, not IEEE minNum);
 *   - dot sums in left-to-right order, ((x*x' + y*y') + z*z') + w*w', and len is simply sqrtf(dot(v,v)).
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

/* Per-component clamp: tests the lower bound before the upper one (the same expression as fv_clamp1 in src/flashvec.cpp,
 * and with lo > hi it places the value by the order of those same two ifs; a NaN compares false against both bounds and is returned unchanged). */
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

/* The NaN rule is exactly that of fv_vec2_min / fv_vec3_min: (a < b) ? a : b, not IEEE minNum. */
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

/* Left to right: ((x*x' + y*y') + z*z') + w*w' -- the SoA kernels of the batch entry points must reproduce this association order,
 * otherwise the bit-level comparison fails (see the bit-exactness promise boundary in include/flashvec.h). */
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
    /* denom is a squared quantity, so the squared threshold FV_EPS_SQR must be used (the same fix as in fv_vec2_project). */
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
