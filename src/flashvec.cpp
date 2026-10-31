#define FV_BUILD 1
/* Implementation of flashvec: compiled as C++, using only the C language subset. */
#include "flashvec.h"

#include <math.h>

#if defined(__GNUC__) && !defined(__clang__)
/* Forbid contracting a + b*c into an FMA. The per-element API promises bit-exact reproducibility (docs/optimization-report.md §4.7),
 * and FMA contraction would make the same expression differ by 1 ulp between a consumer with FMA enabled and one without --
 * the batch entry points (src/flashvec_batch.cpp) carry the same pragma, which is what makes the two sides bit-identical. */
#pragma GCC optimize ("fp-contract=off")
#endif

static fv_real fv_sqrt(fv_real v)
{
    return sqrtf(v);
}

static fv_real fv_abs(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

/* Per-component clamp: tests the lower bound before the upper one, matching the expression written in the header. */
static fv_real fv_clamp1(fv_real x, fv_real lo, fv_real hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

/* ---------------- 2D ---------------- */

fv_vec2 fv_vec2_make(fv_real x, fv_real y)
{
    fv_vec2 v;
    v.x = x;
    v.y = y;
    return v;
}

fv_vec2 fv_vec2_add(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x + b.x, a.y + b.y);
}

fv_vec2 fv_vec2_sub(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x - b.x, a.y - b.y);
}

fv_vec2 fv_vec2_neg(fv_vec2 v)
{
    return fv_vec2_make(-v.x, -v.y);
}

fv_vec2 fv_vec2_scale(fv_vec2 v, fv_real s)
{
    return fv_vec2_make(v.x * s, v.y * s);
}

fv_vec2 fv_vec2_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    return fv_vec2_make(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
}

fv_vec2 fv_vec2_min(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y);
}

fv_vec2 fv_vec2_max(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y);
}

fv_vec2 fv_vec2_clamp(fv_vec2 v, fv_vec2 lo, fv_vec2 hi)
{
    return fv_vec2_make(fv_clamp1(v.x, lo.x, hi.x),
                        fv_clamp1(v.y, lo.y, hi.y));
}

fv_vec2 fv_vec2_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    return fv_vec2_make(a.x + b.x * s,
                        a.y + b.y * s);
}

fv_real fv_vec2_dot(fv_vec2 a, fv_vec2 b)
{
    return a.x * b.x + a.y * b.y;
}

fv_real fv_vec2_len_sqr(fv_vec2 v)
{
    return fv_vec2_dot(v, v);
}

fv_real fv_vec2_len(fv_vec2 v)
{
    return fv_sqrt(fv_vec2_len_sqr(v));
}

fv_real fv_vec2_dist_sqr(fv_vec2 a, fv_vec2 b)
{
    return fv_vec2_len_sqr(fv_vec2_sub(a, b));
}

fv_real fv_vec2_dist(fv_vec2 a, fv_vec2 b)
{
    return fv_sqrt(fv_vec2_dist_sqr(a, b));
}

fv_vec2 fv_vec2_normalize(fv_vec2 v)
{
    fv_real len = fv_vec2_len(v);
    if (len <= FV_EPS) {
        return fv_vec2_make((fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec2_scale(v, inv);
    }
}

fv_vec2 fv_vec2_project(fv_vec2 v, fv_vec2 onto)
{
    fv_real denom = fv_vec2_len_sqr(onto);
    /* denom is a squared length, so the squared threshold FV_EPS_SQR must be used (P0 fix: this used to be the linear threshold FV_EPS). */
    if (fv_abs(denom) <= FV_EPS_SQR) {
        return fv_vec2_make((fv_real)0, (fv_real)0);
    }
    return fv_vec2_scale(onto, fv_vec2_dot(v, onto) / denom);
}

int fv_vec2_near(fv_vec2 a, fv_vec2 b, fv_real eps)
{
    if (fv_abs(a.x - b.x) > eps) {
        return 0;
    }
    if (fv_abs(a.y - b.y) > eps) {
        return 0;
    }
    return 1;
}

/* ---------------- 3D ---------------- */

fv_vec3 fv_vec3_make(fv_real x, fv_real y, fv_real z)
{
    fv_vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

fv_vec3 fv_vec3_add(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

fv_vec3 fv_vec3_sub(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

fv_vec3 fv_vec3_neg(fv_vec3 v)
{
    return fv_vec3_make(-v.x, -v.y, -v.z);
}

fv_vec3 fv_vec3_scale(fv_vec3 v, fv_real s)
{
    return fv_vec3_make(v.x * s, v.y * s, v.z * s);
}

fv_vec3 fv_vec3_lerp(fv_vec3 a, fv_vec3 b, fv_real t)
{
    return fv_vec3_make(a.x + (b.x - a.x) * t,
                        a.y + (b.y - a.y) * t,
                        a.z + (b.z - a.z) * t);
}

fv_vec3 fv_vec3_min(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x < b.x ? a.x : b.x,
                        a.y < b.y ? a.y : b.y,
                        a.z < b.z ? a.z : b.z);
}

fv_vec3 fv_vec3_max(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.x > b.x ? a.x : b.x,
                        a.y > b.y ? a.y : b.y,
                        a.z > b.z ? a.z : b.z);
}

fv_vec3 fv_vec3_clamp(fv_vec3 v, fv_vec3 lo, fv_vec3 hi)
{
    return fv_vec3_make(fv_clamp1(v.x, lo.x, hi.x),
                        fv_clamp1(v.y, lo.y, hi.y),
                        fv_clamp1(v.z, lo.z, hi.z));
}

fv_vec3 fv_vec3_madd(fv_vec3 a, fv_vec3 b, fv_real s)
{
    return fv_vec3_make(a.x + b.x * s,
                        a.y + b.y * s,
                        a.z + b.z * s);
}

fv_real fv_vec3_dot(fv_vec3 a, fv_vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

fv_real fv_vec3_len_sqr(fv_vec3 v)
{
    return fv_vec3_dot(v, v);
}

fv_real fv_vec3_len(fv_vec3 v)
{
    return fv_sqrt(fv_vec3_len_sqr(v));
}

fv_real fv_vec3_dist_sqr(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_len_sqr(fv_vec3_sub(a, b));
}

fv_real fv_vec3_dist(fv_vec3 a, fv_vec3 b)
{
    return fv_sqrt(fv_vec3_dist_sqr(a, b));
}

fv_vec3 fv_vec3_normalize(fv_vec3 v)
{
    fv_real len = fv_vec3_len(v);
    if (len <= FV_EPS) {
        return fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)0);
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_vec3_scale(v, inv);
    }
}

fv_vec3 fv_vec3_project(fv_vec3 v, fv_vec3 onto)
{
    fv_real denom = fv_vec3_len_sqr(onto);
    /* As in fv_vec2_project: a squared quantity takes the squared threshold. */
    if (fv_abs(denom) <= FV_EPS_SQR) {
        return fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)0);
    }
    return fv_vec3_scale(onto, fv_vec3_dot(v, onto) / denom);
}

fv_vec3 fv_vec3_cross(fv_vec3 a, fv_vec3 b)
{
    return fv_vec3_make(a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x);
}

int fv_vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps)
{
    if (fv_abs(a.x - b.x) > eps) {
        return 0;
    }
    if (fv_abs(a.y - b.y) > eps) {
        return 0;
    }
    if (fv_abs(a.z - b.z) > eps) {
        return 0;
    }
    return 1;
}

/* ---------------- 2D bounding box (for broad-phase) ---------------- */

fv_aabb2 fv_aabb2_merge(fv_aabb2 box, fv_vec2 p)
{
    fv_aabb2 r;
    r.lo = fv_vec2_min(box.lo, p);
    r.hi = fv_vec2_max(box.hi, p);
    return r;
}

int fv_aabb2_overlaps(fv_aabb2 a, fv_aabb2 b)
{
    if (a.lo.x > b.hi.x || b.lo.x > a.hi.x) {
        return 0;
    }
    if (a.lo.y > b.hi.y || b.lo.y > a.hi.y) {
        return 0;
    }
    return 1;
}
