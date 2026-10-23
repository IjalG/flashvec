/* verify/ref_strict.cpp —— 独立标量参考实现。编译：-O2 -ffp-contract=off -fno-tree-vectorize */
#include "ref_strict.h"

#include <math.h>

extern "C" {

fv_vec2 fvs_add(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    return r;
}

fv_vec2 fvs_sub(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    return r;
}

fv_vec2 fvs_scale(fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = v.x * s;
    r.y = v.y * s;
    return r;
}

fv_vec2 fvs_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    fv_vec2 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    return r;
}

fv_vec2 fvs_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    fv_vec2 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    return r;
}

fv_real fvs_dot(fv_vec2 a, fv_vec2 b)
{
    return a.x * b.x + a.y * b.y;
}

fv_real fvs_len(fv_vec2 v)
{
    return sqrtf(v.x * v.x + v.y * v.y);
}

fv_vec2 fvs_normalize(fv_vec2 v)
{
    fv_real len = sqrtf(v.x * v.x + v.y * v.y);
    fv_vec2 r;
    if (len <= FV_EPS) {
        r.x = (fv_real)0;
        r.y = (fv_real)0;
        return r;
    }
    {
        fv_real inv = (fv_real)1 / len;
        r.x = v.x * inv;
        r.y = v.y * inv;
    }
    return r;
}

} /* extern "C" */
