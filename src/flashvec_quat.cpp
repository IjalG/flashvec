/* The scalar implementation of fv_quat -- the rotation half of the physics primitives.
 *
 * Written the same way as src/flashvec.cpp / src/flashvec_vec4.cpp: compiled as C++, C language subset
 * only, no global state, no allocation, all values by value. It is a separate translation unit for the
 * same two reasons the vec4 file is one: the work is a new type group rather than an edit of the stable
 * 2D/3D baseline, and keeping the sqrtf / -ffp-contract=off unit separate keeps the bit-exactness story
 * reviewable by hand.
 *
 * The expressions are the reference the header's fv_i_s_quat_* helpers and the SoA batch entries mirror
 * word for word (association order included). Do not "simplify" one side alone: the memcmp tests in
 * tests/test_quat.cpp and tests/test_quat_soa.c compare the two paths bit for bit.
 */
#include "flashvec.h"

#include <math.h>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

fv_quat fv_quat_make(fv_real x, fv_real y, fv_real z, fv_real w)
{
    fv_quat q;
    q.x = x;
    q.y = y;
    q.z = z;
    q.w = w;
    return q;
}

fv_quat fv_quat_add(fv_quat a, fv_quat b)
{
    return fv_quat_make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

fv_quat fv_quat_sub(fv_quat a, fv_quat b)
{
    return fv_quat_make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

fv_quat fv_quat_neg(fv_quat q)
{
    return fv_quat_make(-q.x, -q.y, -q.z, -q.w);
}

fv_quat fv_quat_scale(fv_quat q, fv_real s)
{
    return fv_quat_make(q.x * s, q.y * s, q.z * s, q.w * s);
}

/* Left to right: ((x*x' + y*y') + z*z') + w*w'. The same association order the batch normalize uses. */
fv_real fv_quat_dot(fv_quat a, fv_quat b)
{
    return ((a.x * b.x + a.y * b.y) + a.z * b.z) + a.w * b.w;
}

fv_real fv_quat_len_sqr(fv_quat q)
{
    return fv_quat_dot(q, q);
}

fv_real fv_quat_len(fv_quat q)
{
    return sqrtf(fv_quat_len_sqr(q));
}

/* A zero quaternion is not a rotation; the identity keeps the result usable and NaN-free. This is the
 * deliberate difference from fv_vec*_normalize (which returns the zero vector) and is stated in the
 * header contract. */
fv_quat fv_quat_normalize(fv_quat q)
{
    fv_real len = fv_quat_len(q);
    if (len <= FV_EPS) {
        return fv_quat_identity();
    }
    {
        fv_real inv = (fv_real)1 / len;
        return fv_quat_scale(q, inv);
    }
}

fv_quat fv_quat_conjugate(fv_quat q)
{
    return fv_quat_make(-q.x, -q.y, -q.z, q.w);
}

/* Hamilton product a*b: applies b first, then a. Association order fixed as
 * ((w1*x2 + x1*w2) + y1*z2) - z1*y2, and the three analogues -- see the note in include/flashvec.h. */
fv_quat fv_quat_mul(fv_quat a, fv_quat b)
{
    fv_quat r;
    r.x = ((a.w * b.x + a.x * b.w) + a.y * b.z) - a.z * b.y;
    r.y = ((a.w * b.y - a.x * b.z) + a.y * b.w) + a.z * b.x;
    r.z = ((a.w * b.z + a.x * b.y) - a.y * b.x) + a.z * b.w;
    r.w = ((a.w * b.w - a.x * b.x) - a.y * b.y) - a.z * b.z;
    return r;
}

/* v' = q*(v,0)*conjugate(q), evaluated as t = 2*cross(q.xyz,v); d = cross(q.xyz,t); out = (v + q.w*t) + d.
 * q must be unit; this function does not normalize (same precondition as a rotation matrix). */
fv_vec3 fv_quat_rotate_vec3(fv_quat q, fv_vec3 v)
{
    fv_real tx = (fv_real)2 * (q.y * v.z - q.z * v.y);
    fv_real ty = (fv_real)2 * (q.z * v.x - q.x * v.z);
    fv_real tz = (fv_real)2 * (q.x * v.y - q.y * v.x);
    fv_real dx = q.y * tz - q.z * ty;
    fv_real dy = q.z * tx - q.x * tz;
    fv_real dz = q.x * ty - q.y * tx;
    fv_vec3 r;
    r.x = (v.x + q.w * tx) + dx;
    r.y = (v.y + q.w * ty) + dy;
    r.z = (v.z + q.w * tz) + dz;
    return r;
}
