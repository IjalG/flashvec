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

/* General inverse: conjugate(q) / |q|^2. The threshold is on the squared length, so it is
 * dimensionally consistent with FV_EPS (length <= FV_EPS <=> len_sqr <= FV_EPS_SQR). A degenerate
 * input writes the zero quaternion and returns 0; the value form returns that zero. */
int fv_quat_inverse_checked(fv_quat q, fv_quat *out)
{
    fv_real l2 = fv_quat_len_sqr(q);
    if (!(l2 > FV_EPS_SQR)) {
        *out = fv_quat_make((fv_real)0, (fv_real)0, (fv_real)0, (fv_real)0);
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / l2;
        fv_quat c = fv_quat_conjugate(q);
        *out = fv_quat_scale(c, inv);
        return 1;
    }
}

fv_quat fv_quat_inverse(fv_quat q)
{
    fv_quat r;
    (void)fv_quat_inverse_checked(q, &r);
    return r;
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

/* Rotation matrix -> unit quaternion, Shepperd's method: pick the largest of the four diagonal
 * expressions so the divisor is the largest one (the naive trace-only form loses precision when
 * trace is near -1). The matrix is assumed to be a rotation (orthonormal, det +1); the result is
 * normalized. For a non-rotation this still returns a unit quaternion, but not the matrix's rotation,
 * and the by-value API has no error channel to report that. */
fv_quat fv_quat_from_mat3(const fv_mat3 *m)
{
    fv_real trace = m->m[0] + m->m[4] + m->m[8];
    fv_quat q;
    if (trace > (fv_real)0) {
        fv_real s = sqrtf(trace + (fv_real)1) * (fv_real)2;
        q.w = (fv_real)0.25 * s;
        q.x = (m->m[7] - m->m[5]) / s;
        q.y = (m->m[2] - m->m[6]) / s;
        q.z = (m->m[3] - m->m[1]) / s;
    } else if (m->m[0] > m->m[4] && m->m[0] > m->m[8]) {
        fv_real s = sqrtf((fv_real)1 + m->m[0] - m->m[4] - m->m[8]) * (fv_real)2;
        q.w = (m->m[7] - m->m[5]) / s;
        q.x = (fv_real)0.25 * s;
        q.y = (m->m[1] + m->m[3]) / s;
        q.z = (m->m[2] + m->m[6]) / s;
    } else if (m->m[4] > m->m[8]) {
        fv_real s = sqrtf((fv_real)1 + m->m[4] - m->m[0] - m->m[8]) * (fv_real)2;
        q.w = (m->m[2] - m->m[6]) / s;
        q.x = (m->m[1] + m->m[3]) / s;
        q.y = (fv_real)0.25 * s;
        q.z = (m->m[5] + m->m[7]) / s;
    } else {
        fv_real s = sqrtf((fv_real)1 + m->m[8] - m->m[0] - m->m[4]) * (fv_real)2;
        q.w = (m->m[3] - m->m[1]) / s;
        q.x = (m->m[2] + m->m[6]) / s;
        q.y = (m->m[5] + m->m[7]) / s;
        q.z = (fv_real)0.25 * s;
    }
    return fv_quat_normalize(q);
}

/* Shortest-arc normalized lerp: negating the second quaternion when dot < 0 selects the same rotation
 * with the closer sign, then the components are interpolated and normalized. */
fv_quat fv_quat_nlerp(fv_quat a, fv_quat b, fv_real t)
{
    if (fv_quat_dot(a, b) < (fv_real)0) {
        b = fv_quat_neg(b);
    }
    {
        fv_quat r;
        r.x = a.x + (b.x - a.x) * t;
        r.y = a.y + (b.y - a.y) * t;
        r.z = a.z + (b.z - a.z) * t;
        r.w = a.w + (b.w - a.w) * t;
        return fv_quat_normalize(r);
    }
}

/* Shortest-arc slerp. The small-angle threshold (dot > 0.9995) keeps sin(theta) away from underflow and
 * hands the near-identical case to nlerp, which is stable there. */
fv_quat fv_quat_slerp(fv_quat a, fv_quat b, fv_real t)
{
    fv_real d = fv_quat_dot(a, b);
    if (d < (fv_real)0) {
        b = fv_quat_neg(b);
        d = -d;
    }
    if (d > (fv_real)0.9995 || !(d < (fv_real)1)) {
        return fv_quat_nlerp(a, b, t);
    }
    {
        fv_real theta = acosf(d);
        fv_real st = sinf(theta);
        if (st == (fv_real)0) {
            return fv_quat_nlerp(a, b, t);
        }
        {
            fv_real w1 = sinf((fv_real)(1 - t) * theta) / st;
            fv_real w2 = sinf(t * theta) / st;
            fv_quat r = fv_quat_add(fv_quat_scale(a, w1), fv_quat_scale(b, w2));
            return fv_quat_normalize(r);
        }
    }
}

/* Axis-angle: normalize the axis first (a zero axis has no direction and returns the identity), then
 * the half-angle sine/cosine. */
fv_quat fv_quat_from_axis_angle(fv_vec3 axis, fv_real angle)
{
    fv_real len = fv_vec3_len(axis);
    if (len <= FV_EPS) {
        return fv_quat_identity();
    }
    {
        fv_real inv = (fv_real)1 / len;
        fv_real half = angle * (fv_real)0.5;
        fv_real s = sinf(half);
        fv_real c = cosf(half);
        return fv_quat_make(axis.x * inv * s, axis.y * inv * s, axis.z * inv * s, c);
    }
}

/* Euler XYZ, applied X first, then Y, then Z: q = qz * qy * qx. */
fv_quat fv_quat_from_euler_xyz(fv_real x, fv_real y, fv_real z)
{
    fv_quat qx = fv_quat_from_axis_angle(fv_vec3_make((fv_real)1, (fv_real)0, (fv_real)0), x);
    fv_quat qy = fv_quat_from_axis_angle(fv_vec3_make((fv_real)0, (fv_real)1, (fv_real)0), y);
    fv_quat qz = fv_quat_from_axis_angle(fv_vec3_make((fv_real)0, (fv_real)0, (fv_real)1), z);
    return fv_quat_mul(qz, fv_quat_mul(qy, qx));
}
