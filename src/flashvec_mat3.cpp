/* The scalar implementation of fv_mat3 -- the 3x3 linear-algebra half of the physics primitives.
 *
 * Same rules as src/flashvec.cpp / src/flashvec_vec4.cpp / src/flashvec_quat.cpp: compiled as C++, C
 * language subset only, no global state, no allocation, #pragma GCC optimize("fp-contract=off") so the
 * mul+add chains are never fused into an FMA.
 *
 * Layout: row-major, m[row * 3 + col] -- the same convention as fv_mat4. The expressions below are the
 * reference the header's SoA batch tail reproduces exactly:
 *   transform_vec3: out.x = (m00*x + m01*y) + m02*z, and the two analogous rows;
 *   mul:            (a*b)[i][j] = (ai0*b0j + ai1*b1j) + ai2*b2j.
 */
#include "flashvec.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

fv_mat3 fv_mat3_from_rows(fv_vec3 r0, fv_vec3 r1, fv_vec3 r2)
{
    fv_mat3 m;
    m.m[0] = r0.x;
    m.m[1] = r0.y;
    m.m[2] = r0.z;
    m.m[3] = r1.x;
    m.m[4] = r1.y;
    m.m[5] = r1.z;
    m.m[6] = r2.x;
    m.m[7] = r2.y;
    m.m[8] = r2.z;
    return m;
}

fv_vec3 fv_mat3_transform_vec3(const fv_mat3 *m, fv_vec3 v)
{
    fv_vec3 r;
    r.x = (m->m[0] * v.x + m->m[1] * v.y) + m->m[2] * v.z;
    r.y = (m->m[3] * v.x + m->m[4] * v.y) + m->m[5] * v.z;
    r.z = (m->m[6] * v.x + m->m[7] * v.y) + m->m[8] * v.z;
    return r;
}

fv_mat3 fv_mat3_mul(const fv_mat3 *a, const fv_mat3 *b)
{
    fv_mat3 r;
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r.m[i * 3 + j] = (a->m[i * 3 + 0] * b->m[0 * 3 + j] +
                              a->m[i * 3 + 1] * b->m[1 * 3 + j]) +
                             a->m[i * 3 + 2] * b->m[2 * 3 + j];
        }
    }
    return r;
}

fv_mat3 fv_mat3_transpose(const fv_mat3 *m)
{
    fv_mat3 r;
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r.m[i * 3 + j] = m->m[j * 3 + i];
        }
    }
    return r;
}

/* Row-major rotation matrix of a unit quaternion:
 *   R = [[1-2(yy+zz), 2(xy-wz),   2(xz+wy)],
 *        [2(xy+wz),   1-2(xx+zz), 2(yz-wx)],
 *        [2(xz-wy),   2(yz+wx),   1-2(xx+yy)]]
 * The product order inside each entry is fixed (2*(...), then the subtraction from 1) so the function
 * has one documented bit pattern; it is the matrix form of the same rotation as fv_quat_rotate_vec3,
 * but the two are different expressions and are compared approximately, not with memcmp. */
fv_mat3 fv_mat3_from_quat(fv_quat q)
{
    fv_real xx = q.x * q.x;
    fv_real yy = q.y * q.y;
    fv_real zz = q.z * q.z;
    fv_real xy = q.x * q.y;
    fv_real xz = q.x * q.z;
    fv_real yz = q.y * q.z;
    fv_real wx = q.w * q.x;
    fv_real wy = q.w * q.y;
    fv_real wz = q.w * q.z;
    fv_mat3 r;
    r.m[0] = (fv_real)1 - (fv_real)2 * (yy + zz);
    r.m[1] = (fv_real)2 * (xy - wz);
    r.m[2] = (fv_real)2 * (xz + wy);
    r.m[3] = (fv_real)2 * (xy + wz);
    r.m[4] = (fv_real)1 - (fv_real)2 * (xx + zz);
    r.m[5] = (fv_real)2 * (yz - wx);
    r.m[6] = (fv_real)2 * (xz - wy);
    r.m[7] = (fv_real)2 * (yz + wx);
    r.m[8] = (fv_real)1 - (fv_real)2 * (xx + yy);
    return r;
}
