/* General 4x4 matrix operations -- the P7 completion of the mat4 surface.
 *
 * Row-major (m[row*4 + col], the same convention as fv_mat3 and fv_mat4_transform_*), general 4x4
 * (not affine-only), compiled as C++ using only the C language subset, no allocation, no global
 * state. The unit carries #pragma GCC optimize("fp-contract=off"), so the mul+add chains are never
 * fused into an FMA and the results are reproducible across translation units.
 *
 * The inverse is the adjugate/determinant closed form written with six 2x2 minors (s0..s5) and six
 * complementary minors (c0..c5). That is the standard MESA/GLU formulation: 16 cofactors from 12
 * minors instead of 16 independent 3x3 determinants. The singularity test is scale-relative,
 * |det| <= FV_EPS * max|entry|^4, so uniformly scaling an invertible matrix does not change the
 * verdict. A non-finite entry is treated as singular (zero matrix, checked 0), matching fv_mat3.
 */
#include "flashvec.h"

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

static fv_real fv4m_abs(fv_real v)
{
    return v < (fv_real)0 ? -v : v;
}

static void fv4m_zero(fv_mat4 *m)
{
    int k;
    for (k = 0; k < 16; ++k) {
        m->m[k] = (fv_real)0;
    }
}

fv_mat4 fv_mat4_mul(const fv_mat4 *a, const fv_mat4 *b)
{
    fv_mat4 r;
    int i;
    int j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            r.m[i * 4 + j] = ((a->m[i * 4 + 0] * b->m[0 * 4 + j] +
                               a->m[i * 4 + 1] * b->m[1 * 4 + j]) +
                              a->m[i * 4 + 2] * b->m[2 * 4 + j]) +
                             a->m[i * 4 + 3] * b->m[3 * 4 + j];
        }
    }
    return r;
}

static int fv4m_inverse_impl(const fv_mat4 *m, fv_mat4 *out)
{
    const fv_real *v = m->m;
    fv_real s0 = v[0] * v[5] - v[4] * v[1];
    fv_real s1 = v[0] * v[6] - v[4] * v[2];
    fv_real s2 = v[0] * v[7] - v[4] * v[3];
    fv_real s3 = v[1] * v[6] - v[5] * v[2];
    fv_real s4 = v[1] * v[7] - v[5] * v[3];
    fv_real s5 = v[2] * v[7] - v[6] * v[3];
    fv_real c5 = v[10] * v[15] - v[14] * v[11];
    fv_real c4 = v[9] * v[15] - v[13] * v[11];
    fv_real c3 = v[9] * v[14] - v[13] * v[10];
    fv_real c2 = v[8] * v[15] - v[12] * v[11];
    fv_real c1 = v[8] * v[14] - v[12] * v[10];
    fv_real c0 = v[8] * v[13] - v[12] * v[9];
    fv_real det = (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2) - s4 * c1 + s5 * c0;
    fv_real scale = (fv_real)0;
    fv_real scale2;
    int k;
    for (k = 0; k < 16; ++k) {
        fv_real a = fv4m_abs(v[k]);
        if (a > scale) {
            scale = a;
        }
    }
    scale2 = scale * scale;
    if (!(scale > (fv_real)0) ||
        !(fv4m_abs(det) > FV_EPS * (scale2 * scale2))) {
        fv4m_zero(out);
        return 0;
    }
    {
        fv_real inv = (fv_real)1 / det;
        fv_mat4 r;
        r.m[0] = (v[5] * c5 - v[6] * c4 + v[7] * c3) * inv;
        r.m[1] = (-v[1] * c5 + v[2] * c4 - v[3] * c3) * inv;
        r.m[2] = (v[13] * s5 - v[14] * s4 + v[15] * s3) * inv;
        r.m[3] = (-v[9] * s5 + v[10] * s4 - v[11] * s3) * inv;
        r.m[4] = (-v[4] * c5 + v[6] * c2 - v[7] * c1) * inv;
        r.m[5] = (v[0] * c5 - v[2] * c2 + v[3] * c1) * inv;
        r.m[6] = (-v[12] * s5 + v[14] * s2 - v[15] * s1) * inv;
        r.m[7] = (v[8] * s5 - v[10] * s2 + v[11] * s1) * inv;
        r.m[8] = (v[4] * c4 - v[5] * c2 + v[7] * c0) * inv;
        r.m[9] = (-v[0] * c4 + v[1] * c2 - v[3] * c0) * inv;
        r.m[10] = (v[12] * s4 - v[13] * s2 + v[15] * s0) * inv;
        r.m[11] = (-v[8] * s4 + v[9] * s2 - v[11] * s0) * inv;
        r.m[12] = (-v[4] * c3 + v[5] * c1 - v[6] * c0) * inv;
        r.m[13] = (v[0] * c3 - v[1] * c1 + v[2] * c0) * inv;
        r.m[14] = (-v[12] * s3 + v[13] * s1 - v[14] * s0) * inv;
        r.m[15] = (v[8] * s3 - v[9] * s1 + v[10] * s0) * inv;
        *out = r;
        return 1;
    }
}

fv_mat4 fv_mat4_inverse(const fv_mat4 *m)
{
    fv_mat4 r;
    (void)fv4m_inverse_impl(m, &r);
    return r;
}

int fv_mat4_inverse_checked(const fv_mat4 *m, fv_mat4 *out)
{
    return fv4m_inverse_impl(m, out);
}
