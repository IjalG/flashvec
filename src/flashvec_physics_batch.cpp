/* In-library definitions of the quaternion / mat3 SoA batch entry points (P6).
 *
 * This file follows src/flashvec_batch.cpp exactly: FV_BUILD turns the header into its "declarations
 * only" form, the implementation section is still compiled (it is controlled by !FV_USE_EXTERN), and
 * every public symbol below is one line forwarding to the same-named fv_i_* helper. There is no second
 * copy of any kernel here, so the two forms cannot drift apart.
 *
 * Why a separate translation unit rather than extending src/flashvec_batch.cpp: the P6 group is a new
 * type family, and keeping it in its own unit keeps the stable P1..P5 file reviewable by hand (the same
 * reason src/flashvec_vec4.cpp exists). The symbols exist for the FV_USE_EXTERN form and for already-
 * built callers; a new consumer gets the header's static inline definitions and never references them.
 */
#define FV_BUILD 1

#include "flashvec.h"

void fv_quat_mul_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, const float *__restrict aw,
                     const float *__restrict bx, const float *__restrict by,
                     const float *__restrict bz, const float *__restrict bw,
                     float *__restrict ox, float *__restrict oy,
                     float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_quat_mul_soa(ax, ay, az, aw, bx, by, bz, bw, ox, oy, oz, ow, n);
}

void fv_quat_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az, const float *__restrict aw,
                           float *__restrict ox, float *__restrict oy,
                           float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_quat_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
}

void fv_quat_rotate_vec3_soa(const float *__restrict qx, const float *__restrict qy,
                             const float *__restrict qz, const float *__restrict qw,
                             const float *__restrict vx, const float *__restrict vy,
                             const float *__restrict vz,
                             float *__restrict ox, float *__restrict oy,
                             float *__restrict oz, size_t n)
{
    fv_i_quat_rotate_vec3_soa(qx, qy, qz, qw, vx, vy, vz, ox, oy, oz, n);
}

void fv_mat3_transform_vec3_soa(const float *__restrict m00, const float *__restrict m01,
                                const float *__restrict m02, const float *__restrict m10,
                                const float *__restrict m11, const float *__restrict m12,
                                const float *__restrict m20, const float *__restrict m21,
                                const float *__restrict m22,
                                const float *__restrict vx, const float *__restrict vy,
                                const float *__restrict vz,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, size_t n)
{
    fv_i_mat3_transform_vec3_soa(m00, m01, m02, m10, m11, m12, m20, m21, m22,
                                 vx, vy, vz, ox, oy, oz, n);
}

void fv_mat3_mul_soa(const float *__restrict a00,
                     const float *__restrict a01,
                     const float *__restrict a02,
                     const float *__restrict a10,
                     const float *__restrict a11,
                     const float *__restrict a12,
                     const float *__restrict a20,
                     const float *__restrict a21,
                     const float *__restrict a22,
                     const float *__restrict b00,
                     const float *__restrict b01,
                     const float *__restrict b02,
                     const float *__restrict b10,
                     const float *__restrict b11,
                     const float *__restrict b12,
                     const float *__restrict b20,
                     const float *__restrict b21,
                     const float *__restrict b22,
                     float *__restrict o00,
                     float *__restrict o01,
                     float *__restrict o02,
                     float *__restrict o10,
                     float *__restrict o11,
                     float *__restrict o12,
                     float *__restrict o20,
                     float *__restrict o21,
                     float *__restrict o22,
                     size_t n)
{
    fv_i_mat3_mul_soa(a00, a01, a02, a10, a11, a12, a20, a21, a22, b00, b01, b02, b10, b11, b12, b20, b21, b22, o00, o01, o02, o10, o11, o12, o20, o21, o22, n);
}

void fv_mat3_inverse_soa(const float *__restrict m00,
                         const float *__restrict m01,
                         const float *__restrict m02,
                         const float *__restrict m10,
                         const float *__restrict m11,
                         const float *__restrict m12,
                         const float *__restrict m20,
                         const float *__restrict m21,
                         const float *__restrict m22,
                         float *__restrict o00,
                         float *__restrict o01,
                         float *__restrict o02,
                         float *__restrict o10,
                         float *__restrict o11,
                         float *__restrict o12,
                         float *__restrict o20,
                         float *__restrict o21,
                         float *__restrict o22,
                         size_t n)
{
    fv_i_mat3_inverse_soa(m00, m01, m02, m10, m11, m12, m20, m21, m22, o00, o01, o02, o10, o11, o12, o20, o21, o22, n);
}
