/* In-library definitions of the flashvec batch entry points (after the P1 refactor): there is one implementation only, in include/flashvec.h.
 *
 * This file does three things:
 *   1) use FV_BUILD to turn the header into its "declarations only" form: the public names are extern in the header,
 *      while the implementation section still compiles (that section is controlled by !FV_USE_EXTERN, and FV_BUILD does not exclude it);
 *   2) give each batch entry point an extern symbol whose one-line body forwards to the same-named fv_i_*;
 *   3) define the four per-value integrator primitives (fv_vec2_addmul / fv_vec3_addmul / fv_vec4_addmul / fv_vec2_maddv).
 *      Their header inline bodies are part of P0; this file keeps the same-named extern symbols for FV_USE_EXTERN / ABI.
 *      They live at the end of this file because this is the in-library unit that the integrator group touches.
 *
 * Why these symbols remain: binary compatibility. Old consumers (already-built .o/.so files calling the extern names) keep linking;
 * new consumers get the header's static inline by default and never reference these symbols at link time.
 *
 * Bit-exactness, the mask / index / tail-length semantics and the conditions under which intrinsics appear are all written in
 * include/flashvec.h's "batch entry points (P1)" contract and file header; no second copy is made here, so the two cannot drift apart.
 */

#define FV_BUILD 1

#include "flashvec.h"

void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4(a, b, o, n);
}

void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4(a, b, o, n);
}

void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                      fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4(a, s, o, n);
}

void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real t, fv_vec2 *__restrict o, size_t n)
{
    fv_i_lerp_x4(a, b, t, o, n);
}

void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real s, fv_vec2 *__restrict o, size_t n)
{
    fv_i_madd_x4(a, b, s, o, n);
}

void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                         fv_vec2 *__restrict o, const uint8_t *__restrict mask, size_t n)
{
    fv_i_add_x4_mask(a, b, o, mask, n);
}

void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                        fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m)
{
    fv_i_add_x4_idx(a, b, o, idx, m);
}

void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                          fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m)
{
    fv_i_scale_x4_idx(a, s, o, idx, m);
}

void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict bx, const float *__restrict by,
                     float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa(ax, ay, bx, by, ox, oy, n);
}

void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay, fv_real s,
                       float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_scale_soa(ax, ay, s, ox, oy, n);
}

float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict bx, const float *__restrict by, size_t n)
{
    return fv_i_dot_soa(ax, ay, bx, by, n);
}

void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                     float *__restrict o, size_t n)
{
    fv_i_len_soa(ax, ay, o, n);
}

void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_normalize_soa(ax, ay, ox, oy, n);
}

/* ---------------- memory-traffic entry points (P2): forwarding only, as above ----------------
 * The NT store's alignment peel, its sfence and the mask semantics all live in the include/flashvec.h kernels;
 * the extern symbols here let the FV_USE_EXTERN form (and already-built old consumers) link to that same implementation. */
void fv_vec2_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4_stream(a, b, o, n);
}

void fv_vec2_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4_stream(a, b, o, n);
}

void fv_vec2_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                             fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4_stream(a, s, o, n);
}

void fv_vec2_lerp_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                            fv_real t, fv_vec2 *__restrict o, size_t n)
{
    fv_i_lerp_x4_stream(a, b, t, o, n);
}

void fv_vec2_madd_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                            fv_real s, fv_vec2 *__restrict o, size_t n)
{
    fv_i_madd_x4_stream(a, b, s, o, n);
}

void fv_vec2_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                            const float *__restrict bx, const float *__restrict by,
                            float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa_stream(ax, ay, bx, by, ox, oy, n);
}

void fv_vec2_scale_soa_stream(const float *__restrict ax, const float *__restrict ay, fv_real s,
                              float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_scale_soa_stream(ax, ay, s, ox, oy, n);
}

void fv_vec2_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_normalize_soa_stream(ax, ay, ox, oy, n);
}

/* ---------------- fused pipeline (out = normalize(a * s + b) in one pass): forwarding only, as above ----------------
 * The evaluation order, the zero-vector mask and the NT store convention live in include/flashvec.h; this pair exists so the
 * FV_USE_EXTERN form and already-built consumers reach the same implementation. */
void fv_vec2_madd_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict bx, const float *__restrict by,
                                fv_real s, float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_madd_norm_soa(ax, ay, bx, by, s, ox, oy, n);
}

void fv_vec2_madd_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict bx, const float *__restrict by,
                                       fv_real s, float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_madd_norm_soa_stream(ax, ay, bx, by, s, ox, oy, n);
}

void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                 const uint8_t *__restrict mask, size_t n)
{
    fv_i_add_x4_mask_inplace(a, b, mask, n);
}

/* ---------------- 3D SoA / 4D _x2 and _soa (P3): forwarding only, as above ----------------
 * 3D has SoA only (fv_vec3 stays 12 B with no padding added) and 4D AoS uses _x2 (2 x 16 B elements =
 * one ymm). The kernels, tails, bit-exactness and sectioning are all in include/flashvec.h; no copy of the notes is made here. */
void fv_vec3_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az,
                     const float *__restrict bx, const float *__restrict by,
                     const float *__restrict bz,
                     float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n)
{
    fv_i_v3_add_soa(ax, ay, az, bx, by, bz, ox, oy, oz, n);
}

void fv_vec3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                       const float *__restrict az, fv_real s,
                       float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n)
{
    fv_i_v3_scale_soa(ax, ay, az, s, ox, oy, oz, n);
}

void fv_vec3_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                            const float *__restrict az,
                            const float *__restrict bx, const float *__restrict by,
                            const float *__restrict bz,
                            float *__restrict ox, float *__restrict oy,
                            float *__restrict oz, size_t n)
{
    fv_i_v3_add_soa_stream(ax, ay, az, bx, by, bz, ox, oy, oz, n);
}

void fv_vec3_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                              const float *__restrict az, fv_real s,
                              float *__restrict ox, float *__restrict oy,
                              float *__restrict oz, size_t n)
{
    fv_i_v3_scale_soa_stream(ax, ay, az, s, ox, oy, oz, n);
}

float fv_vec3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, size_t n)
{
    return fv_i_v3_dot_soa(ax, ay, az, bx, by, bz, n);
}

void fv_vec3_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, float *__restrict o, size_t n)
{
    fv_i_v3_len_soa(ax, ay, az, o, n);
}

void fv_vec3_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az,
                           float *__restrict ox, float *__restrict oy, float *__restrict oz,
                           size_t n)
{
    fv_i_v3_normalize_soa(ax, ay, az, ox, oy, oz, n);
}

void fv_vec3_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az,
                                  float *__restrict ox, float *__restrict oy,
                                  float *__restrict oz, size_t n)
{
    fv_i_v3_normalize_soa_stream(ax, ay, az, ox, oy, oz, n);
}

void fv_vec4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_add_x2(a, b, o, n);
}

void fv_vec4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_sub_x2(a, b, o, n);
}

void fv_vec4_scale_x2(const fv_vec4 *__restrict a, fv_real s, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_scale_x2(a, s, o, n);
}

void fv_vec4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real t,
                     fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_lerp_x2(a, b, t, o, n);
}

void fv_vec4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real s,
                     fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_madd_x2(a, b, s, o, n);
}

void fv_vec4_add_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                           fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_add_x2_stream(a, b, o, n);
}

void fv_vec4_sub_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                           fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_sub_x2_stream(a, b, o, n);
}

void fv_vec4_scale_x2_stream(const fv_vec4 *__restrict a, fv_real s,
                             fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_scale_x2_stream(a, s, o, n);
}

void fv_vec4_lerp_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real t,
                            fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_lerp_x2_stream(a, b, t, o, n);
}

void fv_vec4_madd_x2_stream(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real s,
                            fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_madd_x2_stream(a, b, s, o, n);
}

float fv_vec4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az, const float *__restrict aw,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, const float *__restrict bw, size_t n)
{
    return fv_i_v4_dot_soa(ax, ay, az, aw, bx, by, bz, bw, n);
}

void fv_vec4_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, const float *__restrict aw,
                     float *__restrict o, size_t n)
{
    fv_i_v4_len_soa(ax, ay, az, aw, o, n);
}

void fv_vec4_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az, const float *__restrict aw,
                           float *__restrict ox, float *__restrict oy,
                           float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_v4_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
}

void fv_vec4_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az, const float *__restrict aw,
                                  float *__restrict ox, float *__restrict oy,
                                  float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_v4_normalize_soa_stream(ax, ay, az, aw, ox, oy, oz, ow, n);
}

void fv_mat4_transform_points(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                              fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_points(m, in, out, n);
}

void fv_mat4_transform_dirs(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                            fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_dirs(m, in, out, n);
}

void fv_mat4_transform_points_stream(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                                     fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_points_stream(m, in, out, n);
}

void fv_mat4_transform_dirs_stream(const fv_mat4 *__restrict m, const fv_vec4 *__restrict in,
                                   fv_vec4 *__restrict out, size_t n)
{
    fv_i_mat4_tf_dirs_stream(m, in, out, n);
}

void fv_vec2_addmul_x4(fv_vec2 *__restrict x, const fv_vec2 *__restrict v, fv_real s, size_t n)
{
    fv_i_addmul_x4(x, v, s, n);
}

void fv_vec2_addmul_soa(float *__restrict x, float *__restrict y,
                        const float *__restrict vx, const float *__restrict vy, fv_real s,
                        size_t n)
{
    fv_i_addmul_soa(x, y, vx, vy, s, n);
}

void fv_vec3_addmul_soa(float *__restrict x, float *__restrict y, float *__restrict z,
                        const float *__restrict vx, const float *__restrict vy,
                        const float *__restrict vz, fv_real s, size_t n)
{
    fv_i_addmul_v3_soa(x, y, z, vx, vy, vz, s, n);
}

void fv_vec4_addmul_x2(fv_vec4 *__restrict x, const fv_vec4 *__restrict v, fv_real s, size_t n)
{
    fv_i_v4_addmul_x2(x, v, s, n);
}

void fv_vec2_maddv_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                      const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n)
{
    fv_i_maddv_x4(a, b, c, o, n);
}

void fv_vec2_maddv_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                             const fv_vec2 *__restrict c, fv_vec2 *__restrict o, size_t n)
{
    fv_i_maddv_x4_stream(a, b, c, o, n);
}

/* The per-value half of the integrator primitives (fv_vec2_addmul / fv_vec3_addmul / fv_vec4_addmul / fv_vec2_maddv).
 * They are scalar API and therefore in-library symbols, exactly like fv_vec2_add / fv_vec2_madd; they are defined in this
 * unit because it is the one in-library unit this group touches, and they carry the same fp-contract=off pragma the other
 * scalar definitions carry -- a consumer built with -mfma -ffp-contract=fast still links the unfused bit pattern. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif

fv_vec2 fv_vec2_addmul(fv_vec2 x, fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    return r;
}

fv_vec3 fv_vec3_addmul(fv_vec3 x, fv_vec3 v, fv_real s)
{
    fv_vec3 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    r.z = x.z + v.z * s;
    return r;
}

fv_vec4 fv_vec4_addmul(fv_vec4 x, fv_vec4 v, fv_real s)
{
    fv_vec4 r;
    r.x = x.x + v.x * s;
    r.y = x.y + v.y * s;
    r.z = x.z + v.z * s;
    r.w = x.w + v.w * s;
    return r;
}

/* a*b + c, product first, then the sum -- the third operand is a vector here, which is the whole point (fv_vec2_madd's
 * third operand is a scalar). */
fv_vec2 fv_vec2_maddv(fv_vec2 a, fv_vec2 b, fv_vec2 c)
{
    fv_vec2 r;
    r.x = a.x * b.x + c.x;
    r.y = a.y * b.y + c.y;
    return r;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
