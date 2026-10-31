/* In-library definitions of the _fast approximate SoA normalize entries (P7).
 *
 * Same shape as src/flashvec_physics_batch.cpp: FV_BUILD turns the header into its declaration form,
 * the implementation section still compiles, and every public symbol forwards to the same-named fv_i_*
 * helper. No kernel is copied here. The symbols exist for the FV_USE_EXTERN form and already-built
 * callers; new consumers get the header's static inline definitions.
 */
#define FV_BUILD 1

#include "flashvec.h"

void fv_vec2_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_vec2_normalize_soa_fast(ax, ay, ox, oy, n);
}

void fv_vec3_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, size_t n)
{
    fv_i_vec3_normalize_soa_fast(ax, ay, az, ox, oy, oz, n);
}

void fv_vec4_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az, const float *__restrict aw,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_vec4_normalize_soa_fast(ax, ay, az, aw, ox, oy, oz, ow, n);
}

void fv_quat_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict az, const float *__restrict aw,
                                float *__restrict ox, float *__restrict oy,
                                float *__restrict oz, float *__restrict ow, size_t n)
{
    fv_i_quat_normalize_soa_fast(ax, ay, az, aw, ox, oy, oz, ow, n);
}

void fv_vec2_madd_normalize_soa_fast(const float *__restrict ax, const float *__restrict ay,
                                     const float *__restrict bx, const float *__restrict by,
                                     fv_real s, float *__restrict ox, float *__restrict oy,
                                     size_t n)
{
    fv_i_vec2_madd_normalize_soa_fast(ax, ay, bx, by, s, ox, oy, n);
}
