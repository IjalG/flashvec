/* bench/ab/ab_scalar_loop.cpp -- the "before the refactor" anchor: element-wise scalar
 * API calls.
 *
 * This TU only includes the header; the scalar functions are defined elsewhere, so every
 * call is an ordinary cross-TU call. That is both why ab_scalar_lib.cpp exists and the one
 * property this fixture must defend: the anchor has to be dumber than the batch entry
 * point, otherwise what gets measured is the compiler's inlining.
 *
 * Coverage is only the 6 ops the small-N regime needs (add/scale/lerp/dot/len/normalize):
 * the large-N regime compares two batch implementations (A/B); the scalar anchor is
 * verify/perf.cpp's job.
 */
#include "snapshot/include/flashvec.h"

#include "ab_abi.h"

#define AB_S_ ((fv_real)-0.25)
#define AB_T_ ((fv_real)0.375)

extern "C" unsigned long long ab_s_abi(void)
{
    return ((unsigned long long)sizeof(fv_vec2) << 32) | (unsigned)sizeof(fv_real);
}

/* the scalar anchor has no NT store: every write is a plain store. */
extern "C" unsigned ab_s_opflags(int op)
{
    (void)op;
    return 0u;
}

extern "C" void ab_s_loop(int op, size_t n, long iters, const AbArgs *p, float *dot)
{
    const fv_vec2 *a = (const fv_vec2 *)p->a;
    const fv_vec2 *b = (const fv_vec2 *)p->b;
    fv_vec2 *o = (fv_vec2 *)p->o;

    switch (op) {
    case AB_OP_ADD_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_add(a[i], b[i]);
        break;
    case AB_OP_SCALE_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_scale(a[i], AB_S_);
        break;
    case AB_OP_LERP_X4:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) o[i] = fv_vec2_lerp(a[i], b[i], AB_T_);
        break;
    case AB_OP_DOT_SOA: {
        fv_real acc = (fv_real)0;
        for (long k = 0; k < iters; ++k) {
            fv_real s = (fv_real)0;
            for (size_t i = 0; i < n; ++i)
                s = s + fv_vec2_dot(fv_vec2_make(p->ax[i], p->ay[i]),
                                    fv_vec2_make(p->bx[i], p->by[i]));
            acc += s;
        }
        *dot = (float)acc;
        break;
    }
    case AB_OP_LEN_SOA:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) p->lenv[i] = fv_vec2_len(fv_vec2_make(p->ax[i], p->ay[i]));
        break;
    case AB_OP_NORM_SOA:
        for (long k = 0; k < iters; ++k)
            for (size_t i = 0; i < n; ++i) {
                const fv_vec2 r = fv_vec2_normalize(fv_vec2_make(p->ax[i], p->ay[i]));
                p->ox[i] = r.x;
                p->oy[i] = r.y;
            }
        break;
    default:
        break;
    }
}
