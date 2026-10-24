/* bench/ab/ab_abi.h -- the only ABI between the variant TUs and the driver
 *
 * Why this layer exists: include/flashvec.h is being refactored into a "header-inline"
 * form (static inline by default; -DFV_USE_EXTERN switches back to library definitions).
 * The batch entry points have **different linkage** under the two forms (TU-private static
 * when inlined, external symbols in library form), so the driver cannot call them through
 * "declaration + function pointer".
 *
 * The fix: put the timing loop inside each variant's own TU (see ab_variant.cpp) and let it
 * export one whole "sample" at a time. Hence:
 *   - library form -> a plain cross-function call inside the TU (exactly what a consumer
 *                     linking a static library gets);
 *   - inline form  -> inlined inside the TU (exactly what the refactor wants).
 * Both forms are measurable, and neither requires the driver to know the symbol names.
 *
 * This header **includes no library types**: every pointer is void * / float * / uint32_t *,
 * so A and B cannot read garbage just because the two headers disagree on struct layout.
 * Layout consistency has its own runtime check (the sizeof fingerprint exported under
 * AB_ABI_NAME).
 */
#ifndef AB_ABI_H
#define AB_ABI_H

#include <stddef.h>
#include <stdint.h>

#include "ab_ops.h"

struct AbArgs {
    const void *a;   /* const fv_vec2* */
    const void *b;   /* const fv_vec2* */
    void *o;         /* fv_vec2*      */
    const uint8_t *mask;
    const uint32_t *idx;
    const float *ax, *ay, *bx, *by;
    float *ox, *oy, *lenv;
    float *ox2, *oy2;         /* collision: second output set (reflected velocity) */
    float *px, *py, *qx, *qy; /* pipeline scratch (collision normals/velocity; chain3 t/u) */
    /* Component arrays for wide elements (vec3 12 B / vec4 16 B); non-null only when
     * such an op is selected. oz/ow alias vx/vy (the two halves of base2) and are never
     * used at the same time. */
    const float *az, *bz, *aw, *bw;
    float *oz, *ow;
    /* Matrix transforms (16 B AoS elements): mat is 16 floats, row-major, last row
     * (0,0,0,1); vin/vout are AoS arrays of fv_vec4 (one element = 4 floats). */
    const float *mat;
    const void *vin;
    void *vout;
};

/* op ids: derived from the ab_ops.h registry (single source of truth). The variant TUs
 * and the driver share the same set. */
enum {
#define AB_ABI_ENUM_ITEM(ID, ...) AB_OP_##ID,
    AB_LIB_OPS(AB_ABI_ENUM_ITEM)
    AB_LIB_OPS_EXTRA(AB_ABI_ENUM_ITEM)
#undef AB_ABI_ENUM_ITEM
    AB_OP_COUNT
};

/* The variant side tells the driver how this op's stores are really done.
 * The driver uses it for RFO accounting of the output array (plain store = 2x, NT store = 1x).
 * This is **runtime** information because the same ab_variant.cpp is compiled twice (the two
 * headers differ in capability): e.g. when the baseline header has no _stream, the thunk
 * degrades to a plain store and this must report 0. */
#define AB_SIDE_NT 0x1u

#ifdef __cplusplus
extern "C" {
#endif

/* each variant exports three symbols; the names come from -DAB_LOOP_NAME / -DAB_ABI_NAME /
 * -DAB_FLAGS_NAME */
typedef void (*AbLoopFn)(int op, size_t n, long iters, const AbArgs *p, float *dot);

void ab_a_loop(int op, size_t n, long iters, const AbArgs *p, float *dot);
void ab_b_loop(int op, size_t n, long iters, const AbArgs *p, float *dot);
void ab_s_loop(int op, size_t n, long iters, const AbArgs *p, float *dot);

unsigned long long ab_a_abi(void);
unsigned long long ab_b_abi(void);
unsigned long long ab_s_abi(void);

unsigned ab_a_opflags(int op);
unsigned ab_b_opflags(int op);
unsigned ab_s_opflags(int op);

#ifdef __cplusplus
}
#endif

#endif /* AB_ABI_H */
