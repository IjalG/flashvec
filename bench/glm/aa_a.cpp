#include "flashvec.h"
void aa_add4_a(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o, size_t n)
{
    fv_vec4_add_x2(a, b, o, n);
}
