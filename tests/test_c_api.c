/* 用 C 编译器编译本文件，验证头文件的纯 C 兼容性。 */
#include "flashvec.h"

#include <stdio.h>

int main(void)
{
    fv_vec2 a = fv_vec2_make(3, 4);
    fv_vec3 b = fv_vec3_make(1, 2, 3);

    if (!fv_vec2_near(a, fv_vec2_make(3, 4), 0)) {
        return 1;
    }
    if (fv_vec2_len(a) < 4.999f || fv_vec2_len(a) > 5.001f) {
        return 1;
    }
    if (fv_vec3_dot(b, b) != 14) {
        return 1;
    }

    printf("C API ok: |(3,4)| = %f, |(1,2,3)|^2 = %f\n",
           (double)fv_vec2_len(a),
           (double)fv_vec3_len_sqr(b));
    return 0;
}
