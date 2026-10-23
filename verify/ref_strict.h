/* verify/ref_strict.h —— 独立标量参考（验收者 F 自建，不依赖被验证实现的标量路径）
 * 公式与 flashvec 标量路径逐字对应；本 TU 单独编译且带 -ffp-contract=off，
 * 因此不会被编译器收缩成 FMA，可与批量路径做位级比较。 */
#ifndef FV_VERIFY_REF_STRICT_H
#define FV_VERIFY_REF_STRICT_H

#include "flashvec.h"

#ifdef __cplusplus
extern "C" {
#endif

fv_vec2 fvs_add(fv_vec2 a, fv_vec2 b);
fv_vec2 fvs_sub(fv_vec2 a, fv_vec2 b);
fv_vec2 fvs_scale(fv_vec2 v, fv_real s);
fv_vec2 fvs_lerp(fv_vec2 a, fv_vec2 b, fv_real t);
fv_vec2 fvs_madd(fv_vec2 a, fv_vec2 b, fv_real s);
fv_real fvs_dot(fv_vec2 a, fv_vec2 b);
fv_real fvs_len(fv_vec2 v);
fv_vec2 fvs_normalize(fv_vec2 v);

#ifdef __cplusplus
}
#endif
#endif /* FV_VERIFY_REF_STRICT_H */
