/* flashvec 批量入口的库内定义（P1 重构后）：实现只有一份，在 include/flashvec.h。
 *
 * 本文件只做两件事：
 *   1) 用 FV_BUILD 把头文件切成「只有声明」的形态：公开名在头文件里是 extern，
 *      实现小节仍然编译（那一节由 !FV_USE_EXTERN 控制，FV_BUILD 不排除它）；
 *   2) 为 33 个批量入口各给一个 extern 符号，函数体一行转发到同名 fv_i_*。
 *
 * 为什么还留这些符号：二进制兼容。旧消费方（已经编译好、调用 extern 的 .o/.so）继续能链；
 * 新消费方默认拿头文件里的 static inline，链接时不会引用这些符号。
 *
 * 位级一致性、掩码/索引/尾长语义、intrinsics 的出现条件都写在 include/flashvec.h 的
 * 「批量入口（P1）」契约注释与文件头里；这里不复制第二份说明，避免两处说法漂移。
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

/* ---------------- 内存流量优化入口（P2）：同样只做转发 ----------------
 * NT store 的对齐 peel、sfence、掩码语义全在 include/flashvec.h 的内核里；
 * 这里给 extern 符号是为了让 FV_USE_EXTERN 形态（以及已编译好的旧消费方）能链到同一份实现。 */
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

void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                 const uint8_t *__restrict mask, size_t n)
{
    fv_i_add_x4_mask_inplace(a, b, mask, n);
}

/* ---------------- 三维 SoA / 四维 _x2、_soa（P3）：同样只做转发 ----------------
 * 三维只有 SoA 形态（fv_vec3 保持 12 B、不加 padding），四维 AoS 用 _x2（2 个 16 B 元素 =
 * 一个 ymm）。内核、尾部、位级一致性、分段归档全在 include/flashvec.h 里，这里不复制说明。 */
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
