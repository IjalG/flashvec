/* 消费者 TU 的内联效果探针 —— 验收 B1/B2 的证据源。只编译成 .o，不运行；
 * 断言由 make test-inline-effects 用 objdump + nm 做（两类入口分开断言）。
 *
 * n 取编译期常量：
 *   n = 3（< 4）：_x4 的 bulk 为 0，只剩尾巴 —— 用来证明尾巴循环被**完全内联**、没有 per-element call；
 *   n = 8        ：bulk = 8，批量部分调用一次带 target("avx2") 的内核（每批一次，可接受）。
 *
 * 编译档故意用消费方最常见的选择：-O2，**不加** -ffp-contract=off。
 */
#include "flashvec.h"

/* 探针符号在两个语言下都必须同名（C++ 里否则会被 mangle，objdump --disassemble 找不到）。 */
#ifdef __cplusplus
extern "C" {
#endif

/* ---- 第一段：无 mul+add 收缩风险，期望在消费者里完整内联 ---- */
void probe_add_x4_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_add_x4(a, b, o, 3);
}

void probe_add_x4_8(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_add_x4(a, b, o, 8);
}

void probe_sub_x4_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_sub_x4(a, b, o, 3);
}

void probe_scale_x4_3(const fv_vec2 *a, fv_vec2 *o)
{
    fv_vec2_scale_x4(a, 0.5f, o, 3);
}

void probe_add_x4_mask_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o, const uint8_t *m)
{
    fv_vec2_add_x4_mask(a, b, o, m, 3);
}

void probe_add_x4_idx_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o, const uint32_t *i)
{
    fv_vec2_add_x4_idx(a, b, o, i, 3);
}

void probe_scale_x4_idx_3(const fv_vec2 *a, fv_vec2 *o, const uint32_t *i)
{
    fv_vec2_scale_x4_idx(a, 0.5f, o, i, 3);
}

void probe_add_soa_3(const float *ax, const float *ay, const float *bx, const float *by,
                     float *ox, float *oy)
{
    fv_vec2_add_soa(ax, ay, bx, by, ox, oy, 3);
}

void probe_scale_soa_3(const float *ax, const float *ay, float *ox, float *oy)
{
    fv_vec2_scale_soa(ax, ay, 0.5f, ox, oy, 3);
}

/* ---- P2 新增：_stream 与 _mask_inplace（第一段的五个 + 第二段的一个） ---- */
void probe_add_x4_stream_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_add_x4_stream(a, b, o, 3);
}

void probe_sub_x4_stream_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_sub_x4_stream(a, b, o, 3);
}

void probe_scale_x4_stream_3(const fv_vec2 *a, fv_vec2 *o)
{
    fv_vec2_scale_x4_stream(a, 0.5f, o, 3);
}

void probe_add_soa_stream_3(const float *ax, const float *ay, const float *bx, const float *by,
                            float *ox, float *oy)
{
    fv_vec2_add_soa_stream(ax, ay, bx, by, ox, oy, 3);
}

void probe_scale_soa_stream_3(const float *ax, const float *ay, float *ox, float *oy)
{
    fv_vec2_scale_soa_stream(ax, ay, 0.5f, ox, oy, 3);
}

void probe_add_x4_mask_inplace_3(fv_vec2 *a, const fv_vec2 *b, const uint8_t *m)
{
    fv_vec2_add_x4_mask_inplace(a, b, m, 3);
}

/* ---- 第二段：有 mul+add 收缩风险，期望每批一次 TU 内调用（入口不进消费者） ---- */
void probe_lerp_x4_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_lerp_x4(a, b, 0.5f, o, 3);
}

void probe_madd_x4_3(const fv_vec2 *a, const fv_vec2 *b, fv_vec2 *o)
{
    fv_vec2_madd_x4(a, b, 0.5f, o, 3);
}

float probe_dot_soa_3(const float *ax, const float *ay, const float *bx, const float *by)
{
    return fv_vec2_dot_soa(ax, ay, bx, by, 3);
}

void probe_len_soa_3(const float *ax, const float *ay, float *o)
{
    fv_vec2_len_soa(ax, ay, o, 3);
}

void probe_normalize_soa_3(const float *ax, const float *ay, float *ox, float *oy)
{
    fv_vec2_normalize_soa(ax, ay, ox, oy, 3);
}

void probe_normalize_soa_stream_3(const float *ax, const float *ay, float *ox, float *oy)
{
    fv_vec2_normalize_soa_stream(ax, ay, ox, oy, 3);
}

/* ---- P3：三维 SoA 与四维 _x2 / _soa ----
 * n 的取法与上面两类完全同一个判据：
 *   第一段的 _x2 取 n = 1（bulk = n & ~1 = 0）—— 这样才能证明**整段**（含尾巴循环）被完整内联，
 *     探针里一次 fv_i_ 调用都不剩。vec2 的 _x4 取 n = 3 是同一个道理（3 & ~3 = 0）。
 *   第一段的 _soa 取 n = 3（3 & ~7 = 0），同上。
 *   第二段的入口取 n = 3：入口本身不进消费者（pragma 造成的优化属性差异），探针里恰好一次调用。
 */
void probe_v4_add_x2_1(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o)
{
    fv_vec4_add_x2(a, b, o, 1);
}

void probe_v4_sub_x2_1(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o)
{
    fv_vec4_sub_x2(a, b, o, 1);
}

void probe_v4_scale_x2_1(const fv_vec4 *a, fv_vec4 *o)
{
    fv_vec4_scale_x2(a, 0.5f, o, 1);
}

void probe_v3_add_soa_3(const float *ax, const float *ay, const float *az,
                        const float *bx, const float *by, const float *bz,
                        float *ox, float *oy, float *oz)
{
    fv_vec3_add_soa(ax, ay, az, bx, by, bz, ox, oy, oz, 3);
}

void probe_v3_scale_soa_3(const float *ax, const float *ay, const float *az,
                          float *ox, float *oy, float *oz)
{
    fv_vec3_scale_soa(ax, ay, az, 0.5f, ox, oy, oz, 3);
}

void probe_v4_lerp_x2_3(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o)
{
    fv_vec4_lerp_x2(a, b, 0.5f, o, 3);
}

void probe_v4_madd_x2_3(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o)
{
    fv_vec4_madd_x2(a, b, 0.5f, o, 3);
}

float probe_v4_dot_soa_3(const float *ax, const float *ay, const float *az, const float *aw,
                         const float *bx, const float *by, const float *bz, const float *bw)
{
    return fv_vec4_dot_soa(ax, ay, az, aw, bx, by, bz, bw, 3);
}

void probe_v4_len_soa_3(const float *ax, const float *ay, const float *az, const float *aw,
                        float *o)
{
    fv_vec4_len_soa(ax, ay, az, aw, o, 3);
}

void probe_v4_normalize_soa_3(const float *ax, const float *ay, const float *az, const float *aw,
                              float *ox, float *oy, float *oz, float *ow)
{
    fv_vec4_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, 3);
}

float probe_v3_dot_soa_3(const float *ax, const float *ay, const float *az,
                         const float *bx, const float *by, const float *bz)
{
    return fv_vec3_dot_soa(ax, ay, az, bx, by, bz, 3);
}

void probe_v3_len_soa_3(const float *ax, const float *ay, const float *az, float *o)
{
    fv_vec3_len_soa(ax, ay, az, o, 3);
}

void probe_v3_normalize_soa_3(const float *ax, const float *ay, const float *az,
                              float *ox, float *oy, float *oz)
{
    fv_vec3_normalize_soa(ax, ay, az, ox, oy, oz, 3);
}

#ifdef __cplusplus
} /* extern "C" */
#endif