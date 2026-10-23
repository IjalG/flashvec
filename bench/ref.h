/* bench/ref.h -- 严格 IEEE 语义参考实现（与三档被测实现对照用）
 *
 * 单独一个编译单元，固定用 -O2 -ffp-contract=off -fno-tree-vectorize -fno-fast-math 编译，
 * 因此参考结果不受被测实现的编译选项影响，可以当作精度基准。
 */
#ifndef BENCH_REF_H
#define BENCH_REF_H

#ifdef __cplusplus
extern "C" {
#endif

void ref_add(const float *ax, const float *ay, const float *bx, const float *by,
             float *ox, float *oy, unsigned n);
void ref_scale(const float *ax, const float *ay, float s, float *ox, float *oy, unsigned n);
float ref_dot_sum(const float *ax, const float *ay, const float *bx, const float *by, unsigned n);
/* 逐元素点积 a_i·b_i，结果同时写进 ox[i] 与 oy[i]（三档共用同一语义）。 */
void ref_dot_elem(const float *ax, const float *ay, const float *bx, const float *by,
                  float *ox, float *oy, unsigned n);
void ref_len(const float *ax, const float *ay, float *ox, float *oy, unsigned n);
void ref_normalize(const float *ax, const float *ay, float *ox, float *oy, unsigned n);
void ref_collision(const float *px, const float *py, const float *qx, const float *qy,
                   const float *vx, const float *vy, float *ox, float *oy, unsigned n);

#ifdef __cplusplus
}
#endif

#endif /* BENCH_REF_H */
