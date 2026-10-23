/*
 * flashvec -- 轻量向量基础库
 *
 * 设计约定：
 *   - 头文件是纯 C（C11）接口，C 与 C++ 都能直接包含；
 *   - 实现文件用 C++ 编译，但只写 C 语言子集，将来要用 C++ 特性时
 *     不必更换工具链，也不必改 ABI；
 *   - 所有量都是按值传递的小结构体，库内不分配内存、无全局状态。
 *   - 头文件里不出现任何 intrinsics 或目标相关的编译开关：批量入口在这里只是普通 C 声明，
 *     AVX2 内核在 src/flashvec_batch.cpp 里由运行时派发启用（无 AVX2 的机器走标量路径）。
 */
#ifndef FLASHVEC_H
#define FLASHVEC_H

#include <stddef.h> /* size_t：批量入口用 */
#include <stdint.h> /* uint8_t / uint32_t：掩码与索引变体用 */

#ifdef __cplusplus
extern "C" {
#endif

/* 标量精度：整体切到 double 只需要改这一行（实现里的 sqrtf 也要一起改）。 */
typedef float fv_real;

/* 两个零判定阈值，量纲不同，不要混用：
 *   FV_EPS     —— 线性量（长度）阈值，用于归一化与 near 比较；
 *   FV_EPS_SQR —— 平方量阈值 = FV_EPS*FV_EPS，手里已经是长度平方或点积时用它。
 * 用 FV_EPS 去比长度平方会把阈值放宽 1000 倍（project 曾经的量纲 bug，P0 已修）。 */
#define FV_EPS ((fv_real)1e-6)
#define FV_EPS_SQR ((fv_real)1e-12)

typedef struct fv_vec2 {
    fv_real x;
    fv_real y;
} fv_vec2;

typedef struct fv_vec3 {
    fv_real x;
    fv_real y;
    fv_real z;
} fv_vec3;

/* 轴对齐包围盒：lo 是逐分量下界，hi 是逐分量上界（不要求 lo <= hi，见 clamp 的说明）。 */
typedef struct fv_aabb2 {
    fv_vec2 lo;
    fv_vec2 hi;
} fv_aabb2;

/* ---------------- 二维 ---------------- */

fv_vec2 fv_vec2_make(fv_real x, fv_real y);
fv_vec2 fv_vec2_add(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_sub(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_neg(fv_vec2 v);
fv_vec2 fv_vec2_scale(fv_vec2 v, fv_real s);
fv_vec2 fv_vec2_lerp(fv_vec2 a, fv_vec2 b, fv_real t);

/* 逐分量最小 / 最大值。定义是 min(a,b) = (a < b) ? a : b，max(a,b) = (a > b) ? a : b，
 * 不是 IEEE minNum/maxNum：一侧为 NaN 时结果由比较决定
 * （min(NaN, x) == x，min(x, NaN) == NaN；两个都是 NaN 时返回后者）。 */
fv_vec2 fv_vec2_min(fv_vec2 a, fv_vec2 b);
fv_vec2 fv_vec2_max(fv_vec2 a, fv_vec2 b);

/* 逐分量夹取：x < lo ? lo : (x > hi ? hi : x)。
 * lo > hi 时按同一表达式取值（x < lo 给 lo，否则 x > hi 给 hi），不另作承诺；
 * NaN 输入按上面同一条比较规则落位（不与边界比较为真，原样返回）。 */
fv_vec2 fv_vec2_clamp(fv_vec2 v, fv_vec2 lo, fv_vec2 hi);

/* a + b*s（逐分量）。实现不做 FMA 收缩，逐位等于标量表达式 a.x + b.x*s。 */
fv_vec2 fv_vec2_madd(fv_vec2 a, fv_vec2 b, fv_real s);

fv_real fv_vec2_dot(fv_vec2 a, fv_vec2 b);
fv_real fv_vec2_len_sqr(fv_vec2 v);
fv_real fv_vec2_len(fv_vec2 v);
fv_real fv_vec2_dist_sqr(fv_vec2 a, fv_vec2 b);
fv_real fv_vec2_dist(fv_vec2 a, fv_vec2 b);

/* 零向量（长度 <= FV_EPS）返回 (0, 0)，不会放大出 NaN。 */
fv_vec2 fv_vec2_normalize(fv_vec2 v);

/* v 在 onto 方向上的投影分量。
 * onto 的长度平方 <= FV_EPS_SQR（等价长度 <= FV_EPS）时返回 (0, 0)。
 * P0 修正：这里曾经用线性阈值 FV_EPS 去比长度平方，等价长度阈值是 1e-3，
 * 比 normalize 松 1000 倍；现在改用平方阈值 FV_EPS_SQR，与 normalize 对齐。 */
fv_vec2 fv_vec2_project(fv_vec2 v, fv_vec2 onto);

/* 把点 p 并入包围盒：lo 取逐分量较小者，hi 取逐分量较大者。 */
fv_aabb2 fv_aabb2_merge(fv_aabb2 box, fv_vec2 p);

/* 相交判定：相切（a.hi == b.lo 之类）算相交，返回 1；否则返回 0。 */
int fv_aabb2_overlaps(fv_aabb2 a, fv_aabb2 b);

/* 逐分量近似相等，返回 1 / 0。 */
int fv_vec2_near(fv_vec2 a, fv_vec2 b, fv_real eps);

/* ---------------- 三维 ---------------- */

fv_vec3 fv_vec3_make(fv_real x, fv_real y, fv_real z);
fv_vec3 fv_vec3_add(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_sub(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_neg(fv_vec3 v);
fv_vec3 fv_vec3_scale(fv_vec3 v, fv_real s);
fv_vec3 fv_vec3_lerp(fv_vec3 a, fv_vec3 b, fv_real t);

/* 与 fv_vec2_min / max / clamp / madd 同样的语义与 NaN 规则，逐分量进行。 */
fv_vec3 fv_vec3_min(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_max(fv_vec3 a, fv_vec3 b);
fv_vec3 fv_vec3_clamp(fv_vec3 v, fv_vec3 lo, fv_vec3 hi);
fv_vec3 fv_vec3_madd(fv_vec3 a, fv_vec3 b, fv_real s);

fv_real fv_vec3_dot(fv_vec3 a, fv_vec3 b);
fv_real fv_vec3_len_sqr(fv_vec3 v);
fv_real fv_vec3_len(fv_vec3 v);
fv_real fv_vec3_dist_sqr(fv_vec3 a, fv_vec3 b);
fv_real fv_vec3_dist(fv_vec3 a, fv_vec3 b);

fv_vec3 fv_vec3_normalize(fv_vec3 v);
/* 与 fv_vec2_project 同一阈值规则（用 FV_EPS_SQR 比长度平方）。 */
fv_vec3 fv_vec3_project(fv_vec3 v, fv_vec3 onto);
/* 叉积 a x b。 */
fv_vec3 fv_vec3_cross(fv_vec3 a, fv_vec3 b);

int fv_vec3_near(fv_vec3 a, fv_vec3 b, fv_real eps);

/* ============================ 批量入口（P1） ============================
 *
 * 为什么要有这一组函数：逐元素 API 是跨编译单元的函数调用，GCC 直接拒绝向量化
 * （missed: statement clobbers memory），实测边界成本约 4.2 cycle/元素，是打包收益的 3 倍。
 * 批量入口把「每个元素一次 call」换成「一次 call 处理一批」，库内部再决定用不用 SIMD。
 *
 * 布局
 *   _x4   ：输入输出都是 AoS 的 fv_vec2 数组。fv_vec2 是 8 B，4 个正好 32 B，
 *           一个 256-bit 寄存器装 4 个元素，**零重排**，调用方现有数据直接传。
 *   _soa  ：输入输出是分离的分量数组（x[]、y[]）。归约型运算（dot/len/normalize）用它最自然。
 *   _mask ：整批读改写，掩码外的 lane 输出 a[i] 的原值（详见下面的掩码语义）。
 *   _idx  ：只处理 idx[0..m) 指定的下标，真的跳过其余元素（代价是随机访存）。
 *
 * 位级一致性的承诺边界（这是调用方最需要知道的一条）
 *   保证 - 与对应标量 API 逐位一致（同输入下输出的 float 位模式 memcmp 相等）：
 *     fv_vec2_add_x4 / sub_x4 / scale_x4 / lerp_x4 / madd_x4；
 *     fv_vec2_add_x4_mask（含掩码外写回 a[i] 的那一半）；
 *     fv_vec2_add_x4_idx / scale_x4_idx；
 *     fv_vec2_add_soa / scale_soa / len_soa / normalize_soa。
 *   不保证 - fv_vec2_dot_soa：它用 4 个独立累加器 + 末端一次归约打断依赖链
 *     （实测 x2.25），归约顺序与逐元素路径不同，因此**结果依赖归约顺序**：
 *     同输入同输出（确定性），但不是逐位一致。
 *     误差没有闭式上界，下面全是实测值（语料不同可差两个数量级，引用请带条件）：
 *       · n <= 1000 随机数据、按 scale 归一化比较：相对误差 <= 4e-8（tests/test_batch.cpp）
 *       · n = 8192 单次点积、对比逐元素路径：绝对误差 <= 4.8e-7，相对误差 <= 6.2e-8
 *       · 结果接近抵消（相对误差的分母很小）时可达 ~3e-6 相对误差
 *     需要逐位复现时请用标量 fv_vec2_dot。
 *     n < 8 时它退化为纯顺序累加，这时与逐元素路径逐位一致。
 *   保证是怎么来的：这些路径只用 vaddps / vsubps / vmulps / vdivps / vsqrtps，
 *   不用 FMA（vfmadd）、不用倒数近似（rcpps/rsqrtps），归约顺序（除 dot_soa）与标量一致。
 *   tests/test_batch.cpp 对每个函数、每种尾长做 memcmp 级别的对拍。
 *
 * 尾长
 *   n 可以是任意值，包括 0。不足一个向量宽度的尾巴（_x4 是 0-3 个，_soa 是 0-7 个）
 *   走标量路径，语义与主循环完全相同。n == 0 时不读不写任何指针（可以传 NULL）。
 *
 * 掩码语义（容易误读，必须写清）
 *   mask[i] != 0 时 o[i] = 运算结果；mask[i] == 0 时 o[i] 被**显式写回 a[i] 的原值**。
 *   不是「不写」——256-bit 寄存器只能整条写回，跳过存储做不到，而且跳过会在 o 里
 *   留下与 a 不一致的半行数据。想要「保留 o 的旧值」的语义，请用 _idx 变体，
 *   或在调用前让 a 指向那份旧值。
 *
 * 索引语义
 *   对 k ∈ [0, m)：i = idx[k]，然后 o[i] = 运算结果（读 a[i]、b[i]）。
 *   idx 可以重复、可以乱序、不必单调；重复下标写入的值相同，因此与处理顺序无关。
 *   idx 之外的下标既不读也不写，o 的其他位置保持调用前的值。
 *
 * 指针与别名
 *   a / b / o 带 __restrict，按语义不得互相重叠（o 与 a 或 b 同址是 UB）；
 *   idx / mask 也不得与它们重叠。
 *
 * 其他
 *   - n 的量纲是「元素个数」，不是字节数。
 *   - 零向量保护在 normalize_soa 里用掩码而不是分支，语义与 fv_vec2_normalize 相同
 *     （len <= FV_EPS 输出 (0,0)）；实现上先把长度替换成 1 再做除法，所以零向量 lane
 *     不会触发除零（MXCSR 标志也不会被置起）。
 * ====================================================================== */

void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n);
void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                    fv_vec2 *__restrict o, size_t n);
void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                      fv_vec2 *__restrict o, size_t n);
void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real t, fv_vec2 *__restrict o, size_t n);
void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                     fv_real s, fv_vec2 *__restrict o, size_t n);

/* 稀疏变体：mask[i] != 0 才做加法，否则 o[i] = a[i]（写回 a 的原值）。 */
void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                         fv_vec2 *__restrict o, const uint8_t *__restrict mask, size_t n);
/* 稀疏变体：对 k ∈ [0,m) 取 i = idx[k]，做 o[i] = a[i] + b[i]；idx 可重复、可乱序。 */
void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                        fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m);
void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                          fv_vec2 *__restrict o, const uint32_t *__restrict idx, size_t m);

/* SoA 分量数组入口：6 个指针 + n，各自 __restrict。 */
void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict bx, const float *__restrict by,
                     float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay, fv_real s,
                       float *__restrict ox, float *__restrict oy, size_t n);
/* 唯一不保证逐位一致的入口：4 累加器 + 末端一次归约，结果依赖归约顺序（见上）。 */
float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict bx, const float *__restrict by, size_t n);
void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                     float *__restrict o, size_t n);
/* 语义与逐元素 fv_vec2_normalize 一致：len <= FV_EPS 输出 (0, 0)。 */
void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           float *__restrict ox, float *__restrict oy, size_t n);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FLASHVEC_H */
