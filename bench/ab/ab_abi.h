/* bench/ab/ab_abi.h —— 变体 TU 与驱动之间的唯一 ABI
 *
 * 为什么要有这层：include/flashvec.h 正在被重构成"头文件内联"形态（默认
 * static inline；-DFV_USE_EXTERN 才回到库内定义）。两种形态下批量入口的**链接属性
 * 不同**（内联形态是 TU 私有的 static，库形态是外部符号），因此驱动不可能用
 * "声明 + 函数指针"去调它们。
 *
 * 解决办法：把计时循环放进每个变体自己的 TU（见 ab_variant.cpp），由它把
 * "一次取样"整段导出。于是：
 *   - 库形态  -> TU 内是跨函数的普通 call（正是链接静态库的调用方拿到的形态）；
 *   - 内联形态 -> TU 内直接被内联（正是重构想要的形态）。
 * 两种形态都能测，且都不需要驱动知道符号名。
 *
 * 本头文件**不包含任何库类型**：指针一律 void * / float * / uint32_t *，A、B 两侧即使
 * 头文件版本不同也不会因为结构体布局不同而读出垃圾。布局一致性另有一道运行期检查
 * （AB_ABI_NAME 导出的 sizeof 指纹）。
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
    float *ox2, *oy2;         /* collision 的第二组输出（反射后的速度） */
    float *px, *py, *qx, *qy; /* 管线 scratch（collision 的输入法线/速度；chain3 的 t/u） */
    /* 宽元素（vec3 12 B / vec4 16 B）的分量数组；只在这些 op 被选中时非空。
     * oz/ow 与 vx/vy 同址（base2 的前后两半），两者不会同时使用。 */
    const float *az, *bz, *aw, *bw;
    float *oz, *ow;
};

/* op 编号：由 ab_ops.h 的注册表派生（唯一真源）。两侧变体 TU 与驱动共用同一份。 */
enum {
#define AB_ABI_ENUM_ITEM(ID, ...) AB_OP_##ID,
    AB_LIB_OPS(AB_ABI_ENUM_ITEM)
    AB_LIB_OPS_EXTRA(AB_ABI_ENUM_ITEM)
#undef AB_ABI_ENUM_ITEM
    AB_OP_COUNT
};

/* 变体侧向驱动报告"这一个 op 的 store 真的怎么写的"。
 * 驱动用它在 RFO 记账里决定输出数组的真实流量（普通 store = 2x、NT store = 1x）。
 * 这是**运行期**信息，因为同一份 ab_variant.cpp 会被编两次（两侧头文件能力不同）：
 * 例：baseline 头文件没有 _stream 时，thunk 退化成普通 store，这里就必须报 0。 */
#define AB_SIDE_NT 0x1u

#ifdef __cplusplus
extern "C" {
#endif

/* 每个变体导出三个符号；名字由 -DAB_LOOP_NAME / -DAB_ABI_NAME / -DAB_FLAGS_NAME 指定 */
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
