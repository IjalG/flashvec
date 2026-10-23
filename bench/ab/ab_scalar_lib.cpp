/* bench/ab/ab_scalar_lib.cpp —— "重构前"锚点的**库侧**：把 HEAD 快照的 src/flashvec.cpp
 * 原样编成一个独立 TU，逐元素标量 API 在这里有真实的外部符号。
 *
 * 为什么要拆成两个 TU（本文件 + ab_scalar_loop.cpp）：只有在**另一个** TU 里调用，
 * GCC 才拿不到函数体、不能内联、不能向量化——这才是"消费者链接预编译静态库"的真实
 * 形态（README 的档 1）。若写在同一 TU，它会内联展开成档 2（scalar_flat），
 * 那个数字不能当锚点用。
 */
#include "snapshot/include/flashvec.h"
#include "snapshot/src/flashvec.cpp"
