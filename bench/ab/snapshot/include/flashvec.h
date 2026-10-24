/*
 * flashvec -- 轻量向量基础库
 *
 * 设计约定：
 *   - 头文件是纯 C（C11）接口，C 与 C++ 都能直接包含；
 *   - 实现文件用 C++ 编译，但只写 C 语言子集，将来要用 C++ 特性时
 *     不必更换工具链，也不必改 ABI；
 *   - 所有量都是按值传递的小结构体，库内不分配内存、无全局状态。
 *
 * 批量入口的实现形态（P1 重构：从「库内实现」改为「头文件 inline 实现」）
 *   - 默认形态：33 个批量入口在本头文件的「批量入口的实现」一节里以 static inline 定义，
 *     调用方的编译单元直接拿到实现，不跨编译单元调用；
 *   - 库内同名 extern 符号继续导出（src/flashvec_batch.cpp），旧消费方照旧能链；
 *   - intrinsics 的出现条件与语言无关，只与「编译器 + 架构」有关：x86/x86_64 且编译器是
 *     GCC 或 Clang —— **C 编译单元与 C++ 编译单元走同一条规则**。
 *     其余情况（MSVC 等不支持 target 属性的编译器、非 x86 架构）走「可自动向量化的纯 C」
 *     路径：那种编译单元里不会展开出 immintrin.h / __m256 / _mm_* 任何一行。
 *   - 有 intrinsics 不等于会执行：AVX2 内核是否运行由运行时
 *     __builtin_cpu_supports("avx2") != 0 决定；没有 AVX2 的机器与所有尾长走标量路径，
 *     标量表达式与 src/flashvec.cpp 的逐元素 API 逐字相同。
 *   - 上面四条都可以直接验证（在仓库根目录）：
 *       cc  -std=c11   -Iinclude -E tests/test_header_c.c    | grep -c immintrin   # > 0
 *       c++ -std=c++17 -Iinclude -E tests/test_header_cpp.cpp | grep -c immintrin  # > 0
 *       cc  -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -c tests/test_header_c.c  # 0 告警
 *       nm -g --defined-only build/libflashvec.a | grep -c '_x4\|_x2\|_soa'        # 33
 *       make test-inline           # 默认形态 vs FV_USE_EXTERN vs 纯 C 回退，逐字节对拍
 *       make test-inline-effects   # 消费者默认 -O2 下：第一段完整内联、第二段每批一次调用
 *       make test-fma              # 消费者带 -march=native -mfma -ffp-contract=fast 仍位级一致
 *
 * 浮点契约与编译选项的责任边界
 *   - 实现小节里**有收缩风险的那一半**由 #pragma GCC optimize("fp-contract=off")（GCC，
 *     push/pop 成对）罩住，禁止把 a + b*s 收缩成 FMA。它的作用范围是「本头文件内、那几个
 *     函数的代码生成」，对调用方自己的代码、以及调用方的编译选项没有任何作用；
 *   - target("avx2") 挡不住调用方命令行带来的 FMA：实测 -march=native -mfma
 *     -ffp-contract=fast 下，带 target("avx2") 的内核里照样出现 vfmadd。
 *     所以在消费者场景下，**pragma 是唯一的防线**，不是 target 属性；
 *   - 因此：调用方若使用 -ffast-math（含 -fassociative-math / -ffinite-math-only 之类的
 *     重结合与无 NaN 假设），或在不支持该 pragma 的编译器上让 FMA 收缩生效
 *     （-march=native / -mfma 配 -ffp-contract=fast），**逐位一致的承诺由调用方自己负责**；
 *   - 想让结果由库自己的编译选项（-O3 -fno-math-errno -ffp-contract=off）保证，
 *     就定义 FV_USE_EXTERN 调用库内的实现。
 */
#ifndef FLASHVEC_H
#define FLASHVEC_H

#include <stddef.h> /* size_t：批量入口用 */
#include <stdint.h> /* uint8_t / uint32_t：掩码与索引变体用 */

/* 只有「本头文件提供实现」时才需要这两个头：FV_USE_EXTERN 形态的包含集与重构前完全一致。 */
#if !defined(FV_USE_EXTERN)
#include <math.h>   /* sqrtf：批量入口的标量尾路径（含纯 C 路径） */
#include <string.h> /* memcpy：_mask 变体装载 4 个掩码字节 */
#endif

/* ---------------- 编译器/平台能力探针 ----------------
 * FV_I_X86（0/1）      ：是否 x86 / x86_64；
 * FV_I_INTRINSICS（0/1）：批量入口能否使用 AVX2 intrinsics —— x86 + GCC/Clang，
 *                        与语言无关（C 与 C++ 编译单元同规则）；为 0 时整组入口走纯 C 路径。
 * 这两个宏是本头文件的实现细节，但为了可验证而保留下来（tests/ 里直接断言它们的取值）。 */
#if defined(__x86_64__) || defined(__i386__)
#define FV_I_X86 1
#else
#define FV_I_X86 0
#endif

/* FV_I_NO_INTRINSICS：自测/降级用开关，定义它就强制走纯 C 路径。
 * 存在的理由：在 x86 + GCC 的机器上，纯 C 路径否则没有任何测试覆盖，而它正是 MSVC /
 * 非 x86 用户真正会走的那条路（tests/test_parity.c 的 c_pure 变体用它做位级对拍）。 */
#if FV_I_X86 && (defined(__GNUC__) || defined(__clang__)) && !defined(FV_I_NO_INTRINSICS)
#define FV_I_INTRINSICS 1
#else
#define FV_I_INTRINSICS 0
#endif

#if FV_I_INTRINSICS && !defined(FV_USE_EXTERN)
#include <immintrin.h>
#endif

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

/* 四维：4 × fv_real = 16 B，**没有填充字节**，分量按 x, y, z, w 连续存放。
 * 16 B 这个数字是批量入口 fv_vec4_*_x2 的前提（2 个元素 = 32 B = 一个 256-bit 寄存器），
 * 所以用编译期断言把它钉住：要改这个布局，必须先一起改 _x2 内核。 */
typedef struct fv_vec4 {
    fv_real x;
    fv_real y;
    fv_real z;
    fv_real w;
} fv_vec4;

#if defined(__cplusplus)
static_assert(sizeof(fv_vec4) == 16, "fv_vec4 must be exactly 16 B (no padding): _x2 batches rely on it");
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
_Static_assert(sizeof(fv_vec4) == 16, "fv_vec4 must be exactly 16 B (no padding): _x2 batches rely on it");
#endif

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

/* ---------------- 四维 ----------------
 *
 * 布局：4 × fv_real = 16 B，没有填充字节（见上面的 typedef 与编译期断言）。
 *
 * 语义与 fv_vec2 / fv_vec3 逐条对齐 —— 同一套名字、同一组阈值、同一条 NaN 规则：
 *   min(a,b) = (a < b) ? a : b，max(a,b) = (a > b) ? a : b，**不是** IEEE minNum/maxNum：
 *   一侧为 NaN 时结果由比较决定（min(NaN, x) == x，min(x, NaN) == NaN，两个都是 NaN 时
 *   返回后者）；clamp 里与边界比较都不为真的分量原样返回。这条规则与 fv_vec2_min/max/clamp
 *   完全一致，可以直接互相对拍。
 *   madd(a,b,s) = a + b*s，不做 FMA 收缩（逐位等于标量表达式）。
 *   dot 的求和顺序是左结合的 ((x*x' + y*y') + z*z') + w*w'；
 *   len = sqrtf(dot(v, v))；normalize 在 len <= FV_EPS 时返回 (0, 0, 0, 0)；
 *   project 用平方量阈值 FV_EPS_SQR 去比 onto 的长度平方（与 fv_vec2_project 同一条规则）；
 *   near 是逐分量 |a-b| <= eps。
 *   不提供 cross：四维叉积不是标准运算，没有明确需求（YAGNI）。
 *
 * 实现位置：src/flashvec_vec4.cpp —— 与 fv_vec2/fv_vec3 一样是库内 extern 符号，
 * 不是头文件 inline。两个理由：
 *   1) len / normalize 需要 sqrtf；放头文件就得让 FV_USE_EXTERN 形态也包含 <math.h>，
 *      那会破坏「FV_USE_EXTERN 的包含集与重构前一致」这条既有约定；
 *   2) 标量 API 与批量入口的标量尾巴共用同一个 -ffp-contract=off 编译单元，
 *      位级一致更好讲、也更好查。 */
fv_vec4 fv_vec4_make(fv_real x, fv_real y, fv_real z, fv_real w);
fv_vec4 fv_vec4_add(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_sub(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_neg(fv_vec4 v);
fv_vec4 fv_vec4_scale(fv_vec4 v, fv_real s);
fv_vec4 fv_vec4_lerp(fv_vec4 a, fv_vec4 b, fv_real t);
fv_vec4 fv_vec4_min(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_max(fv_vec4 a, fv_vec4 b);
fv_vec4 fv_vec4_clamp(fv_vec4 v, fv_vec4 lo, fv_vec4 hi);
fv_vec4 fv_vec4_madd(fv_vec4 a, fv_vec4 b, fv_real s);

fv_real fv_vec4_dot(fv_vec4 a, fv_vec4 b);
fv_real fv_vec4_len_sqr(fv_vec4 v);
fv_real fv_vec4_len(fv_vec4 v);
fv_real fv_vec4_dist_sqr(fv_vec4 a, fv_vec4 b);
fv_real fv_vec4_dist(fv_vec4 a, fv_vec4 b);

fv_vec4 fv_vec4_normalize(fv_vec4 v);
fv_vec4 fv_vec4_project(fv_vec4 v, fv_vec4 onto);

int fv_vec4_near(fv_vec4 a, fv_vec4 b, fv_real eps);

/* ============================ 批量入口（P1） ============================
 *
 * 为什么要有这一组函数：逐元素 API 是跨编译单元的函数调用，GCC 直接拒绝向量化
 * （missed: statement clobbers memory），**早期探针**实测边界成本约 4.2 cycle/元素，是打包收益的
 * 3 倍（早期微基准；本轮加硬夹具只复测批量入口本身，未复测这条边界成本，只能按「量级」读）。
 * 批量入口把「每个元素一次 call」换成「一次 call 处理一批」，库内部再决定用不用 SIMD。
 *
 * 实现形态（三个宏就是全部开关）
 *   默认            本头文件给出 33 个批量入口的 static inline 实现（见文末「批量入口的实现」），
 *                   调用方的编译单元直接拿到函数体；
 *   FV_USE_EXTERN   头文件只给下面这组 extern 声明，调用进 libflashvec.a —— 重构前的行为；
 *                   也是 C 调用方想用 AVX2 内核时的入口（C 编译单元不走 intrinsics，见文件头）；
 *   FV_BUILD        库自身构建用（src/flashvec_batch.cpp 定义它），对外表现同 FV_USE_EXTERN：
 *                   只有声明；定义在 .cpp 里，函数体一行转发到 fv_i_*。
 *   两种形态是同一份实现：.cpp 不复制任何内核，只把 fv_i_* 接到公开名上。
 *
 * 布局
 *   _x4   ：输入输出都是 AoS 的 fv_vec2 数组。fv_vec2 是 8 B，4 个正好 32 B，
 *           一个 256-bit 寄存器装 4 个元素，**零重排**，调用方现有数据直接传。
 *   _soa  ：输入输出是分离的分量数组（x[]、y[]）。归约型运算（dot/len/normalize）用它最自然。
 *   _mask ：整批读改写，掩码外的 lane 输出 a[i] 的原值（详见下面的掩码语义）。
 *   _idx  ：只处理 idx[0..m) 指定的下标，真的跳过其余元素（代价是随机访存）。
 *   _x2   ：fv_vec4 专用（16 B/元素）。2 个元素 = 32 B = 一个 ymm，**正好装满，一个 lane 都
 *           不浪费**。命名与 vec2 的 _x4 同构：数字说的是「一个 256-bit 寄存器装几个元素」
 *           （vec2 是 8 B/元素 -> 4 个；vec4 是 16 B/元素 -> 2 个），不是「一次处理几个向量」。
 *           拼装问题在这里不存在：fv_vec4 无填充，数组里相邻两个元素就是 32 B 连续内存，
 *           所以两条输入各是一条 256-bit loadu、输出是一条 256-bit storeu —— 与 _x4 一样
 *           零重排。要用两次 128-bit load + _mm256_insertf128_ps 拼装的是另一种布局
 *           （元素步长不等于 16 B，比如 AoSoA 或带 stride 的数组），不是这里的 AoS。
 *   三维没有 AoS 批量，只有 _soa：
 *           fv_vec3 保持 12 B，**不加 padding**（填到 16 B 会破坏调用方已有的 12 B 布局，
 *           文档里有明确结论）。12 B 不是 32 B 的整数分之一：AoS 批量要么让每个 ymm 空转
 *           1/4 lane（3/4 利用率），要么先 shuffle 拆包再打包，省下的宽度又赔回去。
 *           所以三维批量的唯一形态是 SoA（fv_vec3_*_soa，签名 4-9 个指针 + n）——
 *           一条 ymm 里装的是**同一分量的 8 个不同元素**，lane 利用率 100%。
 *           fv_vec3_cross 没有批量入口：3 个 yz 交换需要大量 shuffle，收益不明（YAGNI）。
 *
 * 位级一致性的承诺边界（这是调用方最需要知道的一条）
 *   保证 - 与对应标量 API 逐位一致（同输入下输出的 float 位模式 memcmp 相等）：
 *     fv_vec2_add_x4 / sub_x4 / scale_x4 / lerp_x4 / madd_x4；
 *     fv_vec2_add_x4_mask（含掩码外写回 a[i] 的那一半）；
 *     fv_vec2_add_x4_idx / scale_x4_idx；
 *     fv_vec2_add_soa / scale_soa / len_soa / normalize_soa；
 *     fv_vec3_add_soa / scale_soa / len_soa / normalize_soa（三维只有 SoA 形态）；
 *     fv_vec4_add_x2 / sub_x2 / scale_x2 / lerp_x2 / madd_x2；
 *     fv_vec4_len_soa / normalize_soa。
 *   不保证 - fv_vec2_dot_soa：它用 4 个独立累加器 + 末端一次归约打断依赖链
 *     （早期微基准 x2.25，本轮未复测），归约顺序与逐元素路径不同，因此**结果依赖归约顺序**：
 *     同输入同输出（确定性），但不是逐位一致。
 *     误差没有闭式上界，下面全是实测值（语料不同可差两个数量级，引用请带条件）：
 *       · n <= 1000 随机数据、按 scale 归一化比较：相对误差 <= 4e-8（tests/test_batch.cpp）
 *       · n = 8192 单次点积、对比逐元素路径：绝对误差 <= 4.8e-7，相对误差 <= 6.2e-8
 *       · 结果接近抵消（相对误差的分母很小）时可达 ~3e-6 相对误差
 *     需要逐位复现时请用标量 fv_vec2_dot。
 *     n < 8 时它退化为纯顺序累加，这时与逐元素路径逐位一致（+0.0 与 -0.0 的边界见下）。
 *   同样的例外适用于 fv_vec3_dot_soa / fv_vec4_dot_soa：它们与 fv_vec2_dot_soa 用同一个
 *     归约套路（多个独立累加器 + 末端一次树归约），所以同输入同输出、但不是逐位一致；
 *     容差口径也一样：按 scale（各项绝对值之和，double 累加）归一化后相对误差 < 1e-5
 *     （tests/test_batch.cpp 实测 <= 4e-8；本轮移植到 tests/test_parity.c 的语料上见报告）。
 *     需要逐位复现时请用标量 fv_vec3_dot / fv_vec4_dot。
 *     实现细节：n 不足一个 8 lane 时整段只走尾巴，但尾巴是「sum 从 +0.0 开始逐项累加」，
 *     而标量 API 的第一项不与 +0.0 相加；因此当所有项都是 -0.0 时两者会差一个符号位
 *     （+0.0 vs -0.0）。这条与 vec2 的 dot_soa 相同，属于上面那条例外的范畴。
 *   分段不变性（另一条性质，容易漏）：
 *     逐元素型批量入口（add/sub/scale/lerp/madd/len/normalize 及其 _x4/_soa/_mask/_idx 变体）
 *     具备该性质 —— 同一数组切成任意多段分别调用，输出与一次调用逐位相同
 *     （实测 add_x4 12/12 逐位一致）。
 *     归约型入口（dot_soa，含 vec3/vec4 版本）**不具备**该性质：归约树随切分改变。
 *     实测 n=8192「整段 vs 切 2 段」相对差 2.4e-5、n=1000 约 3.9e-6，12/12 全部不一致。
 *     需要「结果与切分无关」时请用标量 fv_vec2_dot 顺序累加。
 *   保证是怎么来的：这些路径只用 vaddps / vsubps / vmulps / vdivps / vsqrtps，
 *   不用 FMA（vfmadd）、不用倒数近似（rcpps/rsqrtps），归约顺序（除 dot_soa）与标量一致。
 *   实现小节按「有没有 mul+add 收缩风险」分成两段（见下面「批量入口的实现」那一节的说明）：
 *   第一段（add/sub/scale 及 mask/idx/soa 同构入口）没有任何可收缩的表达式，不需要 pragma，
 *   因此在消费者默认编译选项下能完整内联（含 n < 4 / n < 8 的尾巴循环）；
 *   第二段（lerp/madd/len/normalize/dot）整体被 #pragma GCC optimize("fp-contract=off")
 *   罩住（push/pop 成对，只覆盖本头文件里那几个函数），代价是每批一次 TU 内调用。
 *   其它编译器与调用方自己的编译选项不在这个 pragma 管辖范围（见文件头「责任边界」）。
 *   验证手段：objdump -d build/libflashvec.a 中 vfmadd / vrcpps / vrsqrtps / vhaddps 计数为 0；
 *   tests/test_batch.cpp 对每个函数、每种尾长做 memcmp 级别的对拍；
 *   make test-inline 把默认形态与 FV_USE_EXTERN 形态的输出逐字节对拍。
 *
 * 尾长
 *   n 可以是任意值，包括 0。不足一个向量宽度的尾巴（_x4 是 0-3 个，_x2 是 0-1 个，
 *   _soa 是 0-7 个）
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
 *
 * ---- _stream（NT store）变体：只在工作集明显大于 L3 时用 ----
 * 做了什么：主体循环用 _mm256_stream_ps（非临时存储）写出结果，绕过 cache，从而**不产生
 * 「先把输出行读进来」的那一遍 RFO 读**。这是它全部的物理意义 —— 上限完全由字节账决定，
 * 而字节账必须**按入口形状逐个算，不能合并**：一个元素搬几个 float，只有该入口的签名说得了。
 *
 *   入口                        读 B/e  写 B/e  名义  含 RFO  上限（真实比）
 *   add_x4 / sub_x4 / add_soa     16      8      24     32    32/24 = 1.333
 *   scale_x4 / scale_soa           8      8      16     24    24/16 = 1.500
 *   normalize_soa                  8      8      16     24    24/16 = 1.500
 *   dot_soa / len_soa（纯读）    16 / 8    0       -      -    无收益（NT 只影响写）
 *
 * 「读」那一列由签名决定：_x4 读 a / b 两条 fv_vec2 数组（16 B）；add_soa 读 ax / ay / bx / by
 * **四条** float 数组，同样是 16 B（不是 8 B）—— 所以它与 _x4 同形，上限同为 1.333。
 * 「写」那一列同理：add_soa 写 ox / oy 两条（8 B）；scale 类只读 a（8 B）、只写 o（8 B）。
 * 早先那一版把 _soa 写成「2 读 2 写、24/16 = 1.5x」，漏掉了 add_soa 读 4 条数组这一层，是错的；
 * 上限只能逐入口算，别把一个数字套到全部变体上。
 *
 * 这些是字节账给出的**期望上限**，不是物理保证。实测 PASS 点落在 1.29-1.58x，而 4M 档的配对
 * IQR 有 66-98%、1M/8M 档有 20-35%，因此**单点超过上限 20% 以内一律视为噪声**（例如 add_soa
 * 实测 1.582 对上限 1.333，超出 19%，属噪声带内），不得据此声称「实测超过了理论上限」。
 * 报 GB/s 时请同时给这两个口径（名义 / 含 RFO）。
 * 什么时候**不要**用（这条比收益更重要）：结果马上要被读、或工作集落在 L2/L3 里 ——
 * 那里的惩罚是 **2-8 倍**，不是一两成。
 * 实测（Ryzen 7 5700U，L3 = 4 MiB/CCX，taskset 绑核；同进程交错 A/B；每点 4 进程 × 2 构建
 * × 2 顺序，含同轮 A/A 地板）：
 *   24 MB 工作集（n=1M）  add 1.350 / sub 1.293 / scale 1.313
 *   192 MB 工作集（n=8M） add 1.422 / sub 1.519 / add_soa 1.582
 *   L2 内（n=8192）       add 0.217 / scale 0.130 / add_soa 0.229 / normalize 0.506（慢 2-8 倍）
 * 注意：1M/8M 档的配对 IQR 仍有 20-35%，上面的比值只能读作「1.3-1.6 倍量级」，单点只读作
 * 「越过没越过同轮 A/A 地板」。本轮判定（docs/research/raw/ab_campaign_verdict.txt）：n=8192
 * 六个 _stream 点全部 FAIL(变慢)；n=4M 六个点全部 INCONCLUSIVE(地板>15%) —— 该档同轮 A/A 地板
 * 1.34-1.61、配对 IQR 66-98%，测不准，不得当判据；n=131072 是黑名单档（见验收协议）；
 * PASS 只出在 1M 与 8M 档（add 1.350 / sub 1.293 / scale 1.313；1.422 / 1.519 / add_soa 1.582）。
 * 不要用单次运行的数字声称超过某个上限，也不要用 4M 档的单点反推上限。
 * 结论：**必须由调用方显式选择，库永远不自动切换**。判据是工作集字节数（元素数 × 24 B）
 * 明显大于 L3；本机 4 MB 的 L3 下 n <= 131072（3 MB）就该用普通版。
 * 对齐：调用方**不需要**保证任何对齐。内核自己把输出推到 32 B 边界（头部最多 7 个元素用
 * 普通 store 处理），主体走 NT，剩下的尾巴继续用普通 store；fv_vec2 的步长是 8 B，若输出
 * 数组连 8 B 对齐都不满足（32 B 边界不可达），内核整段退回普通 storeu —— 结果逐位不变，
 * 只是没有 NT 的收益。SoA 的两条输出由两条独立的内核处理（各自 peel、各自决定能不能 NT），
 * 不会因为两块输出的相位差不同而互相拖累；normalize 因为两条输出共用同一次计算，要求
 * ox / oy 能**同时**推到 32 B 边界，推不到就整段退回普通 storeu（建议两块输出都按 32 B 分配）。
 * 顺序：每个 stream 内核返回前执行 _mm_sfence()。NT store 是**弱序**写，数据先在写合并缓冲里；
 * 不 sfence 就返回，别的核（以及随后的 DMA / 设备读）可能看不到这些写。
 * 这不是理论担忧：把 sfence 去掉之后，用「写线程 NT store -> release 标志 -> 读线程 acquire
 * 后逐元素核对」的交接实验跑 12 次，**6 次读到旧值**；带 sfence 时 12 次全对。
 * 同时注意：单次读取的检出率约 50%，所以任何「读一遍对不对」的测试都不能声称稳定咬住这个
 * bug —— 正规测法是多次重复 + 双线程交接（tests/test_stream.c 就是这么写的）。
 * 结果：与非 stream 版**逐位相同**（同一批 vaddps / vsubps / vmulps / vdivps / vsqrtps，
 * 只换了存储指令）。tests/test_parity.c 对每个入口、n = 0..17 全尾长 + 随机大 n + 输出指针
 * 逐个相位平移都做 memcmp，并把输出记进跨形态对拍流。
 * 分段：normalize_soa_stream 有 x*x + y*y，按判据属于**第二段**（pragma 之内，消费者里每批
 * 一次调用）；其余 5 个没有 mul+add，属于第一段（消费者默认选项下可完整内联）。
 *
 * ---- _mask_inplace：就地稀疏加（语义入口；实测无吞吐优势） ----
 * 为什么是**新函数**而不是改 _mask：
 *   _mask 的契约是「整批读改写」—— 掩码外的 lane **显式写回 a[i] 的原值**（见「掩码语义」），
 *   而且要求输出指针 o 与 a/b 不重叠（o 与 a 同址是 UB）。就地写的调用方恰恰不能这样：
 *   他没有第二块缓冲区。同一组参数在 _mask 里既不能表示「新输出」又表示「被原地改写的旧
 *   数组」，所以按语义拆成两个入口，而不是放宽 _mask 的契约。
 *
 * 承诺的是**值语义**，只有三条：
 *   1) mask[i] != 0 的元素逐位等于标量 a[i] + b[i]（也与非 inplace 版 _mask 的命中 lane 逐位
 *      相同 —— 两者用的是同一条 vaddps 路径）；
 *   2) mask[i] == 0 的元素**值**保持不变（可验证：与调用前的字节比较，tests/test_stream.c 就是
 *      这么查的）；
 *   3) 其余（未命中元素是否被物理写回、写回几次、cache 行状态如何）**不做任何承诺**。
 * 特别地：**本函数不是同步原语**，不构成任何并发读写之间的顺序保证；它与 fv_vec2_add 一样，
 * 对同一块内存的并发读写仍然是数据竞争（UB）。需要「物理上不写」的场景不是向量库的职责。
 * 为什么不做「只写命中 lane」：那需要 VMASKMOVPS（vmaskmovps maskstore）。早期探针在本机
 * （Zen 2）测到这一族是微码级慢指令，同密度下比整块 storeu 慢 1.2-3.8x，低密度下又比
 * 「块级预扫 + 标量」慢（1% 命中、64 MB 工作集时约 2.5x）—— 这组数字出自加硬夹具之前的
 * 探针，本轮**未复测**，只能按「量级」读。真正与数字无关的那条理由是：若保留它，本函数
 * 「物理上只写命中 lane」就会变成**随掩码密度变化**的性质 —— 没有任何调用方能依赖一个随
 * 输入分布变化的存储保证，却要为它付性能代价。所以这里选择两种写法都不用。
 *
 * 别名：a 是 read-write，所以不存在 _mask 那种 o/a 别名问题；唯一约束是 b 不得与 a 重叠
 *   （mask 也不得与二者重叠）。
 * 实现是**两段密度自适应**的，两段结果逐位相同，选路只影响性能：
 *   采样：只看前 FV_I_MIP_SAMPLE（64）个掩码字节、数非零个数（几十个周期，相对整批可忽略），
 *         得到量级正确的密度估计；估错了只损失性能，不影响语义。
 *   <= NUM/DEN（默认 1/8）命中：块级预扫 + 逐元素标量。整块四个掩码字节全零就整块跳过
 *         （连 a 的 cache 行都不碰；5% 命中时这类块占 95%^4 ≈ 81%）。
 *   >  NUM/DEN 命中：整块 32 B load a / load b + blend + storeu 写回 a。没有数据相关分支，
 *         未命中的 lane 写回 a 的原值（值不变，但物理上写了）。
 *   分界点由三个宏控制，改宏即可调，不必碰函数体：
 *     FV_I_MIP_DENSE_NUM / FV_I_MIP_DENSE_DEN  —— 密度阈值（默认 1/8）
 *     FV_I_MIP_SAMPLE                          —— 采样字节数（默认 64）
 *   这段选路只决定「走哪条循环」，**不带来任何可测的吞吐差异**（见下面的能力边界）。
 *
 * 能力边界（实测，必读）：本入口归档为「试过但无效」，**不是**「已优化的入口」。
 *   实测（加硬夹具，5 种密度 × 4 档规模，含同轮 A/A 地板）：本入口相对 _mask（非就地）
 *   **没有任何可测的吞吐优势**，增益在 -0.10…+0.00 之间，全部落在 A/A 地板以下
 *   （A/A 地板 1.06-1.56；1%/5%/20%/50%/95% × 4 档规模，判定全部 INCONCLUSIVE）。
 *   原因：就地写不会产生 RFO —— a 本来就要先读后写 —— 所以「就地」省不掉流量。
 *   它的价值是**形态**（省掉一块输出缓冲、调用方少一次分配），不是速度。需要吞吐请用
 *   非就地版本并让编译器看到更宽的循环。
 * 附：_mask 的块级预扫判据是「四个掩码字节搬进一个 uint32 后整值为零」，
 *   **不能**写成 (w & 0x01010101u) == 0 —— 那判的是「每个字节的最低位为 0」，掩码字节取
 *   0x02 / 0x80 / 0xFE 这类合法命中值时会误判成没命中（tests/test_stream.c 有回归用例）。
 * ====================================================================== */

/* 公开名字的两种形态：extern 声明（FV_USE_EXTERN / FV_BUILD）或 static inline 定义（默认）。
 * 二者不能同时出现（同一个 TU 里 extern 声明后再 static 定义是硬错误），所以这里用条件编译。 */
#if defined(FV_USE_EXTERN) || defined(FV_BUILD)

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

/* ---------------- 内存流量优化入口（P2） ----------------
 * _stream：结果用 NT store 写出（不经过 cache，不产生 RFO 的读）。**由调用方显式选择**，
 *   库永远不会自动切换 —— 工作集落在 L2/L3 里时它比普通版慢 2-8 倍（本机 L2 规模 0.13-0.51x）。
 *   输出指针的 32 B 对齐由实现内部做 peel 处理，调用方**不需要**保证任何对齐；
 *   函数返回前已经 sfence，后续读取一定能看到结果。
 * _mask_inplace：就地稀疏加，mask[i] == 0 的元素的**值**保持原值（与 _mask 的「写回 a[i]」
 *   是同一个值语义，但不需要第二块缓冲区）。承诺的是值语义，**不承诺**未命中的元素未被物理
 *   写回，也不是同步原语。**实测相对 _mask 没有任何可测的吞吐优势**（增益 -0.10…+0.00，全部
 *   落在 A/A 地板以下）—— 就地写不产生 RFO，省不掉流量；它的价值是形态，不是速度。
 * 详细契约（含实测数字与能力边界）见本头文件「批量入口（P1）」一节的末尾两小节。 */
void fv_vec2_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n);
void fv_vec2_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                           fv_vec2 *__restrict o, size_t n);
void fv_vec2_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                             fv_vec2 *__restrict o, size_t n);
void fv_vec2_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                            const float *__restrict bx, const float *__restrict by,
                            float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_scale_soa_stream(const float *__restrict ax, const float *__restrict ay, fv_real s,
                              float *__restrict ox, float *__restrict oy, size_t n);
void fv_vec2_normalize_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                  float *__restrict ox, float *__restrict oy, size_t n);
/* 就地稀疏加：mask[i] != 0 时 a[i] += b[i]，否则 a[i] 的值不变。a 是 read-write，
 * 唯一约束是 b 不得与 a 重叠（mask 也不得与二者重叠）。 */
void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                 const uint8_t *__restrict mask, size_t n);

/* ---------------- 三维 / 四维批量入口（P3） ----------------
 * 形态选择（详见上面契约的「布局」一节）：
 *   三维只有 SoA（fv_vec3 保持 12 B，不加 padding），签名 4-9 个指针 + n
 *   （add_soa 9 个、scale_soa / dot_soa / normalize_soa 各 6 个、len_soa 4 个）；
 *   四维 AoS 用 _x2（2 个 16 B 元素 = 32 B = 一个 ymm，装满不浪费 lane），SoA 只给
 *   dot / len / normalize 三个归约型入口（add/sub/scale 在 AoS 里已经是零重排，不需要 SoA 版）。
 * 位级承诺与 vec2 的同类入口完全相同：除 dot_soa（归约顺序不同，见上）外，
 * 全部与对应标量 API 逐位一致。指针带 __restrict，按语义不得互相重叠。 */

void fv_vec3_add_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az,
                     const float *__restrict bx, const float *__restrict by,
                     const float *__restrict bz,
                     float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n);
void fv_vec3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                       const float *__restrict az, fv_real s,
                       float *__restrict ox, float *__restrict oy, float *__restrict oz, size_t n);
/* 唯一不保证逐位一致的 vec3 入口：多累加器 + 末端一次归约（与 fv_vec2_dot_soa 同套路）。 */
float fv_vec3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, size_t n);
void fv_vec3_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, float *__restrict o, size_t n);
/* 语义与 fv_vec3_normalize 一致：len <= FV_EPS 输出 (0, 0, 0)。 */
void fv_vec3_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az,
                           float *__restrict ox, float *__restrict oy, float *__restrict oz,
                           size_t n);

/* 四维 AoS：2 个元素一次（元素 16 B、无填充，两条 loadu 各读 32 B = 2 个连续元素）。 */
void fv_vec4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n);
void fv_vec4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                    fv_vec4 *__restrict o, size_t n);
void fv_vec4_scale_x2(const fv_vec4 *__restrict a, fv_real s, fv_vec4 *__restrict o, size_t n);
void fv_vec4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real t,
                     fv_vec4 *__restrict o, size_t n);
void fv_vec4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b, fv_real s,
                     fv_vec4 *__restrict o, size_t n);

/* 四维 SoA：4 个分量数组（dot 要 8 个）+ n。 */
float fv_vec4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                      const float *__restrict az, const float *__restrict aw,
                      const float *__restrict bx, const float *__restrict by,
                      const float *__restrict bz, const float *__restrict bw, size_t n);
void fv_vec4_len_soa(const float *__restrict ax, const float *__restrict ay,
                     const float *__restrict az, const float *__restrict aw,
                     float *__restrict o, size_t n);
/* 语义与 fv_vec4_normalize 一致：len <= FV_EPS 输出 (0, 0, 0, 0)。 */
void fv_vec4_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                           const float *__restrict az, const float *__restrict aw,
                           float *__restrict ox, float *__restrict oy,
                           float *__restrict oz, float *__restrict ow, size_t n);

#endif /* FV_USE_EXTERN || FV_BUILD：公开名只有声明 */

/* ============================ 批量入口的实现 ============================
 * 只在「本头文件提供实现」时编译：
 *   FV_USE_EXTERN —— 整段消失（实现来自 libflashvec.a，头文件保持声明-only）；
 *   FV_BUILD      —— 编译（库需要它），但不含下面的公开名（它们定义在 src/flashvec_batch.cpp）；
 *   默认          —— 编译，再由最后一层 static inline 转发接到公开名上。
 * 实现只有一份：函数体全在 fv_i_* 里，两种形态的差别只是「谁把它接到公开名上」。
 *
 * 为什么分成两段——唯一的分界线是「这个函数里有没有『乘法结果直接进加法』的表达式」：
 *   第一段（无收缩风险，不在 pragma 里）
 *     add / sub / scale / add_mask / add_idx / scale_idx / add_soa / scale_soa：
 *     表达式只有单个加法或单个乘法，FMA 收缩对它们没有意义。放在 pragma 外，消费者用
 *     默认编译选项就能把入口（含 n < 4 / n < 8 的尾巴循环）完整内联进自己的代码。
 *     实测（把等价内核放到 -march=native -mfma -ffp-contract=fast 下）：vfmadd = 0。
 *   第二段（有收缩风险，整体在 pragma 里）
 *     lerp = a + (b-a)*t、madd = a + b*s、len / normalize 的 x*x + y*y、
 *     dot 的 ax*bx + ay*by 累加：实测同样选项下，不带 pragma 的等价内核每处出现 1 条
 *     vfmadd；带上 pragma 后全部为 0。代价是这 5 个入口在消费者默认选项下不内联
 *     （每批一次 TU 内调用）—— 这是有意用「一次调用」换「位级承诺」。
 *
 * 为什么不能只靠 target("avx2")：target 属性只管本函数自己的指令集基线，管不住调用方
 * 命令行带来的 -mfma / -march=native。实测：-march=native -mfma -ffp-contract=fast 下，
 * 带 target("avx2") 的内核照样收缩出 vfmadd。消费者场景下 pragma 是唯一的防线。
 */
#if !defined(FV_USE_EXTERN)

#if FV_I_INTRINSICS
#define FV_I_TARGET_AVX2 __attribute__((target("avx2")))
#endif

/* ======== 第一段：无 mul+add 收缩风险（pragma 之外 → 消费者默认选项下可完整内联） ======== */

/* 标量表达式（非 AVX2 机器与所有尾巴共用）：表达式与 src/flashvec.cpp 的逐元素 API
 * 逐字对应；只有单个加法/乘法，不存在可收缩的 mul+add。 */
static inline fv_vec2 fv_i_s_add(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    return r;
}


static inline fv_vec2 fv_i_s_sub(fv_vec2 a, fv_vec2 b)
{
    fv_vec2 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    return r;
}


static inline fv_vec2 fv_i_s_scale(fv_vec2 v, fv_real s)
{
    fv_vec2 r;
    r.x = v.x * s;
    r.y = v.y * s;
    return r;
}

static inline fv_vec4 fv_i_s_add4(fv_vec4 a, fv_vec4 b)
{
    fv_vec4 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    r.w = a.w + b.w;
    return r;
}


static inline fv_vec4 fv_i_s_sub4(fv_vec4 a, fv_vec4 b)
{
    fv_vec4 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
    r.w = a.w - b.w;
    return r;
}


static inline fv_vec4 fv_i_s_scale4(fv_vec4 v, fv_real s)
{
    fv_vec4 r;
    r.x = v.x * s;
    r.y = v.y * s;
    r.z = v.z * s;
    r.w = v.w * s;
    return r;
}

#if FV_I_INTRINSICS

/* AVX2 内核：只做 add / sub / mul 三种单操作，没有可收缩的 mul+add。
 * 4 个 fv_vec2 = 32 B = 一个 ymm：loadu/storeu 直接吃连续元素，零重排。
 * fv_vec2 只有 4 B 对齐，所以全程用 loadu/storeu（不能用 _mm256_load_ps）。 */

/* 4 个 fv_vec2 = 32 B = 一个 ymm：loadu/storeu 直接吃连续元素，零重排。
 * fv_vec2 只有 4 B 对齐，所以全程用 loadu/storeu（不能用 _mm256_load_ps）。
 * 目标里只写 avx2、不写 fma：FMA 指令根本不可用，收缩在指令层面就被排除。 */

FV_I_TARGET_AVX2 static inline void fv_i_add4_avx2(const fv_vec2 *__restrict a,
                                                   const fv_vec2 *__restrict b,
                                                   fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_sub4_avx2(const fv_vec2 *__restrict a,
                                                   const fv_vec2 *__restrict b,
                                                   fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_sub_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_scale4_avx2(const fv_vec2 *__restrict a, fv_real s,
                                                     fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_mul_ps(va, vs));
    }
}


/* 掩码变体：4 个元素一次（一个 ymm 装 4 个 fv_vec2 = 8 个 float lane）。
 * 关键点：掩码是按**元素**给的，而一个元素占两个 float lane，所以每个掩码字节
 * 必须铺成两个 32-bit lane（x 与 y 共用同一个掩码位），不能按字节平铺。
 * 之后 blendvps 选 a（掩码为 0）或 a+b（掩码非 0），整条写回 o。 */
FV_I_TARGET_AVX2 static inline void fv_i_add_mask4_avx2(const fv_vec2 *__restrict a,
                                                        const fv_vec2 *__restrict b,
                                                        fv_vec2 *__restrict o,
                                                        const uint8_t *__restrict mask, size_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        memcpy(&mb, mask + i, sizeof(mb));
        if (mb == 0u) {
            /* 块级预扫：四个掩码字节全零 → 整块输出就是 a 的原值（掩码外写回 a[i] 的语义），
             * 于是这一块不需要 b 的 load，也不需要 blend。
             * 判据必须是「整个 u32 为零」，不能用 (mb & 0x01010101u) == 0：那只判每个字节的
             * 最低位，掩码字节为 0x02 / 0x80 / 0xFE 这类**偶数或高位非零**的值时会被误判成
             * 「没命中」，而掩码语义是「字节非零就命中」。u32 的零值判断与字节序无关。 */
            _mm256_storeu_ps((float *)(void *)(o + i),
                             _mm256_loadu_ps((const float *)(const void *)(a + i)));
            continue;
        }
        __m128i b4 = _mm_cvtsi32_si128((int)mb);  /* 低 4 字节 = 4 个掩码字节 */
        __m128i w = _mm_cvtepu8_epi16(b4);        /* 4 个字节 -> 4 个 word（0..255） */
        __m128i wd = _mm_unpacklo_epi16(w, w);    /* 每个 word 复制一份 */
        __m256i d = _mm256_cvtepi16_epi32(wd);    /* 8 个 dword：m0 m0 m1 m1 m2 m2 m3 m3 */
        __m256 m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(d, zero));
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_blendv_ps(va, _mm256_add_ps(va, vb), m));
    }
}


/* SoA：8 个 float 一次。 */
FV_I_TARGET_AVX2 static inline void fv_i_add_soa_avx2(const float *__restrict ax,
                                                      const float *__restrict ay,
                                                      const float *__restrict bx,
                                                      const float *__restrict by,
                                                      float *__restrict ox,
                                                      float *__restrict oy, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_add_ps(_mm256_loadu_ps(ax + i),
                                               _mm256_loadu_ps(bx + i)));
        _mm256_storeu_ps(oy + i, _mm256_add_ps(_mm256_loadu_ps(ay + i),
                                               _mm256_loadu_ps(by + i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_scale_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay, fv_real s,
                                                        float *__restrict ox,
                                                        float *__restrict oy, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_mul_ps(_mm256_loadu_ps(ax + i), vs));
        _mm256_storeu_ps(oy + i, _mm256_mul_ps(_mm256_loadu_ps(ay + i), vs));
    }
}

/* ---------------- 流式（NT store）内核：第一段 ----------------
 * 三段结构，对任意指针对齐都正确：
 *   头部 peel —— 用**普通 store** 逐元素推到 32 B 边界（这样任意指针对齐都能进主体）；
 *   主体     —— _mm256_stream_ps（不经 cache、不产生 RFO 的写组合存储），要求 32 B 对齐；
 *   尾巴     —— 不足一个向量的部分继续用普通 store。
 * _mm256_stream_ps 的地址必须 32 B 对齐，否则 #GP；所以主体之前的那次对齐判断是硬约束，
 * 不是优化。输入侧仍然 loadu：只有**输出**需要对齐。
 * 结尾 _mm_sfence()：NT store 是弱序（WC buffer 里的数据对其它核可能还不可见），
 * 不 sfence 就返回会让调用方随后的读线程看到旧值。 */
FV_I_TARGET_AVX2 static inline void fv_i_add4_stream_avx2(const fv_vec2 *__restrict a,
                                                          const fv_vec2 *__restrict b,
                                                          fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    /* 头部：8 B 步长下把 32 B 对齐最多推 3 步 */
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        /* 32 B 边界不可达：只可能发生在 o 不是 8 B 对齐的 fv_vec2 数组上（步长 8 B 无法
         * 消除 4 B 的相位差）。整段交给普通 storeu 内核 —— 头部那几个元素被重算，
         * 结果逐位相同，只是没有 NT 的收益。 */
        fv_i_add4_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_add_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_sub4_stream_avx2(const fv_vec2 *__restrict a,
                                                          const fv_vec2 *__restrict b,
                                                          fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_sub4_avx2(a, b, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_sub_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)),
                                       _mm256_loadu_ps((const float *)(const void *)(b + i))));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_scale4_stream_avx2(const fv_vec2 *__restrict a, fv_real s,
                                                            fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 3u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        fv_i_scale4_avx2(a, s, o, n);
        _mm_sfence();
        return;
    }
    for (; i + 4 <= n; i += 4) {
        _mm256_stream_ps((float *)(void *)(o + i),
                         _mm256_mul_ps(_mm256_loadu_ps((const float *)(const void *)(a + i)), vs));
    }
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
    _mm_sfence();
}


/* SoA 的流式内核按「一条输出分量数组」写：_soa 的流式入口各调用它两次（x 一次、y 一次）。
 * 这样 ox 与 oy 各自独立 peel、各自独立决定能不能 NT，不会因为两块输出的相位差不同而
 * 互相拖累（共用同一个 i 时，两者必须能**同时**对齐）。 */
FV_I_TARGET_AVX2 static inline void fv_i_add1_stream_avx2(const float *__restrict a,
                                                          const float *__restrict b,
                                                          float *__restrict o, size_t n)
{
    size_t i = 0;
    /* float 步长 4 B：7 步就能覆盖 32 B 内的全部相位，正常指针对齐必然可达 */
    for (; i < n && i < 7u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = a[i] + b[i];
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        /* 不可达：只可能是 o 连 4 B 对齐都不满足的病态输入 */
        for (; i + 8 <= n; i += 8) {
            _mm256_storeu_ps(o + i, _mm256_add_ps(_mm256_loadu_ps(a + i),
                                                  _mm256_loadu_ps(b + i)));
        }
    } else {
        for (; i + 8 <= n; i += 8) {
            _mm256_stream_ps(o + i, _mm256_add_ps(_mm256_loadu_ps(a + i),
                                                  _mm256_loadu_ps(b + i)));
        }
    }
    for (; i < n; ++i) {
        o[i] = a[i] + b[i];
    }
    _mm_sfence();
}


FV_I_TARGET_AVX2 static inline void fv_i_scale1_stream_avx2(const float *__restrict a, fv_real s,
                                                            float *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i < n && i < 7u && (((uintptr_t)(const void *)(o + i)) & 31u) != 0u; ++i) {
        o[i] = a[i] * s;
    }
    if ((((uintptr_t)(const void *)(o + i)) & 31u) != 0u) {
        for (; i + 8 <= n; i += 8) {
            _mm256_storeu_ps(o + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), vs));
        }
    } else {
        for (; i + 8 <= n; i += 8) {
            _mm256_stream_ps(o + i, _mm256_mul_ps(_mm256_loadu_ps(a + i), vs));
        }
    }
    for (; i < n; ++i) {
        o[i] = a[i] * s;
    }
    _mm_sfence();
}


/* 就地稀疏加的两条内核 —— 结果**逐位相同**（命中 lane 都是标量 a[i] + b[i]），
 * 差别只在「未命中的 lane 要不要发生访存」。选路由 fv_i_mask_inplace_use_dense() 每批决定一次。
 *
 * 稀疏支（采样密度 <= FV_I_MIP_DENSE_NUM/DEN）：块级预扫 + 逐元素标量。
 *   整块四个掩码字节全零就整块跳过 —— 连 a 的 cache 行都不碰（5% 命中时 4 元素块有 95%^4
 *   ≈ 81% 是全零块），有命中的块里只对命中的元素做 8 B 读改写。
 *   这里不设「四块全命中就走向量」的快路径：在 <= 1/8 的密度下全命中块的概率约 0.02%，
 *   为它加一条判断不值。
 *   **注意**：两条支之间「哪条更划算」没有实测支撑 —— 本入口相对 _mask 是**全档 INCONCLUSIVE**
 *   （增益 -0.10…+0.00，全在同轮 A/A 地板以下，见下面的能力边界）。所以选路只影响代码形态，
 *   两条路径**无可测差异**；上面那句只是拒绝为 0.02% 的事件加一条分支，不是性能论断。 */
FV_I_TARGET_AVX2 static inline void fv_i_branch_add4_prescan_avx2(fv_vec2 *__restrict a,
                                                                  const fv_vec2 *__restrict b,
                                                                  const uint8_t *__restrict mask,
                                                                  size_t n)
{
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        size_t k;
        memcpy(&mb, mask + i, sizeof(mb));
        if (mb == 0u) {
            continue; /* 四个元素都没命中：a / b 的这两条 cache 行整个不碰 */
        }
        for (k = 0; k < 4u; ++k) {
            if (mask[i + k] != 0u) {
                a[i + k] = fv_i_s_add(a[i + k], b[i + k]);
            }
        }
    }
}


/* 密集支（采样密度 > FV_I_MIP_DENSE_NUM/DEN）：整块 32 B load a / load b，blend 出结果，
 * 整条 storeu 写回 a。没有数据相关分支，也没有 VMASKMOVPS（那一族在 Zen 2 上是微码级慢指令；
 * **早期探针**测到同密度下 maskstore 比这里的整条 storeu 慢 1.2-3.8x，本轮加硬夹具**未复测**，
 * 只能按「量级」读 —— 这一族已经从库里删掉，这句只解释「为什么不用它」）。
 * **未命中的 lane 会被写回 a 的原值**（blend 选 va）：值语义上完全等价，但物理上确实写了。
 * 因此本函数只承诺「未命中元素的值不变」，不承诺「未被物理写回」——见契约。 */
FV_I_TARGET_AVX2 static inline void fv_i_blend_add4_inplace_avx2(fv_vec2 *__restrict a,
                                                                 const fv_vec2 *__restrict b,
                                                                 const uint8_t *__restrict mask,
                                                                 size_t n)
{
    const __m256i zero = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        uint32_t mb;
        __m128i b4;
        __m128i w;
        __m128i wd;
        __m256i d;
        __m256 m;
        __m256 va;
        __m256 vb;
        memcpy(&mb, mask + i, sizeof(mb));
        b4 = _mm_cvtsi32_si128((int)mb);
        w = _mm_cvtepu8_epi16(b4);
        wd = _mm_unpacklo_epi16(w, w);
        d = _mm256_cvtepi16_epi32(wd);
        m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(d, zero));
        va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(a + i),
                         _mm256_blendv_ps(va, _mm256_add_ps(va, vb), m));
    }
}


/* ---- fv_vec4 的 AoS 内核（_x2）：2 个 16 B 元素 = 32 B = 一个 ymm ----
 * 这里**不存在**「把两个 16 B 元素拼成一个 ymm」这回事：fv_vec4 没有填充字节，
 * 数组里相邻两个元素本来就是 32 B 连续内存，一条 256-bit loadu 读进来的 8 个 lane 就是
 * (a[i].x, a[i].y, a[i].z, a[i].w, a[i+1].x, ...)，顺序天生正确；输出同理一条 256-bit storeu。
 * 所以 _x2 与 vec2 的 _x4 一样是**零重排**（不需要 insertf128 之类的拼装指令）。
 * 两条 128-bit load + _mm256_insertf128_ps 的写法，只在「元素步长不是 16 B」的布局里才需要
 * （AoSoA、带 stride 的数组），那种情况下才谈得上拼装成本。
 * 元素 16 B、ymm 32 B：正好 2 个元素装满一个寄存器，**一个 lane 都不浪费** ——
 * 这就是 _x2 这个名字的含义（与 vec2 的 _x4 同构：数字 = 一个 256-bit 寄存器装几个元素）。 */
FV_I_TARGET_AVX2 static inline void fv_i_v4_add2_avx2(const fv_vec4 *__restrict a,
                                                    const fv_vec4 *__restrict b,
                                                    fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_sub2_avx2(const fv_vec4 *__restrict a,
                                                    const fv_vec4 *__restrict b,
                                                    fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_sub_ps(va, vb));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_scale2_avx2(const fv_vec4 *__restrict a, fv_real s,
                                                      fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        _mm256_storeu_ps((float *)(void *)(o + i), _mm256_mul_ps(va, vs));
    }
}


/* ---- 三维 SoA 内核：每个分量一条独立的 8 lane 循环 ----
 * 「3 个分量跑 8 lane 会不会浪费 lane」——不会：一条 ymm 里装的永远是**同一分量的 8 个
 * 不同元素**（8 个 ax[i..i+7]），不是「一个元素的 3 个分量 + 1 个填充 lane」。
 * 每个分量的 lane 利用率都是 100%，代价只落在调用方：数据得按分量分开存，
 * 且 n < 8 的尾巴（最多 7 个元素）走标量。AoS 的 12 B 布局则相反：装不进 ymm 的整数倍，
 * 要么空转 1/4 lane、要么 shuffle 拆包，所以三维批量只提供 SoA。 */
FV_I_TARGET_AVX2 static inline void fv_i_v3_add_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       const float *__restrict bx,
                                                       const float *__restrict by,
                                                       const float *__restrict bz,
                                                       float *__restrict ox,
                                                       float *__restrict oy,
                                                       float *__restrict oz, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_add_ps(_mm256_loadu_ps(ax + i),
                                               _mm256_loadu_ps(bx + i)));
        _mm256_storeu_ps(oy + i, _mm256_add_ps(_mm256_loadu_ps(ay + i),
                                               _mm256_loadu_ps(by + i)));
        _mm256_storeu_ps(oz + i, _mm256_add_ps(_mm256_loadu_ps(az + i),
                                               _mm256_loadu_ps(bz + i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v3_scale_soa_avx2(const float *__restrict ax,
                                                         const float *__restrict ay,
                                                         const float *__restrict az, fv_real s,
                                                         float *__restrict ox,
                                                         float *__restrict oy,
                                                         float *__restrict oz, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(ox + i, _mm256_mul_ps(_mm256_loadu_ps(ax + i), vs));
        _mm256_storeu_ps(oy + i, _mm256_mul_ps(_mm256_loadu_ps(ay + i), vs));
        _mm256_storeu_ps(oz + i, _mm256_mul_ps(_mm256_loadu_ps(az + i), vs));
    }
}


#endif /* FV_I_INTRINSICS（第一段的内核） */

/* 第一段的 8 个入口：运行时派发 + 尾巴。
 * 整段没有可收缩的浮点表达式，所以能完整内联进调用方；两个 _idx 变体故意走标量
 * （随机访存下 gather 在本机 Zen 2 上是微码循环，通常不如逐元素；语义上 idx 可重复、
 * 可乱序，重复下标的写入值相同，与处理顺序无关）。
 * 内核带 target("avx2")、与调用方目标不一致，所以批量部分仍是每批一次 TU 内调用。 */
static inline void fv_i_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                               fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                               fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_sub4_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                                 fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_scale4_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                    fv_vec2 *__restrict o, const uint8_t *__restrict mask,
                                    size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add_mask4_avx2(a, b, o, mask, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        if (mask[i] != 0) {
            o[i] = fv_i_s_add(a[i], b[i]);
        } else {
            o[i] = a[i]; /* 掩码外写回 a 的原值（不是「不写」） */
        }
    }
}


static inline void fv_i_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                   size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                                     fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                     size_t m)
{
    size_t k = 0;
    for (; k < m; ++k) {
        const uint32_t i = idx[k];
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_soa(const float *__restrict ax, const float *__restrict ay,
                                const float *__restrict bx, const float *__restrict by,
                                float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_add_soa_avx2(ax, ay, bx, by, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}


static inline void fv_i_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                  fv_real s, float *__restrict ox, float *__restrict oy,
                                  size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_scale_soa_avx2(ax, ay, s, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}

/* ---- 第一段的流式入口：只有存储密集的运算有 stream 变体 ----
 * lerp / madd 没有 stream 变体：它们的结果一般马上被下一步消费，盯着 cache 里的副本更划算，
 * 而且 NT store 的收益来自「写出去的这一大片不会再被读」，那是 add/sub/scale/normalize
 * 这类批处理的特征，不是插值中间结果的特征。
 * 入口形状与对应的非 stream 入口逐字相同（bulk + 标量尾巴），差别只在调哪个内核；
 * 尾长 / n == 0 的语义因此与非 stream 版完全一致（n == 0 时不碰任何指针，也不需要 sfence）。 */
static inline void fv_i_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_add4_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add(a[i], b[i]);
    }
}


static inline void fv_i_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_sub4_stream_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub(a[i], b[i]);
    }
}


static inline void fv_i_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                                        fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_scale4_stream_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale(a[i], s);
    }
}


static inline void fv_i_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict bx, const float *__restrict by,
                                       float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_add1_stream_avx2(ax, bx, ox, bulk);
            fv_i_add1_stream_avx2(ay, by, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}


static inline void fv_i_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                         fv_real s, float *__restrict ox, float *__restrict oy,
                                         size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_scale1_stream_avx2(ax, s, ox, bulk);
            fv_i_scale1_stream_avx2(ay, s, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}


/* 就地稀疏加的选路：低密度走「块级预扫 + 逐元素」（不碰未命中的块），
 * 中高密度走「整块 load + blend + storeu」（没有数据相关分支）。
 * 判据只采样前 64 个掩码字节（成本几十个周期，相对整批可忽略），只要量级正确；
 * 两条路径的结果逐位相同，估算偏了只损失性能、不影响语义。
 * 分界点由下面三个宏决定（改宏即可调，不必碰函数体）：采样长度与密度阈值。 */
#define FV_I_MIP_DENSE_NUM 1u /* 采样密度 > NUM/DEN 时走 dense 路径 */
#define FV_I_MIP_DENSE_DEN 8u
#define FV_I_MIP_SAMPLE 64u

static inline int fv_i_mask_inplace_use_dense(const uint8_t *mask, size_t n)
{
    const size_t samp = (n < FV_I_MIP_SAMPLE) ? n : FV_I_MIP_SAMPLE;
    size_t hits = 0;
    size_t k;
    for (k = 0; k < samp; ++k) {
        hits += (mask[k] != 0u) ? (size_t)1 : (size_t)0;
    }
    return (hits * FV_I_MIP_DENSE_DEN > samp * FV_I_MIP_DENSE_NUM) ? 1 : 0;
}


/* 就地稀疏加：mask[i] != 0 时 a[i] += b[i]，否则 a[i] 的值保持不变。
 * 没有任何可收缩的 mul+add，所以留在第一段。
 * 实现按采样密度在两条内核之间二选一（见上面的 fv_i_mask_inplace_use_dense）：
 * 结果与选哪条路径无关 —— 命中 lane 都逐位等于标量 a[i] + b[i]。
 * 注意两段的**存储行为**不同：低密度不碰未命中块，高密度会把未命中的 lane 写回原值。
 * 这是值语义的承诺，不是存储行为的承诺（见契约「 _mask_inplace 」小节）。 */
static inline void fv_i_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                            const uint8_t *__restrict mask, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            if (fv_i_mask_inplace_use_dense(mask, bulk) != 0) {
                fv_i_blend_add4_inplace_avx2(a, b, mask, bulk);
            } else {
                fv_i_branch_add4_prescan_avx2(a, b, mask, bulk);
            }
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        if (mask[i] != 0) {
            a[i] = fv_i_s_add(a[i], b[i]);
        }
    }
}

/* ---- 三维 / 四维的第一段入口 ----
 * 与上面 vec2 的入口同构：bulk 段交给带 target("avx2") 的内核，剩下的尾巴走标量表达式。
 * 这里没有任何「乘法结果直接进加法」的表达式（只有单个 +、-、×），FMA 收缩对它们无意义，
 * 所以**故意留在 pragma 之外**：消费者用默认编译选项就能把整段（含尾巴循环）完整内联。
 * 尾长：_x2 的 bulk 是 n & ~1（尾巴 0-1 个），_soa 的 bulk 是 n & ~7（尾巴 0-7 个）。 */
static inline void fv_i_v4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_add2_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_add4(a[i], b[i]);
    }
}


static inline void fv_i_v4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_sub2_avx2(a, b, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_sub4(a[i], b[i]);
    }
}


static inline void fv_i_v4_scale_x2(const fv_vec4 *__restrict a, fv_real s, fv_vec4 *__restrict o,
                                  size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_scale2_avx2(a, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_scale4(a[i], s);
    }
}


static inline void fv_i_v3_add_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az,
                                 const float *__restrict bx, const float *__restrict by,
                                 const float *__restrict bz,
                                 float *__restrict ox, float *__restrict oy,
                                 float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_add_soa_avx2(ax, ay, az, bx, by, bz, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
        oz[i] = az[i] + bz[i];
    }
}


static inline void fv_i_v3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, fv_real s,
                                   float *__restrict ox, float *__restrict oy,
                                   float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_scale_soa_avx2(ax, ay, az, s, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
        oz[i] = az[i] * s;
    }
}


/* ======== 第二段：有 mul+add 收缩风险（pragma 之内 → 每批一次 TU 内调用） ======== */

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("fp-contract=off")
#endif


static inline fv_vec2 fv_i_s_lerp(fv_vec2 a, fv_vec2 b, fv_real t)
{
    fv_vec2 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    return r;
}


static inline fv_vec2 fv_i_s_madd(fv_vec2 a, fv_vec2 b, fv_real s)
{
    fv_vec2 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    return r;
}


/* 与 fv_vec2_len 同表达式：sqrtf(x*x + y*y)。 */
static inline fv_real fv_i_s_len(fv_real x, fv_real y)
{
    return sqrtf(x * x + y * y);
}


/* 与 fv_vec2_normalize 同语义、同表达式、同分支。 */
static inline void fv_i_s_normalize(fv_real x, fv_real y, fv_real *ox, fv_real *oy)
{
    fv_real len = fv_i_s_len(x, y);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
    }
}

static inline fv_vec4 fv_i_s_lerp4(fv_vec4 a, fv_vec4 b, fv_real t)
{
    fv_vec4 r;
    r.x = a.x + (b.x - a.x) * t;
    r.y = a.y + (b.y - a.y) * t;
    r.z = a.z + (b.z - a.z) * t;
    r.w = a.w + (b.w - a.w) * t;
    return r;
}


static inline fv_vec4 fv_i_s_madd4(fv_vec4 a, fv_vec4 b, fv_real s)
{
    fv_vec4 r;
    r.x = a.x + b.x * s;
    r.y = a.y + b.y * s;
    r.z = a.z + b.z * s;
    r.w = a.w + b.w * s;
    return r;
}


/* 三维长度：sqrtf((x*x + y*y) + z*z) —— 与 fv_vec3_len_sqr（= fv_vec3_dot(v,v)）的左结合顺序一致。
 * 括号写全是为了把结合顺序固定在源码里，不是给编译器看的。 */
static inline fv_real fv_i_s_len3(fv_real x, fv_real y, fv_real z)
{
    return sqrtf((x * x + y * y) + z * z);
}


static inline void fv_i_s_normalize3(fv_real x, fv_real y, fv_real z,
                                     fv_real *ox, fv_real *oy, fv_real *oz)
{
    fv_real len = fv_i_s_len3(x, y, z);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        *oz = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
        *oz = z * inv;
    }
}


/* 四维长度：sqrtf(((x*x + y*y) + z*z) + w*w) —— 与 fv_vec4_len_sqr 的左结合顺序一致。 */
static inline fv_real fv_i_s_len4(fv_real x, fv_real y, fv_real z, fv_real w)
{
    return sqrtf(((x * x + y * y) + z * z) + w * w);
}


static inline void fv_i_s_normalize4(fv_real x, fv_real y, fv_real z, fv_real w,
                                     fv_real *ox, fv_real *oy, fv_real *oz, fv_real *ow)
{
    fv_real len = fv_i_s_len4(x, y, z, w);
    if (len <= FV_EPS) {
        *ox = (fv_real)0;
        *oy = (fv_real)0;
        *oz = (fv_real)0;
        *ow = (fv_real)0;
        return;
    }
    {
        fv_real inv = (fv_real)1 / len;
        *ox = x * inv;
        *oy = y * inv;
        *oz = z * inv;
        *ow = w * inv;
    }
}

#if FV_I_INTRINSICS

/* lerp = a + (b - a) * t，乘加顺序与标量一致；fp-contract=off 保证不会变成 FMA。 */
FV_I_TARGET_AVX2 static inline void fv_i_lerp4_avx2(const fv_vec2 *__restrict a,
                                                    const fv_vec2 *__restrict b, fv_real t,
                                                    fv_vec2 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
}


/* madd = a + b * s，同上。 */
FV_I_TARGET_AVX2 static inline void fv_i_madd4_avx2(const fv_vec2 *__restrict a,
                                                    const fv_vec2 *__restrict b, fv_real s,
                                                    fv_vec2 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_len_soa_avx2(const float *__restrict ax,
                                                      const float *__restrict ay,
                                                      float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(_mm256_add_ps(
                                    _mm256_mul_ps(x, x), _mm256_mul_ps(y, y))));
    }
}


/* 归一化：零向量保护写成掩码而不是分支。
 * keep = !(len <= FV_EPS)，与标量 if (len <= FV_EPS) 严格互补（NaN 时走 else 分支，
 * 用 _CMP_NLE_UQ 才能对上：unordered 也算「不满足 <=」）。
 * 除法的被除数用 blend 换成 1.0，所以零向量 lane 不产生 inf/NaN，也不置起除零标志；
 * 最后 and(keep) 把该 lane 清成 +0.0，与标量返回 (0,0) 的位模式一致。 */
FV_I_TARGET_AVX2 static inline void fv_i_normalize_soa_avx2(const float *__restrict ax,
                                                            const float *__restrict ay,
                                                            float *__restrict ox,
                                                            float *__restrict oy, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
    }
}


/* dot：4 个独立累加器打断依赖链（早期微基准：单累加器 0.3213 -> 4 累加器 0.1428 ns/元素，
 * 本轮加硬夹具未复测；这条只说「依赖链是瓶颈」，不当作 dot_soa 的增益承诺），
 * 末端一次归约。x 分量与 y 分量累加进同一个累加器，所以归约顺序与逐元素路径不同——
 * 这是 fv_vec2_dot_soa 唯一允许的偏差（头文件已写明）。 */
FV_I_TARGET_AVX2 static inline float fv_i_dot_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict bx,
                                                       const float *__restrict by, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    __m256 s2 = _mm256_setzero_ps();
    __m256 s3 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 32 <= n; i += 32) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 0), _mm256_loadu_ps(bx + i + 0)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 0), _mm256_loadu_ps(by + i + 0))));
        s1 = _mm256_add_ps(s1, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 8), _mm256_loadu_ps(bx + i + 8)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 8), _mm256_loadu_ps(by + i + 8))));
        s2 = _mm256_add_ps(s2, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 16), _mm256_loadu_ps(bx + i + 16)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 16), _mm256_loadu_ps(by + i + 16))));
        s3 = _mm256_add_ps(s3, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i + 24), _mm256_loadu_ps(bx + i + 24)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i + 24), _mm256_loadu_ps(by + i + 24))));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, _mm256_add_ps(
                 _mm256_mul_ps(_mm256_loadu_ps(ax + i), _mm256_loadu_ps(bx + i)),
                 _mm256_mul_ps(_mm256_loadu_ps(ay + i), _mm256_loadu_ps(by + i))));
    }

    /* 末端归约：先把 4 个 ymm 累加器两两相加，再做一次 4 lane 的树。
     * 用 movehdup + shuffle 而不是 hadd（本机 hadd 慢约 2 倍）。
     * q = [q0,q1,q2,q3] -> r = [q0+q1, ..., q2+q3, ...] -> 取 lane0 + lane2。 */
    {
        __m256 t = _mm256_add_ps(_mm256_add_ps(s0, s1), _mm256_add_ps(s2, s3));
        __m128 lo = _mm256_castps256_ps128(t);
        __m128 hi = _mm256_extractf128_ps(t, 1);
        __m128 q = _mm_add_ps(lo, hi);
        __m128 r = _mm_add_ps(q, _mm_movehdup_ps(q));
        sum = _mm_cvtss_f32(_mm_add_ss(r, _mm_shuffle_ps(r, r, 0x02)));
    }

    /* 尾巴按顺序累加（n < 8 时整段都走这里，结果与逐元素路径逐位一致）。 */
    for (; i < n; ++i) {
        sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
    }
    return sum;
}
/* normalize 的流式内核：算术路径与 fv_i_normalize_soa_avx2 逐字相同（含零向量掩码保护），
 * 只把两条输出换成 NT store。ox 与 oy 共用同一个 i，所以 peel 的目标是「两者同时」到 32 B
 * 边界；推不到（两块输出在 32 B 内的相位差不是 8 的倍数）时退回普通 storeu，结果逐位相同。
 * 建议把 ox / oy 都按 32 B 对齐分配，那样一次都不用退。 */
FV_I_TARGET_AVX2 static inline void fv_i_normalize_stream_soa_avx2(const float *__restrict ax,
                                                                  const float *__restrict ay,
                                                                  float *__restrict ox,
                                                                  float *__restrict oy, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    int nt;
    size_t i = 0;

    for (; i < n && i < 7u
           && ((((uintptr_t)(const void *)(ox + i)) & 31u) != 0u
               || (((uintptr_t)(const void *)(oy + i)) & 31u) != 0u); ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
    nt = ((((uintptr_t)(const void *)(ox + i)) & 31u) == 0u
          && (((uintptr_t)(const void *)(oy + i)) & 31u) == 0u) ? 1 : 0;

    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 len = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(x, x),
                                                  _mm256_mul_ps(y, y)));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        __m256 rx = _mm256_and_ps(_mm256_mul_ps(x, inv), keep);
        __m256 ry = _mm256_and_ps(_mm256_mul_ps(y, inv), keep);
        if (nt != 0) {
            _mm256_stream_ps(ox + i, rx);
            _mm256_stream_ps(oy + i, ry);
        } else {
            _mm256_storeu_ps(ox + i, rx);
            _mm256_storeu_ps(oy + i, ry);
        }
    }
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
    _mm_sfence();
}


/* ---- 四维 AoS 的 lerp / madd 内核：2 个元素一次，零重排（同 _x2 的 add/sub/scale） ---- */
FV_I_TARGET_AVX2 static inline void fv_i_v4_lerp2_avx2(const fv_vec4 *__restrict a,
                                                     const fv_vec4 *__restrict b, fv_real t,
                                                     fv_vec4 *__restrict o, size_t n)
{
    const __m256 vt = _mm256_set1_ps(t);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        __m256 vd = _mm256_sub_ps(vb, va);
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vd, vt)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_madd2_avx2(const fv_vec4 *__restrict a,
                                                     const fv_vec4 *__restrict b, fv_real s,
                                                     fv_vec4 *__restrict o, size_t n)
{
    const __m256 vs = _mm256_set1_ps(s);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m256 va = _mm256_loadu_ps((const float *)(const void *)(a + i));
        __m256 vb = _mm256_loadu_ps((const float *)(const void *)(b + i));
        _mm256_storeu_ps((float *)(void *)(o + i),
                         _mm256_add_ps(va, _mm256_mul_ps(vb, vs)));
    }
}


/* ---- SoA 的三个逐 lane 原语：结合顺序的唯一决定点 ----
 * fv_i_v3_mul8 / fv_i_v4_mul8 算的是「一条 ymm 上 8 个元素的逐 lane 内积项」，结合顺序与标量
 * API 逐字对应：
 *   vec3: (x*x' + y*y') + z*z'
 *   vec4: ((x*x' + y*y') + z*z') + w*w'
 * dot / len / normalize 三个内核全部通过它们取值，所以结合顺序只需要在这里对一次。
 * 参数**故意不带** __restrict：len 的调用点是 fv_i_v3_mul8(ax, ay, az, ax, ay, az, i)，
 * 两边的指针指向同一块内存。
 * 这两个函数在 mul 与 add 之间没有可收缩的空间——mul 的结果直接进 add，所以它们必须在
 * pragma（fp-contract=off）之内；否则消费者的 -mfma -ffp-contract=fast 会把它们变成 vfmadd。 */
FV_I_TARGET_AVX2 static inline __m256 fv_i_v3_mul8(const float *ax, const float *ay,
                                                 const float *az,
                                                 const float *bx, const float *by,
                                                 const float *bz, size_t off)
{
    __m256 px = _mm256_mul_ps(_mm256_loadu_ps(ax + off), _mm256_loadu_ps(bx + off));
    __m256 py = _mm256_mul_ps(_mm256_loadu_ps(ay + off), _mm256_loadu_ps(by + off));
    __m256 pz = _mm256_mul_ps(_mm256_loadu_ps(az + off), _mm256_loadu_ps(bz + off));
    return _mm256_add_ps(_mm256_add_ps(px, py), pz);
}


FV_I_TARGET_AVX2 static inline __m256 fv_i_v4_mul8(const float *ax, const float *ay,
                                                 const float *az, const float *aw,
                                                 const float *bx, const float *by,
                                                 const float *bz, const float *bw, size_t off)
{
    __m256 px = _mm256_mul_ps(_mm256_loadu_ps(ax + off), _mm256_loadu_ps(bx + off));
    __m256 py = _mm256_mul_ps(_mm256_loadu_ps(ay + off), _mm256_loadu_ps(by + off));
    __m256 pz = _mm256_mul_ps(_mm256_loadu_ps(az + off), _mm256_loadu_ps(bz + off));
    __m256 pw = _mm256_mul_ps(_mm256_loadu_ps(aw + off), _mm256_loadu_ps(bw + off));
    return _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(px, py), pz), pw);
}


/* 两个 ymm 累加器 -> 标量。写法与 fv_i_dot_soa_avx2 的末端归约相同（movehdup + shuffle，
 * 不用 hadd；那条注释里记了原因）。归约顺序与逐元素路径不同，所以 dot_soa 只承诺确定性
 * 与小误差，不承诺逐位一致（契约里写明）。 */
FV_I_TARGET_AVX2 static inline float fv_i_v34_reduce2(__m256 s0, __m256 s1)
{
    __m256 t = _mm256_add_ps(s0, s1);
    __m128 lo = _mm256_castps256_ps128(t);
    __m128 hi = _mm256_extractf128_ps(t, 1);
    __m128 q = _mm_add_ps(lo, hi);
    __m128 r = _mm_add_ps(q, _mm_movehdup_ps(q));
    return _mm_cvtss_f32(_mm_add_ss(r, _mm_shuffle_ps(r, r, 0x02)));
}


FV_I_TARGET_AVX2 static inline void fv_i_v3_len_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(fv_i_v3_mul8(ax, ay, az, ax, ay, az, i)));
    }
}


/* 归一化：零向量保护写成掩码而不是分支，语义与 fv_vec3_normalize 相同。
 * keep = !(len <= FV_EPS)（_CMP_NLE_UQ：unordered 也算「不满足 <=」，于是 NaN 走 else 分支，
 * 与标量的 if (len <= FV_EPS) 严格互补）；除法的被除数在零向量 lane 被 blend 换成 1.0，
 * 所以不会产生 inf/NaN、也不会置起 MXCSR 的除零标志；最后 and(keep) 把该 lane 清成 +0.0。 */
FV_I_TARGET_AVX2 static inline void fv_i_v3_normalize_soa_avx2(const float *__restrict ax,
                                                             const float *__restrict ay,
                                                             const float *__restrict az,
                                                             float *__restrict ox,
                                                             float *__restrict oy,
                                                             float *__restrict oz, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v3_mul8(ax, ay, az, ax, ay, az, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
        _mm256_storeu_ps(oz + i, _mm256_and_ps(_mm256_mul_ps(z, inv), keep));
    }
}


/* dot：两个独立累加器打断依赖链（每轮 16 个元素），末端一次归约。
 * 3 个分量每轮有 3 条 vmulps + 2 条 vaddps 可以并行，累加链每轮只有 1 条 vaddps，
 * 两个累加器已经足够让这条链不成为瓶颈（本轮不做计时，故不给实测数字）。
 * 归约顺序与逐元素路径不同 —— 这是 fv_vec3_dot_soa 唯一允许的偏差（头文件契约已写明）。 */
FV_I_TARGET_AVX2 static inline float fv_i_v3_dot_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict bx,
                                                        const float *__restrict by,
                                                        const float *__restrict bz, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        s0 = _mm256_add_ps(s0, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i));
        s1 = _mm256_add_ps(s1, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i + 8));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, fv_i_v3_mul8(ax, ay, az, bx, by, bz, i));
    }

    sum = fv_i_v34_reduce2(s0, s1);

    /* 尾巴按顺序累加，表达式与 fv_vec3_dot 的左结合顺序逐字相同。 */
    for (; i < n; ++i) {
        sum = sum + ((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]);
    }
    return sum;
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_len_soa_avx2(const float *__restrict ax,
                                                       const float *__restrict ay,
                                                       const float *__restrict az,
                                                       const float *__restrict aw,
                                                       float *__restrict o, size_t n)
{
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        _mm256_storeu_ps(o + i, _mm256_sqrt_ps(fv_i_v4_mul8(ax, ay, az, aw, ax, ay, az, aw, i)));
    }
}


FV_I_TARGET_AVX2 static inline void fv_i_v4_normalize_soa_avx2(const float *__restrict ax,
                                                             const float *__restrict ay,
                                                             const float *__restrict az,
                                                             const float *__restrict aw,
                                                             float *__restrict ox,
                                                             float *__restrict oy,
                                                             float *__restrict oz,
                                                             float *__restrict ow, size_t n)
{
    const __m256 eps = _mm256_set1_ps(FV_EPS);
    const __m256 one = _mm256_set1_ps((fv_real)1);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 x = _mm256_loadu_ps(ax + i);
        __m256 y = _mm256_loadu_ps(ay + i);
        __m256 z = _mm256_loadu_ps(az + i);
        __m256 w = _mm256_loadu_ps(aw + i);
        __m256 len = _mm256_sqrt_ps(fv_i_v4_mul8(ax, ay, az, aw, ax, ay, az, aw, i));
        __m256 keep = _mm256_cmp_ps(len, eps, _CMP_NLE_UQ);
        __m256 safe = _mm256_blendv_ps(one, len, keep);
        __m256 inv = _mm256_div_ps(one, safe);
        _mm256_storeu_ps(ox + i, _mm256_and_ps(_mm256_mul_ps(x, inv), keep));
        _mm256_storeu_ps(oy + i, _mm256_and_ps(_mm256_mul_ps(y, inv), keep));
        _mm256_storeu_ps(oz + i, _mm256_and_ps(_mm256_mul_ps(z, inv), keep));
        _mm256_storeu_ps(ow + i, _mm256_and_ps(_mm256_mul_ps(w, inv), keep));
    }
}


FV_I_TARGET_AVX2 static inline float fv_i_v4_dot_soa_avx2(const float *__restrict ax,
                                                        const float *__restrict ay,
                                                        const float *__restrict az,
                                                        const float *__restrict aw,
                                                        const float *__restrict bx,
                                                        const float *__restrict by,
                                                        const float *__restrict bz,
                                                        const float *__restrict bw, size_t n)
{
    __m256 s0 = _mm256_setzero_ps();
    __m256 s1 = _mm256_setzero_ps();
    float sum;
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        s0 = _mm256_add_ps(s0, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i));
        s1 = _mm256_add_ps(s1, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i + 8));
    }
    for (; i + 8 <= n; i += 8) {
        s0 = _mm256_add_ps(s0, fv_i_v4_mul8(ax, ay, az, aw, bx, by, bz, bw, i));
    }

    sum = fv_i_v34_reduce2(s0, s1);

    for (; i < n; ++i) {
        sum = sum + (((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]) + aw[i] * bw[i]);
    }
    return sum;
}


#endif /* FV_I_INTRINSICS（第二段的内核） */


static inline void fv_i_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                fv_real t, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_lerp4_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp(a[i], b[i], t);
    }
}


static inline void fv_i_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                fv_real s, fv_vec2 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)3;
        if (bulk != 0) {
            fv_i_madd4_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd(a[i], b[i], s);
    }
}


static inline float fv_i_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict bx, const float *__restrict by,
                                 size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_dot_soa_avx2(ax, ay, bx, by, n);
    }
#endif
    /* 没有 AVX2：单累加器顺序求和，与逐元素路径逐位一致。 */
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + (ax[i] * bx[i] + ay[i] * by[i]);
        }
        return sum;
    }
}


static inline void fv_i_len_soa(const float *__restrict ax, const float *__restrict ay,
                                float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_len_soa_avx2(ax, ay, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len(ax[i], ay[i]);
    }
}


static inline void fv_i_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                      float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_normalize_soa_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
}


/* 流式入口必须留在第二段：它的标量尾巴里有 x*x + y*y（via fv_i_s_normalize），
 * 换到 pragma 外面就会被调用方的 -ffp-contract=fast 收缩成 FMA。 */
static inline void fv_i_normalize_soa_stream(const float *__restrict ax,
                                             const float *__restrict ay,
                                             float *__restrict ox, float *__restrict oy, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_normalize_stream_soa_avx2(ax, ay, ox, oy, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize(ax[i], ay[i], &ox[i], &oy[i]);
    }
}

static inline void fv_i_v4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                 fv_real t, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_lerp2_avx2(a, b, t, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_lerp4(a[i], b[i], t);
    }
}


static inline void fv_i_v4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                 fv_real s, fv_vec4 *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)1;
        if (bulk != 0) {
            fv_i_v4_madd2_avx2(a, b, s, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_madd4(a[i], b[i], s);
    }
}


static inline float fv_i_v3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az,
                                  const float *__restrict bx, const float *__restrict by,
                                  const float *__restrict bz, size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_v3_dot_soa_avx2(ax, ay, az, bx, by, bz, n);
    }
#endif
    /* 没有 AVX2：单累加器顺序求和，与逐元素路径逐位一致（+0.0/-0.0 的边界见契约）。 */
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + ((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]);
        }
        return sum;
    }
}


static inline void fv_i_v3_len_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az, float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_len_soa_avx2(ax, ay, az, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len3(ax[i], ay[i], az[i]);
    }
}


static inline void fv_i_v3_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict az,
                                       float *__restrict ox, float *__restrict oy,
                                       float *__restrict oz, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v3_normalize_soa_avx2(ax, ay, az, ox, oy, oz, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize3(ax[i], ay[i], az[i], &ox[i], &oy[i], &oz[i]);
    }
}


static inline float fv_i_v4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                  const float *__restrict az, const float *__restrict aw,
                                  const float *__restrict bx, const float *__restrict by,
                                  const float *__restrict bz, const float *__restrict bw,
                                  size_t n)
{
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        return fv_i_v4_dot_soa_avx2(ax, ay, az, aw, bx, by, bz, bw, n);
    }
#endif
    {
        fv_real sum = (fv_real)0;
        size_t i = 0;
        for (; i < n; ++i) {
            sum = sum + (((ax[i] * bx[i] + ay[i] * by[i]) + az[i] * bz[i]) + aw[i] * bw[i]);
        }
        return sum;
    }
}


static inline void fv_i_v4_len_soa(const float *__restrict ax, const float *__restrict ay,
                                 const float *__restrict az, const float *__restrict aw,
                                 float *__restrict o, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v4_len_soa_avx2(ax, ay, az, aw, o, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        o[i] = fv_i_s_len4(ax[i], ay[i], az[i], aw[i]);
    }
}


static inline void fv_i_v4_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                       const float *__restrict az, const float *__restrict aw,
                                       float *__restrict ox, float *__restrict oy,
                                       float *__restrict oz, float *__restrict ow, size_t n)
{
    size_t i = 0;
#if FV_I_INTRINSICS
    if (__builtin_cpu_supports("avx2") != 0) {
        size_t bulk = n & ~(size_t)7;
        if (bulk != 0) {
            fv_i_v4_normalize_soa_avx2(ax, ay, az, aw, ox, oy, oz, ow, bulk);
        }
        i = bulk;
    }
#endif
    for (; i < n; ++i) {
        fv_i_s_normalize4(ax[i], ay[i], az[i], aw[i], &ox[i], &oy[i], &oz[i], &ow[i]);
    }
}


#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

/* ======== 公开名：默认形态下的 static inline 定义 ========
 * 函数体只有一行转发，不含任何浮点表达式，所以全部留在 pragma 之外（作用域越小越好）；
 * 转发没有语义，数值全部由 fv_i_* 决定，因此两种形态（这里 / .cpp 里的 extern 定义）
 * 是同一份实现的两种接线方式。 */
#if !defined(FV_BUILD)

static inline void fv_vec2_add_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                  fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4(a, b, o, n);
}


static inline void fv_vec2_sub_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                  fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4(a, b, o, n);
}


static inline void fv_vec2_scale_x4(const fv_vec2 *__restrict a, fv_real s,
                                    fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4(a, s, o, n);
}


static inline void fv_vec2_lerp_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_real t, fv_vec2 *__restrict o, size_t n)
{
    fv_i_lerp_x4(a, b, t, o, n);
}


static inline void fv_vec2_madd_x4(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                   fv_real s, fv_vec2 *__restrict o, size_t n)
{
    fv_i_madd_x4(a, b, s, o, n);
}


static inline void fv_vec2_add_x4_mask(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                       fv_vec2 *__restrict o, const uint8_t *__restrict mask,
                                       size_t n)
{
    fv_i_add_x4_mask(a, b, o, mask, n);
}


static inline void fv_vec2_add_x4_idx(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                      fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                      size_t m)
{
    fv_i_add_x4_idx(a, b, o, idx, m);
}


static inline void fv_vec2_scale_x4_idx(const fv_vec2 *__restrict a, fv_real s,
                                        fv_vec2 *__restrict o, const uint32_t *__restrict idx,
                                        size_t m)
{
    fv_i_scale_x4_idx(a, s, o, idx, m);
}


static inline void fv_vec2_add_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict bx, const float *__restrict by,
                                   float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa(ax, ay, bx, by, ox, oy, n);
}


static inline void fv_vec2_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                     fv_real s, float *__restrict ox, float *__restrict oy,
                                     size_t n)
{
    fv_i_scale_soa(ax, ay, s, ox, oy, n);
}


static inline float fv_vec2_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict bx, const float *__restrict by,
                                    size_t n)
{
    return fv_i_dot_soa(ax, ay, bx, by, n);
}


static inline void fv_vec2_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   float *__restrict o, size_t n)
{
    fv_i_len_soa(ax, ay, o, n);
}


static inline void fv_vec2_normalize_soa(const float *__restrict ax, const float *__restrict ay,
                                         float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_normalize_soa(ax, ay, ox, oy, n);
}


static inline void fv_vec2_add_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                         fv_vec2 *__restrict o, size_t n)
{
    fv_i_add_x4_stream(a, b, o, n);
}


static inline void fv_vec2_sub_x4_stream(const fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                         fv_vec2 *__restrict o, size_t n)
{
    fv_i_sub_x4_stream(a, b, o, n);
}


static inline void fv_vec2_scale_x4_stream(const fv_vec2 *__restrict a, fv_real s,
                                           fv_vec2 *__restrict o, size_t n)
{
    fv_i_scale_x4_stream(a, s, o, n);
}


static inline void fv_vec2_add_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                          const float *__restrict bx, const float *__restrict by,
                                          float *__restrict ox, float *__restrict oy, size_t n)
{
    fv_i_add_soa_stream(ax, ay, bx, by, ox, oy, n);
}


static inline void fv_vec2_scale_soa_stream(const float *__restrict ax, const float *__restrict ay,
                                            fv_real s, float *__restrict ox, float *__restrict oy,
                                            size_t n)
{
    fv_i_scale_soa_stream(ax, ay, s, ox, oy, n);
}


static inline void fv_vec2_normalize_soa_stream(const float *__restrict ax,
                                                const float *__restrict ay,
                                                float *__restrict ox, float *__restrict oy,
                                                size_t n)
{
    fv_i_normalize_soa_stream(ax, ay, ox, oy, n);
}


static inline void fv_vec2_add_x4_mask_inplace(fv_vec2 *__restrict a, const fv_vec2 *__restrict b,
                                               const uint8_t *__restrict mask, size_t n)
{
    fv_i_add_x4_mask_inplace(a, b, mask, n);
}

/* 公开名：默认形态下的 static inline 转发（与上面 vec2 的 20 个转发同构）。
 * 函数体只有一行转发，没有任何浮点表达式，所以全部留在 pragma 之外。 */
static inline void fv_vec3_add_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, const float *__restrict bx,
                                   const float *__restrict by, const float *__restrict bz,
                                   float *__restrict ox, float *__restrict oy,
                                   float *__restrict oz, size_t n)
{
    fv_i_v3_add_soa(ax, ay, az, bx, by, bz, ox, oy, oz, n);
}


static inline void fv_vec3_scale_soa(const float *__restrict ax, const float *__restrict ay,
                                     const float *__restrict az, fv_real s,
                                     float *__restrict ox, float *__restrict oy,
                                     float *__restrict oz, size_t n)
{
    fv_i_v3_scale_soa(ax, ay, az, s, ox, oy, oz, n);
}


static inline void fv_vec4_add_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                  fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_add_x2(a, b, o, n);
}


static inline void fv_vec4_sub_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                  fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_sub_x2(a, b, o, n);
}


static inline void fv_vec4_scale_x2(const fv_vec4 *__restrict a, fv_real s,
                                    fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_scale_x2(a, s, o, n);
}


static inline void fv_vec4_lerp_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                   fv_real t, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_lerp_x2(a, b, t, o, n);
}


static inline void fv_vec4_madd_x2(const fv_vec4 *__restrict a, const fv_vec4 *__restrict b,
                                   fv_real s, fv_vec4 *__restrict o, size_t n)
{
    fv_i_v4_madd_x2(a, b, s, o, n);
}


static inline float fv_vec3_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict az, const float *__restrict bx,
                                    const float *__restrict by, const float *__restrict bz,
                                    size_t n)
{
    return fv_i_v3_dot_soa(ax, ay, az, bx, by, bz, n);
}


static inline void fv_vec3_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, float *__restrict o, size_t n)
{
    fv_i_v3_len_soa(ax, ay, az, o, n);
}


static inline void fv_vec3_normalize_soa(const float *__restrict ax,
                                         const float *__restrict ay,
                                         const float *__restrict az, float *__restrict ox,
                                         float *__restrict oy, float *__restrict oz, size_t n)
{
    fv_i_v3_normalize_soa(ax, ay, az, ox, oy, oz, n);
}


static inline float fv_vec4_dot_soa(const float *__restrict ax, const float *__restrict ay,
                                    const float *__restrict az, const float *__restrict aw,
                                    const float *__restrict bx, const float *__restrict by,
                                    const float *__restrict bz, const float *__restrict bw,
                                    size_t n)
{
    return fv_i_v4_dot_soa(ax, ay, az, aw, bx, by, bz, bw, n);
}


static inline void fv_vec4_len_soa(const float *__restrict ax, const float *__restrict ay,
                                   const float *__restrict az, const float *__restrict aw,
                                   float *__restrict o, size_t n)
{
    fv_i_v4_len_soa(ax, ay, az, aw, o, n);
}


static inline void fv_vec4_normalize_soa(const float *__restrict ax,
                                         const float *__restrict ay,
                                         const float *__restrict az,
                                         const float *__restrict aw, float *__restrict ox,
                                         float *__restrict oy, float *__restrict oz,
                                         float *__restrict ow, size_t n)
{
    fv_i_v4_normalize_soa(ax, ay, az, aw, ox, oy, oz, ow, n);
}

#endif /* !FV_BUILD：默认形态的公开名 */

#endif /* !FV_USE_EXTERN：实现小节 */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FLASHVEC_H */