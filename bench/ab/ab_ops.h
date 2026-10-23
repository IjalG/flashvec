/* bench/ab/ab_ops.h —— op 注册表（本夹具唯一需要"加东西"的地方）
 *
 * 加一个 op 只有两处：
 *   1) 在下面 AB_LIB_OPS(X) 里加一行（或 AB_MEM_OPS：不经过变体 TU 的上限参考）
 *   2) 在 ab_variant.cpp 的 kRows 前面加一个 thunk 函数（那段调用代码）
 * 其余全部由此表派生，驱动逻辑一行都不用改：
 *   - AB_OP_* 枚举                （ab_abi.h 展开）
 *   - fv_* 符号改名表             （Makefile 从**所选头文件**自动生成，见 mkrename.sh）
 *   - 校验和门禁 / 小 N 支持判定  （AB_OPF_SMALL 标志）
 *   - 吞吐矩阵 / TSV / 报告列     （NAME 与字节数）
 *   - RFO 真实流量记账            （RD/WR/RFO 三列）
 *
 * AB_LIB_OPS 行格式（10 列）：
 *   X(ID, NAME, KIND, POL, REF, FLAGS, RD, WR, RFO, NOTE)
 *     ID     C 标识符；枚举名 = AB_OP_<ID>（两侧变体 TU 用同一个编号）
 *     NAME   TSV/报告里的名字（字符串字面量；管线算子没有对应的库符号）
 *     KIND   输出形态：VEC2 / SOA2 / SOA4 / FLOAT / NONE
 *     POL    校验策略：EXACT / REL / NOCHECK
 *     REF    驱动内独立标量参考实现（函数名）
 *     FLAGS  AB_OPF_*（见下）
 *     RD/WR  **名义**每元素读/写字节（逐函数人工核对，不是"数组大小×常数"）
 *     RFO    普通 store 的写分配行为：AB_RFO_MISS / HIT / MIXED / LIBC / NT
 *     NOTE   人读备注（不能含逗号——便于外部脚本按列解析）
 *
 * AB_MEM_OPS 行格式（8 列）：
 *   X(NAME, KIND, LOCAL, FLAGS, RD, WR, RFO, NOTE)
 *     NAME   同上（字符串字面量）
 *     LOCAL  驱动内的本地实现（上限参考：memcpy/memset/纯读/triad）
 *
 * RFO 分类（真实内存流量的关键）：
 *   AB_RFO_MISS   普通 store 写分配 miss：写一次 8 B，RFO 再把该 cache line 读进来一次
 *                 → 该数组真实流量 = 2 × 名义写字节
 *   AB_RFO_HIT    写目标在写之前**已经被读过**（就地更新），行已在 cache → 无额外 RFO
 *   AB_RFO_MIXED  管线里 pass 之间的中间数组可能被 L2/L3 吸收；2× 是**上界**，标注
 *   AB_RFO_LIBC   libc 自己按规模选 NT store（memcpy/memset），不硬填，标注不确定
 *   AB_RFO_NT     NT store（_mm256_stream_ps）：不 RFO → 真实流量 = 名义写字节
 *                 **是否真的走了 NT 由变体 TU 在运行期报告**（见 ab_abi.h 的 *_opflags）
 */
#ifndef AB_OPS_H
#define AB_OPS_H

/* ---- 驱动标志 ---- */
#define AB_OPF_PIPE 0x1u    /* 管线组合算子（--no-pipe 过滤、需要 init_pipe） */
#define AB_OPF_SMALL 0x2u   /* 进小 N 区制：标量锚点 ab_s_loop 已实现该算子 */
#define AB_OPF_INPLACE 0x4u /* 就地更新：写的是自己的第一个输入数组 */
/* 基线里不存在这个入口（vec3/vec4 是新增 API）：A/B 模式无意义，只能走 --intra 与
 * "candidate vs 独立标量参考"的门禁。绝不能把它报成"相对 HEAD 的改进"。 */
#define AB_OPF_NO_BASELINE 0x8u
/* ---- 上限参考的角色（驱动按角色取用，不按下标） ---- */
#define AB_MF_COPY 0x10u
#define AB_MF_SET 0x20u
#define AB_MF_READ 0x40u
#define AB_MF_TRIAD 0x80u

/* ---- RFO 分类 ---- */
#define AB_RFO_MISS 0
#define AB_RFO_HIT 1
#define AB_RFO_MIXED 2
#define AB_RFO_LIBC 3
#define AB_RFO_NT 4

/* =====================================================================
 * 库内批量入口（两侧变体 TU 各自实现，驱动只认 ab_*_loop）
 * 顺序 = 历史顺序：TSV / 基线 / A_A 地板表按名字对齐，别插队，新 op 加到末尾
 * ===================================================================== */
#define AB_LIB_OPS(X)                                                                          \
    X(ADD_X4, "fv_vec2_add_x4", VEC2, 8, EXACT, ref_add_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,     \
      "a/b 读 16；o 写 8 触发 RFO 读 8 → 真实 32 B/e")                                         \
    X(SUB_X4, "fv_vec2_sub_x4", VEC2, 8, EXACT, ref_sub_x4, 0u, 16, 8, AB_RFO_MISS,               \
      "同 add：16 读 + 8 写 + 8 RFO = 32 B/e")                                                 \
    X(SCALE_X4, "fv_vec2_scale_x4", VEC2, 8, EXACT, ref_scale_x4, AB_OPF_SMALL, 8, 8,             \
      AB_RFO_MISS, "只有 a 与 o：8 读 + 8 写 + 8 RFO = 24 B/e")                                \
    X(LERP_X4, "fv_vec2_lerp_x4", VEC2, 8, EXACT, ref_lerp_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,  \
      "16 读 + 8 写 + 8 RFO = 32 B/e")                                                         \
    X(MADD_X4, "fv_vec2_madd_x4", VEC2, 8, EXACT, ref_madd_x4, 0u, 16, 8, AB_RFO_MISS,            \
      "16 读 + 8 写 + 8 RFO = 32 B/e")                                                         \
    X(ADD_MASK, "fv_vec2_add_x4_mask", VEC2, 8, EXACT, ref_add_mask, 0u, 17, 8, AB_RFO_MISS,      \
      "多读 1 B 掩码：17 读 + 8 写 + 8 RFO = 33 B/e")                                         \
    X(ADD_IDX, "fv_vec2_add_x4_idx", VEC2, 8, EXACT, ref_add_idx, 0u, 20, 8, AB_RFO_MISS,         \
      "多读 4 B 下标；散列写仍逐行 RFO：20 读 + 8 写 + 8 RFO = 36 B/e")                         \
    X(SCALE_IDX, "fv_vec2_scale_x4_idx", VEC2, 8, EXACT, ref_scale_idx, 0u, 12, 8, AB_RFO_MISS,   \
      "无 b：12 读 + 8 写 + 8 RFO = 28 B/e")                                                   \
    X(ADD_SOA, "fv_vec2_add_soa", SOA2, 8, EXACT, ref_add_soa, 0u, 16, 8, AB_RFO_MISS,            \
      "ax/ay/bx/by 四条分量数组：16 读 + 8 写 + 8 RFO = 32 B/e")                                \
    X(SCALE_SOA, "fv_vec2_scale_soa", SOA2, 8, EXACT, ref_scale_soa, 0u, 8, 8, AB_RFO_MISS,       \
      "8 读 + 8 写 + 8 RFO = 24 B/e")                                                          \
    X(DOT_SOA, "fv_vec2_dot_soa", NONE, 8, REL, ref_dot_soa, AB_OPF_SMALL, 16, 0, AB_RFO_MISS,    \
      "纯读：16 B/e；与读带宽有关 写流量为 0")                                                 \
    X(LEN_SOA, "fv_vec2_len_soa", FLOAT, 8, EXACT, ref_len_soa, AB_OPF_SMALL, 8, 4, AB_RFO_MISS,  \
      "8 读 + 4 写 + 4 RFO = 16 B/e")                                                          \
    X(NORM_SOA, "fv_vec2_normalize_soa", SOA2, 8, EXACT, ref_norm_soa, AB_OPF_SMALL, 8, 8,        \
      AB_RFO_MISS, "8 读 + 8 写 + 8 RFO = 24 B/e")                                             \
    X(COLLISION, "pipe_collision", SOA4, 8, EXACT, ref_collision, AB_OPF_PIPE, 24, 16,             \
      AB_RFO_MIXED, "两 pass 24 读 + 16 写；中间 ox/oy 可能留在 L3 → 56 B/e 是上界")            \
    X(CHAIN3, "pipe_scale_add_normalize", SOA2, 8, EXACT, ref_chain3, AB_OPF_PIPE, 32, 24,        \
      AB_RFO_MIXED, "三 pass 32 读 + 24 写；px/qx 中间数组可能被吸收 → 80 B/e 是上界")          \
    /* 新入口按上面的格式追加到这一行之上 */

/* =====================================================================
 * 上限参考：不经过变体 TU，A/B 跑同一段代码（比值即纯噪声探针）
 * ===================================================================== */
#define AB_MEM_OPS(X)                                                                          \
    X("ref_memcpy", VEC2, mem_memcpy, AB_MF_COPY, 8, 8, AB_RFO_LIBC,                            \
      "glibc 在大规模上自己选 NT store：真实写流量介于 8 与 16 B/e 之间 标注不确定")            \
    X("ref_memset", VEC2, mem_set, AB_MF_SET, 0, 8, AB_RFO_LIBC,                               \
      "同上：NT store 与否由 libc 决定 名义 8 B/e 真实 8-16 B/e")                               \
    X("ref_pure_read", NONE, mem_read, AB_MF_READ, 8, 0, AB_RFO_MISS,                           \
      "纯读 8 B/e 写 0；小规模上受 load-port 限制而非带宽")                                     \
    X("ref_triad", VEC2, mem_triad, AB_MF_TRIAD, 16, 8, AB_RFO_MISS,                            \
      "与 AoS 批量入口同形态：16 读 + 8 写 + 8 RFO = 32 B/e")

/* =====================================================================
 * 实验区：新落地的入口（_stream / _mask_inplace / …）
 * 默认关闭（AB_OPS_EXTRA=0），这样默认门禁仍是历史上那份 182 项。
 * 打开方式：make ... AB_OPS_EXTRA=1（verify 会自动开）。
 * _stream 系列用 AB_HAVE_STREAM 兜底：baseline = git HEAD 快照里还没有这些入口时，
 * thunk 退化成普通 store —— 于是这一行恰好变成"普通 store vs NT store"的对照。
 * ===================================================================== */
#ifndef AB_OPS_EXTRA
#define AB_OPS_EXTRA 0
#endif

#if AB_OPS_EXTRA
#define AB_LIB_OPS_EXTRA(X)                                                                    \
    X(ADD_X4_STREAM, "fv_vec2_add_x4_stream", VEC2, 8, EXACT, ref_add_x4, 0u, 16, 8, AB_RFO_NT,   \
      "NT store：16 读 + 8 写 = 24 B/e（普通版 32 B/e 的那 8 B RFO 没了）")                     \
    X(SUB_X4_STREAM, "fv_vec2_sub_x4_stream", VEC2, 8, EXACT, ref_sub_x4, 0u, 16, 8, AB_RFO_NT,    \
      "NT store：16 读 + 8 写 = 24 B/e")                                                      \
    X(SCALE_X4_STREAM, "fv_vec2_scale_x4_stream", VEC2, 8, EXACT, ref_scale_x4, 0u, 8, 8,          \
      AB_RFO_NT, "NT store：8 读 + 8 写 = 16 B/e（普通版 24 B/e）")                             \
    X(ADD_SOA_STREAM, "fv_vec2_add_soa_stream", SOA2, 8, EXACT, ref_add_soa, 0u, 16, 8, AB_RFO_NT, \
      "NT store：16 读 + 8 写 = 24 B/e")                                                      \
    X(SCALE_SOA_STREAM, "fv_vec2_scale_soa_stream", SOA2, 8, EXACT, ref_scale_soa, 0u, 8, 8,       \
      AB_RFO_NT, "NT store：8 读 + 8 写 = 16 B/e")                                             \
    X(NORM_SOA_STREAM, "fv_vec2_normalize_soa_stream", SOA2, 8, EXACT, ref_norm_soa, 0u, 8, 8,     \
      AB_RFO_NT, "NT store：8 读 + 8 写 = 16 B/e")                                             \
    X(ADD_MASK_INPLACE, "fv_vec2_add_x4_mask_inplace", VEC2, 8, EXACT, ref_add_mask_inplace,       \
      AB_OPF_INPLACE, 17, 8, AB_RFO_HIT,                                                        \
      "就地更新：a 先读后写 行已在 cache → 无 RFO 真实 25 B/e = 名义")                          \
    /* ---- 三维 SoA（元素 12 B）---- */                                                       \
    X(V3_ADD_SOA, "fv_vec3_add_soa", SOA3, 12, EXACT, ref_v3_add, AB_OPF_NO_BASELINE, 24, 12,       \
      AB_RFO_MISS, "三对分量数组：24 读 + 12 写 + 12 RFO = 48 B/e（元素 12 B）")                \
    X(V3_SCALE_SOA, "fv_vec3_scale_soa", SOA3, 12, EXACT, ref_v3_scale, AB_OPF_NO_BASELINE, 12, 12, \
      AB_RFO_MISS, "12 读 + 12 写 + 12 RFO = 36 B/e")                                          \
    X(V3_DOT_SOA, "fv_vec3_dot_soa", NONE, 12, REL, ref_v3_dot, AB_OPF_NO_BASELINE, 24, 0,          \
      AB_RFO_MISS, "纯读 24 B/e；归约顺序不保证逐位 → REL")                                     \
    X(V3_LEN_SOA, "fv_vec3_len_soa", FLOAT, 12, EXACT, ref_v3_len, AB_OPF_NO_BASELINE, 12, 4,       \
      AB_RFO_MISS, "12 读 + 4 写 + 4 RFO = 20 B/e")                                            \
    X(V3_NORM_SOA, "fv_vec3_normalize_soa", SOA3, 12, EXACT, ref_v3_norm, AB_OPF_NO_BASELINE, 12,   \
      12, AB_RFO_MISS, "12 读 + 12 写 + 12 RFO = 36 B/e；零向量保护与 vec2 同语义")             \
    /* ---- 四维 AoS _x2（元素 16 B）---- */                                                   \
    X(V4_ADD_X2, "fv_vec4_add_x2", VEC4, 16, EXACT, ref_v4_add_x2, AB_OPF_NO_BASELINE, 32, 16,      \
      AB_RFO_MISS, "a/b 各 16 读 + o 16 写 + 16 RFO = 64 B/e（元素 16 B）")                     \
    X(V4_SUB_X2, "fv_vec4_sub_x2", VEC4, 16, EXACT, ref_v4_sub_x2, AB_OPF_NO_BASELINE, 32, 16,      \
      AB_RFO_MISS, "同 add：64 B/e")                                                            \
    X(V4_SCALE_X2, "fv_vec4_scale_x2", VEC4, 16, EXACT, ref_v4_scale_x2, AB_OPF_NO_BASELINE, 16,    \
      16, AB_RFO_MISS, "16 读 + 16 写 + 16 RFO = 48 B/e")                                       \
    X(V4_LERP_X2, "fv_vec4_lerp_x2", VEC4, 16, EXACT, ref_v4_lerp_x2, AB_OPF_NO_BASELINE, 32, 16,   \
      AB_RFO_MISS, "同 add：64 B/e")                                                            \
    X(V4_MADD_X2, "fv_vec4_madd_x2", VEC4, 16, EXACT, ref_v4_madd_x2, AB_OPF_NO_BASELINE, 32, 16,   \
      AB_RFO_MISS, "同 add：64 B/e")                                                            \
    /* ---- 四维 SoA（元素 16 B）---- */                                                       \
    X(V4_DOT_SOA, "fv_vec4_dot_soa", NONE, 16, REL, ref_v4_dot, AB_OPF_NO_BASELINE, 32, 0,          \
      AB_RFO_MISS, "纯读 32 B/e（8 条分量数组）；REL")                                           \
    X(V4_LEN_SOA, "fv_vec4_len_soa", FLOAT, 16, EXACT, ref_v4_len, AB_OPF_NO_BASELINE, 32, 4,       \
      AB_RFO_MISS, "32 读 + 4 写 + 4 RFO = 40 B/e")                                            \
    X(V4_NORM_SOA, "fv_vec4_normalize_soa", SOA4V, 16, EXACT, ref_v4_norm, AB_OPF_NO_BASELINE, 32,  \
      16, AB_RFO_MISS, "32 读 + 16 写 + 16 RFO = 64 B/e")                                       \
    /* 新入口加在这里（格式同 AB_LIB_OPS） */

/* ---- INTRA 对照对：新入口(ID) -> candidate 内部同字节量的参照入口(ID) ----
 * 为什么不是"相对基线"：这些入口在 git HEAD 里根本不存在，拿它跟"不存在"比是
 * "从无到有"，不是"快了多少"。所以 INTRA 模式只在 candidate 内部比，并按字节归一
 * （同一档扫描里两个 op 处理相同字节量：元素数 = n * 8 / 元素宽度）。
 * 报告里必须写明：intra-candidate 对照，**不是**相对 HEAD 的改进。 */
#define AB_INTRA_PAIRS(X)                                                                      \
    X(V4_ADD_X2, ADD_X4)                                                                       \
    X(V4_SUB_X2, SUB_X4)                                                                       \
    X(V4_SCALE_X2, SCALE_X4)                                                                   \
    X(V4_LERP_X2, LERP_X4)                                                                     \
    X(V4_MADD_X2, MADD_X4)                                                                     \
    X(V4_DOT_SOA, DOT_SOA)                                                                     \
    X(V4_LEN_SOA, LEN_SOA)                                                                     \
    X(V4_NORM_SOA, NORM_SOA)                                                                   \
    X(V3_ADD_SOA, ADD_SOA)                                                                     \
    X(V3_SCALE_SOA, SCALE_SOA)                                                                 \
    X(V3_DOT_SOA, DOT_SOA)                                                                     \
    X(V3_LEN_SOA, LEN_SOA)                                                                     \
    X(V3_NORM_SOA, NORM_SOA)
#else
#define AB_INTRA_PAIRS(X)
#define AB_LIB_OPS_EXTRA(X)
#endif

#endif /* AB_OPS_H */
