/* bench/ab/ab_ops.h -- the op registry (the only place this fixture adds anything)
 *
 * Adding an op takes exactly two edits:
 *   1) one line in AB_LIB_OPS(X) below (or AB_MEM_OPS: a ceiling reference that bypasses
 *      the variant TU)
 *   2) one thunk function above kRows in ab_variant.cpp (the call code)
 * Everything else derives from this table; driver logic never changes:
 *   - the AB_OP_* enum                (expanded by ab_abi.h)
 *   - the fv_* symbol rename table    (Makefile generates it from the **selected header**,
 *                                      see mkrename.sh)
 *   - checksum gate / small-N support (AB_OPF_SMALL flag)
 *   - throughput matrix / TSV / report (NAME and the byte counts)
 *   - RFO real-traffic accounting     (the RD/WR/RFO columns)
 *
 * AB_LIB_OPS row format (10 columns):
 *   X(ID, NAME, KIND, POL, REF, FLAGS, RD, WR, RFO, NOTE)
 *     ID     C identifier; enum name = AB_OP_<ID> (both variant TUs use the same number)
 *     NAME   name used in the TSV/report (string literal; pipeline ops have no lib symbol)
 *     KIND   output form: VEC2 / SOA2 / SOA4 / FLOAT / NONE
 *     POL    check policy: EXACT / REL / NOCHECK
 *     REF    independent scalar reference implementation in the driver (function name)
 *     FLAGS  AB_OPF_* (see below)
 *     RD/WR  **nominal** read/write bytes per element (hand-checked per function, not
 *            "array size x constant")
 *     RFO    write-allocate behaviour of the plain store: AB_RFO_MISS / HIT / MIXED /
 *            LIBC / NT
 *     NOTE   human-readable note (no commas -- keeps column parsing easy for scripts)
 *
 * AB_MEM_OPS row format (8 columns):
 *   X(NAME, KIND, LOCAL, FLAGS, RD, WR, RFO, NOTE)
 *     NAME   same as above (string literal)
 *     LOCAL  local implementation inside the driver (ceiling reference: memcpy/memset/
 *            pure read/triad)
 *
 * RFO classes (the key to real memory traffic):
 *   AB_RFO_MISS   plain-store write-allocate miss: write 8 B once, then RFO pulls that
 *                 cache line in once more -> real traffic = 2 x nominal write bytes
 *   AB_RFO_HIT    the write target **was already read** before the write (in-place
 *                 update) so the line is in cache -> no extra RFO
 *   AB_RFO_MIXED  intermediate arrays between pipeline passes may be absorbed by L2/L3;
 *                 2x is an **upper bound** and is flagged as such
 *   AB_RFO_LIBC   libc picks NT store by size (memcpy/memset) so nothing is filled in;
 *                 flagged as unsure
 *   AB_RFO_NT     NT store (_mm256_stream_ps): no RFO -> real traffic = nominal write bytes
 *                 **whether NT was really used is reported at runtime by the variant TU**
 *                 (see *_opflags in ab_abi.h)
 */
#ifndef AB_OPS_H
#define AB_OPS_H

/* ---- driver flags ---- */
#define AB_OPF_PIPE 0x1u    /* pipeline composite op (--no-pipe filters it; needs init_pipe) */
#define AB_OPF_SMALL 0x2u   /* enters the small-N regime: ab_s_loop implements it */
#define AB_OPF_INPLACE 0x4u /* in-place: writes its own first input array */
/* This entry does not exist in the baseline (vec3/vec4 are new API): A/B is meaningless;
 * only --intra and the "candidate vs independent scalar reference" gate apply. It must
 * never be reported as an "improvement over HEAD". */
#define AB_OPF_NO_BASELINE 0x8u
/* ---- ceiling-reference roles (the driver picks by role, not by index) ---- */
#define AB_MF_COPY 0x10u
#define AB_MF_SET 0x20u
#define AB_MF_READ 0x40u
#define AB_MF_TRIAD 0x80u

/* ---- RFO classes ---- */
#define AB_RFO_MISS 0
#define AB_RFO_HIT 1
#define AB_RFO_MIXED 2
#define AB_RFO_LIBC 3
#define AB_RFO_NT 4

/* ---- Matrix transforms (16 B AoS elements, one matrix operand) ----
 * fv_mat4_transform_points / _dirs move 16 B in and 16 B out per element; the kernel stores
 * with a plain store (_mm256_storeu_ps in the AVX2 helper, the scalar helper for the tail), so
 * the write-allocate miss adds another 16 B: 16 + 16 + 16 RFO = 48 B/e real. Object-code
 * check: no _mm256_stream_ps anywhere in the two kernels, hence AB_RFO_MISS, not AB_RFO_NT.
 * n is the element count (each element is one fv_vec4 = 4 floats), and the driver derives
 * elements = n * 8 / 16 like every other 16 B op, so the byte volume stays identical.
 *
 * AB_MAT4_IN_BASELINE: the entries landed after the current git HEAD snapshot, so the driver
 * is told at compile time whether the baseline has them (the Makefile greps the snapshot
 * header). While it does not, the rows carry AB_OPF_NO_BASELINE and INTRA mode is the only
 * way to measure them -- the fixture never fabricates an A/B ratio for an op the baseline
 * cannot run. Once the API is in HEAD the macro flips to 1 and the same rows become normal
 * A/B entries. */
#ifndef AB_MAT4_IN_BASELINE
#define AB_MAT4_IN_BASELINE 0
#endif
#define AB_OPF_MAT4 (AB_MAT4_IN_BASELINE ? 0u : AB_OPF_NO_BASELINE)

/* =====================================================================
 * Library batch entry points (each variant TU implements its own; the driver only ever
 * calls ab_*_loop). Order = historical order: the TSV / baseline / A/A floor table all
 * align by name, so do not insert in the middle -- new ops go at the end.
 * ===================================================================== */
#define AB_LIB_OPS(X)                                                                               \
    X(ADD_X4, "fv_vec2_add_x4", VEC2, 8, EXACT, ref_add_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,       \
      "read 16 B/e from a/b; the 8 B/e write to o RFOs 8 B/e -> 32 B/e real")                       \
    X(SUB_X4, "fv_vec2_sub_x4", VEC2, 8, EXACT, ref_sub_x4, 0u, 16, 8, AB_RFO_MISS,                 \
      "same as add: 16 read + 8 write + 8 RFO = 32 B/e")                                            \
    X(SCALE_X4, "fv_vec2_scale_x4", VEC2, 8, EXACT, ref_scale_x4, AB_OPF_SMALL, 8, 8,               \
      AB_RFO_MISS, "only a and o: 8 read + 8 write + 8 RFO = 24 B/e")                               \
    X(LERP_X4, "fv_vec2_lerp_x4", VEC2, 8, EXACT, ref_lerp_x4, AB_OPF_SMALL, 16, 8, AB_RFO_MISS,    \
      "16 read + 8 write + 8 RFO = 32 B/e")                                                         \
    X(MADD_X4, "fv_vec2_madd_x4", VEC2, 8, EXACT, ref_madd_x4, 0u, 16, 8, AB_RFO_MISS,              \
      "16 read + 8 write + 8 RFO = 32 B/e")                                                         \
    X(ADD_MASK, "fv_vec2_add_x4_mask", VEC2, 8, EXACT, ref_add_mask, 0u, 17, 8, AB_RFO_MISS,        \
      "reads 1 more B for the mask: 17 read + 8 write + 8 RFO = 33 B/e")                            \
    X(ADD_IDX, "fv_vec2_add_x4_idx", VEC2, 8, EXACT, ref_add_idx, 0u, 20, 8, AB_RFO_MISS,           \
      "reads 4 B/e more for indices; scattered writes still RFO per line: 36 B/e real")             \
    X(SCALE_IDX, "fv_vec2_scale_x4_idx", VEC2, 8, EXACT, ref_scale_idx, 0u, 12, 8, AB_RFO_MISS,     \
      "no b: 12 read + 8 write + 8 RFO = 28 B/e")                                                   \
    X(ADD_SOA, "fv_vec2_add_soa", SOA2, 8, EXACT, ref_add_soa, 0u, 16, 8, AB_RFO_MISS,              \
      "ax/ay/bx/by are four component arrays: 16 read + 8 write + 8 RFO = 32 B/e")                  \
    X(SCALE_SOA, "fv_vec2_scale_soa", SOA2, 8, EXACT, ref_scale_soa, 0u, 8, 8, AB_RFO_MISS,         \
      "8 read + 8 write + 8 RFO = 24 B/e")                                                          \
    X(DOT_SOA, "fv_vec2_dot_soa", NONE, 8, REL, ref_dot_soa, AB_OPF_SMALL, 16, 0, AB_RFO_MISS,      \
      "pure read: 16 B/e; read-bandwidth bound and 0 write traffic")                                \
    X(LEN_SOA, "fv_vec2_len_soa", FLOAT, 8, EXACT, ref_len_soa, AB_OPF_SMALL, 8, 4, AB_RFO_MISS,    \
      "8 read + 4 write + 4 RFO = 16 B/e")                                                          \
    X(NORM_SOA, "fv_vec2_normalize_soa", SOA2, 8, EXACT, ref_norm_soa, AB_OPF_SMALL, 8, 8,          \
      AB_RFO_MISS, "8 read + 8 write + 8 RFO = 24 B/e")                                             \
    X(COLLISION, "pipe_collision", SOA4, 8, EXACT, ref_collision, AB_OPF_PIPE, 24, 16,              \
      AB_RFO_MIXED, "two passes 24 read + 16 write; the ox/oy intermediates may stay in L3 -> 56 B/e upper bound") \
    X(CHAIN3, "pipe_scale_add_normalize", SOA2, 8, EXACT, ref_chain3, AB_OPF_PIPE, 32, 24,          \
      AB_RFO_MIXED, "three passes 32 read + 24 write; px/qx may be absorbed -> 80 B/e upper bound") \
    X(MAT4_POINTS, "fv_mat4_transform_points", VEC4M, 16, EXACT, ref_mat4_points,              \
      AB_OPF_MAT4, 16, 16, AB_RFO_MISS,                                                          \
      "in 16 read + out 16 write + 16 RFO = 48 B/e; in.w is the constant 1 (translation on)")    \
    X(MAT4_DIRS, "fv_mat4_transform_dirs", VEC4M, 16, EXACT, ref_mat4_dirs,                      \
      AB_OPF_MAT4, 16, 16, AB_RFO_MISS,                                                          \
      "in 16 read + out 16 write + 16 RFO = 48 B/e; in.w is the constant 0 (translation off)")   \
    /* new entry points are appended above this line, in the format shown */

/* =====================================================================
 * Ceiling references: bypass the variant TU, so A/B runs the very same code (the ratio is
 * a pure noise probe).
 * ===================================================================== */
#define AB_MEM_OPS(X)                                                                               \
    X("ref_memcpy", VEC2, mem_memcpy, AB_MF_COPY, 8, 8, AB_RFO_LIBC,                                \
      "glibc picks NT store by size: 8-16 B/e real write traffic (flagged unsure)")                 \
    X("ref_memset", VEC2, mem_set, AB_MF_SET, 0, 8, AB_RFO_LIBC,                                    \
      "same: libc decides on NT store - 8 B/e nominal and 8-16 B/e real")                           \
    X("ref_pure_read", NONE, mem_read, AB_MF_READ, 8, 0, AB_RFO_MISS,                               \
      "pure read 8 B/e and 0 write; load-port bound when small")                                    \
    X("ref_triad", VEC2, mem_triad, AB_MF_TRIAD, 16, 8, AB_RFO_MISS,                                \
      "same shape as the AoS batch entry: 16 read + 8 write + 8 RFO = 32")

/* =====================================================================
 * Experimental area: newly landed entry points (_stream / _mask_inplace / ...)
 * Off by default (AB_OPS_EXTRA=0) so the default gate is still the historical 182 items.
 * Enable with: make ... AB_OPS_EXTRA=1 (verify turns it on automatically).
 * The _stream family is guarded by AB_HAVE_STREAM: when the baseline (= git HEAD snapshot)
 * has none of these entry points yet, the thunk degrades to a plain store -- and then this
 * very row becomes the "plain store vs NT store" control.
 * ===================================================================== */
#ifndef AB_OPS_EXTRA
#define AB_OPS_EXTRA 0
#endif

#if AB_OPS_EXTRA
#define AB_LIB_OPS_EXTRA(X)                                                                         \
    X(ADD_X4_STREAM, "fv_vec2_add_x4_stream", VEC2, 8, EXACT, ref_add_x4, 0u, 16, 8, AB_RFO_NT,     \
      "NT store: 16 read + 8 write = 24 B/e (the plain 32 B/e version RFOs 8 B)")                   \
    X(SUB_X4_STREAM, "fv_vec2_sub_x4_stream", VEC2, 8, EXACT, ref_sub_x4, 0u, 16, 8, AB_RFO_NT,     \
      "NT store: 16 read + 8 write = 24 B/e")                                                       \
    X(SCALE_X4_STREAM, "fv_vec2_scale_x4_stream", VEC2, 8, EXACT, ref_scale_x4, 0u, 8, 8,           \
      AB_RFO_NT, "NT store: 8 read + 8 write = 16 B/e (the plain version is 24 B/e)")               \
    X(ADD_SOA_STREAM, "fv_vec2_add_soa_stream", SOA2, 8, EXACT, ref_add_soa, 0u, 16, 8, AB_RFO_NT,  \
      "NT store: 16 read + 8 write = 24 B/e")                                                       \
    X(SCALE_SOA_STREAM, "fv_vec2_scale_soa_stream", SOA2, 8, EXACT, ref_scale_soa, 0u, 8, 8,        \
      AB_RFO_NT, "NT store: 8 read + 8 write = 16 B/e")                                             \
    X(NORM_SOA_STREAM, "fv_vec2_normalize_soa_stream", SOA2, 8, EXACT, ref_norm_soa, 0u, 8, 8,      \
      AB_RFO_NT, "NT store: 8 read + 8 write = 16 B/e")                                             \
    X(ADD_MASK_INPLACE, "fv_vec2_add_x4_mask_inplace", VEC2, 8, EXACT, ref_add_mask_inplace,        \
      AB_OPF_INPLACE, 17, 8, AB_RFO_HIT,                                                            \
      "in-place: a is read then written so the line is cached -> no RFO; 25 B/e real = nominal")    \
    /* ---- 3D SoA (12 B elements) ---- */                                                          \
    X(V3_ADD_SOA, "fv_vec3_add_soa", SOA3, 12, EXACT, ref_v3_add, AB_OPF_NO_BASELINE, 24, 12,       \
      AB_RFO_MISS, "three component-array pairs: 24 read + 12 write + 12 RFO = 48 B/e")             \
    X(V3_SCALE_SOA, "fv_vec3_scale_soa", SOA3, 12, EXACT, ref_v3_scale, AB_OPF_NO_BASELINE, 12, 12, \
      AB_RFO_MISS, "12 read + 12 write + 12 RFO = 36 B/e")                                          \
    X(V3_DOT_SOA, "fv_vec3_dot_soa", NONE, 12, REL, ref_v3_dot, AB_OPF_NO_BASELINE, 24, 0,          \
      AB_RFO_MISS, "pure read 24 B/e; the reduction order is not guaranteed bit-identical -> REL")  \
    X(V3_LEN_SOA, "fv_vec3_len_soa", FLOAT, 12, EXACT, ref_v3_len, AB_OPF_NO_BASELINE, 12, 4,       \
      AB_RFO_MISS, "12 read + 4 write + 4 RFO = 20 B/e")                                            \
    X(V3_NORM_SOA, "fv_vec3_normalize_soa", SOA3, 12, EXACT, ref_v3_norm, AB_OPF_NO_BASELINE, 12,   \
      12, AB_RFO_MISS, "12 read + 12 write + 12 RFO = 36 B/e; zero guard as in vec2")               \
    /* ---- 4D AoS _x2 (16 B elements) ---- */                                                      \
    X(V4_ADD_X2, "fv_vec4_add_x2", VEC4, 16, EXACT, ref_v4_add_x2, AB_OPF_NO_BASELINE, 32, 16,      \
      AB_RFO_MISS, "16 read each from a/b + 16 write to o + 16 RFO = 64 B/e")                       \
    X(V4_SUB_X2, "fv_vec4_sub_x2", VEC4, 16, EXACT, ref_v4_sub_x2, AB_OPF_NO_BASELINE, 32, 16,      \
      AB_RFO_MISS, "same as add: 64 B/e")                                                           \
    X(V4_SCALE_X2, "fv_vec4_scale_x2", VEC4, 16, EXACT, ref_v4_scale_x2, AB_OPF_NO_BASELINE, 16,    \
      16, AB_RFO_MISS, "16 read + 16 write + 16 RFO = 48 B/e")                                      \
    X(V4_LERP_X2, "fv_vec4_lerp_x2", VEC4, 16, EXACT, ref_v4_lerp_x2, AB_OPF_NO_BASELINE, 32, 16,   \
      AB_RFO_MISS, "same as add: 64 B/e")                                                           \
    X(V4_MADD_X2, "fv_vec4_madd_x2", VEC4, 16, EXACT, ref_v4_madd_x2, AB_OPF_NO_BASELINE, 32, 16,   \
      AB_RFO_MISS, "same as add: 64 B/e")                                                           \
    /* ---- 4D SoA (16 B elements) ---- */                                                          \
    X(V4_DOT_SOA, "fv_vec4_dot_soa", NONE, 16, REL, ref_v4_dot, AB_OPF_NO_BASELINE, 32, 0,          \
      AB_RFO_MISS, "pure read 32 B/e (8 component arrays); REL")                                    \
    X(V4_LEN_SOA, "fv_vec4_len_soa", FLOAT, 16, EXACT, ref_v4_len, AB_OPF_NO_BASELINE, 32, 4,       \
      AB_RFO_MISS, "32 read + 4 write + 4 RFO = 40 B/e")                                            \
    X(V4_NORM_SOA, "fv_vec4_normalize_soa", SOA4V, 16, EXACT, ref_v4_norm, AB_OPF_NO_BASELINE, 32,  \
      16, AB_RFO_MISS, "32 read + 16 write + 16 RFO = 64 B/e")                                      \
    /* new entry points go here (same format as AB_LIB_OPS) */


/* ---- INTRA pairs: new entry (ID) -> same-byte-volume reference entry (ID) inside the
 * candidate. Why this is not "relative to baseline": these entries do not exist in git HEAD
 * at all, so comparing them against "nothing" measures "something from nothing", not
 * "how much faster". INTRA mode therefore compares inside the candidate only, normalized by
 * bytes (in one tier both ops move the same number of bytes: elements = n * 8 / element
 * width). The report must say so: intra-candidate comparison, **not** an improvement over
 * HEAD. */
#define AB_INTRA_PAIRS(X)                                                                           \
    X(V4_ADD_X2, ADD_X4)                                                                            \
    X(V4_SUB_X2, SUB_X4)                                                                            \
    X(V4_SCALE_X2, SCALE_X4)                                                                        \
    X(V4_LERP_X2, LERP_X4)                                                                          \
    X(V4_MADD_X2, MADD_X4)                                                                          \
    X(V4_DOT_SOA, DOT_SOA)                                                                          \
    X(V4_LEN_SOA, LEN_SOA)                                                                          \
    X(V4_NORM_SOA, NORM_SOA)                                                                        \
    X(V3_ADD_SOA, ADD_SOA)                                                                          \
    X(V3_SCALE_SOA, SCALE_SOA)                                                                      \
    X(V3_DOT_SOA, DOT_SOA)                                                                          \
    X(V3_LEN_SOA, LEN_SOA)                                                                          \
    X(V3_NORM_SOA, NORM_SOA)
#else
#define AB_INTRA_PAIRS(X)
#define AB_LIB_OPS_EXTRA(X)
#endif

#endif /* AB_OPS_H */
