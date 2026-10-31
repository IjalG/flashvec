# flashvec vs GLM: local re-run of the 115-item matrix

This is a local re-run of the external GLM report (`flashvec全量API对比GLM测试报告.md`,
2026-10-31), using the report's own `bench_all.cpp` and `summarize_all.py`. It is **not** the same
machine as the report and not the same GLM version, so the two matrices must not be compared cell by
cell; the useful question is whether the structural negatives that the report found still exist.

## Setup

| Item | Value |
|---|---|
| CPU | AMD Ryzen 7 5700U (Zen 2), 8 cores, AVX2 only, no AVX-512 |
| Compiler | gcc/g++ 12.4 |
| GLM | 0.9.9.8 (Ubuntu `libglm-dev`); the report used GLM 1.1.0 at `6f14f47` |
| Sizes | 8192, 262144, 2097152; 3 rounds x 11 paired reps, `taskset -c 0` |
| Corpus | 115 items: P = 78 scalar, B = 43 batch/SoA/stream/fast, S = 3 small-n, AA = 1 |
| Builds | `-O2`; `-O3 -march=native`; `-O3 -march=native -flto` (report's tiers) |
| Additional | the same three tiers with `-ffp-contract=off` (the library's documented exact build) |

The A/A band (same implementation, two TUs) was 0.73-1.29 at 262k and 0.98-1.07 at 2M on this
machine; the 262k bands are wide because the machine carries background load. Read a single cell as
"clears the same-round band", not as a precise multiplier.

## P group: scalar API

Ratio = t_GLM / t_fv; > 1 means flashvec faster, < 1 means GLM faster.

| Build | Median cell | Geometric mean | flashvec / tie / GLM |
|---|---|---|---|
| `-O2` (report flag) | 0.860 | 0.66x | 33 / 65 / 136 |
| `-O3 -march=native` | 0.817 | 0.61x | 4 / 72 / 158 |
| `-O3 -march=native -flto` | 0.828 | 0.62x | 16 / 79 / 139 |
| `-O2 -ffp-contract=off` | **0.995** | **0.98x** | 27 / 104 / 103 |
| `-O3 -march=native -ffp-contract=off` | **1.000** | **0.97x** | 55 / 133 / 46 |

The report measured a 0.40x median for the pre-P0 `.cpp` scalar API. P0 moved the scalar bodies into
the header, and the scalar median is now at parity when the library is built as documented
(`-ffp-contract=off`). The default-flag rows are lower because a scalar body that contains a
contractible `a + b*c` cannot be forced inline under a consumer's `-ffp-contract=fast` without
letting the caller contract it into an FMA; it therefore stays an out-of-line call to keep the
bit-exact contract. This is a deliberate trade, not a call-boundary bug.

Remaining GLM wins in the `-ffp-contract=off` P group (median of 3 rounds):

| Item | O2 8192 / 2M | O3n 8192 / 2M | Why |
|---|---|---|---|
| `mat3_inverse` | 0.44 / 0.73 | 0.37 / 0.76 | scalar inverse with a singular guard; GLM inlines a branch-light version |
| `mat4_inverse` | 0.74 / 0.78 | 0.72 / 0.76 | general 4x4 adjugate; GLM has a shorter specialized path |
| `vec2_project` | 0.84 / 0.98 | 0.24 / 0.91 | the zero-onto guard blocks vectorization at O3 (`vec3_project` 0.81) |
| `quat_slerp` | 0.89 / 0.90 | 0.88 / 0.87 | exact sin/acos path; GLM's slerp is cheaper |
| `quat_from_mat3` | 1.00 / 1.02 | 0.49 / 0.75 | branch-heavy Shepperd conversion |

These are expensive scalar functions; the batch/SoA entries are the performance path for them
(`mat3_inverse_soa` below is 4.5-4.7x faster than GLM's per-element loop at 8192).

## B group: batch / SoA / stream / fast

`B` at n=8192 is the L2-resident tier; per-size medians:

| Build | n=8192 | n=262144 | n=2097152 | wins at 8192 |
|---|---|---|---|---|
| `-O2` | **2.94x** | 1.09x | 1.05x | 37 / 43 |
| `-O3 -march=native` | 1.23x | 1.06x | 1.01x | 24 / 43 |
| `-O3 -march=native -flto` | 1.25x | 1.08x | 1.02x | 23 / 43 |
| `-O2 -ffp-contract=off` | **2.64x** | 1.05x | 1.05x | 37 / 43 |
| `-O3n -ffp-contract=off` | 1.11x | 1.11x | 1.02x | 23 / 43 |

Geometric means over the whole B group: 1.58x (`-O2`), 1.39x (`-O3n`), 1.40x (`-LTO`);
1.56x / 1.39x for the contract builds.

The report's three batch negatives are gone on this machine:

| Entry | Report (cloud AVX-512) | Local, GLM ratio at 8192 | Local, 2M |
|---|---|---|---|
| `m3_mul_soa` | 0.48x at `-O2` | **5.90x** (`-O2`) / **4.18x** (`-O3n`) | 3.88x / 1.27x |
| `m3_inverse_soa` | 0.68-0.98x | **4.59x** / **4.59x** | 2.78x / 2.62x |
| `v2_madd_normalize_soa` | 0.79x at `-O3n` | **6.55x** (`-O2`) / **7.50x** (`-O3n`) | 1.09x / 1.16x |

Other stable 8192-tier wins (default `-O2` / `-O3n`): `v2_len_soa` 7.73/7.78, `v2_normalize_soa`
7.50/7.50, `v2_dot_soa` 5.86/5.87, `q_mul_soa` 7.77/5.49, `q_rotate_soa` 5.84/5.51,
`m3_transform_soa` 4.47/4.35. At 2M the working set is DRAM-bound and many entries converge to
1.0-1.4x, except `m3_mul_soa` (3.88/1.27), `m3_inverse_soa` (2.78/2.62) and the quaternion entries
(1.2-1.4x).

`_stream` is the one place where the machine flips the sign: on this Ryzen it wins 1.1-1.7x at
262k/2M, while the report's bandwidth-rich cloud Xeon measured 0.47-1.03x. That confirms the
portability warning; the library still does not switch automatically.

## AVX-512 and the report's cloud

The report's cloud Xeon has AVX-512F and GLM auto-vectorizes 512-bit; this machine does not. The
optional AVX-512 backend for the fused madd-normalize entry is compile-checked (`zmm` /
`vrsqrt14ps`, and their removal under `FV_I_NO_AVX512`), but its runtime path has not executed on
AVX-512 hardware here. The cloud matrix therefore still needs a run on that machine; this local
matrix cannot confirm or refute the report's AVX-512 cells.

## Reproduction

The external `bench_all.cpp` and `summarize_all.py` are kept under `bench/glm/` (the summarizer
gained a CLI filename argument). GLM 0.9.9.8 headers came from Ubuntu `libglm-dev`; the report's GLM
1.1.0 is not installed here. The local raw file is `results_final.tsv` (default flags) plus
`results_final_contract.tsv` (`-ffp-contract=off`), not committed because it is machine-specific.
