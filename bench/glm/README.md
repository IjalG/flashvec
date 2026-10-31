# GLM comparison harness

This directory carries the external 115-item GLM benchmark that produced
`docs/glm-rebench.md`. The source was supplied with the report; it is kept here so the local numbers
can be reproduced. It is not part of `make test-all`.

- `bench_all.cpp` — the P / B / S / AA cases and the timing loop.
- `aa_a.cpp`, `aa_b.cpp` — the A/A control: the same `fv_vec4_add_x2` in two TUs.
- `summarize_all.py` — the report's summarizer (accepts the TSV filename as `argv[1]`).
- `run_all.sh` — 3 tiers x 3 sizes x 3 rounds, `taskset -c 0`.

GLM is **not vendored**. The local run used Ubuntu's `libglm-dev` 0.9.9.8; the external report used
GLM 1.1.0 (`6f14f47`). Point `GLM_INCLUDE` at a GLM include directory:

```bash
GLM_INCLUDE=/path/to/glm/include bash bench/glm/run_all.sh
python3 bench/glm/summarize_all.py bench/glm/results_local.tsv
```

The script builds `all_O2`, `all_O3` and `all_LTO` from the current `include/` + `src/`. For the
documented exact scalar build, pass `EXTRA_FLAGS=-ffp-contract=off`; that is the fair P-group
comparison after P0. The raw TSV is machine-specific and is not committed.
