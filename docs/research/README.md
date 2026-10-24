# Research notes

These four documents are the raw research that preceded the library design. They were
written in Chinese and have not been translated: they are working notes, kept for
provenance rather than as user documentation. The English documents that summarize their
conclusions are [`../optimization-report.md`](../optimization-report.md) and
[`../implementation-log.md`](../implementation-log.md).

| File | What it covers |
|---|---|
| `simd-hardware.md` | Anatomy of reference SIMD libraries (GLM, DirectXMath, highway), what GCC will and will not auto-vectorize, and the specific artefacts found in their disassembly. |
| `algorithm-tricks.md` | Approximation strategies, numerical-stability notes, and an assessment of which competitive-programming tricks transfer to a real simulation workload. |
| `layout-memory.md` | AoS / SoA / AoSoA trade-offs against real engine code, memory-hierarchy measurements, and the three candidate batch API routes. |
| `measured.md` | The three-implementation benchmark: scalar AoS API, compiler auto-vectorization, hand-written AVX2, across three working-set sizes. |

The raw logs these notes cite live in [`raw/`](raw/). They are the output of `make run`
and are **overwritten on every run**, so quote them together with a timestamp.

Anything in this directory reflects measurements taken at the time of writing and may not
match the current tree. For current numbers, run `make bench-ab`.
