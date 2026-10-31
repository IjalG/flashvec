#!/bin/bash
# Build and run the 115-item GLM matrix: 3 tiers x 3 sizes x 3 rounds.
# GLM_INCLUDE must point at a GLM include directory (glm/glm.hpp underneath).
# EXTRA_FLAGS can carry e.g. -ffp-contract=off for the documented exact scalar build.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
GLM_INCLUDE=${GLM_INCLUDE:?set GLM_INCLUDE to the GLM include directory}
CXX=${CXX:-g++}
EXTRA=${EXTRA_FLAGS:-}
COMMON="-std=c++17 -DNDEBUG -I$ROOT/include -I$GLM_INCLUDE $EXTRA"

"$CXX" $COMMON -O2 -DBUILD_TAG='"all_O2"' -o "$HERE/all_O2" \
    "$HERE/bench_all.cpp" "$HERE/aa_a.cpp" "$HERE/aa_b.cpp" "$ROOT"/src/*.cpp
"$CXX" $COMMON -O3 -march=native -DBUILD_TAG='"all_O3"' -o "$HERE/all_O3" \
    "$HERE/bench_all.cpp" "$HERE/aa_a.cpp" "$HERE/aa_b.cpp" "$ROOT"/src/*.cpp
"$CXX" $COMMON -O3 -march=native -flto -DBUILD_TAG='"all_LTO"' -o "$HERE/all_LTO" \
    "$HERE/bench_all.cpp" "$HERE/aa_a.cpp" "$HERE/aa_b.cpp" "$ROOT"/src/*.cpp

OUT=${1:-$HERE/results_local.tsv}
: > "$OUT"
for r in 1 2 3; do
  for b in all_O2 all_O3 all_LTO; do
    for N in 8192 262144 2097152; do
      echo "### round=$r build=$b N=$N" >> "$OUT"
      taskset -c 0 "$HERE/$b" "$N" 11 >> "$OUT" 2>> "$HERE/run_all.err"
    done
  done
done
