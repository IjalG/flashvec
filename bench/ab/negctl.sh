#!/bin/sh
# bench/ab/negctl.sh -- negative control: prove the fixture actually bites
#
# How: copy a **private** source tree out of the git HEAD snapshot, break exactly one
# kernel, then point CAND_SRC/CAND_INC at it to build the candidate, run the fixture and
# check that:
#   1) the exit code is 4 (the checksum gate fails hard);
#   2) the output **names** the broken function;
#   3) there is no timing/throughput output at all (gate first, hard failure first, no
#      numbers allowed).
# Note: the repository's src/ and include/ are never touched, so there is no "restore after
# the run" step -- the worktree was never modified. That is also more robust than sed-ing
# the worktree.
set -u

ABDIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$ABDIR/../.." && pwd)
BUILD=$ABDIR/build
SNAP=$ABDIR/snapshot
BIN=$BUILD/ab_bench
TREE=$BUILD/negctl_tree
# Header of the per-op timing table; see the check below for why this token.
TIMING_TOKEN=realB/e

fail=0
case_no=0

run_case() {
    case_no=$((case_no + 1))
    name=$1        # human-readable description of the breakage
    expect=$2      # the function that must be named
    shift 2
    # the remaining arguments are sed expressions (one string each)

    rm -rf "$TREE"
    mkdir -p "$TREE/src" "$TREE/include"
    cp "$SNAP/include/flashvec.h" "$TREE/include/flashvec.h"
    cp "$SNAP/src/flashvec_batch.cpp" "$TREE/src/flashvec_batch.cpp"
    cp "$SNAP/src/flashvec.cpp" "$TREE/src/flashvec.cpp"

    # The kernels live in the header after the header-inline refactor and in
    # src/flashvec_batch.cpp before it, so every expression is applied to both: the one that
    # matches is the one that mutates, and the other is a no-op.
    before=$(cat "$TREE/include/flashvec.h" "$TREE/src/flashvec_batch.cpp" | cksum)
    for e in "$@"; do
        sed -i "$e" "$TREE/include/flashvec.h"
        sed -i "$e" "$TREE/src/flashvec_batch.cpp"
    done
    after=$(cat "$TREE/include/flashvec.h" "$TREE/src/flashvec_batch.cpp" | cksum)

    # first confirm the breakage really landed in the source (otherwise the case is void)
    if [ "$before" = "$after" ]; then
        echo "[FAIL] case$case_no ($name): sed changed no byte, negative control is void"
        fail=$((fail + 1))
        return
    fi

    if ! make -s -f "$ABDIR/Makefile" "$BIN" \
            CAND_INC="$TREE/include" CAND_HDR="$TREE/include/flashvec.h" \
            CAND_SRC="$TREE/src/flashvec_batch.cpp" CAND_WITH_SRC=1 \
            CAND_LABEL=negctl >"$BUILD/negctl_build.log" 2>&1; then
        echo "[FAIL] case$case_no ($name): the mutated tree does not build, gate unproven"
        tail -5 "$BUILD/negctl_build.log"
        fail=$((fail + 1))
        return
    fi

    "$BIN" --label=negctl --sizes=65536,131072 --small=8 --no-cold --no-pipe \
           --trials-min=5 --trials-max=5 --target-ms=2 >"$BUILD/negctl_run.log" 2>&1
    rc=$?

    ok=1
    [ "$rc" -eq 4 ] || { echo "  ! expected exit code 4, got $rc"; ok=0; }
    grep -q "GATE FAIL" "$BUILD/negctl_run.log" || { echo "  ! no GATE FAIL line"; ok=0; }
    grep -q "$expect" "$BUILD/negctl_run.log" || { echo "  ! $expect is not named"; ok=0; }
    # (3) no timing output after a gate failure. The marker is the header of the per-op
    # timing table ("realB/e" = real bytes per element). It is emitted only when measured
    # rows are about to be printed, never on the gate-failure path: the fixture prints the
    # hard-failure banner and returns 4 before that table (ab_bench.cpp, "if (any_gate_fail)").
    # The old marker here was a phrase the fixture never printed at all, so this check was
    # vacuous; TIMING_TOKEN is verified against a real timing run below to stay non-vacuous.
    grep -q "$TIMING_TOKEN" "$BUILD/negctl_run.log" && { echo "  ! timing results were printed anyway"; ok=0; }

    if [ "$ok" -eq 1 ]; then
        echo "[PASS] case$case_no ($name): exit 4, named $expect, no timing output"
        grep -E "GATE FAIL|mismatch|out of tolerance" "$BUILD/negctl_run.log" | head -3 | sed 's/^/       /'
    else
        echo "[FAIL] case$case_no ($name)"
        fail=$((fail + 1))
        grep -E "GATE FAIL|mismatch" "$BUILD/negctl_run.log" | head -5 | sed 's/^/       /'
    fi
}

echo "=== negative control: candidate broken on purpose; the fixture must fail hard and name the function ==="

# The "no timing output" assertion is only meaningful if TIMING_TOKEN really shows up when
# timing output happens. Prove it once per run on a fast, all-good configuration: if the
# token ever disappears (a renamed column, a moved table) this control would silently go
# vacuous again, so that must fail the control itself.
"$BIN" --label=negctl_probe --sizes=4096 --small=4 --no-cold --no-pipe \
       --trials-min=5 --trials-max=5 --target-ms=1 >"$BUILD/negctl_probe.log" 2>&1
if grep -q "$TIMING_TOKEN" "$BUILD/negctl_probe.log"; then
    echo "  probe: timing token '$TIMING_TOKEN' is present in a real timing run (check is live)"
else
    echo "[FAIL] timing token '$TIMING_TOKEN' never appears in a real timing run:"
    echo "       the 'no timing output after a gate failure' assertion would be vacuous"
    exit 1
fi

# 1) mask blend selects the wrong way: lanes outside the mask get a+b, lanes inside get a
run_case "add_x4_mask picks the wrong blend" "fv_vec2_add_x4_mask" \
    's/_mm256_blendv_ps(va, _mm256_add_ps(va, vb), m)/_mm256_blendv_ps(_mm256_add_ps(va, vb), va, m)/'

# 2) the scale kernel adds instead of multiplying (the classic one-token operator bug)
run_case "scale_x4 kernel mul becomes add" "fv_vec2_scale_x4" \
    's/_mm256_storeu_ps((float \*)(void \*)(o + i), _mm256_mul_ps(va, vs));/_mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vs));/'

# 3) normalize inverts its zero-vector guard: only len<=EPS lanes change (the corpus holds
#    one exact zero vector every 512 elements)
run_case "normalize_soa inverts its keep mask" "fv_vec2_normalize_soa" \
    's/_CMP_NLE_UQ/_CMP_LE_OQ/'

rm -rf "$TREE"
make -s -f "$ABDIR/Makefile" "$BIN" CAND_LABEL=head >/dev/null 2>&1 || true

if [ "$fail" -ne 0 ]; then
    echo "=== negative control: $fail case(s) did not fail as expected ==="
    exit 1
fi
echo "=== negative control: all $case_no cases failed hard as expected ==="
