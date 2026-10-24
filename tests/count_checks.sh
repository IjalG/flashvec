#!/bin/sh
# Totals the assertions reported by "make test-all", so the count can be compared before/after a change with one command:
#   bash tests/count_checks.sh              run make test-all and total its output
#   bash tests/count_checks.sh saved.log    total an already captured run (the "before" baseline)
#
# Three kinds of number are summed, because the fixtures report in three shapes and none of them is redundant:
#   "N checks"   -- assertions inside a test program (test_flashvec, test_batch, test_stream, ...)
#   "cases=N"    -- cases inside a parity stream (each case is a recorded input corpus)
#   "PASS" lines -- verdicts printed by the Makefile-level cross-form comparisons (cmp of two record streams)
# FAIL lines are counted separately and make the script exit non-zero.
set -u

if [ "$#" -ge 1 ]; then
    LOG=$1
    OUT=$(cat "$LOG")
    RC=0
else
    ROOT=$(cd "$(dirname "$0")/.." && pwd)
    OUT=$(cd "$ROOT" && make test-all 2>&1)
    RC=$?
    printf '%s\n' "$OUT" | tail -25
fi

CHECKS=$(printf '%s\n' "$OUT" | grep -oE '[0-9]+ checks' | grep -oE '[0-9]+' | awk '{s += $1} END {print s + 0}')
CASES=$(printf '%s\n' "$OUT" | grep -oE 'cases=[0-9]+' | grep -oE '[0-9]+' | awk '{s += $1} END {print s + 0}')
PASS=$(printf '%s\n' "$OUT" | grep -cE '^[[:space:]]*PASS')
FAIL=$(printf '%s\n' "$OUT" | grep -cE '^[[:space:]]*FAIL')
TOTAL=$((CHECKS + CASES + PASS))

echo "== make test-all totals =="
echo "  program assertions (\"N checks\") = $CHECKS"
echo "  parity cases (\"cases=N\")          = $CASES"
echo "  cross-form PASS verdicts          = $PASS"
echo "  FAIL verdicts                     = $FAIL"
echo "  TOTAL checks                      = $TOTAL"

if [ "$RC" -ne 0 ] || [ "$FAIL" -ne 0 ]; then
    echo "  VERDICT: FAILED"
    exit 1
fi
echo "  VERDICT: all green"
exit 0
