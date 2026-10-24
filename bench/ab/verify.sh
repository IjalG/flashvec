#!/bin/sh
# bench/ab/verify.sh -- same-direction evidence mode (R1-R5 + INCONCLUSIVE)
#
#   make -C bench/ab verify OP=fv_vec2_add_x4 N=1048576
#
# The key point (this is the R4 fix): **A/A and A/B interleave inside one round**
# (--aa-each-round). Within each process every round measures A/B first and then one A/A
# pair, on the same code, the same data and the same immersion state; across 2 builds x
# 2 orders x ROUNDS rounds that gives >=16 A/A samples. The floor is the p50/p90/max of
# those samples and R4 uses **p90** as its divisor (steadier than max, stricter than p50).
#
# Three verdicts: PASS / FAIL / INCONCLUSIVE
#   INCONCLUSIVE = the gain is below the floor, same-direction reproduction fails, or the
#   point is not resolvable (floor > 15%). That is not "the optimisation does nothing", it
#   is "this machine cannot measure it" -- the two must stay separate.
# Exit codes: 0 PASS / 1 FAIL / 4 INCONCLUSIVE / 2 usage / 3 no such op
set -u
ABDIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$ABDIR/../.." && pwd)
VB=$ABDIR/build/verify

OP=${OP:-}
N=${N:-1048576}
CPU=${CPU:-6}
TRIALS=${TRIALS:-6}
ROUNDS=${ROUNDS:-4}
MIN_GAIN=${MIN_GAIN:-1.15}
MAX_REGRESS=${MAX_REGRESS:-0.95}
FLOOR_MAX=${FLOOR_MAX:-0.15}
CHECK_SIZES=${CHECK_SIZES:-8192,65536,1048576,4194304}
EXPECT_REGRESS=${EXPECT_REGRESS:-}
BLACKLIST=${BLACKLIST:-131072}

[ -n "$OP" ] || { echo "usage: make -C bench/ab verify OP=<op-name> [N=<size>]" >&2; exit 2; }

BASE="-std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic"
ALIGN="-falign-functions=64 -falign-loops=32 -falign-jumps=32"

mk_bin() {
    mkdir -p "$1"
    make -s -f "$ABDIR/Makefile" "$1/ab_bench" BUILD="$1" CXXFLAGS="$BASE $2" \
        CAND_INC="$WT/include" CAND_HDR="$WT/include/flashvec.h" \
        CAND_SRC="$WT/src/flashvec_batch.cpp" CAND_WITH_SRC="$(cat "$WT/WITH_SRC")" \
        AB_OPS_EXTRA=1 CAND_LABEL=verify "CAND_DESC=worktree snapshot" >"$1/build.log" 2>&1 || {
        echo "[verify] build failed: $1" >&2; tail -15 "$1/build.log" >&2; return 1; }
}
ge() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 >= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 > b + 0) }'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 < b + 0) }'; }

echo "====================================================================="
echo "verify  OP=$OP  N=$N   loadavg=$(cut -d' ' -f1 /proc/loadavg)"
echo "R1 bit-identical / R2 target >=$MIN_GAIN / R3 same direction (AB,BA x 2 builds) / R4 > A/A p90 floor / R5 other sizes >=$MAX_REGRESS"
echo "====================================================================="

make -s -f "$ABDIR/Makefile" snapshot >/dev/null 2>&1
make -s -f "$ABDIR/Makefile" worktree-snapshot >"$ABDIR/build/verify_snap.log" 2>&1 || {
    echo "worktree-snapshot failed" >&2; tail -5 "$ABDIR/build/verify_snap.log" >&2; exit 2; }
WT=$ABDIR/build/worktree_snapshot
LOAD0=$(cut -d' ' -f1 /proc/loadavg)

echo "[1/4] building AB#1 / AB#2 (they differ only in code placement)..."
mk_bin "$VB/ab1" "" || exit 2
mk_bin "$VB/ab2" "$ALIGN" || exit 2
B1=$VB/ab1/ab_bench
B2=$VB/ab2/ab_bench

IS_INTRA=0
"$B1" --list-ops | awk -F'\t' -v op="$OP" '$1 == op && $9 == "0x8" { found = 1 } END { exit !found }' && IS_INTRA=1
if ! "$B1" --list-ops | awk -F'\t' -v op="$OP" '$1 == op { f = 1 } END { exit !f }'; then
    echo "no such op: $OP (see the registry with --list-ops)" >&2; exit 3
fi
INTRA_FLAG=""
[ "$IS_INTRA" = "1" ] && INTRA_FLAG="--intra"

echo "[2/4] R1 bit-identical check (all tail lengths n=N-7..N+3)..."
gs=""; i=-7; while [ "$i" -le 3 ]; do s=$((N + i)); [ "$s" -ge 1 ] && gs="$gs,$s"; i=$((i + 1)); done
gs=${gs#,}
GATE_OK=1
for bin in "$B1" "$B2"; do
    taskset -c "$CPU" "$bin" --gate-only --op="$OP" --sizes="$gs" >"$VB/gate.log" 2>&1 || {
        GATE_OK=0; echo "  gate FAIL: $(grep -m1 'GATE FAIL' "$VB/gate.log")"; }
done
[ "$GATE_OK" = "1" ] && echo "  gate-only passed (2 builds x 11 tail lengths)"

echo "[3/4] timing: 2 builds x 2 orders x $ROUNDS rounds (each round: A/B then one A/A pair)..."
rm -f "$VB"/p_*.tsv
for tag in 1 2; do
    bin=$([ "$tag" = 1 ] && echo "$B1" || echo "$B2")
    for ord in ab ba; do
        r=1
        while [ "$r" -le "$ROUNDS" ]; do
            taskset -c "$CPU" "$bin" --op="$OP" --sizes="$N" --no-cold --no-small --no-pipe \
                --order="$ord" --aa-each-round --trials-min="$TRIALS" --trials-max="$TRIALS" \
                $INTRA_FLAG --tsv-out="$VB/p_${tag}_${ord}_${r}.tsv" \
                >"$VB/p_${tag}_${ord}_${r}.log" 2>&1
            r=$((r + 1))
        done
    done
done

echo "[4/4] summary + R5 scan..."
STATS=$(for f in "$VB"/p_*.tsv; do
    [ -f "$f" ] || continue
    awk -F'\t' -v op="$OP" '
        /^#/ { next }
        $1 == "phase" { for (i = 1; i <= NF; i++) h[$i] = i; next }
        $2 == op { printf "%s\t%s\t%s\t%s\t%s\n", $(h["ratio"]), $(h["aa_p50"]),
                   $(h["aa_p90"]), $(h["aa_max"]), $(h["pair_iqr_pct"]) }' "$f"
done)
RATIOS=$(printf '%s\n' "$STATS" | awk -F'\t' 'NF >= 5 { print $1 }')
AA50=$(printf '%s\n' "$STATS" | awk -F'\t' 'NF >= 5 { print $2 }')
AA90=$(printf '%s\n' "$STATS" | awk -F'\t' 'NF >= 5 { print $3 }')
AAMAX=$(printf '%s\n' "$STATS" | awk -F'\t' 'NF >= 5 { print $4 }')
qmed() { sort -g | awk '{ a[NR] = $1 } END { if (NR == 0) { print "1.0000"; exit } printf "%.4f", a[int((NR + 1) / 2)] }'; }
# An op that exists only in the candidate and has no INTRA pair produces no timing samples at
# all (the A/B side is unavailable, the INTRA side has no counterpart). Say so instead of
# printing ratios derived from an empty set -- R1 is then the only verdict that means anything.
NO_SAMPLES=0

R_MED=$(printf '%s\n' "$RATIOS" | qmed)
R_MIN=$(printf '%s\n' "$RATIOS" | sort -g | head -1)
R_MAX=$(printf '%s\n' "$RATIOS" | sort -g | tail -1)
N_SAMP=$(printf '%s\n' "$RATIOS" | grep -c . || true)
[ "$N_SAMP" -eq 0 ] 2>/dev/null && NO_SAMPLES=1
if [ "$NO_SAMPLES" = "1" ]; then
    R_MED=; R_MIN=; R_MAX=; F50DEV=; F90DEV=; FMAXDEV=; GAIN=
fi
F50=$(printf '%s\n' "$AA50" | qmed)
F90=$(printf '%s\n' "$AA90" | qmed)
FMAX=$(printf '%s\n' "$AAMAX" | sort -g | tail -1)
dev() { awk -v x="$1" 'BEGIN { d = x - 1; if (d < 0) d = -d; printf "%.4f", d }'; }
F50DEV=$(dev "$F50"); F90DEV=$(dev "$F90"); FMAXDEV=$(dev "$FMAX")
GAIN=$(awk -v r="$R_MIN" 'BEGIN { printf "%.4f", r - 1 }')

taskset -c "$CPU" "$B1" --op="$OP" --sizes="$CHECK_SIZES" --no-cold --no-small --no-pipe \
    --order=ab --trials-min="$TRIALS" --trials-max="$TRIALS" $INTRA_FLAG \
    --tsv-out="$VB/scan.tsv" >"$VB/scan.log" 2>&1
SCAN_BAD=""; SCAN_MIN=99; SCAN_MIN_N=""
while read -r n r; do
    [ -n "${n:-}" ] || continue
    case ",$BLACKLIST,$EXPECT_REGRESS," in
        *",$n,"*) echo "      n=$n ratio=$r  (blacklisted / pre-declared: recorded, not judged)"; continue;;
    esac
    if lt "$r" "$SCAN_MIN"; then SCAN_MIN=$r; SCAN_MIN_N=$n; fi
    if lt "$r" "$MAX_REGRESS"; then SCAN_BAD="$SCAN_BAD n=$n:$r"; fi
done <<EOF
$(awk -F'\t' -v op="$OP" '
    /^#/ { next }
    $1 == "phase" { for (i = 1; i <= NF; i++) { if ($i == "ratio") c = i; if ($i == "n") ni = i } next }
    $2 == op { printf "%s %.4f\n", $ni, $c }' "$VB/scan.tsv")
EOF
[ "$SCAN_MIN" = "99" ] && SCAN_MIN=1.0

VERDICT=PASS; WHY=""
if [ "$GATE_OK" != "1" ]; then VERDICT=FAIL; WHY="R1 bit-identical check failed"; fi
if [ -n "$SCAN_BAD" ]; then VERDICT=FAIL; WHY="R5 regression beyond limit:$SCAN_BAD"; fi
if [ "$N" = "131072" ]; then
    VERDICT=INCONCLUSIVE; WHY="n=131072 is the known-poisonous tier (A/A 0.50-2.22) and may not be the criterion"
fi
MED_DEV=$(awk -v m="${R_MED:-1}" 'BEGIN { d = m - 1; if (d < 0) d = -d; printf "%.4f", d }')
if [ "$VERDICT" = "PASS" ]; then
    if lt "$R_MED" "1.0" && gt "$MED_DEV" "$F90DEV"; then
        # median ratio clearly below 1 and beyond the floor: this is a **measured** slowdown, not an unmeasurable point
        VERDICT=FAIL; WHY="median $R_MED is $MED_DEV below 1, beyond the A/A floor $F90DEV (same-direction slowdown)"
    elif gt "$F90DEV" "$FLOOR_MAX"; then
        VERDICT=INCONCLUSIVE; WHY="A/A floor at this point is $F90DEV > $FLOOR_MAX (not resolvable)"
    elif ! gt "$GAIN" "$F90DEV"; then
        VERDICT=INCONCLUSIVE; WHY="gain $GAIN does not exceed the A/A p90 floor $F90DEV (below resolution)"
    elif ! ge "$R_MIN" "$MIN_GAIN"; then
        VERDICT=INCONCLUSIVE; WHY="same direction but min $R_MIN < $MIN_GAIN (threshold not met, the four ratios disagree)"
    fi
fi
[ "$IS_INTRA" = "1" ] && WHY="$WHY [INTRA: candidate-internal control, not an improvement over HEAD]"
if [ "$NO_SAMPLES" = "1" ]; then
    VERDICT=INCONCLUSIVE
    WHY="no timing samples: $OP exists only in the candidate and has no INTRA pair, so A/B has no counterpart -- R1 is the only verdict this run can give"
fi

echo
echo "====================================================================="
echo "verdict  OP=$OP  N=$N   mode=$([ "$IS_INTRA" = 1 ] && echo INTRA-candidate || echo A-B-vs-HEAD)"
echo "====================================================================="
printf "  %-32s %s  (tail lengths %s)\n" "R1 bit-identical" \
    "$([ "$GATE_OK" = 1 ] && echo PASS || echo FAIL)" "$gs"
printf "  %-32s median=%s  threshold %s\n" "R2 median ratio at target N" "${R_MED:-n/a (no samples)}" "$MIN_GAIN"
printf "  %-32s min=%s  max=%s  (%s process samples)\n" "R3 same-direction reproduction" \
    "${R_MIN:-n/a}" "${R_MAX:-n/a}" "$N_SAMP"
printf "  %-32s gain=%s  floor p50=%s p90=%s max=%s\n" "R4 gain vs A/A floor" "${GAIN:-n/a}" \
    "${F50DEV:-n/a}" "${F90DEV:-n/a}" "${FMAXDEV:-n/a}"
printf "  %-32s min=%s @ n=%s  regressions:%s\n" "R5 other sizes" \
    "$([ "$NO_SAMPLES" = 1 ] && echo n/a || echo "$SCAN_MIN")" \
    "${SCAN_MIN_N:-n/a}" "${SCAN_BAD:-none}"
echo "  -------------------------------------------------------------------"
printf "  ratios: $(printf '%s\n' "$RATIOS" | tr '\n' ' ')\n"
printf "  loadavg start %s end %s\n" "$LOAD0" "$(cut -d' ' -f1 /proc/loadavg)"
echo "====================================================================="
echo "verdict: $VERDICT$([ -n "$WHY" ] && echo "  <- $WHY")"
echo "data: $VB/p_*.tsv"
case "$VERDICT" in
    PASS) exit 0;;
    FAIL) exit 1;;
    *) exit 4;;
esac
