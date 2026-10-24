#!/bin/sh
# bench/ab/campaign.sh -- one measurement campaign covering every new entry point
# (_stream / _mask_inplace / vec3 / vec4)
#
#   sh bench/ab/campaign.sh
#
# Quiet-window criterion (the external load on this machine is **permanent**, so loadavg
# cannot be the gate):
#   pgrep -f 'cc1plus|cc1 |make -C bench/ab|g\+\+ .*flashvec' empty for 60 s in a row.
#   loadavg is only recorded, never used to reject; only a rise > 1.0 during the campaign
#   is reported.
#   No clean window within 10 minutes -> run anyway, but tag every verdict POLLUTED (never
#   use those numbers to decide whether a change stays).
#
# Verdict per (op, size): R1 bit-identical / R2 median ratio / R3 same direction (2 builds
# x 2 orders) / R4 gain vs the same-round interleaved A/A floor / R5 other sizes; three
# verdicts PASS / FAIL / INCONCLUSIVE. 131072 is recorded but never judged (A/A there runs
# 0.50-2.22).
set -u
ABDIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$ABDIR/../.." && pwd)
RAW=$ROOT/docs/research/raw
VB=$ABDIR/build/campaign
CPU=${CPU:-6}
TRIALS=${TRIALS:-5}
WAIT_MAX=${WAIT_MAX:-600}
POLLUTED=0

BASE="-std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic"
ALIGN="-falign-functions=64 -falign-loops=32 -falign-jumps=32"
BLACKLIST="131072"

busy() { pgrep -f 'cc1plus|cc1 |make -C bench/ab|g\+\+ .*flashvec' >/dev/null 2>&1; }

echo "====================================================================="
echo "campaign  $(date '+%%F %%H:%%M:%%S')  loadavg=$(cut -d' ' -f1 /proc/loadavg)  nproc=$(nproc)"
echo "window criterion: none of our compile processes (cc1plus/cc1/make -C bench/ab/g++ *flashvec) for 60 s"
echo "note: loadavg sits at 4-9 while cc1plus=0 -- the load is a resident external process, not our build."
echo "====================================================================="

wait_clean() {
    waited=0
    while [ "$waited" -lt "$WAIT_MAX" ]; do
        if busy; then sleep 5; waited=$((waited + 5)); continue; fi
        clean=1; k=0
        while [ "$k" -lt 12 ]; do
            sleep 5; k=$((k + 1))
            if busy; then clean=0; break; fi
        done
        if [ "$clean" = "1" ]; then echo "  window is clean (waited $waited s)"; return 0; fi
        waited=$((waited + 60))
    done
    echo "  no clean window within $WAIT_MAX s -> running anyway, every verdict tagged POLLUTED"
    POLLUTED=1
    return 1
}

mk_bin() {
    mkdir -p "$1"
    make -s -f "$ABDIR/Makefile" "$1/ab_bench" BUILD="$1" CXXFLAGS="$BASE $2" \
        CAND_INC="$WT/include" CAND_HDR="$WT/include/flashvec.h" \
        CAND_SRC="$WT/src/flashvec_batch.cpp" CAND_WITH_SRC="$(cat "$WT/WITH_SRC")" \
        AB_OPS_EXTRA=1 CAND_LABEL=campaign "CAND_DESC=worktree snapshot" >"$1/build.log" 2>&1 || {
        echo "build failed: $1" >&2; tail -10 "$1/build.log" >&2; return 1; }
}

wait_clean || true
LOAD0=$(cut -d' ' -f1 /proc/loadavg)

echo "[0/5] freeze the sources + build both sides (they differ only in code placement)..."
make -s -f "$ABDIR/Makefile" snapshot >/dev/null 2>&1
make -s -f "$ABDIR/Makefile" worktree-snapshot >"$ABDIR/build/campaign_snap.log" 2>&1
WT=$ABDIR/build/worktree_snapshot
mk_bin "$VB/ab1" "" || exit 2
mk_bin "$VB/ab2" "$ALIGN" || exit 2
B1=$VB/ab1/ab_bench
B2=$VB/ab2/ab_bench
mkdir -p "$RAW" "$VB"

echo "[1/5] R1 bit-identical check (once per new entry point, all tail lengths)..."
GATE_TSV=$VB/gate.tsv
: > "$GATE_TSV"
gate_sizes() {
    g=""; i=-7
    while [ "$i" -le 3 ]; do s=$(($1 + i)); [ "$s" -ge 1 ] && g="$g,$s"; i=$((i + 1)); done
    echo "$g" | sed 's/^,//'
}
OPS_NEW=$("$B1" --list-ops | awk -F'\t' '!/^#/ && ($10 != "8" || $1 ~ /_stream|_inplace/) { print $1 }')
for op in $OPS_NEW; do
    ok=1
    for bin in "$B1" "$B2"; do
        taskset -c "$CPU" "$bin" --gate-only --op="$op" --sizes="$(gate_sizes 65536)" \
            >"$VB/gate_$op.log" 2>&1 || ok=0
    done
    printf '%s\t%s\n' "$op" "$ok" >> "$GATE_TSV"
done
echo "  $(grep -c '\t1$' "$GATE_TSV" || echo 0)/$(wc -l < "$GATE_TSV") ops passed"

run_group() { # $1=tag $2=size $3=extra args
    for tag in 1 2; do
        bin=$([ "$tag" = 1 ] && echo "$B1" || echo "$B2")
        for ord in ab ba; do
            taskset -c "$CPU" "$bin" --sizes="$2" --no-cold --no-small --no-pipe \
                --order="$ord" --aa-each-round --trials-min="$TRIALS" --trials-max="$TRIALS" \
                $3 --tsv-out="$VB/r_$1_$2_$tag$ord.tsv" >"$VB/r_$1_$2_$tag$ord.log" 2>&1
        done
    done
}

echo "[2/5] six _stream entries x 5 sizes (the plain-store version is the A side)..."
for n in 8192 131072 1048576 4194304 8388608; do run_group stream "$n" "--ops=_stream"; done

echo "[3/5] _mask_inplace x 5 mask densities x 4 sizes..."
for d in 0.01 0.05 0.20 0.50 0.95; do
    for n in 8192 65536 1048576 4194304; do
        run_group "mask_$d" "$n" "--ops=_inplace --mask-density=$d"
    done
done

echo "[4/5] vec4 / vec3 (INTRA: same-byte-volume control inside the candidate)..."
for n in 8192 131072 1048576 4194304 8388608; do
    run_group v4 "$n" "--intra --ops=fv_vec4"
    run_group v3 "$n" "--intra --ops=fv_vec3"
done

echo "[5/5] summary and verdicts..."
LOAD1=$(cut -d' ' -f1 /proc/loadavg)
cat "$VB"/r_*.tsv > "$VB/all.tsv" 2>/dev/null

awk -F'\t' -v black="$BLACKLIST" -v polluted="$POLLUTED" -v load0="$LOAD0" -v load1="$LOAD1" '
function dev(x) { d = x - 1; if (d < 0) d = -d; return d }
/^#/ { next }
$1 == "phase" { for (i = 1; i <= NF; i++) h[$i] = i; next }
$1 == "warm" {
    key = $2 "\t" $3 "\t" $(h["mode"])
    r[key] = r[key] " " $(h["ratio"])
    aa[key] = aa[key] " " $(h["aa_p90"])
    seen[key] = 1
    if ($(h["pair_iqr_pct"]) + 0 > iqr[key]) iqr[key] = $(h["pair_iqr_pct"])
}
END {
    printf "%-24s %-8s %-7s %6s %6s %6s %7s %6s  %s\n", "op", "n", "mode", "r_med", "r_min", "floor", "gain", "iqr%", "verdict"
    for (k in seen) {
        m = split(r[k], v, " "); cnt = 0
        for (i = 1; i <= m; i++) if (v[i] != "") { cnt++; val[cnt] = v[i] + 0 }
        for (i = 1; i <= cnt; i++) for (j = i + 1; j <= cnt; j++) if (val[j] < val[i]) { t = val[i]; val[i] = val[j]; val[j] = t }
        med = val[int((cnt + 1) / 2)]
        amin = val[1]
        m2 = split(aa[k], w, " "); c2 = 0
        for (i = 1; i <= m2; i++) if (w[i] != "") { c2++; av[c2] = w[i] + 0 }
        for (i = 1; i <= c2; i++) for (j = i + 1; j <= c2; j++) if (av[j] < av[i]) { t = av[i]; av[i] = av[j]; av[j] = t }
        fl = (c2 > 0) ? av[int((c2 + 1) / 2)] : 1.0
        fdev = dev(fl)
        gain = amin - 1
        split(k, kk, "\t")
        if (polluted) vd = "POLLUTED"
        else if (kk[2] == black) vd = "INCONCLUSIVE(blacklisted tier)"
        else if (med < 1.0 && dev(med) > fdev) vd = "FAIL(slower)"
        else if (fdev > 0.15) vd = "INCONCLUSIVE(floor>15%)"
        else if (gain <= fdev) vd = "INCONCLUSIVE(gain<=floor)"
        else if (amin < 1.15) vd = "INCONCLUSIVE(below 1.15)"
        else vd = "PASS"
        printf "%-24s %-8s %-7s %6.3f %6.3f %6.3f %+7.3f %6.1f  %s\n", kk[1], kk[2], kk[3], med, amin, fl, gain, iqr[k], vd
    }
}' "$VB/all.tsv" | sort -k1,1 -k2,2n | tee "$RAW/ab_campaign_verdict.txt"

echo
echo "loadavg start $LOAD0 end $LOAD1    POLLUTED=$POLLUTED"
if [ "$POLLUTED" = "1" ]; then
    echo "!!! no clean window this round: every verdict above is POLLUTED, **do not act on it** !!!"
fi
echo "raw data: $VB/*.tsv   verdict table: $RAW/ab_campaign_verdict.txt"
