#!/bin/sh
# bench/ab/verify.sh —— 同向证据模式（R1–R5 + INCONCLUSIVE）
#
#   make -C bench/ab verify OP=fv_vec2_add_x4 N=1048576
#
# 关键点（R4 的修法）：**A/A 与 A/B 在同一轮里交错**（--aa-each-round）。
# 每个进程内，每轮先测 A/B 再测一对 A/A，用同一段代码、同一份数据、同一浸入状态；
# 跨 2 套构建 × 2 个顺序 × ROUNDS 轮，得到 >=16 个 A/A 样本。地板取这些样本的
# p50/p90/max，R4 用 **p90** 当分母（比 max 稳，比 p50 严）。
#
# 结论三档：PASS / FAIL / INCONCLUSIVE
#   INCONCLUSIVE = 增益在地板以下、同向性不成立、或该点不可分辨（地板 > 15%）。
#   这不是"优化无效"，是"这台机器测不出来"——两者必须分开。
# 退出码：0 PASS / 1 FAIL / 4 INCONCLUSIVE / 2 用法 / 3 没有该 op
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

[ -n "$OP" ] || { echo "用法: make -C bench/ab verify OP=<op名> [N=<规模>]" >&2; exit 2; }

BASE="-std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic"
ALIGN="-falign-functions=64 -falign-loops=32 -falign-jumps=32"

mk_bin() {
    mkdir -p "$1"
    make -s -f "$ABDIR/Makefile" "$1/ab_bench" BUILD="$1" CXXFLAGS="$BASE $2" \
        CAND_INC="$WT/include" CAND_HDR="$WT/include/flashvec.h" \
        CAND_SRC="$WT/src/flashvec_batch.cpp" CAND_WITH_SRC="$(cat "$WT/WITH_SRC")" \
        AB_OPS_EXTRA=1 CAND_LABEL=verify "CAND_DESC=工作树快照" >"$1/build.log" 2>&1 || {
        echo "[verify] 构建失败 $1" >&2; tail -15 "$1/build.log" >&2; return 1; }
}
ge() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 >= b + 0) }'; }
gt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 > b + 0) }'; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a + 0 < b + 0) }'; }

echo "====================================================================="
echo "verify  OP=$OP  N=$N   loadavg=$(cut -d' ' -f1 /proc/loadavg)"
echo "R1 位级 / R2 目标 >=$MIN_GAIN / R3 同向(AB,BA × 2 构建) / R4 > A/A p90 地板 / R5 其它规模 >=$MAX_REGRESS"
echo "====================================================================="

make -s -f "$ABDIR/Makefile" snapshot >/dev/null 2>&1
make -s -f "$ABDIR/Makefile" worktree-snapshot >"$ABDIR/build/verify_snap.log" 2>&1 || {
    echo "worktree-snapshot 失败" >&2; tail -5 "$ABDIR/build/verify_snap.log" >&2; exit 2; }
WT=$ABDIR/build/worktree_snapshot
LOAD0=$(cut -d' ' -f1 /proc/loadavg)

echo "[1/4] 构建 AB#1 / AB#2（只差代码放置）……"
mk_bin "$VB/ab1" "" || exit 2
mk_bin "$VB/ab2" "$ALIGN" || exit 2
B1=$VB/ab1/ab_bench
B2=$VB/ab2/ab_bench

IS_INTRA=0
"$B1" --list-ops | awk -F'\t' -v op="$OP" '$1 == op && $9 == "0x8" { found = 1 } END { exit !found }' && IS_INTRA=1
if ! "$B1" --list-ops | awk -F'\t' -v op="$OP" '$1 == op { f = 1 } END { exit !f }'; then
    echo "没有这个 op：$OP（--list-ops 看注册表）" >&2; exit 3
fi
INTRA_FLAG=""
[ "$IS_INTRA" = "1" ] && INTRA_FLAG="--intra"

echo "[2/4] R1 位级对拍（全尾长 n=N-7…N+3）……"
gs=""; i=-7; while [ "$i" -le 3 ]; do s=$((N + i)); [ "$s" -ge 1 ] && gs="$gs,$s"; i=$((i + 1)); done
gs=${gs#,}
GATE_OK=1
for bin in "$B1" "$B2"; do
    taskset -c "$CPU" "$bin" --gate-only --op="$OP" --sizes="$gs" >"$VB/gate.log" 2>&1 || {
        GATE_OK=0; echo "  gate FAIL: $(grep -m1 'GATE FAIL' "$VB/gate.log")"; }
done
[ "$GATE_OK" = "1" ] && echo "  gate-only 通过（两套构建 × 11 个尾长）"

echo "[3/4] 计时：2 构建 × 2 顺序 × $ROUNDS 轮（每轮 A/B 后紧跟一对 A/A）……"
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

echo "[4/4] 汇总 + R5 扫描……"
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

R_MED=$(printf '%s\n' "$RATIOS" | qmed)
R_MIN=$(printf '%s\n' "$RATIOS" | sort -g | head -1)
R_MAX=$(printf '%s\n' "$RATIOS" | sort -g | tail -1)
N_SAMP=$(printf '%s\n' "$RATIOS" | grep -c . || true)
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
        *",$n,"*) echo "      n=$n ratio=$r  （黑名单/已声明：只记录，不作判据）"; continue;;
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
if [ "$GATE_OK" != "1" ]; then VERDICT=FAIL; WHY="R1 位级对拍失败"; fi
if [ -n "$SCAN_BAD" ]; then VERDICT=FAIL; WHY="R5 回退超限:$SCAN_BAD"; fi
if [ "$N" = "131072" ]; then
    VERDICT=INCONCLUSIVE; WHY="n=131072 是已知最毒的一档（A/A 0.50–2.22），不得当判据"
fi
MED_DEV=$(awk -v m="$R_MED" 'BEGIN { d = m - 1; if (d < 0) d = -d; printf "%.4f", d }')
if [ "$VERDICT" = "PASS" ]; then
    if lt "$R_MED" "1.0" && gt "$MED_DEV" "$F90DEV"; then
        # 中位比显著低于 1 且超出地板：这是**测出来的变慢**，不是测不出来
        VERDICT=FAIL; WHY="中位比 $R_MED 比 1 低 $MED_DEV，超过 A/A 地板 $F90DEV（同向变慢）"
    elif gt "$F90DEV" "$FLOOR_MAX"; then
        VERDICT=INCONCLUSIVE; WHY="该点 A/A 地板 $F90DEV > $FLOOR_MAX（不可分辨）"
    elif ! gt "$GAIN" "$F90DEV"; then
        VERDICT=INCONCLUSIVE; WHY="增益 $GAIN 未超过 A/A p90 地板 $F90DEV（分辨率以下）"
    elif ! ge "$R_MIN" "$MIN_GAIN"; then
        VERDICT=INCONCLUSIVE; WHY="同向但最小 $R_MIN < $MIN_GAIN（未达门槛，四个比值不一致）"
    fi
fi
[ "$IS_INTRA" = "1" ] && WHY="$WHY [INTRA：candidate 内部对照，不是相对 HEAD 的改进]"

echo
echo "====================================================================="
echo "判定  OP=$OP  N=$N   mode=$([ "$IS_INTRA" = 1 ] && echo INTRA-candidate || echo A-B-vs-HEAD)"
echo "====================================================================="
printf "  R1 位级对拍（全尾长 %s）  %s\n" "$gs" "$([ "$GATE_OK" = 1 ] && echo PASS || echo FAIL)"
printf "  R2 目标规模中位比        median=%s  门槛 %s\n" "$R_MED" "$MIN_GAIN"
printf "  R3 同向复现 4 比值       min=%s  max=%s  (%s 个进程样本)\n" "$R_MIN" "$R_MAX" "$N_SAMP"
printf "  R4 增益 vs A/A 地板      gain=%s  floor p50=%s p90=%s max=%s\n" "$GAIN" "$F50DEV" "$F90DEV" "$FMAXDEV"
printf "  R5 其它规模              min=%s @ n=%s  回退:%s\n" "$SCAN_MIN" "${SCAN_MIN_N:-—}" "${SCAN_BAD:-无}"
echo "  -------------------------------------------------------------------"
printf "  四个比值: $(printf '%s\n' "$RATIOS" | tr '\n' ' ')\n"
printf "  loadavg 起 %s 止 %s\n" "$LOAD0" "$(cut -d' ' -f1 /proc/loadavg)"
echo "====================================================================="
echo "结论: $VERDICT$([ -n "$WHY" ] && echo "  <- $WHY")"
echo "数据: $VB/p_*.tsv"
case "$VERDICT" in
    PASS) exit 0;;
    FAIL) exit 1;;
    *) exit 4;;
esac
