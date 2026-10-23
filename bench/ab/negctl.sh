#!/bin/sh
# bench/ab/negctl.sh —— 负对照：证明夹具会"咬人"
#
# 做法：从 git HEAD 快照复制出一棵**私有**源码树，只改一处内核，
# 然后用 CAND_SRC/CAND_INC 指向它编 candidate，跑夹具，检查：
#   1) 退出码必须是 4（校验和门禁硬失败）；
#   2) 输出必须**指名**出错的函数；
#   3) 必须没有任何计时/吞吐结果（门禁在前，硬失败在前，不许出数字）。
# 注意：整个过程不碰仓库的 src/ 与 include/，所以没有"跑完还原"的动作——
# 压根没改过工作树。这也是比 sed 改工作树更稳的做法。
set -u

ABDIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$ABDIR/../.." && pwd)
BUILD=$ABDIR/build
SNAP=$ABDIR/snapshot
BIN=$BUILD/ab_bench
TREE=$BUILD/negctl_tree

fail=0
case_no=0

run_case() {
    case_no=$((case_no + 1))
    name=$1        # 人类可读的破坏说明
    expect=$2      # 必须被指名的函数
    shift 2
    # 余下参数是 sed 表达式（串）

    rm -rf "$TREE"
    mkdir -p "$TREE/src" "$TREE/include"
    cp "$SNAP/include/flashvec.h" "$TREE/include/flashvec.h"
    cp "$SNAP/src/flashvec_batch.cpp" "$TREE/src/flashvec_batch.cpp"
    cp "$SNAP/src/flashvec.cpp" "$TREE/src/flashvec.cpp"

    for e in "$@"; do
        sed -i "$e" "$TREE/src/flashvec_batch.cpp"
    done

    # 先确认破坏真的落到源码上了（否则这个负对照是空的）
    if cmp -s "$SNAP/src/flashvec_batch.cpp" "$TREE/src/flashvec_batch.cpp"; then
        echo "[FAIL] case$case_no ($name): sed 没有改动任何字节，负对照无效"
        fail=$((fail + 1))
        return
    fi

    if ! make -s -f "$ABDIR/Makefile" "$BIN" \
            CAND_INC="$TREE/include" CAND_HDR="$TREE/include/flashvec.h" \
            CAND_SRC="$TREE/src/flashvec_batch.cpp" CAND_WITH_SRC=1 \
            CAND_LABEL=negctl >"$BUILD/negctl_build.log" 2>&1; then
        echo "[FAIL] case$case_no ($name): 变异树编不过，无法验证门禁"
        tail -5 "$BUILD/negctl_build.log"
        fail=$((fail + 1))
        return
    fi

    "$BIN" --label=negctl --sizes=65536,131072 --small=8 --no-cold --no-pipe \
           --trials-min=5 --trials-max=5 --target-ms=2 >"$BUILD/negctl_run.log" 2>&1
    rc=$?

    ok=1
    [ "$rc" -eq 4 ] || { echo "  ! 期待退出码 4，实得 $rc"; ok=0; }
    grep -q "GATE FAIL" "$BUILD/negctl_run.log" || { echo "  ! 没有 GATE FAIL"; ok=0; }
    grep -q "$expect" "$BUILD/negctl_run.log" || { echo "  ! 没有指名 $expect"; ok=0; }
    grep -q "阶段二" "$BUILD/negctl_run.log" && { echo "  ! 竟然还打印了计时结果"; ok=0; }

    if [ "$ok" -eq 1 ]; then
        echo "[PASS] case$case_no ($name): 退出码 4，指名 $expect，未输出任何计时结果"
        grep -E "GATE FAIL|不一致|超出容差|超出容差" "$BUILD/negctl_run.log" | head -3 | sed 's/^/       /'
    else
        echo "[FAIL] case$case_no ($name)"
        fail=$((fail + 1))
        grep -E "GATE FAIL|不一致" "$BUILD/negctl_run.log" | head -5 | sed 's/^/       /'
    fi
}

echo "=== 负对照：候选实现被故意改错，夹具必须硬失败并指名函数 ==="

# 1) 掩码 blend 方向选反：掩码外的 lane 取 a+b，掩码内的取 a（语义完全反了）
run_case "add_x4_mask 的 blend 选反" "fv_vec2_add_x4_mask" \
    's/_mm256_blendv_ps(va, _mm256_add_ps(va, vb), m)/_mm256_blendv_ps(_mm256_add_ps(va, vb), va, m)/'

# 2) scale 内核算成加法（改错一位算子的经典形态）
run_case "scale_x4 内核 mul 改成 add" "fv_vec2_scale_x4" \
    's/_mm256_storeu_ps((float \*)(void \*)(o + i), _mm256_mul_ps(va, vs));/_mm256_storeu_ps((float *)(void *)(o + i), _mm256_add_ps(va, vs));/'

# 3) normalize 的零向量保护掩码取反：只影响 len<=EPS 的 lane（数据里每 512 个埋了一个零向量）
run_case "normalize_soa 的 keep 掩码取反" "fv_vec2_normalize_soa" \
    's/_CMP_NLE_UQ/_CMP_LE_OQ/'

rm -rf "$TREE"
make -s -f "$ABDIR/Makefile" "$BIN" CAND_LABEL=head >/dev/null 2>&1 || true

if [ "$fail" -ne 0 ]; then
    echo "=== 负对照：$fail 个 case 未按预期失败 ==="
    exit 1
fi
echo "=== 负对照：全部 $case_no 个 case 都按预期硬失败 ==="
