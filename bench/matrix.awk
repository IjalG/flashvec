# bench/matrix.awk -- 从原始日志生成「紧凑数字矩阵」（报告里的表格以此为准）
# 用法: awk -F'[ =]+' -f matrix.awk ../docs/research/raw/bench_O3_native.txt \
#          ../docs/research/raw/bench_O3_native_ffastmath.txt \
#          ../docs/research/raw/bench_O3_native_nomatherrno.txt \
#          ../docs/research/raw/bench_O2_shipped.txt
# 口径: 单元格是 5 次中位数的 ns/elem（ns_elem），加速比一律相对 tier1 @ -O3 -march=native。
/^CASE/ {
  for (i = 1; i <= NF; i++) {
    if ($i == "cfg") c = $(i + 1);
    if ($i == "tier") t = $(i + 1);
    if ($i == "op") o = $(i + 1);
    if ($i == "n") n = $(i + 1);
    if ($i == "ns_elem") v = $(i + 1);
    if ($i == "ns_elem_min") m = $(i + 1);
    if ($i == "gbps") g = $(i + 1);
  }
  key = c "|" t "|" o "|" n;
  if (!(key in seen)) { seen[key] = 1; med[key] = v; mn[key] = m; gb[key] = g; }
}

function cell(cfg, tier) {
  const_ = med[cfg "|" tier "|" o "|" n];
  return (const_ == "" ? "-" : const_);
}

function ratio(cfg, tier,   k) {
  k = med[cfg "|" tier "|" o "|" n];
  return (k > 0 && k != "" ? sprintf("x%.2f", b / k) : "-");
}

END {
  split("1048576 65536 8192", sizes, " ");
  split("add scale dot dot_elem len normalize collision", ops, " ");
  for (s = 1; s <= 3; s++) {
    n = sizes[s];
    printf "\n-- N=%s --\n", n;
    for (j = 1; j <= 7; j++) {
      o = ops[j];
      b = med["O3_native|tier1|" o "|" n] + 0;
      printf "%-9s t1=%-8s t2=%-8s t2ff=%-8s t2noerr=%-8s t3=%-8s t3hadd=%-8s t3acc4=%-8s t3rsqrt=%-8s t3_gbps=%-5s | %s %s %s %s %s\n",
        o,
        cell("O3_native", "tier1"), cell("O3_native", "tier2"),
        cell("O3_native_ffastmath", "tier2"), cell("O3_native_nomatherrno", "tier2"),
        cell("O3_native", "tier3"), cell("O3_native", "tier3hadd"),
        cell("O3_native", "tier3_4acc"), cell("O3_native", "tier3_rsqrt"),
        (gb["O3_native|tier3|" o "|" n] == "" ? "-" : gb["O3_native|tier3|" o "|" n]),
        ratio("O3_native", "tier2"), ratio("O3_native_ffastmath", "tier2"),
        ratio("O3_native_nomatherrno", "tier2"), ratio("O3_native", "tier3"),
        ratio("O3_native", "tier3acc4");
    }
  }
}
