# bench/table.awk -- 从 bench 原始日志的 CASE 行生成 Markdown 对比表
# 用法: awk -f table.awk ../docs/research/raw/bench_*.txt
# 说明: 同一 (cfg, tier, op, n) 只保留第一次出现的值；加速比一律相对 tier1 @ O3_native。

{
  if ($1 != "CASE") next;
  for (i = 2; i <= NF; i++) {
    if (split($i, a, "=") == 2) kv[a[1]] = a[2];
  }
  k = kv["cfg"] "|" kv["tier"] "|" kv["op"] "|" kv["n"];
  if (!(k in ns_elem)) {
    ns_elem[k] = kv["ns_elem"] + 0;
    gbps[k] = kv["gbps"] + 0;
  }
  order[kv["n"]] = 1;
}

END {
  split("add scale dot dot_elem len normalize collision", opv, " ");
  split("t1_O3native|tier1|O3_native t1_ffast|tier1|O3_native_ffastmath t1_O2shipped|tier1|O2_shipped " \
        "t2_O3native|tier2|O3_native t2_ffast|tier2|O3_native_ffastmath t2_noerrno|tier2|O3_native_nomatherrno " \
        "t3_avx2|tier3|O3_native t3_ffast|tier3|O3_native_ffastmath t3_noerrno|tier3|O3_native_nomatherrno " \
        "t3_hadd|tier3hadd|O3_native t3_4acc|tier3_4acc|O3_native t3_rsqrt|tier3_rsqrt|O3_native", colv, " ");

  nn = asorti(order, sizes);
  for (s = 1; s <= nn; s++) {
    n = sizes[s];
    printf "\n### N = %s（元素/次，100 次重复取均值，5 次取中位数）\n\n", n;
    printf "| op |";
    for (c = 1; c <= 13; c++) { split(colv[c], cc, "|"); printf " %s ns/elem |", cc[1]; }
    printf " t3 GB/s |\n|---|";
    for (c = 1; c <= 13; c++) printf "---|";
    printf "---|\n";
    for (o = 1; o <= 7; o++) {
      op = opv[o];
      base = ns_elem["O3_native|tier1|" op "|" n];
      printf "| %s |", op;
      for (c = 1; c <= 13; c++) {
        split(colv[c], cc, "|");
        v = ns_elem[cc[3] "|" cc[2] "|" op "|" n];
        if (v == "" ) { printf " - |"; continue; }
        if (base != "" && base > 0) printf " %.3f (x%.2f) |", v, base / v;
        else printf " %.3f |", v;
      }
      g = gbps["O3_native|tier3|" op "|" n];
      printf " %s |\n", (g == "" ? "-" : sprintf("%.1f", g));
    }
  }
}
