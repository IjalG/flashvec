#!/bin/sh
# verify/checks.sh —— 结构检查（验收者 F）：把「不能只靠跑测试证明」的规格条目变成可复现证据。
# 只读仓库文件，不修改任何东西。输出同时写 build/checks.txt。
set -u
cd "$(dirname "$0")" || exit 1
mkdir -p build
{
echo "== 结构检查 (只读)  $(date -Iseconds) =="
echo
echo "-- [1] include/flashvec.h 必须是纯 C：不得出现 intrinsics / C++ 语法 --"
if grep -n -E 'immintrin|xmmintrin|_mm_|__m256|__m128|namespace|template *<|std::|reinterpret_cast|static_cast|class [A-Za-z]' ../include/flashvec.h; then
  echo "  ✗ 上面是命中项"
else
  echo "  PASS: 无任何 intrinsics / C++ 语法"
fi
echo
echo "-- [2] 头文件不得含 target 属性或 builtin 依赖（应为纯声明） --"
grep -n -E 'target\(|__builtin|__attribute__' ../include/flashvec.h || echo "  PASS: 无匹配"
echo
echo "-- [3] 默认构建（根 Makefile）不得使用 -march=native / -mavx2 / -ffast-math --"
grep -n -E '^CXXFLAGS|^CFLAGS' ../Makefile
if grep -n -E '^[^#]*(march=native|mavx|ffast-math|Ofast)' ../Makefile; then
  echo "  ✗ 上面是有效行里的激进编译选项"
else
  echo "  PASS: 默认构建没有 -march=native / -mavx2 / -ffast-math（注释里提到的不算）"
fi
echo
echo "-- [4] AVX2 内核：target 属性单函数 + 运行时 __builtin_cpu_supports 派发（扫描 src/ 下所有 .cpp）--"
found=0
for f in ../src/*.cpp; do
  if grep -q -E 'target\(|__builtin_cpu_supports' "$f"; then
    found=1
    echo "  $f:"
    grep -n -E 'target\(|__builtin_cpu_supports' "$f" | sed 's/^/    /'
  fi
done
[ "$found" = 0 ] && echo "  ✗ 未找到 target 属性/运行时派发"
echo
echo "-- [5] 位级一致路径不得用 FMA / 倒数近似（扫描 src/ 下所有 .cpp）--"
if grep -n -E '_mm256_fmadd|_mm256_fmsub|_mm256_fnmadd|_mm256_fmaddsub|_mm256_rcp|_mm256_rsqrt|_mm_rcp|_mm_rsqrt' ../src/*.cpp; then
  echo "  ↑ 命中：需核对这些指令是否只出现在 dot_soa（允许）里"
else
  echo "  PASS: 无 FMA / 倒数近似 intrinsics"
fi
echo
echo "-- [6] 头文件是否写明 dot_soa 归约顺序免责 + 掩码语义（规格硬要求）--"
grep -n -E '归约顺序|掩码|lane|dot_soa' ../include/flashvec.h | sed 's/^/  /' | head -20
echo
echo "-- [7] 机器/编译器事实 --"
grep -m1 'model name' /proc/cpuinfo | sed 's/^/  /'
printf '  avx2/fma: '; grep -o -m1 -E 'avx2|fma' /proc/cpuinfo | sort -u | tr '\n' ' '; echo
lscpu 2>/dev/null | grep -E '^L[123]|^CPU\(s\)' | sed 's/^/  /'
c++ --version | head -1 | sed 's/^/  /'
echo
echo "-- [8] 编译产物指令普查（-O3 -march=native，两个 TU）--"
for OBJ in build/flashvec_O3n.o build/flashvec_batch_O3n.o; do
  if [ ! -f "$OBJ" ]; then echo "  ($OBJ 尚未构建，跳过)"; continue; fi
  echo "  --- $OBJ ---"
  for sym in $(nm --defined-only "$OBJ" | awk '$2=="T"{print $3}'); do
    cnt=$(objdump -d --no-show-raw-insn --disassemble="$sym" "$OBJ" 2>/dev/null | grep -oE '\b(vfmadd[0-9]*[sp][ds]|vfmsub[0-9]*[sp][ds]|vfnmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vdivps|vsqrtps|vaddps|vsubps|vmulps|vpermps|vpermilps|vpgatherdd|vpgatherdq|vmovups|vblendvps|vblendps|vandps|vorps|vcmpps)\b' | sort | uniq -c | tr '\n' ' ')
    [ -n "$cnt" ] && printf '  %-30s %s\n' "$sym" "$cnt"
  done
done
echo "  （未列出 = 该函数没有这些向量指令；FMA=0 是位级一致性的必要条件之一）"
echo
echo "-- [9] 不带 -march=native 的构建里 AVX2 内核是否仍在（内核是文件内静态符号，小写 t，入口函数只是派发）--"
for OBJ in build/flashvec_batch_O2.o build/flashvec_batch_O3.o; do
  if [ ! -f "$OBJ" ]; then echo "  ($OBJ 尚未构建，跳过)"; continue; fi
  echo "  --- $OBJ 的静态内核符号 ---"
  for sym in $(nm --defined-only "$OBJ" | awk '$2=="t"{print $3}' | c++filt); do
    cnt=$(objdump -d --no-show-raw-insn --disassemble="$sym" "$OBJ" 2>/dev/null | grep -oE '\b(vaddps|vsubps|vmulps|vdivps|vsqrtps|vfmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vmovups|vblendvps|vpermilps)\b' | sort | uniq -c | tr '\n' ' ')
    if [ -z "$cnt" ]; then cnt="（无向量指令）"; fi
    printf '    %-46s %s\n' "$sym" "$cnt"
  done
  echo "  --- 派发证据（入口函数 call 到哪个内核） ---"
  for fn in fv_vec2_add_x4 fv_vec2_add_soa fv_vec2_len_soa fv_vec2_normalize_soa fv_vec2_dot_soa; do
    printf '    %-26s -> %s\n' "$fn" "$(objdump -d --no-show-raw-insn --disassemble="$fn" "$OBJ" | grep -oE 'call.*<[^>]*>' | sed 's/call.*<//; s/>//' | c++filt | tr '\n' ' ')"
  done
done
echo
echo "-- [10] 默认档（-O3，不带 -march=native）实测加速比摘要（来自 make perf 的第三个配置）--"
if [ -f build/perf_O3_plain.txt ]; then sed -n '/=== N=8192/,/^$/p' build/perf_O3_plain.txt | head -14; else echo "  (build/perf_O3_plain.txt 尚未生成)"; fi
} 2>&1 | tee build/checks.txt
