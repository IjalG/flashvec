#!/bin/sh
# verify/checks.sh —— 结构检查（验收者 F）：把「不能只靠跑测试证明」的规格条目变成可复现证据。
# 只读被测文件；唯一的写入是 verify/build/ 下的中间产物与本报告。输出同时写 build/checks.txt。
#
# ============================================================================
# 本脚本在「批量入口搬进头文件」重构之后重写过。每条规则上方都写着「旧定义 → 新定义」
# 以及为什么必须改。重构带来的三个事实，规则全都围着它们改：
#
#   1) 批量入口的**实现**（含 AVX2 内核、target("avx2")、__builtin_cpu_supports 派发）
#      现在是 include/flashvec.h 里的 static inline（默认形态）；src/*.cpp 只剩 extern
#      回退符号（FV_USE_EXTERN / FV_BUILD 形态）。→ 原来扫 src/*.cpp 的规则扫不到了。
#   2) 头文件因此**允许**出现条件 intrinsics（FV_I_INTRINSICS && !FV_USE_EXTERN 时
#      include <immintrin.h>）。→「头文件里不得出现 immintrin」这条旧判据已不成立，
#      但「FV_USE_EXTERN 形态的包含集与重构前一致」这条更强的判据仍然成立且更可判定。
#   3) 旧规则对一个**已经过期的 build/*.o** 做指令普查。那些对象是重构前编的，
#      普查它们等于验收一个不存在的产物。→ 脚本自己按**最严配置**重新编两个形态再普查。
#
# 判定口径：11 条规则里 [7] 只打印机器事实（不设判据），其余每条都有 PASS / FAIL。
# 任何一条 FAIL → 退出码非 0。
# ============================================================================
set -u
cd "$(dirname "$0")" || exit 1
mkdir -p build
TMP=build/checks_tu
rm -rf "$TMP"
rm -f build/.checks_fails
mkdir -p "$TMP"
INC=../include
SRC=../src
RULES=0
FAILS=0

rule() { RULES=$((RULES + 1)); printf '\n-- [%s] %s --\n' "$1" "$2"; }
pass() { printf '  PASS: %s\n' "$*"; }
fail() { FAILS=$((FAILS + 1)); printf '  FAIL: %s\n' "$*"; }
info() { printf '  INFO: %s\n' "$*"; }

{
echo "== 结构检查 (只读)  $(date -Iseconds) =="
info "被测树：$(cd .. && pwd)"
info "规则一览：[1][2] 头文件纯净度 / [3] 编译选项 / [4] 派发位置 / [5] 源码级禁用内建 /"
info "  [6] 头文件契约文字 / [7] 机器事实 / [8] 产物级指令普查 / [9] 默认档内核+派发 /"
info "  [10] 默认档实测摘要 / [11] 已废弃内部名不得复活"

# ---------------------------------------------------------------- [1]
rule 1 "头文件纯净度 A：默认（头文件自带实现）形态在纯 C 编译器下可用，且 intrinsics 确实展开"
# 旧定义：include/flashvec.h 里不得出现 immintrin/_mm_/__m256/namespace/std:: 等（当时头文件是纯声明）。
# 新定义：实现搬进头文件后，头文件在 FV_I_INTRINSICS && !FV_USE_EXTERN 时**必须** include
#   <immintrin.h>。所以旧判据的前半句反转成正对照：
#     * immintrin 展开行数必须 > 0 —— 为 0 反而说明 inline 形态没生效（封装坏了）；
#     * 「纯 C 可用」改由更强的证据承担：cc -std=c11 -pedantic-errors 能编译通过。
if cc -std=c11 -pedantic-errors -Wall -Wextra -I"$INC" -c ../tests/test_header_c.c -o "$TMP/hdr_c.o" 2>"$TMP/hdr_c.err"; then
  pass "cc -std=c11 -pedantic-errors 可编译纯 C 消费者（头文件对 C 可用）"
else
  fail "纯 C 编译失败，前几行："
  sed -n '1,12p' "$TMP/hdr_c.err"
fi
EXPAND_C=$(cc -std=c11 -I"$INC" -E ../tests/test_header_c.c 2>/dev/null | grep -c immintrin)
EXPAND_PP=$(c++ -std=c++17 -I"$INC" -E ../tests/test_header_cpp.cpp 2>/dev/null | grep -c immintrin)
info "cc -std=c11 -E 里 immintrin 行数 = $EXPAND_C（>0 = 条件 intrinsics 生效，与重构后的预期一致）"
info "c++ -std=c++17 -E 里 immintrin 行数 = $EXPAND_PP"
if [ "$EXPAND_C" -gt 0 ] && [ "$EXPAND_PP" -gt 0 ]; then
  pass "默认形态确实展开 intrinsics（正对照；旧规则在这里要求 0，重构后该要求已失效）"
else
  fail "immintrin 展开行数不是正数：头文件的 inline 实现没有生效"
fi
if grep -n -E 'namespace|template *<|std::|reinterpret_cast|static_cast|class [A-Za-z]' "$INC/flashvec.h"; then
  fail "头文件出现 C++ 语法（上面是命中项）"
else
  pass "头文件无 C++ 语法（namespace/template</std::/cast/class）"
fi

# ---------------------------------------------------------------- [2]
rule 2 "头文件纯净度 B：FV_USE_EXTERN 形态的包含集必须与重构前一致"
# 旧定义：头文件不得含 target 属性 / __builtin / extern 实现（都应为纯声明）。
# 新定义：FV_USE_EXTERN 形态**就是**那个「纯声明」形态，把旧规则翻译成可判定的形式：
#   它经预处理器拉进来的头文件里不得出现 math.h / string.h / immintrin.h / stdlib.h。
# 为什么用包含集而不是 grep：头文件文本里的 intrinsics 是条件编译的，只有列出预处理器
#   实际拉进来的头文件，才能证明「回退形态的包含集与重构前逐字相同」。
DEPS=$(printf '#include "flashvec.h"\n' | cc -std=c11 -DFV_USE_EXTERN -I"$INC" -M -x c - 2>/dev/null)
HDRS=$(printf '%s\n' "$DEPS" | tr '[:space:]' '\n' | sed 's/\\$//' | grep -oE '[^/]+\.h$' | sort -u)
info "FV_USE_EXTERN 包含集：$(printf '%s' "$HDRS" | tr '\n' ' ')"
HIT=$(printf '%s\n' "$HDRS" | grep -xE 'math\.h|string\.h|immintrin\.h|stdlib\.h' | tr '\n' ' ')
if [ -n "$HIT" ]; then
  fail "FV_USE_EXTERN 形态多包含了：$HIT（重构前的包含集里没有这些）"
else
  pass "FV_USE_EXTERN 形态不含 math.h / string.h / immintrin.h / stdlib.h"
fi
if printf '#include "flashvec.h"\n' | cc -std=c11 -DFV_USE_EXTERN -pedantic-errors -I"$INC" -c -x c -o "$TMP/ext.o" - 2>"$TMP/ext.err"; then
  pass "FV_USE_EXTERN 形态可被纯 C 以 -pedantic-errors 编译"
else
  fail "FV_USE_EXTERN 形态编译失败，前几行："
  sed -n '1,12p' "$TMP/ext.err"
fi

# ---------------------------------------------------------------- [3]
rule 3 "默认构建不得使用 -march=native / -mavx2 / -ffast-math（判据落在**实际编译命令**上）"
# 旧定义：根 Makefile 里任何非注释行出现 -march=native/-mavx/-ffast-math/-Ofast 就 FAIL。
# 新定义：这条规则的**意图**是「默认产物保持 x86-64 基线」，但那个 grep 是过度近似 ——
#   它会连 FMA_FLAGS 一起判死，而 FMA_FLAGS 只被 tests/test_parity.c 与 tests/test_vec34_parity.c
#   的 *_cpp_fma.o 使用（消费者带 -march=native -mfma -ffp-contract=fast 的**敌意选项对拍**夹具），
#   恰恰是「证明默认档不带这些选项也安全」的那组测试。判它们失败是把夹具当成违规。
#   所以把判据从「扫文本」换成「问 make 默认目标到底会执行什么」：make -n -B all 打印的
#   每一条编译 src/ 的命令都必须干净。这比 grep 更严 —— 它包含变量展开后的真实结果。
if command -v make >/dev/null 2>&1; then
  # 注意 -C .. ：本脚本的工作目录是 verify/，而 -march=native 的判据在根 Makefile 上
  # （verify/Makefile 自己就有 -O3 -march=native 档，那是对拍用的，不是默认产物）。
  DEFAULT_CMDS=$(make -C .. -n -B all 2>/dev/null | grep -E -- '-c src/')
  if [ -z "$DEFAULT_CMDS" ]; then
    info "make -C .. -n -B all 没有打印出编译 src/ 的命令（构建形态可能变了），退化为只看变量行"
  else
    printf '%s\n' "$DEFAULT_CMDS" | sed 's/^/  /'
    if printf '%s\n' "$DEFAULT_CMDS" | grep -E 'march=native|mavx|ffast-math|Ofast|mfma'; then
      fail "默认构建命令里出现激进编译选项（上面是命中行）"
    else
      pass "根 Makefile 默认目标的每条编译命令都是 x86-64 基线（无 -march=native/-mavx/-ffast-math/-mfma）"
    fi
  fi
else
  info "没有 make，跳过默认构建命令检查"
fi
info "CXXFLAGS/CFLAGS 兜底行："
grep -n -E '^(CXXFLAGS|CFLAGS)[[:space:]]*[?:]?=' ../Makefile | sed 's/^/  /'
if grep -n -E '^(CXXFLAGS|CFLAGS)[[:space:]]*[?:]?=.*(march=native|mavx|ffast-math|Ofast)' ../Makefile; then
  fail "CXXFLAGS/CFLAGS 默认值里含激进编译选项"
else
  pass "CXXFLAGS/CFLAGS 默认值干净"
fi
info "Makefile 里 -march=native 的全部出现位置（应只在敌意选项对拍夹具里）："
grep -n 'march=native' ../Makefile | sed 's/^/  /'
# FMA_FLAGS 的每个使用者（建筑在它的 recipe 之上那行规则头）都必须是对拍夹具 *_fma。
# 用 awk 取「上一行」而不是 grep 相邻行：规则头与 recipe 是两行，grep 分不清归属。
FMA_HEADS=$(awk '/FMA_FLAGS\)/ {print prev} {prev=$0}' ../Makefile)
info "引用 FMA_FLAGS 的规则头：$(printf '%s' "$FMA_HEADS" | tr '\n' ' ')"
FMA_BAD=$(printf '%s\n' "$FMA_HEADS" | grep -c -v '_fma')
if [ "$FMA_BAD" -eq 0 ] && [ -n "$FMA_HEADS" ]; then
  pass "-march=native 只出现在敌意选项对拍夹具（*_cpp_fma），不进入任何默认产物"
else
  fail "FMA_FLAGS 有非对拍用途（$FMA_BAD 个规则头不含 _fma）：默认产物可能被污染"
fi

# ---------------------------------------------------------------- [4]
rule 4 "派发现在在 include/flashvec.h（target 属性 + 运行时 __builtin_cpu_supports）"
# 旧定义：扫描 src/ 下所有 .cpp 找 target 属性与 __builtin_cpu_supports。
# 新定义：实现（含内核与派发）搬进头文件后，src/*.cpp 里这两样一个都没有；判据改成在
#   include/flashvec.h 里找，并要求每一处都写成 != 0 —— 返回值是位掩码（avx2 = 0x400），
#   统一写法才能防止有人把它当布尔用。
CS_ALL=$(grep -c '__builtin_cpu_supports' "$INC/flashvec.h")
CS_NE0=$(grep -c '__builtin_cpu_supports("avx2") != 0' "$INC/flashvec.h")
TGT=$(grep -c 'target(' "$INC/flashvec.h")
SRC_CS=$(cat "$SRC"/*.cpp 2>/dev/null | grep -c '__builtin_cpu_supports')
info "include/flashvec.h：__builtin_cpu_supports = $CS_ALL（写成 != 0 的 = $CS_NE0）；target( = $TGT"
info "src/*.cpp：__builtin_cpu_supports = $SRC_CS（0 = 派发确实已全部搬进头文件，这是重构后的架构事实）"
if [ "$CS_ALL" -gt 0 ] && [ "$CS_NE0" -eq "$CS_ALL" ] && [ "$TGT" -gt 0 ]; then
  pass "头文件里有 target 属性与 $CS_ALL 处运行时派发，且全部写成 != 0"
else
  fail "头文件里的派发/属性不完整或写法不统一（total=$CS_ALL ne0=$CS_NE0 target=$TGT）"
fi
info "派发点位置（前 6 处，行号来自 include/flashvec.h）："
grep -n '__builtin_cpu_supports' "$INC/flashvec.h" | head -6 | sed 's/^/    /'

# ---------------------------------------------------------------- [5]
rule 5 "位级一致路径不得用 FMA / 倒数近似 / VMASKMOVPS（源码级，含头文件）"
# 旧定义：扫描 src/ 下所有 .cpp。
# 新定义：实现现在主要在 include/flashvec.h，所以扫描 include/ + src/ 两处。
#   扫的是**内建名**（_mm*_fmadd_ps 这类），不是汇编助记符：头文件的契约注释里会写
#   「vfmadd 计数为 0」「不用 rcpps/rsqrtps」这类句子，那是文档，不该被判失败。
#   其中 _mm*_maskload/maskstore（= VMASKMOVPS）是被实测否决的稀疏策略：Zen 2 上它对
#   未命中的 lane 仍然发生访存，只写命中 lane 的收益拿不到，所以连"允许"都不给。
BAN='_mm(128|256)?_(fmadd|fmsub|fnmadd|fmaddsub|fmsubadd|rcp|rsqrt|maskload|maskstore)'
if grep -rnE "$BAN" ../include ../src; then
  fail "上面是禁用内建名的命中项"
else
  pass "include/ + src/ 里 FMA / rcpps / rsqrtps / maskload / maskstore 内建名 = 0"
fi
DOC=$(grep -cE '\b(vfmadd|vrcpps|vrsqrtps|haddps|vmaskmovps)\b' "$INC/flashvec.h")
info "头文件里作为**注释文字**出现的助记符行数 = $DOC（允许：契约写「计数为 0」必须能提这些名字）"

# ---------------------------------------------------------------- [6]
rule 6 "头文件是否写明 dot_soa 归约顺序免责 + 掩码语义（规格硬要求，原样保留）"
# 先收进变量再打印前 20 行：`grep | head` 会让 grep 收到 SIGPIPE，在输出里留一条无意义的报错。
KW=$(grep -n -E '归约顺序|掩码|lane|dot_soa' "$INC/flashvec.h")
printf '%s\n' "$KW" | sed -n '1,20p' | sed 's/^/  /'
info "（上面只印前 20 行；全文 $(printf '%s\n' "$KW" | wc -l) 行）"
C_RED=$(grep -c '归约顺序' "$INC/flashvec.h")
C_MSK=$(grep -c '掩码' "$INC/flashvec.h")
C_DOT=$(grep -c 'dot_soa' "$INC/flashvec.h")
info "关键词计数：归约顺序=$C_RED 掩码=$C_MSK dot_soa=$C_DOT"
if [ "$C_RED" -gt 0 ] && [ "$C_MSK" -gt 0 ] && [ "$C_DOT" -gt 0 ]; then
  pass "三项契约文字均在头文件里"
else
  fail "头文件缺契约文字（归约顺序/掩码/dot_soa 之一）"
fi

# ---------------------------------------------------------------- [7]
rule 7 "机器/编译器事实（只打印，不设判据）"
grep -m1 'model name' /proc/cpuinfo | sed 's/^/  /'
printf '  avx2/fma: '; grep -o -m1 -E 'avx2|fma' /proc/cpuinfo | sort -u | tr '\n' ' '; echo
lscpu 2>/dev/null | grep -E '^L[123]|^CPU\(s\)' | sed 's/^/  /'
# 缓存层级必须按实测标：lscpu 的「8 MiB (2 instances)」是芯片总量，单线程只吃得到
# 自己那个 CCX 的 4 MiB。这条曾在项目文档里被写成「L3 8 MiB 共享」，并据此算错了拐点。
for i in 2 3; do
  if [ -r "/sys/devices/system/cpu/cpu0/cache/index$i/size" ]; then
    printf '  index%s: level=%s size=%s shared_cpu_list=%s\n' "$i" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/level)" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/size)" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/shared_cpu_list)"
  fi
done
info "口径：L3 4 MiB/CCX x 2 CCX = 8 MiB 芯片总量；单线程可见 4 MiB（本机拐点约 175K 元素 @24 B/元素）"
c++ --version | head -1 | sed 's/^/  /'
cc --version | head -1 | sed 's/^/  /'

# ---------------------------------------------------------------- [8]
rule 8 "产物级指令普查：7 项禁用助记符必须为 0，vmovntps 必须 > 0（NT store 正对照）"
# 旧定义：对 build/flashvec_O3n.o / build/flashvec_batch_O3n.o 逐符号普查。
# 新定义：内核在头文件里，而 build/ 下那些对象是重构前编的（普查它们等于验收不存在的产物）。
#   脚本自己按**最严配置**重新编两个形态：
#     形态 A 库内 extern：src/flashvec.cpp + flashvec_batch.cpp + flashvec_vec4.cpp
#            （vec3 的 SoA 与 vec4 的 _x2/_soa 的 extern 定义都在 flashvec_batch.cpp 里）
#     形态 B 头文件 inline 消费者：tests/test_stream.c（纯 C；它的热点只调库、不带自己的
#            算术，所以整 TU 普查合法。反例：tests/test_vec4.c 自带标量参考实现，
#            -ffp-contract=fast 下它自己的 a+b*c 会被收缩成 vfmadd132ss —— 那不是库的
#            指令，不能算进这条判据。）
#   配置 -O3 -march=native -mfma -ffp-contract=fast：最狠的一档。头文件契约声称即使
#   消费者的选项会把标量表达式收缩成 FMA，内核里也不许出现 vfmadd（靠 pragma 防）。
#   只有在这一档下普查，判据才有判别力。
BAD_RE='\b(vfmadd[0-9]*[sp][ds]|vfmsub[0-9]*[sp][ds]|vfnmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vhaddps|haddps|vmaskmovps)\b'
ALL_RE='\b(vfmadd[0-9]*[sp][ds]|vfmsub[0-9]*[sp][ds]|vfnmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vhaddps|haddps|vmaskmovps|vmovntps|vaddps|vmovups|vblendvps)\b'
HARD='-O3 -march=native -mfma -ffp-contract=fast'
for s in flashvec flashvec_batch flashvec_vec4; do
  c++ -std=c++17 $HARD -I"$INC" -c "$SRC/$s.cpp" -o "$TMP/lib_$s.o" 2>"$TMP/lib_$s.err" \
    || { fail "编译 $SRC/$s.cpp 失败："; sed -n '1,8p' "$TMP/lib_$s.err"; }
done
cc -std=c11 $HARD -I"$INC" -c ../tests/test_stream.c -o "$TMP/inline_stream.o" 2>"$TMP/inline_stream.err" \
  || { fail "编译 tests/test_stream.c 失败："; sed -n '1,8p' "$TMP/inline_stream.err"; }
TOTAL_BAD=0
TOTAL_NT=0
for obj in "$TMP/lib_flashvec.o" "$TMP/lib_flashvec_batch.o" "$TMP/lib_flashvec_vec4.o" "$TMP/inline_stream.o"; do
  if [ ! -f "$obj" ]; then continue; fi
  bad=$(objdump -d --no-show-raw-insn "$obj" | grep -oE "$BAD_RE" | wc -l)
  nt=$(objdump -d --no-show-raw-insn "$obj" | grep -oE 'vmovntps' | wc -l)
  sig=$(objdump -d --no-show-raw-insn "$obj" | grep -oE "$ALL_RE" | sort | uniq -c | tr '\n' ' ')
  TOTAL_BAD=$((TOTAL_BAD + bad))
  TOTAL_NT=$((TOTAL_NT + nt))
  printf '  %-24s 禁用=%s vmovntps=%s | %s\n' "$(basename "$obj")" "$bad" "$nt" "$sig"
done
info "（上面只列这几类指令；未列出 = 该对象没有。haddps 也算禁用项：归约的杠杆在写法，不在 hadd）"
if [ "$TOTAL_BAD" -eq 0 ]; then
  pass "7 项禁用助记符合计 = 0（vfmadd/vfmsub/vfnmadd/vrcpps/vrsqrtps/haddps/vmaskmovps）"
else
  fail "禁用助记符合计 = $TOTAL_BAD，必须为 0"
fi
if [ "$TOTAL_NT" -gt 0 ]; then
  pass "vmovntps = $TOTAL_NT > 0：NT store 真的在用（_stream 族不是空壳）"
else
  fail "vmovntps = 0：_stream（NT store）入口没有生成 NT store 指令"
fi
info "vmovntps 的来源（逐符号，库内 extern 形态；逐符号之和 <= 对象总计 ="
info "  某些内核被 GCC 复制/内联到别的符号里，判据看上面的对象总计，不看单行）："
for sym in fv_vec2_add_x4_stream fv_vec2_sub_x4_stream fv_vec2_scale_x4_stream \
           fv_vec2_add_soa_stream fv_vec2_scale_soa_stream fv_vec2_normalize_soa_stream; do
  n=$(objdump -d --no-show-raw-insn --disassemble="$sym" "$TMP/lib_flashvec_batch.o" 2>/dev/null | grep -oE 'vmovntps' | wc -l)
  printf '    %-32s vmovntps=%s\n' "$sym" "$n"
done

# ---------------------------------------------------------------- [9]
rule 9 "不带 -march=native 的默认档：AVX2 内核仍然存在，且派发是运行时的（不是编译期展开成标量）"
# 旧定义：在 build/flashvec_batch_O2.o / _O3.o 里找文件内静态内核符号 fv_*_avx2 与
#   入口函数的 call 关系（当时内核在 src/flashvec_batch.cpp 里，是文件内静态函数）。
# 新定义：内核是头文件里的 static inline（带 target("avx2")），「内核还在」对每个消费者
#   TU 都成立，所以改成直接编一个**纯 C、无 -march=native** 的消费者来判：
#     (a) fv_i_*_avx2 局部符号存在（内核没被优化掉，也没被条件编译关掉）；
#     (b) 里面有 AVX2 指令（vaddps 等）——基线下能出向量指令只可能来自 target 属性；
#     (c) 派发是运行时的：对象引用 __cpu_model，并在 avx2 位 0x400 上做判断。
if cc -std=c11 -O2 -I"$INC" -c ../tests/test_batch_c.c -o "$TMP/plain_o2.o" 2>"$TMP/plain_o2.err"; then
  info "纯 C 消费者（cc -std=c11 -O2，无 -march=native）编译成功"
else
  fail "默认档编译失败："; sed -n '1,8p' "$TMP/plain_o2.err"
fi
KS=$(nm --defined-only "$TMP/plain_o2.o" 2>/dev/null | awk '$2=="t"{print $3}' | grep -cE 'avx2')
VEC=$(objdump -d --no-show-raw-insn "$TMP/plain_o2.o" 2>/dev/null | grep -oE '\b(vaddps|vsubps|vmulps|vdivps|vsqrtps|vblendvps|vmovups)\b' | sort | uniq -c | tr '\n' ' ')
CM=$(nm -u "$TMP/plain_o2.o" 2>/dev/null | grep -c '__cpu_model')
BIT=$(objdump -d --no-show-raw-insn "$TMP/plain_o2.o" 2>/dev/null | grep -cE 'and[ql]?[[:space:]]+\$0x400')
printf '  %-24s avx2内核局部符号=%s  __cpu_model引用=%s  avx2位(0x400)判断=%s | %s\n' \
  "plain_o2.o" "$KS" "$CM" "$BIT" "$VEC"
for s in flashvec_batch flashvec_vec4; do
  c++ -std=c++17 -O3 -fno-math-errno -ffp-contract=off -I"$INC" -c "$SRC/$s.cpp" -o "$TMP/def_$s.o" 2>/dev/null
done
KS2=$(nm --defined-only "$TMP/def_flashvec_batch.o" 2>/dev/null | awk '$2=="t"{print $3}' | grep -cE 'avx2')
VEC2=$(objdump -d --no-show-raw-insn "$TMP/def_flashvec_batch.o" 2>/dev/null | grep -oE '\b(vaddps|vsubps|vmulps|vdivps|vsqrtps|vmovntps|vblendvps)\b' | sort | uniq -c | tr '\n' ' ')
printf '  %-24s avx2内核局部符号=%s | %s\n' "def_flashvec_batch.o" "$KS2" "$VEC2"
info "默认档编译选项取根 Makefile 的 -O3 -fno-math-errno -ffp-contract=off（无 -march=native）"
if [ "$KS" -gt 0 ] && [ "$KS2" -gt 0 ] && [ "$CM" -gt 0 ] && [ "$BIT" -gt 0 ]; then
  pass "基线下 AVX2 内核存在（纯 C 消费者 $KS 个 / 库 TU $KS2 个），派发走 __cpu_model 的 avx2 位"
else
  fail "基线下内核缺失或派发不是运行时的（KS=$KS KS2=$KS2 cpu_model=$CM bit=$BIT）"
fi

# ---------------------------------------------------------------- [10]
rule 10 "默认档（-O3，不带 -march=native）实测加速比摘要（旧规则，原样保留）"
# 本脚本**不做性能测量**（测量属于 bench/ab 与 verify/perf 的活），只引已有产物。
if [ ! -f build/perf_O3_plain.txt ]; then
  info "build/perf_O3_plain.txt 尚未生成（本轮不重跑 perf：不做性能测量）"
elif [ "$(wc -c < build/perf_O3_plain.txt)" -lt 200 ]; then
  info "build/perf_O3_plain.txt 只有 $(wc -c < build/perf_O3_plain.txt) 字节 —— 是 10-17 那次"
  info "  Makefile 变量转义 bug 留下的过期垃圾（P 变量被 make 吃掉），不是测量结果，见 verify/README.md §7"
else
  sed -n '/=== N=8192/,/^$/p' build/perf_O3_plain.txt | head -14 | sed 's/^/  /'
fi

# ---------------------------------------------------------------- [11]
rule 11 "已废弃内部名不得复活：VMASKMOVPS 策略的两个内核名在整仓出现次数必须为 0"
# 这两个名字绑定的是「只写命中 lane」的 VMASKMOVPS 方案，在本机被实测否决（Zen 2 的
#   maskmov 对未命中的 lane 仍然发生访存，稀疏收益拿不到 → 改用整块读改写 + 显式写回
#   a[i] 的 _mask 语义，以及按密度二选一的 _mask_inplace）。权威说法只留头文件契约一处。
# 所以它一旦在任何文件里复活（旧分支被 rebase 回来、旧的 bench 夹具被复制回来），这条
#   规则要立刻响。模式末位写成 [2]，否则本脚本自己的规则文本会命中自己。
# 扫描范围 = 代码与构建文件（*.h *.c *.cpp *.hpp *.cc + Makefile）。
# 为什么不含 *.md 与本脚本自己：它们**讨论**这条禁令，必然会写出被禁的名字；把说明也算命中，
# 这条规则就只能靠删掉说明来通过 —— 那不是检查，是自欺。名字在代码里复活才是要抓的事。
OLDNAME='fv_i_maskload_add4_avx[2]|fv_i_maskstore_add4_avx[2]'
SCAN='..'
HITS=$(grep -rInE "$OLDNAME" "$SCAN" --exclude-dir=.git --exclude-dir=build \
        --include='*.h' --include='*.c' --include='*.cpp' --include='*.hpp' --include='*.cc' \
        --include='Makefile' 2>/dev/null | wc -l)
info "扫描范围：仓库内所有 .h/.c/.cpp/.hpp/.cc 与 Makefile（跳过 .git 与 build/；-I 跳过二进制）"
if [ "$HITS" -eq 0 ]; then
  pass "已废弃内部名在代码里出现 0 次"
else
  fail "已废弃内部名在代码里出现 $HITS 次："
  grep -rInE "$OLDNAME" "$SCAN" --exclude-dir=.git --exclude-dir=build \
       --include='*.h' --include='*.c' --include='*.cpp' --include='*.hpp' --include='*.cc' \
       --include='Makefile' 2>/dev/null | sed 's/^/    /'
fi

# ---------------------------------------------------------------- 汇总
echo
printf '== 汇总：%d 条规则，FAIL=%d ==\n' "$RULES" "$FAILS"
if [ "$FAILS" -eq 0 ]; then
  echo "== 全部规则 PASS =="
else
  echo "== 有 $FAILS 条规则未通过（逐条见上）=="
fi
echo "$FAILS" > build/.checks_fails
} 2>&1 | tee build/checks.txt

RC=$(cat build/.checks_fails 2>/dev/null || echo 1)
rm -f build/.checks_fails
[ "$RC" = 0 ] && exit 0 || exit 1
