#!/bin/sh
# verify/checks.sh -- structural checks (reviewer F): turn the spec items that "running tests cannot prove" into reproducible evidence.
# Reads the tree under test only; the only writes are the intermediate objects under verify/build/ and this report.
# The report is also written to build/checks.txt.
# ============================================================================
# This script was rewritten after the "batch entry points move into the header" refactor. Each rule is
# preceded by "old definition -> new definition" and by why it had to change. Three facts came out of
# the refactor, and every rule is built around them:
#   1) the **implementation** of the batch entry points (AVX2 kernels, target("avx2"), the
#      __builtin_cpu_supports("avx2") != 0 dispatch) now lives in include/flashvec.h as static inline
#      (the default form); src/*.cpp keeps only the extern fallback symbols (FV_USE_EXTERN / FV_BUILD
#      forms). -> the old rules that scanned src/*.cpp no longer see it.
#   2) the header may therefore contain conditional intrinsics (it includes <immintrin.h> when
#      FV_I_INTRINSICS && !FV_USE_EXTERN). -> "the header must not contain immintrin" no longer holds,
#      but the stronger statement "the FV_USE_EXTERN include set matches the pre-refactor one" holds
#      and is more decidable.
#   3) the old rules disassembled a **stale build/*.o** -- objects compiled before the refactor, so that
#      census accepts an artifact that does not exist. -> the script recompiles both forms itself first.
# ============================================================================
# Verdict: 11 rules; [7] prints machine facts only (no criterion); all others PASS / FAIL; any FAIL exits non-zero.
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
echo "== structural checks (read-only)  $(date -Iseconds) =="
info "tree under test: $(cd .. && pwd)"
info "rules: [1][2] header purity / [3] compiler flags / [4] dispatch location / [5] banned builtins in source /"
info "  [6] header contract text / [7] machine facts / [8] artifact-level instruction census / [9] default-tier kernels + dispatch /"
info "  [10] default-tier measured summary / [11] retired internal names must not come back"

# ---------------------------------------------------------------- [1]
rule 1 "header purity A: the default form (implementation in the header) works under a plain C compiler, and intrinsics really expand"
# Old definition: include/flashvec.h must contain none of immintrin/_mm_/__m256/namespace/std:: (the header was pure declarations then).
# New definition: now that the implementation lives in the header, the header **must** include
#   <immintrin.h> when FV_I_INTRINSICS && !FV_USE_EXTERN. So the first half of the old criterion flips into a positive control:
#     * the number of expanded immintrin lines must be > 0 -- 0 would mean the inline form is not in effect (the wrapping is broken);
#     * "usable from plain C" is carried by stronger evidence: cc -std=c11 -pedantic-errors compiles it.
if cc -std=c11 -pedantic-errors -Wall -Wextra -I"$INC" -c ../tests/test_header_c.c -o "$TMP/hdr_c.o" 2>"$TMP/hdr_c.err"; then
  pass "cc -std=c11 -pedantic-errors compiles a plain C consumer (the header is usable from C)"
else
  fail "plain C compilation failed, first lines:"
  sed -n '1,12p' "$TMP/hdr_c.err"
fi
EXPAND_C=$(cc -std=c11 -I"$INC" -E ../tests/test_header_c.c 2>/dev/null | grep -c immintrin)
EXPAND_PP=$(c++ -std=c++17 -I"$INC" -E ../tests/test_header_cpp.cpp 2>/dev/null | grep -c immintrin)
info "immintrin lines under cc -std=c11 -E = $EXPAND_C (>0 = conditional intrinsics are in effect, as the refactor expects)"
info "immintrin lines under c++ -std=c++17 -E = $EXPAND_PP"
if [ "$EXPAND_C" -gt 0 ] && [ "$EXPAND_PP" -gt 0 ]; then
  pass "the default form really expands intrinsics (positive control; the old rule required 0 here, which no longer applies)"
else
  fail "the immintrin expansion count is not positive: the header's inline implementation is not in effect"
fi
if grep -n -E 'namespace|template *<|std::|reinterpret_cast|static_cast|class [A-Za-z]' "$INC/flashvec.h"; then
  fail "C++ syntax in the header (the hits are above)"
else
  pass "the header has no C++ syntax (namespace/template</std::/cast/class)"
fi

# ---------------------------------------------------------------- [2]
rule 2 "header purity B: the include set of the FV_USE_EXTERN form must match the pre-refactor one"
# Old definition: the header must contain no target attribute / __builtin / extern implementation (all declarations only).
# New definition: the FV_USE_EXTERN form **is** that "declarations only" form, so the old rule is restated decidably:
#   the headers the preprocessor pulls in must contain no math.h / string.h / immintrin.h / stdlib.h.
# Why the include set instead of grep: the intrinsics in the header text are conditionally compiled, so only listing
#   what the preprocessor actually pulls in can prove "the fallback form's include set is word-for-word the pre-refactor one".
DEPS=$(printf '#include "flashvec.h"\n' | cc -std=c11 -DFV_USE_EXTERN -I"$INC" -M -x c - 2>/dev/null)
HDRS=$(printf '%s\n' "$DEPS" | tr '[:space:]' '\n' | sed 's/\\$//' | grep -oE '[^/]+\.h$' | sort -u)
info "FV_USE_EXTERN include set: $(printf '%s' "$HDRS" | tr '\n' ' ')"
HIT=$(printf '%s\n' "$HDRS" | grep -xE 'math\.h|string\.h|immintrin\.h|stdlib\.h' | tr '\n' ' ')
if [ -n "$HIT" ]; then
  fail "the FV_USE_EXTERN form pulls in extra headers: $HIT (the pre-refactor include set had none of them)"
else
  pass "the FV_USE_EXTERN form contains no math.h / string.h / immintrin.h / stdlib.h"
fi
if printf '#include "flashvec.h"\n' | cc -std=c11 -DFV_USE_EXTERN -pedantic-errors -I"$INC" -c -x c -o "$TMP/ext.o" - 2>"$TMP/ext.err"; then
  pass "the FV_USE_EXTERN form compiles as plain C under -pedantic-errors"
else
  fail "the FV_USE_EXTERN form failed to compile, first lines:"
  sed -n '1,12p' "$TMP/ext.err"
fi

# ---------------------------------------------------------------- [3]
rule 3 "the default build must not use -march=native / -mavx2 / -ffast-math (the criterion sits on the **actual compile commands**)"
# Old definition: any non-comment line of the root Makefile containing -march=native/-mavx/-ffast-math/-Ofast fails.
# New definition: the **intent** of this rule is "the default artifacts keep the x86-64 baseline", but that grep is a
#   gross approximation -- it also condemns FMA_FLAGS, which is used only by tests/test_parity.c and tests/test_vec34_parity.c
#   through their *_cpp_fma.o (the **hostile-flags comparison** fixture, a consumer built with -march=native -mfma -ffp-contract=fast),
#   which is exactly the test group that proves the default tier is safe without those flags. Failing it treats a fixture as a violation.
#   So the criterion moves from "scan the text" to "ask what the default make target would actually execute": every compile command
#   for src/ printed by make -n -B all must be clean. That is stricter than grep -- it contains the real, expanded result.
if command -v make >/dev/null 2>&1; then
  # Note the -C .. : this script's working directory is verify/, while the -march=native criterion is about the root Makefile
  # (verify/Makefile has a -O3 -march=native tier of its own; that is for the comparison runs, not the default artifacts).
  DEFAULT_CMDS=$(make -C .. -n -B all 2>/dev/null | grep -E -- '-c src/')
  if [ -z "$DEFAULT_CMDS" ]; then
    info "make -C .. -n -B all printed no command compiling src/ (the build shape may have changed); falling back to the variable lines only"
  else
    printf '%s\n' "$DEFAULT_CMDS" | sed 's/^/  /'
    if printf '%s\n' "$DEFAULT_CMDS" | grep -E 'march=native|mavx|ffast-math|Ofast|mfma'; then
      fail "aggressive compiler flags in the default build commands (the hits are above)"
    else
      pass "every compile command of the root Makefile's default target is x86-64 baseline (no -march=native/-mavx/-ffast-math/-mfma)"
    fi
  fi
else
  info "no make available, skipping the default build command check"
fi
info "CXXFLAGS/CFLAGS fallback lines:"
grep -n -E '^(CXXFLAGS|CFLAGS)[[:space:]]*[?:]?=' ../Makefile | sed 's/^/  /'
if grep -n -E '^(CXXFLAGS|CFLAGS)[[:space:]]*[?:]?=.*(march=native|mavx|ffast-math|Ofast)' ../Makefile; then
  fail "the CXXFLAGS/CFLAGS defaults contain aggressive compiler flags"
else
  pass "the CXXFLAGS/CFLAGS defaults are clean"
fi
info "every occurrence of -march=native in the Makefile (expected only in the hostile-flag comparison fixtures):"
grep -n 'march=native' ../Makefile | sed 's/^/  /'
# Every user of FMA_FLAGS (the rule head sitting above its recipe) must be a comparison fixture, *_fma.
# awk takes "the previous line" rather than grep taking neighbouring lines: the rule head and the recipe are two lines, and grep cannot tell which owns the flag.
FMA_HEADS=$(awk '/FMA_FLAGS\)/ {print prev} {prev=$0}' ../Makefile)
info "rule heads referencing FMA_FLAGS: $(printf '%s' "$FMA_HEADS" | tr '\n' ' ')"
FMA_BAD=$(printf '%s\n' "$FMA_HEADS" | grep -c -v '_fma')
if [ "$FMA_BAD" -eq 0 ] && [ -n "$FMA_HEADS" ]; then
  pass "-march=native appears only in the hostile-flag comparison fixtures (*_cpp_fma) and never enters a default artifact"
else
  fail "FMA_FLAGS has a non-comparison use ($FMA_BAD rule heads do not contain _fma): a default artifact may be contaminated"
fi

# ---------------------------------------------------------------- [4]
rule 4 "dispatch now lives in include/flashvec.h (target attribute + run-time __builtin_cpu_supports)"
# Old definition: scan every .cpp under src/ for the target attribute and __builtin_cpu_supports.
# New definition: once the implementation (kernels and dispatch) moved into the header, neither appears in src/*.cpp;
#   the criterion moves to include/flashvec.h and requires every site to be written as != 0 -- the return value is a bit mask
#   (avx2 = 0x400), and one uniform spelling is what keeps anyone from using it as a boolean.
CS_ALL=$(grep -c '__builtin_cpu_supports' "$INC/flashvec.h")
CS_NE0=$(grep -c '__builtin_cpu_supports("avx2") != 0' "$INC/flashvec.h")
TGT=$(grep -c 'target(' "$INC/flashvec.h")
SRC_CS=$(cat "$SRC"/*.cpp 2>/dev/null | grep -c '__builtin_cpu_supports')
info "include/flashvec.h: __builtin_cpu_supports = $CS_ALL (written as != 0: $CS_NE0); target( = $TGT"
info "src/*.cpp: __builtin_cpu_supports = $SRC_CS (0 = dispatch really did move into the header, an architectural fact after the refactor)"
if [ "$CS_ALL" -gt 0 ] && [ "$CS_NE0" -eq "$CS_ALL" ] && [ "$TGT" -gt 0 ]; then
  pass "the header has the target attribute and $CS_ALL run-time dispatch sites, all written as != 0"
else
  fail "the header's dispatch/attribute is incomplete or not uniformly written (total=$CS_ALL ne0=$CS_NE0 target=$TGT)"
fi
info "dispatch sites (first 6, line numbers from include/flashvec.h):"
grep -n '__builtin_cpu_supports' "$INC/flashvec.h" | head -6 | sed 's/^/    /'

# ---------------------------------------------------------------- [5]
rule 5 "the bit-exact paths must not use FMA / reciprocal approximations / VMASKMOVPS (source level, header included)"
# Old definition: scan every .cpp under src/.
# New definition: the implementation now mostly lives in include/flashvec.h, so scan include/ + src/.
#   What is scanned are **builtin names** (_mm*_fmadd_ps and the like), not assembly mnemonics: the header's contract
#   comments legitimately say things like 'the vfmadd count is 0' or 'no rcpps/rsqrtps', and that is documentation, not a failure.
#   Among them _mm*_maskload/maskstore (= VMASKMOVPS) is the sparse strategy that measurement rejected: on Zen 2 it still
#   touches memory for lanes that missed, so the gain from writing only the hit lanes never materialises -- it is not even allowed.
BAN='_mm(128|256)?_(fmadd|fmsub|fnmadd|fmaddsub|fmsubadd|rcp|rsqrt|maskload|maskstore)'
if grep -rnE "$BAN" ../include ../src; then
  fail "the banned builtin names are hit above"
else
  pass "banned builtin names (FMA / rcpps / rsqrtps / maskload / maskstore) in include/ + src/ = 0"
fi
DOC=$(grep -cE '\b(vfmadd|vrcpps|vrsqrtps|haddps|vmaskmovps)\b' "$INC/flashvec.h")
info "mnemonic lines appearing as **comment text** in the header = $DOC (allowed: a contract saying 'the count is 0' must be able to name them)"

# ---------------------------------------------------------------- [6]
rule 6 "the header writes out the dot_soa reduction-order disclaimer and the mask semantics (hard spec requirement)"
# Old definition: three keyword counts > 0. One stray mention of each word was enough, so the rule could not
#   tell a contract from a footnote; and it counted Chinese keywords that an English header no longer has.
# New definition: a mention is not a contract. Seven assertions, each attached to the entry it governs:
#   (a) the disclaimer sentence itself exists ("depends on the reduction order") and sits inside the
#       "Not guaranteed" statement, not merely somewhere in the file;
#   (b) the extern declaration of fv_vec2_dot_soa carries its own note too ("reduction order" within 5
#       lines above it), so the exception is stated where the API is read, not only in the contract above;
#   (c) fv_vec3_dot_soa is named where the same exception is stated;
#   (d) fv_vec4_dot_soa is named where the same exception is stated;
#   (e) the mask contract has a heading of its own ("Mask semantics");
#   (f) the mask rule is spelled out under that heading: a lane outside the mask gets a[i]'s original
#       value *written back*, which is the part callers misread as "not written";
#   (g) the deviation is never promoted into the promise: no line of the "Guaranteed - ..." list may name
#       a *_dot_soa entry.
# Thresholds are deliberately not the criterion. "mask" occurs about 85 times in the header, so any
#   count-based threshold would be met by unrelated sentences (mask[i] != 0 inside a kernel body, _mask in
#   a signature). (f) is what makes the wording load-bearing, and (a)-(d) are what make deleting one
#   explanation fail this rule instead of passing on a leftover mention elsewhere.
# Evidence first, verdict second: the keyword lines are printed before the assertions run.
KW=$(grep -n -E 'reduction order|[Mm]ask semantics|lane|dot_soa' "$INC/flashvec.h")
printf '%s\n' "$KW" | sed -n '1,20p' | sed 's/^/  /'
info "(only the first 20 lines are printed above; $(printf '%s\n' "$KW" | wc -l) lines in total)"
CONTRACT=$(awk '
function near(s, ln, last, win,   j) {
    for (j = ln - win; j <= ln + win; j++)
        if (j >= 1 && j <= last && index(L[j], s) > 0) return 1
    return 0
}
function above(s, ln, win,   j) {
    for (j = ln - win; j < ln; j++)
        if (j >= 1 && index(L[j], s) > 0) return 1
    return 0
}
{ L[NR] = $0 }
END {
    a = 0; b = 0; c = 0; d = 0; e = 0; f = 0; g = 0
    inpromise = 0
    for (i = 1; i <= NR; i++) {
        if (index(L[i], "Guaranteed - ") > 0) inpromise = 1
        if (index(L[i], "Not guaranteed - ") > 0) inpromise = 0
        if (inpromise && index(L[i], "dot_soa") > 0) g = 1
        if (index(L[i], "depends on the reduction order") > 0 && near("Not guaranteed", i, NR, 5)) a = 1
        if (index(L[i], "float fv_vec2_dot_soa(") == 1 && above("reduction order", i, 5)) b = 1
        if (index(L[i], "fv_vec3_dot_soa") > 0)
            if (near("exception", i, NR, 5) || near("reduction", i, NR, 5) ||
                near("not promised", i, NR, 5) || near("not guaranteed", i, NR, 5)) c = 1
        if (index(L[i], "fv_vec4_dot_soa") > 0)
            if (near("exception", i, NR, 5) || near("reduction", i, NR, 5) ||
                near("not promised", i, NR, 5) || near("not guaranteed", i, NR, 5)) d = 1
        if (index(L[i], "Mask semantics") > 0) {
            e = 1
            for (j = i; j <= i + 10 && j <= NR; j++)
                if (index(L[j], "mask") > 0 && index(L[j], "written back") > 0 &&
                    index(L[j], "a[i]") > 0) f = 1
        }
    }
    printf "%d %d %d %d %d %d %d\n", a, b, c, d, e, f, g
}' "$INC/flashvec.h")
[ -n "$CONTRACT" ] || CONTRACT="0 0 0 0 0 0 0"
set -- $CONTRACT
info "assertions a=disclaimer b=attached-to-fv_vec2_dot_soa c=vec3-sibling d=vec4-sibling e=mask-heading f=mask-write-back g=promoted-must-be-0 -> $1 $2 $3 $4 $5 $6 $7"
if [ "$1" = 1 ] && [ "$2" = 1 ] && [ "$3" = 1 ] && [ "$4" = 1 ] && [ "$5" = 1 ] && [ "$6" = 1 ] && [ "$7" = 0 ]; then
  pass "both contract points are written out and attached to the entries they govern (dot_soa reduction order, mask write-back)"
else
  fail "a contract point is missing or detached from its entry: a=$1 b=$2 c=$3 d=$4 e=$5 f=$6 g=$7 (g must be 0)"
fi


# ---------------------------------------------------------------- [7]
rule 7 "machine / compiler facts (printed only, no criterion)"
grep -m1 'model name' /proc/cpuinfo | sed 's/^/  /'
printf '  avx2/fma: '; grep -o -m1 -E 'avx2|fma' /proc/cpuinfo | sort -u | tr '\n' ' '; echo
lscpu 2>/dev/null | grep -E '^L[123]|^CPU\(s\)' | sed 's/^/  /'
# The cache hierarchy must be labelled from measurement: lscpu's "8 MiB (2 instances)" is the whole chip, while a single thread only gets
# its own CCX's 4 MiB. The project docs once called this "8 MiB of shared L3" and computed the knee point wrongly because of it.
for i in 2 3; do
  if [ -r "/sys/devices/system/cpu/cpu0/cache/index$i/size" ]; then
    printf '  index%s: level=%s size=%s shared_cpu_list=%s\n' "$i" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/level)" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/size)" \
      "$(cat /sys/devices/system/cpu/cpu0/cache/index$i/shared_cpu_list)"
  fi
done
info "convention: L3 4 MiB/CCX x 2 CCX = 8 MiB per chip; a single thread sees 4 MiB (the knee on this machine is about 175K elements @24 B/element)"
c++ --version | head -1 | sed 's/^/  /'
cc --version | head -1 | sed 's/^/  /'

# ---------------------------------------------------------------- [8]
rule 8 "artifact-level instruction census: the 7 banned mnemonics must be 0 and vmovntps must be > 0 (the NT store positive control)"
# Old definition: census per symbol over build/flashvec_O3n.o / build/flashvec_batch_O3n.o.
# New definition: the kernels live in the header, while those objects under build/ were compiled before the refactor (censusing them means accepting an artifact that does not exist).
#   The script recompiles both forms itself under the **strictest configuration**:
#     form A, in-library extern: src/flashvec.cpp + flashvec_batch.cpp + flashvec_vec4.cpp
#            (the extern definitions of vec3's SoA and of vec4's _x2/_soa are all in flashvec_batch.cpp)
#     form B, a header-inline consumer: tests/test_stream.c (plain C; its hot loops only call the library and carry no arithmetic
#            of their own, so censusing the whole TU is legitimate. Counterexample: tests/test_vec4.c carries its own scalar reference
#            implementation, and under -ffp-contract=fast its own a+b*c contracts into vfmadd132ss -- those are not the library's
#            instructions and cannot count towards this criterion.)
#   Configuration -O3 -march=native -mfma -ffp-contract=fast: the harshest tier. The header contract claims that even when
#   the consumer's options would contract scalar expressions into FMA, no vfmadd may appear inside a kernel (the pragma prevents it).
#   Only a census at this tier gives the criterion any discriminating power.
BAD_RE='\b(vfmadd[0-9]*[sp][ds]|vfmsub[0-9]*[sp][ds]|vfnmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vhaddps|haddps|vmaskmovps)\b'
ALL_RE='\b(vfmadd[0-9]*[sp][ds]|vfmsub[0-9]*[sp][ds]|vfnmadd[0-9]*[sp][ds]|vrcpps|vrsqrtps|vhaddps|haddps|vmaskmovps|vmovntps|vaddps|vmovups|vblendvps)\b'
HARD='-O3 -march=native -mfma -ffp-contract=fast'
for s in flashvec flashvec_batch flashvec_vec4; do
  c++ -std=c++17 $HARD -I"$INC" -c "$SRC/$s.cpp" -o "$TMP/lib_$s.o" 2>"$TMP/lib_$s.err" \
    || { fail "compiling $SRC/$s.cpp failed:"; sed -n '1,8p' "$TMP/lib_$s.err"; }
done
cc -std=c11 $HARD -I"$INC" -c ../tests/test_stream.c -o "$TMP/inline_stream.o" 2>"$TMP/inline_stream.err" \
  || { fail "compiling tests/test_stream.c failed:"; sed -n '1,8p' "$TMP/inline_stream.err"; }
TOTAL_BAD=0
TOTAL_NT=0
for obj in "$TMP/lib_flashvec.o" "$TMP/lib_flashvec_batch.o" "$TMP/lib_flashvec_vec4.o" "$TMP/inline_stream.o"; do
  if [ ! -f "$obj" ]; then continue; fi
  bad=$(objdump -d --no-show-raw-insn "$obj" | grep -oE "$BAD_RE" | wc -l)
  nt=$(objdump -d --no-show-raw-insn "$obj" | grep -oE 'vmovntps' | wc -l)
  sig=$(objdump -d --no-show-raw-insn "$obj" | grep -oE "$ALL_RE" | sort | uniq -c | tr '\n' ' ')
  TOTAL_BAD=$((TOTAL_BAD + bad))
  TOTAL_NT=$((TOTAL_NT + nt))
  printf '  %-24s banned=%s vmovntps=%s | %s\n' "$(basename "$obj")" "$bad" "$nt" "$sig"
done
info "(only these instruction classes are listed; not listed = absent from that object. haddps counts as banned too: the lever for a reduction is how it is written, not hadd)"
if [ "$TOTAL_BAD" -eq 0 ]; then
  pass "the 7 banned mnemonics add up to 0 (vfmadd/vfmsub/vfnmadd/vrcpps/vrsqrtps/haddps/vmaskmovps)"
else
  fail "banned mnemonics add up to $TOTAL_BAD, must be 0"
fi
if [ "$TOTAL_NT" -gt 0 ]; then
  pass "vmovntps = $TOTAL_NT > 0: NT stores are really used (the _stream family is not an empty shell)"
else
  fail "vmovntps = 0: the _stream (NT store) entry points emit no NT store instruction"
fi
info "where vmovntps comes from (per symbol, in-library extern form; the per-symbol sum <= the object total ="
info "  some kernels are copied or inlined by GCC into other symbols, so the criterion reads the object total above, not a single line):"
for sym in fv_vec2_add_x4_stream fv_vec2_sub_x4_stream fv_vec2_scale_x4_stream \
           fv_vec2_add_soa_stream fv_vec2_scale_soa_stream fv_vec2_normalize_soa_stream; do
  n=$(objdump -d --no-show-raw-insn --disassemble="$sym" "$TMP/lib_flashvec_batch.o" 2>/dev/null | grep -oE 'vmovntps' | wc -l)
  printf '    %-32s vmovntps=%s\n' "$sym" "$n"
done

# ---------------------------------------------------------------- [9]
rule 9 "the default tier without -march=native: the AVX2 kernels still exist and dispatch is run-time (not expanded to scalar at compile time)"
# Old definition: find the file-local kernel symbols fv_*_avx2 and the call relation to the entry functions
#   in build/flashvec_batch_O2.o / _O3.o (the kernels were file-local static functions in src/flashvec_batch.cpp then).
# New definition: the kernels are static inline in the header (carrying target("avx2")), so "the kernel is still there" is true for every
#   consumer TU; the criterion is therefore a consumer compiled as **plain C with no -march=native**:
#     (a) the fv_i_*_avx2 local symbols exist (the kernel was neither optimised away nor compiled out);
#     (b) they contain AVX2 instructions (vaddps and the like) -- at baseline only the target attribute can produce vector instructions;
#     (c) dispatch is run-time: the object references __cpu_model and tests the avx2 bit 0x400.
if cc -std=c11 -O2 -I"$INC" -c ../tests/test_batch_c.c -o "$TMP/plain_o2.o" 2>"$TMP/plain_o2.err"; then
  info "the plain C consumer (cc -std=c11 -O2, no -march=native) compiled successfully"
else
  fail "the default tier failed to compile:"; sed -n '1,8p' "$TMP/plain_o2.err"
fi
KS=$(nm --defined-only "$TMP/plain_o2.o" 2>/dev/null | awk '$2=="t"{print $3}' | grep -cE 'avx2')
VEC=$(objdump -d --no-show-raw-insn "$TMP/plain_o2.o" 2>/dev/null | grep -oE '\b(vaddps|vsubps|vmulps|vdivps|vsqrtps|vblendvps|vmovups)\b' | sort | uniq -c | tr '\n' ' ')
CM=$(nm -u "$TMP/plain_o2.o" 2>/dev/null | grep -c '__cpu_model')
BIT=$(objdump -d --no-show-raw-insn "$TMP/plain_o2.o" 2>/dev/null | grep -cE 'and[ql]?[[:space:]]+\$0x400')
printf '  %-24s avx2 kernel local symbols=%s  __cpu_model refs=%s  avx2 bit (0x400) test=%s | %s\n' \
  "plain_o2.o" "$KS" "$CM" "$BIT" "$VEC"
for s in flashvec_batch flashvec_vec4; do
  c++ -std=c++17 -O3 -fno-math-errno -ffp-contract=off -I"$INC" -c "$SRC/$s.cpp" -o "$TMP/def_$s.o" 2>/dev/null
done
KS2=$(nm --defined-only "$TMP/def_flashvec_batch.o" 2>/dev/null | awk '$2=="t"{print $3}' | grep -cE 'avx2')
VEC2=$(objdump -d --no-show-raw-insn "$TMP/def_flashvec_batch.o" 2>/dev/null | grep -oE '\b(vaddps|vsubps|vmulps|vdivps|vsqrtps|vmovntps|vblendvps)\b' | sort | uniq -c | tr '\n' ' ')
printf '  %-24s avx2 kernel local symbols=%s | %s\n' "def_flashvec_batch.o" "$KS2" "$VEC2"
info "the default tier uses the root Makefile's -O3 -fno-math-errno -ffp-contract=off (no -march=native)"
if [ "$KS" -gt 0 ] && [ "$KS2" -gt 0 ] && [ "$CM" -gt 0 ] && [ "$BIT" -gt 0 ]; then
  pass "AVX2 kernels exist at baseline ($KS in the plain C consumer / $KS2 in the library TU) and dispatch goes through __cpu_model's avx2 bit"
else
  fail "at baseline the kernels are missing or dispatch is not run-time (KS=$KS KS2=$KS2 cpu_model=$CM bit=$BIT)"
fi

# ---------------------------------------------------------------- [10]
rule 10 "default tier (-O3, no -march=native) measured speedup summary (old rule, kept as is)"
# This script **does no performance measurement** (that is the job of bench/ab and verify/perf); it only cites existing artifacts.
if [ ! -f build/perf_O3_plain.txt ]; then
  info "build/perf_O3_plain.txt does not exist yet (no perf rerun this round: no performance measurement)"
elif [ "$(wc -c < build/perf_O3_plain.txt)" -lt 200 ]; then
  info "build/perf_O3_plain.txt holds only $(wc -c < build/perf_O3_plain.txt) bytes -- the stale leftovers of that 10-17"
  info "  Makefile variable-escaping bug (make swallowed the P variable), not a measurement; see verify/README.md §7"
else
  sed -n '/=== N=8192/,/^$/p' build/perf_O3_plain.txt | head -14 | sed 's/^/  /'
fi

# ---------------------------------------------------------------- [11]
rule 11 "retired internal names must not come back: the two kernel names of the VMASKMOVPS strategy must appear 0 times in the whole repo"
# Those two names belong to the "write only the hit lanes" VMASKMOVPS scheme, which measurement rejected on this machine (on Zen 2 a
#   maskmov still touches memory for lanes that missed, so the sparse gain never materialises -> hence the whole-block read-modify-write
#   with an explicit write-back of a[i] for _mask, and the density-based two-way choice of _mask_inplace). The header contract is the one authoritative statement.
# So if it comes back in any file (an old branch rebased in, an old bench fixture copied back), this
#   rule must fire at once. The last digit of the pattern is written as [2], or this script's own rule text would match itself.
# Scan scope = code and build files (*.h *.c *.cpp *.hpp *.cc + Makefile).
# Why *.md files and this script are excluded: they **discuss** the ban and must write the banned names; if explanations counted as hits,
#   the rule could only be passed by deleting the explanation -- that is not a check, it is self-deception. A name coming back in code is what matters.
OLDNAME='fv_i_maskload_add4_avx[2]|fv_i_maskstore_add4_avx[2]'
SCAN='..'
HITS=$(grep -rInE "$OLDNAME" "$SCAN" --exclude-dir=.git --exclude-dir=build \
        --include='*.h' --include='*.c' --include='*.cpp' --include='*.hpp' --include='*.cc' \
        --include='Makefile' 2>/dev/null | wc -l)
info "scan scope: every .h/.c/.cpp/.hpp/.cc and Makefile in the repo (skipping .git and build/; -I skips binaries)"
if [ "$HITS" -eq 0 ]; then
  pass "retired internal names appear 0 times in code"
else
  fail "retired internal names appear $HITS times in code:"
  grep -rInE "$OLDNAME" "$SCAN" --exclude-dir=.git --exclude-dir=build \
       --include='*.h' --include='*.c' --include='*.cpp' --include='*.hpp' --include='*.cc' \
       --include='Makefile' 2>/dev/null | sed 's/^/    /'
fi

# ---------------------------------------------------------------- summary
echo
printf '== summary: %d rules, FAIL=%d ==\n' "$RULES" "$FAILS"
if [ "$FAILS" -eq 0 ]; then
  echo "== all rules PASS =="
else
  echo "== $FAILS rule(s) did not pass (see each one above) =="
fi
echo "$FAILS" > build/.checks_fails
} 2>&1 | tee build/checks.txt

RC=$(cat build/.checks_fails 2>/dev/null || echo 1)
rm -f build/.checks_fails
[ "$RC" = 0 ] && exit 0 || exit 1
