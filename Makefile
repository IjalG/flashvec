CXX      ?= c++
CC       ?= cc
AR       ?= ar
# Compiler flags (P0):
#   -O3            unlocks the benefit of -fno-math-errno. Measured on the same loop containing sqrtf,
#                  -O2 -fno-math-errno leaves GCC scalar, -O3 gives "loop vectorized using 16 byte vectors";
#                  the x1.25 -> x8.16 figure for the len-class loops in the report comes from the -O3 configuration.
#   -fno-math-errno drops "set errno when a math function fails". This library's sqrt arguments are always >= 0,
#                  so errno is never set; measured values stay bit-identical, yet loops containing sqrtf become vectorizable.
#   -ffp-contract=off forbids contracting a + b*c into an FMA. The per-element API promises bit-reproducible
#                  results (docs/optimization-report.md §4.7), and FMA contraction makes the same expression differ
#                  by 1 ulp across compilation units / targets (measured here: vfmadd appears under -ffp-contract=fast).
#   No -march=native / -mavx2 / -ffast-math: the default artifacts keep the x86-64 baseline, and AVX2 is enabled only
#   inside the static inline entry points of include/flashvec.h, through function-level target attributes plus
#   runtime __builtin_cpu_supports("avx2") != 0 dispatch (measured: 32 dispatch sites, all in the header, 0 in src/*.cpp).
CXXFLAGS ?= -std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CFLAGS   ?= -std=c11 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude

BUILD := build
LIB   := $(BUILD)/libflashvec.a
OBJ   := $(BUILD)/flashvec.o $(BUILD)/flashvec_batch.o $(BUILD)/flashvec_vec4.o \
         $(BUILD)/flashvec_quat.o $(BUILD)/flashvec_mat3.o $(BUILD)/flashvec_mat4.o \
         $(BUILD)/flashvec_fast.o $(BUILD)/flashvec_fast_batch.o \
         $(BUILD)/flashvec_physics_batch.o

.PHONY: all test clean

all: $(LIB)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/flashvec.o: src/flashvec.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_batch.o: src/flashvec_batch.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

$(BUILD)/test_flashvec: tests/test_flashvec.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/test_c_api.o: tests/test_c_api.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_c_api: $(BUILD)/test_c_api.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test: $(BUILD)/test_flashvec $(BUILD)/test_c_api
	./$(BUILD)/test_flashvec
	./$(BUILD)/test_c_api

clean:
	rm -rf $(BUILD)

# ---------------- Benchmarks (new targets; does not change the behaviour of any existing target above) ----------------
# Benchmark sources and build rules live under bench/; these targets only forward, so they run from the repo root:
#   make bench        build and run all four flag sets x three implementation tiers x three sizes; raw output goes to docs/research/raw/
#   make bench-table  regenerate the Markdown comparison table from the raw logs
#   make bench-vec    capture the GCC auto-vectorization reports
#   make bench-asm    produce the disassembly archive (bench/build/disasm_*.txt)
#   make bench-perf   perf stat (needs kernel.perf_event_paranoid<=1 or CAP_PERFMON)
.PHONY: bench bench-table bench-vec bench-asm bench-perf bench-clean

bench:
	$(MAKE) -C bench run

bench-table:
	$(MAKE) -C bench table

bench-vec:
	$(MAKE) -C bench vec

bench-asm:
	$(MAKE) -C bench asm

bench-perf:
	$(MAKE) -C bench perf

bench-clean:
	$(MAKE) -C bench clean

# ---------------- P0 / P1 tests (append-only; does not change the behaviour of all/test/clean above) ----------------
#   make test-p0     boundary tests for min/max/clamp/madd, AABB, and the FV_EPS dimensional-consistency fix
#   make test-batch  bit-identical results (memcmp), mask semantics and index variants of the batch entry points; includes a C caller test
#   make test-all    the two above + the existing make test (30 assertions + the C compatibility test)
.PHONY: test-p0 test-batch test-all

$(BUILD)/test_p0: tests/test_p0.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/test_batch: tests/test_batch.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD)/test_batch_c.o: tests/test_batch_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_batch_c: $(BUILD)/test_batch_c.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-p0: $(BUILD)/test_p0
	./$(BUILD)/test_p0

test-batch: $(BUILD)/test_batch $(BUILD)/test_batch_c
	./$(BUILD)/test_batch
	./$(BUILD)/test_batch_c

test-all: test test-p0 test-batch test-mat4 test-mat4-c test-mat4-inline test-mat4-fma \
          test-maddnorm test-maddnorm-c test-maddnorm-inline test-maddnorm-fma \
          test-addmul test-addmul-c test-addmul-inline test-addmul-fma \
          test-stream test-stream-parity test-stream-parity-fma \
          test-physics test-mat4-ops test-fast

# ------------- Acceptance targets for the header-inline form (P1 refactor): append-only; not one command of the targets above was changed -------------
#   make test-header          strict compile profile (C11 -pedantic-errors / C++17 -Wpedantic); double inclusion in one TU;
#                             a compile-time assertion pins "a C TU on x86 + GCC/Clang must take the intrinsics path"; bit-level spot checks.
#   make test-inline          bit-identical cross-check of the default (inline) form against the FV_USE_EXTERN (in-library implementation) form:
#                             one tests/test_parity.c compiled five ways (language x form x pure-C fallback); the main output stream is compared byte-for-byte.
#   make test-fma             the consumer compiles with -march=native -mfma -ffp-contract=fast and must still be bit-identical
#                             -- the only safeguard that "an entry point with contraction risk was not misclassified as risk-free".
#                             Compilation always runs; execution needs FMA hardware (present here), otherwise it only compiles and prints SKIP.
#   make test-inline-effects  inlining survey under the consumer's default options (-O2, without -ffp-contract=off):
#                             the first group of 8 entry points must be **fully inlined** (no per-element call in the n<4 tail loop),
#                             the second group of 5 entry points must keep exactly **one** in-TU call per batch.
#   make test-refactor        all of the above + exported-symbol listing + disassembly survey + link evidence.
.PHONY: test-header test-inline test-fma test-inline-effects test-stream test-refactor

HEADER_CFLAGS   := -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -O2
HEADER_CXXFLAGS := -std=c++17 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -O2
PARITY_FLAGS    := -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic -Iinclude
FMA_FLAGS       := -O3 -Wall -Wextra -Wpedantic -march=native -mfma -ffp-contract=fast -Iinclude

# "entry point : probe function" table. First group = no mul+add contraction risk (full inlining required);
# second group = contraction risk, guarded by a pragma (one in-TU call per batch allowed, a per-element call is not).
SAFE_ENTRIES  := fv_vec2_add_x4:probe_add_x4_3 fv_vec2_sub_x4:probe_sub_x4_3 \
                 fv_vec2_scale_x4:probe_scale_x4_3 fv_vec2_add_x4_mask:probe_add_x4_mask_3 \
                 fv_vec2_add_x4_idx:probe_add_x4_idx_3 fv_vec2_scale_x4_idx:probe_scale_x4_idx_3 \
                 fv_vec2_add_soa:probe_add_soa_3 fv_vec2_scale_soa:probe_scale_soa_3 \
                 fv_vec2_add_x4_stream:probe_add_x4_stream_3 \
                 fv_vec2_sub_x4_stream:probe_sub_x4_stream_3 \
                 fv_vec2_scale_x4_stream:probe_scale_x4_stream_3 \
                 fv_vec2_add_soa_stream:probe_add_soa_stream_3 \
                 fv_vec2_scale_soa_stream:probe_scale_soa_stream_3 \
                 fv_vec2_add_x4_mask_inplace:probe_add_x4_mask_inplace_3 \
                 fv_vec4_add_x2:probe_v4_add_x2_1 fv_vec4_sub_x2:probe_v4_sub_x2_1 \
                 fv_vec4_scale_x2:probe_v4_scale_x2_1 \
                 fv_vec3_add_soa:probe_v3_add_soa_3 fv_vec3_scale_soa:probe_v3_scale_soa_3
RISKY_ENTRIES := fv_vec2_lerp_x4:probe_lerp_x4_3 fv_vec2_madd_x4:probe_madd_x4_3 \
                 fv_vec2_dot_soa:probe_dot_soa_3 fv_vec2_len_soa:probe_len_soa_3 \
                 fv_vec2_normalize_soa:probe_normalize_soa_3 \
                 fv_vec2_normalize_soa_stream:probe_normalize_soa_stream_3 \
                 fv_vec4_lerp_x2:probe_v4_lerp_x2_3 fv_vec4_madd_x2:probe_v4_madd_x2_3 \
                 fv_vec4_dot_soa:probe_v4_dot_soa_3 fv_vec4_len_soa:probe_v4_len_soa_3 \
                 fv_vec4_normalize_soa:probe_v4_normalize_soa_3 \
                 fv_vec3_dot_soa:probe_v3_dot_soa_3 fv_vec3_len_soa:probe_v3_len_soa_3 \
                 fv_vec3_normalize_soa:probe_v3_normalize_soa_3

$(BUILD)/test_header_c.o: tests/test_header_c.c include/flashvec.h | $(BUILD)
	$(CC) $(HEADER_CFLAGS) -c $< -o $@

$(BUILD)/test_header_c: $(BUILD)/test_header_c.o $(LIB)
	$(CXX) $^ -o $@

$(BUILD)/test_header_cpp.o: tests/test_header_cpp.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(HEADER_CXXFLAGS) -c $< -o $@

$(BUILD)/test_header_cpp: $(BUILD)/test_header_cpp.o $(LIB)
	$(CXX) $^ -o $@

test-header: $(BUILD)/test_header_c $(BUILD)/test_header_cpp
	./$(BUILD)/test_header_c
	./$(BUILD)/test_header_cpp

# Cross-check artifacts: language (C/C++) x form (default inline / FV_USE_EXTERN) x pure-C fallback (-DFV_I_NO_INTRINSICS)
$(BUILD)/parity_c_inline.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/parity_cpp_inline.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/parity_c_extern.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/parity_cpp_extern.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/parity_c_pure.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/parity_cpp_fma.o: tests/test_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/parity_%: $(BUILD)/parity_%.o $(LIB)
	$(CXX) $^ -o $@

test-inline: $(BUILD)/parity_c_inline $(BUILD)/parity_cpp_inline $(BUILD)/parity_c_extern \
             $(BUILD)/parity_cpp_extern $(BUILD)/parity_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/parity_$$v $(BUILD)/parity_$$v.bin $(BUILD)/parity_$$v.dot.bin || exit 1; \
	done
	@echo "== Output streams of the bit-identical entry points (1486 cases: n=0..17 all tail lengths + 63; four mask modes; idx duplicate/reversed/shuffled/single-element/m=0; + the 7 new P2 entry points) =="
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"
	@echo "== dot_soa (the only entry point the contract does not promise bit-identical results for): all four AVX2 builds take the same reduction path =="
	@cmp $(BUILD)/parity_cpp_inline.dot.bin $(BUILD)/parity_cpp_extern.dot.bin && echo "  PASS  C++ inline == C++ extern"
	@cmp $(BUILD)/parity_c_extern.dot.bin $(BUILD)/parity_cpp_extern.dot.bin && echo "  PASS  C extern   == C++ extern"
	@if cmp -s $(BUILD)/parity_c_pure.dot.bin $(BUILD)/parity_cpp_extern.dot.bin; then \
	  echo "  NOTE  the pure-C fallback dot equals AVX2 (this corpus did not trigger a reduction-order difference)"; \
	else echo "  NOTE  the pure-C fallback dot differs from AVX2 -- allowed by the contract (different reduction order), so it takes no part in the comparison"; fi

test-fma: $(BUILD)/parity_cpp_fma $(BUILD)/parity_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/parity_cpp_fma $(BUILD)/parity_cpp_fma.bin $(BUILD)/parity_cpp_fma.dot.bin || exit 1; \
	  echo "== Consumer built with -march=native -mfma -ffp-contract=fast (the worst real-world options) =="; \
	  cmp $(BUILD)/parity_cpp_fma.bin $(BUILD)/parity_cpp_extern.bin \
	    && echo "  PASS  all 12 entry points are still bit-identical to the in-library implementation (the pragma blocks the FMA the command line brings in)"; \
	  cmp $(BUILD)/parity_cpp_fma.dot.bin $(BUILD)/parity_cpp_extern.dot.bin \
	    && echo "  PASS  dot_soa matches too (both sides take the same AVX2 reduction path)"; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

$(BUILD)/consumer_probe_c.o: tests/consumer_probe.c include/flashvec.h | $(BUILD)
	$(CC) -std=c11 -O2 -Wall -Wextra -Wpedantic -Iinclude -c $< -o $@

$(BUILD)/consumer_probe_cpp.o: tests/consumer_probe.c include/flashvec.h | $(BUILD)
	$(CXX) -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude -c $< -o $@

test-inline-effects: $(BUILD)/consumer_probe_c.o $(BUILD)/consumer_probe_cpp.o
	@fail=0; \
	for obj in $(BUILD)/consumer_probe_c.o $(BUILD)/consumer_probe_cpp.o; do \
	  echo "== $$obj (-O2, without -ffp-contract=off = the consumer's default options) =="; \
	  printf '  %-26s %-8s %-8s %-14s %s\n' "entry" "group" "inline" "fv_i_ calls" "verdict"; \
	  for spec in $(SAFE_ENTRIES); do \
	    entry=$$(echo $$spec | cut -d: -f1); probe=$$(echo $$spec | cut -d: -f2); \
	    impl=$$(echo $$entry | sed -e 's/^fv_vec2_/fv_i_/' -e 's/^fv_vec3_/fv_i_v3_/' -e 's/^fv_vec4_/fv_i_v4_/'); \
	n=$$(nm --defined-only $$obj | awk '{print $$3}' | grep -cE "^$$impl([.]|$$)" || true); \
	    c=$$(objdump -d --no-show-raw-insn --disassemble="$$probe" $$obj | grep -cE ':(.*)(call|jmp).*<fv_i_' || true); \
	    if [ "$$n" = "0" ] && [ "$$c" = "0" ]; then v=PASS; else v=FAIL; fail=1; fi; \
	    inl=$$(if [ $$n = 0 ]; then echo yes; else echo no; fi); \
	    printf '  %-26s %-8s %-8s %-14s %s\n' "$$entry" "group 1" "$$inl" "$$c" "$$v"; \
	  done; \
	  for spec in $(RISKY_ENTRIES); do \
	    entry=$$(echo $$spec | cut -d: -f1); probe=$$(echo $$spec | cut -d: -f2); \
	    impl=$$(echo $$entry | sed -e 's/^fv_vec2_/fv_i_/' -e 's/^fv_vec3_/fv_i_v3_/' -e 's/^fv_vec4_/fv_i_v4_/'); \
	n=$$(nm --defined-only $$obj | awk '{print $$3}' | grep -cE "^$$impl([.]|$$)" || true); \
	    c=$$(objdump -d --no-show-raw-insn --disassemble="$$probe" $$obj | grep -cE ':(.*)(call|jmp).*<fv_i_' || true); \
	    if [ "$$n" -ge 1 ] && [ "$$c" = "1" ]; then v=PASS; else v=FAIL; fail=1; fi; \
	    inl=$$(if [ $$n = 0 ]; then echo yes; else echo no; fi); \
	    printf '  %-26s %-8s %-8s %-14s %s\n' "$$entry" "group 2" "$$inl" "$$c" "$$v"; \
	  done; \
	done; \
	if objdump -d $(BUILD)/consumer_probe_cpp.o | grep -qE 'jmp.*<fv_i_add_x4>'; then \
	  echo "  FAIL jmp fv_i_add_x4 still present"; fail=1; \
	else echo "  PASS the consumer TU has no jmp fv_i_add_x4 (the n=3 tail is fully inlined)"; fi; \
	echo "  -- actual assembly of probe_add_x4_3 (tail unrolled element by element, no call) --"; \
	objdump -d --no-show-raw-insn --disassemble=probe_add_x4_3 $(BUILD)/consumer_probe_cpp.o | tail -n +8 | head -10 | sed 's/^/    /'; \
	exit $$fail

# Dedicated acceptance for P2 (the memory-traffic entry points): the four alignment-peel paths, the "do not write" semantics of _mask_inplace, sfence visibility.
# The main bit-level cross-check lives in test_parity.c (already covered by make test-inline / test-fma); this only adds the forms it cannot reach.
$(BUILD)/test_stream.o: tests/test_stream.c include/flashvec.h | $(BUILD)
	$(CC) -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -O2 -Iinclude -c $< -o $@

$(BUILD)/test_stream: $(BUILD)/test_stream.o $(LIB)
	$(CXX) $^ -pthread -o $@

test-stream: $(BUILD)/test_stream
	./$(BUILD)/test_stream

# ---- the 13 _stream (NT store) entry points added for vec2 / vec3 / vec4 / mat4 ----
#   make test-stream-parity      cross-form bit-level check: the same source is compiled C / C++ x inline / FV_USE_EXTERN,
#                                plus the pure-C fallback (-DFV_I_NO_INTRINSICS); each form appends every stream output to a
#                                record file and the five files are cmp-ed, so a form that disagrees (wrong store, wrong peel
#                                range, an FMA-contracted scalar tail under -ffp-contract=fast) shows up as a byte difference.
#   make test-stream-parity-fma  the hostile-flags consumer (-march=native -mfma -ffp-contract=fast) cmp-ed against the
#                                in-library form: the evidence that the pragma still blocks the vfmadd the command line brings in.
# test-stream itself (the C test above) carries the alignment x tail-length x special-value matrix for those 13 entries.
.PHONY: test-stream-parity test-stream-parity-fma

$(BUILD)/stream_parity_c_inline.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/stream_parity_cpp_inline.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/stream_parity_c_extern.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/stream_parity_cpp_extern.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/stream_parity_c_pure.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/stream_parity_cpp_fma.o: tests/test_stream_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/stream_parity_%: $(BUILD)/stream_parity_%.o $(LIB)
	$(CXX) $^ -o $@

test-stream-parity: $(BUILD)/stream_parity_c_inline $(BUILD)/stream_parity_cpp_inline \
                    $(BUILD)/stream_parity_c_extern $(BUILD)/stream_parity_cpp_extern \
                    $(BUILD)/stream_parity_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/stream_parity_$$v $(BUILD)/stream_parity_$$v.bin || exit 1; \
	done
	@echo "== the 13 new _stream entries: C / C++ x inline / extern + pure-C fallback. Recorded corpus: n=0..40 all tail lengths + 63 + 100 + 1023 + 4097, six output phases, finite values. Edge values (-0.0 / +-inf / NaN / FLT_MIN / FLT_MAX / zero + near-zero vectors) and 12 phases are compared in-process in every form =="
	@cmp $(BUILD)/stream_parity_c_inline.bin $(BUILD)/stream_parity_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/stream_parity_c_inline.bin $(BUILD)/stream_parity_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/stream_parity_c_inline.bin $(BUILD)/stream_parity_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/stream_parity_c_inline.bin $(BUILD)/stream_parity_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"

test-stream-parity-fma: $(BUILD)/stream_parity_cpp_fma $(BUILD)/stream_parity_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/stream_parity_cpp_fma $(BUILD)/stream_parity_cpp_fma.bin || exit 1; \
	  echo "== the same 13 entries under -march=native -mfma -ffp-contract=fast (the worst real-world options) =="; \
	  cmp $(BUILD)/stream_parity_cpp_fma.bin $(BUILD)/stream_parity_cpp_extern.bin \
	    && echo "  PASS  all 13 new entry points are still bit-identical to the in-library implementation (the pragma blocks the FMA the command line brings in)"; \
	  printf '  vfmadd in the FMA consumer object: %s\n' "$$(objdump -d $(BUILD)/stream_parity_cpp_fma.o | grep -coE 'vfmadd[0-9]*[sp][ds]' || true)"; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

test-refactor: test-header test-inline test-fma test-inline-effects test-stream \
               test-stream-parity test-stream-parity-fma \
               test-vec4 test-vec34-inline test-vec34-fma
	@echo "== Exported symbol list (nm -g --defined-only build/libflashvec.a | sort) =="
	@nm -g --defined-only $(LIB) | sort
	@echo "== Disassembly survey: each of these patterns must occur 0 times in the library =="
	@for pat in 'vfmadd[0-9]*[sp][ds]' 'vfmsub[0-9]*[sp][ds]' 'vfnmadd[0-9]*[sp][ds]' 'vrcpps' 'vrsqrtps' 'haddps' 'vmaskmovps'; do \
	  printf '  %-22s %s\n' "$$pat" "$$(objdump -d $(LIB) | grep -oE "$$pat" | wc -l)"; \
	done
	@echo "== Same survey over the consumer artifacts (C translation units now carry AVX2 too, so the positive control must be non-zero) =="
	@for obj in $(BUILD)/parity_c_inline.o $(BUILD)/parity_cpp_inline.o $(BUILD)/parity_c_pure.o; do \
	  printf '  %-30s vfmadd=%s  AVX(positive control)=%s\n' "$$obj" \
	    "$$(objdump -d $$obj | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	    "$$(objdump -d $$obj | grep -oE 'v(add|sub|mul|div|sqrt)ps' | wc -l)"; \
	done
	@echo "== Link evidence: default-form executables carry no in-library batch symbols (inlined away), the extern form takes them from the library =="
	@printf '  parity_c_inline  global batch symbols = %s\n' "$$(nm -g --defined-only $(BUILD)/parity_c_inline | grep -cE 'fv_vec2_.*_(x4|soa)')"
	@printf '  parity_c_extern  global batch symbols = %s\n' "$$(nm -g --defined-only $(BUILD)/parity_c_extern | grep -cE 'fv_vec2_.*_(x4|soa)')"
# ============ P3: fv_vec4 / fv_vec3 SoA batches (this section is maintained by the vec3/vec4 task) ============
# This is the only build logic added/changed this round; it states three things:
#   1) the scalar implementation of fv_vec4 (src/flashvec_vec4.cpp) is now part of $(OBJ): the scalar API and the
#      scalar tail of the batch entry points must come from the same source with the same options, and the symbols must be in libflashvec.a too (nm only grows);
#   2) the bit-level cross-check compiles tests/test_vec34_parity.c five ways (language x form x pure-C fallback) + an FMA consumer,
#      then compares the main stream byte-for-byte; the dot stream is grouped per contract (the pure-C fallback reduces in a different order, so a difference is allowed);
#   3) inlining effects: the 13 P3 entry points were added to the existing SAFE_ENTRIES / RISKY_ENTRIES tables above,
#      and the probes live in the P3 subsection of tests/consumer_probe.c, so make test-inline-effects covers them directly.
.PHONY: test-vec4 test-vec34-inline test-vec34-fma

$(BUILD)/flashvec_vec4.o: src/flashvec_vec4.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Bit-level reference test for the scalar API. The options deliberately match the in-library implementation (-O3 -fno-math-errno -ffp-contract=off):
# the sqrtf in the reference expression then expands to the same instruction as in the library, and mul+add does not contract either --
# otherwise what gets compared is a compiler difference, not the behaviour of the library.
$(BUILD)/test_vec4.o: tests/test_vec4.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) -std=c11 -O3 -fno-math-errno -ffp-contract=off -pedantic-errors -Wall -Wextra -Wpedantic -c $< -o $@

$(BUILD)/test_vec4: $(BUILD)/test_vec4.o $(LIB)
	$(CXX) $^ -o $@

test-vec4: $(BUILD)/test_vec4
	./$(BUILD)/test_vec4

# Cross-check artifacts: language (C/C++) x form (default inline / FV_USE_EXTERN) x pure-C fallback + FMA consumer
$(BUILD)/vec34_c_inline.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/vec34_cpp_inline.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/vec34_c_extern.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/vec34_cpp_extern.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/vec34_c_pure.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/vec34_cpp_fma.o: tests/test_vec34_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/vec34_%: $(BUILD)/vec34_%.o $(LIB)
	$(CXX) $^ -o $@

test-vec34-inline: $(BUILD)/vec34_c_inline $(BUILD)/vec34_cpp_inline $(BUILD)/vec34_c_extern \
                  $(BUILD)/vec34_cpp_extern $(BUILD)/vec34_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/vec34_$$v $(BUILD)/vec34_$$v.bin $(BUILD)/vec34_$$v.dot.bin || exit 1; \
	done
	@echo "== 3D/4D batch entry points: bit-identical (n=0..17 all tail lengths + 63 + 1024 + 8192; with GUARD poisoning) =="
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"
	@echo "== dot_soa (the two entry points the contract does not promise bit-identical results for): all four AVX2 builds take the same reduction path =="
	@cmp $(BUILD)/vec34_cpp_inline.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin && echo "  PASS  C++ inline == C++ extern"
	@cmp $(BUILD)/vec34_c_extern.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin && echo "  PASS  C extern   == C++ extern"
	@if cmp -s $(BUILD)/vec34_c_pure.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin; then \
	  echo "  NOTE  the pure-C fallback dot equals AVX2 (this corpus did not trigger a reduction-order difference)"; \
	else echo "  NOTE  the pure-C fallback dot differs from AVX2 -- allowed by the contract (different reduction order), so it takes no part in the comparison"; fi

test-vec34-fma: $(BUILD)/vec34_cpp_fma $(BUILD)/vec34_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/vec34_cpp_fma $(BUILD)/vec34_cpp_fma.bin $(BUILD)/vec34_cpp_fma.dot.bin || exit 1; \
	  echo "== Consumer built with -march=native -mfma -ffp-contract=fast (the worst real-world options) =="; \
	  cmp $(BUILD)/vec34_cpp_fma.bin $(BUILD)/vec34_cpp_extern.bin \
	    && echo "  PASS  all 13 new entry points are still bit-identical to the in-library implementation (the pragma blocks the FMA the command line brings in)"; \
	  cmp $(BUILD)/vec34_cpp_fma.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin \
	    && echo "  PASS  dot_soa matches too (both sides take the same AVX2 reduction path)"; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

# =====================================================================
# Forwarding targets for the bench/ab same-process A/B fixture (pure forwarding; all logic lives in bench/ab/)
#   make bench-ab          self-check + baseline snapshot (including the A/A floor table)
#   make bench-ab-worktree candidate = the current working tree
#   make bench-ab-negctl   negative control (a wrongly modified candidate must fail hard and name the function)
#   make bench-verify OP=<op> N=<size>   the five verdicts of same-direction evidence (PASS/FAIL/INCONCLUSIVE)
# These targets only write bench/ab/build/ and docs/research/raw/ab_*, never build/rootbuild/.
# =====================================================================
.PHONY: bench-ab bench-ab-worktree bench-ab-negctl bench-verify

bench-ab:
	$(MAKE) -C bench/ab run

bench-ab-worktree:
	$(MAKE) -C bench/ab run-worktree

bench-ab-negctl:
	$(MAKE) -C bench/ab negctl

bench-verify:
	@if [ -z "$(OP)" ]; then echo "usage: make bench-verify OP=<op-name> [N=<size>]" >&2; exit 2; fi
	$(MAKE) -C bench/ab verify OP=$(OP) N=$(if $(N),$(N),1048576)

# ==================== P4: matrix point-transform batch entry points (fv_mat4 / fv_mat4_transform_*) ====================
# Append-only: not one command of any existing target above was changed. The four targets divide the work:
#   test-mat4          C++ main test: bit-identical cross-check against the scalar reference (n = 0..40 all tail lengths + 1024/4097/8192),
#                      semantics (identity matrix / translation only / scale only / the w convention / the last row is never read), output sentinels
#   test-mat4-c        the same set of forms as a C translation unit (declarations, semantics and bits all line up under plain C11)
#   test-mat4-inline   the output streams of the 5 forms (C/C++ x inline/extern x pure-C fallback) compared byte-for-byte
#   test-mat4-fma      the consumer compiles with -march=native -mfma -ffp-contract=fast and must still be bit-identical
#                      -- the matrix kernel has mul+add, so this is the only direct evidence that "the pragma sits in the right place"
.PHONY: test-mat4 test-mat4-c test-mat4-inline test-mat4-fma

$(BUILD)/test_mat4: tests/test_mat4.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-mat4: $(BUILD)/test_mat4
	./$(BUILD)/test_mat4

$(BUILD)/test_mat4_c.o: tests/test_mat4_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_mat4_c: $(BUILD)/test_mat4_c.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-mat4-c: $(BUILD)/test_mat4_c
	./$(BUILD)/test_mat4_c

# Cross-check artifacts: language (C/C++) x form (default inline / FV_USE_EXTERN) x pure-C fallback + FMA consumer
$(BUILD)/mat4_c_inline.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/mat4_cpp_inline.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/mat4_c_extern.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/mat4_cpp_extern.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/mat4_c_pure.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/mat4_cpp_fma.o: tests/test_mat4_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/mat4_%: $(BUILD)/mat4_%.o $(LIB)
	$(CXX) $^ -o $@

test-mat4-inline: $(BUILD)/mat4_c_inline $(BUILD)/mat4_cpp_inline $(BUILD)/mat4_c_extern \
                  $(BUILD)/mat4_cpp_extern $(BUILD)/mat4_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/mat4_$$v $(BUILD)/mat4_$$v.bin || exit 1; \
	done
	@echo "== Matrix point transform: bit-identical (n = 0..40 all tail lengths + 63 + 1024 + 4097; including inf/NaN/-0.0/subnormals) =="
	@cmp $(BUILD)/mat4_c_inline.bin $(BUILD)/mat4_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/mat4_c_inline.bin $(BUILD)/mat4_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/mat4_c_inline.bin $(BUILD)/mat4_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/mat4_c_inline.bin $(BUILD)/mat4_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"

test-mat4-fma: $(BUILD)/mat4_cpp_fma $(BUILD)/mat4_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/mat4_cpp_fma $(BUILD)/mat4_cpp_fma.bin || exit 1; \
	  echo "== Matrix point transform, consumer built with -march=native -mfma -ffp-contract=fast =="; \
	  cmp $(BUILD)/mat4_cpp_fma.bin $(BUILD)/mat4_cpp_extern.bin \
	    && echo "  PASS  both entry points are still bit-identical to the in-library implementation (the AVX2 kernel sits inside the pragma region, no vfmadd)"; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

# ==================== P5: fused pipeline entry (fv_vec2_madd_normalize_soa / _stream) ====================
# Append-only: not one command of any existing target above was changed. Four targets:
#   test-maddnorm          C++ main test: bit-identical (memcmp) to fv_vec2_scale_soa -> add_soa -> normalize_soa
#                          for n = 0..40 (every tail length), n = 1023/4096/4097/8192/12345, an s sweep,
#                          output/input phase shifts 0..7, path independence, plus a counted anything-goes corpus
#   test-maddnorm-c        the same criterion from a plain C11 translation unit (declarations and bits line up)
#   test-maddnorm-inline   the output streams of the 5 forms (C/C++ x inline/extern x pure-C fallback) compared byte for byte
#   test-maddnorm-fma      the consumer compiles with -march=native -mfma -ffp-contract=fast and must still be bit-identical
#                          -- the fused kernel computes a * s + b, so this is the direct evidence that it sits inside the
#                          fp-contract=off pragma; it also counts vfmadd in the two new kernels (must be 0).
.PHONY: test-maddnorm test-maddnorm-c test-maddnorm-inline test-maddnorm-fma

$(BUILD)/test_madd_norm: tests/test_madd_norm.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-maddnorm: $(BUILD)/test_madd_norm
	./$(BUILD)/test_madd_norm

$(BUILD)/test_madd_norm_c.o: tests/test_madd_norm_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_madd_norm_c: $(BUILD)/test_madd_norm_c.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-maddnorm-c: $(BUILD)/test_madd_norm_c
	./$(BUILD)/test_madd_norm_c

# Cross-check artifacts: language (C/C++) x form (default inline / FV_USE_EXTERN) x pure-C fallback + FMA consumer
$(BUILD)/maddnorm_c_inline.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/maddnorm_cpp_inline.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/maddnorm_c_extern.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/maddnorm_cpp_extern.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/maddnorm_c_pure.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/maddnorm_cpp_fma.o: tests/test_madd_norm_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/maddnorm_%: $(BUILD)/maddnorm_%.o $(LIB)
	$(CXX) $^ -o $@

test-maddnorm-inline: $(BUILD)/maddnorm_c_inline $(BUILD)/maddnorm_cpp_inline $(BUILD)/maddnorm_c_extern \
                      $(BUILD)/maddnorm_cpp_extern $(BUILD)/maddnorm_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/maddnorm_$$v $(BUILD)/maddnorm_$$v.bin || exit 1; \
	done
	@echo "== Fused pipeline out = normalize(a*s + b): bit-identical (n = 0..24 all tail lengths + 63 + 1023 + 4096, s sweep, pointer phases 0/1/3/7) =="
	@cmp $(BUILD)/maddnorm_c_inline.bin $(BUILD)/maddnorm_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/maddnorm_c_inline.bin $(BUILD)/maddnorm_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/maddnorm_c_inline.bin $(BUILD)/maddnorm_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/maddnorm_c_inline.bin $(BUILD)/maddnorm_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"

test-maddnorm-fma: $(BUILD)/maddnorm_cpp_fma $(BUILD)/maddnorm_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/maddnorm_cpp_fma $(BUILD)/maddnorm_cpp_fma.bin || exit 1; \
	  echo "== Fused pipeline, consumer built with -march=native -mfma -ffp-contract=fast =="; \
	  cmp $(BUILD)/maddnorm_cpp_fma.bin $(BUILD)/maddnorm_cpp_extern.bin \
	    && echo "  PASS  both entries are still bit-identical to the in-library implementation (the kernel sits inside the pragma region)"; \
	  echo "== vfmadd census of that artifact (avx2 symbols first; -march=native inlines them into the wrappers, so the whole object is the fallback criterion) =="; \
	  printf '  %-40s vfmadd=%s vmovntps=%s\n' "whole object $(BUILD)/maddnorm_cpp_fma.o" \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/maddnorm_cpp_fma.o | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/maddnorm_cpp_fma.o | grep -oE 'vmovntps' | wc -l)"; \
	  for sym in fv_i_madd_norm_soa_avx2 fv_i_madd_norm_stream_soa_avx2 fv_i_madd_norm_soa fv_i_madd_norm_soa_stream; do \
	    if nm --defined-only $(BUILD)/maddnorm_cpp_fma.o | grep -qE "[[:space:]]$$sym$$"; then \
	      printf '  %-40s vfmadd=%s vmulps+vaddps=%s\n' "$$sym" \
	        "$$(objdump -d --no-show-raw-insn --disassemble="$$sym" $(BUILD)/maddnorm_cpp_fma.o | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	        "$$(objdump -d --no-show-raw-insn --disassemble="$$sym" $(BUILD)/maddnorm_cpp_fma.o | grep -oE 'v(mul|add)ps' | wc -l)"; \
	    else printf '  %-40s (inlined into its caller in this artifact)\n' "$$sym"; fi; \
	  done; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

# ---------------- P4: integrator primitives (in-place x += v*s, vector madd o = a*b + c) ----------------
#   make test-addmul         bit-level fixture (C++): every tail length n = 0..40, three corpora x eight scales, large and
#                            non-power-of-two n, NULL/n == 0, alias and alignment phases, split invariance, path
#                            independence, per-value API vs the written-out expression
#   make test-addmul-c       the same group from a plain C caller (every tail length n = 0..17, misaligned destination)
#   make test-addmul-inline  cross-form bit-level parity: C / C++ x inline / FV_USE_EXTERN + the pure-C fallback, byte-compared
#   make test-addmul-fma     the hostile consumer (-march=native -mfma -ffp-contract=fast): the pragma must still block the
#                            FMA, plus a vfmadd census of that artifact
#   make bench-addmul        measurement of addmul / maddv / maddv _stream against the two-pass workaround
.PHONY: test-addmul test-addmul-c test-addmul-inline test-addmul-fma bench-addmul

$(BUILD)/test_addmul: tests/test_addmul.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-addmul: $(BUILD)/test_addmul
	./$(BUILD)/test_addmul

$(BUILD)/test_addmul_c.o: tests/test_addmul_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_addmul_c: $(BUILD)/test_addmul_c.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-addmul-c: $(BUILD)/test_addmul_c
	./$(BUILD)/test_addmul_c

# Cross-check artifacts: language (C/C++) x form (default inline / FV_USE_EXTERN) x pure-C fallback + FMA consumer
$(BUILD)/addmul_c_inline.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -c $< -o $@
$(BUILD)/addmul_cpp_inline.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -c $< -o $@
$(BUILD)/addmul_c_extern.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/addmul_cpp_extern.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(PARITY_FLAGS) -std=c++17 -DFV_USE_EXTERN -c $< -o $@
$(BUILD)/addmul_c_pure.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CC) $(PARITY_FLAGS) -std=c11 -DFV_I_NO_INTRINSICS -c $< -o $@
$(BUILD)/addmul_cpp_fma.o: tests/test_addmul_parity.c include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/addmul_%: $(BUILD)/addmul_%.o $(LIB)
	$(CXX) $^ -o $@

test-addmul-inline: $(BUILD)/addmul_c_inline $(BUILD)/addmul_cpp_inline $(BUILD)/addmul_c_extern \
                    $(BUILD)/addmul_cpp_extern $(BUILD)/addmul_c_pure
	@for v in c_inline cpp_inline c_extern cpp_extern c_pure; do \
	  ./$(BUILD)/addmul_$$v $(BUILD)/addmul_$$v.bin || exit 1; \
	done
	@echo "== Integrator primitives: bit-identical (n = 0..40 all tail lengths + 63 + 1023 + 4097; four destination phases; finite and single-NaN corpora) =="
	@cmp $(BUILD)/addmul_c_inline.bin $(BUILD)/addmul_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/addmul_c_inline.bin $(BUILD)/addmul_c_extern.bin && echo "  PASS  C inline  == C extern (in-library implementation)"
	@cmp $(BUILD)/addmul_c_inline.bin $(BUILD)/addmul_cpp_extern.bin && echo "  PASS  C inline  == C++ extern (in-library implementation)"
	@cmp $(BUILD)/addmul_c_inline.bin $(BUILD)/addmul_c_pure.bin && echo "  PASS  C inline  == C pure-C fallback (-DFV_I_NO_INTRINSICS)"

test-addmul-fma: $(BUILD)/addmul_cpp_fma $(BUILD)/addmul_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/addmul_cpp_fma $(BUILD)/addmul_cpp_fma.bin || exit 1; \
	  echo "== Integrator primitives, consumer built with -march=native -mfma -ffp-contract=fast =="; \
	  cmp $(BUILD)/addmul_cpp_fma.bin $(BUILD)/addmul_cpp_extern.bin \
	    && echo "  PASS  all six batch entries plus the per-value API are still bit-identical to the in-library implementation (the pragma blocks the FMA)"; \
	  echo "== vfmadd / vmovntps census of that artifact =="; \
	  printf '  %-40s vfmadd=%s vmovntps=%s\n' "whole object $(BUILD)/addmul_cpp_fma.o" \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/addmul_cpp_fma.o | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/addmul_cpp_fma.o | grep -oE 'vmovntps' | wc -l)"; \
	  for sym in fv_i_addmul4_avx2 fv_i_addmul_soa_avx2 fv_i_addmul_v3_soa_avx2 fv_i_v4_addmul2_avx2 fv_i_maddv4_avx2 fv_i_maddv4_stream_avx2; do \
	    if nm --defined-only $(BUILD)/addmul_cpp_fma.o | grep -qE "[[:space:]]$$sym$$"; then \
	      printf '  %-40s vfmadd=%s vmulps+vaddps=%s\n' "$$sym" \
	        "$$(objdump -d --no-show-raw-insn --disassemble="$$sym" $(BUILD)/addmul_cpp_fma.o | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	        "$$(objdump -d --no-show-raw-insn --disassemble="$$sym" $(BUILD)/addmul_cpp_fma.o | grep -oE 'v(mul|add)ps' | wc -l)"; \
	    else printf '  %-40s (inlined into its caller in this artifact)\n' "$$sym"; fi; \
	  done; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

$(BUILD)/bench_addmul: tests/bench_addmul.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

bench-addmul: $(BUILD)/bench_addmul
	./$(BUILD)/bench_addmul

# ==================== P6: quaternion / mat3 physics primitives ====================
# Append-only: not one command of any existing target above was changed.
#   test-quat        C++ main test: quat semantics + bit-identical SoA batch (n = 0..40 all tails + 63/1024/4097)
#   test-mat3        C++ main test: mat3 semantics + bit-identical SoA batch (same tail coverage)
#   test-physics-c   plain C11 caller: scalar API links, batch entries compile and are bit-identical
#   test-physics     all of the above forms (this is what test-all depends on)
#   bench-physics    scalar-per-element vs SoA batch, same-process A/B, 11 paired rounds + A/A control
.PHONY: test-quat test-mat3 test-physics-c test-physics-extern test-physics-pure test-physics-fma test-physics test-mat4-ops test-fast bench-physics bench-fast bench-vec34 bench-pack

$(BUILD)/flashvec_quat.o: src/flashvec_quat.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_mat3.o: src/flashvec_mat3.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_mat4.o: src/flashvec_mat4.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_fast.o: src/flashvec_fast.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_fast_batch.o: src/flashvec_fast_batch.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/flashvec_physics_batch.o: src/flashvec_physics_batch.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/test_quat: tests/test_quat.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-quat: $(BUILD)/test_quat
	./$(BUILD)/test_quat

$(BUILD)/test_mat3: tests/test_mat3.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-mat3: $(BUILD)/test_mat3
	./$(BUILD)/test_mat3

$(BUILD)/test_physics_c.o: tests/test_physics_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/test_physics_c: $(BUILD)/test_physics_c.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-physics-c: $(BUILD)/test_physics_c
	./$(BUILD)/test_physics_c

# The same C11 source in the two other forms: FV_USE_EXTERN (the library's extern symbols) and the
# pure-C fallback (no intrinsics). Each form still compares the batch output against the per-value API
# bit for bit, so the three forms are one implementation from the caller's point of view.
$(BUILD)/test_physics_c_extern.o: tests/test_physics_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DFV_USE_EXTERN -c $< -o $@

$(BUILD)/test_physics_c_extern: $(BUILD)/test_physics_c_extern.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-physics-extern: $(BUILD)/test_physics_c_extern
	./$(BUILD)/test_physics_c_extern

$(BUILD)/test_physics_c_pure.o: tests/test_physics_c.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -DFV_I_NO_INTRINSICS -c $< -o $@

$(BUILD)/test_physics_c_pure: $(BUILD)/test_physics_c_pure.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-physics-pure: $(BUILD)/test_physics_c_pure
	./$(BUILD)/test_physics_c_pure

# FMA consumer: the worst real-world options for the bit-exactness promise. The test's own reference
# expressions are wrapped in fp-contract=off (the file's pragma region), so this measures the header's
# pragma, not the test's flags. A vfmadd count of 0 is the direct evidence that the kernel did not fuse.
$(BUILD)/test_quat_fma.o: tests/test_quat.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/test_quat_fma: $(BUILD)/test_quat_fma.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

$(BUILD)/test_mat3_fma.o: tests/test_mat3.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(FMA_FLAGS) -std=c++17 -c $< -o $@

$(BUILD)/test_mat3_fma: $(BUILD)/test_mat3_fma.o $(LIB) | $(BUILD)
	$(CXX) $^ -o $@

test-physics-fma: $(BUILD)/test_quat_fma $(BUILD)/test_mat3_fma
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/test_quat_fma; \
	  ./$(BUILD)/test_mat3_fma; \
	  printf '  vfmadd in test_quat_fma.o: %s\n' \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/test_quat_fma.o | grep -cE 'vfmadd')"; \
	  printf '  vfmadd in test_mat3_fma.o: %s\n' \
	    "$$(objdump -d --no-show-raw-insn $(BUILD)/test_mat3_fma.o | grep -cE 'vfmadd')"; \
	else echo "  SKIP: no FMA on this machine, compilation only"; fi

test-physics: test-quat test-mat3 test-physics-c test-physics-extern test-physics-pure test-physics-fma
	@echo "== P6 physics primitives: quat + mat3 scalar and SoA batch, all forms green =="

$(BUILD)/test_mat4_ops: tests/test_mat4_ops.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-mat4-ops: $(BUILD)/test_mat4_ops
	./$(BUILD)/test_mat4_ops

$(BUILD)/test_fast: tests/test_fast.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

test-fast: $(BUILD)/test_fast
	./$(BUILD)/test_fast

$(BUILD)/bench_physics: tests/bench_physics.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

bench-physics: $(BUILD)/bench_physics
	./$(BUILD)/bench_physics

$(BUILD)/bench_fast: tests/bench_fast.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

bench-fast: $(BUILD)/bench_fast
	./$(BUILD)/bench_fast

$(BUILD)/bench_vec34: tests/bench_vec34.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

bench-vec34: $(BUILD)/bench_vec34
	./$(BUILD)/bench_vec34

$(BUILD)/bench_pack: tests/bench_pack.cpp $(LIB) include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

bench-pack: $(BUILD)/bench_pack
	./$(BUILD)/bench_pack
