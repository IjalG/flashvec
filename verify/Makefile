# verify/Makefile -- verifier F's independent build (writes only under verify/; never touches include/ src/ tests/ bench/ or the root Makefile)
#
#   make parity   -- bit-identical verification: build the library once at -O2 and once at -O3 -march=native, run each, plus ASan/UBSan
#   make perf     -- benchmarks: the -O3 -march=native and default -O2 tiers, raw output written to build/
#   make check    -- pure-C11 header/ABI check + the repository's existing tests (recompiled here) + structure checks
#   make extra    -- implementer E's own test suite (compiled and run only, never modified)
#   make all      -- all of the above

CXX ?= c++
CC  ?= cc
INC := -I../include
WARN := -Wall -Wextra
SRC1 := ../src/flashvec.cpp
SRC2 := ../src/flashvec_batch.cpp
B := build

# The verification side (parity.cpp / ref_strict.cpp) always uses -ffp-contract=off: my own reference values and expected values
# must not be contracted into FMA by the compiler, which would produce false positives.
PARITY_FLAGS := $(WARN) $(INC) -ffp-contract=off
LIB_O2  := -std=c++17 -O2
LIB_O3N := -std=c++17 -O3 -march=native
OBJ_O2  := $(B)/flashvec_O2.o $(B)/flashvec_batch_O2.o
OBJ_O3N := $(B)/flashvec_O3n.o $(B)/flashvec_batch_O3n.o
# -O3 without -march=native: aligned with the root Makefile's current release default (-O3 -fno-math-errno -ffp-contract=off)
DEF_O3  := -std=c++17 -O3 -fno-math-errno -ffp-contract=off
OBJ_O3  := $(B)/flashvec_O3.o $(B)/flashvec_batch_O3.o

.PHONY: all parity perf check extra clean run-parity run-perf

all: parity perf check extra

$(B):
	mkdir -p $(B)

# ---------------- The library under test, built at two tiers (two TUs) ----------------
$(B)/flashvec_O2.o: $(SRC1) ../include/flashvec.h | $(B)
	$(CXX) $(LIB_O2) $(WARN) $(INC) -c $< -o $@
$(B)/flashvec_batch_O2.o: $(SRC2) ../include/flashvec.h | $(B)
	$(CXX) $(LIB_O2) $(WARN) $(INC) -c $< -o $@
$(B)/flashvec_O3n.o: $(SRC1) ../include/flashvec.h | $(B)
	$(CXX) $(LIB_O3N) $(WARN) $(INC) -c $< -o $@
$(B)/flashvec_batch_O3n.o: $(SRC2) ../include/flashvec.h | $(B)
	$(CXX) $(LIB_O3N) $(WARN) $(INC) -c $< -o $@

$(B)/flashvec_O3.o: $(SRC1) ../include/flashvec.h | $(B)
	$(CXX) $(DEF_O3) $(WARN) $(INC) -c $< -o $@
$(B)/flashvec_batch_O3.o: $(SRC2) ../include/flashvec.h | $(B)
	$(CXX) $(DEF_O3) $(WARN) $(INC) -c $< -o $@

$(B)/ref_strict.o: ref_strict.cpp ref_strict.h ../include/flashvec.h | $(B)
	$(CXX) -std=c++17 -O2 -ffp-contract=off -fno-tree-vectorize $(WARN) $(INC) -c $< -o $@

# ---------------- Bit-identical verification ----------------
$(B)/parity_O2_plain: parity.cpp ref_strict.h ../include/flashvec.h $(OBJ_O2) $(B)/ref_strict.o | $(B)
	$(CXX) -std=c++17 -O2 $(PARITY_FLAGS) parity.cpp $(B)/ref_strict.o $(OBJ_O2) -lm -o $@

$(B)/parity_O3_native: parity.cpp ref_strict.h ../include/flashvec.h $(OBJ_O3N) $(B)/ref_strict.o | $(B)
	$(CXX) -std=c++17 -O3 -march=native $(PARITY_FLAGS) parity.cpp $(B)/ref_strict.o $(OBJ_O3N) -lm -o $@

$(B)/parity_asan: parity.cpp ref_strict.h ../include/flashvec.h $(SRC1) $(SRC2) ref_strict.cpp | $(B)
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(PARITY_FLAGS) \
		parity.cpp ref_strict.cpp $(SRC1) $(SRC2) -lm -o $@

parity: $(B)/parity_O2_plain $(B)/parity_O3_native $(B)/parity_asan
	@echo '================ parity (-O2, default build tier: no -march=native) ================'
	@$(B)/parity_O2_plain > $(B)/parity_O2_plain.txt 2>&1; rc=$$?; cat $(B)/parity_O2_plain.txt; echo "exit=$$rc"
	@echo '================ parity (-O3 -march=native) ================'
	@$(B)/parity_O3_native > $(B)/parity_O3_native.txt 2>&1; rc=$$?; cat $(B)/parity_O3_native.txt; echo "exit=$$rc"
	@echo '================ parity (ASan+UBSan) ================'
	@$(B)/parity_asan > $(B)/parity_asan.txt 2>&1; rc=$$?; echo "exit=$$rc"; grep -E 'FAIL|ERROR|SUMMARY|runtime error' $(B)/parity_asan.txt | head -10; tail -2 $(B)/parity_asan.txt

run-parity: parity

# ---------------- Benchmarks ----------------
$(B)/perf_O3_native: perf.cpp ../include/flashvec.h $(OBJ_O3N) | $(B)
	$(CXX) -std=c++17 -O3 -march=native -ffp-contract=off $(WARN) $(INC) perf.cpp $(OBJ_O3N) -lm -o $@

$(B)/perf_O2_plain: perf.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CXX) -std=c++17 -O2 -ffp-contract=off $(WARN) $(INC) perf.cpp $(OBJ_O2) -lm -o $@

# Measurement of the release default (without -march=native) -- answers "does the default build still get an AVX2 win"
$(B)/perf_O3_plain: perf.cpp ../include/flashvec.h $(OBJ_O3) | $(B)
	$(CXX) $(DEF_O3) $(WARN) $(INC) perf.cpp $(OBJ_O3) -lm -o $@

# Pinned core + unbuffered output: avoids both the migration artifact and the "a crash loses the buffer" artifact
perf: $(B)/perf_O3_native $(B)/perf_O3_plain $(B)/perf_O2_plain
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (-O3 -march=native) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O3_native > $(B)/perf_O3_native.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O3_native.txt
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (release default -O3, no -march=native) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O3_plain > $(B)/perf_O3_plain.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O3_plain.txt
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (-O2 without -march=native, the old default tier) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O2_plain > $(B)/perf_O2_plain.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O2_plain.txt

run-perf: perf

# ---------------- Structure / pure C / existing tests ----------------
$(B)/c_api_check_O2: c_api_check.c ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CC) -std=c11 -O3 -pedantic-errors $(WARN) $(INC) -ffp-contract=off c_api_check.c $(OBJ_O2) -lm -o $@

$(B)/test_flashvec_F: ../tests/test_flashvec.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@

$(B)/test_c_api_F.o: ../tests/test_c_api.c ../include/flashvec.h | $(B)
	$(CC) -std=c11 -O2 -pedantic-errors $(WARN) $(INC) -c $< -o $@

$(B)/test_c_api_F: $(B)/test_c_api_F.o $(OBJ_O2)
	$(CXX) $^ -lm -o $@

check: $(B)/c_api_check_O2 $(B)/test_flashvec_F $(B)/test_c_api_F
	@echo '================ Pure-C11 header + ABI check (-pedantic-errors) ================'
	@$(B)/c_api_check_O2 > $(B)/c_api_check.txt 2>&1; rc=$$?; cat $(B)/c_api_check.txt; echo "exit=$$rc"
	@echo '================ Existing tests from the repository (recompiled here, build/ untouched) ================'
	@$(B)/test_flashvec_F > $(B)/test_flashvec_F.txt 2>&1; rc=$$?; cat $(B)/test_flashvec_F.txt; echo "exit=$$rc"
	@$(B)/test_c_api_F > $(B)/test_c_api_F.txt 2>&1; rc=$$?; cat $(B)/test_c_api_F.txt; echo "exit=$$rc"
	@./checks.sh

# ---------------- Implementer E own test suite (compiled and run only, never modified) ----------------
$(B)/test_batch_F: ../tests/test_batch.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	-$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@
$(B)/test_batch_c_F.o: ../tests/test_batch_c.c ../include/flashvec.h | $(B)
	-$(CC) -std=c11 -O2 -pedantic-errors $(WARN) $(INC) -c $< -o $@
$(B)/test_batch_c_F: $(B)/test_batch_c_F.o $(OBJ_O2)
	-$(CXX) $^ -lm -o $@
$(B)/test_p0_F: ../tests/test_p0.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	-$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@
extra: $(B)/test_batch_F $(B)/test_batch_c_F $(B)/test_p0_F
	@echo '================ Implementer E test suite (recompiled here and run independently) ================'
	@for t in test_batch_F test_batch_c_F test_p0_F; do if [ -x $(B)/$$t ]; then $(B)/$$t > $(B)/$$t.txt 2>&1; echo "$$t exit=$$?"; tail -6 $(B)/$$t.txt; fi; done

clean:
	rm -rf $(B)
