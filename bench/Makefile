# bench/Makefile -- flashvec three-implementation benchmark (new file; does not change any existing target of the main Makefile)
#
# Usage:
#   make -C bench run    build + run every tier; raw output goes to ../docs/research/raw/
#   make -C bench table  generate the Markdown comparison table from the raw logs (bench/table.awk)
#   make -C bench vec    capture the GCC auto-vectorization reports (-fopt-info-vec)
#   make -C bench asm    disassembly archive (to inspect call overhead / ABI argument passing)
#   make -C bench perf   perf stat for IPC / cache-miss / branch-miss (limited by perf_event_paranoid)
#   make -C bench clean  remove bench/build

ROOT   := ..
BUILD  := build
RAW    := $(ROOT)/docs/research/raw
CXX    ?= g++
CC     ?= cc
AWK    ?= awk
PERF   ?= perf
CPU    ?= 2

CPPFLAGS += -I$(ROOT)/include

# The three sets of flags under test (-march=native expands to AVX2 + FMA on this machine; see the report)
NATIVE   := -O3 -march=native
FFAST    := -O3 -march=native -ffast-math
SHIPPED  := -O2
# Exactly one strictly semantics-safe switch, to isolate "was it fast-math or math-errno that blocked vectorization"
NOERRNO  := -O3 -march=native -fno-math-errno
COMMON   := -std=c++17 -Wall -Wextra -Wpedantic
SIZEFL   := -DFV_N=1048576 -DFV_NMID=65536 -DFV_NSMALL=8192 -DFV_REPS=100 -DFV_REPS_MID=400 -DFV_REPS_SMALL=4000 -DFV_TRIALS=5

LIB_NATIVE   := $(BUILD)/flashvec_O3_native.o
LIB_FFAST    := $(BUILD)/flashvec_O3_native_ffastmath.o
LIB_SHIPPED  := $(BUILD)/flashvec_O2_shipped.o
LIB_NOERRNO  := $(BUILD)/flashvec_O3_native_nomatherrno.o
REF          := $(BUILD)/ref_strict.o
BIN_NATIVE   := $(BUILD)/bench_O3_native
BIN_FFAST    := $(BUILD)/bench_O3_native_ffastmath
BIN_SHIPPED  := $(BUILD)/bench_O2_shipped
BIN_PERF     := $(BUILD)/bench_perf
BIN_NOERRNO  := $(BUILD)/bench_O3_native_nomatherrno
BIN_MEMBW    := $(BUILD)/membw
BINS         := $(BIN_NATIVE) $(BIN_FFAST) $(BIN_SHIPPED) $(BIN_NOERRNO)

.PHONY: all run checks table matrix vec asm perf clean

all: $(BINS) $(BIN_MEMBW)

$(BUILD):
	mkdir -p $(BUILD)

# Library under test: one object per flag set; tier 1 links them directly (across translation units, so they are not inlined)
$(LIB_NATIVE): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(NATIVE) -c $< -o $@

$(LIB_FFAST): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(FFAST) -c $< -o $@

$(LIB_SHIPPED): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(SHIPPED) -c $< -o $@

$(LIB_NOERRNO): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(NOERRNO) -c $< -o $@

# Accuracy reference implementation: strictly defined IEEE semantics, independent of the flags under test
$(REF): ref.cpp ref.h | $(BUILD)
	$(CXX) $(CPPFLAGS) -std=c++17 -O2 -fno-fast-math -ffp-contract=off -fno-tree-vectorize -c $< -o $@

$(BIN_NATIVE): bench.cpp ref.h $(REF) $(LIB_NATIVE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(NATIVE) $(SIZEFL) -DFV_CFG_NAME='"O3_native"' bench.cpp $(REF) $(LIB_NATIVE) -o $@

$(BIN_FFAST): bench.cpp ref.h $(REF) $(LIB_FFAST) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(FFAST) $(SIZEFL) -DFV_CFG_NAME='"O3_native_ffastmath"' bench.cpp $(REF) $(LIB_FFAST) -o $@

$(BIN_SHIPPED): bench.cpp ref.h $(REF) $(LIB_SHIPPED) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(SHIPPED) $(SIZEFL) -DFV_CFG_NAME='"O2_shipped"' bench.cpp $(REF) $(LIB_SHIPPED) -o $@

$(BIN_NOERRNO): bench.cpp ref.h $(REF) $(LIB_NOERRNO) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(NOERRNO) $(SIZEFL) -DFV_CFG_NAME='"O3_native_nomatherrno"' bench.cpp $(REF) $(LIB_NOERRNO) -o $@

# perf only: shortened repeat counts so perf stat can finish
$(BIN_PERF): bench.cpp ref.h $(REF) $(LIB_NATIVE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(NATIVE) -DFV_N=1048576 -DFV_NMID=65536 -DFV_NSMALL=8192 -DFV_REPS=10 -DFV_REPS_MID=50 -DFV_REPS_SMALL=200 -DFV_TRIALS=1 -DFV_CFG_NAME='"perf"' bench.cpp $(REF) $(LIB_NATIVE) -o $@

$(BIN_MEMBW): membw.c | $(BUILD)
	$(CC) -O3 -march=native -Wall -Wextra $< -o $@

run: $(BINS) $(BIN_MEMBW)
	@mkdir -p $(RAW)
	taskset -c $(CPU) $(BIN_MEMBW) | tee $(RAW)/membw.txt
	$(BIN_NATIVE)  all $(CPU) | tee $(RAW)/bench_O3_native.txt
	$(BIN_FFAST)   all $(CPU) | tee $(RAW)/bench_O3_native_ffastmath.txt
	$(BIN_SHIPPED) all $(CPU) | tee $(RAW)/bench_O2_shipped.txt
	$(BIN_NOERRNO) all $(CPU) | tee $(RAW)/bench_O3_native_nomatherrno.txt
	@echo "raw logs -> $(RAW)"

# checks only: repeat counts squeezed to 1, just to get DIFF/IDENT (accuracy and bit-identical) evidence quickly
$(BUILD)/bench_checks: bench.cpp ref.h $(REF) $(LIB_NATIVE) | $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(NATIVE) -DFV_N=1048576 -DFV_NMID=65536 -DFV_NSMALL=8192 -DFV_REPS=1 -DFV_REPS_MID=1 -DFV_REPS_SMALL=1 -DFV_TRIALS=1 -DFV_CFG_NAME='"checks"' bench.cpp $(REF) $(LIB_NATIVE) -o $@

checks: $(BUILD)/bench_checks
	@mkdir -p $(RAW)
	$(BUILD)/bench_checks all $(CPU) | grep -E '^DIFF|^IDENT' | tee $(RAW)/ident_check.txt

matrix:
	@$(AWK) -F'[ =]+' -f matrix.awk $(RAW)/bench_O3_native.txt $(RAW)/bench_O3_native_ffastmath.txt 	  $(RAW)/bench_O3_native_nomatherrno.txt $(RAW)/bench_O2_shipped.txt | tee $(RAW)/matrix.txt

table:
	@$(AWK) -f table.awk $(RAW)/bench_O3_native.txt $(RAW)/bench_O3_native_ffastmath.txt $(RAW)/bench_O2_shipped.txt | tee $(RAW)/table.md

vec:
	@mkdir -p $(BUILD)
	$(CXX) $(CPPFLAGS) $(COMMON) $(NATIVE) $(SIZEFL) -DFV_CFG_NAME='"vecreport"' \
	  -fopt-info-vec-optimized=$(BUILD)/vec_O3_native.txt \
	  -fopt-info-vec-missed=$(BUILD)/vecmiss_O3_native.txt -c bench.cpp -o $(BUILD)/bench_vecreport_native.o
	$(CXX) $(CPPFLAGS) $(COMMON) $(FFAST) $(SIZEFL) -DFV_CFG_NAME='"vecreport"' \
	  -fopt-info-vec-optimized=$(BUILD)/vec_O3_native_ffastmath.txt \
	  -fopt-info-vec-missed=$(BUILD)/vecmiss_O3_native_ffastmath.txt -c bench.cpp -o $(BUILD)/bench_vecreport_ffastmath.o
	@echo "== -O3 -march=native: vectorized t1_/t2_ loops =="
	@grep -E '^bench.cpp:[0-9]+:[0-9]+: optimized: loop vectorized' $(BUILD)/vec_O3_native.txt | head -40
	@echo "== -O3 -march=native -ffast-math: vectorized t1_/t2_ loops =="
	@grep -E '^bench.cpp:[0-9]+:[0-9]+: optimized: loop vectorized' $(BUILD)/vec_O3_native_ffastmath.txt | head -40

asm:
	@mkdir -p $(BUILD)
	objdump -d --no-show-raw-insn $(BIN_NATIVE) > $(BUILD)/disasm_O3_native.txt
	objdump -d --no-show-raw-insn $(BIN_SHIPPED) > $(BUILD)/disasm_O2_shipped.txt
	objdump -d --no-show-raw-insn $(BIN_FFAST) > $(BUILD)/disasm_O3_native_ffastmath.txt
	@echo "objdump -> $(BUILD)/disasm_*.txt"

perf: $(BIN_PERF)
	@echo "# kernel.perf_event_paranoid = $$(cat /proc/sys/kernel/perf_event_paranoid) (unprivileged perf stat is refused when >1)"
	@echo "# how to relax it (needs root): sudo sysctl -w kernel.perf_event_paranoid=1"
	taskset -c $(CPU) $(PERF) stat -e cycles,instructions,cache-references,cache-misses,branch-instructions,branch-misses \
	  $(BIN_PERF) all $(CPU)

clean:
	rm -rf $(BUILD)
