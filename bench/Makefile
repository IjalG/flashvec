# bench/Makefile -- flashvec 三档实现基准（新增文件；不改动主 Makefile 的既有目标）
#
# 用法：
#   make -C bench run    编译 + 跑全部档位，原始输出写入 ../docs/research/raw/
#   make -C bench table  从原始日志生成 Markdown 对比表（bench/table.awk）
#   make -C bench vec    抓 GCC 自动向量化报告（-fopt-info-vec）
#   make -C bench asm    反汇编归档（看调用开销 / ABI 传参）
#   make -C bench perf   perf stat 抓 IPC / cache-miss / branch-miss（受 perf_event_paranoid 限制）
#   make -C bench clean  清掉 bench/build

ROOT   := ..
BUILD  := build
RAW    := $(ROOT)/docs/research/raw
CXX    ?= g++
CC     ?= cc
AWK    ?= awk
PERF   ?= perf
CPU    ?= 2

CPPFLAGS += -I$(ROOT)/include

# 三组被测编译选项（-march=native 在本机展开为 AVX2 + FMA，见报告）
NATIVE   := -O3 -march=native
FFAST    := -O3 -march=native -ffast-math
SHIPPED  := -O2
# 只加一个严格语义安全的开关，用来隔离「到底是 fast-math 还是 math-errno 挡住了向量化」
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

# 被测库：三组编译选项各编一份，档1 直接链接它们（跨编译单元，不会被内联）
$(LIB_NATIVE): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(NATIVE) -c $< -o $@

$(LIB_FFAST): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(FFAST) -c $< -o $@

$(LIB_SHIPPED): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(SHIPPED) -c $< -o $@

$(LIB_NOERRNO): $(ROOT)/src/flashvec.cpp $(ROOT)/include/flashvec.h | $(BUILD)
	$(CXX) -I$(ROOT)/include $(COMMON) $(NOERRNO) -c $< -o $@

# 精度参考实现：固定严格 IEEE 语义，和被测编译选项无关
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

# perf 专用：缩短重复次数，让 perf stat 跑得完
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

# 校验段专用：重复次数压到 1，只为快速拿到 DIFF/IDENT（精度与逐位一致性）证据
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
	@echo "== -O3 -march=native：被向量化的 t1_/t2_ 循环 =="
	@grep -E '^bench.cpp:[0-9]+:[0-9]+: optimized: loop vectorized' $(BUILD)/vec_O3_native.txt | head -40
	@echo "== -O3 -march=native -ffast-math：被向量化的 t1_/t2_ 循环 =="
	@grep -E '^bench.cpp:[0-9]+:[0-9]+: optimized: loop vectorized' $(BUILD)/vec_O3_native_ffastmath.txt | head -40

asm:
	@mkdir -p $(BUILD)
	objdump -d --no-show-raw-insn $(BIN_NATIVE) > $(BUILD)/disasm_O3_native.txt
	objdump -d --no-show-raw-insn $(BIN_SHIPPED) > $(BUILD)/disasm_O2_shipped.txt
	objdump -d --no-show-raw-insn $(BIN_FFAST) > $(BUILD)/disasm_O3_native_ffastmath.txt
	@echo "objdump -> $(BUILD)/disasm_*.txt"

perf: $(BIN_PERF)
	@echo "# kernel.perf_event_paranoid = $$(cat /proc/sys/kernel/perf_event_paranoid)（>1 时非特权 perf stat 会被拒绝）"
	@echo "# 放开方式（需 root）：sudo sysctl -w kernel.perf_event_paranoid=1"
	taskset -c $(CPU) $(PERF) stat -e cycles,instructions,cache-references,cache-misses,branch-instructions,branch-misses \
	  $(BIN_PERF) all $(CPU)

clean:
	rm -rf $(BUILD)
