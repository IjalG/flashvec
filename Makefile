CXX      ?= c++
CC       ?= cc
AR       ?= ar
# 编译选项（P0）：
#   -O3            解锁 -fno-math-errno 的收益。实测同一段含 sqrtf 的循环，
#                  -O2 -fno-math-errno 下 GCC 不向量化，-O3 下 "loop vectorized using 16 byte vectors"；
#                  报告里 len 类循环 x1.25 -> x8.16 的那个数字来自 -O3 配置。
#   -fno-math-errno 放弃「数学函数失败时设置 errno」。本库的 sqrt 参数恒 >= 0，
#                  errno 从不被设置；实测数值逐位不变，却让含 sqrtf 的循环可向量化。
#   -ffp-contract=off 禁止把 a + b*c 收缩成 FMA。逐元素 API 承诺逐位可复现
#                  （docs/optimization-report.md §4.7），而 FMA 收缩会让同一表达式在不同
#                  编译单元 / 不同目标下差 1 ulp（本机实测：-ffp-contract=fast 下出现 vfmadd）。
#   不加 -march=native / -mavx2 / -ffast-math：默认产物保持 x86-64 基线，AVX2 只在
#   src/flashvec_batch.cpp 里用函数级 target 属性 + 运行时派发启用。
CXXFLAGS ?= -std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CFLAGS   ?= -std=c11 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude

BUILD := build
LIB   := $(BUILD)/libflashvec.a
OBJ   := $(BUILD)/flashvec.o $(BUILD)/flashvec_batch.o

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

# ---------------- 基准（新增目标，不改动上面任何既有目标的行为） ----------------
# 基准源码与构建规则都在 bench/ 下；这里只是转发，方便在仓库根目录直接跑：
#   make bench        编译并跑完四组编译选项 × 三档实现 × 三种规模，原始输出写 docs/research/raw/
#   make bench-table  由原始日志重新生成 Markdown 对比表
#   make bench-vec    抓 GCC 自动向量化报告
#   make bench-asm    生成反汇编归档（bench/build/disasm_*.txt）
#   make bench-perf   perf stat（需要 kernel.perf_event_paranoid<=1 或有 CAP_PERFMON）
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

# ---------------- P0 / P1 测试（只追加，不改动上面 all/test/clean 的任何行为） ----------------
#   make test-p0     min/max/clamp/madd、AABB、FV_EPS 量纲修正的边界测试
#   make test-batch  批量入口的位级一致性（memcmp）、掩码语义、索引变体；含一个 C 语言的调用测试
#   make test-all    上面两个 + 既有 make test（30 条断言 + C 兼容性测试）
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

test-all: test test-p0 test-batch
