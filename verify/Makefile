# verify/Makefile —— 验收者 F 的独立构建（只写 verify/ 下的文件，不碰 include/ src/ tests/ bench/ 根 Makefile）
#
#   make parity   —— 位级一致性验证：库在 -O2 与 -O3 -march=native 下各构建一份并各跑一遍 + ASan/UBSan
#   make perf     —— 基准：-O3 -march=native 与默认 -O2 两档，原始输出写 build/
#   make check    —— 纯 C11 头文件/ABI 检查 + 仓库既有测试（我自己重编译）+ 结构检查
#   make extra    —— 实现者 E 自己的测试套件（我只编译运行，不修改）
#   make all      —— 以上全部

CXX ?= c++
CC  ?= cc
INC := -I../include
WARN := -Wall -Wextra
SRC1 := ../src/flashvec.cpp
SRC2 := ../src/flashvec_batch.cpp
B := build

# 校验方（parity.cpp / ref_strict.cpp）一律 -ffp-contract=off：保证我自己的参考值与期望值
# 不会被编译器收缩成 FMA，避免误报。
PARITY_FLAGS := $(WARN) $(INC) -ffp-contract=off
LIB_O2  := -std=c++17 -O2
LIB_O3N := -std=c++17 -O3 -march=native
OBJ_O2  := $(B)/flashvec_O2.o $(B)/flashvec_batch_O2.o
OBJ_O3N := $(B)/flashvec_O3n.o $(B)/flashvec_batch_O3n.o
# -O3 不带 -march=native：与根 Makefile 现在的发布默认档（-O3 -fno-math-errno -ffp-contract=off）对齐
DEF_O3  := -std=c++17 -O3 -fno-math-errno -ffp-contract=off
OBJ_O3  := $(B)/flashvec_O3.o $(B)/flashvec_batch_O3.o

.PHONY: all parity perf check extra clean run-parity run-perf

all: parity perf check extra

$(B):
	mkdir -p $(B)

# ---------------- 被验证库的两档构建（两个 TU） ----------------
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

# ---------------- 位级一致性验证 ----------------
$(B)/parity_O2_plain: parity.cpp ref_strict.h ../include/flashvec.h $(OBJ_O2) $(B)/ref_strict.o | $(B)
	$(CXX) -std=c++17 -O2 $(PARITY_FLAGS) parity.cpp $(B)/ref_strict.o $(OBJ_O2) -lm -o $@

$(B)/parity_O3_native: parity.cpp ref_strict.h ../include/flashvec.h $(OBJ_O3N) $(B)/ref_strict.o | $(B)
	$(CXX) -std=c++17 -O3 -march=native $(PARITY_FLAGS) parity.cpp $(B)/ref_strict.o $(OBJ_O3N) -lm -o $@

$(B)/parity_asan: parity.cpp ref_strict.h ../include/flashvec.h $(SRC1) $(SRC2) ref_strict.cpp | $(B)
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(PARITY_FLAGS) \
		parity.cpp ref_strict.cpp $(SRC1) $(SRC2) -lm -o $@

parity: $(B)/parity_O2_plain $(B)/parity_O3_native $(B)/parity_asan
	@echo '================ parity (-O2，默认构建档：无 -march=native) ================'
	@$(B)/parity_O2_plain > $(B)/parity_O2_plain.txt 2>&1; rc=$$?; cat $(B)/parity_O2_plain.txt; echo "exit=$$rc"
	@echo '================ parity (-O3 -march=native) ================'
	@$(B)/parity_O3_native > $(B)/parity_O3_native.txt 2>&1; rc=$$?; cat $(B)/parity_O3_native.txt; echo "exit=$$rc"
	@echo '================ parity (ASan+UBSan) ================'
	@$(B)/parity_asan > $(B)/parity_asan.txt 2>&1; rc=$$?; echo "exit=$$rc"; grep -E 'FAIL|ERROR|SUMMARY|runtime error' $(B)/parity_asan.txt | head -10; tail -2 $(B)/parity_asan.txt

run-parity: parity

# ---------------- 基准 ----------------
$(B)/perf_O3_native: perf.cpp ../include/flashvec.h $(OBJ_O3N) | $(B)
	$(CXX) -std=c++17 -O3 -march=native -ffp-contract=off $(WARN) $(INC) perf.cpp $(OBJ_O3N) -lm -o $@

$(B)/perf_O2_plain: perf.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CXX) -std=c++17 -O2 -ffp-contract=off $(WARN) $(INC) perf.cpp $(OBJ_O2) -lm -o $@

# 发布默认档实测（不带 -march=native）——用来回答「默认构建还有没有 AVX2 收益」
$(B)/perf_O3_plain: perf.cpp ../include/flashvec.h $(OBJ_O3) | $(B)
	$(CXX) $(DEF_O3) $(WARN) $(INC) perf.cpp $(OBJ_O3) -lm -o $@

# 固定核 + 无缓冲输出：避免迁移与「崩溃丢失缓冲」两种伪影
perf: $(B)/perf_O3_native $(B)/perf_O3_plain $(B)/perf_O2_plain
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (-O3 -march=native) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O3_native > $(B)/perf_O3_native.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O3_native.txt
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (发布默认档 -O3，无 -march=native) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O3_plain > $(B)/perf_O3_plain.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O3_plain.txt
	@PIN=""; command -v taskset >/dev/null 2>&1 && PIN="taskset -c 2"; \
	echo "================ perf (-O2 无 -march=native，旧默认档) ================"; \
	$$PIN stdbuf -o0 -e0 $(B)/perf_O2_plain > $(B)/perf_O2_plain.txt 2>&1; echo "exit=$$?"; cat $(B)/perf_O2_plain.txt

run-perf: perf

# ---------------- 结构 / 纯 C / 既有测试 ----------------
$(B)/c_api_check_O2: c_api_check.c ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CC) -std=c11 -O3 -pedantic-errors $(WARN) $(INC) -ffp-contract=off c_api_check.c $(OBJ_O2) -lm -o $@

$(B)/test_flashvec_F: ../tests/test_flashvec.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@

$(B)/test_c_api_F.o: ../tests/test_c_api.c ../include/flashvec.h | $(B)
	$(CC) -std=c11 -O2 -pedantic-errors $(WARN) $(INC) -c $< -o $@

$(B)/test_c_api_F: $(B)/test_c_api_F.o $(OBJ_O2)
	$(CXX) $^ -lm -o $@

check: $(B)/c_api_check_O2 $(B)/test_flashvec_F $(B)/test_c_api_F
	@echo '================ 纯 C11 头文件 + ABI 检查 (-pedantic-errors) ================'
	@$(B)/c_api_check_O2 > $(B)/c_api_check.txt 2>&1; rc=$$?; cat $(B)/c_api_check.txt; echo "exit=$$rc"
	@echo '================ 仓库既有测试（我自己重编译，不动 build/） ================'
	@$(B)/test_flashvec_F > $(B)/test_flashvec_F.txt 2>&1; rc=$$?; cat $(B)/test_flashvec_F.txt; echo "exit=$$rc"
	@$(B)/test_c_api_F > $(B)/test_c_api_F.txt 2>&1; rc=$$?; cat $(B)/test_c_api_F.txt; echo "exit=$$rc"
	@./checks.sh

# ---------------- E 自己的测试套件（只编译运行，不修改） ----------------
$(B)/test_batch_F: ../tests/test_batch.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	-$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@
$(B)/test_batch_c_F.o: ../tests/test_batch_c.c ../include/flashvec.h | $(B)
	-$(CC) -std=c11 -O2 -pedantic-errors $(WARN) $(INC) -c $< -o $@
$(B)/test_batch_c_F: $(B)/test_batch_c_F.o $(OBJ_O2)
	-$(CXX) $^ -lm -o $@
$(B)/test_p0_F: ../tests/test_p0.cpp ../include/flashvec.h $(OBJ_O2) | $(B)
	-$(CXX) -std=c++17 -O2 $(WARN) $(INC) $< $(OBJ_O2) -lm -o $@
extra: $(B)/test_batch_F $(B)/test_batch_c_F $(B)/test_p0_F
	@echo '================ 实现者 E 自己的测试（我重编译后独立跑） ================'
	@for t in test_batch_F test_batch_c_F test_p0_F; do if [ -x $(B)/$$t ]; then $(B)/$$t > $(B)/$$t.txt 2>&1; echo "$$t exit=$$?"; tail -6 $(B)/$$t.txt; fi; done

clean:
	rm -rf $(B)
