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
#   include/flashvec.h 的 static inline 入口里用函数级 target 属性 + 运行时
#   __builtin_cpu_supports("avx2") != 0 派发启用（实测 32 处派发全在头文件，src/*.cpp 里 0 处）。
CXXFLAGS ?= -std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CFLAGS   ?= -std=c11 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude

BUILD := build
LIB   := $(BUILD)/libflashvec.a
OBJ   := $(BUILD)/flashvec.o $(BUILD)/flashvec_batch.o $(BUILD)/flashvec_vec4.o

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

# ------------- 头文件 inline 形态（P1 重构）的验收目标（只新增；上面所有目标的命令一字未动） -------------
#   make test-header          严格编译档（C11 -pedantic-errors / C++17 -Wpedantic）；同一 TU 双重包含；
#                             编译期断言锁住「x86 + GCC/Clang 的 C TU 必须走 intrinsics 路径」；位级抽查。
#   make test-inline          默认（inline）形态 vs FV_USE_EXTERN（库内实现）形态的位级对拍：
#                             同一份 tests/test_parity.c 编五份（语言 × 形态 × 纯 C 回退），主输出流逐字节 cmp。
#   make test-fma             消费者带 -march=native -mfma -ffp-contract=fast 编译，结果仍必须位级一致
#                             —— 这是「有收缩风险的入口没被误判成无风险」的唯一保险。
#                             编译始终执行；运行需要 FMA 硬件（本机有），没有时只编译并打印 SKIP。
#   make test-inline-effects  消费者默认选项（-O2，不加 -ffp-contract=off）下的内联效果普查：
#                             第一段 8 个入口必须**完整内联**（n<4 的尾巴循环无 per-element call），
#                             第二段 5 个入口必须只保留**每批一次** TU 内调用。
#   make test-refactor        上面全部 + 导出符号清单 + 反汇编普查 + 链接证据。
.PHONY: test-header test-inline test-fma test-inline-effects test-stream test-refactor

HEADER_CFLAGS   := -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -O2
HEADER_CXXFLAGS := -std=c++17 -pedantic-errors -Wall -Wextra -Wpedantic -Iinclude -O2
PARITY_FLAGS    := -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic -Iinclude
FMA_FLAGS       := -O3 -Wall -Wextra -Wpedantic -march=native -mfma -ffp-contract=fast -Iinclude

# 「入口:探针函数」表。第一段 = 无 mul+add 收缩风险（要求完整内联）；
# 第二段 = 有收缩风险、由 pragma 守着（允许每批一次 TU 内调用，不允许 per-element call）。
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

# 对拍产物：语言（C/C++）× 形态（默认 inline / FV_USE_EXTERN）× 纯 C 回退（-DFV_I_NO_INTRINSICS）
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
	@echo "== 逐位一致入口的输出流（1486 例：n=0..17 全尾长 + 63；mask 四种模式；idx 重复/逆序/乱序/单元素/m=0；+ P2 的 7 个新入口）=="
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_c_extern.bin && echo "  PASS  C inline  == C extern（库内实现）"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_cpp_extern.bin && echo "  PASS  C inline  == C++ extern（库内实现）"
	@cmp $(BUILD)/parity_c_inline.bin $(BUILD)/parity_c_pure.bin && echo "  PASS  C inline  == C 纯 C 回退（-DFV_I_NO_INTRINSICS）"
	@echo "== dot_soa（契约里唯一不保证逐位一致的入口）：四个 AVX2 进程走同一条归约路径 =="
	@cmp $(BUILD)/parity_cpp_inline.dot.bin $(BUILD)/parity_cpp_extern.dot.bin && echo "  PASS  C++ inline == C++ extern"
	@cmp $(BUILD)/parity_c_extern.dot.bin $(BUILD)/parity_cpp_extern.dot.bin && echo "  PASS  C extern   == C++ extern"
	@if cmp -s $(BUILD)/parity_c_pure.dot.bin $(BUILD)/parity_cpp_extern.dot.bin; then \
	  echo "  NOTE  纯 C 回退的 dot 与 AVX2 相同（本次语料未触发归约顺序差异）"; \
	else echo "  NOTE  纯 C 回退的 dot 与 AVX2 不同 —— 契约允许（归约顺序不同），故不参与比较"; fi

test-fma: $(BUILD)/parity_cpp_fma $(BUILD)/parity_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/parity_cpp_fma $(BUILD)/parity_cpp_fma.bin $(BUILD)/parity_cpp_fma.dot.bin || exit 1; \
	  echo "== 消费者带 -march=native -mfma -ffp-contract=fast（最恶劣的真实选项）=="; \
	  cmp $(BUILD)/parity_cpp_fma.bin $(BUILD)/parity_cpp_extern.bin \
	    && echo "  PASS  12 个入口仍与库内实现位级一致（pragma 挡住了命令行带来的 FMA）"; \
	  cmp $(BUILD)/parity_cpp_fma.dot.bin $(BUILD)/parity_cpp_extern.dot.bin \
	    && echo "  PASS  dot_soa 也一致（两侧走同一条 AVX2 归约路径）"; \
	else echo "  SKIP: 本机无 FMA，只完成编译"; fi

$(BUILD)/consumer_probe_c.o: tests/consumer_probe.c include/flashvec.h | $(BUILD)
	$(CC) -std=c11 -O2 -Wall -Wextra -Wpedantic -Iinclude -c $< -o $@

$(BUILD)/consumer_probe_cpp.o: tests/consumer_probe.c include/flashvec.h | $(BUILD)
	$(CXX) -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude -c $< -o $@

test-inline-effects: $(BUILD)/consumer_probe_c.o $(BUILD)/consumer_probe_cpp.o
	@fail=0; \
	for obj in $(BUILD)/consumer_probe_c.o $(BUILD)/consumer_probe_cpp.o; do \
	  echo "== $$obj（-O2，未加 -ffp-contract=off = 消费方默认选项）=="; \
	  printf '  %-26s %-8s %-8s %-14s %s\n' "入口" "类别" "内联" "调 fv_i_ 次数" "判定"; \
	  for spec in $(SAFE_ENTRIES); do \
	    entry=$$(echo $$spec | cut -d: -f1); probe=$$(echo $$spec | cut -d: -f2); \
	    impl=$$(echo $$entry | sed -e 's/^fv_vec2_/fv_i_/' -e 's/^fv_vec3_/fv_i_v3_/' -e 's/^fv_vec4_/fv_i_v4_/'); \
	n=$$(nm --defined-only $$obj | awk '{print $$3}' | grep -cE "^$$impl([.]|$$)" || true); \
	    c=$$(objdump -d --no-show-raw-insn --disassemble="$$probe" $$obj | grep -cE ':(.*)(call|jmp).*<fv_i_' || true); \
	    if [ "$$n" = "0" ] && [ "$$c" = "0" ]; then v=PASS; else v=FAIL; fail=1; fi; \
	    inl=$$(if [ $$n = 0 ]; then echo 是; else echo 否; fi); \
	    printf '  %-26s %-8s %-8s %-14s %s\n' "$$entry" "第一段" "$$inl" "$$c" "$$v"; \
	  done; \
	  for spec in $(RISKY_ENTRIES); do \
	    entry=$$(echo $$spec | cut -d: -f1); probe=$$(echo $$spec | cut -d: -f2); \
	    impl=$$(echo $$entry | sed -e 's/^fv_vec2_/fv_i_/' -e 's/^fv_vec3_/fv_i_v3_/' -e 's/^fv_vec4_/fv_i_v4_/'); \
	n=$$(nm --defined-only $$obj | awk '{print $$3}' | grep -cE "^$$impl([.]|$$)" || true); \
	    c=$$(objdump -d --no-show-raw-insn --disassemble="$$probe" $$obj | grep -cE ':(.*)(call|jmp).*<fv_i_' || true); \
	    if [ "$$n" -ge 1 ] && [ "$$c" = "1" ]; then v=PASS; else v=FAIL; fail=1; fi; \
	    inl=$$(if [ $$n = 0 ]; then echo 是; else echo 否; fi); \
	    printf '  %-26s %-8s %-8s %-14s %s\n' "$$entry" "第二段" "$$inl" "$$c" "$$v"; \
	  done; \
	done; \
	if objdump -d $(BUILD)/consumer_probe_cpp.o | grep -qE 'jmp.*<fv_i_add_x4>'; then \
	  echo "  FAIL 仍有 jmp fv_i_add_x4"; fail=1; \
	else echo "  PASS 消费者 TU 里没有 jmp fv_i_add_x4（n=3 的尾巴已完整内联）"; fi; \
	echo "  -- probe_add_x4_3 的实际汇编（尾巴逐元素展开，无 call）--"; \
	objdump -d --no-show-raw-insn --disassemble=probe_add_x4_3 $(BUILD)/consumer_probe_cpp.o | tail -n +8 | head -10 | sed 's/^/    /'; \
	exit $$fail

# P2（内存流量入口）的专项验收：对齐 peel 的四条路径、_mask_inplace 的「不写」语义、sfence 可见性。
# 位级主对拍在 test_parity.c 里（make test-inline / test-fma 已经覆盖），这里只补它覆盖不到的形态。
$(BUILD)/test_stream.o: tests/test_stream.c include/flashvec.h | $(BUILD)
	$(CC) -std=c11 -pedantic-errors -Wall -Wextra -Wpedantic -O2 -Iinclude -c $< -o $@

$(BUILD)/test_stream: $(BUILD)/test_stream.o $(LIB)
	$(CXX) $^ -pthread -o $@

test-stream: $(BUILD)/test_stream
	./$(BUILD)/test_stream

test-refactor: test-header test-inline test-fma test-inline-effects test-stream \
               test-vec4 test-vec34-inline test-vec34-fma
	@echo "== 导出符号清单（nm -g --defined-only build/libflashvec.a | sort）=="
	@nm -g --defined-only $(LIB) | sort
	@echo "== 反汇编普查：这些模式在库内出现次数必须都是 0 =="
	@for pat in 'vfmadd[0-9]*[sp][ds]' 'vfmsub[0-9]*[sp][ds]' 'vfnmadd[0-9]*[sp][ds]' 'vrcpps' 'vrsqrtps' 'haddps' 'vmaskmovps'; do \
	  printf '  %-22s %s\n' "$$pat" "$$(objdump -d $(LIB) | grep -oE "$$pat" | wc -l)"; \
	done
	@echo "== 同样普查消费者产物（C 编译单元现在也带 AVX2，正对照应当非 0）=="
	@for obj in $(BUILD)/parity_c_inline.o $(BUILD)/parity_cpp_inline.o $(BUILD)/parity_c_pure.o; do \
	  printf '  %-30s vfmadd=%s  AVX(正对照)=%s\n' "$$obj" \
	    "$$(objdump -d $$obj | grep -oE 'vfmadd[0-9]*[sp][ds]' | wc -l)" \
	    "$$(objdump -d $$obj | grep -oE 'v(add|sub|mul|div|sqrt)ps' | wc -l)"; \
	done
	@echo "== 链接证据：默认形态的可执行文件里没有库内批量符号（inline 展开），extern 形态的从库里来 =="
	@printf '  parity_c_inline  全局批量符号 = %s\n' "$$(nm -g --defined-only $(BUILD)/parity_c_inline | grep -cE 'fv_vec2_.*_(x4|soa)')"
	@printf '  parity_c_extern  全局批量符号 = %s\n' "$$(nm -g --defined-only $(BUILD)/parity_c_extern | grep -cE 'fv_vec2_.*_(x4|soa)')"
# ============ P3：fv_vec4 / fv_vec3 SoA 批量（本段由 vec3/vec4 任务维护） ============
# 这一段是本轮唯一新增/改动的构建逻辑，说明三件事：
#   1) fv_vec4 的标量实现（src/flashvec_vec4.cpp）已并入 $(OBJ)：标量 API 与批量入口的标量
#      尾巴必须同源同选项，符号也必须在 libflashvec.a 里（nm 只增不减）；
#   2) 位级对拍用 tests/test_vec34_parity.c 编 5 种形态（语言 × 形态 × 纯 C 回退）+ FMA 消费者，
#      主流程逐字节 cmp；dot 流按契约分组（纯 C 回退的归约顺序不同，允许不同）；
#   3) 内联效应：P3 的 13 个入口已经加进上面既有的 SAFE_ENTRIES / RISKY_ENTRIES 表，
#      探针在 tests/consumer_probe.c 的 P3 小节里，所以 make test-inline-effects 直接覆盖它们。
.PHONY: test-vec4 test-vec34-inline test-vec34-fma

$(BUILD)/flashvec_vec4.o: src/flashvec_vec4.cpp include/flashvec.h | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# 标量 API 的位级参考测试。选项故意与库内实现对齐（-O3 -fno-math-errno -ffp-contract=off）：
# 参考表达式里的 sqrtf 于是与库里展开成同一条指令，mul+add 也不会收缩 —— 否则比的是编译器
# 差异，不是库的行为。
$(BUILD)/test_vec4.o: tests/test_vec4.c include/flashvec.h | $(BUILD)
	$(CC) $(CPPFLAGS) -std=c11 -O3 -fno-math-errno -ffp-contract=off -pedantic-errors -Wall -Wextra -Wpedantic -c $< -o $@

$(BUILD)/test_vec4: $(BUILD)/test_vec4.o $(LIB)
	$(CXX) $^ -o $@

test-vec4: $(BUILD)/test_vec4
	./$(BUILD)/test_vec4

# 对拍产物：语言（C/C++）× 形态（默认 inline / FV_USE_EXTERN）× 纯 C 回退 + FMA 消费者
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
	@echo "== 三维/四维批量入口：位级一致（n=0..17 全尾长 + 63 + 1024 + 8192；含 GUARD 毒化） =="
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_cpp_inline.bin && echo "  PASS  C inline  == C++ inline"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_c_extern.bin && echo "  PASS  C inline  == C extern（库内实现）"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_cpp_extern.bin && echo "  PASS  C inline  == C++ extern（库内实现）"
	@cmp $(BUILD)/vec34_c_inline.bin $(BUILD)/vec34_c_pure.bin && echo "  PASS  C inline  == C 纯 C 回退（-DFV_I_NO_INTRINSICS）"
	@echo "== dot_soa（契约里不保证逐位一致的两个入口）：四个 AVX2 进程走同一条归约路径 =="
	@cmp $(BUILD)/vec34_cpp_inline.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin && echo "  PASS  C++ inline == C++ extern"
	@cmp $(BUILD)/vec34_c_extern.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin && echo "  PASS  C extern   == C++ extern"
	@if cmp -s $(BUILD)/vec34_c_pure.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin; then \
	  echo "  NOTE  纯 C 回退的 dot 与 AVX2 相同（本次语料未触发归约顺序差异）"; \
	else echo "  NOTE  纯 C 回退的 dot 与 AVX2 不同 —— 契约允许（归约顺序不同），故不参与比较"; fi

test-vec34-fma: $(BUILD)/vec34_cpp_fma $(BUILD)/vec34_cpp_extern
	@if grep -qm1 '^flags.* fma' /proc/cpuinfo 2>/dev/null; then \
	  ./$(BUILD)/vec34_cpp_fma $(BUILD)/vec34_cpp_fma.bin $(BUILD)/vec34_cpp_fma.dot.bin || exit 1; \
	  echo "== 消费者带 -march=native -mfma -ffp-contract=fast（最恶劣的真实选项）=="; \
	  cmp $(BUILD)/vec34_cpp_fma.bin $(BUILD)/vec34_cpp_extern.bin \
	    && echo "  PASS  13 个新入口仍与库内实现位级一致（pragma 挡住了命令行带来的 FMA）"; \
	  cmp $(BUILD)/vec34_cpp_fma.dot.bin $(BUILD)/vec34_cpp_extern.dot.bin \
	    && echo "  PASS  dot_soa 也一致（两侧走同一条 AVX2 归约路径）"; \
	else echo "  SKIP: 本机无 FMA，只完成编译"; fi

# =====================================================================
# bench/ab 同进程 A/B 夹具的转发目标（只做转发，逻辑全在 bench/ab/）
#   make bench-ab          自检 + 基线快照（含 A/A 地板表）
#   make bench-ab-worktree candidate = 当前工作树
#   make bench-ab-negctl   负对照（候选被改错时必须硬失败并指名函数）
#   make bench-verify OP=<op> N=<规模>   同向证据五条判定（PASS/FAIL/INCONCLUSIVE）
# 这些目标只写 bench/ab/build/ 与 docs/research/raw/ab_*，不碰 build/rootbuild/。
# =====================================================================
.PHONY: bench-ab bench-ab-worktree bench-ab-negctl bench-verify

bench-ab:
	$(MAKE) -C bench/ab run

bench-ab-worktree:
	$(MAKE) -C bench/ab run-worktree

bench-ab-negctl:
	$(MAKE) -C bench/ab negctl

bench-verify:
	@if [ -z "$(OP)" ]; then echo "用法: make bench-verify OP=<op名> [N=<规模>]" >&2; exit 2; fi
	$(MAKE) -C bench/ab verify OP=$(OP) N=$(if $(N),$(N),1048576)
