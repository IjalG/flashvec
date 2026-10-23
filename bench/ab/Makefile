# bench/ab/Makefile —— 同进程 A/B 交错夹具（自包含；不碰顶层 Makefile 与 bench/ 既有文件）
#
# 用法（用 -C 调用；路径全部按本文件位置推导，可从任何目录调用）：
#   make -C bench/ab run            自检 + 基线快照：candidate == git HEAD 快照
#                                   -> docs/research/raw/ab_baseline_head.{txt,tsv}
#                                      docs/research/raw/ab_baseline_head_matrix.tsv
#                                      docs/research/raw/ab_baseline_head_smallN.tsv
#                                      docs/research/raw/ab_aa_floor.tsv（A/A 地板表）
#   make -C bench/ab run-worktree   candidate = 当前工作树（真正比改动）
#   make -C bench/ab verify OP=fv_vec2_add_x4 N=1048576
#                                   同向证据：AB/BA × 两套独立构建 + A/A 地板 + 五条判定
#   make -C bench/ab aa-floor       单独重跑自检并落盘 A/A 地板表
#   make -C bench/ab negctl         负对照：把 candidate 故意改错，必须硬失败并指名函数
#   make -C bench/ab snapshot       重新抓取 git HEAD 快照（HEAD 变了会自动重抓）
#   make -C bench/ab clean          只删 bench/ab/build（不碰任何别人的 build/）
#
# 重要：baseline 永远编译 snapshot/ 里的源码（git show HEAD:...），绝不编译工作树 src/，
# 因为别的代理正在改它。candidate 默认也用快照（== 自检），跑工作树要显式 run-worktree。
#
# 两种库形态都支持（见 ab_variant.cpp 顶部的说明）：
#   CAND_WITH_SRC=1  批量入口在 src/flashvec_batch.cpp 里（HEAD 的形态）
#   CAND_WITH_SRC=0  批量入口是头文件里的 static inline（重构后的形态，不编 .cpp）
# 换形态只需 make ... CAND_HDR=… CAND_WITH_SRC=0，不需要改任何源码。
#
# 加一个 op：改 ab_ops.h 的注册表一行 + ab_variant.cpp 一个 thunk。见 README §9。

ABDIR   := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT    := $(abspath $(ABDIR)/../..)
BUILD   := $(ABDIR)/build
SNAP    := $(ABDIR)/snapshot
RAW     := $(ROOT)/docs/research/raw

CXX      ?= g++
# 与顶层 Makefile 默认**逐字一致**：baseline 与 candidate 必须同选项，否则比值没意义。
CXXFLAGS ?= -std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic

# 绑一个核：cpu0-7 同属一个 CCX（L3 = 4 MB）；跨 CCX 会换 L3，破坏缓存层级标注。
CPU         ?= 6
SIZES       ?= 4096,16384,65536,131072,262144,1048576,4194304,16777216
SMALL       ?= 4,8,16,32,64
TRIALS_MIN  ?= 5
TRIALS_MAX  ?= 10   # 与已验收的那次基线快照一致；12/15 只是把噪声平得更久，不减偏置
TARGET_MS   ?= 8
SMALL_MS    ?= 20

# 实验区：新落地的入口（_stream / _mask_inplace …）。默认 0 = 门禁仍是历史上那份 182 项。
AB_OPS_EXTRA ?= 0

# ---- candidate 的三个旋钮（默认 = HEAD 快照，自检模式）----
CAND_INC      ?= $(SNAP)/include
CAND_HDR      ?= $(CAND_INC)/flashvec.h
CAND_SRC      ?= $(SNAP)/src/flashvec_batch.cpp
CAND_WITH_SRC ?= 1
CAND_LABEL    ?= head
CAND_DESC     ?= git HEAD 快照 ($(CAND_SRC)) [自检模式：与 baseline 同一份源码]
LABEL         ?= baseline_head
SELFCHECK     ?= 1
SELF_FLAG     := $(if $(filter 1,$(SELFCHECK)),--selfcheck)
# A/A 地板表只在"两侧同一份源码"的构建里才有意义
FLOOR_FLAG    := $(if $(filter 1,$(SELFCHECK)),--floor-out=$(RAW)/ab_aa_floor.tsv,)

GIT_REV := $(shell git -C $(ROOT) rev-parse --short HEAD 2>/dev/null || echo unknown)

# 本侧头文件里有没有新入口（递归变量：每次使用时重算，避免 snapshot 刚刷新就读到旧值）。
# 没有就退化成普通 store（见 ab_variant.cpp 的 thunk），不报错。
A_HAVE_STREAM   = $(shell grep -qs 'fv_vec2_add_x4_stream' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_STREAM   = $(shell grep -qs 'fv_vec2_add_x4_stream' $(CAND_HDR) && echo 1 || echo 0)
A_HAVE_INPLACE  = $(shell grep -qs 'fv_vec2_add_x4_mask_inplace' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_INPLACE  = $(shell grep -qs 'fv_vec2_add_x4_mask_inplace' $(CAND_HDR) && echo 1 || echo 0)
# vec3/vec4 入口（元素 12 B / 16 B）：HEAD 里没有，只有 candidate 有
A_HAVE_WIDE     = $(shell grep -qs 'fv_vec4_add_x2' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_WIDE     = $(shell grep -qs 'fv_vec4_add_x2' $(CAND_HDR) && echo 1 || echo 0)

BIN  := $(BUILD)/ab_bench
OBJS := $(BUILD)/ab_bench.o $(BUILD)/var_a.o $(BUILD)/var_b.o $(BUILD)/scalar_lib.o \
        $(BUILD)/scalar_loop.o

.PHONY: all run run-worktree build negctl verify aa-floor snapshot clean help FORCE
all: $(BIN)
build: $(BIN)

$(BUILD):
	@mkdir -p $@

$(SNAP)/src $(SNAP)/include:
	@mkdir -p $@

# ---- baseline 源码快照：git show HEAD，内容不变就不动 mtime（避免无谓重编译）----
snapshot: | $(SNAP)/src $(SNAP)/include
	@git -C $(ROOT) show HEAD:src/flashvec_batch.cpp > $(SNAP)/src/.tmp.cpp
	@cmp -s $(SNAP)/src/.tmp.cpp $(SNAP)/src/flashvec_batch.cpp && rm -f $(SNAP)/src/.tmp.cpp || mv $(SNAP)/src/.tmp.cpp $(SNAP)/src/flashvec_batch.cpp
	@git -C $(ROOT) show HEAD:src/flashvec.cpp > $(SNAP)/src/.tmp2.cpp
	@cmp -s $(SNAP)/src/.tmp2.cpp $(SNAP)/src/flashvec.cpp && rm -f $(SNAP)/src/.tmp2.cpp || mv $(SNAP)/src/.tmp2.cpp $(SNAP)/src/flashvec.cpp
	@git -C $(ROOT) show HEAD:include/flashvec.h > $(SNAP)/include/.tmp.h
	@cmp -s $(SNAP)/include/.tmp.h $(SNAP)/include/flashvec.h && rm -f $(SNAP)/include/.tmp.h || mv $(SNAP)/include/.tmp.h $(SNAP)/include/flashvec.h
	@git -C $(ROOT) rev-parse HEAD > $(SNAP)/.tmp.rev
	@cmp -s $(SNAP)/.tmp.rev $(SNAP)/HEAD.txt && rm -f $(SNAP)/.tmp.rev || mv $(SNAP)/.tmp.rev $(SNAP)/HEAD.txt
	@echo "[snapshot] baseline <- git HEAD $$(cat $(SNAP)/HEAD.txt)"

$(SNAP)/src/flashvec_batch.cpp $(SNAP)/src/flashvec.cpp $(SNAP)/include/flashvec.h: | snapshot
	@:

# ---- 符号改名表：从**本侧头文件**自动生成（新入口落地即自动隔离，不用手工维护）----
$(BUILD)/rename_a.h: $(SNAP)/include/flashvec.h $(ABDIR)/mkrename.sh | $(BUILD)
	@sh $(ABDIR)/mkrename.sh $(SNAP)/include/flashvec.h ab_a_ $@
	@echo "[rename] A <- $@ ($$(grep -c '^#define' $@) 个符号)"

$(BUILD)/rename_b.h: $(CAND_HDR) $(ABDIR)/mkrename.sh | $(BUILD)
	@sh $(ABDIR)/mkrename.sh $(CAND_HDR) ab_b_ $@
	@echo "[rename] B <- $@ ($$(grep -c '^#define' $@) 个符号)"

# ---- 变体 A：baseline（HEAD 快照，库内定义形态）----
$(BUILD)/var_a.o: $(ABDIR)/ab_variant.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                  $(SNAP)/include/flashvec.h $(SNAP)/src/flashvec_batch.cpp \
                  $(BUILD)/rename_a.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -include $(BUILD)/rename_a.h \
	  -DAB_HDR='"$(SNAP)/include/flashvec.h"' \
	  -DAB_SRC='"$(SNAP)/src/flashvec_batch.cpp"' -DAB_WITH_SRC=1 \
	  -DAB_HAVE_STREAM=$(A_HAVE_STREAM) -DAB_HAVE_INPLACE=$(A_HAVE_INPLACE) \
	  -DAB_HAVE_WIDE=$(A_HAVE_WIDE) \
	  -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_LOOP_NAME=ab_a_loop -DAB_ABI_NAME=ab_a_abi -DAB_FLAGS_NAME=ab_a_opflags \
	  -c $(ABDIR)/ab_variant.cpp -o $@

# ---- 变体 B：candidate（形态与源码树都由变量决定）----
$(BUILD)/candidate.cfg: FORCE | $(BUILD)
	@printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' '$(CAND_HDR)' '$(CAND_SRC)' \
	  '$(CAND_WITH_SRC)' '$(CXXFLAGS)' '$(AB_OPS_EXTRA)' '$(A_HAVE_STREAM)' \
	  '$(B_HAVE_STREAM)' '$(B_HAVE_INPLACE)' > $@.tmp
	@cmp -s $@.tmp $@ && rm -f $@.tmp || mv $@.tmp $@

# $(CAND_SRC) / $(CAND_HDR) 必须作为**普通**依赖：否则在同一棵树里改了几行源码
# （负对照就是这么干的，用户改工作树也一样）时 make 认不出变化，会默默用旧对象。
$(BUILD)/var_b.o: $(ABDIR)/ab_variant.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                  $(CAND_HDR) $(CAND_SRC) $(BUILD)/candidate.cfg $(BUILD)/rename_b.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(CAND_INC) -include $(BUILD)/rename_b.h \
	  -DAB_HDR='"$(CAND_HDR)"' \
	  -DAB_SRC='"$(CAND_SRC)"' -DAB_WITH_SRC=$(CAND_WITH_SRC) \
	  -DAB_HAVE_STREAM=$(B_HAVE_STREAM) -DAB_HAVE_INPLACE=$(B_HAVE_INPLACE) \
	  -DAB_HAVE_WIDE=$(B_HAVE_WIDE) \
	  -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_LOOP_NAME=ab_b_loop -DAB_ABI_NAME=ab_b_abi -DAB_FLAGS_NAME=ab_b_opflags \
	  -c $(ABDIR)/ab_variant.cpp -o $@

# ---- 标量锚点（"重构前"，逐元素调标量 API）：库侧与循环侧必须分两个 TU ----
$(BUILD)/scalar_lib.o: $(ABDIR)/ab_scalar_lib.cpp $(SNAP)/src/flashvec.cpp \
                       $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -c $< -o $@

$(BUILD)/scalar_loop.o: $(ABDIR)/ab_scalar_loop.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                        $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -c $< -o $@

$(BUILD)/ab_bench.o: $(ABDIR)/ab_bench.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                     $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(ABDIR) -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_SNAP_HAS_VEC4=$(A_HAVE_WIDE) \
	  -DAB_LABEL='"$(LABEL)"' -DAB_CAND_DESC='"$(CAND_DESC)"' -DAB_GIT_REV='"$(GIT_REV)"' \
	  -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "[link] 双实现符号: $$(nm -g --defined-only $@ | grep -cE ' T ab_[abs]_(loop|abi|opflags)')  个 loop/abi/opflags 导出"

$(RAW):
	@mkdir -p $@

run: $(BIN) | $(RAW)
	@taskset -c $(CPU) $(BIN) $(SELF_FLAG) --aa-each-round --label=$(LABEL) --sizes=$(SIZES) --small=$(SMALL) \
	  --trials-min=$(TRIALS_MIN) --trials-max=$(TRIALS_MAX) \
	  --target-ms=$(TARGET_MS) --small-target-ms=$(SMALL_MS) \
	  --tsv-out=$(RAW)/ab_$(LABEL).tsv \
	  --matrix-out=$(RAW)/ab_$(LABEL)_matrix.tsv \
	  --small-out=$(RAW)/ab_$(LABEL)_smallN.tsv $(FLOOR_FLAG) \
	  > $(RAW)/ab_$(LABEL).txt 2>&1 ; \
	rc=$$? ; cat $(RAW)/ab_$(LABEL).txt ; exit $$rc

# A/A 地板表单独重跑（只跑自检、只落地板表；A/B 都是 HEAD 快照）
aa-floor: $(BIN) | $(RAW)
	@taskset -c $(CPU) $(BIN) --label=aa_floor --aa-each-round --sizes=$(SIZES) --no-small --no-pipe \
	  --trials-min=$(TRIALS_MIN) --trials-max=$(TRIALS_MAX) --target-ms=$(TARGET_MS) \
	  --floor-out=$(RAW)/ab_aa_floor.tsv > /dev/null
	@echo "A/A floor -> $(RAW)/ab_aa_floor.tsv"

# ---- candidate = 当前工作树：**先拍快照再编** ----
# 别的代理正在改 include/ 与 src/。直接编活文件会偶发编译失败，或者更糟：编到半改状态的
# 头文件而看不出来。所以先 cp -a 一份到 build/worktree_snapshot/，把 git status 与 sha256
# 一起记下来，再从这份快照编。candidate 的"身份"因此是可追溯的。
WT_SNAP := $(BUILD)/worktree_snapshot

.PHONY: worktree-snapshot
worktree-snapshot: | $(BUILD)
	@rm -rf $(WT_SNAP)
	@mkdir -p $(WT_SNAP)/include $(WT_SNAP)/src
	@cp -a $(ROOT)/include/flashvec.h $(WT_SNAP)/include/
	@cp -a $(ROOT)/src/flashvec_batch.cpp $(ROOT)/src/flashvec.cpp $(WT_SNAP)/src/
	@( cd $(ROOT) && git status --porcelain -- include src ) > $(WT_SNAP)/GIT_STATUS.txt 2>&1 || true
	@( cd $(WT_SNAP) && sha256sum include/flashvec.h src/flashvec_batch.cpp src/flashvec.cpp ) \
	  > $(WT_SNAP)/SHA256.txt
	@if grep -q 'FV_BUILD' $(WT_SNAP)/src/flashvec_batch.cpp; then echo 0 > $(WT_SNAP)/WITH_SRC; \
	 else echo 1 > $(WT_SNAP)/WITH_SRC; fi
	@echo "[worktree-snapshot] 源文件指纹：" && sed 's/^/  /' $(WT_SNAP)/SHA256.txt
	@echo "[worktree-snapshot] 工作树改动：" && sed 's/^/  /' $(WT_SNAP)/GIT_STATUS.txt
	@echo "[worktree-snapshot] 形态判定 CAND_WITH_SRC=$$(cat $(WT_SNAP)/WITH_SRC) （1=库内定义 0=头文件内联）"

run-worktree: worktree-snapshot
	@$(MAKE) --no-print-directory -f $(ABDIR)/Makefile run \
	  CAND_INC=$(WT_SNAP)/include CAND_HDR=$(WT_SNAP)/include/flashvec.h \
	  CAND_SRC=$(WT_SNAP)/src/flashvec_batch.cpp \
	  CAND_WITH_SRC=$$(cat $(WT_SNAP)/WITH_SRC) \
	  CAND_LABEL=worktree LABEL=vs_worktree SELFCHECK=0 \
	  'CAND_DESC=当前工作树快照 (build/worktree_snapshot，指纹与 git status 见该目录)'

# ---- 同向证据：OP / N 由命令行给 ----
OP ?=
N  ?= 1048576
verify: $(BIN)
	@env OP='$(OP)' N='$(N)' CPU='$(CPU)' TRIALS_MIN='$(TRIALS_MIN)' \
	  TRIALS_MAX='$(TRIALS_MAX)' TARGET_MS='$(TARGET_MS)' \
	  CHECK_SIZES='$(CHECK_SIZES)' EXPECT_REGRESS='$(EXPECT_REGRESS)' \
	  MIN_GAIN='$(MIN_GAIN)' MAX_REGRESS='$(MAX_REGRESS)' \
	  sh $(ABDIR)/verify.sh

negctl: $(BIN)
	@sh $(ABDIR)/negctl.sh

clean:
	rm -rf $(BUILD)

help:
	@sed -n '2,35p' $(ABDIR)/Makefile

FORCE:
