# bench/ab/Makefile -- in-process A/B interleaving fixture (self-contained: it touches
# neither the top-level Makefile nor any existing file under bench/)
#
# Usage (always call with -C; every path derives from this file's location, so it works from
# any directory):
#   make -C bench/ab run            self-check + baseline snapshot: candidate == git HEAD
#                                   snapshot
#                                   -> docs/research/raw/ab_baseline_head.{txt,tsv}
#                                      docs/research/raw/ab_baseline_head_matrix.tsv
#                                      docs/research/raw/ab_baseline_head_smallN.tsv
#                                      docs/research/raw/ab_aa_floor.tsv (A/A floor table)
#   make -C bench/ab run-worktree   candidate = current worktree (the real comparison)
#   make -C bench/ab verify OP=fv_vec2_add_x4 N=1048576
#                                   same-direction evidence: AB/BA x two independent builds +
#                                   A/A floor + five verdict criteria
#   make -C bench/ab aa-floor       rerun the self-check alone and write the A/A floor table
#   make -C bench/ab negctl         negative control: break the candidate on purpose, the
#                                   fixture must fail hard and name the function
#   make -C bench/ab snapshot       re-take the git HEAD snapshot (automatic when HEAD moves)
#   make -C bench/ab clean          remove only bench/ab/build (never anyone else's build/)
#
# Important: baseline always compiles the sources in snapshot/ (git show HEAD:...), never the
# worktree src/, because other agents are editing it. The candidate defaults to the snapshot
# too (== the self-check); use run-worktree explicitly to measure the worktree.
#
# Both library forms are supported (see the notes at the top of ab_variant.cpp):
#   CAND_WITH_SRC=1  batch entry points live in src/flashvec_batch.cpp (the HEAD form)
#   CAND_WITH_SRC=0  batch entry points are static inline in the header (post-refactor form,
#                    no .cpp compiled)
# Switching form only needs make ... CAND_HDR=... CAND_WITH_SRC=0; no source changes.
#
# Adding an op: one line in the ab_ops.h registry + one thunk in ab_variant.cpp. README §9.

ABDIR   := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ROOT    := $(abspath $(ABDIR)/../..)
BUILD   := $(ABDIR)/build
SNAP    := $(ABDIR)/snapshot
RAW     := $(ROOT)/docs/research/raw

CXX      ?= g++
# **literally identical** to the top-level Makefile defaults: baseline and candidate must
# use the same options, otherwise the ratio means nothing.
CXXFLAGS ?= -std=c++17 -O3 -fno-math-errno -ffp-contract=off -Wall -Wextra -Wpedantic

# Pin to one core: cpu0-7 share one CCX (L3 = 4 MB); crossing a CCX switches L3 and
# invalidates the cache-hierarchy annotation.
CPU         ?= 6
SIZES       ?= 4096,16384,65536,131072,262144,1048576,4194304,16777216
SMALL       ?= 4,8,16,32,64
TRIALS_MIN  ?= 5
TRIALS_MAX  ?= 10   # matches the accepted baseline snapshot; 12/15 only average longer
TARGET_MS   ?= 8
SMALL_MS    ?= 20

# Experimental area: newly landed entry points (_stream / _mask_inplace ...). Default 0 keeps
# the gate at the historical 182 items.
AB_OPS_EXTRA ?= 0

# ---- the three candidate knobs (default = HEAD snapshot, i.e. self-check mode) ----
CAND_INC      ?= $(SNAP)/include
CAND_HDR      ?= $(CAND_INC)/flashvec.h
CAND_SRC      ?= $(SNAP)/src/flashvec_batch.cpp
CAND_WITH_SRC ?= 1
CAND_LABEL    ?= head
# Repo-relative on purpose: the value is recorded in the console banner and in every TSV, so an
# absolute path would bake this checkout's location (which may be non-ASCII) into the artifacts.
CAND_DESC     ?= git HEAD snapshot ($(patsubst $(ROOT)/%,%,$(CAND_SRC))) [self-check mode: same source as baseline]
LABEL         ?= baseline_head
SELFCHECK     ?= 1
SELF_FLAG     := $(if $(filter 1,$(SELFCHECK)),--selfcheck)
# The A/A floor table is only meaningful for builds where both sides are the same source
FLOOR_FLAG    := $(if $(filter 1,$(SELFCHECK)),--floor-out=$(RAW)/ab_aa_floor.tsv,)

GIT_REV := $(shell git -C $(ROOT) rev-parse --short HEAD 2>/dev/null || echo unknown)

# Does this side header have the new entry points (recursive variable: re-evaluated on every
# use so a freshly refreshed snapshot is never read stale). Without them the thunk degrades to
# a plain store (see ab_variant.cpp) and nothing fails.
A_HAVE_STREAM   = $(shell grep -qs 'fv_vec2_add_x4_stream' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_STREAM   = $(shell grep -qs 'fv_vec2_add_x4_stream' $(CAND_HDR) && echo 1 || echo 0)
A_HAVE_INPLACE  = $(shell grep -qs 'fv_vec2_add_x4_mask_inplace' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_INPLACE  = $(shell grep -qs 'fv_vec2_add_x4_mask_inplace' $(CAND_HDR) && echo 1 || echo 0)
# vec3/vec4 entry points (12 B / 16 B elements): absent from HEAD, candidate only
A_HAVE_WIDE     = $(shell grep -qs 'fv_vec4_add_x2' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_WIDE     = $(shell grep -qs 'fv_vec4_add_x2' $(CAND_HDR) && echo 1 || echo 0)
# matrix transforms: the driver also needs to know whether the BASELINE has them, because that
# decides whether they are A/B ops or INTRA-only (see AB_MAT4_IN_BASELINE in ab_ops.h).
A_HAVE_MAT4     = $(shell grep -qs 'fv_mat4_transform_points' $(SNAP)/include/flashvec.h && echo 1 || echo 0)
B_HAVE_MAT4     = $(shell grep -qs 'fv_mat4_transform_points' $(CAND_HDR) && echo 1 || echo 0)

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

# ---- baseline source snapshot: git show HEAD; leave mtime alone when nothing changed
# ---- (so nothing is needlessly rebuilt) ----
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

# ---- symbol rename table: generated from **this side header** (a new entry point is
# ---- isolated as soon as it lands; nothing to maintain by hand) ----
$(BUILD)/rename_a.h: $(SNAP)/include/flashvec.h $(ABDIR)/mkrename.sh | $(BUILD)
	@sh $(ABDIR)/mkrename.sh $(SNAP)/include/flashvec.h ab_a_ $@
	@echo "[rename] A <- $@ ($$(grep -c '^#define' $@) symbols)"

$(BUILD)/rename_b.h: $(CAND_HDR) $(ABDIR)/mkrename.sh | $(BUILD)
	@sh $(ABDIR)/mkrename.sh $(CAND_HDR) ab_b_ $@
	@echo "[rename] B <- $@ ($$(grep -c '^#define' $@) symbols)"

# ---- variant A: baseline (HEAD snapshot, library form) ----
$(BUILD)/var_a.o: $(ABDIR)/ab_variant.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                  $(SNAP)/include/flashvec.h $(SNAP)/src/flashvec_batch.cpp \
                  $(BUILD)/rename_a.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -include $(BUILD)/rename_a.h \
	  -DAB_HDR='"$(SNAP)/include/flashvec.h"' \
	  -DAB_SRC='"$(SNAP)/src/flashvec_batch.cpp"' -DAB_WITH_SRC=1 \
	  -DAB_HAVE_STREAM=$(A_HAVE_STREAM) -DAB_HAVE_INPLACE=$(A_HAVE_INPLACE) \
	  -DAB_HAVE_WIDE=$(A_HAVE_WIDE) \
	  -DAB_HAVE_MAT4=$(A_HAVE_MAT4) -DAB_MAT4_IN_BASELINE=$(A_HAVE_MAT4) \
	  -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_LOOP_NAME=ab_a_loop -DAB_ABI_NAME=ab_a_abi -DAB_FLAGS_NAME=ab_a_opflags \
	  -c $(ABDIR)/ab_variant.cpp -o $@

# ---- variant B: candidate (form and source tree both come from variables) ----
$(BUILD)/candidate.cfg: FORCE | $(BUILD)
	@printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' '$(CAND_HDR)' '$(CAND_SRC)' \
	  '$(CAND_WITH_SRC)' '$(CXXFLAGS)' '$(AB_OPS_EXTRA)' '$(A_HAVE_STREAM)' \
	  '$(B_HAVE_STREAM)' '$(B_HAVE_INPLACE)' '$(A_HAVE_MAT4)' '$(B_HAVE_MAT4)' > $@.tmp
	@cmp -s $@.tmp $@ && rm -f $@.tmp || mv $@.tmp $@

# $(CAND_SRC) / $(CAND_HDR) must be **plain** prerequisites: otherwise, when a few source
# lines change inside the same tree (that is exactly what the negative control does, and so
# does the user editing the worktree), make cannot tell and silently reuses the stale object.
$(BUILD)/var_b.o: $(ABDIR)/ab_variant.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                  $(CAND_HDR) $(CAND_SRC) $(BUILD)/candidate.cfg $(BUILD)/rename_b.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(CAND_INC) -include $(BUILD)/rename_b.h \
	  -DAB_HDR='"$(CAND_HDR)"' \
	  -DAB_SRC='"$(CAND_SRC)"' -DAB_WITH_SRC=$(CAND_WITH_SRC) \
	  -DAB_HAVE_STREAM=$(B_HAVE_STREAM) -DAB_HAVE_INPLACE=$(B_HAVE_INPLACE) \
	  -DAB_HAVE_WIDE=$(B_HAVE_WIDE) \
	  -DAB_HAVE_MAT4=$(B_HAVE_MAT4) -DAB_MAT4_IN_BASELINE=$(A_HAVE_MAT4) \
	  -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_LOOP_NAME=ab_b_loop -DAB_ABI_NAME=ab_b_abi -DAB_FLAGS_NAME=ab_b_opflags \
	  -c $(ABDIR)/ab_variant.cpp -o $@

# ---- scalar anchor ("before the refactor", element-wise scalar API): library side and loop
# ---- side must live in two TUs ----
$(BUILD)/scalar_lib.o: $(ABDIR)/ab_scalar_lib.cpp $(SNAP)/src/flashvec.cpp \
                       $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -c $< -o $@

$(BUILD)/scalar_loop.o: $(ABDIR)/ab_scalar_loop.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                        $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(SNAP)/include -c $< -o $@

$(BUILD)/ab_bench.o: $(ABDIR)/ab_bench.cpp $(ABDIR)/ab_abi.h $(ABDIR)/ab_ops.h \
                     $(SNAP)/include/flashvec.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -I$(ABDIR) -DAB_OPS_EXTRA=$(AB_OPS_EXTRA) \
	  -DAB_MAT4_IN_BASELINE=$(A_HAVE_MAT4) \
	  -DAB_SNAP_HAS_VEC4=$(A_HAVE_WIDE) \
	  -DAB_LABEL='"$(LABEL)"' -DAB_CAND_DESC='"$(CAND_DESC)"' -DAB_GIT_REV='"$(GIT_REV)"' \
	  -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "[link] dual-implementation symbols: $$(nm -g --defined-only $@ | grep -cE ' T ab_[abs]_(loop|abi|opflags)')  loop/abi/opflags exports"

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

# Rerun for the A/A floor table alone (self-check only, floor table only; A/B are both the
# HEAD snapshot)
aa-floor: $(BIN) | $(RAW)
	@taskset -c $(CPU) $(BIN) --label=aa_floor --aa-each-round --sizes=$(SIZES) --no-small --no-pipe \
	  --trials-min=$(TRIALS_MIN) --trials-max=$(TRIALS_MAX) --target-ms=$(TARGET_MS) \
	  --floor-out=$(RAW)/ab_aa_floor.tsv > /dev/null
	@echo "A/A floor -> $(RAW)/ab_aa_floor.tsv"

# ---- candidate = the current worktree: **snapshot first, compile second** ----
# Other agents are editing include/ and src/. Compiling the live files fails now and then or,
# worse, compiles a half-edited header without showing it. So cp -a into
# build/worktree_snapshot/ first, record git status and sha256 alongside, and compile from
# that snapshot. The candidate's "identity" is therefore traceable.
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
	@echo "[worktree-snapshot] source fingerprints:" && sed 's/^/  /' $(WT_SNAP)/SHA256.txt
	@echo "[worktree-snapshot] worktree changes:" && sed 's/^/  /' $(WT_SNAP)/GIT_STATUS.txt
	@echo "[worktree-snapshot] form CAND_WITH_SRC=$$(cat $(WT_SNAP)/WITH_SRC) (1=library 0=header-inline)"

run-worktree: worktree-snapshot
	@$(MAKE) --no-print-directory -f $(ABDIR)/Makefile run \
	  CAND_INC=$(WT_SNAP)/include CAND_HDR=$(WT_SNAP)/include/flashvec.h \
	  CAND_SRC=$(WT_SNAP)/src/flashvec_batch.cpp \
	  CAND_WITH_SRC=$$(cat $(WT_SNAP)/WITH_SRC) \
	  CAND_LABEL=worktree LABEL=vs_worktree SELFCHECK=0 \
	  'CAND_DESC=current worktree snapshot (build/worktree_snapshot; fingerprints and git status in that dir)'

# ---- same-direction evidence: OP / N come from the command line ----
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
	@sed -n '2,32p' $(ABDIR)/Makefile

FORCE:
