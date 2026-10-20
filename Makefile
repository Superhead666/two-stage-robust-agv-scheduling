CXX ?= c++
PYTHON ?= python3
GUROBI_ROOT ?= /Library/gurobi1302/macos_universal2
GUROBI_LIB ?= gurobi130
CPLEX_ROOT ?= /Applications/CPLEX_Studio2211
BUILD ?= build
OPTFLAGS ?= -O3 -DNDEBUG
LDFLAGS ?=
CPPFLAGS += -Isrc/common -Isrc/hig -Isrc/lbbd -Isrc/ccg -I$(GUROBI_ROOT)/include
CXXFLAGS += -std=c++17 $(OPTFLAGS) -Wall -Wextra -Wpedantic -MMD -MP
LDLIBS := -L$(GUROBI_ROOT)/lib -Wl,-rpath,$(GUROBI_ROOT)/lib -lgurobi_c++ -l$(GUROBI_LIB)

ifeq ($(shell uname -s),Darwin)
  ifeq ($(shell uname -m),arm64)
    CPLEX_PLATFORM ?= arm64_osx
  else
    CPLEX_PLATFORM ?= x86-64_osx
  endif
  CP_SYSTEM_LIBS := -lpthread -framework CoreFoundation -framework IOKit
else
  CPLEX_PLATFORM ?= x86-64_linux
  CP_SYSTEM_LIBS := -lpthread -ldl -lm
endif
CPOFLAGS := -isystem $(CPLEX_ROOT)/cpoptimizer/include -isystem $(CPLEX_ROOT)/concert/include -fno-strict-aliasing
CPOLIBS := -L$(CPLEX_ROOT)/cpoptimizer/lib/$(CPLEX_PLATFORM)/static_pic -L$(CPLEX_ROOT)/cplex/lib/$(CPLEX_PLATFORM)/static_pic -L$(CPLEX_ROOT)/concert/lib/$(CPLEX_PLATFORM)/static_pic -lcp -lcplex -lconcert $(CP_SYSTEM_LIBS)

CORE := $(BUILD)/common/core_dp.o $(BUILD)/common/robust.o
NATIVE := $(BUILD)/common/core_native.o
BLOCK := $(BUILD)/common/block_mip.o
EVAL := $(BUILD)/hig/evaluator.o $(BLOCK)
HIG := $(BUILD)/hig/long_search.o $(BUILD)/hig/tight.o $(BUILD)/hig/proposals.o $(EVAL) $(CORE)
CCG := $(BUILD)/ccg/master.o $(BUILD)/ccg/oracle.o $(BUILD)/ccg/solve.o $(BLOCK) $(CORE)
CPCCG := $(BUILD)/ccg/master_with_cp.o $(BUILD)/ccg/cp_master.o $(BUILD)/ccg/oracle.o $(BUILD)/ccg/solve.o $(BLOCK) $(CORE)
PROGRAMS := agv_hig agv_lbbd agv_ccg agv_initialize agv_schedule_metrics run_case agv_model_compare
VALIDATION := common_validation hig_validation hig_alignment lbbd_packing lbbd_structural ccg_validation model_comparison

.PHONY: all cp check check-cp check-data manifest clean
all: $(addprefix $(BUILD)/,$(PROGRAMS))
cp: $(BUILD)/agv_ccg_cp

$(BUILD)/agv_hig: $(BUILD)/hig/main.o $(HIG)
	$(CXX) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/agv_lbbd: $(BUILD)/lbbd/branch_check.o $(EVAL) $(CORE)
	$(CXX) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/agv_ccg: $(BUILD)/ccg/main.o $(CCG)
	$(CXX) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/agv_ccg_cp: $(BUILD)/ccg/main.o $(CPCCG)
	$(CXX) $^ $(LDLIBS) $(CPOLIBS) $(LDFLAGS) -o $@
$(BUILD)/agv_initialize: $(BUILD)/runner/initialize.o $(NATIVE)
	$(CXX) $^ $(LDFLAGS) -o $@
$(BUILD)/agv_schedule_metrics: $(BUILD)/runner/schedule_metrics.o $(NATIVE)
	$(CXX) $^ $(LDFLAGS) -o $@
$(BUILD)/run_case: $(BUILD)/runner/run_case.o
	$(CXX) $^ $(LDFLAGS) -o $@
$(BUILD)/agv_model_compare: $(BUILD)/runner/model_compare.o $(HIG)
	$(CXX) $^ $(LDLIBS) $(LDFLAGS) -o $@

$(BUILD)/common/core_dp.o: src/common/core.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Ddeterministic=deterministic_dp -c $< -o $@
$(BUILD)/common/core_native.o: src/common/core.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
$(BUILD)/ccg/master_with_cp.o: src/ccg/master.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DAGV_WITH_CP -c $< -o $@
$(BUILD)/ccg/cp_master.o: src/ccg/cp_master.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CPOFLAGS) -DAGV_WITH_CP -c $< -o $@
$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/tests/common_validation: tests/common/validation.cpp tests/common/fleet.cpp $(NATIVE) $(BUILD)/common/robust.o
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/hig_validation: tests/hig/check.cpp $(HIG)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/hig_alignment: tests/hig/alignment.cpp $(HIG)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/lbbd_packing: tests/lbbd/check.cpp $(NATIVE)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
$(BUILD)/tests/lbbd_structural: tests/lbbd/check_structural.cpp $(NATIVE)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/ccg_validation: tests/ccg/main.cpp tests/ccg/validation.cpp $(CCG)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/ccg_cp_validation: tests/ccg/main.cpp tests/ccg/validation.cpp $(CPCCG)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CPOFLAGS) -DAGV_WITH_CP $^ $(LDLIBS) $(CPOLIBS) $(LDFLAGS) -o $@
$(BUILD)/tests/model_comparison: tests/model_comparison.cpp $(HIG)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDLIBS) $(LDFLAGS) -o $@

check: all $(addprefix $(BUILD)/tests/,$(VALIDATION))
	LC_ALL=C $(BUILD)/tests/common_validation validate
	LC_ALL=C $(BUILD)/tests/hig_validation
	LC_ALL=C $(BUILD)/tests/hig_alignment
	LC_ALL=C $(BUILD)/tests/lbbd_packing tests/data/lbbd
	LC_ALL=C $(BUILD)/tests/lbbd_structural
	LC_ALL=C $(BUILD)/tests/ccg_validation
	LC_ALL=C $(BUILD)/tests/model_comparison
	LC_ALL=C $(PYTHON) tests/integration_smoke.py --build $(BUILD)
	LC_ALL=C $(PYTHON) tests/experiment_driver.py --build $(BUILD)
	LC_ALL=C $(PYTHON) tests/sensitivity_driver.py --build $(BUILD)
	LC_ALL=C $(PYTHON) tests/overlap_verifier.py --build $(BUILD)
	LC_ALL=C $(PYTHON) tests/model_comparison_driver.py --build $(BUILD)
check-cp: all cp $(BUILD)/tests/ccg_cp_validation
	LC_ALL=C $(BUILD)/tests/ccg_cp_validation
	LC_ALL=C $(PYTHON) tests/integration_smoke.py --build $(BUILD) --cp
check-data:
	$(PYTHON) scripts/verify_archive.py
manifest:
	$(PYTHON) scripts/verify_archive.py --write-manifest
clean:
	@$(PYTHON) -c "from pathlib import Path; import shutil; root=Path.cwd().resolve(); path=Path('$(BUILD)').resolve(); assert path.parent==root and path.name.startswith('build'); shutil.rmtree(path,ignore_errors=True)"

-include $(shell test ! -d $(BUILD) || find $(BUILD) -name '*.d')
