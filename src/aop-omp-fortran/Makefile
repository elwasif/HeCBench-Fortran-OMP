# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
CXX       = amdclang++
OPTIMIZE  ?= yes
DEBUG     ?= no
DEVICE    ?= gpu
ARCH      ?= gfx90a
LAUNCHER  ?=

program = main
source  = main.f90
cpp_source = aop_rng.cpp
obj     = $(source:.f90=.o) $(cpp_source:.cpp=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall -cpp

ifeq ($(DEBUG),yes)
  FFLAGS += -g -DDEBUG
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS += -O3
endif

ifeq ($(DEVICE),gpu)
  FFLAGS += -fopenmp --offload-arch=$(ARCH)
else
  FFLAGS += -fopenmp
endif

LDFLAGS =
BENCHMARK_ARGS ?= -runs 100

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -lstdc++ -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(EXTRA_CXXFLAGS) -O3 -c $< -o $@

clean:
	rm -rf $(program) $(obj) *.mod

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
