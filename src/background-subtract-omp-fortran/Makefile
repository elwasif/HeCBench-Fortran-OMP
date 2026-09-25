# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
CXX       = amdclang++
OPTIMIZE  = yes
DEBUG     = no
DEVICE    = gpu
ARCH      = gfx90a
LAUNCHER  =

program = main
source  = main.f90
cxxsource = mt19937_rng.cpp
obj     = $(source:.f90=.o) $(cxxsource:.cpp=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall

ifeq ($(DEBUG),yes)
  FFLAGS += -g
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
BENCHMARK_ARGS ?= 64 64 0 3

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS) -lstdc++

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -O3 -c $< -o $@

clean:
	rm -rf $(program) $(obj)

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
