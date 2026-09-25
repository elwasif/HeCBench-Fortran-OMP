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
cpp_source = input_generator.cpp
obj     = $(source:.f90=.o)
cpp_obj = $(cpp_source:.cpp=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall
CXXFLAGS := $(EXTRA_CXXFLAGS) -Wall

ifeq ($(DEBUG),yes)
  FFLAGS += -g
  CXXFLAGS += -g
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS += -O3
  CXXFLAGS += -O3
endif

ifeq ($(DEVICE),gpu)
  FFLAGS += -fopenmp --offload-arch=$(ARCH)
else
  FFLAGS += -fopenmp
endif

LDFLAGS =
BENCHMARK_ARGS ?= 1 4 1

$(program): $(obj) $(cpp_obj)
	$(FC) $(FFLAGS) $(obj) $(cpp_obj) -o $@ $(LDFLAGS) -lstdc++

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) $(cpp_obj)

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
