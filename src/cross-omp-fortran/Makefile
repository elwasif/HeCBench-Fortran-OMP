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
cpp_source = cross_rng.cpp
obj     = $(source:.f90=.o) $(cpp_source:.cpp=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall
CXXFLAGS := $(EXTRA_CXXFLAGS)

ifeq ($(DEBUG),yes)
  FFLAGS += -g
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

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS) -lstdc++

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

run: $(program)
	$(LAUNCHER) ./$(program) 10000000 100
