# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
CXX       = amdclang++
OPTIMIZE  ?= yes
DEBUG     ?= no
DEVICE    ?= gpu
ARCH      ?= gfx90a
LAUNCHER  ?=

program = main
f_source  = main.f90
cxx_source = nlll_rng.cpp
obj       = $(f_source:.f90=.o) $(cxx_source:.cpp=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall
CXXFLAGS := $(EXTRA_CXXFLAGS) -std=c++17 -Wall

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

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS) -lstdc++

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

run: $(program)
	$(LAUNCHER) ./$(program) 2048 1000 1000
