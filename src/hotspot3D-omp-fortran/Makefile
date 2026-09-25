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
cpp_source = hotspot3d_reference.cpp
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
  CXXFLAGS += -fopenmp --offload-arch=$(ARCH)
else
  FFLAGS += -fopenmp
  CXXFLAGS += -fopenmp
endif

LDFLAGS =

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) output.out

run: $(program)
	$(LAUNCHER) ./$(program) 4 2 2 fixtures/power_4x2 fixtures/temp_4x2 output.out
