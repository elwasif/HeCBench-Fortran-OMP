# SPDX-License-Identifier: CC0-1.0
#===============================================================================
# User Options
#===============================================================================

FC        = amdflang
CXX       = amdclang++
OPTIMIZE  = yes
DEBUG     = no
DEVICE    = gpu
ARCH      = gfx90a
LAUNCHER  =

#===============================================================================
# Program name & source code list
#===============================================================================

program = main
source = main.f90
cxx_source = mt19937_helper.cpp
obj = $(source:.f90=.o) $(cxx_source:.cpp=.o)

#===============================================================================
# Sets Flags
#===============================================================================

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

LDFLAGS = -lstdc++

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *.mod

run: $(program)
	$(LAUNCHER) ./$(program) 1024 1
