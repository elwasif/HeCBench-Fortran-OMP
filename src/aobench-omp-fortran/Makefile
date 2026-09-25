# SPDX-License-Identifier: CC0-1.0
#===============================================================================
# User Options
#===============================================================================

FC        = amdflang
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

obj = $(source:.f90=.o)

#===============================================================================
# Sets Flags
#===============================================================================

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
BENCHMARK_ARGS ?= 1

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) ao.ppm

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
