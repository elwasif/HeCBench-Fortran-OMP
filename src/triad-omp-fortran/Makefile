# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
OPTIMIZE  ?= yes
DEBUG     ?= no
DEVICE    ?= gpu
ARCH      ?= gfx90a
LAUNCHER  ?=

program = main
source  = main.f90
obj     = $(source:.f90=.o)

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

BENCHMARK_ARGS ?= --passes 1

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *.mod

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
