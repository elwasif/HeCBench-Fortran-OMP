# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
CC        = amdclang
OPTIMIZE  = yes
DEBUG     = no
DEVICE    = gpu
ARCH      = gfx90a
LAUNCHER  =

program = main
source  = main.f90
csource = rand_fill.c
obj     = $(source:.f90=.o) $(csource:.c=.o)

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

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) -O3 -c $< -o $@

clean:
	rm -rf $(program) $(obj)

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
