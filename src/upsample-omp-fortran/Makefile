# SPDX-License-Identifier: CC0-1.0
FC        = amdflang
CC        = amdclang
OPTIMIZE  = yes
DEBUG     = no
DEVICE    = gpu
ARCH      = gfx90a
LAUNCHER  =

program = main
fsource = main.f90
csource = random_helper.c
obj     = $(fsource:.f90=.o) $(csource:.c=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall
CFLAGS := $(EXTRA_CFLAGS) -Wall -fopenmp

ifeq ($(DEBUG),yes)
  FFLAGS += -g
  CFLAGS += -g
endif

ifeq ($(OPTIMIZE),yes)
  FFLAGS += -O3
  CFLAGS += -O3
endif

ifeq ($(DEVICE),gpu)
  FFLAGS += -fopenmp --offload-arch=$(ARCH)
else
  FFLAGS += -fopenmp
endif

LDFLAGS =
BENCHMARK_ARGS ?= 2 4 8 8 1

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
