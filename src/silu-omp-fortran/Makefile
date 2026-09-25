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
csources = random_helper.c
obj     = $(source:.f90=.o) $(csources:.c=.o)

FFLAGS := $(EXTRA_FFLAGS) -Wall
CFLAGS := $(EXTRA_CFLAGS) -Wall

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
  CFLAGS += -fopenmp
else
  FFLAGS += -fopenmp
  CFLAGS += -fopenmp
endif

LDFLAGS =
BENCHMARK_ARGS ?= 1 1 16 16 1

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
