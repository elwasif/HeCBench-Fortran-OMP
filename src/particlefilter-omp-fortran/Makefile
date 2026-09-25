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

BLOCK_SIZE ?= 256
FFLAGS := $(EXTRA_FFLAGS) -Wall -cpp -DBLOCK_SIZE=$(BLOCK_SIZE)

ifeq ($(DEBUG),yes)
  FFLAGS += -g -DDEBUG
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
BENCHMARK_ARGS ?= -x 128 -y 128 -z 10 -np 400000

$(program): $(obj)
	$(FC) $(FFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) *.mod output.txt

run: $(program)
	$(LAUNCHER) ./$(program) $(BENCHMARK_ARGS)
