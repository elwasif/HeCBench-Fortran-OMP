# HeCBench-OMP-Fortran

This repository contains a curated subset of OpenMP benchmarks derived from the
HeCBench benchmark suite. It keeps the original C/C++ OpenMP benchmark sources
alongside corresponding Fortran OpenMP offloading ports.

This repository accompanies the paper "HeCBench-OMP-Fortran: An OpenMP
Offloading Benchmark Suite for Fortran Codes", submitted to IWOMP 2026 - The
22nd International Workshop on OpenMP.

## Contents

- `src/*-omp/`: original C/C++ OpenMP benchmark implementations from HeCBench.
- `src/*-omp-fortran/`: Fortran OpenMP offloading ports for the same benchmark
  set.
- `STRUCTURE_REPAIR_CHANGES.md`: benchmark-by-benchmark summary of source and
  build-file changes between the initial port state and the structure-repair
  state.

The repository contains 188 benchmark pairs. Each retained benchmark has both
an original C/C++ OpenMP implementation and a matching Fortran OpenMP
offloading implementation.  The Fortran ports were produced by ChatGPT-5.5.

## Licensing

`SPDX-License-Identifier: BSD-3-Clause AND CC0-1.0`

- The original benchmarks (`src/*-omp/`) are licensed under the [BSD 3-Clause License](LICENSE.bsd).
- The AI-generated translations (`src/*-omp-fortran/`) are in the public domain under the [Creative Commons CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](LICENSE.cc0).

## Provenance

The original benchmark sources come from
[HeCBench](https://github.com/ORNL/HeCBench), a heterogeneous computing
benchmark suite containing CUDA, HIP, SYCL, and OpenMP implementations. This
repository was produced by filtering the development repository to retain the
selected C/C++ OpenMP and Fortran OpenMP benchmark pairs.

## Tags

- `initial-port`: the initial retained Fortran OpenMP port state.
- `structure-repair`: the retained state after source-level structure repairs
  to improve fidelity with the original C/C++ OpenMP benchmark implementations.
