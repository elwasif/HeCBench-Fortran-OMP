// SPDX-License-Identifier: CC0-1.0
#include <stdint.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void silu_make_random_float(float *array, int64_t n)
{
#pragma omp parallel for
  for (int64_t i = 0; i < n; i++) {
    array[i] = 2.0f * ((float)rand() / (float)RAND_MAX) - 1.0f;
  }
}
