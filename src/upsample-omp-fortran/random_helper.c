// SPDX-License-Identifier: CC0-1.0
#include <stdint.h>
#include <stdlib.h>

void fill_random_float(float *arr, int64_t n) {
  #pragma omp parallel for
  for (int64_t i = 0; i < n; ++i) {
    arr[i] = rand() / (float)RAND_MAX * 2.0f - 1.0f;
  }
}
