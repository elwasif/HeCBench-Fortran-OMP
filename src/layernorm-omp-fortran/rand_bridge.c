// SPDX-License-Identifier: CC0-1.0
#include <stdlib.h>

void hecbench_srand(unsigned int seed) {
  srand(seed);
}

float hecbench_rand_float(void) {
  return rand() / (float)RAND_MAX * 2.0f - 1.0f;
}
