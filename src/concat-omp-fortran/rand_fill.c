// SPDX-License-Identifier: CC0-1.0
#include <stdlib.h>

int concat_rand_mod(int modulus) {
  return rand() % modulus;
}

void concat_fill_rand_float(float *values, long long n) {
  for (long long i = 0; i < n; ++i) {
    values[i] = (float)(rand() % n);
  }
}
