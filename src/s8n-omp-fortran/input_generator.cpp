// SPDX-License-Identifier: CC0-1.0
#include <random>

extern "C" void s8n_fill_input_cpp(int *values, int input_size) {
  std::default_random_engine g(123);
  std::uniform_int_distribution<> distr(-256, 255);
  for (int i = 0; i < input_size; ++i) {
    values[i] = distr(g);
  }
}
