// SPDX-License-Identifier: CC0-1.0
#include <random>

extern "C" void dense_embedding_random_fill(float *dense, float *input,
                                             int dense_size, int input_size) {
  std::default_random_engine g(123);
  std::uniform_real_distribution<float> distr(-1.0f, 1.0f);

  for (int i = 0; i < dense_size; ++i) {
    dense[i] = distr(g);
  }

  for (int i = 0; i < input_size; ++i) {
    input[i] = distr(g);
  }
}
