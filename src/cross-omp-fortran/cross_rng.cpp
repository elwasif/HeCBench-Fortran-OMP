// SPDX-License-Identifier: CC0-1.0
#include <random>

extern "C" void fill_cross_inputs_float(int num_elems, float* a, float* b)
{
  std::default_random_engine g(123);
  std::uniform_real_distribution<float> distr(-2.f, 2.f);
  for (int i = 0; i < num_elems; ++i) {
    a[i] = distr(g);
    b[i] = distr(g);
  }
}

extern "C" void fill_cross_inputs_double(int num_elems, double* a, double* b)
{
  std::default_random_engine g(123);
  std::uniform_real_distribution<double> distr(-2.f, 2.f);
  for (int i = 0; i < num_elems; ++i) {
    a[i] = distr(g);
    b[i] = distr(g);
  }
}
