// SPDX-License-Identifier: CC0-1.0
#include <cstdint>
#include <random>

extern "C" void nlll_generate_inputs(
    std::int64_t input_size,
    std::int64_t weights_size,
    std::int64_t target_size,
    std::int64_t n_classes,
    float* input,
    float* weights,
    std::int32_t* target)
{
  std::default_random_engine g(123);
  std::uniform_real_distribution<float> d1(-1.f, 1.f);
  std::uniform_int_distribution<int> d2(0, static_cast<int>(n_classes) - 1);

  for (std::int64_t i = 0; i < input_size; ++i) {
    input[i] = d1(g);
  }

  for (std::int64_t i = 0; i < weights_size; ++i) {
    weights[i] = d1(g);
  }

  for (std::int64_t i = 0; i < target_size; ++i) {
    target[i] = static_cast<std::int32_t>(d2(g) + 1);
  }
}
