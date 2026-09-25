// SPDX-License-Identifier: CC0-1.0
#include <cstdint>
#include <random>

extern "C" void moe_sum_initialize_input(float* input, int64_t input_size, int topk)
{
  std::mt19937 gen(topk);
  std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

  for (int64_t i = 0; i < input_size; i++) {
    input[i] = dis(gen);
  }
}
