// SPDX-License-Identifier: CC0-1.0
#include <cmath>
#include <random>

extern "C" void moe_initialize_inputs(
    float* gating_output,
    int* topk_indices,
    int* topk_indices_ref,
    int num_tokens,
    int num_experts,
    int topk)
{
  const int output_size = num_tokens * num_experts;

  std::mt19937 gen(19937);
  std::uniform_int_distribution<> distrib(-100000, 100000);

  for (int i = 0; i < output_size; i++) {
    gating_output[i] = distrib(gen);
  }

  for (int i = 0; i < topk; i++) {
    for (int j = 0; j < num_tokens; j++) {
      topk_indices_ref[i * num_tokens + j] =
          topk_indices[i * num_tokens + j] = std::abs(distrib(gen)) % num_experts;
    }
  }
}
