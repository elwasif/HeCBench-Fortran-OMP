// SPDX-License-Identifier: CC0-1.0
#include <cstddef>
#include <random>

extern "C" void fill_logprob_logits(float* logits, std::size_t logits_size)
{
  std::default_random_engine g(123);
  std::uniform_real_distribution<float> distr(-6.f, 6.f);

  for (std::size_t i = 0; i < logits_size; ++i) {
    logits[i] = distr(g);
  }
}
