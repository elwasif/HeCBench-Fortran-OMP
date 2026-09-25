// SPDX-License-Identifier: CC0-1.0
#include <cstddef>
#include <random>

namespace {
std::mt19937 gen;
std::uniform_real_distribution<float> dis(0.0f, 4.0f);
}

extern "C" void michalewicz_mt19937_reset(int seed)
{
  gen.seed(static_cast<std::mt19937::result_type>(seed));
  dis.reset();
}

extern "C" void michalewicz_mt19937_fill(float *values, std::size_t count)
{
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = dis(gen);
  }
}
