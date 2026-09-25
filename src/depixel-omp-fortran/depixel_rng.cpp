// SPDX-License-Identifier: CC0-1.0
#include <random>

extern "C" {

void depixel_seed_rng(int seed);
float depixel_next_random();

}

namespace {
std::mt19937 gen;
std::uniform_real_distribution<float> dis(0.f, 0.4f);
}

void depixel_seed_rng(int seed)
{
  gen.seed(seed);
}

float depixel_next_random()
{
  return dis(gen);
}
