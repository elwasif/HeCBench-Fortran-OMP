// SPDX-License-Identifier: CC0-1.0
#include <cstdint>
#include <random>

namespace {
std::mt19937 generator;
std::uniform_int_distribution<int> distribute(0, 255);
}

extern "C" void bs_seed_rng(std::int64_t seed) {
  generator.seed(static_cast<std::mt19937::result_type>(seed));
}

extern "C" int bs_next_rand_byte() {
  return distribute(generator);
}
