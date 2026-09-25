// SPDX-License-Identifier: CC0-1.0
#include <random>

namespace {
std::default_random_engine rng;
std::uniform_real_distribution<double> dist(0.0, 1.0);
}

extern "C" void phmm_rng_seed(int seed) {
  rng.seed(seed);
  dist.reset();
}

extern "C" double phmm_rng_next() {
  return dist(rng);
}
