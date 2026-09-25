// SPDX-License-Identifier: CC0-1.0
#include <random>
#include <sys/resource.h>

namespace {
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-1.0, 1.0);
}

extern "C" void dslash_rng_seed_random_device()
{
  gen.seed(rd());
}

extern "C" double dslash_rng_uniform_real()
{
  return dis(gen);
}

extern "C" double dslash_maxrss_mb()
{
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    return usage.ru_maxrss / 1024.0;
  }
  return 0.0;
}
