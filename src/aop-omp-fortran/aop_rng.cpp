// SPDX-License-Identifier: CC0-1.0
#include <cstddef>
#include <random>

namespace {
std::default_random_engine rng;
std::normal_distribution<double> norm_dist(0.0, 1.0);
}

extern "C" void aop_reset_rng()
{
  rng = std::default_random_engine();
  norm_dist.reset();
}

extern "C" void aop_fill_samples(double *samples, std::size_t count)
{
  for (std::size_t i = 0; i < count; ++i)
    samples[i] = norm_dist(rng);
}
