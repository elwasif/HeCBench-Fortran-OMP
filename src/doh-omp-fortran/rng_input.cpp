// SPDX-License-Identifier: CC0-1.0
#include <random>

extern "C" void doh_fill_input(float *input_img, int img_size)
{
  std::default_random_engine rng(123);
  std::normal_distribution<float> norm_dist(0.f, 1.f);

  for (int i = 0; i < img_size; i++) {
    input_img[i] = norm_dist(rng);
  }
}
