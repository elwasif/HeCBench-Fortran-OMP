// SPDX-License-Identifier: CC0-1.0
#include <stdint.h>
#include <string.h>

#include <random>

using half = _Float16;

static std::mt19937 engine;

extern "C" void relu_rng_init(int seed)
{
  engine.seed(seed);
}

extern "C" void relu_fill_gradient_feature(int count, uint16_t* gradient, uint16_t* feature)
{
  std::uniform_real_distribution<float> real_dist(-1.f, 1.f);

  for (int i = 0; i < count; i++) {
    half feature_h = half(real_dist(engine));
    half gradient_h = half(1.f);
    memcpy(&feature[i], &feature_h, sizeof(feature[i]));
    memcpy(&gradient[i], &gradient_h, sizeof(gradient[i]));
  }
}

extern "C" void relu_fill_int_inputs(int count, int32_t* input)
{
  std::uniform_int_distribution<unsigned char> int_dist(0, 255);

  for (int i = 0; i < count; i++) {
    uint32_t packed = (uint32_t)int_dist(engine)       |
                      (uint32_t)int_dist(engine) <<  8 |
                      (uint32_t)int_dist(engine) << 16 |
                      (uint32_t)int_dist(engine) << 24;
    memcpy(&input[i], &packed, sizeof(input[i]));
  }
}
