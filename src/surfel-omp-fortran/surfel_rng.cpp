// SPDX-License-Identifier: CC0-1.0
#include <cmath>
#include <random>

namespace {
constexpr int COL_P_X = 0;
constexpr int COL_P_Y = 1;
constexpr int COL_P_Z = 2;
constexpr int COL_N_X = 3;
constexpr int COL_N_Y = 4;
constexpr int COL_N_Z = 5;
constexpr int COL_RSq = 6;
constexpr int COL_DIM = 7;
}

extern "C" void surfel_fill_src(float *h_src, int n)
{
  std::mt19937 gen(19937);
  std::uniform_real_distribution<float> dis1(-5.0f, 5.0f);
  std::uniform_real_distribution<float> dis2(0.3f, 5.0f);
  std::uniform_real_distribution<float> dis3(-1.0f, 1.0f);
  std::uniform_real_distribution<float> dis4(4.0e-4f, 2.5e-3f);

  for (int i = 0; i < n; i++) {
    h_src[i * COL_DIM + COL_P_X] = dis1(gen);
    h_src[i * COL_DIM + COL_P_Y] = dis1(gen);
    h_src[i * COL_DIM + COL_P_Z] = dis2(gen);
    float nx = dis3(gen);
    float ny = dis3(gen);
    float nz = dis3(gen);
    float s = std::sqrt(nx * nx + ny * ny + nz * nz);
    h_src[i * COL_DIM + COL_N_X] = nx / s;
    h_src[i * COL_DIM + COL_N_Y] = ny / s;
    h_src[i * COL_DIM + COL_N_Z] = nz / s;
    h_src[i * COL_DIM + COL_RSq] = dis4(gen);
  }
}
