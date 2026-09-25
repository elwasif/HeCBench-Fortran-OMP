// SPDX-License-Identifier: CC0-1.0
#include <cstddef>
#include <random>

extern "C" void dp_fill_real32(float *srcA, float *srcB,
                               long long iNumElements,
                               long long src_size,
                               float *dst_ref) {
  std::mt19937 engine(19937);
  std::uniform_int_distribution<int> dis(-32, 32);

  *dst_ref = 0.0f;
  for (long long i = 0; i < iNumElements; ++i) {
    srcA[i] = static_cast<float>(dis(engine));
    srcB[i] = static_cast<float>(dis(engine));
    *dst_ref += srcA[i] * srcB[i];
  }
  for (long long i = iNumElements; i < src_size; ++i) {
    srcA[i] = 0.0f;
    srcB[i] = 0.0f;
  }
}

extern "C" void dp_fill_real64(double *srcA, double *srcB,
                               long long iNumElements,
                               long long src_size,
                               double *dst_ref) {
  std::mt19937 engine(19937);
  std::uniform_int_distribution<int> dis(-32, 32);

  *dst_ref = 0.0;
  for (long long i = 0; i < iNumElements; ++i) {
    srcA[i] = static_cast<double>(dis(engine));
    srcB[i] = static_cast<double>(dis(engine));
    *dst_ref += srcA[i] * srcB[i];
  }
  for (long long i = iNumElements; i < src_size; ++i) {
    srcA[i] = 0.0;
    srcB[i] = 0.0;
  }
}
