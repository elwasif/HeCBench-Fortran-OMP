// SPDX-License-Identifier: CC0-1.0
#include <cstdint>
#include <random>

extern "C" void adamw_initialize_inputs(
    int64_t vector_size,
    int64_t float_size,
    float *g,
    float *p,
    float *p_ref,
    float *m_qscale,
    float *v_qscale,
    float *m_qscale_ref,
    float *v_qscale_ref,
    int8_t *m,
    int8_t *v,
    int8_t *m_ref,
    int8_t *v_ref)
{
  std::mt19937 gen(19937);
  std::uniform_real_distribution<float> dist(0, 1);

  for (int64_t i = 0; i < float_size; i++) {
    m_qscale[i] = dist(gen);
    v_qscale[i] = dist(gen);
    g[i] = dist(gen);
    p[i] = dist(gen);
    p_ref[i] = p[i];
    m_qscale_ref[i] = m_qscale[i];
    v_qscale_ref[i] = v_qscale[i];
  }

  for (int64_t i = 0; i < vector_size; i++) {
    m[i] = static_cast<int8_t>(256 * dist(gen));
    m_ref[i] = m[i];
    v[i] = static_cast<int8_t>(256 * dist(gen));
    v_ref[i] = v[i];
  }
}
