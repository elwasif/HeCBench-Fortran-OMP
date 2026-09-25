// SPDX-License-Identifier: CC0-1.0
#include <stdint.h>

double damage_lcg_random_double(uint64_t *seed)
{
  const uint64_t mask = (1ULL << 63) - 1ULL;
  const uint64_t multiplier = 2806196910506780709ULL;

  *seed = (multiplier * (*seed) + 1ULL) & mask;
  return (double)(*seed) / 9223372036854775808.0;
}
