! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_int8_t
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(c_int8_t), save :: sbox(0:15) = [ &
    -64_c_int8_t, 80_c_int8_t, 96_c_int8_t, -80_c_int8_t, -112_c_int8_t, 0_c_int8_t, -96_c_int8_t, -48_c_int8_t, 48_c_int8_t, -32_c_int8_t, -16_c_int8_t, -128_c_int8_t, 64_c_int8_t, 112_c_int8_t, 16_c_int8_t, 32_c_int8_t &
  ]
  integer(c_int8_t), save :: sbox_pmt_3(0:255) = [ &
    -16_c_int8_t, -79_c_int8_t, -76_c_int8_t, -27_c_int8_t, -31_c_int8_t, -96_c_int8_t, -28_c_int8_t, -15_c_int8_t, -91_c_int8_t, -12_c_int8_t, -11_c_int8_t, -32_c_int8_t, -80_c_int8_t, -75_c_int8_t, -95_c_int8_t, -92_c_int8_t, &
    114_c_int8_t, 51_c_int8_t, 54_c_int8_t, 103_c_int8_t, 99_c_int8_t, 34_c_int8_t, 102_c_int8_t, 115_c_int8_t, 39_c_int8_t, 118_c_int8_t, 119_c_int8_t, 98_c_int8_t, 50_c_int8_t, 55_c_int8_t, 35_c_int8_t, 38_c_int8_t, &
    120_c_int8_t, 57_c_int8_t, 60_c_int8_t, 109_c_int8_t, 105_c_int8_t, 40_c_int8_t, 108_c_int8_t, 121_c_int8_t, 45_c_int8_t, 124_c_int8_t, 125_c_int8_t, 104_c_int8_t, 56_c_int8_t, 61_c_int8_t, 41_c_int8_t, 44_c_int8_t, &
    -38_c_int8_t, -101_c_int8_t, -98_c_int8_t, -49_c_int8_t, -53_c_int8_t, -118_c_int8_t, -50_c_int8_t, -37_c_int8_t, -113_c_int8_t, -34_c_int8_t, -33_c_int8_t, -54_c_int8_t, -102_c_int8_t, -97_c_int8_t, -117_c_int8_t, -114_c_int8_t, &
    -46_c_int8_t, -109_c_int8_t, -106_c_int8_t, -57_c_int8_t, -61_c_int8_t, -126_c_int8_t, -58_c_int8_t, -45_c_int8_t, -121_c_int8_t, -42_c_int8_t, -41_c_int8_t, -62_c_int8_t, -110_c_int8_t, -105_c_int8_t, -125_c_int8_t, -122_c_int8_t, &
    80_c_int8_t, 17_c_int8_t, 20_c_int8_t, 69_c_int8_t, 65_c_int8_t, 0_c_int8_t, 68_c_int8_t, 81_c_int8_t, 5_c_int8_t, 84_c_int8_t, 85_c_int8_t, 64_c_int8_t, 16_c_int8_t, 21_c_int8_t, 1_c_int8_t, 4_c_int8_t, &
    -40_c_int8_t, -103_c_int8_t, -100_c_int8_t, -51_c_int8_t, -55_c_int8_t, -120_c_int8_t, -52_c_int8_t, -39_c_int8_t, -115_c_int8_t, -36_c_int8_t, -35_c_int8_t, -56_c_int8_t, -104_c_int8_t, -99_c_int8_t, -119_c_int8_t, -116_c_int8_t, &
    -14_c_int8_t, -77_c_int8_t, -74_c_int8_t, -25_c_int8_t, -29_c_int8_t, -94_c_int8_t, -26_c_int8_t, -13_c_int8_t, -89_c_int8_t, -10_c_int8_t, -9_c_int8_t, -30_c_int8_t, -78_c_int8_t, -73_c_int8_t, -93_c_int8_t, -90_c_int8_t, &
    90_c_int8_t, 27_c_int8_t, 30_c_int8_t, 79_c_int8_t, 75_c_int8_t, 10_c_int8_t, 78_c_int8_t, 91_c_int8_t, 15_c_int8_t, 94_c_int8_t, 95_c_int8_t, 74_c_int8_t, 26_c_int8_t, 31_c_int8_t, 11_c_int8_t, 14_c_int8_t, &
    -8_c_int8_t, -71_c_int8_t, -68_c_int8_t, -19_c_int8_t, -23_c_int8_t, -88_c_int8_t, -20_c_int8_t, -7_c_int8_t, -83_c_int8_t, -4_c_int8_t, -3_c_int8_t, -24_c_int8_t, -72_c_int8_t, -67_c_int8_t, -87_c_int8_t, -84_c_int8_t, &
    -6_c_int8_t, -69_c_int8_t, -66_c_int8_t, -17_c_int8_t, -21_c_int8_t, -86_c_int8_t, -18_c_int8_t, -5_c_int8_t, -81_c_int8_t, -2_c_int8_t, -1_c_int8_t, -22_c_int8_t, -70_c_int8_t, -65_c_int8_t, -85_c_int8_t, -82_c_int8_t, &
    -48_c_int8_t, -111_c_int8_t, -108_c_int8_t, -59_c_int8_t, -63_c_int8_t, -128_c_int8_t, -60_c_int8_t, -47_c_int8_t, -123_c_int8_t, -44_c_int8_t, -43_c_int8_t, -64_c_int8_t, -112_c_int8_t, -107_c_int8_t, -127_c_int8_t, -124_c_int8_t, &
    112_c_int8_t, 49_c_int8_t, 52_c_int8_t, 101_c_int8_t, 97_c_int8_t, 32_c_int8_t, 100_c_int8_t, 113_c_int8_t, 37_c_int8_t, 116_c_int8_t, 117_c_int8_t, 96_c_int8_t, 48_c_int8_t, 53_c_int8_t, 33_c_int8_t, 36_c_int8_t, &
    122_c_int8_t, 59_c_int8_t, 62_c_int8_t, 111_c_int8_t, 107_c_int8_t, 42_c_int8_t, 110_c_int8_t, 123_c_int8_t, 47_c_int8_t, 126_c_int8_t, 127_c_int8_t, 106_c_int8_t, 58_c_int8_t, 63_c_int8_t, 43_c_int8_t, 46_c_int8_t, &
    82_c_int8_t, 19_c_int8_t, 22_c_int8_t, 71_c_int8_t, 67_c_int8_t, 2_c_int8_t, 70_c_int8_t, 83_c_int8_t, 7_c_int8_t, 86_c_int8_t, 87_c_int8_t, 66_c_int8_t, 18_c_int8_t, 23_c_int8_t, 3_c_int8_t, 6_c_int8_t, &
    88_c_int8_t, 25_c_int8_t, 28_c_int8_t, 77_c_int8_t, 73_c_int8_t, 8_c_int8_t, 76_c_int8_t, 89_c_int8_t, 13_c_int8_t, 92_c_int8_t, 93_c_int8_t, 72_c_int8_t, 24_c_int8_t, 29_c_int8_t, 9_c_int8_t, 12_c_int8_t &
  ]
  integer(c_int8_t), save :: sbox_pmt_2(0:255) = [ &
    60_c_int8_t, 108_c_int8_t, 45_c_int8_t, 121_c_int8_t, 120_c_int8_t, 40_c_int8_t, 57_c_int8_t, 124_c_int8_t, 105_c_int8_t, 61_c_int8_t, 125_c_int8_t, 56_c_int8_t, 44_c_int8_t, 109_c_int8_t, 104_c_int8_t, 41_c_int8_t, &
    -100_c_int8_t, -52_c_int8_t, -115_c_int8_t, -39_c_int8_t, -40_c_int8_t, -120_c_int8_t, -103_c_int8_t, -36_c_int8_t, -55_c_int8_t, -99_c_int8_t, -35_c_int8_t, -104_c_int8_t, -116_c_int8_t, -51_c_int8_t, -56_c_int8_t, -119_c_int8_t, &
    30_c_int8_t, 78_c_int8_t, 15_c_int8_t, 91_c_int8_t, 90_c_int8_t, 10_c_int8_t, 27_c_int8_t, 94_c_int8_t, 75_c_int8_t, 31_c_int8_t, 95_c_int8_t, 26_c_int8_t, 14_c_int8_t, 79_c_int8_t, 74_c_int8_t, 11_c_int8_t, &
    -74_c_int8_t, -26_c_int8_t, -89_c_int8_t, -13_c_int8_t, -14_c_int8_t, -94_c_int8_t, -77_c_int8_t, -10_c_int8_t, -29_c_int8_t, -73_c_int8_t, -9_c_int8_t, -78_c_int8_t, -90_c_int8_t, -25_c_int8_t, -30_c_int8_t, -93_c_int8_t, &
    -76_c_int8_t, -28_c_int8_t, -91_c_int8_t, -15_c_int8_t, -16_c_int8_t, -96_c_int8_t, -79_c_int8_t, -12_c_int8_t, -31_c_int8_t, -75_c_int8_t, -11_c_int8_t, -80_c_int8_t, -92_c_int8_t, -27_c_int8_t, -32_c_int8_t, -95_c_int8_t, &
    20_c_int8_t, 68_c_int8_t, 5_c_int8_t, 81_c_int8_t, 80_c_int8_t, 0_c_int8_t, 17_c_int8_t, 84_c_int8_t, 65_c_int8_t, 21_c_int8_t, 85_c_int8_t, 16_c_int8_t, 4_c_int8_t, 69_c_int8_t, 64_c_int8_t, 1_c_int8_t, &
    54_c_int8_t, 102_c_int8_t, 39_c_int8_t, 115_c_int8_t, 114_c_int8_t, 34_c_int8_t, 51_c_int8_t, 118_c_int8_t, 99_c_int8_t, 55_c_int8_t, 119_c_int8_t, 50_c_int8_t, 38_c_int8_t, 103_c_int8_t, 98_c_int8_t, 35_c_int8_t, &
    -68_c_int8_t, -20_c_int8_t, -83_c_int8_t, -7_c_int8_t, -8_c_int8_t, -88_c_int8_t, -71_c_int8_t, -4_c_int8_t, -23_c_int8_t, -67_c_int8_t, -3_c_int8_t, -72_c_int8_t, -84_c_int8_t, -19_c_int8_t, -24_c_int8_t, -87_c_int8_t, &
    -106_c_int8_t, -58_c_int8_t, -121_c_int8_t, -45_c_int8_t, -46_c_int8_t, -126_c_int8_t, -109_c_int8_t, -42_c_int8_t, -61_c_int8_t, -105_c_int8_t, -41_c_int8_t, -110_c_int8_t, -122_c_int8_t, -57_c_int8_t, -62_c_int8_t, -125_c_int8_t, &
    62_c_int8_t, 110_c_int8_t, 47_c_int8_t, 123_c_int8_t, 122_c_int8_t, 42_c_int8_t, 59_c_int8_t, 126_c_int8_t, 107_c_int8_t, 63_c_int8_t, 127_c_int8_t, 58_c_int8_t, 46_c_int8_t, 111_c_int8_t, 106_c_int8_t, 43_c_int8_t, &
    -66_c_int8_t, -18_c_int8_t, -81_c_int8_t, -5_c_int8_t, -6_c_int8_t, -86_c_int8_t, -69_c_int8_t, -2_c_int8_t, -21_c_int8_t, -65_c_int8_t, -1_c_int8_t, -70_c_int8_t, -82_c_int8_t, -17_c_int8_t, -22_c_int8_t, -85_c_int8_t, &
    52_c_int8_t, 100_c_int8_t, 37_c_int8_t, 113_c_int8_t, 112_c_int8_t, 32_c_int8_t, 49_c_int8_t, 116_c_int8_t, 97_c_int8_t, 53_c_int8_t, 117_c_int8_t, 48_c_int8_t, 36_c_int8_t, 101_c_int8_t, 96_c_int8_t, 33_c_int8_t, &
    28_c_int8_t, 76_c_int8_t, 13_c_int8_t, 89_c_int8_t, 88_c_int8_t, 8_c_int8_t, 25_c_int8_t, 92_c_int8_t, 73_c_int8_t, 29_c_int8_t, 93_c_int8_t, 24_c_int8_t, 12_c_int8_t, 77_c_int8_t, 72_c_int8_t, 9_c_int8_t, &
    -98_c_int8_t, -50_c_int8_t, -113_c_int8_t, -37_c_int8_t, -38_c_int8_t, -118_c_int8_t, -101_c_int8_t, -34_c_int8_t, -53_c_int8_t, -97_c_int8_t, -33_c_int8_t, -102_c_int8_t, -114_c_int8_t, -49_c_int8_t, -54_c_int8_t, -117_c_int8_t, &
    -108_c_int8_t, -60_c_int8_t, -123_c_int8_t, -47_c_int8_t, -48_c_int8_t, -128_c_int8_t, -111_c_int8_t, -44_c_int8_t, -63_c_int8_t, -107_c_int8_t, -43_c_int8_t, -112_c_int8_t, -124_c_int8_t, -59_c_int8_t, -64_c_int8_t, -127_c_int8_t, &
    22_c_int8_t, 70_c_int8_t, 7_c_int8_t, 83_c_int8_t, 82_c_int8_t, 2_c_int8_t, 19_c_int8_t, 86_c_int8_t, 67_c_int8_t, 23_c_int8_t, 87_c_int8_t, 18_c_int8_t, 6_c_int8_t, 71_c_int8_t, 66_c_int8_t, 3_c_int8_t &
  ]
  integer(c_int8_t), save :: sbox_pmt_1(0:255) = [ &
    15_c_int8_t, 27_c_int8_t, 75_c_int8_t, 94_c_int8_t, 30_c_int8_t, 10_c_int8_t, 78_c_int8_t, 31_c_int8_t, 90_c_int8_t, 79_c_int8_t, 95_c_int8_t, 14_c_int8_t, 11_c_int8_t, 91_c_int8_t, 26_c_int8_t, 74_c_int8_t, &
    39_c_int8_t, 51_c_int8_t, 99_c_int8_t, 118_c_int8_t, 54_c_int8_t, 34_c_int8_t, 102_c_int8_t, 55_c_int8_t, 114_c_int8_t, 103_c_int8_t, 119_c_int8_t, 38_c_int8_t, 35_c_int8_t, 115_c_int8_t, 50_c_int8_t, 98_c_int8_t, &
    -121_c_int8_t, -109_c_int8_t, -61_c_int8_t, -42_c_int8_t, -106_c_int8_t, -126_c_int8_t, -58_c_int8_t, -105_c_int8_t, -46_c_int8_t, -57_c_int8_t, -41_c_int8_t, -122_c_int8_t, -125_c_int8_t, -45_c_int8_t, -110_c_int8_t, -62_c_int8_t, &
    -83_c_int8_t, -71_c_int8_t, -23_c_int8_t, -4_c_int8_t, -68_c_int8_t, -88_c_int8_t, -20_c_int8_t, -67_c_int8_t, -8_c_int8_t, -19_c_int8_t, -3_c_int8_t, -84_c_int8_t, -87_c_int8_t, -7_c_int8_t, -72_c_int8_t, -24_c_int8_t, &
    45_c_int8_t, 57_c_int8_t, 105_c_int8_t, 124_c_int8_t, 60_c_int8_t, 40_c_int8_t, 108_c_int8_t, 61_c_int8_t, 120_c_int8_t, 109_c_int8_t, 125_c_int8_t, 44_c_int8_t, 41_c_int8_t, 121_c_int8_t, 56_c_int8_t, 104_c_int8_t, &
    5_c_int8_t, 17_c_int8_t, 65_c_int8_t, 84_c_int8_t, 20_c_int8_t, 0_c_int8_t, 68_c_int8_t, 21_c_int8_t, 80_c_int8_t, 69_c_int8_t, 85_c_int8_t, 4_c_int8_t, 1_c_int8_t, 81_c_int8_t, 16_c_int8_t, 64_c_int8_t, &
    -115_c_int8_t, -103_c_int8_t, -55_c_int8_t, -36_c_int8_t, -100_c_int8_t, -120_c_int8_t, -52_c_int8_t, -99_c_int8_t, -40_c_int8_t, -51_c_int8_t, -35_c_int8_t, -116_c_int8_t, -119_c_int8_t, -39_c_int8_t, -104_c_int8_t, -56_c_int8_t, &
    47_c_int8_t, 59_c_int8_t, 107_c_int8_t, 126_c_int8_t, 62_c_int8_t, 42_c_int8_t, 110_c_int8_t, 63_c_int8_t, 122_c_int8_t, 111_c_int8_t, 127_c_int8_t, 46_c_int8_t, 43_c_int8_t, 123_c_int8_t, 58_c_int8_t, 106_c_int8_t, &
    -91_c_int8_t, -79_c_int8_t, -31_c_int8_t, -12_c_int8_t, -76_c_int8_t, -96_c_int8_t, -28_c_int8_t, -75_c_int8_t, -16_c_int8_t, -27_c_int8_t, -11_c_int8_t, -92_c_int8_t, -95_c_int8_t, -15_c_int8_t, -80_c_int8_t, -32_c_int8_t, &
    -113_c_int8_t, -101_c_int8_t, -53_c_int8_t, -34_c_int8_t, -98_c_int8_t, -118_c_int8_t, -50_c_int8_t, -97_c_int8_t, -38_c_int8_t, -49_c_int8_t, -33_c_int8_t, -114_c_int8_t, -117_c_int8_t, -37_c_int8_t, -102_c_int8_t, -54_c_int8_t, &
    -81_c_int8_t, -69_c_int8_t, -21_c_int8_t, -2_c_int8_t, -66_c_int8_t, -86_c_int8_t, -18_c_int8_t, -65_c_int8_t, -6_c_int8_t, -17_c_int8_t, -1_c_int8_t, -82_c_int8_t, -85_c_int8_t, -5_c_int8_t, -70_c_int8_t, -22_c_int8_t, &
    13_c_int8_t, 25_c_int8_t, 73_c_int8_t, 92_c_int8_t, 28_c_int8_t, 8_c_int8_t, 76_c_int8_t, 29_c_int8_t, 88_c_int8_t, 77_c_int8_t, 93_c_int8_t, 12_c_int8_t, 9_c_int8_t, 89_c_int8_t, 24_c_int8_t, 72_c_int8_t, &
    7_c_int8_t, 19_c_int8_t, 67_c_int8_t, 86_c_int8_t, 22_c_int8_t, 2_c_int8_t, 70_c_int8_t, 23_c_int8_t, 82_c_int8_t, 71_c_int8_t, 87_c_int8_t, 6_c_int8_t, 3_c_int8_t, 83_c_int8_t, 18_c_int8_t, 66_c_int8_t, &
    -89_c_int8_t, -77_c_int8_t, -29_c_int8_t, -10_c_int8_t, -74_c_int8_t, -94_c_int8_t, -26_c_int8_t, -73_c_int8_t, -14_c_int8_t, -25_c_int8_t, -9_c_int8_t, -90_c_int8_t, -93_c_int8_t, -13_c_int8_t, -78_c_int8_t, -30_c_int8_t, &
    37_c_int8_t, 49_c_int8_t, 97_c_int8_t, 116_c_int8_t, 52_c_int8_t, 32_c_int8_t, 100_c_int8_t, 53_c_int8_t, 112_c_int8_t, 101_c_int8_t, 117_c_int8_t, 36_c_int8_t, 33_c_int8_t, 113_c_int8_t, 48_c_int8_t, 96_c_int8_t, &
    -123_c_int8_t, -111_c_int8_t, -63_c_int8_t, -44_c_int8_t, -108_c_int8_t, -128_c_int8_t, -60_c_int8_t, -107_c_int8_t, -48_c_int8_t, -59_c_int8_t, -43_c_int8_t, -124_c_int8_t, -127_c_int8_t, -47_c_int8_t, -112_c_int8_t, -64_c_int8_t &
  ]
  integer(c_int8_t), save :: sbox_pmt_0(0:255) = [ &
    -61_c_int8_t, -58_c_int8_t, -46_c_int8_t, -105_c_int8_t, -121_c_int8_t, -126_c_int8_t, -109_c_int8_t, -57_c_int8_t, -106_c_int8_t, -45_c_int8_t, -41_c_int8_t, -125_c_int8_t, -62_c_int8_t, -42_c_int8_t, -122_c_int8_t, -110_c_int8_t, &
    -55_c_int8_t, -52_c_int8_t, -40_c_int8_t, -99_c_int8_t, -115_c_int8_t, -120_c_int8_t, -103_c_int8_t, -51_c_int8_t, -100_c_int8_t, -39_c_int8_t, -35_c_int8_t, -119_c_int8_t, -56_c_int8_t, -36_c_int8_t, -116_c_int8_t, -104_c_int8_t, &
    -31_c_int8_t, -28_c_int8_t, -16_c_int8_t, -75_c_int8_t, -91_c_int8_t, -96_c_int8_t, -79_c_int8_t, -27_c_int8_t, -76_c_int8_t, -15_c_int8_t, -11_c_int8_t, -95_c_int8_t, -32_c_int8_t, -12_c_int8_t, -92_c_int8_t, -80_c_int8_t, &
    107_c_int8_t, 110_c_int8_t, 122_c_int8_t, 63_c_int8_t, 47_c_int8_t, 42_c_int8_t, 59_c_int8_t, 111_c_int8_t, 62_c_int8_t, 123_c_int8_t, 127_c_int8_t, 43_c_int8_t, 106_c_int8_t, 126_c_int8_t, 46_c_int8_t, 58_c_int8_t, &
    75_c_int8_t, 78_c_int8_t, 90_c_int8_t, 31_c_int8_t, 15_c_int8_t, 10_c_int8_t, 27_c_int8_t, 79_c_int8_t, 30_c_int8_t, 91_c_int8_t, 95_c_int8_t, 11_c_int8_t, 74_c_int8_t, 94_c_int8_t, 14_c_int8_t, 26_c_int8_t, &
    65_c_int8_t, 68_c_int8_t, 80_c_int8_t, 21_c_int8_t, 5_c_int8_t, 0_c_int8_t, 17_c_int8_t, 69_c_int8_t, 20_c_int8_t, 81_c_int8_t, 85_c_int8_t, 1_c_int8_t, 64_c_int8_t, 84_c_int8_t, 4_c_int8_t, 16_c_int8_t, &
    99_c_int8_t, 102_c_int8_t, 114_c_int8_t, 55_c_int8_t, 39_c_int8_t, 34_c_int8_t, 51_c_int8_t, 103_c_int8_t, 54_c_int8_t, 115_c_int8_t, 119_c_int8_t, 35_c_int8_t, 98_c_int8_t, 118_c_int8_t, 38_c_int8_t, 50_c_int8_t, &
    -53_c_int8_t, -50_c_int8_t, -38_c_int8_t, -97_c_int8_t, -113_c_int8_t, -118_c_int8_t, -101_c_int8_t, -49_c_int8_t, -98_c_int8_t, -37_c_int8_t, -33_c_int8_t, -117_c_int8_t, -54_c_int8_t, -34_c_int8_t, -114_c_int8_t, -102_c_int8_t, &
    105_c_int8_t, 108_c_int8_t, 120_c_int8_t, 61_c_int8_t, 45_c_int8_t, 40_c_int8_t, 57_c_int8_t, 109_c_int8_t, 60_c_int8_t, 121_c_int8_t, 125_c_int8_t, 41_c_int8_t, 104_c_int8_t, 124_c_int8_t, 44_c_int8_t, 56_c_int8_t, &
    -29_c_int8_t, -26_c_int8_t, -14_c_int8_t, -73_c_int8_t, -89_c_int8_t, -94_c_int8_t, -77_c_int8_t, -25_c_int8_t, -74_c_int8_t, -13_c_int8_t, -9_c_int8_t, -93_c_int8_t, -30_c_int8_t, -10_c_int8_t, -90_c_int8_t, -78_c_int8_t, &
    -21_c_int8_t, -18_c_int8_t, -6_c_int8_t, -65_c_int8_t, -81_c_int8_t, -86_c_int8_t, -69_c_int8_t, -17_c_int8_t, -66_c_int8_t, -5_c_int8_t, -1_c_int8_t, -85_c_int8_t, -22_c_int8_t, -2_c_int8_t, -82_c_int8_t, -70_c_int8_t, &
    67_c_int8_t, 70_c_int8_t, 82_c_int8_t, 23_c_int8_t, 7_c_int8_t, 2_c_int8_t, 19_c_int8_t, 71_c_int8_t, 22_c_int8_t, 83_c_int8_t, 87_c_int8_t, 3_c_int8_t, 66_c_int8_t, 86_c_int8_t, 6_c_int8_t, 18_c_int8_t, &
    -63_c_int8_t, -60_c_int8_t, -48_c_int8_t, -107_c_int8_t, -123_c_int8_t, -128_c_int8_t, -111_c_int8_t, -59_c_int8_t, -108_c_int8_t, -47_c_int8_t, -43_c_int8_t, -127_c_int8_t, -64_c_int8_t, -44_c_int8_t, -124_c_int8_t, -112_c_int8_t, &
    -23_c_int8_t, -20_c_int8_t, -8_c_int8_t, -67_c_int8_t, -83_c_int8_t, -88_c_int8_t, -71_c_int8_t, -19_c_int8_t, -68_c_int8_t, -7_c_int8_t, -3_c_int8_t, -87_c_int8_t, -24_c_int8_t, -4_c_int8_t, -84_c_int8_t, -72_c_int8_t, &
    73_c_int8_t, 76_c_int8_t, 88_c_int8_t, 29_c_int8_t, 13_c_int8_t, 8_c_int8_t, 25_c_int8_t, 77_c_int8_t, 28_c_int8_t, 89_c_int8_t, 93_c_int8_t, 9_c_int8_t, 72_c_int8_t, 92_c_int8_t, 12_c_int8_t, 24_c_int8_t, &
    97_c_int8_t, 100_c_int8_t, 112_c_int8_t, 53_c_int8_t, 37_c_int8_t, 32_c_int8_t, 49_c_int8_t, 101_c_int8_t, 52_c_int8_t, 113_c_int8_t, 117_c_int8_t, 33_c_int8_t, 96_c_int8_t, 116_c_int8_t, 36_c_int8_t, 48_c_int8_t &
  ]
  !$omp declare target(sbox, sbox_pmt_3, sbox_pmt_2, sbox_pmt_1, sbox_pmt_0)


  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer :: num, repeat, argc
  character(len=256) :: arg
  integer(c_int8_t), allocatable :: h_plain(:), h_key(:), h_cipher(:), ciphers(:)
  integer(int64) :: h_checksum, d_checksum
  integer(c_int8_t) :: plain(0:7), key(0:9)
  integer :: i, k, n
  integer(int32), parameter :: rounds = 31_int32
  integer(int32), parameter :: seed = 8_int32
  real(real64) :: start_time, end_time, elapsed

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, arg)
    write(*,'(A,A,A)') 'Usage: ', trim(arg), ' <number of plain texts> <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) num
  call get_command_argument(2, arg)
  read(arg, *) repeat
  if (num <= 0 .or. repeat <= 0) stop 1

  allocate(h_plain(0:8 * num - 1), h_key(0:10 * num - 1), h_cipher(0:8 * num - 1))
  allocate(ciphers(0:8 * num - 1))

  call c_srand(int(seed, c_int))
  plain = [int(iachar('P'), c_int8_t), int(iachar('R'), c_int8_t), int(iachar('E'), c_int8_t), &
    int(iachar('S'), c_int8_t), int(iachar('E'), c_int8_t), int(iachar('N'), c_int8_t), &
    int(iachar('T'), c_int8_t), 0_c_int8_t]

  do i = 0, num - 1
    do k = 0, 9
      key(k) = byte_value(mod(c_rand(), 256))
      h_key(i * 10 + k) = key(k)
    end do
    do k = 0, 7
      h_plain(i * 8 + k) = plain(k)
    end do
    call deterministic_shuffle(plain)
  end do

  h_checksum = 0_int64
  do n = 0, repeat
    do i = 0, num - 1
      call present_rounds(h_plain(i * 8:), h_key(i * 10:), rounds, h_cipher(i * 8:))
      do k = 0, 7
        h_checksum = h_checksum + int(byte_index(h_cipher(i * 8 + k)), int64)
      end do
    end do
  end do

  ciphers = 0_c_int8_t
  d_checksum = 0_int64
  elapsed = 0.0_real64

  !$omp target data map(to: h_plain(0:8 * num - 1), h_key(0:10 * num - 1)) &
  !$omp& map(alloc: ciphers(0:8 * num - 1))
  do n = 0, repeat
    start_time = omp_get_wtime()
    call present_kernel(num, h_plain, h_key, rounds, ciphers)
    end_time = omp_get_wtime()
    if (n > 0) elapsed = elapsed + (end_time - start_time)

    !$omp target update from(ciphers(0:8 * num - 1))
    do i = 0, 8 * num - 1
      d_checksum = d_checksum + int(byte_index(ciphers(i)), int64)
    end do
  end do
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', (elapsed * 1.0e6_real64) / real(repeat, real64), ' (us)'
  if (h_checksum /= d_checksum) then
    write(*,'(A)') 'FAIL'
  else
    write(*,'(A)') 'PASS'
  end if

  deallocate(ciphers, h_cipher, h_key, h_plain)

contains

  subroutine deterministic_shuffle(values)
    integer(c_int8_t), intent(inout) :: values(0:7)
    integer(int32), parameter :: order(0:7) = [4_int32, 5_int32, 1_int32, 2_int32, 6_int32, 7_int32, 3_int32, 0_int32]
    integer(c_int8_t) :: tmp(0:7)
    integer :: j

    tmp = values
    do j = 0, 7
      values(j) = tmp(order(j))
    end do
  end subroutine deterministic_shuffle

  subroutine present_kernel(num_items, plains, keys, rounds, ciphers)
    integer, intent(in) :: num_items
    integer(c_int8_t), intent(in) :: plains(0:), keys(0:)
    integer(int32), intent(in) :: rounds
    integer(c_int8_t), intent(inout) :: ciphers(0:)
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256)
    do idx = 0, num_items - 1
      call present_rounds(plains(idx * 8:), keys(idx * 10:), rounds, ciphers(idx * 8:))
    end do
  end subroutine present_kernel

  subroutine present_rounds(plain, key, rounds, cipher)
    !$omp declare target
    integer(c_int8_t), intent(in) :: plain(0:), key(0:)
    integer(int32), intent(in) :: rounds
    integer(c_int8_t), intent(inout) :: cipher(0:)
    integer(int32) :: round_counter
    integer(c_int8_t) :: state(0:7), round_key(0:9)

    state(0) = bxor8(plain(0), key(0))
    state(1) = bxor8(plain(1), key(1))
    state(2) = bxor8(plain(2), key(2))
    state(3) = bxor8(plain(3), key(3))
    state(4) = bxor8(plain(4), key(4))
    state(5) = bxor8(plain(5), key(5))
    state(6) = bxor8(plain(6), key(6))
    state(7) = bxor8(plain(7), key(7))

    round_key(9) = byte_value(bor(shl8(key(6), 5), shr8(key(7), 3)))
    round_key(8) = byte_value(bor(shl8(key(5), 5), shr8(key(6), 3)))
    round_key(7) = byte_value(bor(shl8(key(4), 5), shr8(key(5), 3)))
    round_key(6) = byte_value(bor(shl8(key(3), 5), shr8(key(4), 3)))
    round_key(5) = byte_value(bor(shl8(key(2), 5), shr8(key(3), 3)))
    round_key(4) = byte_value(bor(shl8(key(1), 5), shr8(key(2), 3)))
    round_key(3) = byte_value(bor(shl8(key(0), 5), shr8(key(1), 3)))
    round_key(2) = byte_value(bor(shl8(key(9), 5), shr8(key(0), 3)))
    round_key(1) = byte_value(bor(shl8(key(8), 5), shr8(key(9), 3)))
    round_key(0) = byte_value(bor(shl8(key(7), 5), shr8(key(8), 3)))
    round_key(0) = byte_value(bor(band(byte_index(round_key(0)), z'0F'), byte_index(sbox(shr8(round_key(0), 4)))))
    round_key(7) = bxor8(round_key(7), byte_value(shr8(byte_value(1), 1)))
    round_key(8) = bxor8(round_key(8), byte_value(shl8(byte_value(1), 7)))

    call substitute_permute(state, cipher)

    do round_counter = 2_int32, rounds
      state(0) = bxor8(cipher(0), round_key(0))
      state(1) = bxor8(cipher(1), round_key(1))
      state(2) = bxor8(cipher(2), round_key(2))
      state(3) = bxor8(cipher(3), round_key(3))
      state(4) = bxor8(cipher(4), round_key(4))
      state(5) = bxor8(cipher(5), round_key(5))
      state(6) = bxor8(cipher(6), round_key(6))
      state(7) = bxor8(cipher(7), round_key(7))

      call substitute_permute(state, cipher)
      call update_round_key(round_key, state, round_counter)
    end do

    if (rounds == 31_int32) then
      cipher(0) = bxor8(cipher(0), round_key(0))
      cipher(1) = bxor8(cipher(1), round_key(1))
      cipher(2) = bxor8(cipher(2), round_key(2))
      cipher(3) = bxor8(cipher(3), round_key(3))
      cipher(4) = bxor8(cipher(4), round_key(4))
      cipher(5) = bxor8(cipher(5), round_key(5))
      cipher(6) = bxor8(cipher(6), round_key(6))
      cipher(7) = bxor8(cipher(7), round_key(7))
    end if
  end subroutine present_rounds

  subroutine substitute_permute(state, cipher)
    !$omp declare target
    integer(c_int8_t), intent(in) :: state(0:7)
    integer(c_int8_t), intent(inout) :: cipher(0:7)

    cipher(0) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_3(byte_index(state(0)))), z'C0'), band(byte_index(sbox_pmt_2(byte_index(state(1)))), z'30')), &
      band(byte_index(sbox_pmt_1(byte_index(state(2)))), z'0C')), band(byte_index(sbox_pmt_0(byte_index(state(3)))), z'03')))
    cipher(1) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_3(byte_index(state(4)))), z'C0'), band(byte_index(sbox_pmt_2(byte_index(state(5)))), z'30')), &
      band(byte_index(sbox_pmt_1(byte_index(state(6)))), z'0C')), band(byte_index(sbox_pmt_0(byte_index(state(7)))), z'03')))
    cipher(2) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_0(byte_index(state(0)))), z'C0'), band(byte_index(sbox_pmt_3(byte_index(state(1)))), z'30')), &
      band(byte_index(sbox_pmt_2(byte_index(state(2)))), z'0C')), band(byte_index(sbox_pmt_1(byte_index(state(3)))), z'03')))
    cipher(3) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_0(byte_index(state(4)))), z'C0'), band(byte_index(sbox_pmt_3(byte_index(state(5)))), z'30')), &
      band(byte_index(sbox_pmt_2(byte_index(state(6)))), z'0C')), band(byte_index(sbox_pmt_1(byte_index(state(7)))), z'03')))
    cipher(4) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_1(byte_index(state(0)))), z'C0'), band(byte_index(sbox_pmt_0(byte_index(state(1)))), z'30')), &
      band(byte_index(sbox_pmt_3(byte_index(state(2)))), z'0C')), band(byte_index(sbox_pmt_2(byte_index(state(3)))), z'03')))
    cipher(5) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_1(byte_index(state(4)))), z'C0'), band(byte_index(sbox_pmt_0(byte_index(state(5)))), z'30')), &
      band(byte_index(sbox_pmt_3(byte_index(state(6)))), z'0C')), band(byte_index(sbox_pmt_2(byte_index(state(7)))), z'03')))
    cipher(6) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_2(byte_index(state(0)))), z'C0'), band(byte_index(sbox_pmt_1(byte_index(state(1)))), z'30')), &
      band(byte_index(sbox_pmt_0(byte_index(state(2)))), z'0C')), band(byte_index(sbox_pmt_3(byte_index(state(3)))), z'03')))
    cipher(7) = byte_value(bor(bor(bor(band(byte_index(sbox_pmt_2(byte_index(state(4)))), z'C0'), band(byte_index(sbox_pmt_1(byte_index(state(5)))), z'30')), &
      band(byte_index(sbox_pmt_0(byte_index(state(6)))), z'0C')), band(byte_index(sbox_pmt_3(byte_index(state(7)))), z'03')))
  end subroutine substitute_permute

  subroutine update_round_key(round_key, scratch, round_counter)
    !$omp declare target
    integer(c_int8_t), intent(inout) :: round_key(0:9)
    integer(c_int8_t), intent(inout) :: scratch(0:7)
    integer(int32), intent(in) :: round_counter

    round_key(5) = bxor8(round_key(5), byte_value(shl8(byte_value(round_counter), 2)))
    scratch(2) = round_key(9)
    scratch(1) = round_key(8)
    scratch(0) = round_key(7)
    round_key(9) = byte_value(bor(shl8(round_key(6), 5), shr8(round_key(7), 3)))
    round_key(8) = byte_value(bor(shl8(round_key(5), 5), shr8(round_key(6), 3)))
    round_key(7) = byte_value(bor(shl8(round_key(4), 5), shr8(round_key(5), 3)))
    round_key(6) = byte_value(bor(shl8(round_key(3), 5), shr8(round_key(4), 3)))
    round_key(5) = byte_value(bor(shl8(round_key(2), 5), shr8(round_key(3), 3)))
    round_key(4) = byte_value(bor(shl8(round_key(1), 5), shr8(round_key(2), 3)))
    round_key(3) = byte_value(bor(shl8(round_key(0), 5), shr8(round_key(1), 3)))
    round_key(2) = byte_value(bor(shl8(scratch(2), 5), shr8(round_key(0), 3)))
    round_key(1) = byte_value(bor(shl8(scratch(1), 5), shr8(scratch(2), 3)))
    round_key(0) = byte_value(bor(shl8(scratch(0), 5), shr8(scratch(1), 3)))
    round_key(0) = byte_value(bor(band(byte_index(round_key(0)), z'0F'), byte_index(sbox(shr8(round_key(0), 4)))))
  end subroutine update_round_key

  integer(int32) function byte_index(value) result(unsigned_value)
    !$omp declare target
    integer(c_int8_t), intent(in) :: value
    unsigned_value = iand(int(value, int32), 255_int32)
  end function byte_index

  integer(c_int8_t) function byte_value(value) result(byte)
    !$omp declare target
    integer(int32), intent(in) :: value
    byte = int(iand(value, 255_int32), c_int8_t)
  end function byte_value

  integer(int32) function band(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(lhs, rhs)
  end function band

  integer(int32) function bor(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(ior(lhs, rhs), 255_int32)
  end function bor

  integer(int32) function bxor(lhs, rhs) result(value)
    !$omp declare target
    integer(int32), intent(in) :: lhs, rhs
    value = iand(ieor(lhs, rhs), 255_int32)
  end function bxor

  integer(c_int8_t) function bxor8(lhs, rhs) result(value)
    !$omp declare target
    integer(c_int8_t), intent(in) :: lhs, rhs
    value = byte_value(bxor(byte_index(lhs), byte_index(rhs)))
  end function bxor8

  integer(int32) function shl8(lhs, amount) result(value)
    !$omp declare target
    integer(c_int8_t), intent(in) :: lhs
    integer(int32), intent(in) :: amount
    value = iand(ishft(byte_index(lhs), amount), 255_int32)
  end function shl8

  integer(int32) function shr8(lhs, amount) result(value)
    !$omp declare target
    integer(c_int8_t), intent(in) :: lhs
    integer(int32), intent(in) :: amount
    value = ishft(byte_index(lhs), -amount)
  end function shr8

end program main
