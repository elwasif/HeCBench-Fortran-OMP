! SPDX-License-Identifier: CC0-1.0
program rodrigues_main
  use, intrinsic :: iso_fortran_env, only: real32
  use, intrinsic :: iso_c_binding, only: c_char, c_float, c_int, c_null_char
  use omp_lib
  implicit none

  type, bind(C) :: float3
    real(c_float) :: x, y, z, pad
  end type float3

  type, bind(C) :: float4
    real(c_float) :: x, y, z, w
  end type float4

  integer :: n, repeat, i
  type(float3), allocatable :: h(:), h_ref(:)
  type(float4), allocatable :: h2(:), h2_ref(:)
  real(real32) :: wx, wy, wz, norm, angle, max_error
  real(real32), parameter :: tolerance = 5.0e-4_real32
  real(real32) :: a, b, c, d
  real(8) :: start_time, end_time

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name="rand")
      import :: c_int
      integer(c_int) :: value
    end function c_rand

    function c_atoi(str) result(value) bind(C, name="atoi")
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: str(*)
      integer(c_int) :: value
    end function c_atoi
  end interface

  call parse_args(n, repeat)

  wx = -0.3_real32
  wy = -0.6_real32
  wz = 0.15_real32
  norm = 1.0_real32 / sqrt(wx * wx + wy * wy + wz * wz)
  wx = wx * norm
  wy = wy * norm
  wz = wz * norm
  angle = 0.5_real32

  allocate(h(n), h_ref(n))
  allocate(h2(n), h2_ref(n))

  call c_srand(123_c_int)
  do i = 1, n
    a = real(c_rand(), real32)
    b = real(c_rand(), real32)
    c = real(c_rand(), real32)
    d = sqrt(a * a + b * b + c * c)
    h(i) = float3(a / d, b / d, c / d, 0.0_real32)
    h2(i) = float4(a / d, b / d, c / d, 0.0_real32)
  end do

  h_ref = h
  h2_ref = h2

  !$omp target data map(to: h, h2)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call rotate3_device(n, angle, wx, wy, wz, h)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') "Average kernel execution time (float3): ", &
       real((end_time - start_time) * 1.0d6 / dble(repeat), real32), " (us)"

  start_time = omp_get_wtime()
  do i = 1, repeat
    call rotate4_device(n, angle, wx, wy, wz, h2)
  end do
  end_time = omp_get_wtime()
  !$omp target update from(h, h2)
  !$omp end target data
  write(*,'(A,F0.6,A)') "Average kernel execution time (float4): ", &
       real((end_time - start_time) * 1.0d6 / dble(repeat), real32), " (us)"

  do i = 1, repeat
    call rotate3_cpu(n, angle, wx, wy, wz, h_ref)
    call rotate4_cpu(n, angle, wx, wy, wz, h2_ref)
  end do

  max_error = 0.0_real32
  do i = 1, n
    max_error = max(max_error, abs(h(i)%x - h_ref(i)%x))
    max_error = max(max_error, abs(h(i)%y - h_ref(i)%y))
    max_error = max(max_error, abs(h(i)%z - h_ref(i)%z))
    max_error = max(max_error, abs(h2(i)%x - h2_ref(i)%x))
    max_error = max(max_error, abs(h2(i)%y - h2_ref(i)%y))
    max_error = max(max_error, abs(h2(i)%z - h2_ref(i)%z))
    max_error = max(max_error, abs(h2(i)%w - h2_ref(i)%w))
  end do
  if (max_error > tolerance) then
    write(*,'(A,ES12.4)') "FAIL: max rotation error = ", max_error
    error stop 1
  end if

contains

  subroutine parse_args(n, repeat)
    integer, intent(out) :: n, repeat
    character(len=64) :: arg
    character(kind=c_char, len=65) :: c_arg
    integer :: argc, arg_len

    argc = command_argument_count()
    if (argc /= 2) then
      call get_command_argument(0, arg)
      write(*,'(A,A,A)') "Usage: ", trim(arg), " <number of points> <repeat>"
      stop 1
    end if

    call get_command_argument(1, arg)
    arg_len = min(len_trim(arg), len(c_arg) - 1)
    c_arg = c_null_char
    c_arg(1:arg_len) = arg(1:arg_len)
    n = int(c_atoi(c_arg), kind(n))

    call get_command_argument(2, arg)
    arg_len = min(len_trim(arg), len(c_arg) - 1)
    c_arg = c_null_char
    c_arg(1:arg_len) = arg(1:arg_len)
    repeat = int(c_atoi(c_arg), kind(repeat))
  end subroutine parse_args

  subroutine rotate3_device(n, angle, wx, wy, wz, d)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    type(float3), intent(inout) :: d(n)
    integer :: i
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(s, c, mc, px, py, pz, m1, m2, m3, m4, m5, m6, m7, m8, m9)
    do i = 1, n
      s = sin(angle)
      c = cos(angle)
      mc = 1.0_real32 - c
      px = d(i)%x
      py = d(i)%y
      pz = d(i)%z
      m1 = c + wx * wx * mc
      m2 = wz * s + wx * wy * mc
      m3 = -wy * s + wx * wz * mc
      m4 = -wz * s + wx * wy * mc
      m5 = c + wy * wy * mc
      m6 = wx * s + wy * wz * mc
      m7 = wy * s + wx * wz * mc
      m8 = -wx * s + wy * wz * mc
      m9 = c + wz * wz * mc
      d(i)%x = px * m1 + py * m2 + pz * m3
      d(i)%y = px * m4 + py * m5 + pz * m6
      d(i)%z = px * m7 + py * m8 + pz * m9
    end do
    !$omp end target teams distribute parallel do
  end subroutine rotate3_device

  subroutine rotate4_device(n, angle, wx, wy, wz, d)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    type(float4), intent(inout) :: d(n)
    integer :: i
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(s, c, mc, px, py, pz, m1, m2, m3, m4, m5, m6, m7, m8, m9)
    do i = 1, n
      s = sin(angle)
      c = cos(angle)
      mc = 1.0_real32 - c
      px = d(i)%x
      py = d(i)%y
      pz = d(i)%z
      m1 = c + wx * wx * mc
      m2 = wz * s + wx * wy * mc
      m3 = -wy * s + wx * wz * mc
      m4 = -wz * s + wx * wy * mc
      m5 = c + wy * wy * mc
      m6 = wx * s + wy * wz * mc
      m7 = wy * s + wx * wz * mc
      m8 = -wx * s + wy * wz * mc
      m9 = c + wz * wz * mc
      d(i) = float4(px * m1 + py * m2 + pz * m3, &
                    px * m4 + py * m5 + pz * m6, &
                    px * m7 + py * m8 + pz * m9, 0.0_real32)
    end do
    !$omp end target teams distribute parallel do
  end subroutine rotate4_device

  subroutine rotate3_cpu(n, angle, wx, wy, wz, d)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    type(float3), intent(inout) :: d(n)
    integer :: i

    do i = 1, n
      call rotate_one(angle, wx, wy, wz, d(i)%x, d(i)%y, d(i)%z)
    end do
  end subroutine rotate3_cpu

  subroutine rotate4_cpu(n, angle, wx, wy, wz, d)
    integer, intent(in) :: n
    real(real32), intent(in) :: angle, wx, wy, wz
    type(float4), intent(inout) :: d(n)
    integer :: i

    do i = 1, n
      call rotate_one(angle, wx, wy, wz, d(i)%x, d(i)%y, d(i)%z)
      d(i)%w = 0.0_real32
    end do
  end subroutine rotate4_cpu

  subroutine rotate_one(angle, wx, wy, wz, x, y, z)
    real(real32), intent(in) :: angle, wx, wy, wz
    real(real32), intent(inout) :: x, y, z
    real(real32) :: s, c, mc, px, py, pz
    real(real32) :: m1, m2, m3, m4, m5, m6, m7, m8, m9

    s = sin(angle)
    c = cos(angle)
    mc = 1.0_real32 - c
    px = x
    py = y
    pz = z
    m1 = c + wx * wx * mc
    m2 = wz * s + wx * wy * mc
    m3 = -wy * s + wx * wz * mc
    m4 = -wz * s + wx * wy * mc
    m5 = c + wy * wy * mc
    m6 = wx * s + wy * wz * mc
    m7 = wy * s + wx * wz * mc
    m8 = -wx * s + wy * wz * mc
    m9 = c + wz * wz * mc
    x = px * m1 + py * m2 + pz * m3
    y = px * m4 + py * m5 + pz * m6
    z = px * m7 + py * m8 + pz * m9
  end subroutine rotate_one

end program rodrigues_main
