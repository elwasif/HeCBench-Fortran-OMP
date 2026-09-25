! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: par_t = 13_int32
  real(real64), parameter :: par_l = 2.65_real64
  real(real64), parameter :: theta = 45.0_real64
  real(real64), parameter :: tolerance = 1.0d-3

  character(len=256) :: arg0, arg1, arg2, arg3
  integer(int32) :: height, width, repeat
  integer(int64) :: n
  real(real64), allocatable :: h_filter(:), d_filter(:)
  real(real64) :: avg_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <height> <width> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) height
  read(arg2, *) width
  read(arg3, *) repeat

  if (height <= 0_int32 .or. width <= 0_int32 .or. repeat <= 0_int32) then
    print '(A)', 'FAIL'
    stop 1
  end if

  n = int(height, int64) * int(width, int64)
  allocate(h_filter(0:n - 1_int64), d_filter(0:n - 1_int64))

  call generate_gabor_kernel_host(height, width, h_filter)
  call generate_gabor_kernel_device(repeat, height, width, d_filter, avg_us)

  ok = arrays_match(h_filter, d_filter, tolerance)
  print '(A,F0.6,A)', 'Average kernel execution time: ', avg_us, ' (us)'
  print '(A)', merge('PASS', 'FAIL', ok)

  deallocate(h_filter, d_filter)

contains

  subroutine generate_gabor_kernel_host(height, width, gabor_spatial)
    integer(int32), intent(in) :: height, width
    real(real64), intent(out) :: gabor_spatial(0:)
    integer(int32) :: x, y
    integer(int64) :: idx
    real(real64) :: sx, sy, sx_2, sy_2, fx
    real(real64) :: ctheta, stheta, center_y, center_x, scale
    real(real64) :: centered_y, centered_x, u, v, pi

    pi = acos(-1.0_real64)
    sx = real(par_t, real64) / (2.0_real64 * sqrt(2.0_real64 * log(2.0_real64)))
    sy = par_l * sx
    sx_2 = sx * sx
    sy_2 = sy * sy
    fx = 1.0_real64 / real(par_t, real64)
    ctheta = cos(theta)
    stheta = sin(theta)
    center_y = real(height, real64) / 2.0_real64
    center_x = real(width, real64) / 2.0_real64
    scale = 1.0_real64 / (2.0_real64 * pi * sx * sy)

    do y = 0_int32, height - 1_int32
      centered_y = real(y, real64) - center_y
      do x = 0_int32, width - 1_int32
        idx = int(y, int64) * int(width, int64) + int(x, int64)
        centered_x = real(x, real64) - center_x
        u = ctheta * centered_x - stheta * centered_y
        v = ctheta * centered_y + stheta * centered_x
        gabor_spatial(idx) = scale * exp(-0.5_real64 * (u * u / sx_2 + v * v / sy_2)) * &
                             cos(2.0_real64 * pi * fx * u)
      end do
    end do
  end subroutine generate_gabor_kernel_host

  subroutine generate_gabor_kernel_device(repeat, height, width, gabor_spatial, avg_us)
    integer(int32), intent(in) :: repeat, height, width
    real(real64), intent(out) :: gabor_spatial(0:)
    real(real64), intent(out) :: avg_us
    integer(int32) :: iter
    integer(int64) :: idx, n
    integer(int32) :: x, y
    real(real64) :: sx, sy, sx_2, sy_2, fx
    real(real64) :: ctheta, stheta, center_y, center_x, scale
    real(real64) :: centered_y, centered_x, u, v, pi
    real(real64) :: start_time, end_time

    n = int(height, int64) * int(width, int64)
    pi = acos(-1.0_real64)
    sx = real(par_t, real64) / (2.0_real64 * sqrt(2.0_real64 * log(2.0_real64)))
    sy = par_l * sx
    sx_2 = sx * sx
    sy_2 = sy * sy
    fx = 1.0_real64 / real(par_t, real64)
    ctheta = cos(theta)
    stheta = sin(theta)
    center_y = real(height, real64) / 2.0_real64
    center_x = real(width, real64) / 2.0_real64
    scale = 1.0_real64 / (2.0_real64 * pi * sx * sy)

    !$omp target data map(from: gabor_spatial(0:n - 1_int64))
    start_time = omp_get_wtime()
    do iter = 1_int32, repeat
      !$omp target teams distribute parallel do thread_limit(256) &
      !$omp& private(x, y, centered_y, centered_x, u, v)
      do idx = 0_int64, n - 1_int64
        y = int(idx / int(width, int64), int32)
        x = int(mod(idx, int(width, int64)), int32)
        centered_y = real(y, real64) - center_y
        centered_x = real(x, real64) - center_x
        u = ctheta * centered_x - stheta * centered_y
        v = ctheta * centered_y + stheta * centered_x
        gabor_spatial(idx) = scale * exp(-0.5_real64 * (u * u / sx_2 + v * v / sy_2)) * &
                             cos(2.0_real64 * pi * fx * u)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
  end subroutine generate_gabor_kernel_device

  logical function arrays_match(lhs, rhs, tol) result(ok)
    real(real64), intent(in) :: lhs(0:), rhs(0:)
    real(real64), intent(in) :: tol
    integer(int64) :: i, n

    n = ubound(lhs, 1)
    ok = .true.
    do i = 0_int64, n
      if (abs(lhs(i) - rhs(i)) > tol) then
        ok = .false.
        exit
      end if
    end do
  end function arrays_match

end program main
