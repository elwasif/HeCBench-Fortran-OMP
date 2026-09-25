! SPDX-License-Identifier: CC0-1.0
program lombscargle_main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: x_shape = 1000
  integer, parameter :: freqs_shape = 100000
  real(real32), parameter :: a = 2.0_real32
  real(real32), parameter :: w = 1.0_real32
  real(real32), parameter :: phi = 1.57_real32
  real(real32), parameter :: y_dot = 2.0_real32 / 1.5_real32
  integer :: repeat, n, i
  real(real32), allocatable :: x(:), y(:), freqs(:), p(:), p2(:)
  real(real64) :: start_time, elapsed_us
  character(len=64) :: arg
  logical :: error

  if (command_argument_count() /= 1) then
    write(*,'("Usage: ./main <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) then
    write(*,'("Usage: ./main <repeat>")')
    stop 1
  end if

  allocate(x(0:x_shape-1), y(0:x_shape-1), freqs(0:freqs_shape-1))
  allocate(p(0:freqs_shape-1), p2(0:freqs_shape-1))

  do i = 0, x_shape - 1
    x(i) = 0.01_real32 + real(i, real32) * (31.4_real32 - 0.01_real32) / real(x_shape, real32)
  end do

  do i = 0, x_shape - 1
    y(i) = a * sin(w * x(i) + phi)
  end do

  do i = 0, freqs_shape - 1
    freqs(i) = 0.01_real32 + real(i, real32) * (10.0_real32 - 0.01_real32) / real(freqs_shape, real32)
  end do

  !$omp target data map(to: x, y, freqs) map(from: p)
  start_time = omp_get_wtime()
  do n = 1, repeat
    call lombscargle_kernel(x_shape, freqs_shape, x, y, freqs, p, y_dot)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  !$omp end target data

  write(*,'("Average kernel execution time ",F0.6," (us)")') elapsed_us

  call lombscargle_cpu(x_shape, freqs_shape, x, y, freqs, p2, y_dot)

  error = .false.
  do i = 0, freqs_shape - 1
    if (abs(p(i) - p2(i)) > 1.0e-1_real32) then
      write(*,'(F0.3,1X,F0.3)') p(i), p2(i)
      error = .true.
      exit
    end if
  end do

  if (error) then
    write(*,'("FAIL")')
    stop 1
  else
    write(*,'("PASS")')
  end if

contains

  subroutine lombscargle_kernel(x_count, freq_count, x, y, freqs, pgram, y_dot_value)
    integer, intent(in) :: x_count, freq_count
    real(real32), intent(in) :: x(0:), y(0:), freqs(0:), y_dot_value
    real(real32), intent(out) :: pgram(0:)
    integer :: tid, j
    real(real32) :: freq, xc, xs, cc, ss, cs
    real(real32) :: c, s, c_tau, s_tau, tau, c_tau2, s_tau2, cs_tau
    real(real32) :: c_term, s_term, c_denom, s_denom

    !$omp target teams distribute parallel do thread_limit(256) private(j, freq, xc, xs, cc, ss, cs, c, s) &
    !$omp& private(c_tau, s_tau, tau, c_tau2, s_tau2, cs_tau, c_term, s_term, c_denom, s_denom)
    do tid = 0, freq_count - 1
      freq = freqs(tid)
      xc = 0.0_real32
      xs = 0.0_real32
      cc = 0.0_real32
      ss = 0.0_real32
      cs = 0.0_real32

      do j = 0, x_count - 1
        s = sin(freq * x(j))
        c = cos(freq * x(j))
        xc = xc + y(j) * c
        xs = xs + y(j) * s
        cc = cc + c * c
        ss = ss + s * s
        cs = cs + c * s
      end do

      tau = atan2(2.0_real32 * cs, cc - ss) / (2.0_real32 * freq)
      s_tau = sin(freq * tau)
      c_tau = cos(freq * tau)
      c_tau2 = c_tau * c_tau
      s_tau2 = s_tau * s_tau
      cs_tau = 2.0_real32 * c_tau * s_tau
      c_term = c_tau * xc + s_tau * xs
      s_term = c_tau * xs - s_tau * xc
      c_denom = c_tau2 * cc + cs_tau * cs + s_tau2 * ss
      s_denom = c_tau2 * ss - cs_tau * cs + s_tau2 * cc
      pgram(tid) = 0.5_real32 * ((c_term * c_term / c_denom) + (s_term * s_term / s_denom)) * y_dot_value
    end do
    !$omp end target teams distribute parallel do
  end subroutine lombscargle_kernel

  subroutine lombscargle_cpu(x_count, freq_count, x, y, freqs, pgram, y_dot_value)
    integer, intent(in) :: x_count, freq_count
    real(real32), intent(in) :: x(0:), y(0:), freqs(0:), y_dot_value
    real(real32), intent(out) :: pgram(0:)
    integer :: tid, j
    real(real32) :: freq, xc, xs, cc, ss, cs
    real(real32) :: c, s, c_tau, s_tau, tau, c_tau2, s_tau2, cs_tau
    real(real32) :: c_term, s_term, c_denom, s_denom

    do tid = 0, freq_count - 1
      freq = freqs(tid)
      xc = 0.0_real32
      xs = 0.0_real32
      cc = 0.0_real32
      ss = 0.0_real32
      cs = 0.0_real32

      do j = 0, x_count - 1
        s = sin(freq * x(j))
        c = cos(freq * x(j))
        xc = xc + y(j) * c
        xs = xs + y(j) * s
        cc = cc + c * c
        ss = ss + s * s
        cs = cs + c * s
      end do

      tau = atan2(2.0_real32 * cs, cc - ss) / (2.0_real32 * freq)
      s_tau = sin(freq * tau)
      c_tau = cos(freq * tau)
      c_tau2 = c_tau * c_tau
      s_tau2 = s_tau * s_tau
      cs_tau = 2.0_real32 * c_tau * s_tau
      c_term = c_tau * xc + s_tau * xs
      s_term = c_tau * xs - s_tau * xc
      c_denom = c_tau2 * cc + cs_tau * cs + s_tau2 * ss
      s_denom = c_tau2 * ss - cs_tau * cs + s_tau2 * cc
      pgram(tid) = 0.5_real32 * ((c_term * c_term / c_denom) + (s_term * s_term / s_denom)) * y_dot_value
    end do
  end subroutine lombscargle_cpu

end program lombscargle_main
