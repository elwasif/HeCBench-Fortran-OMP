! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: repeat = 1000
  integer, parameter :: n = 370
  integer, parameter :: ldat = n
  integer, parameter :: incx = 1
  integer, parameter :: incy = 1
  integer, parameter :: at_size = n * ldat
  integer, parameter :: x_size = n * incx
  integer, parameter :: y_size = n * incy
  real(real32), parameter :: alpha_re = 3.14_real32
  real(real32), parameter :: alpha_im = 1.59_real32
  real(real32), parameter :: beta_re = 2.71_real32
  real(real32), parameter :: beta_im = 8.28_real32

  real(real32) :: AT(0:2*at_size-1)
  real(real32) :: X(0:2*x_size-1)
  real(real32) :: Y_cpu(0:2*y_size-1)
  real(real32) :: Y_gpu(0:2*y_size-1)
  integer :: i, j
  real(real64) :: start_time, elapsed_us
  logical :: ok

  do i = 0, n - 1
    X(re(i * incx + 0)) = real(i + 5, real32)
    X(im(i * incx + 0)) = real(i * 2, real32)
    Y_cpu(re(i * incy + 0)) = real(i * 3, real32)
    Y_cpu(im(i * incy + 0)) = real(i + 7, real32)
    Y_gpu(re(i * incy + 0)) = real(i * 3, real32)
    Y_gpu(im(i * incy + 0)) = real(i + 7, real32)
    do j = 0, ldat - 1
      AT(re(i * ldat + j)) = real(i + j, real32)
      AT(im(i * ldat + j)) = real(i + 3, real32)
    end do
  end do

  call chemv_cpu(alpha_re, alpha_im, beta_re, beta_im, AT, X, Y_cpu)

  !$omp target data map(to: AT, X) map(tofrom: Y_gpu)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call chemv_kernel0(AT, X, Y_gpu, alpha_im, alpha_re, beta_im, beta_re)
    call chemv_kernel1(AT, X, Y_gpu, alpha_im, alpha_re)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average execution time of chemv kernels: ', elapsed_us, ' (us)'
  !$omp end target data

  ok = .true.
  do i = 0, n - 1
    if ((abs(Y_cpu(re(i * incy + 0)) - Y_gpu(re(i * incy + 0))) > 1.0e-3_real32) .or. &
        (abs(Y_cpu(im(i * incy + 0)) - Y_gpu(im(i * incy + 0))) > 1.0e-3_real32)) then
      ok = .false.
      write(*,'(I0,1X,F0.6,1X,F0.6)') i, Y_cpu(re(i * incy + 0)), Y_gpu(re(i * incy + 0))
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
    stop 1
  end if

contains

  pure integer function re(idx)
    integer, intent(in) :: idx
    re = 2 * idx
  end function re

  pure integer function im(idx)
    integer, intent(in) :: idx
    im = 2 * idx + 1
  end function im

  subroutine chemv_cpu(alpha_re, alpha_im, beta_re, beta_im, AT, X, Y)
    real(real32), intent(in) :: alpha_re, alpha_im, beta_re, beta_im
    real(real32), intent(in) :: AT(0:2*at_size-1), X(0:2*x_size-1)
    real(real32), intent(inout) :: Y(0:2*y_size-1)
    integer :: i0, i1, i2, i3
    real(real32) :: var5_Re, var5_Im
    real(real32) :: var2_Re, var3_Im, var2_Im, var4_Im, var4_Re, var3_Re
    real(real32) :: var99_Re, var96_Re, var98_Im, var96_Im, var94_Im, var95_Im
    real(real32) :: var94_Re, var95_Re, var97_Im, var99_Im, var97_Re, var98_Re

    do i0 = 0, n - 1
      var5_Re = (Y(re(i0 * incy + 0)) * beta_re) - (Y(im(i0 * incy + 0)) * beta_im)
      var5_Im = (Y(im(i0 * incy + 0)) * beta_re) + (Y(re(i0 * incy + 0)) * beta_im)
      Y(re(i0 * incy + 0)) = var5_Re
      Y(im(i0 * incy + 0)) = var5_Im
    end do

    do i1 = 0, n - 1
      var2_Re = alpha_re * AT(re(i1 * ldat + i1))
      var2_Im = alpha_im * AT(re(i1 * ldat + i1))
      var3_Re = (var2_Re * X(re(i1 * incx + 0))) - (var2_Im * X(im(i1 * incx + 0)))
      var3_Im = (var2_Im * X(re(i1 * incx + 0))) + (var2_Re * X(im(i1 * incx + 0)))
      var4_Re = Y(re(i1 * incy + 0)) + var3_Re
      var4_Im = Y(im(i1 * incy + 0)) + var3_Im
      Y(re(i1 * incy + 0)) = var4_Re
      Y(im(i1 * incy + 0)) = var4_Im
    end do

    do i2 = 0, n - 2
      do i3 = 0, (n - 1) - (1 + i2)
        var94_Re = (alpha_re * AT(re(i2 * ldat + ((1 + i2) + i3)))) - &
            (alpha_im * (-AT(im(i2 * ldat + ((1 + i2) + i3)))))
        var94_Im = (alpha_im * AT(re(i2 * ldat + ((1 + i2) + i3)))) + &
            (alpha_re * (-AT(im(i2 * ldat + ((1 + i2) + i3)))))
        var95_Re = (var94_Re * X(re(((i3 + i2) + 1) * incx + 0))) - &
            (var94_Im * X(im(((i3 + i2) + 1) * incx + 0)))
        var95_Im = (var94_Im * X(re(((i3 + i2) + 1) * incx + 0))) + &
            (var94_Re * X(im(((i3 + i2) + 1) * incx + 0)))
        var96_Re = Y(re(i2 * incy + 0)) + var95_Re
        var96_Im = Y(im(i2 * incy + 0)) + var95_Im
        Y(re(i2 * incy + 0)) = var96_Re
        Y(im(i2 * incy + 0)) = var96_Im

        var97_Re = (alpha_re * AT(re(i2 * ldat + ((1 + i2) + i3)))) - &
            (alpha_im * AT(im(i2 * ldat + ((1 + i2) + i3))))
        var97_Im = (alpha_im * AT(re(i2 * ldat + ((1 + i2) + i3)))) + &
            (alpha_re * AT(im(i2 * ldat + ((1 + i2) + i3))))
        var98_Re = (var97_Re * X(re(i2 * incx + 0))) - (var97_Im * X(im(i2 * incx + 0)))
        var98_Im = (var97_Im * X(re(i2 * incx + 0))) + (var97_Re * X(im(i2 * incx + 0)))
        var99_Re = Y(re(((i3 + i2) + 1) * incy + 0)) + var98_Re
        var99_Im = Y(im(((i3 + i2) + 1) * incy + 0)) + var98_Im
        Y(re(((i3 + i2) + 1) * incy + 0)) = var99_Re
        Y(im(((i3 + i2) + 1) * incy + 0)) = var99_Im
      end do
    end do
  end subroutine chemv_cpu

  subroutine chemv_kernel0(AT, X, Y, alpha_im, alpha_re, beta_im, beta_re)
    real(real32), intent(in) :: AT(0:2*at_size-1), X(0:2*x_size-1)
    real(real32), intent(inout) :: Y(0:2*y_size-1)
    real(real32), intent(in) :: alpha_im, alpha_re, beta_im, beta_re
    integer :: gid, b0, t0, c1, c3
    real(real32) :: private_var5_Re, private_var5_Im, private_var2_Re, private_var3_Im
    real(real32) :: private_var2_Im, private_var4_Im, private_var4_Re, private_var3_Re
    real(real32) :: private_var99_Re, private_var98_Im, private_var97_Im
    real(real32) :: private_var99_Im, private_var97_Re, private_var98_Re

    !$omp target teams distribute parallel do num_teams(12) thread_limit(32) private(b0,t0,c1,c3) &
    !$omp& private(private_var5_Re,private_var5_Im,private_var2_Re,private_var3_Im,private_var2_Im) &
    !$omp& private(private_var4_Im,private_var4_Re,private_var3_Re,private_var99_Re,private_var98_Im) &
    !$omp& private(private_var97_Im,private_var99_Im,private_var97_Re,private_var98_Re)
    do gid = 0, 383
    b0 = gid / 32
    t0 = mod(gid, 32)
    do c1 = 0, min(368, 32 * b0 + 30), 32
      if (32 * b0 + t0 <= 369 .and. c1 == 0) then
        private_var5_Re = (Y(2 * (32 * b0 + t0)) * beta_re) - &
            (Y(2 * (32 * b0 + t0) + 1) * beta_im)
        private_var5_Im = (Y(2 * (32 * b0 + t0) + 1) * beta_re) + &
            (Y(2 * (32 * b0 + t0)) * beta_im)
        Y(2 * (32 * b0 + t0)) = private_var5_Re
        Y(2 * (32 * b0 + t0) + 1) = private_var5_Im
        private_var2_Re = alpha_re * AT(2 * (11872 * b0 + 371 * t0))
        private_var2_Im = alpha_im * AT(2 * (11872 * b0 + 371 * t0))
        private_var3_Re = (private_var2_Re * X(2 * (32 * b0 + t0))) - &
            (private_var2_Im * X(2 * (32 * b0 + t0) + 1))
        private_var3_Im = (private_var2_Im * X(2 * (32 * b0 + t0))) + &
            (private_var2_Re * X(2 * (32 * b0 + t0) + 1))
        private_var4_Re = Y(2 * (32 * b0 + t0)) + private_var3_Re
        private_var4_Im = Y(2 * (32 * b0 + t0) + 1) + private_var3_Im
        Y(2 * (32 * b0 + t0)) = private_var4_Re
        Y(2 * (32 * b0 + t0) + 1) = private_var4_Im
      end if
      if (32 * b0 + t0 <= 369) then
        do c3 = 0, min(31, 32 * b0 + t0 - c1 - 1)
          private_var97_Re = (alpha_re * AT(2 * (32 * b0 + t0 + 370 * c1 + 370 * c3))) - &
              (alpha_im * AT(2 * (32 * b0 + t0 + 370 * c1 + 370 * c3) + 1))
          private_var97_Im = (alpha_im * AT(2 * (32 * b0 + t0 + 370 * c1 + 370 * c3))) + &
              (alpha_re * AT(2 * (32 * b0 + t0 + 370 * c1 + 370 * c3) + 1))
          private_var98_Re = (private_var97_Re * X(2 * (c1 + c3))) - &
              (private_var97_Im * X(2 * (c1 + c3) + 1))
          private_var98_Im = (private_var97_Im * X(2 * (c1 + c3))) + &
              (private_var97_Re * X(2 * (c1 + c3) + 1))
          private_var99_Re = Y(2 * (32 * b0 + t0)) + private_var98_Re
          private_var99_Im = Y(2 * (32 * b0 + t0) + 1) + private_var98_Im
          Y(2 * (32 * b0 + t0)) = private_var99_Re
          Y(2 * (32 * b0 + t0) + 1) = private_var99_Im
        end do
      end if
    end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine chemv_kernel0

  subroutine chemv_kernel1(AT, X, Y, alpha_im, alpha_re)
    real(real32), intent(in) :: AT(0:2*at_size-1), X(0:2*x_size-1)
    real(real32), intent(inout) :: Y(0:2*y_size-1)
    real(real32), intent(in) :: alpha_im, alpha_re
    integer :: gid, b0, t0, c1, c3
    real(real32) :: private_var96_Re, private_var96_Im, private_var94_Im
    real(real32) :: private_var95_Im, private_var94_Re, private_var95_Re

    !$omp target teams distribute parallel do num_teams(12) thread_limit(32) private(b0,t0,c1,c3) &
    !$omp& private(private_var96_Re,private_var96_Im,private_var94_Im,private_var95_Im) &
    !$omp& private(private_var94_Re,private_var95_Re)
    do gid = 0, 383
    b0 = gid / 32
    t0 = mod(gid, 32)
    do c1 = 5888 * b0, min(67712, 5856 * b0 + 6016), 32
      do c3 = max(0, 5888 * b0 + 184 * t0 - c1), &
          min(31, 5856 * b0 + 183 * t0 - c1 + 368)
        private_var94_Re = (alpha_re * AT(2 * (5984 * b0 + 187 * t0 + c1 + c3 + 1))) - &
            (alpha_im * (-AT(2 * (5984 * b0 + 187 * t0 + c1 + c3 + 1) + 1)))
        private_var94_Im = (alpha_im * AT(2 * (5984 * b0 + 187 * t0 + c1 + c3 + 1))) + &
            (alpha_re * (-AT(2 * (5984 * b0 + 187 * t0 + c1 + c3 + 1) + 1)))
        private_var95_Re = (private_var94_Re * X(2 * (-5856 * b0 - 183 * t0 + c1 + c3 + 1))) - &
            (private_var94_Im * X(2 * (-5856 * b0 - 183 * t0 + c1 + c3 + 1) + 1))
        private_var95_Im = (private_var94_Im * X(2 * (-5856 * b0 - 183 * t0 + c1 + c3 + 1))) + &
            (private_var94_Re * X(2 * (-5856 * b0 - 183 * t0 + c1 + c3 + 1) + 1))
        private_var96_Re = Y(2 * (32 * b0 + t0)) + private_var95_Re
        private_var96_Im = Y(2 * (32 * b0 + t0) + 1) + private_var95_Im
        Y(2 * (32 * b0 + t0)) = private_var96_Re
        Y(2 * (32 * b0 + t0) + 1) = private_var96_Im
      end do
    end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine chemv_kernel1

end program main
