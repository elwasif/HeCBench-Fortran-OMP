! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand

    subroutine c_exit(status) bind(C, name="exit")
      import :: c_int
      integer(c_int), value :: status
    end subroutine c_exit
  end interface

  character(len=256) :: arg0, arg1, arg2
  integer :: length, repeat, i, ios
  integer, allocatable :: y(:)
  real(real32), allocatable :: x1(:), x2(:), dout(:), dx1(:), dx2(:), rdx1(:), rdx2(:)
  real(real32), parameter :: margin = 0.01_real32
  real(real32), parameter :: rand_max = 2147483647.0_real32
  real(real64) :: start_time, end_time, elapsed_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    call c_exit(1_c_int)
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *, iostat=ios) length
  if (ios /= 0) length = 0
  read(arg2, *, iostat=ios) repeat
  if (ios /= 0) repeat = 0

  allocate(y(length), x1(length), x2(length), dout(length), dx1(length), dx2(length), rdx1(length), rdx2(length))

  call c_srand(123_c_int)
  do i = 1, length
    x1(i) = random_uniform_signed()
    x2(i) = random_uniform_signed()
    dout(i) = random_uniform_signed()
    if (random_uniform_signed() < 0.0_real32) then
      y(i) = -1
    else
      y(i) = 1
    end if
  end do
  dx1 = 0.0_real32
  dx2 = 0.0_real32

  !$omp target data map(to: x1(1:length), x2(1:length), dout(1:length), y(1:length)) &
  !$omp& map(from: dx1(1:length), dx2(1:length))
  do i = 1, repeat
    call mrc_gradient(length, y, x1, x2, dout, margin, dx1, dx2)
    call mrc_gradient2(length, y, x1, x2, dout, margin, dx1, dx2)
  end do

  start_time = omp_get_wtime()
  do i = 1, repeat
    call mrc_gradient(length, y, x1, x2, dout, margin, dx1, dx2)
  end do
  end_time = omp_get_wtime()
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time of MRC kernel: ', elapsed_us, ' (us)'

  start_time = omp_get_wtime()
  do i = 1, repeat
    call mrc_gradient2(length, y, x1, x2, dout, margin, dx1, dx2)
  end do
  end_time = omp_get_wtime()
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time of MRC2 kernel: ', elapsed_us, ' (us)'
  !$omp end target data

  call reference(length, y, x1, x2, dout, margin, rdx1, rdx2)

  ok = .true.
  do i = 1, length
    if (abs(dx1(i) - rdx1(i)) > 1.0e-3_real32 .or. abs(dx2(i) - rdx2(i)) > 1.0e-3_real32) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(y, x1, x2, dout, dx1, dx2, rdx1, rdx2)

contains

  function random_uniform_signed() result(value)
    real(real32) :: value

    value = real(c_rand(), real32) / rand_max
    value = value * 4.0_real32 - 2.0_real32
  end function random_uniform_signed

  subroutine mrc_gradient(n, y, x1, x2, dout, margin, dx1, dx2)
    integer, intent(in) :: n
    integer, intent(in) :: y(:)
    real(real32), intent(in) :: x1(:), x2(:), dout(:), margin
    real(real32), intent(out) :: dx1(:), dx2(:)
    integer :: idx
    real(real32) :: dist

    !$omp target teams distribute parallel do thread_limit(256) private(dist)
    do idx = 1, n
      dist = -real(y(idx), real32) * (x1(idx) - x2(idx)) + margin
      if (dist < 0.0_real32) then
        dx1(idx) = 0.0_real32
        dx2(idx) = 0.0_real32
      else
        dx1(idx) = -real(y(idx), real32) * dout(idx)
        dx2(idx) = real(y(idx), real32) * dout(idx)
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine mrc_gradient

  subroutine mrc_gradient2(n, y, x1, x2, dout, margin, dx1, dx2)
    integer, intent(in) :: n
    integer, intent(in) :: y(:)
    real(real32), intent(in) :: x1(:), x2(:), dout(:), margin
    real(real32), intent(out) :: dx1(:), dx2(:)
    integer :: idx
    real(real32) :: yval, oval, dist

    !$omp target teams distribute parallel do thread_limit(256) private(yval, oval, dist)
    do idx = 1, n
      yval = real(y(idx), real32)
      oval = dout(idx)
      dist = -yval * (x1(idx) - x2(idx)) + margin
      if (dist < 0.0_real32) then
        dx1(idx) = 0.0_real32
        dx2(idx) = 0.0_real32
      else
        dx1(idx) = -yval * oval
        dx2(idx) = yval * oval
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine mrc_gradient2

  subroutine reference(n, y, x1, x2, dout, margin, dx1, dx2)
    integer, intent(in) :: n
    integer, intent(in) :: y(:)
    real(real32), intent(in) :: x1(:), x2(:), dout(:), margin
    real(real32), intent(out) :: dx1(:), dx2(:)
    integer :: idx
    real(real32) :: dist

    do idx = 1, n
      dist = -real(y(idx), real32) * (x1(idx) - x2(idx)) + margin
      if (dist < 0.0_real32) then
        dx1(idx) = 0.0_real32
        dx2(idx) = 0.0_real32
      else
        dx1(idx) = -real(y(idx), real32) * dout(idx)
        dx2(idx) = real(y(idx), real32) * dout(idx)
      end if
    end do
  end subroutine reference

end program main
