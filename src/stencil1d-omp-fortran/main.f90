! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: radius = 7_int32
  integer(int32), parameter :: block_size = 256_int32

  character(len=256) :: arg0, arg1, arg2
  integer(int32) :: length, repeat, pad_size
  integer(int32), allocatable :: a(:), b(:)
  integer(int32) :: i, j, expected
  logical :: ok
  real(real64) :: start_time, end_time, avg_s

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <length> <repeat>'
    print '(A,I0)', 'length is a multiple of ', block_size
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) length
  read(arg2, *) repeat
  if (length <= 0_int32 .or. repeat <= 0_int32) stop 1

  pad_size = length + radius
  allocate(a(0:pad_size - 1), b(0:length - 1))

  do i = 0, pad_size - 1
    a(i) = i
  end do
  b = 0_int32

  start_time = omp_get_wtime()
  do i = 1, repeat
    call stencil_kernel(length, pad_size, a, b)
  end do
  end_time = omp_get_wtime()

  avg_s = (end_time - start_time) / real(repeat, real64)
  print '(A,F8.6,A)', 'Average kernel execution time: ', avg_s, ' (s)'

  ok = .true.
  do i = 0, (2_int32 * radius) - 1
    expected = 0_int32
    do j = i, i + 2_int32 * radius
      if (j >= radius) expected = expected + (a(j) - radius)
    end do
    if (expected /= b(i)) then
      print '(A,I0,A,I0,A,I0,A)', 'Error at ', i, ': ', expected, ' (host) != ', b(i), ' (device)'
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    do i = 2_int32 * radius, length - 1
      expected = 0_int32
      do j = i - radius, i + radius
        expected = expected + a(j)
      end do
      if (expected /= b(i)) then
        print '(A,I0,A,I0,A,I0,A)', 'Error at ', i, ': ', expected, ' (host) != ', b(i), ' (device)'
        ok = .false.
        exit
      end if
    end do
  end if

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(a, b)

contains

  subroutine stencil_kernel(n, padded_n, a, b)
    integer(int32), intent(in) :: n, padded_n
    integer(int32), intent(in) :: a(0:padded_n - 1)
    integer(int32), intent(out) :: b(0:n - 1)
    integer(int32) :: i, j, gindex, offset, result
    integer(int32) :: temp(0:block_size + 2_int32 * radius - 1)

    !$omp target teams distribute map(to: a(0:padded_n - 1)) map(from: b(0:n - 1)) &
    !$omp& private(temp, j, gindex, offset, result)
    do i = 0, n - 1, block_size
      !$omp parallel do schedule(static,1) private(gindex)
      do j = 0, block_size - 1
        gindex = i + j
        temp(j + radius) = a(gindex)
        if (j < radius) then
          if (gindex < radius) then
            temp(j) = 0_int32
          else
            temp(j) = a(gindex - radius)
          end if
          temp(j + radius + block_size) = a(gindex + block_size)
        end if
      end do
      !$omp end parallel do

      !$omp parallel do schedule(static,1) private(result, offset)
      do j = 0, block_size - 1
        result = 0_int32
        do offset = -radius, radius
          result = result + temp(j + radius + offset)
        end do
        b(i + j) = result
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine stencil_kernel

end program main
