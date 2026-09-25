! SPDX-License-Identifier: CC0-1.0
program popcount_main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
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
  end interface

  integer, parameter :: block_size = 256
  integer(int64), parameter :: m1 = int(z'5555555555555555', int64)
  integer(int64), parameter :: m2 = int(z'3333333333333333', int64)
  integer(int64), parameter :: m4 = int(z'0f0f0f0f0f0f0f0f', int64)
  integer(int64), parameter :: h01 = int(z'0101010101010101', int64)
  integer(int32), parameter :: lut(0:255) = [integer(int32) :: &
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, &
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, &
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, &
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7, &
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, &
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7, &
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7, &
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8]
  integer :: length, repeat
  integer(int64), allocatable :: data(:)
  integer(int32), allocatable :: result(:)
  integer :: i, n, bit_index
  real(real64) :: start_time, elapsed_us
  character(len=64) :: arg
  integer(int64) :: x
  integer(int32) :: count_value, i1, i2, i3, i4, i5, i6, i7, i8

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <length> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) length
  call get_command_argument(2, arg)
  read(arg, *) repeat

  allocate(data(0:length-1), result(0:length-1))

  call c_srand(2_c_int)
  do i = 0, length - 1
    data(i) = ior(ishft(int(c_rand(), int64), 32), int(c_rand(), int64))
  end do

  !$omp target data map(to: data) map(alloc: result)
  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size) private(x)
    do i = 0, length - 1
      x = data(i)
      x = x - iand(shiftr(x, 1), m1)
      x = iand(x, m2) + iand(shiftr(x, 2), m2)
      x = iand(x + shiftr(x, 4), m4)
      x = x + shiftr(x, 8)
      x = x + shiftr(x, 16)
      x = x + shiftr(x, 32)
      result(i) = int(iand(x, 127_int64), int32)
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc1): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)

  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size) private(x)
    do i = 0, length - 1
      x = data(i)
      x = x - iand(shiftr(x, 1), m1)
      x = iand(x, m2) + iand(shiftr(x, 2), m2)
      x = iand(x + shiftr(x, 4), m4)
      result(i) = int(shiftr(x * h01, 56), int32)
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc2): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)

  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size) private(x, count_value)
    do i = 0, length - 1
      count_value = 0_int32
      x = data(i)
      do while (x /= 0_int64)
        count_value = count_value + 1_int32
        x = iand(x, x - 1_int64)
      end do
      result(i) = count_value
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc3): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)

  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size) private(x, bit_index, count_value)
    do i = 0, length - 1
      x = data(i)
      count_value = 0_int32
      do bit_index = 0, 63
        count_value = count_value + int(iand(x, 1_int64), int32)
        x = shiftr(x, 1)
      end do
      result(i) = count_value
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc4): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)

  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size) private(x, i1, i2, i3, i4, i5, i6, i7, i8)
    do i = 0, length - 1
      x = data(i)
      i1 = lut(int(iand(x, 255_int64)))
      i2 = lut(int(iand(shiftr(x, 8), 255_int64)))
      i3 = lut(int(iand(shiftr(x, 16), 255_int64)))
      i4 = lut(int(iand(shiftr(x, 24), 255_int64)))
      i5 = lut(int(iand(shiftr(x, 32), 255_int64)))
      i6 = lut(int(iand(shiftr(x, 40), 255_int64)))
      i7 = lut(int(iand(shiftr(x, 48), 255_int64)))
      i8 = lut(int(iand(shiftr(x, 56), 255_int64)))
      result(i) = (i1 + i2) + (i3 + i4) + (i5 + i6) + (i7 + i8)
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc5): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)

  start_time = omp_get_wtime()
  do n = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, length - 1
      result(i) = int(popcnt(data(i)), int32)
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (pc6): ",F0.6," (us)")') elapsed_us
  !$omp target update from(result)
  call check_results(data, result, length)
  !$omp end target data

contains

  subroutine check_results(data, result, length)
    integer(int64), intent(in) :: data(0:)
    integer(int32), intent(in) :: result(0:)
    integer, intent(in) :: length
    integer :: i
    logical :: error

    error = .false.
    do i = 0, length - 1
      if (popcnt(data(i)) /= result(i)) then
        error = .true.
        exit
      end if
    end do

    if (error) then
      write(*,'("Fail")')
    else
      write(*,'("Success")')
    end if
  end subroutine check_results

end program popcount_main
