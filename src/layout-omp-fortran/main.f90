! SPDX-License-Identifier: CC0-1.0
program layout_main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: tree_num = 4096
  integer, parameter :: tree_size = 4096
  integer, parameter :: group_size = 256
  integer :: iterations, elements, i, j, n
  integer(int32), allocatable :: data(:), output(:), reference(:)
  integer(int64) :: expected
  real(real64) :: start_time, elapsed_us
  character(len=64) :: arg
  logical :: fail

  if (command_argument_count() /= 1) then
    write(*,'("Usage: ./main <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) iterations
  if (iterations < 1) then
    write(*,'("Iterations cannot be 0 or negative. Exiting..")')
    stop 1
  end if

  elements = tree_size * tree_num
  allocate(data(0:elements-1), output(0:tree_num-1), reference(0:tree_num-1))

  do i = 0, tree_num - 1
    expected = 0_int64
    do j = 0, tree_size - 1
      expected = expected + int(i * tree_size + j, int64)
    end do
    reference(i) = int(expected, int32)
  end do

  !$omp target data map(alloc: data, output)
  do i = 0, tree_num - 1
    do j = 0, tree_size - 1
      data(j + i * tree_size) = int(j + i * tree_size, int32)
    end do
  end do
  !$omp target update to(data)

  start_time = omp_get_wtime()
  do n = 1, iterations
    call run_aos(data, output)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(iterations, real64)
  write(*,'("Average kernel execution time (AoS): ",F0.6," (us)")') elapsed_us

  !$omp target update from(output)
  fail = any(output /= reference)
  if (fail) then
    write(*,'("FAIL")')
  else
    write(*,'("PASS")')
  end if

  do i = 0, tree_num - 1
    do j = 0, tree_size - 1
      data(i + j * tree_num) = int(j + i * tree_size, int32)
    end do
  end do
  !$omp target update to(data)

  start_time = omp_get_wtime()
  do n = 1, iterations
    call run_soa(data, output)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(iterations, real64)
  write(*,'("Average kernel execution time (SoA): ",F0.6," (us)")') elapsed_us

  !$omp target update from(output)
  !$omp end target data
  fail = any(output /= reference)
  if (fail) then
    write(*,'("FAIL")')
  else
    write(*,'("PASS")')
  end if

contains

  subroutine run_aos(data, output)
    integer(int32), intent(in) :: data(0:)
    integer(int32), intent(out) :: output(0:)
    integer :: gid, idx
    integer(int32) :: result

    !$omp target teams distribute parallel do thread_limit(group_size) private(idx, result)
    do gid = 0, tree_num - 1
      result = 0_int32
      do idx = 0, tree_size - 1
        result = result + data(idx + gid * tree_size)
      end do
      output(gid) = result
    end do
    !$omp end target teams distribute parallel do
  end subroutine run_aos

  subroutine run_soa(data, output)
    integer(int32), intent(in) :: data(0:)
    integer(int32), intent(out) :: output(0:)
    integer :: gid, idx
    integer(int32) :: result

    !$omp target teams distribute parallel do thread_limit(group_size) private(idx, result)
    do gid = 0, tree_num - 1
      result = 0_int32
      do idx = 0, tree_size - 1
        result = result + data(gid + idx * tree_num)
      end do
      output(gid) = result
    end do
    !$omp end target teams distribute parallel do
  end subroutine run_soa

end program layout_main
