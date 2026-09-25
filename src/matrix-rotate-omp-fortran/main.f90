! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer(int32) :: n, repeat, i, j
  real(real32), allocatable :: serial_res(:), parallel_res(:)
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <matrix size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) n
  read(arg2, *) repeat
  if (n <= 0_int32 .or. repeat <= 0_int32) stop 1

  allocate(serial_res(n * n), parallel_res(n * n))

  do i = 0, n - 1
    do j = 0, n - 1
      serial_res(index_2d(i, j, n)) = real(i * n + j, real32)
      parallel_res(index_2d(i, j, n)) = real(i * n + j, real32)
    end do
  end do

  do i = 1, repeat
    call rotate_matrix_serial(serial_res, n)
  end do

  call rotate_matrix_parallel(parallel_res, n, repeat)

  ok = all(serial_res == parallel_res)
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(serial_res, parallel_res)

contains

  integer pure function index_2d(row, col, nsize) result(idx)
    integer(int32), intent(in) :: row, col, nsize
    idx = row * nsize + col + 1_int32
  end function index_2d

  subroutine rotate_matrix_parallel(matrix, nsize, repeat_count)
    real(real32), intent(inout) :: matrix(:)
    integer(int32), intent(in) :: nsize, repeat_count
    integer(int32) :: iter, layer, first, last, inner, offset
    real(real32) :: saved_top
    real(real64) :: start_time, end_time, avg_s

    !$omp target data map(tofrom: matrix(1:nsize*nsize))
    start_time = omp_get_wtime()
    do iter = 1, repeat_count
      !$omp target teams distribute parallel do thread_limit(256) private(first, last, inner, offset, saved_top)
      do layer = 0, (nsize / 2) - 1
        first = layer
        last = nsize - 1 - layer
        do inner = first, last - 1
          offset = inner - first
          saved_top = matrix(index_2d(first, inner, nsize))
          matrix(index_2d(first, inner, nsize)) = matrix(index_2d(last - offset, first, nsize))
          matrix(index_2d(last - offset, first, nsize)) = matrix(index_2d(last, last - offset, nsize))
          matrix(index_2d(last, last - offset, nsize)) = matrix(index_2d(inner, last, nsize))
          matrix(index_2d(inner, last, nsize)) = saved_top
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    avg_s = (end_time - start_time) / real(repeat_count, real64)
    print '(A,F8.6,A)', 'Average kernel execution time: ', avg_s, ' (s)'
  end subroutine rotate_matrix_parallel

  subroutine rotate_matrix_serial(matrix, nsize)
    real(real32), intent(inout) :: matrix(:)
    integer(int32), intent(in) :: nsize
    integer(int32) :: layer, first, last, inner, offset
    real(real32) :: saved_top

    do layer = 0, (nsize / 2) - 1
      first = layer
      last = nsize - 1 - layer
      do inner = first, last - 1
        offset = inner - first
        saved_top = matrix(index_2d(first, inner, nsize))
        matrix(index_2d(first, inner, nsize)) = matrix(index_2d(last - offset, first, nsize))
        matrix(index_2d(last - offset, first, nsize)) = matrix(index_2d(last, last - offset, nsize))
        matrix(index_2d(last, last - offset, nsize)) = matrix(index_2d(inner, last, nsize))
        matrix(index_2d(inner, last, nsize)) = saved_top
      end do
    end do
  end subroutine rotate_matrix_serial

end program main
