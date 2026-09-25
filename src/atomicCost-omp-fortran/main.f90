! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: block_size = 256_int32
  integer(int32), parameter :: length = 922521600_int32
  character(len=256) :: arg0, arg1, arg2
  integer(int32) :: nelems, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <N> <repeat>'
    print '(A)', 'N: the number of elements to sum per thread (1 - 16)'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) nelems
  read(arg2, *) repeat

  if (nelems <= 0_int32 .or. nelems > 16_int32) stop 1
  if (repeat <= 0_int32) stop 1
  if (mod(length, block_size) /= 0_int32) stop 1

  print '(A)', ''
  print '(A)', 'FP64 atomic add'
  call atomic_cost_real64(length, nelems, repeat)

  print '(A)', ''
  print '(A)', 'INT32 atomic add'
  call atomic_cost_int32(length, nelems, repeat)

  print '(A)', ''
  print '(A)', 'FP32 atomic add'
  call atomic_cost_real32(length, nelems, repeat)

contains

  subroutine atomic_cost_real64(total_length, chunk_size, repeat)
    integer(int32), intent(in) :: total_length, chunk_size, repeat
    integer(int32) :: num_threads, iter
    real(real64), allocatable :: result_wi(:), result_wo(:)
    real(real64) :: start_time, end_time, avg_us

    call validate_shape(total_length, chunk_size, num_threads)
    allocate(result_wi(num_threads), result_wo(num_threads))
    result_wi = 0.0_real64
    result_wo = 0.0_real64

    print '(A)', ''
    print '(A)', ''
    print '(A,I0,A)', 'Each thread sums up ', chunk_size, ' elements'

    !$omp target data map(tofrom: result_wi(1:num_threads), result_wo(1:num_threads))
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call with_atomic_real64(result_wi, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wi(1:num_threads))

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call without_atomic_real64(result_wo, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithoutAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wo(1:num_threads))
    !$omp end target data

    call print_result(all(result_wi == result_wo))
    deallocate(result_wi, result_wo)
  end subroutine atomic_cost_real64

  subroutine atomic_cost_int32(total_length, chunk_size, repeat)
    integer(int32), intent(in) :: total_length, chunk_size, repeat
    integer(int32) :: num_threads, iter
    integer(int32), allocatable :: result_wi(:), result_wo(:)
    real(real64) :: start_time, end_time, avg_us

    call validate_shape(total_length, chunk_size, num_threads)
    allocate(result_wi(num_threads), result_wo(num_threads))
    result_wi = 0_int32
    result_wo = 0_int32

    print '(A)', ''
    print '(A)', ''
    print '(A,I0,A)', 'Each thread sums up ', chunk_size, ' elements'

    !$omp target data map(tofrom: result_wi(1:num_threads), result_wo(1:num_threads))
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call with_atomic_int32(result_wi, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wi(1:num_threads))

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call without_atomic_int32(result_wo, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithoutAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wo(1:num_threads))
    !$omp end target data

    call print_result(all(result_wi == result_wo))
    deallocate(result_wi, result_wo)
  end subroutine atomic_cost_int32

  subroutine atomic_cost_real32(total_length, chunk_size, repeat)
    integer(int32), intent(in) :: total_length, chunk_size, repeat
    integer(int32) :: num_threads, iter
    real(real32), allocatable :: result_wi(:), result_wo(:)
    real(real64) :: start_time, end_time, avg_us

    call validate_shape(total_length, chunk_size, num_threads)
    allocate(result_wi(num_threads), result_wo(num_threads))
    result_wi = 0.0_real32
    result_wo = 0.0_real32

    print '(A)', ''
    print '(A)', ''
    print '(A,I0,A)', 'Each thread sums up ', chunk_size, ' elements'

    !$omp target data map(tofrom: result_wi(1:num_threads), result_wo(1:num_threads))
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call with_atomic_real32(result_wi, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wi(1:num_threads))

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call without_atomic_real32(result_wo, chunk_size, num_threads)
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of WithoutAtomicOnGlobalMem: ', avg_us, ' (us)'
    !$omp target update from(result_wo(1:num_threads))
    !$omp end target data

    call print_result(all(result_wi == result_wo))
    deallocate(result_wi, result_wo)
  end subroutine atomic_cost_real32

  subroutine validate_shape(total_length, chunk_size, num_threads)
    integer(int32), intent(in) :: total_length, chunk_size
    integer(int32), intent(out) :: num_threads

    if (mod(total_length, chunk_size) /= 0_int32) stop 1
    num_threads = total_length / chunk_size
    if (mod(num_threads, block_size) /= 0_int32) stop 1
  end subroutine validate_shape

  subroutine with_atomic_real64(result, chunk_size, num_threads)
    real(real64), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        !$omp atomic update
        result(tid) = result(tid) + real(mod(i, 2_int32), real64)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine with_atomic_real64

  subroutine without_atomic_real64(result, chunk_size, num_threads)
    real(real64), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        result(tid) = result(tid) + real(mod(i, 2_int32), real64)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine without_atomic_real64

  subroutine with_atomic_int32(result, chunk_size, num_threads)
    integer(int32), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        !$omp atomic update
        result(tid) = result(tid) + mod(i, 2_int32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine with_atomic_int32

  subroutine without_atomic_int32(result, chunk_size, num_threads)
    integer(int32), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        result(tid) = result(tid) + mod(i, 2_int32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine without_atomic_int32

  subroutine with_atomic_real32(result, chunk_size, num_threads)
    real(real32), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        !$omp atomic update
        result(tid) = result(tid) + real(mod(i, 2_int32), real32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine with_atomic_real32

  subroutine without_atomic_real32(result, chunk_size, num_threads)
    real(real32), intent(inout) :: result(:)
    integer(int32), intent(in) :: chunk_size, num_threads
    integer(int32) :: tid, i

    !$omp target teams distribute parallel do thread_limit(block_size) private(i)
    do tid = 1, num_threads
      do i = (tid - 1_int32) * chunk_size, tid * chunk_size - 1_int32
        result(tid) = result(tid) + real(mod(i, 2_int32), real32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine without_atomic_real32

  subroutine print_result(ok)
    logical, intent(in) :: ok

    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine print_result

end program main
