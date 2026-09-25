! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: block_size = 256_int32

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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer(int32) :: num_slice, slice_size, repeat, num_elem
  real(real32), allocatable :: input(:), output_gpu(:), output_cpu(:)
  integer(int32) :: i
  logical :: ok
  real(real64) :: start_time, end_time, avg_ms

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of slices> <slice size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) num_slice
  read(arg2, *) slice_size
  read(arg3, *) repeat

  if (num_slice <= 0_int32 .or. slice_size <= 0_int32 .or. repeat <= 0_int32) stop 1
  num_elem = num_slice * slice_size

  allocate(input(num_elem), output_gpu(num_elem), output_cpu(num_elem))

  call c_srand(2_c_int)
  do i = 1, num_elem
    input(i) = real(modulo(c_rand(), 13_c_int), real32)
  end do
  output_gpu = 0.0_real32
  output_cpu = 0.0_real32

  !$omp target data map(to: input(1:num_elem)) map(from: output_gpu(1:num_elem))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call softmax_gpu(num_slice, slice_size, input, output_gpu)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  avg_ms = ((end_time - start_time) * 1.0d3) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time: ', avg_ms, ' (ms)'

  call softmax_cpu(num_slice, slice_size, input, output_cpu)
  ok = all(abs(output_cpu - output_gpu) <= 1.0e-3_real32)
  if (ok) then
    print '(A)', 'PASS'
  else
    call print_first_mismatch(output_cpu, output_gpu)
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(input, output_gpu, output_cpu)

contains

  subroutine softmax_gpu(num_slice, slice_size, input, output)
    integer(int32), intent(in) :: num_slice, slice_size
    real(real32), intent(in) :: input(:)
    real(real32), intent(inout) :: output(:)
    integer(int32) :: i, j, base
    real(real32) :: max_value, sum_value

    !$omp target teams distribute parallel do simd thread_limit(block_size) private(j, base, max_value, sum_value)
    do i = 1, num_slice
      base = (i - 1_int32) * slice_size
      max_value = input(base + 1_int32)
      do j = 2, slice_size
        if (input(base + j) > max_value) max_value = input(base + j)
      end do
      sum_value = 0.0_real32
      do j = 1, slice_size
        sum_value = sum_value + exp(input(base + j) - max_value)
      end do
      do j = 1, slice_size
        output(base + j) = exp(input(base + j) - max_value) / sum_value
      end do
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine softmax_gpu

  subroutine softmax_cpu(num_slice, slice_size, input, output)
    integer(int32), intent(in) :: num_slice, slice_size
    real(real32), intent(in) :: input(:)
    real(real32), intent(out) :: output(:)
    integer(int32) :: i, j, base
    real(real32) :: max_value, sum_value, e_value

    do i = 1, num_slice
      base = (i - 1_int32) * slice_size
      max_value = input(base + 1_int32)
      do j = 2, slice_size
        if (input(base + j) > max_value) max_value = input(base + j)
      end do
      sum_value = 0.0_real32
      do j = 1, slice_size
        e_value = exp(input(base + j) - max_value)
        sum_value = sum_value + e_value
        output(base + j) = e_value
      end do
      do j = 1, slice_size
        output(base + j) = output(base + j) / sum_value
      end do
    end do
  end subroutine softmax_cpu

  subroutine print_first_mismatch(expected, actual)
    real(real32), intent(in) :: expected(:), actual(:)
    integer(int32) :: idx

    do idx = 1, size(expected)
      if (abs(expected(idx) - actual(idx)) > 1.0e-3_real32) then
        print '("@index ",I0," host: ",F0.6," device: ",F0.6)', idx - 1, expected(idx), actual(idx)
        return
      end if
    end do
  end subroutine print_first_mismatch

end program main
