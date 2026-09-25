! SPDX-License-Identifier: CC0-1.0
program flip_main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer(int64) :: num_dims, dim_size, num_flip_dims
  integer(int32) :: repeat
  character(len=64) :: arg

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <number of dimensions> <size of each dimension> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) num_dims
  call get_command_argument(2, arg)
  read(arg, *) dim_size
  call get_command_argument(3, arg)
  read(arg, *) repeat

  num_flip_dims = num_dims

  write(*,'("=========== Data type is FP32 ==========")')
  call flip_real32(num_dims, num_flip_dims, dim_size, repeat)

  write(*,'("=========== Data type is FP64 ==========")')
  call flip_real64(num_dims, num_flip_dims, dim_size, repeat)

contains

  subroutine print_property(name, values)
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: values(0:)
    integer :: i

    write(*,'(A,": ( ")', advance='no') trim(name)
    do i = 0, size(values) - 1
      write(*,'(I0,1X)', advance='no') values(i)
    end do
    write(*,'(")")')
  end subroutine print_property

  subroutine setup_properties(num_dims, num_flip_dims, dim_size, shape, flip_dims, stride)
    integer(int64), intent(in) :: num_dims, num_flip_dims, dim_size
    integer(int64), intent(out) :: shape(0:), flip_dims(0:), stride(0:)
    integer(int64) :: i

    do i = 0, num_dims - 1
      shape(i) = dim_size
    end do
    do i = 0, num_flip_dims - 1
      flip_dims(i) = i
    end do

    stride = 1_int64
    if (num_dims >= 3) then
      stride(0) = shape(1) * shape(2)
      stride(1) = shape(2)
      stride(2) = 1_int64
    else if (num_dims == 2) then
      stride(0) = shape(1)
      stride(1) = 1_int64
    else if (num_dims == 1) then
      stride(0) = 1_int64
    end if
  end subroutine setup_properties

  subroutine flip_real32(num_dims, num_flip_dims, dim_size, repeat)
    integer(int64), intent(in) :: num_dims, num_flip_dims, dim_size
    integer(int32), intent(in) :: repeat
    integer(int64), allocatable :: shape(:), flip_dims(:), stride(:)
    integer(int64) :: n, i
    real(real32), allocatable :: input(:), output(:), output_ref(:)
    real(real64) :: start_time, elapsed_ms
    logical :: error

    allocate(shape(0:num_dims-1), flip_dims(0:num_flip_dims-1), stride(0:num_dims-1))
    call setup_properties(num_dims, num_flip_dims, dim_size, shape, flip_dims, stride)
    n = product(shape)

    call print_property("shape", shape)
    call print_property("flip_dims", flip_dims)
    call print_property("stride", stride)

    allocate(input(0:n-1), output(0:n-1), output_ref(0:n-1))
    do i = 0, n - 1
      input(i) = real(i, real32)
    end do

    !$omp target data map(to: input, shape, flip_dims, stride) map(alloc: output)
    call flip_kernel_real32(input, output, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    !$omp target update from(output)
    call flip_cpu_real32(input, output_ref, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    error = any(output /= output_ref)
    write(*,'(A)') merge("FAIL", "PASS", error)

    start_time = omp_get_wtime()
    do i = 1, repeat
      call flip_kernel_real32(input, output, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    end do
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
    write(*,'("Average execution time of the flip kernel: ",F0.6," (ms)")') elapsed_ms
    !$omp end target data
  end subroutine flip_real32

  subroutine flip_real64(num_dims, num_flip_dims, dim_size, repeat)
    integer(int64), intent(in) :: num_dims, num_flip_dims, dim_size
    integer(int32), intent(in) :: repeat
    integer(int64), allocatable :: shape(:), flip_dims(:), stride(:)
    integer(int64) :: n, i
    real(real64), allocatable :: input(:), output(:), output_ref(:)
    real(real64) :: start_time, elapsed_ms
    logical :: error

    allocate(shape(0:num_dims-1), flip_dims(0:num_flip_dims-1), stride(0:num_dims-1))
    call setup_properties(num_dims, num_flip_dims, dim_size, shape, flip_dims, stride)
    n = product(shape)

    call print_property("shape", shape)
    call print_property("flip_dims", flip_dims)
    call print_property("stride", stride)

    allocate(input(0:n-1), output(0:n-1), output_ref(0:n-1))
    do i = 0, n - 1
      input(i) = real(i, real64)
    end do

    !$omp target data map(to: input, shape, flip_dims, stride) map(alloc: output)
    call flip_kernel_real64(input, output, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    !$omp target update from(output)
    call flip_cpu_real64(input, output_ref, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    error = any(output /= output_ref)
    write(*,'(A)') merge("FAIL", "PASS", error)

    start_time = omp_get_wtime()
    do i = 1, repeat
      call flip_kernel_real64(input, output, n, flip_dims, num_flip_dims, stride, stride, shape, num_dims)
    end do
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
    write(*,'("Average execution time of the flip kernel: ",F0.6," (ms)")') elapsed_ms
    !$omp end target data
  end subroutine flip_real64

  subroutine flip_kernel_real32(input, output, n, flip_dims, flip_dims_size, strides, strides_contiguous, shape, total_dims)
    real(real32), intent(in) :: input(0:)
    real(real32), intent(out) :: output(0:)
    integer(int64), intent(in) :: n, flip_dims(0:), flip_dims_size, strides(0:), strides_contiguous(0:), shape(0:), total_dims
    integer(int64) :: linear_index, cur_indices, rem, dst_offset, i, j, temp

    !$omp target teams distribute parallel do thread_limit(256) private(cur_indices, rem, dst_offset, i, j, temp)
    do linear_index = 0, n - 1
      cur_indices = linear_index
      rem = 0_int64
      dst_offset = 0_int64
      do i = 0, total_dims - 1
        temp = cur_indices
        cur_indices = cur_indices / strides_contiguous(i)
        rem = temp - cur_indices * strides_contiguous(i)
        do j = 0, flip_dims_size - 1
          if (i == flip_dims(j)) then
            cur_indices = shape(i) - 1_int64 - cur_indices
          end if
        end do
        dst_offset = dst_offset + cur_indices * strides(i)
        cur_indices = rem
      end do
      output(linear_index) = input(dst_offset)
    end do
    !$omp end target teams distribute parallel do
  end subroutine flip_kernel_real32

  subroutine flip_kernel_real64(input, output, n, flip_dims, flip_dims_size, strides, strides_contiguous, shape, total_dims)
    real(real64), intent(in) :: input(0:)
    real(real64), intent(out) :: output(0:)
    integer(int64), intent(in) :: n, flip_dims(0:), flip_dims_size, strides(0:), strides_contiguous(0:), shape(0:), total_dims
    integer(int64) :: linear_index, cur_indices, rem, dst_offset, i, j, temp

    !$omp target teams distribute parallel do thread_limit(256) private(cur_indices, rem, dst_offset, i, j, temp)
    do linear_index = 0, n - 1
      cur_indices = linear_index
      rem = 0_int64
      dst_offset = 0_int64
      do i = 0, total_dims - 1
        temp = cur_indices
        cur_indices = cur_indices / strides_contiguous(i)
        rem = temp - cur_indices * strides_contiguous(i)
        do j = 0, flip_dims_size - 1
          if (i == flip_dims(j)) then
            cur_indices = shape(i) - 1_int64 - cur_indices
          end if
        end do
        dst_offset = dst_offset + cur_indices * strides(i)
        cur_indices = rem
      end do
      output(linear_index) = input(dst_offset)
    end do
    !$omp end target teams distribute parallel do
  end subroutine flip_kernel_real64

  subroutine flip_cpu_real32(input, output, n, flip_dims, flip_dims_size, strides, strides_contiguous, shape, total_dims)
    real(real32), intent(in) :: input(0:)
    real(real32), intent(out) :: output(0:)
    integer(int64), intent(in) :: n, flip_dims(0:), flip_dims_size, strides(0:), strides_contiguous(0:), shape(0:), total_dims
    integer(int64) :: linear_index, cur_indices, rem, dst_offset, i, j, temp

    do linear_index = 0, n - 1
      cur_indices = linear_index
      rem = 0_int64
      dst_offset = 0_int64
      do i = 0, total_dims - 1
        temp = cur_indices
        cur_indices = cur_indices / strides_contiguous(i)
        rem = temp - cur_indices * strides_contiguous(i)
        do j = 0, flip_dims_size - 1
          if (i == flip_dims(j)) then
            cur_indices = shape(i) - 1_int64 - cur_indices
          end if
        end do
        dst_offset = dst_offset + cur_indices * strides(i)
        cur_indices = rem
      end do
      output(linear_index) = input(dst_offset)
    end do
  end subroutine flip_cpu_real32

  subroutine flip_cpu_real64(input, output, n, flip_dims, flip_dims_size, strides, strides_contiguous, shape, total_dims)
    real(real64), intent(in) :: input(0:)
    real(real64), intent(out) :: output(0:)
    integer(int64), intent(in) :: n, flip_dims(0:), flip_dims_size, strides(0:), strides_contiguous(0:), shape(0:), total_dims
    integer(int64) :: linear_index, cur_indices, rem, dst_offset, i, j, temp

    do linear_index = 0, n - 1
      cur_indices = linear_index
      rem = 0_int64
      dst_offset = 0_int64
      do i = 0, total_dims - 1
        temp = cur_indices
        cur_indices = cur_indices / strides_contiguous(i)
        rem = temp - cur_indices * strides_contiguous(i)
        do j = 0, flip_dims_size - 1
          if (i == flip_dims(j)) then
            cur_indices = shape(i) - 1_int64 - cur_indices
          end if
        end do
        dst_offset = dst_offset + cur_indices * strides(i)
        cur_indices = rem
      end do
      output(linear_index) = input(dst_offset)
    end do
  end subroutine flip_cpu_real64

end program flip_main
