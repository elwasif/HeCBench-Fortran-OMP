! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer(c_int), parameter :: c_rand_max = 2147483647_c_int

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

  character(len=256) :: arg0
  integer :: batch_size, input_channels, input_height, input_width
  integer :: output_height, output_width, repeat
  integer :: input_numel, output_numel, nthreads, i
  integer, parameter :: ksize_height = 11, ksize_width = 11
  integer, parameter :: stride_height = 4, stride_width = 4
  integer, parameter :: padding_height = 1, padding_width = 1
  real(real32), allocatable :: input(:), output(:), output_grad(:)
  real(real32), allocatable :: input_grad(:), input_grad_ref(:)
  real(real64) :: start_time, end_time, avg_time
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 7) then
    print '(3A)', 'Usage: ', trim(arg0), ' <batch> <input channels> <input height> '
    print '(A)', '<input width> <output height> <output width> <repeat>'
    stop 1
  end if

  batch_size = read_arg(1)
  input_channels = read_arg(2)
  input_height = read_arg(3)
  input_width = read_arg(4)
  output_height = read_arg(5)
  output_width = read_arg(6)
  repeat = read_arg(7)

  if (batch_size <= 0 .or. input_channels <= 0 .or. input_height <= 0 .or. &
      input_width <= 0 .or. output_height <= 0 .or. output_width <= 0 .or. repeat <= 0) then
    stop 1
  end if

  input_numel = batch_size * input_channels * input_height * input_width
  output_numel = batch_size * input_channels * output_height * output_width
  nthreads = input_numel

  allocate(input(input_numel), output(output_numel), output_grad(output_numel))
  allocate(input_grad(input_numel), input_grad_ref(input_numel))

  call c_srand(123_c_int)
  do i = 1, input_numel
    input(i) = real(c_rand(), real32) / real(c_rand_max, real32)
  end do

  do i = 1, output_numel
    output(i) = real(c_rand(), real32) / real(c_rand_max, real32)
    output_grad(i) = real(input_width * input_height, real32)
  end do

  input_grad = 0.0_real32
  input_grad_ref = 0.0_real32

  !$omp target data map(to: input(1:input_numel), output(1:output_numel), output_grad(1:output_numel)) &
  !$omp& map(from: input_grad(1:input_numel))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call pool2d_grad_device(nthreads, input, output, output_grad, input_channels, input_height, input_width, &
                            output_height, output_width, ksize_height, ksize_width, stride_height, stride_width, &
                            padding_height, padding_width, input_grad)
  end do
  end_time = omp_get_wtime()
  avg_time = (end_time - start_time) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time: ', avg_time, ' (s)'
  !$omp end target data

  call pool2d_grad_host(nthreads, input, output, output_grad, input_channels, input_height, input_width, &
                        output_height, output_width, ksize_height, ksize_width, stride_height, stride_width, &
                        padding_height, padding_width, input_grad_ref)

  ok = .true.
  do i = 1, input_numel
    if (abs(input_grad(i) - input_grad_ref(i)) > 1.0e-3_real32) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(input, output, output_grad, input_grad, input_grad_ref)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine pool2d_grad_device(nthreads, input_data, output_data, output_grad, channels, input_height, input_width, &
                                output_height, output_width, ksize_height, ksize_width, stride_height, stride_width, &
                                padding_height, padding_width, input_grad)
    integer, intent(in) :: nthreads, channels, input_height, input_width, output_height, output_width
    integer, intent(in) :: ksize_height, ksize_width, stride_height, stride_width, padding_height, padding_width
    real(real32), intent(in) :: input_data(:), output_data(:), output_grad(:)
    real(real32), intent(out) :: input_grad(:)
    integer :: index, w_offset, h_offset, offset_c, batch_idx, tmp
    integer :: phstart, phend, pwstart, pwend, ph, pw, hstart, hend, wstart, wend
    integer :: pool_size, output_stride, output_sub_idx
    real(real32) :: gradient, scale

    !$omp target teams distribute parallel do thread_limit(256) private(w_offset, h_offset, offset_c, batch_idx, tmp) &
    !$omp& private(phstart, phend, pwstart, pwend, ph, pw, hstart, hend, wstart, wend) &
    !$omp& private(pool_size, output_stride, output_sub_idx, gradient, scale)
    do index = 0, nthreads - 1
      w_offset = mod(index, input_width) + padding_width
      tmp = index / input_width
      h_offset = mod(tmp, input_height) + padding_height
      tmp = tmp / input_height
      offset_c = mod(tmp, channels)
      batch_idx = tmp / channels

      if (h_offset < ksize_height) then
        phstart = 0
      else
        phstart = (h_offset - ksize_height) / stride_height + 1
      end if
      if (w_offset < ksize_width) then
        pwstart = 0
      else
        pwstart = (w_offset - ksize_width) / stride_width + 1
      end if
      phend = min(h_offset / stride_height + 1, output_height)
      pwend = min(w_offset / stride_width + 1, output_width)

      gradient = 0.0_real32
      output_stride = (batch_idx * channels + offset_c) * output_height * output_width

      do ph = phstart, phend - 1
        do pw = pwstart, pwend - 1
          hstart = ph * stride_height - padding_height
          wstart = pw * stride_width - padding_width
          hend = min(hstart + ksize_height, input_height)
          wend = min(wstart + ksize_width, input_width)
          hstart = max(hstart, 0)
          wstart = max(wstart, 0)
          pool_size = (hend - hstart) * (wend - wstart)
          output_sub_idx = ph * output_width + pw
          scale = 1.0_real32 / real(pool_size, real32)
          gradient = gradient + scale * output_grad(output_stride + output_sub_idx + 1)
        end do
      end do
      input_grad(index + 1) = gradient
    end do
    !$omp end target teams distribute parallel do
  end subroutine pool2d_grad_device

  subroutine pool2d_grad_host(nthreads, input_data, output_data, output_grad, channels, input_height, input_width, &
                              output_height, output_width, ksize_height, ksize_width, stride_height, stride_width, &
                              padding_height, padding_width, input_grad)
    integer, intent(in) :: nthreads, channels, input_height, input_width, output_height, output_width
    integer, intent(in) :: ksize_height, ksize_width, stride_height, stride_width, padding_height, padding_width
    real(real32), intent(in) :: input_data(:), output_data(:), output_grad(:)
    real(real32), intent(out) :: input_grad(:)
    integer :: index, w_offset, h_offset, offset_c, batch_idx, tmp
    integer :: phstart, phend, pwstart, pwend, ph, pw, hstart, hend, wstart, wend
    integer :: pool_size, output_stride, output_sub_idx
    real(real32) :: gradient, scale

    do index = 0, nthreads - 1
      w_offset = mod(index, input_width) + padding_width
      tmp = index / input_width
      h_offset = mod(tmp, input_height) + padding_height
      tmp = tmp / input_height
      offset_c = mod(tmp, channels)
      batch_idx = tmp / channels

      if (h_offset < ksize_height) then
        phstart = 0
      else
        phstart = (h_offset - ksize_height) / stride_height + 1
      end if
      if (w_offset < ksize_width) then
        pwstart = 0
      else
        pwstart = (w_offset - ksize_width) / stride_width + 1
      end if
      phend = min(h_offset / stride_height + 1, output_height)
      pwend = min(w_offset / stride_width + 1, output_width)

      gradient = 0.0_real32
      output_stride = (batch_idx * channels + offset_c) * output_height * output_width

      do ph = phstart, phend - 1
        do pw = pwstart, pwend - 1
          hstart = ph * stride_height - padding_height
          wstart = pw * stride_width - padding_width
          hend = min(hstart + ksize_height, input_height)
          wend = min(wstart + ksize_width, input_width)
          hstart = max(hstart, 0)
          wstart = max(wstart, 0)
          pool_size = (hend - hstart) * (wend - wstart)
          output_sub_idx = ph * output_width + pw
          scale = 1.0_real32 / real(pool_size, real32)
          gradient = gradient + scale * output_grad(output_stride + output_sub_idx + 1)
        end do
      end do
      input_grad(index + 1) = gradient
    end do
  end subroutine pool2d_grad_host

end program main
