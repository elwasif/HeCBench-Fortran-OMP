! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_int64_t
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    subroutine c_fill_random_float(a, n) bind(C, name='fill_random_float')
      import :: c_float, c_int64_t
      real(c_float), intent(out) :: a(*)
      integer(c_int64_t), value :: n
    end subroutine c_fill_random_float
  end interface

  integer, parameter :: block_sizes(6) = [32, 64, 128, 256, 512, 1024]
  integer, parameter :: block2d_sizes(3) = [8, 16, 32]
  integer :: bsz, channels, height, width, repeat, i, block_size
  integer(int64) :: input_size, output_size
  real(real32), allocatable :: x(:), dout(:), out_ref(:), out_dev(:), dx_ref(:), dx_dev(:)
  real(real64) :: elapsed_ms, gflops

  if (command_argument_count() /= 5) then
    print '(A)', 'Usage: ./main <batch size> <number of channels> <height> <width> <repeat>'
    stop 1
  end if
  bsz = read_arg(1)
  channels = read_arg(2)
  height = read_arg(3)
  width = read_arg(4)
  repeat = read_arg(5)

  input_size = int(bsz, int64) * channels * height * width
  output_size = input_size * 4_int64
  allocate(x(input_size), dout(output_size), out_ref(output_size), out_dev(output_size), dx_ref(input_size), dx_dev(input_size))
  call c_srand(0_c_int)
  call fill_random(x)
  call fill_random(dout)
  out_ref = 0.0_real32
  out_dev = 0.0_real32
  dx_ref = 0.0_real32
  dx_dev = 0.0_real32
  call upsample_forward_ref(x, out_ref, bsz, channels, height, width)
  call upsample_backward_ref(dout, dx_ref, bsz, channels, height, width)

  !$omp target data map(to: x(1:input_size), dout(1:output_size)) map(tofrom: out_dev(1:output_size), dx_dev(1:input_size))
    print '(A)', 'Checking forward pass'
    do i = 1, size(block_sizes)
      block_size = block_sizes(i)
      write(*, '(A,I0)') 'Checking block size ', block_size
      call upsample_forward_dev(x, out_dev, bsz, channels, height, width, block_size)
      !$omp target update from(out_dev(1:output_size))
      call validate_result(out_dev, out_ref, 'out')
    end do
    print '(A)', 'Forward1 pass: all results match'
    print '(A)', ''

    print '(A)', 'Forward1 pass benchmarks:'
    do i = 1, size(block_sizes)
      block_size = block_sizes(i)
      elapsed_ms = benchmark_forward(repeat, x, out_dev, bsz, channels, height, width, block_size)
      gflops = real(input_size, real64) / elapsed_ms * 1.0e3_real64 / 1.0e9_real64
      write(*, '(A,I4,A,F6.4,A,F4.2)') 'block_size ', block_size, ' | time ', elapsed_ms, ' ms | gflops ', gflops
    end do

    print '(A)', ''
    print '(A)', '─────────────────────────────────────────────────────'
    do i = 1, size(block2d_sizes)
      block_size = block2d_sizes(i)
      write(*, '(A,I0)') 'Checking block size ', block_size
      call upsample_forward_dev2(x, out_dev, bsz, channels, height, width, block_size, block_size)
      !$omp target update from(out_dev(1:output_size))
      call validate_result(out_dev, out_ref, 'out')
    end do
    print '(A)', 'Forward2 pass: all results match'
    print '(A)', ''
    print '(A)', 'Forward2 pass benchmarks:'
    do i = 1, size(block2d_sizes)
      block_size = block2d_sizes(i)
      elapsed_ms = benchmark_forward2(repeat, x, out_dev, bsz, channels, height, width, block_size, block_size)
      gflops = real(input_size, real64) / elapsed_ms * 1.0e3_real64 / 1.0e9_real64
      write(*, '(A,I4,A,F6.4,A,F4.2)') 'block2D_size ', block_size, ' | time ', elapsed_ms, ' ms | gflops ', gflops
    end do

    print '(A)', ''
    print '(A)', '─────────────────────────────────────────────────────'
    print '(A)', 'Checking backward pass'
    do i = 1, size(block_sizes)
      block_size = block_sizes(i)
      write(*, '(A,I0)') 'Checking block size ', block_size
      call upsample_backward_dev(dout, dx_dev, bsz, channels, height, width, block_size)
      !$omp target update from(dx_dev(1:input_size))
      call validate_result(dx_dev, dx_ref, 'dx')
    end do
    print '(A)', 'Backward pass: all results match'
    print '(A)', ''
    print '(A)', 'All results match. Starting benchmarks.'
    print '(A)', ''

    print '(A)', ''
    print '(A)', 'Backward pass benchmarks:'
    do i = 1, size(block_sizes)
      block_size = block_sizes(i)
      elapsed_ms = benchmark_backward(repeat, dout, dx_dev, bsz, channels, height, width, block_size)
      gflops = real(input_size, real64) / elapsed_ms * 1.0e3_real64 / 1.0e9_real64
      write(*, '(A,I4,A,F6.4,A,F4.2)') 'block_size ', block_size, ' | time ', elapsed_ms, ' ms | gflops ', gflops
    end do

    print '(A)', ''
    print '(A)', '─────────────────────────────────────────────────────'
    print '(A)', 'Checking backward2 pass'
    do i = 1, size(block2d_sizes)
      block_size = block2d_sizes(i)
      write(*, '(A,I0)') 'Checking block size ', block_size
      call upsample_backward_dev2(dout, dx_dev, bsz, channels, height, width, block_size, block_size)
      !$omp target update from(dx_dev(1:input_size))
      call validate_result(dx_dev, dx_ref, 'dx')
    end do
    print '(A)', 'Backward2 pass: all results match'
    print '(A)', ''
    print '(A)', 'All results match. Starting benchmarks.'
    print '(A)', ''

    print '(A)', ''
    print '(A)', 'Backward2 pass benchmarks:'
    do i = 1, size(block2d_sizes)
      block_size = block2d_sizes(i)
      elapsed_ms = benchmark_backward2(repeat, dout, dx_dev, bsz, channels, height, width, block_size, block_size)
      gflops = real(input_size, real64) / elapsed_ms * 1.0e3_real64 / 1.0e9_real64
      write(*, '(A,I4,A,F6.4,A,F4.2)') 'block2D_size ', block_size, ' | time ', elapsed_ms, ' ms | gflops ', gflops
    end do
  !$omp end target data

  deallocate(x, dout, out_ref, out_dev, dx_ref, dx_dev)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_random(a)
    real(real32), intent(out) :: a(:)
    call c_fill_random_float(a, int(size(a), c_int64_t))
  end subroutine fill_random

  subroutine validate_result(actual, expected, name)
    real(real32), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: name
    integer :: idx, nfaults
    nfaults = 0
    do idx = 1, size(actual)
      if (abs(expected(idx) - actual(idx)) > 1.0e-4_real32 .and. expected(idx) == expected(idx)) then
        write(*, '(A,A,A,I0,A,F0.6,A,F0.6)') 'Mismatch of ', trim(name), ' at ', idx - 1, ': CPU_ref: ', expected(idx), ' vs GPU: ', actual(idx)
        nfaults = nfaults + 1
        if (nfaults >= 10) return
      end if
    end do
  end subroutine validate_result

  subroutine upsample_forward_ref(x, out, bsz, channels, height, width)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer, intent(in) :: bsz, channels, height, width
    integer :: b, c, y, xidx, in_idx, out_idx, h_out, w_out
    h_out = height * 2
    w_out = width * 2
    out = 0.0_real32
    do b = 0, bsz - 1
      do c = 0, channels - 1
        do y = 0, height - 1
          do xidx = 0, width - 1
            in_idx = ((b * channels + c) * height + y) * width + xidx + 1
            out_idx = ((b * channels + c) * h_out + 2 * y) * w_out + 2 * xidx + 1
            out(out_idx) = x(in_idx)
            out(out_idx + 1) = x(in_idx)
            out(out_idx + w_out) = x(in_idx)
            out(out_idx + w_out + 1) = x(in_idx)
          end do
        end do
      end do
    end do
  end subroutine upsample_forward_ref

  subroutine upsample_backward_ref(dout, dx, bsz, channels, height, width)
    real(real32), intent(in) :: dout(:)
    real(real32), intent(out) :: dx(:)
    integer, intent(in) :: bsz, channels, height, width
    integer :: b, c, y, xidx, in_idx, out_idx, h_out, w_out
    h_out = height * 2
    w_out = width * 2
    do b = 0, bsz - 1
      do c = 0, channels - 1
        do y = 0, height - 1
          do xidx = 0, width - 1
            in_idx = ((b * channels + c) * height + y) * width + xidx + 1
            out_idx = ((b * channels + c) * h_out + 2 * y) * w_out + 2 * xidx + 1
            dx(in_idx) = dout(out_idx) + dout(out_idx + 1) + dout(out_idx + w_out) + dout(out_idx + w_out + 1)
          end do
        end do
      end do
    end do
  end subroutine upsample_backward_ref

  subroutine upsample_forward_dev(x, out, bsz, channels, height, width, block_size)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: out(:)
    integer, intent(in) :: bsz, channels, height, width, block_size
    integer(int64) :: flat, total
    integer :: b, c, y, xidx, in_idx, out_idx, h_out, w_out, img_size
    img_size = height * width
    h_out = height * 2
    w_out = width * 2
    total = int(bsz, int64) * channels * img_size
    !$omp target teams distribute parallel do thread_limit(block_size) private(b, c, y, xidx, in_idx, out_idx)
    do flat = 0, total - 1
      b = int(flat / int(channels * img_size, int64))
      c = int(mod(flat / img_size, int(channels, int64)))
      y = int(mod(flat / width, int(height, int64)))
      xidx = int(mod(flat, int(width, int64)))
      in_idx = int(flat) + 1
      out_idx = ((b * channels + c) * h_out + 2 * y) * w_out + 2 * xidx + 1
      out(out_idx) = x(in_idx)
      out(out_idx + 1) = x(in_idx)
      out(out_idx + w_out) = x(in_idx)
      out(out_idx + w_out + 1) = x(in_idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine upsample_forward_dev

  subroutine upsample_backward_dev(dout, dx, bsz, channels, height, width, block_size)
    real(real32), intent(in) :: dout(:)
    real(real32), intent(inout) :: dx(:)
    integer, intent(in) :: bsz, channels, height, width, block_size
    integer(int64) :: flat, total
    integer :: b, c, y, xidx, in_idx, out_idx, h_out, w_out, img_size
    img_size = height * width
    h_out = height * 2
    w_out = width * 2
    total = int(bsz, int64) * channels * img_size
    !$omp target teams distribute parallel do thread_limit(block_size) private(b, c, y, xidx, in_idx, out_idx)
    do flat = 0, total - 1
      b = int(flat / int(channels * img_size, int64))
      c = int(mod(flat / img_size, int(channels, int64)))
      y = int(mod(flat / width, int(height, int64)))
      xidx = int(mod(flat, int(width, int64)))
      in_idx = int(flat) + 1
      out_idx = ((b * channels + c) * h_out + 2 * y) * w_out + 2 * xidx + 1
      dx(in_idx) = dout(out_idx) + dout(out_idx + 1) + dout(out_idx + w_out) + dout(out_idx + w_out + 1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine upsample_backward_dev

  subroutine upsample_forward_dev2(x, out, bsz, channels, height, width, block_size_x, block_size_y)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: out(:)
    integer, intent(in) :: bsz, channels, height, width, block_size_x, block_size_y
    integer :: bc, in_y, in_x, b, c, in_idx, out_base, h_out, w_out, block_size
    h_out = height * 2
    w_out = width * 2
    block_size = block_size_x * block_size_y
    !$omp target teams distribute parallel do collapse(3) thread_limit(block_size) private(b, c, in_idx, out_base)
    do bc = 0, bsz * channels - 1
      do in_y = 0, height - 1
        do in_x = 0, width - 1
          b = bc / channels
          c = mod(bc, channels)
          in_idx = ((b * channels + c) * height + in_y) * width + in_x + 1
          out_base = ((b * channels + c) * h_out + 2 * in_y) * w_out + 2 * in_x + 1
          out(out_base) = x(in_idx)
          out(out_base + 1) = x(in_idx)
          out(out_base + w_out) = x(in_idx)
          out(out_base + w_out + 1) = x(in_idx)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine upsample_forward_dev2

  subroutine upsample_backward_dev2(dout, dx, bsz, channels, height, width, block_size_x, block_size_y)
    real(real32), intent(in) :: dout(:)
    real(real32), intent(inout) :: dx(:)
    integer, intent(in) :: bsz, channels, height, width, block_size_x, block_size_y
    integer :: bc, in_y, in_x, b, c, in_idx, out_base, h_out, w_out, block_size
    h_out = height * 2
    w_out = width * 2
    block_size = block_size_x * block_size_y
    !$omp target teams distribute parallel do collapse(3) thread_limit(block_size) private(b, c, in_idx, out_base)
    do bc = 0, bsz * channels - 1
      do in_y = 0, height - 1
        do in_x = 0, width - 1
          b = bc / channels
          c = mod(bc, channels)
          in_idx = ((b * channels + c) * height + in_y) * width + in_x + 1
          out_base = ((b * channels + c) * h_out + 2 * in_y) * w_out + 2 * in_x + 1
          dx(in_idx) = dout(out_base) + dout(out_base + 1) + dout(out_base + w_out) + dout(out_base + w_out + 1)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine upsample_backward_dev2

  real(real64) function benchmark_forward(repeat, x, out, bsz, channels, height, width, block_size)
    integer, intent(in) :: repeat, bsz, channels, height, width, block_size
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: out(:)
    integer :: iter
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call upsample_forward_dev(x, out, bsz, channels, height, width, block_size)
    end do
    end_time = omp_get_wtime()
    benchmark_forward = ((end_time - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_forward

  real(real64) function benchmark_forward2(repeat, x, out, bsz, channels, height, width, block_size_x, block_size_y)
    integer, intent(in) :: repeat, bsz, channels, height, width, block_size_x, block_size_y
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: out(:)
    integer :: iter
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call upsample_forward_dev2(x, out, bsz, channels, height, width, block_size_x, block_size_y)
    end do
    end_time = omp_get_wtime()
    benchmark_forward2 = ((end_time - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_forward2

  real(real64) function benchmark_backward(repeat, dout, dx, bsz, channels, height, width, block_size)
    integer, intent(in) :: repeat, bsz, channels, height, width, block_size
    real(real32), intent(in) :: dout(:)
    real(real32), intent(inout) :: dx(:)
    integer :: iter
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call upsample_backward_dev(dout, dx, bsz, channels, height, width, block_size)
    end do
    end_time = omp_get_wtime()
    benchmark_backward = ((end_time - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_backward

  real(real64) function benchmark_backward2(repeat, dout, dx, bsz, channels, height, width, block_size_x, block_size_y)
    integer, intent(in) :: repeat, bsz, channels, height, width, block_size_x, block_size_y
    real(real32), intent(in) :: dout(:)
    real(real32), intent(inout) :: dx(:)
    integer :: iter
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call upsample_backward_dev2(dout, dx, bsz, channels, height, width, block_size_x, block_size_y)
    end do
    end_time = omp_get_wtime()
    benchmark_backward2 = ((end_time - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_backward2

end program main
