! SPDX-License-Identifier: CC0-1.0
program resize_benchmark
  use iso_fortran_env, only: int8, int16, int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: channels_per_iter = 8

  integer :: argc
  character(len=64) :: arg
  integer :: in_width, in_height, out_width, out_height, num_channels, repeat

  argc = command_argument_count()
  if (argc /= 6) then
    call get_command_argument(0, arg)
    print '(A,A,A)', 'Usage: ', trim(arg), ' <input image width> <input image height>'
    print '(A)', '          <output image width> <output image height>'
    print '(A)', '          <image channels> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) in_width
  call get_command_argument(2, arg)
  read(arg, *) in_height
  call get_command_argument(3, arg)
  read(arg, *) out_width
  call get_command_argument(4, arg)
  read(arg, *) out_height
  call get_command_argument(5, arg)
  read(arg, *) num_channels
  call get_command_argument(6, arg)
  read(arg, *) repeat

  print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', 'Resize ', num_channels, ' images from (', &
      in_width, ' x ', in_height, ') to (', out_width, ' x ', out_height, ')'

  print *
  print '(A)', 'The size of each pixel is 1 byte'
  call resize_image_i1(in_width, in_height, out_width, out_height, num_channels, repeat, .false.)
  print *
  print '(A)', 'Bilinear resizing'
  call resize_image_i1(in_width, in_height, out_width, out_height, num_channels, repeat, .true.)

  print *
  print '(A)', 'The size of each pixel is 2 bytes'
  call resize_image_i2(in_width, in_height, out_width, out_height, num_channels, repeat, .false.)
  print *
  print '(A)', 'Bilinear resizing'
  call resize_image_i2(in_width, in_height, out_width, out_height, num_channels, repeat, .true.)

  print *
  print '(A)', 'The size of each pixel is 4 bytes'
  call resize_image_i4(in_width, in_height, out_width, out_height, num_channels, repeat, .false.)
  print *
  print '(A)', 'Bilinear resizing'
  call resize_image_i4(in_width, in_height, out_width, out_height, num_channels, repeat, .true.)

contains

  subroutine resize_image_i1(in_width, in_height, out_width, out_height, num_channels, repeat, bilinear)
    integer, intent(in) :: in_width, in_height, out_width, out_height, num_channels, repeat
    logical, intent(in) :: bilinear
    integer(int8), allocatable :: input(:), output(:), reference(:)
    integer(int64) :: in_image_size, out_image_size, in_size, out_size, i
    real(real32) :: fx, fy
    real(real64) :: start_time, end_time, elapsed_ns, perf

    call resize_sizes(in_width, in_height, out_width, out_height, num_channels, &
                      in_image_size, out_image_size, in_size, out_size, fx, fy)
    allocate(input(in_size), output(out_size), reference(out_size))
    do i = 1_int64, in_size
      input(i) = int(mod(i, 13_int64), int8)
    end do

    !$omp target data map(to: input(1:in_size)) map(from: output(1:out_size))
    start_time = omp_get_wtime()
    do i = 1, repeat
      if (bilinear) then
        call resize_bilinear_i1(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      else
        call resize_nearest_i1(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      end if
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    if (bilinear) then
      call host_bilinear_i1(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    else
      call host_nearest_i1(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    end if
    if (any(output /= reference)) print '(A)', 'Resize validation FAILED'

    elapsed_ns = (end_time - start_time) * 1.0e9_real64
    perf = real((in_size + out_size) * 1_int64, real64) * real(repeat, real64) / elapsed_ns
    print '(A,A,A,A,A)', 'Average kernel execution time: ', &
        trim(format_real64(elapsed_ns * 1.0e-3_real64 / real(repeat, real64))), &
        ' (us)    Perf: ', trim(format_real64(perf)), ' (GB/s)'
    deallocate(input, output, reference)
  end subroutine resize_image_i1

  subroutine resize_image_i2(in_width, in_height, out_width, out_height, num_channels, repeat, bilinear)
    integer, intent(in) :: in_width, in_height, out_width, out_height, num_channels, repeat
    logical, intent(in) :: bilinear
    integer(int16), allocatable :: input(:), output(:), reference(:)
    integer(int64) :: in_image_size, out_image_size, in_size, out_size, i
    real(real32) :: fx, fy
    real(real64) :: start_time, end_time, elapsed_ns, perf

    call resize_sizes(in_width, in_height, out_width, out_height, num_channels, &
                      in_image_size, out_image_size, in_size, out_size, fx, fy)
    allocate(input(in_size), output(out_size), reference(out_size))
    do i = 1_int64, in_size
      input(i) = int(mod(i, 13_int64), int16)
    end do

    !$omp target data map(to: input(1:in_size)) map(from: output(1:out_size))
    start_time = omp_get_wtime()
    do i = 1, repeat
      if (bilinear) then
        call resize_bilinear_i2(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      else
        call resize_nearest_i2(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      end if
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    if (bilinear) then
      call host_bilinear_i2(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    else
      call host_nearest_i2(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    end if
    if (any(output /= reference)) print '(A)', 'Resize validation FAILED'

    elapsed_ns = (end_time - start_time) * 1.0e9_real64
    perf = real((in_size + out_size) * 2_int64, real64) * real(repeat, real64) / elapsed_ns
    print '(A,A,A,A,A)', 'Average kernel execution time: ', &
        trim(format_real64(elapsed_ns * 1.0e-3_real64 / real(repeat, real64))), &
        ' (us)    Perf: ', trim(format_real64(perf)), ' (GB/s)'
    deallocate(input, output, reference)
  end subroutine resize_image_i2

  subroutine resize_image_i4(in_width, in_height, out_width, out_height, num_channels, repeat, bilinear)
    integer, intent(in) :: in_width, in_height, out_width, out_height, num_channels, repeat
    logical, intent(in) :: bilinear
    integer(int32), allocatable :: input(:), output(:), reference(:)
    integer(int64) :: in_image_size, out_image_size, in_size, out_size, i
    real(real32) :: fx, fy
    real(real64) :: start_time, end_time, elapsed_ns, perf

    call resize_sizes(in_width, in_height, out_width, out_height, num_channels, &
                      in_image_size, out_image_size, in_size, out_size, fx, fy)
    allocate(input(in_size), output(out_size), reference(out_size))
    do i = 1_int64, in_size
      input(i) = int(mod(i, 13_int64), int32)
    end do

    !$omp target data map(to: input(1:in_size)) map(from: output(1:out_size))
    start_time = omp_get_wtime()
    do i = 1, repeat
      if (bilinear) then
        call resize_bilinear_i4(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      else
        call resize_nearest_i4(output, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
      end if
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    if (bilinear) then
      call host_bilinear_i4(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    else
      call host_nearest_i4(reference, out_size, out_height, out_width, input, in_height, in_width, fx, fy)
    end if
    if (any(output /= reference)) print '(A)', 'Resize validation FAILED'

    elapsed_ns = (end_time - start_time) * 1.0e9_real64
    perf = real((in_size + out_size) * 4_int64, real64) * real(repeat, real64) / elapsed_ns
    print '(A,A,A,A,A)', 'Average kernel execution time: ', &
        trim(format_real64(elapsed_ns * 1.0e-3_real64 / real(repeat, real64))), &
        ' (us)    Perf: ', trim(format_real64(perf)), ' (GB/s)'
    deallocate(input, output, reference)
  end subroutine resize_image_i4

  function format_real64(value) result(text)
    real(real64), intent(in) :: value
    character(len=32) :: text

    write(text, '(F32.6)') value
    text = adjustl(text)
  end function format_real64

  subroutine resize_sizes(in_width, in_height, out_width, out_height, num_channels, &
                          in_image_size, out_image_size, in_size, out_size, fx, fy)
    integer, intent(in) :: in_width, in_height, out_width, out_height, num_channels
    integer(int64), intent(out) :: in_image_size, out_image_size, in_size, out_size
    real(real32), intent(out) :: fx, fy

    in_image_size = int(in_height, int64) * int(in_width, int64)
    out_image_size = int(out_height, int64) * int(out_width, int64)
    in_size = int(num_channels, int64) * in_image_size
    out_size = int(num_channels, int64) * out_image_size
    fx = real(in_width / out_width, real32)
    fy = real(in_height / out_height, real32)
  end subroutine resize_sizes

  subroutine resize_nearest_i1(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int8), intent(inout) :: output(:)
    integer(int8), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c
    real(real32) :: in_yf, in_xf

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,y,x,in_y,in_x,in_idx,out_idx,c,in_yf,in_xf)
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_nearest_i1

  subroutine resize_bilinear_i1(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int8), intent(inout) :: output(:)
    integer(int8), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int8) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,c_end,y,x,in_x0,in_x1,in_y0,in_y1,in_y2) &
    !$omp& private(in_offset_r0,in_offset_r1,out_idx,c,in_x,in_y,v_00,v_01,v_10,v_11)
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int32) + &
            int(int(in_y - real(in_y0, real32), int8), int32) * int(int(v_10 - v_00, int8), int32) + &
            int(int(in_x - real(in_x0, real32), int8), int32) * int(int(v_01 - v_00, int8), int32) + &
            int(int(in_y - real(in_y0, real32), int8), int32) * &
            int(int(in_x - real(in_x0, real32), int8), int32) * &
            int(int(v_11 - v_01 - v_10 + v_00, int8), int32), int8)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_bilinear_i1

  subroutine resize_nearest_i2(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int16), intent(inout) :: output(:)
    integer(int16), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c
    real(real32) :: in_yf, in_xf

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,y,x,in_y,in_x,in_idx,out_idx,c,in_yf,in_xf)
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_nearest_i2

  subroutine resize_bilinear_i2(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int16), intent(inout) :: output(:)
    integer(int16), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int16) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,c_end,y,x,in_x0,in_x1,in_y0,in_y1,in_y2) &
    !$omp& private(in_offset_r0,in_offset_r1,out_idx,c,in_x,in_y,v_00,v_01,v_10,v_11)
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int32) + &
            int(int(in_y - real(in_y0, real32), int16), int32) * int(int(v_10 - v_00, int16), int32) + &
            int(int(in_x - real(in_x0, real32), int16), int32) * int(int(v_01 - v_00, int16), int32) + &
            int(int(in_y - real(in_y0, real32), int16), int32) * &
            int(int(in_x - real(in_x0, real32), int16), int32) * &
            int(int(v_11 - v_01 - v_10 + v_00, int16), int32), int16)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_bilinear_i2

  subroutine resize_nearest_i4(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int32), intent(inout) :: output(:)
    integer(int32), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c
    real(real32) :: in_yf, in_xf

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,y,x,in_y,in_x,in_idx,out_idx,c,in_yf,in_xf)
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_nearest_i4

  subroutine resize_bilinear_i4(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int32), intent(inout) :: output(:)
    integer(int32), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int32) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    !$omp target teams distribute parallel do num_teams(29184) thread_limit(256) &
    !$omp& private(iter,in_image_size,out_image_size,c_start,c_end,y,x,in_x0,in_x1,in_y0,in_y1,in_y2) &
    !$omp& private(in_offset_r0,in_offset_r1,out_idx,c,in_x,in_y,v_00,v_01,v_10,v_11)
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int64) + &
            int(int(in_y - real(in_y0, real32), int32), int64) * int(int(v_10 - v_00, int32), int64) + &
            int(int(in_x - real(in_x0, real32), int32), int64) * int(int(v_01 - v_00, int32), int64) + &
            int(int(in_y - real(in_y0, real32), int32), int64) * &
            int(int(in_x - real(in_x0, real32), int32), int64) * &
            int(int(v_11 - v_01 - v_10 + v_00, int32), int64), int32)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine resize_bilinear_i4

  subroutine nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                             in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
    integer(int64), intent(in) :: iter
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer, intent(out) :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx
    real(real32) :: in_yf, in_xf

    in_image_size = in_height * in_width
    out_image_size = out_height * out_width
    c_start = int(iter / int(out_image_size, int64)) * channels_per_iter
    y = int(mod(iter, int(out_image_size, int64))) / out_width
    x = mod(int(iter), out_width)
    in_yf = (real(y, real32) + 0.5_real32) * o2i_fy
    in_y = nint(in_yf)
    in_xf = (real(x, real32) + 0.5_real32) * o2i_fx
    in_x = nint(in_xf)
    in_x = min(in_x, in_width - 1)
    in_y = min(in_y, in_height - 1)
    in_idx = c_start * in_image_size + in_y * in_width + in_x
    out_idx = c_start * out_image_size + y * out_width + x
  end subroutine nearest_indices

  subroutine bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                              in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                              in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
    integer(int64), intent(in) :: iter
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer, intent(out) :: in_image_size, out_image_size, c_start, c_end, y, x
    integer, intent(out) :: in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx
    real(real32), intent(out) :: in_x, in_y

    in_image_size = in_height * in_width
    out_image_size = out_height * out_width
    c_start = int(iter / int(out_image_size, int64)) * channels_per_iter
    c_end = c_start + channels_per_iter
    y = int(mod(iter, int(out_image_size, int64))) / out_width
    x = mod(int(iter), out_width)
    in_x = max((real(x, real32) + 0.5_real32) * o2i_fx - 0.5_real32, 0.0_real32)
    in_y = max((real(y, real32) + 0.5_real32) * o2i_fy - 0.5_real32, 0.0_real32)
    in_x0 = int(in_x)
    in_x1 = min(in_x0 + 1, in_width - 1)
    in_y0 = int(in_y)
    in_y1 = min(in_y0, in_height - 1)
    in_y2 = min(in_y0 + 1, in_height - 1)
    in_offset_r0 = c_start * in_image_size + in_y1 * in_width
    in_offset_r1 = c_start * in_image_size + in_y2 * in_width
    out_idx = c_start * out_image_size + y * out_width + x
  end subroutine bilinear_indices

  subroutine host_nearest_i1(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int8), intent(out) :: output(:)
    integer(int8), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_nearest_i1

  subroutine host_bilinear_i1(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int8), intent(out) :: output(:)
    integer(int8), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int8) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int32) + &
            int(int(in_y - real(in_y0, real32), int8), int32) * int(int(v_10 - v_00, int8), int32) + &
            int(int(in_x - real(in_x0, real32), int8), int32) * int(int(v_01 - v_00, int8), int32) + &
            int(int(in_y - real(in_y0, real32), int8), int32) * &
            int(int(in_x - real(in_x0, real32), int8), int32) * &
            int(int(v_11 - v_01 - v_10 + v_00, int8), int32), int8)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_bilinear_i1

  subroutine host_nearest_i2(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int16), intent(out) :: output(:)
    integer(int16), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_nearest_i2

  subroutine host_bilinear_i2(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int16), intent(out) :: output(:)
    integer(int16), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int16) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int32) + &
            int(int(in_y - real(in_y0, real32), int16), int32) * int(int(v_10 - v_00, int16), int32) + &
            int(int(in_x - real(in_x0, real32), int16), int32) * int(int(v_01 - v_00, int16), int32) + &
            int(int(in_y - real(in_y0, real32), int16), int32) * &
            int(int(in_x - real(in_x0, real32), int16), int32) * &
            int(int(v_11 - v_01 - v_10 + v_00, int16), int32), int16)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_bilinear_i2

  subroutine host_nearest_i4(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int32), intent(out) :: output(:)
    integer(int32), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx, c

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call nearest_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                           in_image_size, out_image_size, c_start, y, x, in_y, in_x, in_idx, out_idx)
      do c = 1, channels_per_iter
        output(out_idx + 1) = input(in_idx + 1)
        in_idx = in_idx + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_nearest_i4

  subroutine host_bilinear_i4(output, output_size, out_height, out_width, input, in_height, in_width, o2i_fy, o2i_fx)
    integer(int32), intent(out) :: output(:)
    integer(int32), intent(in) :: input(:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: out_height, out_width, in_height, in_width
    real(real32), intent(in) :: o2i_fy, o2i_fx
    integer(int64) :: iter, iters_required
    integer :: in_image_size, out_image_size, c_start, c_end, y, x, in_x0, in_x1, in_y0, in_y1, in_y2
    integer :: in_offset_r0, in_offset_r1, out_idx, c
    integer(int32) :: v_00, v_01, v_10, v_11
    real(real32) :: in_x, in_y

    iters_required = output_size / channels_per_iter
    do iter = 0_int64, iters_required - 1_int64
      call bilinear_indices(iter, out_height, out_width, in_height, in_width, o2i_fy, o2i_fx, &
                            in_image_size, out_image_size, c_start, c_end, y, x, in_x, in_y, &
                            in_x0, in_x1, in_y0, in_y1, in_y2, in_offset_r0, in_offset_r1, out_idx)
      do c = c_start, c_end - 1
        v_00 = input(in_offset_r0 + in_x0 + 1)
        v_01 = input(in_offset_r0 + in_x1 + 1)
        v_10 = input(in_offset_r1 + in_x0 + 1)
        v_11 = input(in_offset_r1 + in_x1 + 1)
        output(out_idx + 1) = int( &
            int(v_00, int64) + &
            int(int(in_y - real(in_y0, real32), int32), int64) * int(int(v_10 - v_00, int32), int64) + &
            int(int(in_x - real(in_x0, real32), int32), int64) * int(int(v_01 - v_00, int32), int64) + &
            int(int(in_y - real(in_y0, real32), int32), int64) * &
            int(int(in_x - real(in_x0, real32), int32), int64) * &
            int(int(v_11 - v_01 - v_10 + v_00, int32), int64), int32)
        in_offset_r0 = in_offset_r0 + in_image_size
        in_offset_r1 = in_offset_r1 + in_image_size
        out_idx = out_idx + out_image_size
      end do
    end do
  end subroutine host_bilinear_i4

end program resize_benchmark
