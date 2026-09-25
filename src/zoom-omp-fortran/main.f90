! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_long_long
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg
  integer :: input_sizes(4), repeat
  real(real32) :: zf(2)

  interface
    subroutine zoom_fill_input(input_img, img_size) bind(C, name="zoom_fill_input")
      import :: c_float, c_long_long
      real(c_float), intent(out) :: input_img(*)
      integer(c_long_long), value :: img_size
    end subroutine zoom_fill_input
  end interface

  if (command_argument_count() /= 5) then
    call get_command_argument(0, arg)
    write(*,'(A,A,A)') 'Usage: ', trim(arg), ' <batch> <channel> <height> <width> <repeat>'
    stop 1
  end if

  input_sizes(1) = read_int_arg(1)
  input_sizes(2) = read_int_arg(2)
  input_sizes(3) = read_int_arg(3)
  input_sizes(4) = read_int_arg(4)
  repeat = read_int_arg(5)
  if (any(input_sizes <= 0) .or. repeat <= 0) stop 1

  zf = [1.5_real32, 2.5_real32]
  call zoom(repeat, input_sizes, zf)

  zf = [0.6_real32, 0.9_real32]
  call zoom(repeat, input_sizes, zf)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  subroutine zoom(repeat, input_sizes, zoom_factor)
    integer, intent(in) :: repeat, input_sizes(4)
    real(real32), intent(in) :: zoom_factor(2)
    integer :: n, c, h, w, ho, wo, batch_size
    integer :: pad_dims(2, 2), slice_dims(2, 2), diff, half
    integer(int64) :: pitch, img_size
    real(real32), allocatable :: input_img(:), output_img(:), output_img_ref(:)
    integer(int64) :: i
    logical :: is_zoom_out, is_zoom_in, ok
    real(real64) :: start_time, elapsed

    n = input_sizes(1)
    c = input_sizes(2)
    h = input_sizes(3)
    w = input_sizes(4)
    ho = floor_real32(real(h, real32) * zoom_factor(1))
    wo = floor_real32(real(w, real32) * zoom_factor(2))
    is_zoom_out = ho < h .and. wo < w
    is_zoom_in = ho > h .and. wo > w
    if (.not. is_zoom_out .and. .not. is_zoom_in) then
      write(*,'(A)') 'Zoom factors only handle simultaneous expansion(or shrinkage) in both dimensions. Exit'
      stop 1
    end if

    pitch = int(h, int64) * int(w, int64)
    batch_size = c * n
    pad_dims = 0
    slice_dims = 0
    diff = h - ho
    half = abs(diff) / 2
    if (diff > 0) then
      pad_dims(1, 1) = half
      pad_dims(2, 1) = diff - half
    else
      slice_dims(1, 1) = half
      slice_dims(2, 1) = h + half
    end if
    diff = w - wo
    half = abs(diff) / 2
    if (diff > 0) then
      pad_dims(1, 2) = half
      pad_dims(2, 2) = diff - half
    else
      slice_dims(1, 2) = half
      slice_dims(2, 2) = w + half
    end if

    img_size = pitch * int(batch_size, int64)
    allocate(input_img(0:img_size - 1), output_img(0:img_size - 1), output_img_ref(0:img_size - 1))
    call zoom_fill_input(input_img, int(img_size, c_long_long))
    output_img = 0.0_real32
    output_img_ref = 0.0_real32

    !$omp target data map(to: input_img(0:img_size-1)) map(from: output_img(0:img_size-1))
      start_time = omp_get_wtime()
      do i = 1, repeat
        if (is_zoom_in) then
          call zoom_in_kernel(input_img, output_img, h, w, ho, wo, pitch, &
            slice_dims(1, 1), slice_dims(2, 1), slice_dims(1, 2), slice_dims(2, 2), batch_size)
        else
          call zoom_out_kernel(input_img, output_img, h, w, ho, wo, pitch, &
            pad_dims(1, 1), pad_dims(2, 1), pad_dims(1, 2), pad_dims(2, 2), batch_size)
          call zoom_out_edge_pad(output_img, h, w, pitch, pad_dims(1, 1), pad_dims(1, 2), &
            pad_dims(1, 1) + ho, pad_dims(1, 2) + wo, batch_size)
        end if
      end do
      elapsed = omp_get_wtime() - start_time
      if (is_zoom_in) then
        write(*,'(A,F0.6,A)') 'Average execution time of the zoom-in kernel: ', &
          1.0e6_real64 * elapsed / real(repeat, real64), ' (us)'
      else
        write(*,'(A,F0.6,A)') 'Average execution time of the zoom-out kernel: ', &
          1.0e6_real64 * elapsed / real(repeat, real64), ' (us)'
      end if
    !$omp end target data

    if (is_zoom_in) then
      call zoom_in_reference(input_img, output_img_ref, h, w, ho, wo, pitch, &
        slice_dims(1, 1), slice_dims(2, 1), slice_dims(1, 2), slice_dims(2, 2), batch_size)
    else
      call zoom_out_reference(input_img, output_img_ref, h, w, ho, wo, pitch, &
        pad_dims(1, 1), pad_dims(2, 1), pad_dims(1, 2), pad_dims(2, 2), batch_size)
      call zoom_out_edge_pad_reference(output_img_ref, h, w, pitch, pad_dims(1, 1), pad_dims(1, 2), &
        pad_dims(1, 1) + ho, pad_dims(1, 2) + wo, batch_size)
    end if

    ok = .true.
    do i = 0, img_size - 1
      if (abs(output_img(i) - output_img_ref(i)) > 1.0e-4_real32) then
        ok = .false.
        exit
      end if
    end do
    write(*,'(A)') merge('PASS', 'FAIL', ok)
    deallocate(input_img, output_img, output_img_ref)
  end subroutine zoom

  integer function floor_real32(value) result(out)
    real(real32), intent(in) :: value
    out = floor(value)
  end function floor_real32

  integer function ceil_real32(value) result(out)
    real(real32), intent(in) :: value
    out = ceiling(value)
  end function ceil_real32

  subroutine zoom_in_kernel(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, &
      out_h_start, out_h_end, out_w_start, out_w_end, batch_size)
    real(real32), intent(in) :: input_tensor(0:)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: input_h, input_w, output_h, output_w, out_h_start, out_h_end, out_w_start, out_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow, i, j, src_row, src_col, start_h, end_h, start_w, end_w, del_h, del_w
    real(real32) :: ratio_h, ratio_w, sum_value
    integer(int64) :: base, out_idx, in_idx

    ratio_h = real(input_h, real32) / real(output_h, real32)
    ratio_w = real(input_w, real32) / real(output_w, real32)
    !$omp target teams distribute parallel do collapse(3) &
    !$omp& private(base, start_h, end_h, start_w, end_w, del_h, del_w, sum_value, i, j, src_row, src_col, in_idx, out_idx)
    do b = 0, batch_size - 1
      do oh = out_h_start, out_h_end - 1
        do ow = out_w_start, out_w_end - 1
          base = int(b, int64) * pitch
          start_h = floor_real32(real(oh, real32) * ratio_h)
          end_h = ceil_real32(real(oh + 1, real32) * ratio_h)
          start_w = floor_real32(real(ow, real32) * ratio_w)
          end_w = ceil_real32(real(ow + 1, real32) * ratio_w)
          del_h = end_h - start_h
          del_w = end_w - start_w
          sum_value = 0.0_real32
          do i = 0, del_h - 1
            src_row = start_h + i
            if (src_row >= input_h) cycle
            do j = 0, del_w - 1
              src_col = start_w + j
              if (src_col >= input_w) cycle
              in_idx = base + int(src_row * input_w + src_col, int64)
              sum_value = sum_value + input_tensor(in_idx)
            end do
          end do
          out_idx = base + int((oh - out_h_start) * input_w + (ow - out_w_start), int64)
          output_tensor(out_idx) = sum_value / real(del_h * del_w, real32)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine zoom_in_kernel

  subroutine zoom_out_kernel(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, &
      out_h_start, out_h_end, out_w_start, out_w_end, batch_size)
    real(real32), intent(in) :: input_tensor(0:)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: input_h, input_w, output_h, output_w, out_h_start, out_h_end, out_w_start, out_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow, i, j, src_row, src_col, start_h, end_h, start_w, end_w, del_h, del_w
    real(real32) :: ratio_h, ratio_w, sum_value
    integer(int64) :: base, out_idx, in_idx

    ratio_h = real(input_h, real32) / real(output_h, real32)
    ratio_w = real(input_w, real32) / real(output_w, real32)
    !$omp target teams distribute parallel do collapse(3) &
    !$omp& private(base, start_h, end_h, start_w, end_w, del_h, del_w, sum_value, i, j, src_row, src_col, in_idx, out_idx)
    do b = 0, batch_size - 1
      do oh = 0, output_h - 1
        do ow = 0, output_w - 1
          base = int(b, int64) * pitch
          start_h = floor_real32(real(oh, real32) * ratio_h)
          end_h = ceil_real32(real(oh + 1, real32) * ratio_h)
          start_w = floor_real32(real(ow, real32) * ratio_w)
          end_w = ceil_real32(real(ow + 1, real32) * ratio_w)
          del_h = end_h - start_h
          del_w = end_w - start_w
          sum_value = 0.0_real32
          do i = 0, del_h - 1
            src_row = start_h + i
            if (src_row >= input_h) cycle
            do j = 0, del_w - 1
              src_col = start_w + j
              if (src_col >= input_w) cycle
              in_idx = base + int(src_row * input_w + src_col, int64)
              sum_value = sum_value + input_tensor(in_idx)
            end do
          end do
          out_idx = base + int((oh + out_h_start) * input_w + (ow + out_w_start), int64)
          output_tensor(out_idx) = sum_value / real(del_h * del_w, real32)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine zoom_out_kernel

  subroutine zoom_out_edge_pad(output_tensor, height, width, pitch, no_padding_h_start, no_padding_w_start, &
      no_padding_h_end, no_padding_w_end, batch_size)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: height, width, no_padding_h_start, no_padding_w_start, no_padding_h_end, no_padding_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow
    integer(int64) :: base, dst, src

    !$omp target teams distribute parallel do collapse(3) private(base, dst, src)
    do b = 0, batch_size - 1
      do oh = 0, height - 1
        do ow = 0, width - 1
          base = int(b, int64) * pitch
          dst = base + int(oh * width + ow, int64)
          src = -1_int64
          if (oh < no_padding_h_start .and. ow >= no_padding_w_start .and. ow < no_padding_w_end) then
            src = base + int(no_padding_h_start * width + ow, int64)
          else if (oh >= no_padding_h_end .and. ow >= no_padding_w_start .and. ow < no_padding_w_end) then
            src = base + int((no_padding_h_end - 1) * width + ow, int64)
          else if (ow < no_padding_w_start .and. oh >= no_padding_h_start .and. oh < no_padding_h_end) then
            src = base + int(oh * width + no_padding_w_start, int64)
          else if (ow >= no_padding_w_end .and. oh >= no_padding_h_start .and. oh < no_padding_h_end) then
            src = base + int(oh * width + (no_padding_w_end - 1), int64)
          else if (oh < no_padding_h_start .and. ow < no_padding_w_start) then
            src = base + int(no_padding_h_start * width + no_padding_w_start, int64)
          else if (oh < no_padding_h_start .and. ow >= no_padding_w_end) then
            src = base + int(no_padding_h_start * width + (no_padding_w_end - 1), int64)
          else if (oh >= no_padding_h_end .and. ow < no_padding_w_start) then
            src = base + int((no_padding_h_end - 1) * width + no_padding_w_start, int64)
          else if (oh >= no_padding_h_end .and. ow >= no_padding_w_end) then
            src = base + int((no_padding_h_end - 1) * width + (no_padding_w_end - 1), int64)
          end if
          if (src >= 0_int64) output_tensor(dst) = output_tensor(src)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine zoom_out_edge_pad

  subroutine zoom_in_reference(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, &
      out_h_start, out_h_end, out_w_start, out_w_end, batch_size)
    real(real32), intent(in) :: input_tensor(0:)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: input_h, input_w, output_h, output_w, out_h_start, out_h_end, out_w_start, out_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow
    do b = 0, batch_size - 1
      do oh = 0, output_h - 1
        do ow = 0, output_w - 1
          if (oh < out_h_start .or. oh >= out_h_end .or. ow < out_w_start .or. ow >= out_w_end) cycle
          call zoom_in_pixel(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, &
            out_h_start, out_w_start, b, oh, ow)
        end do
      end do
    end do
  end subroutine zoom_in_reference

  subroutine zoom_in_pixel(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, out_h_start, out_w_start, b, oh, ow)
    real(real32), intent(in) :: input_tensor(0:)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: input_h, input_w, output_h, output_w, out_h_start, out_w_start, b, oh, ow
    integer(int64), intent(in) :: pitch
    integer :: i, j, src_row, src_col, start_h, end_h, start_w, end_w, del_h, del_w
    real(real32) :: ratio_h, ratio_w, sum_value
    integer(int64) :: base, out_idx, in_idx
    ratio_h = real(input_h, real32) / real(output_h, real32)
    ratio_w = real(input_w, real32) / real(output_w, real32)
    start_h = floor_real32(real(oh, real32) * ratio_h)
    end_h = ceil_real32(real(oh + 1, real32) * ratio_h)
    start_w = floor_real32(real(ow, real32) * ratio_w)
    end_w = ceil_real32(real(ow + 1, real32) * ratio_w)
    del_h = end_h - start_h
    del_w = end_w - start_w
    base = int(b, int64) * pitch
    sum_value = 0.0_real32
    do i = 0, del_h - 1
      src_row = start_h + i
      if (src_row >= input_h) cycle
      do j = 0, del_w - 1
        src_col = start_w + j
        if (src_col >= input_w) cycle
        in_idx = base + int(src_row * input_w + src_col, int64)
        sum_value = sum_value + input_tensor(in_idx)
      end do
    end do
    out_idx = base + int((oh - out_h_start) * input_w + (ow - out_w_start), int64)
    output_tensor(out_idx) = sum_value / real(del_h * del_w, real32)
  end subroutine zoom_in_pixel

  subroutine zoom_out_reference(input_tensor, output_tensor, input_h, input_w, output_h, output_w, pitch, &
      out_h_start, out_h_end, out_w_start, out_w_end, batch_size)
    real(real32), intent(in) :: input_tensor(0:)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: input_h, input_w, output_h, output_w, out_h_start, out_h_end, out_w_start, out_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow, i, j, src_row, src_col, start_h, end_h, start_w, end_w, del_h, del_w
    real(real32) :: ratio_h, ratio_w, sum_value
    integer(int64) :: base, out_idx, in_idx
    ratio_h = real(input_h, real32) / real(output_h, real32)
    ratio_w = real(input_w, real32) / real(output_w, real32)
    do b = 0, batch_size - 1
      base = int(b, int64) * pitch
      do oh = 0, output_h - 1
        do ow = 0, output_w - 1
          start_h = floor_real32(real(oh, real32) * ratio_h)
          end_h = ceil_real32(real(oh + 1, real32) * ratio_h)
          start_w = floor_real32(real(ow, real32) * ratio_w)
          end_w = ceil_real32(real(ow + 1, real32) * ratio_w)
          del_h = end_h - start_h
          del_w = end_w - start_w
          sum_value = 0.0_real32
          do i = 0, del_h - 1
            src_row = start_h + i
            if (src_row >= input_h) cycle
            do j = 0, del_w - 1
              src_col = start_w + j
              if (src_col >= input_w) cycle
              in_idx = base + int(src_row * input_w + src_col, int64)
              sum_value = sum_value + input_tensor(in_idx)
            end do
          end do
          out_idx = base + int((oh + out_h_start) * input_w + (ow + out_w_start), int64)
          output_tensor(out_idx) = sum_value / real(del_h * del_w, real32)
        end do
      end do
    end do
  end subroutine zoom_out_reference

  subroutine zoom_out_edge_pad_reference(output_tensor, height, width, pitch, no_padding_h_start, no_padding_w_start, &
      no_padding_h_end, no_padding_w_end, batch_size)
    real(real32), intent(inout) :: output_tensor(0:)
    integer, intent(in) :: height, width, no_padding_h_start, no_padding_w_start, no_padding_h_end, no_padding_w_end, batch_size
    integer(int64), intent(in) :: pitch
    integer :: b, oh, ow
    integer(int64) :: base
    do b = 0, batch_size - 1
      base = int(b, int64) * pitch
      do oh = 0, height - 1
        do ow = 0, width - 1
          call pad_one(output_tensor, base, height, width, no_padding_h_start, no_padding_w_start, no_padding_h_end, no_padding_w_end, oh, ow)
        end do
      end do
    end do
  end subroutine zoom_out_edge_pad_reference

  subroutine pad_one(output_tensor, base, height, width, no_padding_h_start, no_padding_w_start, no_padding_h_end, no_padding_w_end, oh, ow)
    real(real32), intent(inout) :: output_tensor(0:)
    integer(int64), intent(in) :: base
    integer, intent(in) :: height, width, no_padding_h_start, no_padding_w_start, no_padding_h_end, no_padding_w_end, oh, ow
    integer(int64) :: dst, src
    integer :: unused_height
    unused_height = height
    dst = base + int(oh * width + ow, int64)
    src = -1_int64
    if (oh < no_padding_h_start .and. ow >= no_padding_w_start .and. ow < no_padding_w_end) then
      src = base + int(no_padding_h_start * width + ow, int64)
    else if (oh >= no_padding_h_end .and. ow >= no_padding_w_start .and. ow < no_padding_w_end) then
      src = base + int((no_padding_h_end - 1) * width + ow, int64)
    else if (ow < no_padding_w_start .and. oh >= no_padding_h_start .and. oh < no_padding_h_end) then
      src = base + int(oh * width + no_padding_w_start, int64)
    else if (ow >= no_padding_w_end .and. oh >= no_padding_h_start .and. oh < no_padding_h_end) then
      src = base + int(oh * width + (no_padding_w_end - 1), int64)
    else if (oh < no_padding_h_start .and. ow < no_padding_w_start) then
      src = base + int(no_padding_h_start * width + no_padding_w_start, int64)
    else if (oh < no_padding_h_start .and. ow >= no_padding_w_end) then
      src = base + int(no_padding_h_start * width + (no_padding_w_end - 1), int64)
    else if (oh >= no_padding_h_end .and. ow < no_padding_w_start) then
      src = base + int((no_padding_h_end - 1) * width + no_padding_w_start, int64)
    else if (oh >= no_padding_h_end .and. ow >= no_padding_w_end) then
      src = base + int((no_padding_h_end - 1) * width + (no_padding_w_end - 1), int64)
    end if
    if (src >= 0_int64) output_tensor(dst) = output_tensor(src)
  end subroutine pad_one

end program main
