! SPDX-License-Identifier: CC0-1.0
module doh_kernels
  use, intrinsic :: iso_fortran_env, only : real32, real64, int64
  use omp_lib
  implicit none

contains

  integer function clip_int(x, low, high)
    integer, intent(in) :: x, low, high
    if (x > high) then
      clip_int = high
    else if (x < low) then
      clip_int = low
    else
      clip_int = x
    end if
  end function clip_int

  real(real32) function integ(img, img_rows, img_cols, r_in, c_in, rl, cl)
    real(real32), intent(in) :: img(:)
    integer, intent(in) :: img_rows, img_cols, r_in, c_in, rl, cl
    integer :: r, c, r2, c2
    r = clip_int(r_in, 0, img_rows - 1)
    c = clip_int(c_in, 0, img_cols - 1)
    r2 = clip_int(r + rl, 0, img_rows - 1)
    c2 = clip_int(c + cl, 0, img_cols - 1)
    integ = img(r * img_cols + c + 1) + img(r2 * img_cols + c2 + 1) - &
        img(r * img_cols + c2 + 1) - img(r2 * img_cols + c + 1)
    integ = max(0.0_real32, integ)
  end function integ

  subroutine hessian_matrix_det(img, img_rows, img_cols, sigma, out)
    real(real32), intent(in) :: img(:), sigma
    integer, intent(in) :: img_rows, img_cols
    real(real32), intent(out) :: out(:)
    integer :: tid, r, c, size, b, l, w
    real(real32) :: w_i, tl, br, bl, tr, dxy, mid, side, dxx, dyy
    !$omp target teams distribute parallel do thread_limit(256) private(r, c, size, b, l, w, w_i, tl, br, bl, tr, dxy, mid, side, dxx, dyy)
    do tid = 0, img_rows * img_cols - 1
      r = tid / img_cols
      c = mod(tid, img_cols)
      size = int(3.0_real32 * sigma)
      b = (size - 1) / 2 + 1
      l = size / 3
      w = size
      w_i = 1.0_real32 / real(size * size, real32)

      tl = integ(img, img_rows, img_cols, r - l, c - l, l, l)
      br = integ(img, img_rows, img_cols, r + 1, c + 1, l, l)
      bl = integ(img, img_rows, img_cols, r - l, c + 1, l, l)
      tr = integ(img, img_rows, img_cols, r + 1, c - l, l, l)
      dxy = -(bl + tr - tl - br) * w_i

      mid = integ(img, img_rows, img_cols, r - l + 1, c - l, 2 * l - 1, w)
      side = integ(img, img_rows, img_cols, r - l + 1, c - l / 2, 2 * l - 1, l)
      dxx = -(mid - 3.0_real32 * side) * w_i

      mid = integ(img, img_rows, img_cols, r - l, c - b + 1, w, 2 * b - 1)
      side = integ(img, img_rows, img_cols, r - b / 2, c - b + 1, b, 2 * b - 1)
      dyy = -(mid - 3.0_real32 * side) * w_i

      out(tid + 1) = dxx * dyy - 0.81_real32 * (dxy * dxy)
    end do
    !$omp end target teams distribute parallel do
  end subroutine hessian_matrix_det

  subroutine hessian_matrix_det_cpu(img, img_rows, img_cols, sigma, out)
    real(real32), intent(in) :: img(:), sigma
    integer, intent(in) :: img_rows, img_cols
    real(real32), intent(out) :: out(:)
    integer :: tid, r, c, size, b, l, w
    real(real32) :: w_i, tl, br, bl, tr, dxy, mid, side, dxx, dyy
    do tid = 0, img_rows * img_cols - 1
      r = tid / img_cols
      c = mod(tid, img_cols)
      size = int(3.0_real32 * sigma)
      b = (size - 1) / 2 + 1
      l = size / 3
      w = size
      w_i = 1.0_real32 / real(size * size, real32)

      tl = integ(img, img_rows, img_cols, r - l, c - l, l, l)
      br = integ(img, img_rows, img_cols, r + 1, c + 1, l, l)
      bl = integ(img, img_rows, img_cols, r - l, c + 1, l, l)
      tr = integ(img, img_rows, img_cols, r + 1, c - l, l, l)
      dxy = -(bl + tr - tl - br) * w_i

      mid = integ(img, img_rows, img_cols, r - l + 1, c - l, 2 * l - 1, w)
      side = integ(img, img_rows, img_cols, r - l + 1, c - l / 2, 2 * l - 1, l)
      dxx = -(mid - 3.0_real32 * side) * w_i

      mid = integ(img, img_rows, img_cols, r - l, c - b + 1, w, 2 * b - 1)
      side = integ(img, img_rows, img_cols, r - b / 2, c - b + 1, b, 2 * b - 1)
      dyy = -(mid - 3.0_real32 * side) * w_i

      out(tid + 1) = dxx * dyy - 0.81_real32 * (dxy * dxy)
    end do
  end subroutine hessian_matrix_det_cpu

end module doh_kernels

program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64, int64
  use omp_lib
  use doh_kernels
  implicit none

  integer :: h, w, repeat, img_size, i, j, y, x
  real(real32), allocatable :: input_img(:), integral_img(:), output_img(:), reference_img(:)
  real(real32) :: sigma, s
  real(real64) :: start_time, end_time, elapsed_us, checksum

  interface
    subroutine doh_fill_input(input_img, img_size) bind(C, name='doh_fill_input')
      import :: c_float, c_int
      real(c_float), intent(out) :: input_img(*)
      integer(c_int), value :: img_size
    end subroutine doh_fill_input
  end interface

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <height> <width> <repeat>'
    stop 1
  end if

  h = read_arg(1)
  w = read_arg(2)
  repeat = read_arg(3)
  img_size = h * w
  sigma = 4.0_real32

  allocate(input_img(img_size), integral_img(img_size), output_img(img_size), reference_img(img_size))
  call fill_input(input_img)

  print '(A)', 'Integrating the input image may take a while...'
  do i = 0, h - 1
    do j = 0, w - 1
      s = 0.0_real32
      do y = 0, i
        do x = 0, j
          s = s + input_img(y * w + x + 1)
        end do
      end do
      integral_img(i * w + j + 1) = s
    end do
  end do

  !$omp target data map(to: integral_img(1:img_size)) map(from: output_img(1:img_size))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call hessian_matrix_det(integral_img, h, w, sigma, output_img)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  call hessian_matrix_det_cpu(integral_img, h, w, sigma, reference_img)
  if (maxval(abs(output_img - reference_img)) > 1.0e-3_real32) then
    print '(A)', 'Internal validation failed'
  end if

  checksum = sum(real(output_img, real64))
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time : ', elapsed_us, ' (us)'
  print '(A,F8.6)', 'Kernel checksum: ', checksum

  deallocate(input_img, integral_img, output_img, reference_img)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_input(values)
    real(real32), intent(out) :: values(:)
    call doh_fill_input(values, int(size(values), c_int))
  end subroutine fill_input

end program main
