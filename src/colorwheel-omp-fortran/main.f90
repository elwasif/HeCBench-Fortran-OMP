! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: ry = 15
  integer, parameter :: yg = 6
  integer, parameter :: gc = 4
  integer, parameter :: cb = 11
  integer, parameter :: bm = 13
  integer, parameter :: mr = 6
  integer, parameter :: maxcols = ry + yg + gc + cb + bm + mr

  character(len=256) :: arg0, arg
  integer :: size, repeat, half_size, img_size
  integer :: x, y, idx, i, fail, max_error, err
  real(real32) :: truerange, range, fx, fy
  real(real64) :: start_time, elapsed
  integer(int8), allocatable :: pix(:), res(:), d_pix(:)

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <range> <size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) truerange
  call get_command_argument(2, arg); read(arg, *) size
  call get_command_argument(3, arg); read(arg, *) repeat
  if (truerange == 0.0_real32 .or. size <= 0 .or. repeat <= 0) stop 1

  range = 1.04_real32 * truerange
  half_size = size / 2
  if (half_size <= 0) stop 1

  img_size = size * size * 3
  allocate(pix(img_size), res(img_size), d_pix(img_size))
  pix = 0
  res = 0
  d_pix = 0

  do y = 0, size - 1
    do x = 0, size - 1
      fx = real(x, real32) / real(half_size, real32) * range - range
      fy = real(y, real32) / real(half_size, real32) * range - range
      if (x == half_size .or. y == half_size) cycle
      idx = (y * size + x) * 3 + 1
      call compute_color(fx / truerange, fy / truerange, pix(idx:idx + 2))
    end do
  end do

  write(*,'(A)') 'Start execution on a device'

  !$omp target data map(tofrom: d_pix(1:img_size))
  start_time = omp_get_wtime()
  do i = 1, repeat
    !$omp target teams distribute parallel do collapse(2) private(fx, fy, idx)
    do y = 0, size - 1
      do x = 0, size - 1
        fx = real(x, real32) / real(half_size, real32) * range - range
        fy = real(y, real32) / real(half_size, real32) * range - range
        if (x /= half_size .and. y /= half_size) then
          idx = (y * size + x) * 3 + 1
          call compute_color(fx / truerange, fy / truerange, d_pix(idx:idx + 2))
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time : ', &
    (elapsed * 1000.0_real64) / real(repeat, real64), ' (ms)'

  fail = 0
  max_error = 0
  do i = 1, img_size
    err = abs(ubyte_value(res(i)) - ubyte_value(pix(i)))
    if (err > 1) then
      fail = 1
      if (err > max_error) max_error = err
    end if
  end do

  if (fail /= 0) then
    write(*,'(A,I0)') 'Verification failed. Maximum error between host and device results: ', max_error
  else
    write(*,'(A)') 'PASS'
  end if

  deallocate(pix, res, d_pix)

contains

  subroutine compute_color(fx, fy, pix)
    real(real32), intent(in) :: fx, fy
    integer(int8), intent(out) :: pix(3)
    integer :: cw(0:maxcols - 1, 0:2)
    integer :: k, j, b, k0, k1
    real(real32) :: rad, a, fk, f, col0, col1, col
    real(real32), parameter :: pi = acos(-1.0_real32)

    k = 0
    do j = 0, ry - 1
      call set_col(cw, 255, 255 * j / ry, 0, k)
      k = k + 1
    end do
    do j = 0, yg - 1
      call set_col(cw, 255 - 255 * j / yg, 255, 0, k)
      k = k + 1
    end do
    do j = 0, gc - 1
      call set_col(cw, 0, 255, 255 * j / gc, k)
      k = k + 1
    end do
    do j = 0, cb - 1
      call set_col(cw, 0, 255 - 255 * j / cb, 255, k)
      k = k + 1
    end do
    do j = 0, bm - 1
      call set_col(cw, 255 * j / bm, 0, 255, k)
      k = k + 1
    end do
    do j = 0, mr - 1
      call set_col(cw, 255, 0, 255 - 255 * j / mr, k)
      k = k + 1
    end do

    rad = sqrt(fx * fx + fy * fy)
    a = atan2(-fy, -fx) / pi
    fk = (a + 1.0_real32) / 2.0_real32 * real(maxcols - 1, real32)
    k0 = int(fk)
    k1 = modulo(k0 + 1, maxcols)
    f = fk - real(k0, real32)

    do b = 0, 2
      col0 = real(cw(k0, b), real32) / 255.0_real32
      col1 = real(cw(k1, b), real32) / 255.0_real32
      col = (1.0_real32 - f) * col0 + f * col1
      if (rad <= 1.0_real32) then
        col = 1.0_real32 - rad * (1.0_real32 - col)
      else
        col = col * 0.75_real32
      end if
      pix(3 - b) = to_ubyte(int(255.0_real32 * col))
    end do
  end subroutine compute_color

  pure integer(int8) function to_ubyte(value) result(byte)
    integer, intent(in) :: value
    if (value <= 127) then
      byte = int(value, int8)
    else
      byte = int(value - 256, int8)
    end if
  end function to_ubyte

  pure integer function ubyte_value(byte) result(value)
    integer(int8), intent(in) :: byte
    value = int(byte)
    if (value < 0) value = value + 256
  end function ubyte_value

  subroutine set_col(cw, r, g, b, k)
    integer, intent(inout) :: cw(0:maxcols - 1, 0:2)
    integer, intent(in) :: r, g, b, k
    cw(k, 0) = r
    cw(k, 1) = g
    cw(k, 2) = b
  end subroutine set_col

end program main
