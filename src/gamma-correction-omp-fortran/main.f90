! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg
  integer :: width, height, block_size, repeat
  integer :: image_size, i, iter
  integer, parameter :: b = 1, g = 2, r = 3, a = 4, channels = 4
  integer, allocatable :: image(:,:), reference(:,:), pixel(:,:)
  real(real64) :: start_time, total_time
  character(len=32) :: avg_time_text

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <image width> <image height> <block size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) block_size
  call get_command_argument(4, arg); read(arg, *) repeat

  if (width <= 0 .or. height <= 0 .or. block_size <= 0 .or. repeat <= 0) stop 1

  image_size = width * height
  allocate(image(channels, image_size), reference(channels, image_size), pixel(channels, image_size))

  call fill_fractal(width, height, image)
  do i = 1, image_size
    call gamma_pixel(reference(:, i), image(:, i))
  end do

  pixel = image
  total_time = 0.0_real64

  !$omp target data map(from: pixel(1:channels, 1:image_size))
  do iter = 1, repeat
    pixel = image
    !$omp target update to(pixel(1:channels, 1:image_size))

    start_time = omp_get_wtime()
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, image_size
      call gamma_pixel(pixel(:, i), pixel(:, i))
    end do
    !$omp end target teams distribute parallel do
    total_time = total_time + (omp_get_wtime() - start_time)
  end do
  !$omp end target data

  write(avg_time_text, '(F20.6)') total_time / real(repeat, real64)
  avg_time_text = adjustl(avg_time_text)
  if (avg_time_text(1:1) == '.') avg_time_text = '0' // trim(avg_time_text)
  write(*,'(A,A,A)') 'Average kernel execution time ', trim(avg_time_text), ' (s)'

  if (all(pixel == reference)) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(image, reference, pixel)

contains

  subroutine fill_fractal(width, height, image)
    integer, intent(in) :: width, height
    integer, intent(out) :: image(:,:)
    integer :: x, y, idx
    integer :: fractal_pixel

    do y = 0, height - 1
      do x = 0, width - 1
        idx = y * width + x + 1
        fractal_pixel = int(fractal_value(x, y, width, height))
        if (fractal_pixel < 0) fractal_pixel = 0
        if (fractal_pixel > 255) fractal_pixel = 255
        image(b, idx) = fractal_pixel
        image(g, idx) = fractal_pixel
        image(r, idx) = fractal_pixel
        image(a, idx) = fractal_pixel
      end do
    end do
  end subroutine fill_fractal

  real(real64) function fractal_value(x, y, width, height)
    integer, intent(in) :: x, y, width, height
    integer :: iter
    real(real64) :: fx, fy, res, nx, ny, val
    real(real64), parameter :: cx = -0.7436_real64
    real(real64), parameter :: cy = 0.1319_real64
    real(real64), parameter :: magn = 2000000.0_real64

    fx = (real(x, real64) - real(width, real64) / 2.0_real64) * (1.0_real64 / magn) + cx
    fy = (real(y, real64) - real(height, real64) / 2.0_real64) * (1.0_real64 / magn) + cy

    res = 0.0_real64
    nx = 0.0_real64
    ny = 0.0_real64
    val = 0.0_real64
    iter = 0
    do while (nx * nx + ny * ny <= 4.0_real64 .and. iter < 1000)
      val = nx * nx - ny * ny + fx
      ny = 2.0_real64 * nx * ny + fy
      nx = val
      res = res + exp(-sqrt(nx * nx + ny * ny))
      iter = iter + 1
    end do

    fractal_value = res
  end function fractal_value

  subroutine gamma_pixel(pixel, input_pixel)
    integer, intent(out) :: pixel(:)
    integer, intent(in) :: input_pixel(:)
    integer :: gamma_value
    real(real32) :: v

    v = (0.3_real32 * real(input_pixel(r), real32) + &
         0.59_real32 * real(input_pixel(g), real32) + &
         0.11_real32 * real(input_pixel(b), real32)) / 255.0_real32
    gamma_value = int(255.0_real32 * v * v)
    if (gamma_value > 255) gamma_value = 255
    pixel(b) = gamma_value
    pixel(g) = gamma_value
    pixel(r) = gamma_value
    pixel(a) = gamma_value
  end subroutine gamma_pixel

end program main
