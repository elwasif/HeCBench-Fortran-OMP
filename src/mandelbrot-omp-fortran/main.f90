! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: row_size = 1080
  integer, parameter :: col_size = 1920
  integer, parameter :: max_iterations = 100
  integer, parameter :: threads_per_block_x = 16
  integer, parameter :: threads_per_block_y = 16
  integer, parameter :: image_size = row_size * col_size

  character(len=256) :: arg0, arg1
  integer :: repetitions, rep
  integer, allocatable :: parallel_data(:), serial_data(:)
  real(real64) :: parallel_start, parallel_elapsed
  real(real64) :: kernel_time, serial_start, serial_elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    call usage(trim(arg0))
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repetitions
  if (repetitions <= 0) call usage(trim(arg0))

  allocate(parallel_data(image_size), serial_data(image_size))

  kernel_time = evaluate_parallel(parallel_data)
  kernel_time = 0.0_real64
  parallel_start = omp_get_wtime()
  do rep = 1, repetitions
    kernel_time = kernel_time + evaluate_parallel(parallel_data)
  end do
  parallel_elapsed = omp_get_wtime() - parallel_start

  call print_mandel()

  serial_start = omp_get_wtime()
  call evaluate_serial(serial_data)
  serial_elapsed = omp_get_wtime() - serial_start

  write(*,'(A20,ES12.6,A)') 'serial time: ', serial_elapsed, 's'
  write(*,'(A,ES12.6,A)') 'Average parallel time: ', &
    (parallel_elapsed / real(repetitions, real64)) * 1.0e3_real64, ' ms'
  write(*,'(A,ES12.6,A)') 'Average kernel execution time: ', &
    (kernel_time / real(repetitions, real64)) * 1.0e3_real64, ' ms'

  call verify_parallel(parallel_data, serial_data)
  write(*,'(A)') 'Success'

  deallocate(parallel_data, serial_data)

contains

  subroutine usage(program_name)
    character(len=*), intent(in) :: program_name

    write(*,'(A)') ' Incorrect parameters'
    write(*,'(A)', advance='no') ' Usage: '
    write(*,'(A,A)') trim(program_name), ' <repeat>'
    write(*,*)
    stop 255
  end subroutine usage

  subroutine print_mandel()
    integer :: i, j

    if (row_size > 128 .or. col_size > 128) then
      write(*,'(A)') 'No Print() output due to size too large'
      return
    end if

    do i = 0, row_size - 1
      do j = 0, col_size - 1
        write(*,'(A)', advance='no') ' '
      end do
      write(*,*)
    end do
  end subroutine print_mandel

  subroutine evaluate_serial(data)
    integer, intent(out) :: data(:)
    integer :: i, j

    do i = 0, row_size - 1
      do j = 0, col_size - 1
        data(i * col_size + j + 1) = mandel_point(scale_row(i), scale_col(j))
      end do
    end do
  end subroutine evaluate_serial

  real(real64) function evaluate_parallel(data)
    integer, intent(out) :: data(:)
    integer :: i, j
    real(real64) :: start_time

    !$omp target data map(from: data(1:image_size))
    start_time = omp_get_wtime()
    !$omp target teams distribute parallel do simd collapse(2) &
    !$omp& thread_limit(threads_per_block_x * threads_per_block_y)
    do i = 0, row_size - 1
      do j = 0, col_size - 1
        data(i * col_size + j + 1) = mandel_point(scale_row(i), scale_col(j))
      end do
    end do
    !$omp end target teams distribute parallel do simd
    evaluate_parallel = omp_get_wtime() - start_time
    !$omp end target data
  end function evaluate_parallel

  subroutine verify_parallel(parallel_data, serial_data)
    integer, intent(in) :: parallel_data(:), serial_data(:)
    integer :: idx, diff
    real(real64) :: ratio

    diff = 0
    do idx = 1, image_size
      if (parallel_data(idx) /= serial_data(idx)) diff = diff + 1
    end do

    ratio = real(diff, real64) / real(image_size, real64)
    if (ratio > 0.05_real64) then
      write(*,'(A)') 'Fail verification - diff larger than tolerance'
      stop 1
    end if
  end subroutine verify_parallel

  real(real32) function scale_row(i)
    integer, intent(in) :: i

    scale_row = -1.5_real32 + real(i, real32) * (2.0_real32 / real(row_size, real32))
  end function scale_row

  real(real32) function scale_col(i)
    integer, intent(in) :: i

    scale_col = -1.0_real32 + real(i, real32) * (2.0_real32 / real(col_size, real32))
  end function scale_col

  integer function mandel_point(c_real, c_imag)
    real(real32), intent(in) :: c_real, c_imag
    integer :: iter
    real(real32) :: z_real, z_imag, r, im

    mandel_point = 0
    z_real = 0.0_real32
    z_imag = 0.0_real32
    do iter = 1, max_iterations
      r = z_real
      im = z_imag
      if ((r * r + im * im) >= 4.0_real32) exit
      z_real = r * r - im * im + c_real
      z_imag = 2.0_real32 * r * im + c_imag
      mandel_point = mandel_point + 1
    end do
  end function mandel_point

end program main
