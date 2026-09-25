! SPDX-License-Identifier: CC0-1.0
program tsa_main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer :: width, height, repeat
  character(len=64) :: arg

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <matrix width> <matrix height> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) width
  call get_command_argument(2, arg)
  read(arg, *) height
  call get_command_argument(3, arg)
  read(arg, *) repeat

  write(*,'("TSA in float32")')
  call tsa32(width, height, repeat)
  write(*,*)
  write(*,'("TSA in float64")')
  call tsa64(width, height, repeat)

contains

  subroutine tsa32(width, height, repeat)
    integer, intent(in) :: width, height, repeat
    integer :: numel, i
    real(real32), allocatable :: p_real(:), p_imag(:), p2_real(:), p2_imag(:), h_real(:), h_imag(:)
    real(real32) :: a, b
    real(real64) :: start_time, elapsed_us
    logical :: ok

    numel = width * height
    allocate(p_real(0:numel-1), p_imag(0:numel-1), p2_real(0:numel-1), p2_imag(0:numel-1), &
             h_real(0:numel-1), h_imag(0:numel-1))
    call init_p32(p_real, p_imag, width, height)
    h_real = p_real
    h_imag = p_imag
    a = cos(0.02_real32)
    b = sin(0.02_real32)

    call reference32(h_real, h_imag, a, b, width, height, repeat)

    !$omp target data map(to: p_real, p_imag) map(alloc: p2_real, p2_imag)
    start_time = omp_get_wtime()
    do i = 1, repeat
      if (mod(i, 2) == 1) then
        call trotter_kernel32(p_real, p_imag, p2_real, p2_imag, a, b, width, height)
      else
        call trotter_kernel32(p2_real, p2_imag, p_real, p_imag, a, b, width, height)
      end if
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    if (mod(repeat, 2) == 1) then
      !$omp target update from(p2_real, p2_imag)
    else
      !$omp target update from(p_real, p_imag)
    end if
    !$omp end target data

    write(*,'("Average kernel execution time: ",F0.6," (us)")') elapsed_us
    if (mod(repeat, 2) == 1) then
      p_real = p2_real
      p_imag = p2_imag
    end if
    ok = maxval(abs(p_real - h_real)) <= 1.0e-3_real32 .and. maxval(abs(p_imag - h_imag)) <= 1.0e-3_real32
    write(*,'(A)') merge("PASS", "FAIL", ok)
    if (.not. ok) stop 1
  end subroutine tsa32

  subroutine tsa64(width, height, repeat)
    integer, intent(in) :: width, height, repeat
    integer :: numel, i
    real(real64), allocatable :: p_real(:), p_imag(:), p2_real(:), p2_imag(:), h_real(:), h_imag(:)
    real(real64) :: a, b
    real(real64) :: start_time, elapsed_us
    logical :: ok

    numel = width * height
    allocate(p_real(0:numel-1), p_imag(0:numel-1), p2_real(0:numel-1), p2_imag(0:numel-1), &
             h_real(0:numel-1), h_imag(0:numel-1))
    call init_p64(p_real, p_imag, width, height)
    h_real = p_real
    h_imag = p_imag
    a = cos(0.02_real64)
    b = sin(0.02_real64)

    call reference64(h_real, h_imag, a, b, width, height, repeat)

    !$omp target data map(to: p_real, p_imag) map(alloc: p2_real, p2_imag)
    start_time = omp_get_wtime()
    do i = 1, repeat
      if (mod(i, 2) == 1) then
        call trotter_kernel64(p_real, p_imag, p2_real, p2_imag, a, b, width, height)
      else
        call trotter_kernel64(p2_real, p2_imag, p_real, p_imag, a, b, width, height)
      end if
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    if (mod(repeat, 2) == 1) then
      !$omp target update from(p2_real, p2_imag)
    else
      !$omp target update from(p_real, p_imag)
    end if
    !$omp end target data

    write(*,'("Average kernel execution time: ",F0.6," (us)")') elapsed_us
    if (mod(repeat, 2) == 1) then
      p_real = p2_real
      p_imag = p2_imag
    end if
    ok = maxval(abs(p_real - h_real)) <= 1.0e-3_real64 .and. maxval(abs(p_imag - h_imag)) <= 1.0e-3_real64
    write(*,'(A)') merge("PASS", "FAIL", ok)
    if (.not. ok) stop 1
  end subroutine tsa64

  subroutine init_p32(p_real, p_imag, width, height)
    real(real32), intent(out) :: p_real(0:), p_imag(0:)
    integer, intent(in) :: width, height
    integer :: i, j, idx
    real(real32) :: s, mag, phase

    s = 64.0_real32
    do j = 1, height
      do i = 1, width
        idx = (j - 1) * width + i - 1
        mag = exp(-(((real(i, real32) - 180.0_real32) ** 2 + &
                     (real(j, real32) - 300.0_real32) ** 2) / (2.0_real32 * s ** 2)))
        phase = 0.4_real32 * (real(i + j, real32) - 480.0_real32)
        p_real(idx) = mag * cos(phase)
        p_imag(idx) = mag * sin(phase)
      end do
    end do
  end subroutine init_p32

  subroutine init_p64(p_real, p_imag, width, height)
    real(real64), intent(out) :: p_real(0:), p_imag(0:)
    integer, intent(in) :: width, height
    integer :: i, j, idx
    real(real64) :: s, mag, phase

    s = 64.0_real64
    do j = 1, height
      do i = 1, width
        idx = (j - 1) * width + i - 1
        mag = exp(-(((real(i, real64) - 180.0_real64) ** 2 + &
                     (real(j, real64) - 300.0_real64) ** 2) / (2.0_real64 * s ** 2)))
        phase = 0.4_real64 * (real(i + j, real64) - 480.0_real64)
        p_real(idx) = mag * cos(phase)
        p_imag(idx) = mag * sin(phase)
      end do
    end do
  end subroutine init_p64

  subroutine reference32(pr, pi, a, b, width, height, repeat)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height, repeat
    integer :: i

    do i = 1, repeat
      call kernel1_32(pr, pi, a, b, width, height)
      call kernel2_32(pr, pi, a, b, width, height)
      call kernel3_32(pr, pi, a, b, width, height)
      call kernel4_32(pr, pi, a, b, width, height)
      call kernel4_32(pr, pi, a, b, width, height)
      call kernel3_32(pr, pi, a, b, width, height)
      call kernel2_32(pr, pi, a, b, width, height)
      call kernel1_32(pr, pi, a, b, width, height)
    end do
  end subroutine reference32

  subroutine reference64(pr, pi, a, b, width, height, repeat)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height, repeat
    integer :: i

    do i = 1, repeat
      call kernel1_64(pr, pi, a, b, width, height)
      call kernel2_64(pr, pi, a, b, width, height)
      call kernel3_64(pr, pi, a, b, width, height)
      call kernel4_64(pr, pi, a, b, width, height)
      call kernel4_64(pr, pi, a, b, width, height)
      call kernel3_64(pr, pi, a, b, width, height)
      call kernel2_64(pr, pi, a, b, width, height)
      call kernel1_64(pr, pi, a, b, width, height)
    end do
  end subroutine reference64

  subroutine trotter_sequence32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    call kernel1_32_device(pr, pi, a, b, width, height)
    call kernel2_32_device(pr, pi, a, b, width, height)
    call kernel3_32_device(pr, pi, a, b, width, height)
    call kernel4_32_device(pr, pi, a, b, width, height)
    call kernel4_32_device(pr, pi, a, b, width, height)
    call kernel3_32_device(pr, pi, a, b, width, height)
    call kernel2_32_device(pr, pi, a, b, width, height)
    call kernel1_32_device(pr, pi, a, b, width, height)
  end subroutine trotter_sequence32

  subroutine trotter_sequence64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    call kernel1_64_device(pr, pi, a, b, width, height)
    call kernel2_64_device(pr, pi, a, b, width, height)
    call kernel3_64_device(pr, pi, a, b, width, height)
    call kernel4_64_device(pr, pi, a, b, width, height)
    call kernel4_64_device(pr, pi, a, b, width, height)
    call kernel3_64_device(pr, pi, a, b, width, height)
    call kernel2_64_device(pr, pi, a, b, width, height)
    call kernel1_64_device(pr, pi, a, b, width, height)
  end subroutine trotter_sequence64

  include 'tsa_kernels.inc'

end program tsa_main
