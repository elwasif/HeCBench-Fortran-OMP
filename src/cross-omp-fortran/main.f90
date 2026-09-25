! SPDX-License-Identifier: CC0-1.0
program cross_main
  use, intrinsic :: iso_c_binding, only : c_double, c_float, c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  interface
    subroutine fill_cross_inputs_float(num_elems, a, b) bind(C, name="fill_cross_inputs_float")
      import :: c_float, c_int
      integer(c_int), value :: num_elems
      real(c_float) :: a(*), b(*)
    end subroutine fill_cross_inputs_float

    subroutine fill_cross_inputs_double(num_elems, a, b) bind(C, name="fill_cross_inputs_double")
      import :: c_double, c_int
      integer(c_int), value :: num_elems
      real(c_double) :: a(*), b(*)
    end subroutine fill_cross_inputs_double
  end interface

  integer :: nrows, repeat
  character(len=64) :: arg

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <number of rows in a 2D tensor> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nrows
  call get_command_argument(2, arg)
  read(arg, *) repeat

  write(*,'("=========== Data type is FP32 ==========")')
  call eval_real32(nrows, repeat)

  write(*,'("=========== Data type is FP64 ==========")')
  call eval_real64(nrows, repeat)

contains

  subroutine eval_real32(nrows, repeat)
    integer, intent(in) :: nrows, repeat
    integer :: num_elems, i
    real(real32), allocatable :: a(:), b(:), o(:), o2(:), o3(:)
    real(real64) :: start_time, elapsed_us
    logical :: ok

    num_elems = nrows * 3
    allocate(a(0:num_elems-1), b(0:num_elems-1), o(0:num_elems-1), o2(0:num_elems-1), o3(0:num_elems-1))
    call fill_cross_inputs_float(num_elems, a, b)

    !$omp target data map(to: a, b) map(from: o, o2, o3)
    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross1_real32(nrows, o, a, b, 1, 1, 1)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross1 kernel: ",F0.6," (us)")') elapsed_us

    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross2_real32(nrows, o2, a, b, 1, 1, 1)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross2 kernel: ",F0.6," (us)")') elapsed_us

    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross3_real32(nrows, o3, a, b)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross3 kernel: ",F0.6," (us)")') elapsed_us
    !$omp end target data

    ok = .true.
    do i = 0, num_elems - 1
      if (abs(o(i) - o2(i)) > 1.0e-3_real32 .or. abs(o(i) - o3(i)) > 1.0e-3_real32) then
        ok = .false.
        exit
      end if
    end do
    write(*,'(A)') merge("PASS", "FAIL", ok)
  end subroutine eval_real32

  subroutine eval_real64(nrows, repeat)
    integer, intent(in) :: nrows, repeat
    integer :: num_elems, i
    real(real64), allocatable :: a(:), b(:), o(:), o2(:), o3(:)
    real(real64) :: start_time, elapsed_us
    logical :: ok

    num_elems = nrows * 3
    allocate(a(0:num_elems-1), b(0:num_elems-1), o(0:num_elems-1), o2(0:num_elems-1), o3(0:num_elems-1))
    call fill_cross_inputs_double(num_elems, a, b)

    !$omp target data map(to: a, b) map(from: o, o2, o3)
    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross1_real64(nrows, o, a, b, 1, 1, 1)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross1 kernel: ",F0.6," (us)")') elapsed_us

    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross2_real64(nrows, o2, a, b, 1, 1, 1)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross2 kernel: ",F0.6," (us)")') elapsed_us

    start_time = omp_get_wtime()
    do i = 1, repeat
      call cross3_real64(nrows, o3, a, b)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'("Average execution time of cross3 kernel: ",F0.6," (us)")') elapsed_us
    !$omp end target data

    ok = .true.
    do i = 0, num_elems - 1
      if (abs(o(i) - o2(i)) > 1.0e-3_real64 .or. abs(o(i) - o3(i)) > 1.0e-3_real64) then
        ok = .false.
        exit
      end if
    end do
    write(*,'(A)') merge("PASS", "FAIL", ok)
  end subroutine eval_real64

  subroutine cross1_real32(numel, out, x1, x2, ostride, x1stride, x2stride)
    integer, intent(in) :: numel, ostride, x1stride, x2stride
    real(real32), intent(out) :: out(0:)
    real(real32), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      out(out_base + 0 * ostride) = x1(x1_base + 1 * x1stride) * x2(x2_base + 2 * x2stride) - &
          x1(x1_base + 2 * x1stride) * x2(x2_base + 1 * x2stride)
      out(out_base + 1 * ostride) = x1(x1_base + 2 * x1stride) * x2(x2_base + 0 * x2stride) - &
          x1(x1_base + 0 * x1stride) * x2(x2_base + 2 * x2stride)
      out(out_base + 2 * ostride) = x1(x1_base + 0 * x1stride) * x2(x2_base + 1 * x2stride) - &
          x1(x1_base + 1 * x1stride) * x2(x2_base + 0 * x2stride)
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross1_real32

  subroutine cross2_real32(numel, out, x1, x2, ostride, x1stride, x2stride)
    integer, intent(in) :: numel, ostride, x1stride, x2stride
    real(real32), intent(out) :: out(0:)
    real(real32), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base
    real(real32) :: x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2
    real(real32) :: val0, val1, val2

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base, x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2, val0, val1, val2)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      x1_c0 = x1(x1_base + 0 * x1stride)
      x1_c1 = x1(x1_base + 1 * x1stride)
      x1_c2 = x1(x1_base + 2 * x1stride)
      x2_c0 = x2(x2_base + 0 * x2stride)
      x2_c1 = x2(x2_base + 1 * x2stride)
      x2_c2 = x2(x2_base + 2 * x2stride)
      val0 = x1_c1 * x2_c2 - x1_c2 * x2_c1
      val1 = x1_c2 * x2_c0 - x1_c0 * x2_c2
      val2 = x1_c0 * x2_c1 - x1_c1 * x2_c0
      out(out_base + 0 * ostride) = val0
      out(out_base + 1 * ostride) = val1
      out(out_base + 2 * ostride) = val2
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross2_real32

  subroutine cross3_real32(numel, out, x1, x2)
    integer, intent(in) :: numel
    real(real32), intent(out) :: out(0:)
    real(real32), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base
    real(real32) :: x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2
    real(real32) :: val0, val1, val2

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base, x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2, val0, val1, val2)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      x1_c0 = x1(x1_base)
      x1_c1 = x1(x1_base + 1)
      x1_c2 = x1(x1_base + 2)
      x2_c0 = x2(x2_base)
      x2_c1 = x2(x2_base + 1)
      x2_c2 = x2(x2_base + 2)
      val0 = x1_c1 * x2_c2 - x1_c2 * x2_c1
      val1 = x1_c2 * x2_c0 - x1_c0 * x2_c2
      val2 = x1_c0 * x2_c1 - x1_c1 * x2_c0
      out(out_base) = val0
      out(out_base + 1) = val1
      out(out_base + 2) = val2
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross3_real32

  subroutine cross1_real64(numel, out, x1, x2, ostride, x1stride, x2stride)
    integer, intent(in) :: numel, ostride, x1stride, x2stride
    real(real64), intent(out) :: out(0:)
    real(real64), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      out(out_base) = x1(x1_base + 1) * x2(x2_base + 2) - x1(x1_base + 2) * x2(x2_base + 1)
      out(out_base + 1) = x1(x1_base + 2) * x2(x2_base) - x1(x1_base) * x2(x2_base + 2)
      out(out_base + 2) = x1(x1_base) * x2(x2_base + 1) - x1(x1_base + 1) * x2(x2_base)
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross1_real64

  subroutine cross2_real64(numel, out, x1, x2, ostride, x1stride, x2stride)
    integer, intent(in) :: numel, ostride, x1stride, x2stride
    real(real64), intent(out) :: out(0:)
    real(real64), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base
    real(real64) :: x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2
    real(real64) :: val0, val1, val2

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base, x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2, val0, val1, val2)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      x1_c0 = x1(x1_base + 0 * x1stride)
      x1_c1 = x1(x1_base + 1 * x1stride)
      x1_c2 = x1(x1_base + 2 * x1stride)
      x2_c0 = x2(x2_base + 0 * x2stride)
      x2_c1 = x2(x2_base + 1 * x2stride)
      x2_c2 = x2(x2_base + 2 * x2stride)
      val0 = x1_c1 * x2_c2 - x1_c2 * x2_c1
      val1 = x1_c2 * x2_c0 - x1_c0 * x2_c2
      val2 = x1_c0 * x2_c1 - x1_c1 * x2_c0
      out(out_base + 0 * ostride) = val0
      out(out_base + 1 * ostride) = val1
      out(out_base + 2 * ostride) = val2
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross2_real64

  subroutine cross3_real64(numel, out, x1, x2)
    integer, intent(in) :: numel
    real(real64), intent(out) :: out(0:)
    real(real64), intent(in) :: x1(0:), x2(0:)
    integer :: i, out_base, x1_base, x2_base
    real(real64) :: x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2
    real(real64) :: val0, val1, val2

    !$omp target teams distribute parallel do thread_limit(256) private(out_base, x1_base, x2_base, x1_c0, x1_c1, x1_c2, x2_c0, x2_c1, x2_c2, val0, val1, val2)
    do i = 0, numel - 1
      out_base = 3 * i
      x1_base = 3 * i
      x2_base = 3 * i
      x1_c0 = x1(x1_base)
      x1_c1 = x1(x1_base + 1)
      x1_c2 = x1(x1_base + 2)
      x2_c0 = x2(x2_base)
      x2_c1 = x2(x2_base + 1)
      x2_c2 = x2(x2_base + 2)
      val0 = x1_c1 * x2_c2 - x1_c2 * x2_c1
      val1 = x1_c2 * x2_c0 - x1_c0 * x2_c2
      val2 = x1_c0 * x2_c1 - x1_c1 * x2_c0
      out(out_base) = val0
      out(out_base + 1) = val1
      out(out_base + 2) = val2
    end do
    !$omp end target teams distribute parallel do
  end subroutine cross3_real64

end program cross_main
