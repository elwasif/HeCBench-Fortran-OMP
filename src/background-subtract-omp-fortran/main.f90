! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int, c_int64_t
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  interface
    subroutine c_seed_rng(seed) bind(C, name="bs_seed_rng")
      import :: c_int64_t
      integer(c_int64_t), value :: seed
    end subroutine c_seed_rng

    function c_next_rand_byte() bind(C, name="bs_next_rand_byte") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_next_rand_byte
  end interface

  character(len=256) :: arg0, arg
  integer :: width, height, merged, repeat, img_size, iter, j
  integer, pointer :: img(:), img1(:), img2(:), tmp(:)
  integer, allocatable :: bn(:), bn_ref(:), mp(:), tn(:), tn_ref(:)
  integer :: max_error
  integer(int64) :: time_ticks
  real(real64) :: start_time, kernel_time_us

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <image width> <image height> <merge> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) merged
  call get_command_argument(4, arg); read(arg, *) repeat
  if (width <= 0 .or. height <= 0 .or. repeat <= 0) stop 1

  img_size = width * height
  allocate(img(img_size), img1(img_size), img2(img_size))
  allocate(bn(img_size), bn_ref(img_size), mp(img_size), tn(img_size), tn_ref(img_size))

  call seed_rng(123_int64)
  do j = 1, img_size
    bn(j) = next_rand_byte()
    bn_ref(j) = bn(j)
    tn(j) = 128
    tn_ref(j) = 128
  end do
  img = 0; img1 = 0; img2 = 0; mp = 0
  time_ticks = 0_int64

  !$omp target data map(tofrom: bn(1:img_size), tn(1:img_size)) &
  !$omp& map(alloc: mp(1:img_size), img(1:img_size), img1(1:img_size), img2(1:img_size))
    do iter = 1, repeat
      do j = 1, img_size
        img(j) = next_rand_byte()
      end do

      !$omp target update to(img(1:img_size))

      tmp => img2
      img2 => img1
      img1 => img
      img => tmp

      if (iter >= 3) then
        start_time = omp_get_wtime()
        if (merged /= 0) then
          call merge_kernel(img_size, img, img1, img2, tn, bn)
        else
          call find_moving_pixels(img_size, img, img1, img2, tn, mp)
          call update_background(img_size, img, mp, bn)
          call update_threshold(img_size, img, mp, bn, tn)
        end if
        time_ticks = time_ticks + int((omp_get_wtime() - start_time) * 1000000000.0_real64, int64)
        call merge_ref(img_size, img, img1, img2, tn_ref, bn_ref)
      end if
    end do

    if (repeat <= 2) then
      kernel_time_us = 0.0_real64
    else
      kernel_time_us = real(time_ticks, real64) * 1.0e-3_real64 / real(repeat - 2, real64)
    end if
    write(*,'(A,F0.6,A)') 'Average kernel execution time: ', kernel_time_us, ' (us)'
  !$omp end target data

  max_error = 0
  do j = 1, img_size
    max_error = max(max_error, abs(tn(j) - tn_ref(j)))
    max_error = max(max_error, abs(bn(j) - bn_ref(j)))
  end do
  write(*,'(A,I0)') 'Max error is ', max_error
  write(*,'(A)') merge('FAIL', 'PASS', max_error /= 0)

  deallocate(img, img1, img2, bn, bn_ref, mp, tn, tn_ref)

contains

  subroutine find_moving_pixels(n, img, img1, img2, tn, mp)
    integer, intent(in) :: n, img(:), img1(:), img2(:), tn(:)
    integer, intent(out) :: mp(:)
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, n
      if (abs(img(i) - img1(i)) > tn(i) .or. abs(img(i) - img2(i)) > tn(i)) then
        mp(i) = 255
      else
        mp(i) = 0
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine find_moving_pixels

  subroutine update_background(n, img, mp, bn)
    integer, intent(in) :: n, img(:), mp(:)
    integer, intent(inout) :: bn(:)
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, n
      if (mp(i) == 0) bn(i) = int(0.92_real32 * real(bn(i), real32) + 0.08_real32 * real(img(i), real32))
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_background

  subroutine update_threshold(n, img, mp, bn, tn)
    integer, intent(in) :: n, img(:), mp(:), bn(:)
    integer, intent(inout) :: tn(:)
    integer :: i
    real(real32) :: th

    !$omp target teams distribute parallel do thread_limit(block_size) private(th)
    do i = 1, n
      if (mp(i) == 0) then
        th = 0.92_real32 * real(tn(i), real32) + 0.24_real32 * real(img(i) - bn(i), real32)
        tn(i) = int(max(th, 20.0_real32))
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_threshold

  subroutine merge_kernel(n, img, img1, img2, tn, bn)
    integer, intent(in) :: n, img(:), img1(:), img2(:)
    integer, intent(inout) :: tn(:), bn(:)
    integer :: i
    real(real32) :: th

    !$omp target teams distribute parallel do thread_limit(block_size) private(th)
    do i = 1, n
      if (abs(img(i) - img1(i)) <= tn(i) .and. abs(img(i) - img2(i)) <= tn(i)) then
        bn(i) = int(0.92_real32 * real(bn(i), real32) + 0.08_real32 * real(img(i), real32))
        th = 0.92_real32 * real(tn(i), real32) + 0.24_real32 * real(img(i) - bn(i), real32)
        tn(i) = int(max(th, 20.0_real32))
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine merge_kernel

  subroutine merge_ref(n, img, img1, img2, tn, bn)
    integer, intent(in) :: n, img(:), img1(:), img2(:)
    integer, intent(inout) :: tn(:), bn(:)
    integer :: i
    real(real32) :: th

    do i = 1, n
      if (abs(img(i) - img1(i)) <= tn(i) .and. abs(img(i) - img2(i)) <= tn(i)) then
        bn(i) = int(0.92_real32 * real(bn(i), real32) + 0.08_real32 * real(img(i), real32))
        th = 0.92_real32 * real(tn(i), real32) + 0.24_real32 * real(img(i) - bn(i), real32)
        tn(i) = int(max(th, 20.0_real32))
      end if
    end do
  end subroutine merge_ref

  subroutine seed_rng(seed)
    integer(int64), intent(in) :: seed
    call c_seed_rng(int(seed, c_int64_t))
  end subroutine seed_rng

  integer function next_rand_byte() result(value)
    value = int(c_next_rand_byte())
  end function next_rand_byte

end program main
