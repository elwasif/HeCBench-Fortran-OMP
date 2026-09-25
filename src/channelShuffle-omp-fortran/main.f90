! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_threads = 256
  character(len=256) :: arg0, arg
  integer :: group_size, width, height, repeat
  integer :: n, channels, numel, i, data_size
  real(real32), allocatable :: x(:), y(:), y_ref(:)
  real(real64) :: elapsed
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <group size> <width> <height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) group_size
  call get_command_argument(2, arg); read(arg, *) width
  call get_command_argument(3, arg); read(arg, *) height
  call get_command_argument(4, arg); read(arg, *) repeat
  if (group_size <= 0 .or. width <= 0 .or. height <= 0 .or. repeat <= 0) stop 1

  n = 1
  do while (n <= 64)
    channels = 32
    do while (channels <= 512)
      write(*,*)
      write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '(N=', n, ' C=', channels, ' W=', width, ' H=', height, ')'

      numel = n * channels * width * height
      data_size = numel
      allocate(x(numel), y(numel), y_ref(numel))
      do i = 1, numel
        x(i) = real(i - 1, real32) / real(numel, real32)
      end do
      y = 0.0_real32
      y_ref = 0.0_real32

      !$omp target data map(to: x(1:numel)) map(alloc: y(1:numel))
      elapsed = channel_shuffle_nhwc(x, n, channels, group_size, data_size, y, repeat)
      !$omp target update from(y(1:numel))
      call channel_shuffle_nhwc_cpu(x, n, channels, group_size, data_size, y_ref)
      ok = all(y == y_ref)
      if (ok) then
        write(*,'(A,F0.6,A)') 'Average time of channel shuffle (NHWC): ', elapsed * 1.0e3_real64, ' (ms)'
      else
        write(*,'(A)') 'Failed to pass channel shuffle (NHWC) check'
      end if

      elapsed = channel_shuffle_nchw(x, n, channels, group_size, data_size, y, repeat)
      !$omp target update from(y(1:numel))
      call channel_shuffle_nchw_cpu(x, n, channels, group_size, data_size, y_ref)
      ok = all(y == y_ref)
      if (ok) then
        write(*,'(A,F0.6,A)') 'Average time of channel shuffle (NCHW): ', elapsed * 1.0e3_real64, ' (ms)'
      else
        write(*,'(A)') 'Failed to pass channel shuffle (NCHW) check'
      end if
      !$omp end target data

      deallocate(x, y, y_ref)
      channels = channels * 4
    end do
    n = n * 4
  end do

contains

  real(real64) function channel_shuffle_nhwc(x, n, channels, groups, numel, y, repeat) result(elapsed)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: y(:)
    integer, intent(in) :: n, channels, groups, numel, repeat
    integer :: iter, k, hxw, outer
    real(real64) :: start_time

    if (mod(channels, groups) /= 0 .or. numel < n * channels) error stop 'invalid channel shuffle shape'
    k = channels / groups
    hxw = numel / (n * channels)
    outer = n * hxw

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call channel_shuffle_nhwc_kernel(outer, groups, k, x, y)
    end do
    elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  end function channel_shuffle_nhwc

  real(real64) function channel_shuffle_nchw(x, n, channels, groups, numel, y, repeat) result(elapsed)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: y(:)
    integer, intent(in) :: n, channels, groups, numel, repeat
    integer :: iter, k, hxw
    real(real64) :: start_time

    if (mod(channels, groups) /= 0 .or. numel < n * channels) error stop 'invalid channel shuffle shape'
    k = channels / groups
    hxw = numel / (n * channels)

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call channel_shuffle_nchw_kernel(n, groups, k, hxw, x, y)
    end do
    elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  end function channel_shuffle_nchw

  subroutine channel_shuffle_nhwc_kernel(outer, groups, k, x, y)
    integer, intent(in) :: outer, groups, k
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: y(:)
    integer :: o, i, channels

    channels = groups * k
    !$omp target teams distribute parallel do collapse(2) num_threads(num_threads) &
    !$omp& map(to: x(1:outer*channels)) map(tofrom: y(1:outer*channels))
    do o = 0, outer - 1
      do i = 0, channels - 1
        y(o * channels + i + 1) = x(o * channels + modulo(i, groups) * k + i / groups + 1)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine channel_shuffle_nhwc_kernel

  subroutine channel_shuffle_nchw_kernel(n, groups, k, hxw, x, y)
    integer, intent(in) :: n, groups, k, hxw
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: y(:)
    integer :: batch, c, s, channels

    channels = groups * k
    !$omp target teams distribute parallel do collapse(3) num_threads(num_threads) &
    !$omp& map(to: x(1:n*channels*hxw)) map(tofrom: y(1:n*channels*hxw))
    do batch = 0, n - 1
      do c = 0, channels - 1
        do s = 0, hxw - 1
          y((batch * channels + c) * hxw + s + 1) = &
            x((batch * channels + modulo(c, groups) * k + c / groups) * hxw + s + 1)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine channel_shuffle_nchw_kernel

  subroutine channel_shuffle_nhwc_cpu(x, n, channels, groups, numel, y)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: y(:)
    integer, intent(in) :: n, channels, groups, numel
    integer :: k, hxw, outer, o, i

    k = channels / groups
    hxw = numel / (n * channels)
    outer = n * hxw
    do o = 0, outer - 1
      do i = 0, channels - 1
        y(o * channels + i + 1) = x(o * channels + modulo(i, groups) * k + i / groups + 1)
      end do
    end do
  end subroutine channel_shuffle_nhwc_cpu

  subroutine channel_shuffle_nchw_cpu(x, n, channels, groups, numel, y)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: y(:)
    integer, intent(in) :: n, channels, groups, numel
    integer :: k, hxw, batch, c, s

    k = channels / groups
    hxw = numel / (n * channels)
    do batch = 0, n - 1
      do c = 0, channels - 1
        do s = 0, hxw - 1
          y((batch * channels + c) * hxw + s + 1) = &
            x((batch * channels + modulo(c, groups) * k + c / groups) * hxw + s + 1)
        end do
      end do
    end do
  end subroutine channel_shuffle_nchw_cpu

end program main
