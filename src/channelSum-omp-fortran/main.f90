! SPDX-License-Identifier: CC0-1.0
program channel_sum_main
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: num_threads = 256
  integer :: w, h, repeat, n, c, numel, hxw
  integer :: i
  integer, allocatable :: x(:), sum_out(:), sumsq(:), ref_sum(:)
  real(real64) :: elapsed
  character(len=64) :: arg
  logical :: ok

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <width> <height> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) w
  call get_command_argument(2, arg)
  read(arg, *) h
  call get_command_argument(3, arg)
  read(arg, *) repeat

  hxw = w * h

  n = 1
  do while (n <= 64)
    c = 32
    do while (c <= 512)
      write(*,*)
      write(*,'("(N=",I0," C=",I0," W=",I0," H=",I0,")")') n, c, w, h

      numel = n * c * hxw
      allocate(x(0:numel-1), sum_out(0:c-1), sumsq(0:c-1), ref_sum(0:c-1))

      call c_srand(int(numel, c_int))
      do i = 0, numel - 1
        x(i) = mod(c_rand(), 256_c_int)
      end do

      !$omp target data map(to: x) map(from: sum_out, sumsq)
      call compute_channel_sum_nhwc(n, c, hxw, x, sum_out, sumsq, elapsed, repeat)
      !$omp target update from(sum_out)
      call ref_nhwc(n, c, hxw, x, ref_sum, sumsq)
      ok = check_values(c, sum_out, ref_sum)
      write(*,'("Average time of channel sum (nhwc): ",F0.6," (ms)")') elapsed * 1.0e3_real64
      write(*,'("Verification ",A," for channel sum (nhwc)")') merge("PASS", "FAIL", ok)

      call compute_channel_sum_nchw(n, c, hxw, x, sum_out, sumsq, elapsed, repeat)
      !$omp target update from(sum_out)
      !$omp end target data
      call ref_nchw(n, c, hxw, x, ref_sum, sumsq)
      ok = check_values(c, sum_out, ref_sum)
      write(*,'("Average time of channel sum (nchw): ",F0.6," (ms)")') elapsed * 1.0e3_real64
      write(*,'("Verification ",A," for channel sum (nchw)")') merge("PASS", "FAIL", ok)

      deallocate(x, sum_out, sumsq, ref_sum)
      c = c * 4
    end do
    n = n * 4
  end do

contains

  subroutine compute_channel_sum_nhwc(n, c_count, hxw_count, x, sum_out, sumsq, elapsed, repeat)
    integer, intent(in) :: n, c_count, hxw_count, repeat
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    real(real64), intent(out) :: elapsed
    integer :: iter
    real(real64) :: start_time

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call channel_sum_nhwc(n, c_count, hxw_count, x, sum_out, sumsq)
    end do
    elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  end subroutine compute_channel_sum_nhwc

  subroutine compute_channel_sum_nchw(n, c_count, hxw_count, x, sum_out, sumsq, elapsed, repeat)
    integer, intent(in) :: n, c_count, hxw_count, repeat
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    real(real64), intent(out) :: elapsed
    integer :: iter
    real(real64) :: start_time

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call channel_sum_nchw(n, c_count, hxw_count, x, sum_out, sumsq)
    end do
    elapsed = (omp_get_wtime() - start_time) / real(repeat, real64)
  end subroutine compute_channel_sum_nchw

  subroutine channel_sum_nhwc(n, c_count, hxw_count, x, sum_out, sumsq)
    integer, intent(in) :: n, c_count, hxw_count
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    integer :: chan, idx, input_index, m_val, v_val

    !$omp target teams distribute private(idx, input_index, m_val, v_val) num_teams(c_count)
    do chan = 0, c_count - 1
      m_val = 0
      v_val = 0
      !$omp parallel do reduction(+:m_val, v_val) num_threads(num_threads)
      do idx = 0, n * hxw_count - 1
        input_index = idx * c_count + chan
        m_val = m_val + x(input_index)
        v_val = v_val + x(input_index) * x(input_index)
      end do
      !$omp end parallel do
      sum_out(chan) = m_val
      sumsq(chan) = v_val
    end do
    !$omp end target teams distribute
  end subroutine channel_sum_nhwc

  subroutine channel_sum_nchw(n, c_count, hxw_count, x, sum_out, sumsq)
    integer, intent(in) :: n, c_count, hxw_count
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    integer :: chan, sample, hw, input_index, m_val, v_val

    !$omp target teams distribute private(sample, hw, input_index, m_val, v_val) num_teams(c_count)
    do chan = 0, c_count - 1
      m_val = 0
      v_val = 0
      !$omp parallel do collapse(2) reduction(+:m_val, v_val) num_threads(num_threads)
      do sample = 0, n - 1
        do hw = 0, hxw_count - 1
          input_index = (sample * c_count + chan) * hxw_count + hw
          m_val = m_val + x(input_index)
          v_val = v_val + x(input_index) * x(input_index)
        end do
      end do
      !$omp end parallel do
      sum_out(chan) = m_val
      sumsq(chan) = v_val
    end do
    !$omp end target teams distribute
  end subroutine channel_sum_nchw

  subroutine ref_nhwc(n, c_count, hxw_count, x, sum_out, sumsq)
    integer, intent(in) :: n, c_count, hxw_count
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    integer :: chan, idx, input_index, m_val, v_val

    do chan = 0, c_count - 1
      m_val = 0
      v_val = 0
      do idx = 0, n * hxw_count - 1
        input_index = idx * c_count + chan
        m_val = m_val + x(input_index)
        v_val = v_val + x(input_index) * x(input_index)
      end do
      sum_out(chan) = m_val
      sumsq(chan) = v_val
    end do
  end subroutine ref_nhwc

  subroutine ref_nchw(n, c_count, hxw_count, x, sum_out, sumsq)
    integer, intent(in) :: n, c_count, hxw_count
    integer, intent(in) :: x(0:)
    integer, intent(out) :: sum_out(0:), sumsq(0:)
    integer :: chan, sample, hw, input_index, m_val, v_val

    do chan = 0, c_count - 1
      m_val = 0
      v_val = 0
      do sample = 0, n - 1
        do hw = 0, hxw_count - 1
          input_index = (sample * c_count + chan) * hxw_count + hw
          m_val = m_val + x(input_index)
          v_val = v_val + x(input_index) * x(input_index)
        end do
      end do
      sum_out(chan) = m_val
      sumsq(chan) = v_val
    end do
  end subroutine ref_nchw

  logical function check_values(size, device_values, reference_values)
    integer, intent(in) :: size
    integer, intent(in) :: device_values(0:), reference_values(0:)
    integer :: idx

    check_values = .true.
    do idx = 0, size - 1
      if (abs(device_values(idx) - reference_values(idx)) > 1) then
        check_values = .false.
        exit
      end if
    end do
  end function check_values

end program channel_sum_main
