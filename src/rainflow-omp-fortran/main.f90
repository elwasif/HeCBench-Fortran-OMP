! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer :: num_history, repeat, n, total_length, i
  integer(int32), allocatable :: history_lengths(:), result_lengths(:), ref_result_lengths(:), points(:)
  real(real64), allocatable :: history(:), extrema(:), results(:, :)
  real(real64) :: start_time, end_time, avg_us
  logical :: ok
  integer(c_int), parameter :: c_rand_max = 2147483647_c_int

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

  if (command_argument_count() /= 2) stop 1
  num_history = read_arg(1)
  repeat = read_arg(2)
  if (num_history <= 0 .or. repeat <= 0) stop 1

  allocate(history_lengths(num_history + 1), result_lengths(num_history), ref_result_lengths(num_history))
  call c_srand(123_c_int)
  total_length = 0
  do n = 1, num_history
    history_lengths(n) = int(total_length, int32)
    total_length = total_length + (mod(c_rand(), 10_c_int) + 1_c_int) * 100
  end do
  history_lengths(num_history + 1) = int(total_length, int32)

  print '(A,I0)', 'Total history length = ', total_length

  allocate(history(total_length), extrema(total_length), points(total_length), results(3, total_length))
  do i = 1, total_length
    history(i) = real(c_rand(), real64) / real(c_rand_max, real64)
  end do

  result_lengths = 0_int32
  ref_result_lengths = 0_int32

  !$omp target data map(to: history_lengths(1:num_history+1), history(1:total_length)) &
  !$omp& map(alloc: extrema(1:total_length), points(1:total_length), results(1:3,1:total_length)) &
  !$omp& map(from: result_lengths(1:num_history))
  start_time = omp_get_wtime()
  do n = 1, repeat
    call rainflow_device(history, history_lengths, extrema, points, results, result_lengths, num_history)
  end do
  end_time = omp_get_wtime()
  avg_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time: ', avg_us, ' (us)'
  !$omp end target data

  call rainflow_reference(history, history_lengths, extrema, points, results, ref_result_lengths, num_history)

  ok = .true.
  do i = 1, num_history
    if (result_lengths(i) /= ref_result_lengths(i)) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(history_lengths, result_lengths, ref_result_lengths, history, extrema, points, results)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine rainflow_device(history, history_lengths, extrema, points, results, result_lengths, num_history)
    real(real64), intent(in) :: history(:)
    integer(int32), intent(in) :: history_lengths(:)
    real(real64), intent(out) :: extrema(:)
    real(real64), intent(out) :: results(:, :)
    integer(int32), intent(out) :: points(:), result_lengths(:)
    integer, intent(in) :: num_history
    integer :: i, offset, history_length

    !$omp target teams distribute parallel do thread_limit(256) private(offset, history_length)
    do i = 1, num_history
      offset = int(history_lengths(i)) + 1
      history_length = int(history_lengths(i + 1) - history_lengths(i))
      call execute_history(history, offset, history_length, extrema, points, results, result_lengths(i))
    end do
    !$omp end target teams distribute parallel do
  end subroutine rainflow_device

  subroutine rainflow_reference(history, history_lengths, extrema, points, results, result_lengths, num_history)
    real(real64), intent(in) :: history(:)
    integer(int32), intent(in) :: history_lengths(:)
    real(real64), intent(out) :: extrema(:)
    real(real64), intent(out) :: results(:, :)
    integer(int32), intent(out) :: points(:), result_lengths(:)
    integer, intent(in) :: num_history
    integer :: i, offset, history_length

    do i = 1, num_history
      offset = int(history_lengths(i)) + 1
      history_length = int(history_lengths(i + 1) - history_lengths(i))
      call execute_history(history, offset, history_length, extrema, points, results, result_lengths(i))
    end do
  end subroutine rainflow_reference

  subroutine execute_history(history, offset, history_length, extrema, points, results, result_length)
    real(real64), intent(in) :: history(:)
    integer, intent(in) :: offset, history_length
    real(real64), intent(out) :: extrema(:)
    real(real64), intent(out) :: results(:, :)
    integer(int32), intent(out) :: points(:)
    integer(int32), intent(out) :: result_length
    integer :: extrema_length, pidx, eidx, ridx, i
    real(real64) :: x_range, y_range, y_mean, range, mean

    call extrema_history(history, offset, history_length, extrema, extrema_length)

    pidx = -1
    eidx = -1
    ridx = -1
    do i = 0, extrema_length - 1
      pidx = pidx + 1
      eidx = eidx + 1
      points(offset + pidx) = int(eidx, int32)
      do while (pidx >= 2)
        x_range = abs(extrema(offset + int(points(offset + pidx - 1))) - extrema(offset + int(points(offset + pidx))))
        y_range = abs(extrema(offset + int(points(offset + pidx - 2))) - extrema(offset + int(points(offset + pidx - 1))))
        if (x_range < y_range) exit
        ridx = ridx + 1
        y_mean = 0.5_real64 * (extrema(offset + int(points(offset + pidx - 2))) + &
                               extrema(offset + int(points(offset + pidx - 1))))
        if (pidx == 2) then
          results(1, offset + ridx) = 0.5_real64
          results(2, offset + ridx) = y_range
          results(3, offset + ridx) = y_mean
          points(offset) = points(offset + 1)
          points(offset + 1) = points(offset + 2)
          pidx = 1
        else
          results(1, offset + ridx) = 1.0_real64
          results(2, offset + ridx) = y_range
          results(3, offset + ridx) = y_mean
          points(offset + pidx - 2) = points(offset + pidx)
          pidx = pidx - 2
        end if
      end do
    end do

    do i = 0, pidx - 1
      range = abs(extrema(offset + int(points(offset + i))) - extrema(offset + int(points(offset + i + 1))))
      mean = 0.5_real64 * (extrema(offset + int(points(offset + i))) + extrema(offset + int(points(offset + i + 1))))
      ridx = ridx + 1
      results(1, offset + ridx) = 0.5_real64
      results(2, offset + ridx) = range
      results(3, offset + ridx) = mean
    end do
    result_length = int(ridx + 1, int32)
  end subroutine execute_history

  subroutine extrema_history(history, offset, history_length, extrema, extrema_length)
    real(real64), intent(in) :: history(:)
    integer, intent(in) :: offset, history_length
    real(real64), intent(out) :: extrema(:)
    integer, intent(out) :: extrema_length
    integer :: i, eidx

    extrema(offset) = history(offset)
    eidx = 0
    do i = 1, history_length - 2
      if ((history(offset + i) > extrema(offset + eidx) .and. history(offset + i) > history(offset + i + 1)) .or. &
          (history(offset + i) < extrema(offset + eidx) .and. history(offset + i) < history(offset + i + 1))) then
        eidx = eidx + 1
        extrema(offset + eidx) = history(offset + i)
      end if
    end do
    eidx = eidx + 1
    extrema(offset + eidx) = history(offset + history_length - 1)
    extrema_length = eidx + 1
  end subroutine extrema_history

end program main
