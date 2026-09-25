! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_threads = 128
  integer(int64), parameter :: thread_work_size = 4_int64
  integer(int64), parameter :: block_work_size = thread_work_size * int(num_threads, int64)
  real(real32), parameter :: tolerance = 1.0e-3_real32

  character(len=256) :: arg0, arg1
  integer :: repeat
  integer(int64) :: numel, num_teams, i
  real(real32), allocatable :: x1(:), x2(:), cos_values(:), sin_values(:), o1(:), o2(:)
  real(real64) :: start_time, end_time, avg_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0) stop 1

  numel = block_work_size * 10000_int64
  num_teams = (numel + block_work_size - 1_int64) / block_work_size

  print '(A,I0)', 'Number of elements: ', numel
  print '(A,I0,A,I0)', 'Number of teams (blocks): ', num_teams, ', threads per team: ', num_threads

  allocate(x1(numel), x2(numel), cos_values(numel), sin_values(numel), o1(numel), o2(numel))

  do i = 1_int64, numel
    x1(i) = real(i, real32) / real(numel, real32)
    x2(i) = real(i, real32) / real(numel, real32)
    cos_values(i) = cos(real(i - 1_int64, real32) / (10000.0_real32 ** &
        (real(i - 1_int64, real32) / real(numel, real32))))
    sin_values(i) = sin(real(i - 1_int64, real32) / (10000.0_real32 ** &
        (real(i - 1_int64, real32) / real(numel, real32))))
  end do

  o1 = 0.0_real32
  o2 = 0.0_real32

  !$omp target data map(to: x1(1:numel), x2(1:numel), cos_values(1:numel), sin_values(1:numel)) &
  !$omp& map(from: o1(1:numel), o2(1:numel))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call launch_kernel(numel, num_teams, o1, o2, x1, x2, cos_values, sin_values)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time: ', avg_us, ' (us)'

  ok = verify_outputs(numel, x1, x2, cos_values, sin_values, o1, o2)
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(x1, x2, cos_values, sin_values, o1, o2)

contains

  subroutine launch_kernel(n, teams, out1, out2, in1, in2, cvals, svals)
    integer(int64), intent(in) :: n, teams
    real(real32), intent(out) :: out1(:), out2(:)
    real(real32), intent(in) :: in1(:), in2(:), cvals(:), svals(:)
    integer(int64) :: idx

    !$omp target teams distribute parallel do num_teams(teams) num_threads(num_threads)
    do idx = 1_int64, n
      out1(idx) = in1(idx) * cvals(idx) - in2(idx) * svals(idx)
      out2(idx) = in1(idx) * svals(idx) + in2(idx) * cvals(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine launch_kernel

  function verify_outputs(n, in1, in2, cvals, svals, out1, out2) result(ok)
    integer(int64), intent(in) :: n
    real(real32), intent(in) :: in1(:), in2(:), cvals(:), svals(:), out1(:), out2(:)
    logical :: ok
    integer(int64) :: idx
    real(real32) :: expected1, expected2

    ok = .true.
    do idx = 1_int64, n
      expected1 = in1(idx) * cvals(idx) - in2(idx) * svals(idx)
      expected2 = in1(idx) * svals(idx) + in2(idx) * cvals(idx)
      if (abs(expected1 - out1(idx)) > tolerance .or. abs(expected2 - out2(idx)) > tolerance) then
        ok = .false.
        return
      end if
    end do
  end function verify_outputs

end program main
