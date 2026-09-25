! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: len = 256_int32
  character(len=256) :: arg0, arg1
  integer(int32) :: iteration
  integer(int32) :: test(len), scratch(len), gold_even(len), gold_odd(len)
  integer(int32) :: i, iter, count
  integer(int64) :: total_count
  logical :: error
  real(real64) :: start_time, end_time, total_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ./', trim(arg0), ' <iterations>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) iteration
  if (iteration <= 0_int32) stop 1

  do i = 1, len
    gold_even(i) = i - 1_int32
    gold_odd(i) = len - i
  end do

  error = .false.
  total_time = 0.0_real64
  total_count = 0_int64

  !$omp target data map(alloc: test(1:len))
  do iter = 1, iteration
    count = reverse_count(iter)
    total_count = total_count + int(count, int64)
    test = gold_even
    !$omp target update to(test(1:len))

    start_time = omp_get_wtime()
    call reverse_repeated(test, count)
    end_time = omp_get_wtime()
    total_time = total_time + end_time - start_time

    !$omp target update from(test(1:len))
    if (mod(count, 2_int32) == 0_int32) then
      error = any(test /= gold_even)
    else
      error = any(test /= gold_odd)
    end if
    if (error) exit
  end do
  !$omp end target data

  print '(A,F8.6,A)', 'Total kernel execution time: ', total_time, ' (s)'
  if (error) then
    print '(A)', 'FAIL'
    stop 1
  end if
  print '(A)', 'PASS'

contains

  integer(int32) function reverse_count(iter) result(count)
    integer(int32), intent(in) :: iter

    count = 100_int32 + modulo(37_int32 * (iter - 1_int32), 9900_int32)
  end function reverse_count

  subroutine reverse_repeated(test, count)
    integer(int32), intent(inout) :: test(:)
    integer(int32), intent(in) :: count
    integer(int32) :: pass, t

    do pass = 1, count
      !$omp target teams num_teams(1) thread_limit(len)
      block
        integer(int32) :: s(len)

        !$omp parallel private(t) shared(s, test)
        t = omp_get_thread_num() + 1_int32
        s(t) = test(t)
        !$omp barrier
        test(t) = s(len - t + 1_int32)
        !$omp end parallel
      end block
      !$omp end target teams
    end do
  end subroutine reverse_repeated

end program main
