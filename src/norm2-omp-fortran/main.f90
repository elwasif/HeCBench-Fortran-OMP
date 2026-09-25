! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: first_n = 512 * 1024
  integer(int32), parameter :: last_n = 512 * 1024 * 1024
  character(len=256) :: arg0, arg1
  integer :: argc
  integer(int32) :: repeat
  logical :: ok
  real(real32), allocatable :: h_result(:)

  argc = command_argument_count()
  call get_command_argument(0, arg0)
  if (argc /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  repeat = max(1_int32, repeat)

  allocate(h_result(repeat))
  ok = .true.

  call run_benchmark(repeat, h_result, ok)

  deallocate(h_result)

  if (ok) print '(A)', 'PASS'

contains

  subroutine run_benchmark(repeat, h_result, ok)
    integer(int32), intent(in) :: repeat
    real(real32), intent(out) :: h_result(repeat)
    logical, intent(inout) :: ok
    integer :: i, j
    integer(int32) :: n
    character(len=16) :: elements_label
    real(real32), allocatable :: a(:)
    real(real64) :: gold, sum, kstart, kend, avg_us, gops
    real(real32) :: t

    n = first_n
    do while (n <= last_n)
      allocate(a(n))

      gold = 0.0_real64
      do i = 1, n
        a(i) = real(mod(i, 7), kind=real32)
        gold = gold + real(a(i), kind=real64) * real(a(i), kind=real64)
      end do
      gold = sqrt(gold)

      !$omp target data map(to: a(1:n))
      kstart = omp_get_wtime()

      do j = 1, repeat
        sum = 0.0_real64
        !$omp target teams distribute parallel do thread_limit(256) reduction(+:sum)
        do i = 1, n
          t = a(i) * a(i)
          sum = sum + real(t, kind=real64)
        end do
        !$omp end target teams distribute parallel do
        h_result(j) = real(sqrt(sum), kind=real32)
      end do

      kend = omp_get_wtime()
      !$omp end target data

      avg_us = ((kend - kstart) * 1.0d6) / real(repeat, kind=real64)
      gops = real((2_int64 * int(n, kind=int64) + 1_int64) * int(repeat, int64), kind=real64) / &
        ((kend - kstart) * 1.0d9)
      write(elements_label, '(F6.2)') real(n, kind=real64) / (1024.0_real64 * 1024.0_real64)
      elements_label = adjustl(elements_label)
      print '("#elements = ",A," M: average omp nrm2 execution time = ",F0.6," (us), performance = ",F0.6," (Gop/s)")', &
        trim(elements_label), avg_us, gops

      do j = 1, repeat
        if (abs(real(gold, kind=real32) - h_result(j)) > 1.0e-3_real32) then
          print '("FAIL at iteration ",I0,": gold=",F0.6," actual=",F0.6," for ",I0," elements")', &
            j - 1, real(gold, kind=real32), h_result(j), n
          ok = .false.
          exit
        end if
      end do

      deallocate(a)
      if (.not. ok) exit
      n = n * 2
    end do
  end subroutine run_benchmark

end program main
