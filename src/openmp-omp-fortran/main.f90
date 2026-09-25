! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg
  integer :: repeat, num_gpus, b, pass, f, num_threads
  integer(int64), parameter :: nwords_per_gpu = 33554432_int64
  integer(int64) :: nwords
  integer(int32), allocatable :: a(:)
  real(real64) :: start_time, end_time, overhead
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() < 1) stop 1
  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat < 0) stop 1

  write(*,'(A,A)') trim(arg0), ' Starting...'
  write(*,*)

  num_gpus = 1
  write(*,'(A,I0)') 'number of host CPUs:' // achar(9), omp_get_num_procs()
  write(*,'(A,I0)') 'number of devices:' // achar(9), num_gpus

  nwords = int(num_gpus, int64) * nwords_per_gpu
  b = 3
  allocate(a(0:nwords - 1))
  overhead = 0.0_real64

  do pass = 0, 1
    f = 1
    do
      if (f > 32) exit
      num_threads = f * num_gpus
      start_time = omp_get_wtime()
      call omp_set_num_threads(num_threads)
      !$omp parallel shared(a, nwords, repeat, b)
      call run_thread_partition(a, nwords, repeat, b)
      !$omp end parallel
      end_time = omp_get_wtime()

      write(*,'(A,F8.6,A,I0,A)') 'Work took ', end_time - start_time, &
        ' seconds with ', num_threads, ' CPU threads'

      if (f == 1) then
        if (pass == 0) then
          overhead = end_time - start_time
        else
          overhead = overhead - (end_time - start_time)
        end if
      end if

      ok = correct_result(a, nwords, b, repeat)
      if (ok) then
        write(*,'(A)') 'PASS'
      else
        write(*,'(A)') 'FAIL'
      end if

      f = f * 2
    end do
  end do

  write(*,'(A,F8.6,A)') 'Runtime overhead of first run is ', overhead, ' seconds'
  deallocate(a)

contains

  subroutine run_thread_partition(a, nwords, repeat, b)
    integer(int32), intent(inout), target :: a(0:)
    integer(int64), intent(in) :: nwords
    integer, intent(in) :: repeat, b
    integer :: cpu_thread_id, num_cpu_threads
    integer :: j
    integer(int64) :: nwords_per_kernel, first, n
    integer(int32), pointer :: sub_a(:)

    cpu_thread_id = omp_get_thread_num()
    num_cpu_threads = omp_get_num_threads()
    nwords_per_kernel = nwords / int(num_cpu_threads, int64)
    first = int(cpu_thread_id, int64) * nwords_per_kernel
    sub_a(0:nwords_per_kernel - 1_int64) => a(first:first + nwords_per_kernel - 1_int64)

    do n = 0_int64, nwords_per_kernel - 1_int64
      sub_a(n) = int(first + n, int32)
    end do

    !$omp target data map(tofrom: sub_a(0:nwords_per_kernel - 1_int64))
    !$omp target teams distribute parallel do thread_limit(256)
    do n = 0_int64, nwords_per_kernel - 1_int64
      do j = 0, repeat - 1
        sub_a(n) = sub_a(n) + mod(j, b)
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp end target data
  end subroutine run_thread_partition

  logical function correct_result(data, nwords, b, repeat)
    integer(int32), intent(in) :: data(0:)
    integer(int64), intent(in) :: nwords
    integer, intent(in) :: b, repeat
    integer :: j, sum
    integer(int64) :: i

    sum = 0
    do j = 0, repeat - 1
      sum = sum + mod(j, b)
    end do

    correct_result = .true.
    do i = 0_int64, nwords - 1_int64
      if (data(i) /= int(i, int32) + sum) then
        write(*,'(A,I0,A,I0)') 'check: ', data(i), ' != ', sum
        correct_result = .false.
        return
      end if
    end do
  end function correct_result

end program main
