! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: num_elems, block_size, repeat, i
  integer :: nres(1)
  integer, allocatable :: input(:), output(:), h_output(:)
  integer :: h_flt_count
  real(real64) :: start_time, end_time, avg_ms
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of elements> <block size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) num_elems
  read(arg2, *) block_size
  read(arg3, *) repeat

  if (num_elems <= 0 .or. block_size <= 0 .or. repeat <= 0) stop 1

  allocate(input(num_elems), output(num_elems), h_output(num_elems))
  do i = 1, num_elems
    input(i) = (i - 1) - num_elems / 2
  end do
  output = 0
  nres = 0

  !$omp target data map(to: input(1:num_elems)) map(tofrom: nres(1:1)) map(from: output(1:num_elems))
  start_time = omp_get_wtime()
  do i = 1, repeat
    nres(1) = 0
    !$omp target update to(nres(1:1))
    call filter_positive(input, output, nres, num_elems, block_size)
  end do
  end_time = omp_get_wtime()
  avg_ms = (end_time - start_time) * 1.0e3_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time ', avg_ms, ' (ms)'
  !$omp end target data

  h_flt_count = 0
  do i = 1, num_elems
    if (input(i) > 0) then
      h_flt_count = h_flt_count + 1
      h_output(h_flt_count) = input(i)
    end if
  end do

  call sort_int(h_output, h_flt_count)
  call sort_int(output, nres(1))

  ok = (h_flt_count == nres(1))
  if (ok) then
    do i = 1, h_flt_count
      if (h_output(i) /= output(i)) then
        ok = .false.
        exit
      end if
    end do
  end if

  print '(A)'
  if (ok) then
    print '(A)', 'Filter using shared memory PASS '
  else
    print '(A)', 'Filter using shared memory FAIL '
    stop 1
  end if

  deallocate(input, output, h_output)

contains

  subroutine filter_positive(input, output, nres, num_elems, block_size)
    integer, intent(in) :: input(:), num_elems, block_size
    integer, intent(out) :: output(:)
    integer, intent(inout) :: nres(:)
    integer :: idx, pos, value, teams, l_n, old

    teams = (num_elems + block_size - 1) / block_size
    !$omp target teams num_teams(teams) thread_limit(block_size) private(l_n)
    !$omp parallel private(idx, pos, value, old)
      idx = omp_get_team_num() * omp_get_num_threads() + omp_get_thread_num() + 1
      if (omp_get_thread_num() == 0) l_n = 0
      !$omp barrier

      if (idx <= num_elems) then
        value = input(idx)
        if (value > 0) then
          !$omp atomic capture
          pos = l_n
          l_n = l_n + 1
          !$omp end atomic
        end if
      end if
      !$omp barrier

      if (omp_get_thread_num() == 0) then
        !$omp atomic capture
        old = nres(1)
        nres(1) = nres(1) + l_n
        !$omp end atomic
        l_n = old
      end if
      !$omp barrier

      if (idx <= num_elems) then
        if (value > 0) output(pos + l_n + 1) = value
      end if
      !$omp barrier
    !$omp end parallel
    !$omp end target teams
  end subroutine filter_positive

  subroutine sort_int(a, n)
    integer, intent(inout) :: a(:)
    integer, intent(in) :: n
    integer :: i, j, key

    do i = 2, n
      key = a(i)
      j = i - 1
      do while (j >= 1 .and. a(j) > key)
        a(j + 1) = a(j)
        j = j - 1
      end do
      a(j + 1) = key
    end do
  end subroutine sort_int

end program main
