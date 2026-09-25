! SPDX-License-Identifier: CC0-1.0
program threadfence_main
  use, intrinsic :: iso_fortran_env, only: real32
  use omp_lib
  implicit none

  integer :: repeat, n, blocks, grids, iter
  integer, allocatable :: h_count(:)
  real(real32), allocatable :: h_array(:), h_result(:)
  real(8) :: elapsed, start_time, end_time
  logical :: ok

  call parse_args(repeat, n)

  blocks = 256
  grids = (n + blocks - 1) / blocks

  allocate(h_array(n))
  allocate(h_result(grids))
  allocate(h_count(1))
  h_array = -1.0_real32
  h_count(1) = 0
  elapsed = 0.0_8
  ok = .true.

  !$omp target data map(to: h_array) map(tofrom: h_count) map(alloc: h_result)
  do iter = 1, repeat
    start_time = omp_get_wtime()
    call sum_device(grids, blocks, h_array, n, h_count, h_result)
    end_time = omp_get_wtime()
    elapsed = elapsed + (end_time - start_time)

    !$omp target update from(h_result(1:1))

    if (h_result(1) /= -1.0_real32 * real(n, real32)) then
      ok = .false.
      exit
    end if
  end do
  !$omp end target data

  if (ok) then
    write(*,'(A,F0.6,A)') "Average kernel execution time: ", &
         real(elapsed * 1.0d3 / dble(repeat), real32), " (ms)"
  end if

  write(*,'(A)') merge("PASS", "FAIL", ok)

contains

  subroutine parse_args(repeat, n)
    integer, intent(out) :: repeat, n
    character(len=64) :: arg
    integer :: status

    if (command_argument_count() /= 2) then
      call get_command_argument(0, arg)
      write(*,'(A,A,A)') "Usage: ", trim(arg), " <repeat> <array length>"
      stop 1
    end if

    call get_command_argument(1, arg)
    read(arg, *, iostat=status) repeat
    if (status /= 0) repeat = 0

    call get_command_argument(2, arg)
    read(arg, *, iostat=status) n
    if (status /= 0) n = 0
  end subroutine parse_args

  subroutine sum_device(teams, blocks, array, n, count, result)
    integer, intent(in) :: teams, blocks
    integer, intent(in) :: n
    integer, intent(inout) :: count(1)
    real(real32), intent(in) :: array(n)
    real(real32), intent(inout) :: result(teams)
    integer :: bid, num_blocks, block_size, lid, gid, i, value
    logical :: isLastBlockDone
    real(real32) :: partialSum

    !$omp target teams num_teams(teams) thread_limit(blocks)
    !$omp parallel private(bid, num_blocks, block_size, lid, gid, i, value) &
    !$omp& shared(array, count, result, partialSum, isLastBlockDone)
    bid = omp_get_team_num()
    num_blocks = teams
    block_size = blocks
    lid = omp_get_thread_num()
    gid = bid * block_size + lid

    if (lid == 0) partialSum = 0.0_real32
    !$omp barrier

    if (gid < n) then
      !$omp atomic update
      partialSum = partialSum + array(gid + 1)
    end if

    !$omp barrier

    if (lid == 0) then
      result(bid + 1) = partialSum

      !$omp atomic capture
      value = count(1)
      count(1) = count(1) + 1
      !$omp end atomic

      isLastBlockDone = (value == (num_blocks - 1))
    end if

    !$omp barrier

    if (isLastBlockDone) then
      if (lid == 0) partialSum = 0.0_real32
      !$omp barrier

      do i = lid, num_blocks - 1, block_size
        !$omp atomic update
        partialSum = partialSum + result(i + 1)
      end do

      !$omp barrier

      if (lid == 0) then
        result(1) = partialSum
        count(1) = 0
      end if
    end if
    !$omp end parallel
    !$omp end target teams
  end subroutine sum_device

end program threadfence_main
