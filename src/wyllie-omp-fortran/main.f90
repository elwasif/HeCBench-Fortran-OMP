! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: nil = -1_int32
  integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)

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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer(int32) :: elems, set_random_list, repeat, iter, i
  integer(int32), allocatable :: next(:), rank(:)
  integer(int64), allocatable :: list(:), original_list(:), d_res(:), h_res(:)
  real(real64) :: time_ms
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./' // trim(arg0) // ' <list size> <0 or 1> <repeat>'
    print '(A)', '0 and 1 indicate an ordered list and a random list, respectively'
    stop 255
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) elems
  read(arg2, *) set_random_list
  read(arg3, *) repeat

  if (elems <= 1_int32 .or. repeat <= 0_int32) then
    print '(A)', 'FAIL'
    stop 1
  end if

  allocate(next(0:elems - 1), rank(0:elems - 1))
  allocate(list(0:elems - 1), original_list(0:elems - 1))
  allocate(d_res(0:elems - 1), h_res(0:elems - 1))

  if (set_random_list /= 0_int32) then
    call random_list(next)
  else
    call ordered_list(next)
  end if

  do i = 0, elems - 1
    if (next(i) == nil) then
      rank(i) = 0_int32
    else
      rank(i) = 1_int32
    end if
    list(i) = pack_node(next(i), rank(i))
  end do
  original_list = list

  call run_device(list, original_list, elems, repeat, time_ms)

  do i = 0, elems - 1
    d_res(i) = iand(list(i), mask32)
  end do

  call build_reference(next, h_res)
  ok = all(h_res == d_res)

  print '(A,F0.6,A)', 'Average kernel execution time: ', time_ms, ' (ms)'
  print '(A)', merge('PASS', 'FAIL', ok)

contains

  subroutine ordered_list(next)
    integer(int32), intent(out) :: next(0:)
    integer(int32) :: i, n

    n = ubound(next, 1)
    do i = 0, n - 1
      next(i) = i + 1
    end do
    next(n) = nil
  end subroutine ordered_list

  subroutine random_list(next)
    integer(int32), intent(out) :: next(0:)
    logical, allocatable :: free_list(:)
    integer(int32) :: n_elems, nil_pos, prev_pos, pos, remaining

    n_elems = size(next)
    allocate(free_list(0:n_elems - 1))
    free_list = .true.

    call c_srand(123_c_int)
    nil_pos = modulo(int(c_rand(), int32), n_elems - 1_int32) + 1_int32
    free_list(nil_pos) = .false.
    next(nil_pos) = nil
    prev_pos = nil_pos
    remaining = n_elems

    do
      remaining = remaining - 1_int32
      if (remaining <= 1_int32) exit
      do
        pos = modulo(int(c_rand(), int32), size(next) - 1) + 1_int32
        if (free_list(pos)) exit
      end do
      free_list(pos) = .false.
      next(pos) = prev_pos
      prev_pos = pos
    end do

    next(0) = prev_pos
    deallocate(free_list)
  end subroutine random_list

  subroutine run_device(list, original_list, elems, repeat, time_ms)
    integer(int64), intent(inout) :: list(0:)
    integer(int64), intent(in) :: original_list(0:)
    integer(int32), intent(in) :: elems, repeat
    real(real64), intent(out) :: time_ms
    integer(int32) :: iter, teams
    real(real64) :: time_accum, start_time, end_time

    teams = (elems + 255_int32) / 256_int32
    time_accum = 0.0_real64

    !$omp target data map(tofrom: list(0:elems - 1))
    do iter = 0, repeat
      list = original_list
      !$omp target update to(list(0:elems - 1))

      start_time = omp_get_wtime()
      !$omp target teams num_teams(teams) thread_limit(256) firstprivate(elems)
      !$omp parallel
      block
        integer(int32) :: index, next_index, next_next
        integer(int64) :: node, next_node, temp

        index = omp_get_team_num() * 256 + omp_get_thread_num()
        if (index < elems) then
          do
            node = list(index)
            next_index = unpack_next(node)
            if (next_index == nil) exit

            next_node = list(next_index)
            next_next = unpack_next(next_node)
            if (next_next == nil) exit

            temp = iand(node, mask32)
            temp = temp + iand(next_node, mask32)
            temp = temp + ishft(int(next_next, int64), 32)
            !$omp barrier
            list(index) = temp
          end do
        end if
      end block
      !$omp end parallel
      !$omp end target teams
      end_time = omp_get_wtime()

      if (iter > 0) then
        time_accum = time_accum + (end_time - start_time)
      end if
    end do
    !$omp end target data

    time_ms = (time_accum * 1.0d3) / real(repeat, real64)
  end subroutine run_device

  subroutine build_reference(next, ranks)
    integer(int32), intent(in) :: next(0:)
    integer(int64), intent(out) :: ranks(0:)
    integer(int32) :: i, r

    ranks(0) = int(size(next) - 1, int64)
    i = 0_int32
    do r = 1, size(next) - 1
      ranks(next(i)) = int(size(next) - 1 - r, int64)
      i = next(i)
    end do
  end subroutine build_reference

  integer(int64) function pack_node(next_index, rank_value) result(value)
    integer(int32), intent(in) :: next_index, rank_value

    value = ishft(int(next_index, int64), 32) + iand(int(rank_value, int64), mask32)
  end function pack_node

  integer(int32) function unpack_next(node) result(next_index)
    integer(int64), intent(in) :: node

    next_index = int(ishft(node, -32), int32)
  end function unpack_next

end program main
