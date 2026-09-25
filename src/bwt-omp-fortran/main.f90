! SPDX-License-Identifier: CC0-1.0
module bwt_mod
  use iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_int, c_signed_char
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: etx = 0

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

contains

  subroutine generate_sequence(sequence, n)
    integer(c_signed_char), intent(out) :: sequence(0:)
    integer, intent(in) :: n
    integer :: i
    integer(c_int) :: pick
    integer(c_signed_char), parameter :: alphabet(0:3) = [ &
      int(iachar('A'), c_signed_char), int(iachar('T'), c_signed_char), &
      int(iachar('C'), c_signed_char), int(iachar('G'), c_signed_char)]

    call c_srand(123_c_int)
    do i = 0, n - 1
      pick = mod(c_rand(), 4_c_int)
      sequence(i) = alphabet(int(pick))
    end do
    sequence(n) = etx
  end subroutine generate_sequence

  subroutine generate_table(table, table_size, n)
    integer(int32), intent(out) :: table(0:)
    integer, intent(in) :: table_size, n
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, table_size - 1
      if (i < n) then
        table(i) = i
      else
        table(i) = -1
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine generate_table

  subroutine bitonic_sort_step(table, table_size, j, k, genome, n)
    integer(int32), intent(inout) :: table(0:)
    integer(c_signed_char), intent(in) :: genome(0:)
    integer, intent(in) :: table_size, j, k, n
    integer :: i, ixj, t1, t2, cmp_a, cmp_b, offset
    logical :: forward, swap_entries

    !$omp target teams distribute parallel do thread_limit(block_size) private(ixj, forward, t1, t2, cmp_a, cmp_b, offset, swap_entries)
    do i = 0, table_size - 1
      ixj = ieor(i, j)
      if (i < ixj) then
        forward = iand(i, k) == 0
        t1 = table(i)
        t2 = table(ixj)
        cmp_a = merge(t2, t1, forward)
        cmp_b = merge(t1, t2, forward)
        if (cmp_a < 0) then
          swap_entries = .false.
        else if (cmp_b < 0) then
          swap_entries = .true.
        else
          swap_entries = .false.
          do offset = 0, n - 1
            if (genome(mod(cmp_a + offset, n)) /= genome(mod(cmp_b + offset, n))) then
              swap_entries = genome(mod(cmp_a + offset, n)) < genome(mod(cmp_b + offset, n))
              exit
            end if
          end do
        end if
        if (swap_entries) then
          table(i) = t2
          table(ixj) = t1
        end if
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine bitonic_sort_step

  logical function compare_rotations_cpu(a, b, genome, n) result(is_less)
    integer, intent(in) :: a, b, n
    integer(c_signed_char), intent(in) :: genome(0:)
    integer :: i

    if (a < 0) then
      is_less = .false.
      return
    end if
    if (b < 0) then
      is_less = .true.
      return
    end if

    do i = 0, n - 1
      if (genome(mod(a + i, n)) /= genome(mod(b + i, n))) then
        is_less = genome(mod(a + i, n)) < genome(mod(b + i, n))
        return
      end if
    end do
    is_less = .false.
  end function compare_rotations_cpu

  subroutine reconstruct_sequence(table, sequence, transformed, n)
    integer(int32), intent(in) :: table(0:)
    integer(c_signed_char), intent(in) :: sequence(0:)
    integer(c_signed_char), intent(out) :: transformed(0:)
    integer, intent(in) :: n
    integer :: i

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, n - 1
      transformed(i) = sequence(mod(n + table(i) - 1, n))
    end do
    !$omp end target teams distribute parallel do
  end subroutine reconstruct_sequence

  subroutine bwt_gpu(sequence, transformed, suffix_table, n)
    integer(c_signed_char), intent(in) :: sequence(0:)
    integer(c_signed_char), intent(out) :: transformed(0:)
    integer(int32), allocatable, intent(out) :: suffix_table(:)
    integer, intent(in) :: n
    integer :: table_size, j, k

    table_size = 1
    do while (table_size < n)
      table_size = table_size * 2
    end do

    allocate(suffix_table(0:table_size - 1))

    !$omp target data map(from: suffix_table(0:table_size - 1), transformed(0:n - 1)) map(to: sequence(0:n - 1))
    call generate_table(suffix_table, table_size, n)

    k = 2
    do while (k <= table_size)
      j = k / 2
      do while (j > 0)
        call bitonic_sort_step(suffix_table, table_size, j, k, sequence, n)
        j = j / 2
      end do
      k = k * 2
    end do

    call reconstruct_sequence(suffix_table, sequence, transformed, n)
    !$omp end target data
  end subroutine bwt_gpu

  subroutine bwt_cpu(sequence, transformed, n)
    integer(c_signed_char), intent(in) :: sequence(0:)
    integer(c_signed_char), intent(out) :: transformed(0:)
    integer, intent(in) :: n
    integer(int32), allocatable :: next(:)
    integer :: i, head

    allocate(next(0:n - 1))
    do i = 0, n - 1
      if (i == n - 1) then
        next(i) = -1
      else
        next(i) = i + 1
      end if
    end do
    head = list_sort(0, next, sequence, n)

    i = 0
    do while (head >= 0)
      transformed(i) = sequence(mod(n + head - 1, n))
      head = next(head)
      i = i + 1
    end do

    deallocate(next)
  end subroutine bwt_cpu

  recursive integer function list_sort(head, next, genome, n) result(sorted)
    integer, intent(in) :: head, n
    integer(int32), intent(inout) :: next(0:)
    integer(c_signed_char), intent(in) :: genome(0:)
    integer :: slow, fast, mid

    if (head < 0 .or. next(head) < 0) then
      sorted = head
      return
    end if

    slow = head
    fast = next(head)
    do while (fast >= 0)
      fast = next(fast)
      if (fast >= 0) then
        slow = next(slow)
        fast = next(fast)
      end if
    end do

    mid = next(slow)
    next(slow) = -1
    sorted = merge_lists(list_sort(head, next, genome, n), list_sort(mid, next, genome, n), next, genome, n)
  end function list_sort

  recursive integer function merge_lists(left, right, next, genome, n) result(head)
    integer, intent(in) :: left, right, n
    integer(int32), intent(inout) :: next(0:)
    integer(c_signed_char), intent(in) :: genome(0:)
    integer :: chosen, rest

    if (left < 0) then
      head = right
    else if (right < 0) then
      head = left
    else if (compare_rotations_cpu(right, left, genome, n)) then
      chosen = right
      rest = next(right)
      next(chosen) = merge_lists(left, rest, next, genome, n)
      head = chosen
    else
      chosen = left
      rest = next(left)
      next(chosen) = merge_lists(rest, right, next, genome, n)
      head = chosen
    end if
  end function merge_lists

  logical function arrays_equal(a, b, n) result(equal)
    integer(c_signed_char), intent(in) :: a(0:), b(0:)
    integer, intent(in) :: n
    integer :: i

    equal = .true.
    do i = 0, n - 1
      if (a(i) /= b(i)) then
        equal = .false.
        return
      end if
    end do
  end function arrays_equal

end module bwt_mod

program main
  use iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_signed_char
  use omp_lib
  use bwt_mod
  implicit none

  integer :: n, arg_count
  character(len=64) :: arg
  integer(c_signed_char), allocatable :: sequence(:), cpu_seq(:), gpu_seq(:)
  integer(int32), allocatable :: suffix_table(:)
  real(real64) :: start_time, cpu_ms, gpu_ms

  arg_count = command_argument_count()
  if (arg_count > 0) then
    call get_command_argument(1, arg)
    read(arg, *) n
  else
    n = 1000000
  end if

  write(*,'(A,I0)') 'running a sample sequence of length ', n

  allocate(sequence(0:n), cpu_seq(0:n), gpu_seq(0:n))
  call generate_sequence(sequence, n)

  start_time = omp_get_wtime()
  call bwt_cpu(sequence, cpu_seq, n)
  cpu_ms = (omp_get_wtime() - start_time) * 1000.0_real64

  start_time = omp_get_wtime()
  call bwt_gpu(sequence, gpu_seq, suffix_table, n)
  gpu_ms = (omp_get_wtime() - start_time) * 1000.0_real64

  write(*,'(A,I0,A)') 'Host time: ', int(cpu_ms), ' ms'
  write(*,'(A,I0,A)') 'Device time: ', int(gpu_ms), ' ms'

  if (arrays_equal(cpu_seq, gpu_seq, n)) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(sequence, cpu_seq, gpu_seq, suffix_table)
end program main
