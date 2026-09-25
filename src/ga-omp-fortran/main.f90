! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int8, int32, real64
  use omp_lib
  implicit none

  integer, parameter :: batch_size = 1024

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

  character(len=128) :: arg1, arg2, arg3, arg4
  integer :: target_size, query_size, coarse_match_length, coarse_match_threshold
  integer(int8), allocatable :: target_sequence(:), query_sequence(:)
  integer(int8) :: device_batch(batch_size), host_batch(batch_size)
  integer :: current_position, max_searchable_length, end_position, length
  real(real64) :: start_time, end_time, total_time
  logical :: ok

  if (command_argument_count() /= 4) then
    print '(A)', 'Usage: ./main <target sequence length> <query sequence length> <coarse match length> <coarse match threshold>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  read(arg1, *) target_size
  read(arg2, *) query_size
  read(arg3, *) coarse_match_length
  read(arg4, *) coarse_match_threshold
  if (target_size <= coarse_match_length .or. query_size < coarse_match_length) stop 1

  allocate(target_sequence(target_size), query_sequence(query_size))
  call initialize_sequences(target_sequence, query_sequence)

  max_searchable_length = target_size - coarse_match_length
  current_position = 0
  total_time = 0.0_real64
  ok = .true.

  !$omp target data map(to: target_sequence(1:target_size), query_sequence(1:query_size)) map(alloc: device_batch(1:batch_size))
  do while (current_position < max_searchable_length)
    call clear_device_batch(device_batch)
    host_batch = 0_int8

    end_position = current_position + batch_size
    if (end_position >= max_searchable_length) end_position = max_searchable_length
    length = end_position - current_position

    start_time = omp_get_wtime()
    call ga_kernel(target_sequence, query_sequence, device_batch, length, query_size, &
        coarse_match_length, coarse_match_threshold, current_position)
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)

    call ga_reference(target_sequence, query_sequence, host_batch, length, query_size, &
        coarse_match_length, coarse_match_threshold, current_position)

    !$omp target update from(device_batch(1:batch_size))
    if (any(device_batch /= host_batch)) then
      ok = .false.
      exit
    end if

    current_position = end_position
  end do
  !$omp end target data

  print '(A,F0.6,A)', 'Total kernel execution time ', total_time, ' (s)'
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(target_sequence, query_sequence)

contains

  subroutine initialize_sequences(target_sequence, query_sequence)
    integer(int8), intent(out) :: target_sequence(:), query_sequence(:)
    integer(int8), parameter :: alphabet(4) = [ &
        int(iachar('A'), int8), int(iachar('C'), int8), &
        int(iachar('T'), int8), int(iachar('G'), int8)]
    integer :: i

    call c_srand(123_c_int)
    do i = 1, size(target_sequence)
      target_sequence(i) = alphabet(mod(c_rand(), 4) + 1)
    end do
    do i = 1, size(query_sequence)
      query_sequence(i) = alphabet(mod(c_rand(), 4) + 1)
    end do
  end subroutine initialize_sequences

  subroutine clear_device_batch(device_batch)
    integer(int8), intent(inout) :: device_batch(:)
    integer :: i

    !$omp target teams distribute parallel do thread_limit(256) private(i)
    do i = 1, batch_size
      device_batch(i) = 0_int8
    end do
    !$omp end target teams distribute parallel do
  end subroutine clear_device_batch

  subroutine ga_kernel(target, query, batch_result, length, query_sequence_length, &
      coarse_match_length, coarse_match_threshold, current_position)
    integer(int8), intent(in) :: target(:), query(:)
    integer(int8), intent(inout) :: batch_result(:)
    integer, intent(in) :: length, query_sequence_length, coarse_match_length
    integer, intent(in) :: coarse_match_threshold, current_position
    integer :: tid, i, j, distance, max_length
    logical :: match

    max_length = query_sequence_length - coarse_match_length
    !$omp target teams distribute parallel do thread_limit(256) private(tid, i, j, distance, match) firstprivate(max_length)
    do tid = 0, length - 1
      match = .false.
      do i = 0, max_length
        distance = 0
        do j = 0, coarse_match_length - 1
          if (target(current_position + tid + j + 1) /= query(i + j + 1)) distance = distance + 1
        end do
        if (distance < coarse_match_threshold) then
          match = .true.
          exit
        end if
      end do
      if (match) batch_result(tid + 1) = 1_int8
    end do
    !$omp end target teams distribute parallel do
  end subroutine ga_kernel

  subroutine ga_reference(target, query, batch_result, length, query_sequence_length, &
      coarse_match_length, coarse_match_threshold, current_position)
    integer(int8), intent(in) :: target(:), query(:)
    integer(int8), intent(inout) :: batch_result(:)
    integer, intent(in) :: length, query_sequence_length, coarse_match_length
    integer, intent(in) :: coarse_match_threshold, current_position
    integer :: tid, i, j, distance, max_length
    logical :: match

    max_length = query_sequence_length - coarse_match_length
    do tid = 0, length - 1
      match = .false.
      do i = 0, max_length
        distance = 0
        do j = 0, coarse_match_length - 1
          if (target(current_position + tid + j + 1) /= query(i + j + 1)) distance = distance + 1
        end do
        if (distance < coarse_match_threshold) then
          match = .true.
          exit
        end if
      end do
      if (match) batch_result(tid + 1) = 1_int8
    end do
  end subroutine ga_reference

end program main
