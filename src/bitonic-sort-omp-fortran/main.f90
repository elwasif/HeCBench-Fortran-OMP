! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0
  integer :: n, seed, size, i
  integer(int32), allocatable :: data_cpu(:), data_gpu(:)
  logical :: ok

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

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    call usage(trim(arg0))
    stop 1
  end if

  n = read_arg(1)
  seed = read_arg(2)
  if (n < 0 .or. n >= 30) then
    call usage(trim(arg0))
    stop 1
  end if

  size = 2 ** n
  print '(A)'
  print '(A,I0,A,I0)', 'Array size: ', size, ', seed: ', seed

  allocate(data_cpu(size), data_gpu(size))
  call c_srand(int(seed, c_int))
  do i = 1, size
    data_cpu(i) = int(mod(c_rand(), 1000_c_int), int32)
    data_gpu(i) = data_cpu(i)
  end do

  print '(A)', 'Bitonic sort (parallel)..'
  call parallel_bitonic_sort(data_gpu, n)

  print '(A)', 'Bitonic sort (serial)..'
  call bitonic_sort(data_cpu, n)

  ok = .true.
  do i = 1, size
    if (data_gpu(i) /= data_cpu(i)) then
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

  deallocate(data_cpu, data_gpu)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine usage(program_name)
    character(len=*), intent(in) :: program_name

    print '(A)', ' Incorrect parameters'
    print '(3A)', ' Usage: ', trim(program_name), ' n k '
    print '(A)', ''
    print '(A)', ' n: Integer exponent presenting the size of the input array. The number of element in'
    print '(A)', '    the array must be power of 2 (e.g., 1, 2, 4, ...). Please enter the corresponding'
    print '(A)', '    exponent between 0 and 29.'
    print '(A)', ' k: Seed used to generate a random sequence.'
  end subroutine usage

  subroutine parallel_bitonic_sort(input, n)
    integer(int32), intent(inout) :: input(:)
    integer, intent(in) :: n
    integer :: size, step, stage, seq_len, two_power
    real(real64) :: start_time, end_time

    size = 2 ** n
    !$omp target data map(tofrom: input(1:size))
    start_time = omp_get_wtime()
    do step = 0, n - 1
      do stage = step, 0, -1
        seq_len = 2 ** (stage + 1)
        two_power = 2 ** (step - stage)
        call bitonic_stage_device(input, size, seq_len, two_power)
      end do
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Total kernel execution time: ', (end_time - start_time) * 1.0e3_real64, ' (ms)'
    !$omp end target data
  end subroutine parallel_bitonic_sort

  subroutine bitonic_stage_device(input, size, seq_len, two_power)
    integer(int32), intent(inout) :: input(:)
    integer, intent(in) :: size, seq_len, two_power
    integer :: idx, seq_num, swapped_ele, h_len, odd, temp
    logical :: increasing

    h_len = seq_len / 2
    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(seq_num, swapped_ele, odd, temp, increasing)
    do idx = 0, size - 1
      seq_num = idx / seq_len
      swapped_ele = -1
      if (idx < (seq_len * seq_num) + h_len) swapped_ele = idx + h_len
      odd = seq_num / two_power
      increasing = mod(odd, 2) == 0

      if (swapped_ele /= -1) then
        if (((input(idx + 1) > input(swapped_ele + 1)) .and. increasing) .or. &
            ((input(idx + 1) < input(swapped_ele + 1)) .and. (.not. increasing))) then
          temp = input(idx + 1)
          input(idx + 1) = input(swapped_ele + 1)
          input(swapped_ele + 1) = temp
        end if
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine bitonic_stage_device

  subroutine bitonic_sort(array, n)
    integer(int32), intent(inout) :: array(:)
    integer, intent(in) :: n
    integer :: step, stage, num_sequence, sequence_len

    do step = 0, n - 1
      do stage = step, 0, -1
        num_sequence = 2 ** (n - stage - 1)
        sequence_len = 2 ** (stage + 1)
        call swap_elements(step, stage, num_sequence, sequence_len, array)
      end do
    end do
  end subroutine bitonic_sort

  subroutine swap_elements(step, stage, num_sequence, seq_len, array)
    integer, intent(in) :: step, stage, num_sequence, seq_len
    integer(int32), intent(inout) :: array(:)
    integer :: seq_num, odd, h_len, idx, swapped_ele, temp
    logical :: increasing

    h_len = seq_len / 2
    do seq_num = 0, num_sequence - 1
      odd = seq_num / (2 ** (step - stage))
      increasing = mod(odd, 2) == 0
      do idx = seq_num * seq_len, seq_num * seq_len + h_len - 1
        swapped_ele = idx + h_len
        if (((array(idx + 1) > array(swapped_ele + 1)) .and. increasing) .or. &
            ((array(idx + 1) < array(swapped_ele + 1)) .and. (.not. increasing))) then
          temp = array(idx + 1)
          array(idx + 1) = array(swapped_ele + 1)
          array(swapped_ele + 1) = temp
        end if
      end do
    end do
  end subroutine swap_elements

end program main
