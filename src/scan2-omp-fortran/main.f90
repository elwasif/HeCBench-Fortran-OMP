! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

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

  integer, parameter :: group_size = 256
  character(len=256) :: arg
  integer :: iterations, length, block_size

  if (command_argument_count() /= 3) then
    call get_command_argument(0, arg)
    write(*,'(A,A,A)') 'Usage: ', trim(arg), ' <repeat> <input length> <block size>'
    stop 1
  end if

  iterations = read_int_arg(1)
  length = read_int_arg(2)
  block_size = read_int_arg(3)

  if (iterations < 1) then
    write(*,'(A)') 'Error, iterations cannot be 0 or negative. Exiting..'
    stop 255
  end if

  if (.not. is_power_of_2(length)) length = round_to_power_of_2(length)

  if ((length / block_size > group_size) .and. (.not. is_power_of_2(length))) then
    write(*,'(A,I0)') 'Invalid length: ', length
    stop 255
  end if

  call run_scan(iterations, length, block_size)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  logical function is_power_of_2(value) result(ok)
    integer, intent(in) :: value
    ok = value > 0 .and. iand(value, -value) == value
  end function is_power_of_2

  integer function round_to_power_of_2(value) result(out)
    integer, intent(in) :: value
    integer :: v, i, bytes
    v = value - 1
    bytes = storage_size(value) / 8
    do i = 0, bytes - 1
      v = ior(v, ishft(v, -(2 ** i)))
    end do
    out = v + 1
  end function round_to_power_of_2

  subroutine run_scan(iterations, length, block_size_in)
    integer, intent(in) :: iterations, length, block_size_in
    integer :: block_size, pass, output_buffer_size, block_sum_buffer_size, temp_length
    integer :: i, n, update_len
    integer, allocatable :: output_offset(:), block_sum_offset(:)
    real(real32), allocatable :: input_buffer(:), output_buffer(:), block_sum_buffer(:), temp_buffer(:)
    real(real32), allocatable :: verification_output(:)
    real(real64) :: t, start_time, elapsed_us

    block_size = min(block_size_in, length / 2)
    t = log(real(length, real64)) / log(real(block_size, real64))
    pass = int(t)
    if (abs(t - real(pass, real64)) < 1.0e-7_real64) pass = pass - 1
    if (pass < 1) pass = 1

    allocate(output_offset(0:pass-1), block_sum_offset(0:pass-1))
    output_buffer_size = 0
    do i = 0, pass - 1
      output_offset(i) = output_buffer_size
      output_buffer_size = output_buffer_size + int(real(length, real64) / real(block_size, real64) ** i)
    end do

    block_sum_buffer_size = 0
    do i = 0, pass - 1
      block_sum_offset(i) = block_sum_buffer_size
      block_sum_buffer_size = block_sum_buffer_size + int(real(length, real64) / real(block_size, real64) ** (i + 1))
    end do

    temp_length = int(real(length, real64) / real(block_size, real64) ** pass)

    allocate(input_buffer(0:length-1), output_buffer(0:output_buffer_size-1), &
             block_sum_buffer(0:block_sum_buffer_size-1), temp_buffer(0:temp_length-1))
    call fill_random(input_buffer, length, 0.0_real32, 255.0_real32)
    output_buffer = 0.0_real32
    block_sum_buffer = 0.0_real32
    temp_buffer = 0.0_real32

    write(*,'(A,I0,A)') 'Executing kernel for ', iterations, ' iterations'
    write(*,'(A)') '-------------------------------------------'

    !$omp target data map(to: input_buffer(0:length-1)) &
    !$omp& map(alloc: temp_buffer(0:temp_length-1), block_sum_buffer(0:block_sum_buffer_size-1), output_buffer(0:output_buffer_size-1))
      start_time = omp_get_wtime()
      do n = 1, iterations
        call b_scan(block_size, length, input_buffer, output_buffer(output_offset(0):), block_sum_buffer(block_sum_offset(0):))

        do i = 1, pass - 1
          call b_scan(block_size, int(real(length, real64) / real(block_size, real64) ** i), &
            block_sum_buffer(block_sum_offset(i - 1):), output_buffer(output_offset(i):), block_sum_buffer(block_sum_offset(i):))
        end do

        call p_scan(block_size, temp_length, block_sum_buffer(block_sum_offset(pass - 1):), temp_buffer)

        call b_addition(block_size, int(real(length, real64) / real(block_size, real64) ** (pass - 1)), &
          temp_buffer, output_buffer(output_offset(pass - 1):))

        do i = pass - 1, 1, -1
          call b_addition(block_size, int(real(length, real64) / real(block_size, real64) ** (i - 1)), &
            output_buffer(output_offset(i):), output_buffer(output_offset(i - 1):))
        end do
      end do
      elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64
      write(*,'(A,F0.6,A)') 'Average execution time of scan kernels: ', elapsed_us / real(iterations, real64), ' (us)'

      if (pass == 1) then
        update_len = output_buffer_size
      else
        update_len = output_offset(1)
      end if
      !$omp target update from(output_buffer(0:update_len-1))
    !$omp end target data

    allocate(verification_output(0:length-1))
    verification_output = 0.0_real32
    call scan_large_arrays_cpu_reference(verification_output, input_buffer, length)

    if (compare(output_buffer, verification_output, length, 0.001_real32)) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    deallocate(verification_output, input_buffer, temp_buffer, block_sum_buffer, block_sum_offset, output_buffer, output_offset)
  end subroutine run_scan

  subroutine b_scan(block_size, len, input, output, sum_buffer)
    integer, intent(in) :: block_size, len
    real(real32), intent(in) :: input(0:)
    real(real32), intent(inout) :: output(0:), sum_buffer(0:)
    integer :: teams
    integer :: tid, bid, gid, stride
    real(real32) :: block(0:255), cache0, cache1

    teams = len / block_size
    !$omp target teams num_teams(teams) thread_limit(block_size / 2) private(block)
      !$omp parallel private(tid, bid, gid, stride, cache0, cache1)
        tid = omp_get_thread_num()
        bid = omp_get_team_num()
        gid = bid * block_size / 2 + tid

        block(2 * tid) = input(2 * gid)
        block(2 * tid + 1) = input(2 * gid + 1)
        !$omp barrier

        cache0 = block(0)
        cache1 = cache0 + block(1)

        stride = 1
        do while (stride < block_size)
          if (2 * tid >= stride) then
            cache0 = block(2 * tid - stride) + block(2 * tid)
            cache1 = block(2 * tid + 1 - stride) + block(2 * tid + 1)
          end if
          !$omp barrier

          block(2 * tid) = cache0
          block(2 * tid + 1) = cache1
          !$omp barrier

          stride = stride * 2
        end do

        sum_buffer(bid) = block(block_size - 1)

        if (tid == 0) then
          output(2 * gid) = 0.0_real32
          output(2 * gid + 1) = block(2 * tid)
        else
          output(2 * gid) = block(2 * tid - 1)
          output(2 * gid + 1) = block(2 * tid)
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine b_scan

  subroutine p_scan(block_size, len, input, output)
    integer, intent(in) :: block_size, len
    real(real32), intent(in) :: input(0:)
    real(real32), intent(inout) :: output(0:)
    integer :: tid, bid, gid, stride
    real(real32) :: block(0:3), cache0, cache1

    !$omp target teams num_teams(1) thread_limit(len / 2) private(block)
      !$omp parallel private(tid, bid, gid, stride, cache0, cache1)
        tid = omp_get_thread_num()
        bid = omp_get_team_num()
        gid = bid * len / 2 + tid

        block(2 * tid) = input(2 * gid)
        block(2 * tid + 1) = input(2 * gid + 1)
        !$omp barrier

        cache0 = block(0)
        cache1 = cache0 + block(1)

        stride = 1
        do while (stride < block_size)
          if (2 * tid >= stride) then
            cache0 = block(2 * tid - stride) + block(2 * tid)
            cache1 = block(2 * tid + 1 - stride) + block(2 * tid + 1)
          end if
          !$omp barrier

          block(2 * tid) = cache0
          block(2 * tid + 1) = cache1
          !$omp barrier

          stride = stride * 2
        end do

        if (tid == 0) then
          output(2 * gid) = 0.0_real32
          output(2 * gid + 1) = block(2 * tid)
        else
          output(2 * gid) = block(2 * tid - 1)
          output(2 * gid + 1) = block(2 * tid)
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine p_scan

  subroutine b_addition(block_size, len, input, output)
    integer, intent(in) :: block_size, len
    real(real32), intent(in) :: input(0:)
    real(real32), intent(inout) :: output(0:)
    integer :: tid, bid, gid
    real(real32) :: value

    !$omp target teams num_teams(len / block_size) thread_limit(block_size) private(value)
      !$omp parallel private(tid, bid, gid)
        tid = omp_get_thread_num()
        bid = omp_get_team_num()
        gid = bid * block_size + tid

        if (tid == 0) value = input(bid)
        !$omp barrier

        output(gid) = output(gid) + value
      !$omp end parallel
    !$omp end target teams
  end subroutine b_addition

  subroutine scan_large_arrays_cpu_reference(output, input, length)
    real(real32), intent(inout) :: output(0:)
    real(real32), intent(in) :: input(0:)
    integer, intent(in) :: length
    integer :: i
    output(0) = 0.0_real32
    do i = 1, length - 1
      output(i) = input(i - 1) + output(i - 1)
    end do
  end subroutine scan_large_arrays_cpu_reference

  subroutine fill_random(array, length, range_min, range_max)
    real(real32), intent(out) :: array(0:)
    integer, intent(in) :: length
    real(real32), intent(in) :: range_min, range_max
    integer :: i
    integer(c_int) :: raw
    real(real64) :: range

    range = real(range_max - range_min, real64) + 1.0_real64
    call c_srand(123_c_int)
    do i = 0, length - 1
      raw = c_rand()
      array(i) = range_min + real(range * real(raw, real64) / 2147483648.0_real64, real32)
    end do
  end subroutine fill_random

  logical function compare(ref_data, data, length, epsilon) result(ok)
    real(real32), intent(in) :: ref_data(0:), data(0:), epsilon
    integer, intent(in) :: length
    integer :: i
    real(real32) :: error, ref, diff, norm_ref, norm_error

    error = 0.0_real32
    ref = 0.0_real32
    do i = 1, length - 1
      diff = ref_data(i) - data(i)
      error = error + diff * diff
      ref = ref + ref_data(i) * ref_data(i)
    end do
    norm_ref = sqrt(ref)
    if (abs(ref) < epsilon) then
      ok = .false.
      return
    end if
    norm_error = sqrt(error)
    error = norm_error / norm_ref
    ok = error < epsilon
  end function compare

end program main
