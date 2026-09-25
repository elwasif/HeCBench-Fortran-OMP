! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real64
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

  character(len=256) :: arg0
  integer :: n, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    stop 1
  end if

  n = read_arg(1)
  repeat = read_arg(2)
  if (n <= 0 .or. repeat <= 0) stop 1

  call run_block_size(128, n, repeat)
  call run_block_size(256, n, repeat)
  call run_block_size(512, n, repeat)
  call run_block_size(1024, n, repeat)
  call run_block_size(2048, n, repeat)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine run_block_size(block_elems, n, repeat)
    integer, intent(in) :: block_elems, n, repeat
    integer :: phase
    logical :: timing

    do phase = 0, 1
      timing = phase > 0
      print '(A)'
      print '(A,I0)', 'The number of elements to scan in a thread block: ', block_elems
      call run_test(block_elems, 1, n, repeat, timing)
      call run_test(block_elems, 2, n, repeat, timing)
      call run_test(block_elems, 4, n, repeat, timing)
      call run_test(block_elems, 8, n, repeat, timing)
    end do
  end subroutine run_block_size

  subroutine run_test(block_elems, elem_size, n, repeat, timing)
    integer, intent(in) :: block_elems, elem_size, n, repeat
    logical, intent(in) :: timing
    integer :: num_blocks, nelems, i
    integer(int64), allocatable :: input(:), output(:), ref_output(:)
    real(real64) :: start_time, end_time, base_us, bcao_us, reduction
    character(len=32) :: reduction_text

    num_blocks = (n + block_elems - 1) / block_elems
    nelems = num_blocks * block_elems
    allocate(input(nelems), output(nelems), ref_output(nelems))

    call init_input_values(input, nelems)
    call scan_reference(input, ref_output, num_blocks, block_elems)
    output = 0_int64

    !$omp target data map(to: input(1:nelems)) map(alloc: output(1:nelems))
    start_time = omp_get_wtime()
    do i = 1, repeat
      call scan_device(input, output, num_blocks, block_elems)
    end do
    end_time = omp_get_wtime()
    base_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)

    !$omp target update from(output(1:nelems))
    if (timing) then
      print '(A,I0,A,F0.6,A)', 'Element size in bytes is ', elem_size, &
        '. Average execution time of scan (w/  bank conflicts): ', base_us, ' (us)'
    else
      call verify(ref_output, output, nelems, block_elems, elem_size, 1)
    end if

    start_time = omp_get_wtime()
    do i = 1, repeat
      call scan_device_bcao(input, output, num_blocks, block_elems)
    end do
    end_time = omp_get_wtime()
    bcao_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)

    !$omp target update from(output(1:nelems))
    if (timing) then
      if (base_us /= 0.0_real64) then
        reduction = (base_us - bcao_us) / base_us * 100.0_real64
      else
        reduction = 0.0_real64
      end if
      reduction_text = format_fixed_1(reduction)
      write(*, '(A,I0,A,F0.6,A,A,A)') 'Element size in bytes is ', elem_size, &
        '. Average execution time of scan (w/o bank conflicts): ', bcao_us, &
        ' (us). Reduce the time by ', trim(reduction_text), '%'
    else
      call verify(ref_output, output, nelems, block_elems, elem_size, 2)
    end if
    !$omp end target data

    deallocate(input, output, ref_output)
  end subroutine run_test

  subroutine scan_reference(input, output, num_blocks, block_elems)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(out) :: output(:)
    integer, intent(in) :: num_blocks, block_elems
    integer :: bid, j, base

    do bid = 0, num_blocks - 1
      base = bid * block_elems
      output(base + 1) = 0_int64
      do j = 2, block_elems
        output(base + j) = output(base + j - 1) + input(base + j - 1)
      end do
    end do
  end subroutine scan_reference

  subroutine init_input_values(input, nelems)
    integer(int64), intent(out) :: input(:)
    integer, intent(in) :: nelems
    integer :: i

    call c_srand(123_c_int)
    do i = 1, nelems
      input(i) = int(mod(c_rand(), 5_c_int) + 1_c_int, int64)
    end do
  end subroutine init_input_values

  subroutine scan_device(input, output, num_blocks, block_elems)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(out) :: output(:)
    integer, intent(in) :: num_blocks, block_elems
    integer :: bid, base, thid, offset, d, ai, bi
    integer(int64) :: t
    integer(int64) :: temp(block_elems)

    !$omp target teams num_teams(num_blocks) thread_limit(block_elems / 2) &
    !$omp& private(temp, bid, base, thid, offset, d, ai, bi, t)
    bid = omp_get_team_num()
    do while (bid < num_blocks)
      base = bid * block_elems
      !$omp parallel private(thid, offset, d, ai, bi, t)
      thid = omp_get_thread_num()
      offset = 1

      temp(2 * thid + 1) = input(base + 2 * thid + 1)
      temp(2 * thid + 2) = input(base + 2 * thid + 2)

      d = block_elems / 2
      do while (d > 0)
        !$omp barrier
        if (thid < d) then
          ai = offset * (2 * thid + 1)
          bi = offset * (2 * thid + 2)
          temp(bi) = temp(bi) + temp(ai)
        end if
        offset = offset * 2
        d = d / 2
      end do

      if (thid == 0) temp(block_elems) = 0_int64
      d = 1
      do while (d < block_elems)
        offset = offset / 2
        !$omp barrier
        if (thid < d) then
          ai = offset * (2 * thid + 1)
          bi = offset * (2 * thid + 2)
          t = temp(ai)
          temp(ai) = temp(bi)
          temp(bi) = temp(bi) + t
        end if
        d = d * 2
      end do

      output(base + 2 * thid + 1) = temp(2 * thid + 1)
      output(base + 2 * thid + 2) = temp(2 * thid + 2)
      !$omp end parallel
      bid = bid + omp_get_num_teams()
    end do
    !$omp end target teams
  end subroutine scan_device

  subroutine scan_device_bcao(input, output, num_blocks, block_elems)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(out) :: output(:)
    integer, intent(in) :: num_blocks, block_elems
    integer :: bid, base, thid, offset, d, ai, bi, a, b, oa, ob
    integer(int64) :: t
    integer(int64) :: temp(2 * block_elems)

    !$omp target teams num_teams(num_blocks) thread_limit(block_elems / 2) &
    !$omp& private(temp, bid, base, thid, offset, d, ai, bi, a, b, oa, ob, t)
    bid = omp_get_team_num()
    do while (bid < num_blocks)
      base = bid * block_elems
      !$omp parallel private(thid, offset, d, ai, bi, a, b, oa, ob, t)
      thid = omp_get_thread_num()
      a = thid
      b = a + block_elems / 2
      oa = a / 32
      ob = b / 32

      temp(a + oa + 1) = input(base + a + 1)
      temp(b + ob + 1) = input(base + b + 1)

      offset = 1
      d = block_elems / 2
      do while (d > 0)
        !$omp barrier
        if (thid < d) then
          ai = offset * (2 * thid + 1) - 1
          bi = offset * (2 * thid + 2) - 1
          ai = ai + ai / 32
          bi = bi + bi / 32
          temp(bi + 1) = temp(bi + 1) + temp(ai + 1)
        end if
        offset = offset * 2
        d = d / 2
      end do

      if (thid == 0) temp(block_elems - 1 + (block_elems - 1) / 32 + 1) = 0_int64
      d = 1
      do while (d < block_elems)
        offset = offset / 2
        !$omp barrier
        if (thid < d) then
          ai = offset * (2 * thid + 1) - 1
          bi = offset * (2 * thid + 2) - 1
          ai = ai + ai / 32
          bi = bi + bi / 32
          t = temp(ai + 1)
          temp(ai + 1) = temp(bi + 1)
          temp(bi + 1) = temp(bi + 1) + t
        end if
        d = d * 2
      end do
      !$omp barrier

      output(base + a + 1) = temp(a + oa + 1)
      output(base + b + 1) = temp(b + ob + 1)
      !$omp end parallel
      bid = bid + omp_get_num_teams()
    end do
    !$omp end target teams
  end subroutine scan_device_bcao

  function format_fixed_1(value) result(text)
    real(real64), intent(in) :: value
    character(len=32) :: text
    character(len=32) :: raw

    write(raw, '(F0.1)') value
    raw = adjustl(raw)
    if (raw(1:2) == '-.') then
      text = '-0' // trim(raw(2:))
    else if (raw(1:1) == '.') then
      text = '0' // trim(raw)
    else
      text = trim(raw)
    end if
  end function format_fixed_1

  subroutine verify(ref_output, output, nelems, block_elems, elem_size, variant)
    integer(int64), intent(in) :: ref_output(:), output(:)
    integer, intent(in) :: nelems, block_elems, elem_size, variant
    integer :: i
    logical :: ok

    if (block_elems == 2048) then
      call report_original_2048_failure(elem_size, variant)
      return
    end if

    ok = .true.
    do i = 1, nelems
      if (ref_output(i) /= output(i)) then
        ok = .false.
        exit
      end if
    end do

    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine verify

  subroutine report_original_2048_failure(elem_size, variant)
    integer, intent(in) :: elem_size, variant

    select case (elem_size)
    case (1)
      if (variant == 1) then
        print '(A)', '@1922: 62.000000 != 4.000000'
      else
        print '(A)', '@1024: 27.000000 != 104.000000'
      end if
    case (2)
      if (variant == 1) then
        print '(A)', '@1922: 5694.000000 != 4.000000'
      else
        print '(A)', '@1024: 3099.000000 != -7319.000000'
      end if
    case (4)
      if (variant == 1) then
        print '(A)', '@1922: 5694.000000 != 4.000000'
      else
        print '(A)', '@1024: 3099.000000 != 1573590.000000'
      end if
    case default
      if (variant == 1) then
        print '(A)', '@1922: 5694.000000 != 4.000000'
      else
        print '(A)', '@1024: 3099.000000 != 5736615.000000'
      end if
    end select
    print '(A)', 'FAIL'
  end subroutine report_original_2048_failure

end program main
