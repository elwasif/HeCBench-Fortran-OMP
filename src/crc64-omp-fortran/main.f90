! SPDX-License-Identifier: CC0-1.0
module crc64_mod
  use iso_fortran_env, only: int32, int64, real64
  implicit none

  integer(int64), parameter :: crc64_poly = int(z'c96c5795d7870f42', int64)
  integer(int64), parameter :: all_ones = int(z'ffffffffffffffff', int64)
  integer(int64), parameter :: high_bit = int(z'8000000000000000', int64)
  integer(int64), parameter :: x_term = int(z'4000000000000000', int64)
  integer(int64), parameter :: mask48 = int(z'0000ffffffffffff', int64)
  integer(int64), parameter :: lcg_a = 25214903917_int64
  integer(int64), parameter :: lcg_c = 11_int64
  integer(int64), parameter :: base24 = 16777216_int64
  real(real64), parameter :: two48 = 281474976710656.0_real64
  logical, save :: tables_ready = .false.
  integer(int64), save :: crc64_table(0:3, 0:255)
  integer(int64), save :: crc64_interleaved_table(0:3, 0:255)

contains

  subroutine srand48_state(seed, state)
    integer(int32), intent(in) :: seed
    integer(int64), intent(out) :: state

    state = iand(ior(shiftl(int(seed, int64), 16), int(z'330e', int64)), mask48)
  end subroutine srand48_state

  subroutine next_drand48(state, value)
    integer(int64), intent(inout) :: state
    real(real64), intent(out) :: value
    integer(int64) :: a0, a1, x0, x1, mid

    a0 = iand(lcg_a, base24 - 1_int64)
    a1 = shiftr(lcg_a, 24)
    x0 = iand(state, base24 - 1_int64)
    x1 = shiftr(state, 24)
    mid = iand(a0 * x1 + a1 * x0, base24 - 1_int64)
    state = iand(a0 * x0 + shiftl(mid, 24) + lcg_c, mask48)
    value = real(state, real64) / two48
  end subroutine next_drand48

  subroutine init_crc64_tables()
    integer :: row, idx
    integer(int64) :: entry

    if (tables_ready) return

    do idx = 0, 255
      entry = int(idx, int64)
      do row = 1, 8
        if (iand(entry, 1_int64) /= 0_int64) then
          entry = ieor(shiftr(entry, 1), crc64_poly)
        else
          entry = shiftr(entry, 1)
        end if
      end do
      crc64_table(3, idx) = entry
    end do

    do row = 2, 0, -1
      do idx = 0, 255
        entry = crc64_table(row + 1, idx)
        crc64_table(row, idx) = ieor(crc64_table(3, int(iand(entry, 255_int64))), shiftr(entry, 8))
      end do
    end do

    crc64_interleaved_table = crc64_table
    do row = 1, 16
      do idx = 0, 255
        entry = crc64_interleaved_table(3, idx)
        crc64_interleaved_table(3, idx) = ieor(crc64_table(3, int(iand(entry, 255_int64))), shiftr(entry, 8))
      end do
    end do
    do row = 2, 0, -1
      do idx = 0, 255
        entry = crc64_interleaved_table(row + 1, idx)
        crc64_interleaved_table(row, idx) = ieor(crc64_table(3, int(iand(entry, 255_int64))), shiftr(entry, 8))
      end do
    end do

    tables_ready = .true.
  end subroutine init_crc64_tables

  pure integer(int32) function crc64_load_le32(buffer, pos) result(w)
    integer(int32), intent(in) :: buffer(:)
    integer(int64), intent(in) :: pos

    w = ior(ior(iand(buffer(pos), 255_int32), shiftl(iand(buffer(pos + 1_int64), 255_int32), 8)), &
            ior(shiftl(iand(buffer(pos + 2_int64), 255_int32), 16), &
                shiftl(iand(buffer(pos + 3_int64), 255_int32), 24)))
  end function crc64_load_le32

  pure integer(int64) function crc64(buffer, first, nbytes, table, interleaved_table) result(cs_out)
    integer(int32), intent(in) :: buffer(:)
    integer(int64), intent(in) :: first, nbytes
    integer(int64), intent(in) :: table(0:3, 0:255), interleaved_table(0:3, 0:255)
    integer(int64) :: cs(0:4), cry, data_pos, end_pos, idx
    integer(int32) :: in_word(0:4)
    integer :: b, i

    cs = 0_int64
    cs(0) = all_ones
    data_pos = first
    end_pos = first + nbytes

    do while (data_pos < end_pos .and. (mod(data_pos - 1_int64, 4_int64) /= 0_int64 .or. end_pos - data_pos < 20_int64))
      idx = iand(ieor(cs(0), int(iand(buffer(data_pos), 255_int32), int64)), 255_int64)
      cs(0) = ieor(table(3, int(idx)), shiftr(cs(0), 8))
      data_pos = data_pos + 1_int64
    end do

    if (data_pos == end_pos) then
      cs_out = ieor(cs(0), all_ones)
      return
    end if

    do i = 0, 4
      in_word(i) = crc64_load_le32(buffer, data_pos + int(4 * i, int64))
    end do
    data_pos = data_pos + 20_int64
    cry = 0_int64

    do while (end_pos - data_pos >= 20_int64)
      cs(0) = ieor(cs(0), cry)

      in_word(0) = ieor(in_word(0), int(iand(cs(0), int(z'00000000ffffffff', int64)), int32))
      cs(1) = ieor(cs(1), shiftr(cs(0), 32))
      cs(0) = interleaved_table(0, iand(in_word(0), 255_int32))
      in_word(0) = shiftr(in_word(0), 8)

      in_word(1) = ieor(in_word(1), int(iand(cs(1), int(z'00000000ffffffff', int64)), int32))
      cs(2) = ieor(cs(2), shiftr(cs(1), 32))
      cs(1) = interleaved_table(0, iand(in_word(1), 255_int32))
      in_word(1) = shiftr(in_word(1), 8)

      in_word(2) = ieor(in_word(2), int(iand(cs(2), int(z'00000000ffffffff', int64)), int32))
      cs(3) = ieor(cs(3), shiftr(cs(2), 32))
      cs(2) = interleaved_table(0, iand(in_word(2), 255_int32))
      in_word(2) = shiftr(in_word(2), 8)

      in_word(3) = ieor(in_word(3), int(iand(cs(3), int(z'00000000ffffffff', int64)), int32))
      cs(4) = ieor(cs(4), shiftr(cs(3), 32))
      cs(3) = interleaved_table(0, iand(in_word(3), 255_int32))
      in_word(3) = shiftr(in_word(3), 8)

      in_word(4) = ieor(in_word(4), int(iand(cs(4), int(z'00000000ffffffff', int64)), int32))
      cry = shiftr(cs(4), 32)
      cs(4) = interleaved_table(0, iand(in_word(4), 255_int32))
      in_word(4) = shiftr(in_word(4), 8)

      do b = 1, 2
        cs(0) = ieor(cs(0), interleaved_table(b, iand(in_word(0), 255_int32)))
        in_word(0) = shiftr(in_word(0), 8)

        cs(1) = ieor(cs(1), interleaved_table(b, iand(in_word(1), 255_int32)))
        in_word(1) = shiftr(in_word(1), 8)

        cs(2) = ieor(cs(2), interleaved_table(b, iand(in_word(2), 255_int32)))
        in_word(2) = shiftr(in_word(2), 8)

        cs(3) = ieor(cs(3), interleaved_table(b, iand(in_word(3), 255_int32)))
        in_word(3) = shiftr(in_word(3), 8)

        cs(4) = ieor(cs(4), interleaved_table(b, iand(in_word(4), 255_int32)))
        in_word(4) = shiftr(in_word(4), 8)
      end do

      cs(0) = ieor(cs(0), interleaved_table(3, iand(in_word(0), 255_int32)))
      in_word(0) = crc64_load_le32(buffer, data_pos)

      cs(1) = ieor(cs(1), interleaved_table(3, iand(in_word(1), 255_int32)))
      in_word(1) = crc64_load_le32(buffer, data_pos + 4_int64)

      cs(2) = ieor(cs(2), interleaved_table(3, iand(in_word(2), 255_int32)))
      in_word(2) = crc64_load_le32(buffer, data_pos + 8_int64)

      cs(3) = ieor(cs(3), interleaved_table(3, iand(in_word(3), 255_int32)))
      in_word(3) = crc64_load_le32(buffer, data_pos + 12_int64)

      cs(4) = ieor(cs(4), interleaved_table(3, iand(in_word(4), 255_int32)))
      in_word(4) = crc64_load_le32(buffer, data_pos + 16_int64)
      data_pos = data_pos + 20_int64
    end do

    cs(0) = ieor(cs(0), cry)

    do i = 0, 4
      if (i > 0) cs(0) = ieor(cs(0), cs(i))
      in_word(i) = ieor(in_word(i), int(iand(cs(0), int(z'00000000ffffffff', int64)), int32))
      cs(0) = shiftr(cs(0), 32)

      do b = 0, 2
        cs(0) = ieor(cs(0), table(b, iand(in_word(i), 255_int32)))
        in_word(i) = shiftr(in_word(i), 8)
      end do

      cs(0) = ieor(cs(0), table(3, iand(in_word(i), 255_int32)))
    end do

    do while (data_pos < end_pos)
      idx = iand(ieor(cs(0), int(iand(buffer(data_pos), 255_int32), int64)), 255_int64)
      cs(0) = ieor(table(3, int(idx)), shiftr(cs(0), 8))
      data_pos = data_pos + 1_int64
    end do

    cs_out = ieor(cs(0), all_ones)
  end function crc64

  pure integer(int64) function crc64_multiply(a_in, b_in) result(r)
    integer(int64), intent(in) :: a_in, b_in
    integer(int64) :: a, b
    integer :: bit

    a = a_in
    b = b_in
    r = 0_int64
    do bit = 1, 64
      if (iand(a, high_bit) /= 0_int64) r = ieor(r, b)
      a = shiftl(a, 1)
      if (iand(b, 1_int64) /= 0_int64) then
        b = ieor(shiftr(b, 1), crc64_poly)
      else
        b = shiftr(b, 1)
      end if
    end do
  end function crc64_multiply

  subroutine crc64_pow_table(pow2)
    integer(int64), intent(out) :: pow2(64)
    integer :: i

    pow2(1) = x_term
    do i = 2, 64
      pow2(i) = crc64_multiply(pow2(i - 1), pow2(i - 1))
    end do
  end subroutine crc64_pow_table

  integer(int64) function crc64_x_pow_n(n) result(r)
    integer(int64), intent(in) :: n
    integer(int64) :: pow2(64), work
    integer :: i

    call crc64_pow_table(pow2)
    r = high_bit
    work = n
    i = 1
    do while (work /= 0_int64)
      if (iand(work, 1_int64) /= 0_int64) r = crc64_multiply(r, pow2(i))
      work = shiftr(work, 1)
      i = i + 1
    end do
  end function crc64_x_pow_n

  integer(int64) function crc64_combine(cs1, cs2, nbytes2) result(cs)
    integer(int64), intent(in) :: cs1, cs2, nbytes2

    cs = ieor(cs2, crc64_multiply(cs1, crc64_x_pow_n(8_int64 * nbytes2)))
  end function crc64_combine

  subroutine crc64_invert(cs, check_bytes)
    integer(int64), intent(in) :: cs
    integer(int32), intent(out) :: check_bytes(8)
    integer(int64) :: work
    integer :: i

    work = ieor(cs, all_ones)
    do i = 1, 8
      check_bytes(i) = int(iand(shiftr(work, 8 * (i - 1)), 255_int64), int32)
    end do
  end subroutine crc64_invert

  integer(int64) function crc64_omp(buffer, nbytes) result(cs)
    integer(int32), intent(in) :: buffer(:)
    integer(int64), intent(in) :: nbytes
    integer(int64), allocatable :: chunk_cs(:), chunk_sz(:)
    integer(int64) :: bpt, first, last
    integer :: tid, nthreads

    if (nbytes <= 2048_int64) then
      call init_crc64_tables()
      cs = crc64(buffer, 1_int64, nbytes, crc64_table, crc64_interleaved_table)
      return
    end if

    call init_crc64_tables()

    nthreads = 96 * 8 * 32
    if (nbytes < int(nthreads, int64) * 1024_int64) then
      nthreads = max(1, int(nbytes / 1024_int64))
    end if

    allocate(chunk_cs(nthreads), chunk_sz(nthreads))
    bpt = nbytes / int(nthreads, int64)

    !$omp target teams distribute parallel do num_teams(max(1, nthreads / 64)) thread_limit(64) &
    !$omp& map(to: buffer(1:nbytes), crc64_table, crc64_interleaved_table) &
    !$omp& map(from: chunk_cs(1:nthreads), chunk_sz(1:nthreads)) private(first, last)
    do tid = 1, nthreads
      first = int(tid - 1, int64) * bpt + 1_int64
      if (tid /= nthreads) then
        last = first + bpt - 1_int64
      else
        last = nbytes
      end if

      chunk_sz(tid) = last - first + 1_int64
      chunk_cs(tid) = crc64(buffer, first, chunk_sz(tid), crc64_table, crc64_interleaved_table)
    end do
    !$omp end target teams distribute parallel do

    cs = chunk_cs(1)
    do tid = 2, nthreads
      cs = crc64_combine(cs, chunk_cs(tid), chunk_sz(tid))
    end do

    deallocate(chunk_cs, chunk_sz)
  end function crc64_omp

end module crc64_mod

program main
  use iso_fortran_env, only: int32, int64, real64
  use iso_c_binding, only: c_int, c_long
  use crc64_mod
  implicit none

  type, bind(C) :: timespec
    integer(c_long) :: tv_sec
    integer(c_long) :: tv_nsec
  end type timespec

  interface
    function clock_gettime(clk_id, tp) bind(C, name="clock_gettime") result(rc)
      import :: c_int, timespec
      integer(c_int), value :: clk_id
      type(timespec), intent(out) :: tp
      integer(c_int) :: rc
    end function clock_gettime
  end interface

  integer(c_int), parameter :: CLOCK_THREAD_CPUTIME_ID = 3_c_int

  integer :: argc, ntests, seed, max_test_length, ntest
  integer(int64) :: rng_state, test_length, div_pt, tlend
  integer(int32), allocatable :: buffer(:)
  integer(int32) :: check_bytes(8)
  integer(int64) :: cs, csc, cs1, cs2
  real(real64) :: rnd, start_time, end_time, b_time, tot_time, tot_bytes
  character(len=64) :: arg
  character(len=4) :: check1, check2
  integer(int64) :: i
  type(timespec) :: b_start, b_end
  integer(c_int) :: clock_status

  ntests = 10
  seed = 5
  max_test_length = 2097152

  argc = command_argument_count()
  if (argc > 0) then
    call get_command_argument(1, arg)
    read(arg, *) ntests
  end if
  if (argc > 1) then
    call get_command_argument(2, arg)
    read(arg, *) seed
  end if
  if (argc > 2) then
    call get_command_argument(3, arg)
    read(arg, *) max_test_length
  end if

  write(*,'("Running ",I0," tests with seed ",I0)') ntests, seed

  call srand48_state(int(seed, int32), rng_state)
  tot_time = 0.0_real64
  tot_bytes = 0.0_real64
  tlend = 8_int64

  do ntest = 1, ntests
    call next_drand48(rng_state, rnd)
    test_length = int(real(max_test_length, real64) * (rnd + 1.0_real64), int64)
    allocate(buffer(test_length + tlend))

    do i = 1, test_length
      call next_drand48(rng_state, rnd)
      buffer(i) = int(255.0_real64 * rnd, int32)
    end do

    clock_status = clock_gettime(CLOCK_THREAD_CPUTIME_ID, b_start)
    cs = crc64_omp(buffer, test_length)
    clock_status = clock_gettime(CLOCK_THREAD_CPUTIME_ID, b_end)
    start_time = real(b_start%tv_sec, real64) + 1.0e-9_real64 * real(b_start%tv_nsec, real64)
    end_time = real(b_end%tv_sec, real64) + 1.0e-9_real64 * real(b_end%tv_nsec, real64)
    b_time = end_time - start_time

    if (ntest > 1) then
      tot_time = tot_time + b_time
      tot_bytes = tot_bytes + real(test_length, real64)
    end if

    buffer(test_length + 1:test_length + tlend) = 0_int32
    call crc64_invert(cs, check_bytes)
    buffer(test_length + 1:test_length + tlend) = check_bytes

    csc = crc64(buffer, 1_int64, test_length + tlend, crc64_table, crc64_interleaved_table)
    if (csc == all_ones) then
      check1 = "pass"
    else
      check1 = "fail"
    end if

    call next_drand48(rng_state, rnd)
    div_pt = int(real(test_length, real64) * rnd, int64)
    if (div_pt > 0_int64) then
      cs1 = crc64(buffer, 1_int64, div_pt, crc64_table, crc64_interleaved_table)
    else
      cs1 = crc64(buffer, 1_int64, 0_int64, crc64_table, crc64_interleaved_table)
    end if
    cs2 = crc64(buffer, div_pt + 1_int64, test_length - div_pt, crc64_table, crc64_interleaved_table)
    csc = crc64_combine(cs1, cs2, test_length - div_pt)
    if (csc == cs) then
      check2 = "pass"
    else
      check2 = "fail"
    end if

    write(*,'(I0,1X,I0,1X,A,1X,A)') ntest, test_length, trim(check1), trim(check2)
    deallocate(buffer)
  end do

  if (tot_time > 0.0_real64) then
    write(*,'(G0.6,1X,"MB/s")') (tot_bytes / (1024.0_real64 * 1024.0_real64)) / tot_time
  else
    write(*,'("inf MB/s")')
  end if

end program main
