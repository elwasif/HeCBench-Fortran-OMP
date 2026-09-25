! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use, intrinsic :: iso_c_binding, only : c_int, c_int8_t
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

  integer, parameter :: block_size = 256
  character(len=256) :: arg0
  integer :: num_keys, repeat, total_length
  integer :: i, c
  integer(int32), allocatable :: length(:), offsets(:)
  integer(c_int8_t), allocatable :: keys(:)
  integer(int64), allocatable :: ref_out(:), device_out(:)
  real(real64) :: start_time, end_time, avg_time
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of keys> <repeat>'
    stop 1
  end if

  num_keys = read_arg(1)
  repeat = read_arg(2)
  if (num_keys <= 0 .or. repeat <= 0) stop 1

  allocate(length(num_keys), offsets(num_keys + 1))
  offsets(1) = 0_int32
  total_length = 0
  call c_srand(3_c_int)
  do i = 1, num_keys
    length(i) = int(mod(c_rand(), 10000_c_int), int32)
    total_length = total_length + int(length(i))
    offsets(i + 1) = int(total_length, int32)
  end do

  allocate(keys(max(total_length, 1)))
  do i = 1, num_keys
    do c = 0, int(length(i)) - 1
      keys(int(offsets(i)) + c + 1) = uint8_storage(c)
    end do
  end do

  allocate(ref_out(2 * num_keys), device_out(2 * num_keys))
  do i = 1, num_keys
    call murmurhash3_x64_128(keys, int(offsets(i)) + 1, int(length(i)), i - 1, ref_out(2 * i - 1), ref_out(2 * i))
  end do

  device_out = 0_int64

  !$omp target data map(to: keys(1:max(total_length, 1)), offsets(1:num_keys+1), length(1:num_keys)) &
  !$omp& map(from: device_out(1:2*num_keys))
  start_time = omp_get_wtime()
  do c = 1, repeat
    call murmurhash_device(num_keys, keys, offsets, length, device_out)
  end do
  end_time = omp_get_wtime()
  avg_time = (end_time - start_time) / real(repeat, real64)
  print '(A,F8.6,A)', 'Average kernel execution time ', avg_time, ' (s)'
  !$omp end target data

  ok = .true.
  do i = 1, 2 * num_keys
    if (device_out(i) /= ref_out(i)) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'SUCCESS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(length, offsets, keys, ref_out, device_out)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine murmurhash_device(num_keys, keys, offsets, length, out)
    integer, intent(in) :: num_keys
    integer(c_int8_t), intent(in) :: keys(:)
    integer(int32), intent(in) :: offsets(:), length(:)
    integer(int64), intent(out) :: out(:)
    integer :: i
    integer(int64) :: h1, h2

    !$omp target teams distribute parallel do thread_limit(block_size) private(h1, h2)
    do i = 1, num_keys
      call murmurhash3_x64_128(keys, int(offsets(i)) + 1, int(length(i)), i - 1, h1, h2)
      out(2 * i - 1) = h1
      out(2 * i) = h2
    end do
    !$omp end target teams distribute parallel do
  end subroutine murmurhash_device

  subroutine murmurhash3_x64_128(data, start_pos, len, seed, out1, out2)
    integer(c_int8_t), intent(in) :: data(:)
    integer, intent(in) :: start_pos, len, seed
    integer(int64), intent(out) :: out1, out2
    integer(int64), parameter :: c1 = int(z'87c37b91114253d5', int64)
    integer(int64), parameter :: c2 = int(z'4cf5ad432745937f', int64)
    integer :: nblocks, block, tail_pos, tail_len
    integer(int64) :: h1, h2, k1, k2

    nblocks = len / 16
    h1 = int(seed, int64)
    h2 = int(seed, int64)

    do block = 0, nblocks - 1
      k1 = getblock64(data, start_pos + 16 * block)
      k2 = getblock64(data, start_pos + 16 * block + 8)

      k1 = k1 * c1
      k1 = rotl64(k1, 31)
      k1 = k1 * c2
      h1 = ieor(h1, k1)

      h1 = rotl64(h1, 27)
      h1 = h1 + h2
      h1 = h1 * 5_int64 + int(z'52dce729', int64)

      k2 = k2 * c2
      k2 = rotl64(k2, 33)
      k2 = k2 * c1
      h2 = ieor(h2, k2)

      h2 = rotl64(h2, 31)
      h2 = h2 + h1
      h2 = h2 * 5_int64 + int(z'38495ab5', int64)
    end do

    tail_pos = start_pos + nblocks * 16
    tail_len = iand(len, 15)
    k1 = 0_int64
    k2 = 0_int64

    if (tail_len >= 15) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 14)), 48))
    if (tail_len >= 14) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 13)), 40))
    if (tail_len >= 13) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 12)), 32))
    if (tail_len >= 12) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 11)), 24))
    if (tail_len >= 11) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 10)), 16))
    if (tail_len >= 10) k2 = ieor(k2, shiftl(byte64(data(tail_pos + 9)), 8))
    if (tail_len >= 9) then
      k2 = ieor(k2, byte64(data(tail_pos + 8)))
      k2 = k2 * c2
      k2 = rotl64(k2, 33)
      k2 = k2 * c1
      h2 = ieor(h2, k2)
    end if

    if (tail_len >= 8) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 7)), 56))
    if (tail_len >= 7) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 6)), 48))
    if (tail_len >= 6) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 5)), 40))
    if (tail_len >= 5) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 4)), 32))
    if (tail_len >= 4) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 3)), 24))
    if (tail_len >= 3) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 2)), 16))
    if (tail_len >= 2) k1 = ieor(k1, shiftl(byte64(data(tail_pos + 1)), 8))
    if (tail_len >= 1) then
      k1 = ieor(k1, byte64(data(tail_pos)))
      k1 = k1 * c1
      k1 = rotl64(k1, 31)
      k1 = k1 * c2
      h1 = ieor(h1, k1)
    end if

    h1 = ieor(h1, int(len, int64))
    h2 = ieor(h2, int(len, int64))

    h1 = h1 + h2
    h2 = h2 + h1

    h1 = fmix64(h1)
    h2 = fmix64(h2)

    h1 = h1 + h2
    h2 = h2 + h1

    out1 = h1
    out2 = h2
  end subroutine murmurhash3_x64_128

  integer(int64) function getblock64(data, start_pos)
    integer(c_int8_t), intent(in) :: data(:)
    integer, intent(in) :: start_pos
    integer :: n

    getblock64 = 0_int64
    do n = 0, 7
      getblock64 = ior(getblock64, shiftl(byte64(data(start_pos + n)), n * 8))
    end do
  end function getblock64

  integer(int64) function fmix64(k_in)
    integer(int64), intent(in) :: k_in
    integer(int64), parameter :: c1 = int(z'ff51afd7ed558ccd', int64)
    integer(int64), parameter :: c2 = int(z'c4ceb9fe1a85ec53', int64)

    fmix64 = k_in
    fmix64 = ieor(fmix64, shiftr(fmix64, 33))
    fmix64 = fmix64 * c1
    fmix64 = ieor(fmix64, shiftr(fmix64, 33))
    fmix64 = fmix64 * c2
    fmix64 = ieor(fmix64, shiftr(fmix64, 33))
  end function fmix64

  integer(int64) function rotl64(x, r)
    integer(int64), intent(in) :: x
    integer, intent(in) :: r

    rotl64 = ior(shiftl(x, r), shiftr(x, 64 - r))
  end function rotl64

  integer(int64) function byte64(value)
    integer(c_int8_t), intent(in) :: value

    byte64 = int(iand(int(value, int32), int(z'000000ff', int32)), int64)
  end function byte64

  integer(c_int8_t) function uint8_storage(value)
    integer, intent(in) :: value
    integer :: byte_value

    byte_value = iand(value, 255)
    if (byte_value < 128) then
      uint8_storage = int(byte_value, c_int8_t)
    else
      uint8_storage = int(byte_value - 256, c_int8_t)
    end if
  end function uint8_storage

end program main
