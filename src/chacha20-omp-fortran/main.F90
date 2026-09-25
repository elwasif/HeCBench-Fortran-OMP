! SPDX-License-Identifier: CC0-1.0
#define ROTL32(V,N) ior(ishft((V), (N)), shiftr((V), 32 - (N)))
#define QR(A,B,C,D) state(A)=add32(state(A),state(B)); tmp=ieor(state(D),state(A)); state(D)=ROTL32(tmp,16); state(C)=add32(state(C),state(D)); tmp=ieor(state(B),state(C)); state(B)=ROTL32(tmp,12); state(A)=add32(state(A),state(B)); tmp=ieor(state(D),state(A)); state(D)=ROTL32(tmp,8); state(C)=add32(state(C),state(D)); tmp=ieor(state(B),state(C)); state(B)=ROTL32(tmp,7)

program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: key_len = 64
  integer, parameter :: nonce_len = 16
  integer, parameter :: keystream_len = 512
  integer, parameter :: result_len = keystream_len / 2
  integer(int64), parameter :: mask32 = int(z'ffffffff', int64)

  character(len=*), parameter :: text_key = &
       '000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f'
  character(len=*), parameter :: text_nonce = '0001020304050607'
  character(len=*), parameter :: text_keystream = &
       'f798a189f195e66982105ffb640bb7757f579da31602fc93ec01ac56f85ac3c1' // &
       '34a4547b733b46413042c9440049176905d3be59ea1c53f15916155c2be8241a' // &
       '38008b9a26bc35941e2444177c8ade6689de95264986d95889fb60e84629c9bd' // &
       '9a5acb1cc118be563eb9b3a4a472f82e09a7e778492b562ef7130e88dfe031' // &
       'c79db9d4f7c7a899151b9a475032b63fc385245fe054e3dd5a97a5f576fe064' // &
       '025d3ce042c566ab2c507b138db853e3d6959660996546cc9c4a6eafdc777c' // &
       '040d70eaf46f76dad3979e5c5360c3317166a1c894c94a371876a94df7628fe' // &
       '4eaaf2ccb27d5aaae0ad7ad0f9d4b6ad3b54098746d4524d38407a6deb3ab78fab78c9'

  character(len=256) :: arg0, arg
  integer :: repeat, iter, i, j, idx, block, word, byte_pos
  integer(int32) :: tmp
  integer(int64) :: counter
  integer(int8) :: hi, lo
  integer(int8) :: char_to_uint(0:255)
  integer(int8) :: key_chars(key_len), nonce_chars(nonce_len), keystream_chars(keystream_len)
  integer(int8) :: raw_key(key_len / 2), raw_nonce(nonce_len / 2)
  integer(int8) :: raw_keystream(result_len), result(result_len)
  integer(int32) :: state0(16), state(16), out_state(16)
  integer(int8) :: keystream8(64)
  real(real64) :: start_time, elapsed_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) repeat

  char_to_uint = 0_int8
  do i = 0, 9
    char_to_uint(iachar('0') + i) = int(i, int8)
  end do
  do i = 0, 25
    char_to_uint(iachar('a') + i) = int(i + 10, int8)
    char_to_uint(iachar('A') + i) = int(i + 10, int8)
  end do

  do i = 1, key_len
    key_chars(i) = int(iachar(text_key(i:i)), int8)
  end do
  do i = 1, nonce_len
    nonce_chars(i) = int(iachar(text_nonce(i:i)), int8)
  end do
  do i = 1, keystream_len
    keystream_chars(i) = int(iachar(text_keystream(i:i)), int8)
  end do

  raw_key = 0_int8
  raw_nonce = 0_int8
  raw_keystream = 0_int8
  result = 0_int8

  !$omp target data map(to: char_to_uint, key_chars, nonce_chars, keystream_chars) &
  !$omp& map(alloc: raw_key, raw_nonce) map(from: raw_keystream, result)
  start_time = omp_get_wtime()
  do iter = 1, repeat
    !$omp target teams distribute parallel do
    do i = 1, result_len
      result(i) = 0_int8
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams num_teams(1) thread_limit(256) private(i, j, idx, hi, lo, tmp, counter, block, word, byte_pos, state0, state, out_state, keystream8)
    !$omp parallel private(i, j, idx, hi, lo, tmp, counter, block, word, byte_pos, state0, state, out_state, keystream8)
    !$omp do
    do i = 1, key_len / 2
      hi = char_to_uint(key_chars(2 * i - 1))
      lo = char_to_uint(key_chars(2 * i))
      raw_key(i) = int(ior(ishft(int(hi, int32), 4), int(lo, int32)), int8)
    end do
    !$omp end do
    !$omp do
    do i = 1, nonce_len / 2
      hi = char_to_uint(nonce_chars(2 * i - 1))
      lo = char_to_uint(nonce_chars(2 * i))
      raw_nonce(i) = int(ior(ishft(int(hi, int32), 4), int(lo, int32)), int8)
    end do
    !$omp end do
    !$omp do
    do i = 1, result_len
      hi = char_to_uint(keystream_chars(2 * i - 1))
      lo = char_to_uint(keystream_chars(2 * i))
      raw_keystream(i) = int(ior(ishft(int(hi, int32), 4), int(lo, int32)), int8)
    end do
    !$omp end do

    !$omp single
    state0 = 0_int32
    state0(1) = int(z'61707865', int32)
    state0(2) = int(z'3320646e', int32)
    state0(3) = int(z'79622d32', int32)
    state0(4) = int(z'6b206574', int32)
    do i = 1, 8
      idx = 4 * (i - 1) + 1
      state0(4 + i) = pack4(raw_key(idx))
    end do
    state0(13) = 0_int32
    state0(14) = 0_int32
    state0(15) = pack4(raw_nonce(1))
    state0(16) = pack4(raw_nonce(5))

    counter = 0_int64
    byte_pos = 1
    do block = 1, (result_len + 63) / 64
      state = state0
      state(13) = u32(counter)
      state(14) = u32(ishft(counter, -32))
      out_state = state
      do j = 1, 10
        QR(1, 5, 9, 13)
        QR(2, 6, 10, 14)
        QR(3, 7, 11, 15)
        QR(4, 8, 12, 16)
        QR(1, 6, 11, 16)
        QR(2, 7, 12, 13)
        QR(3, 8, 9, 14)
        QR(4, 5, 10, 15)
      end do
      do word = 1, 16
        out_state(word) = add32(out_state(word), state(word))
        keystream8(4 * word - 3) = int(iand(out_state(word), int(z'ff', int32)), int8)
        keystream8(4 * word - 2) = int(iand(shiftr(out_state(word), 8), int(z'ff', int32)), int8)
        keystream8(4 * word - 1) = int(iand(shiftr(out_state(word), 16), int(z'ff', int32)), int8)
        keystream8(4 * word) = int(iand(shiftr(out_state(word), 24), int(z'ff', int32)), int8)
      end do
      do i = 1, 64
        if (byte_pos <= result_len) then
          result(byte_pos) = ieor(result(byte_pos), keystream8(i))
          byte_pos = byte_pos + 1
        end if
      end do
      counter = counter + 1_int64
    end do
    !$omp end single
    !$omp end parallel
    !$omp end target teams
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average execution time of kernels: ', elapsed_us, ' (us)'
  !$omp end target data

  ok = .true.
  do i = 1, result_len
    if (result(i) /= raw_keystream(i)) then
      ok = .false.
      exit
    end if
  end do
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
    stop 1
  end if

contains

  elemental integer(int32) function u32(value)
    integer(int64), intent(in) :: value
    u32 = transfer(iand(value, mask32), u32)
  end function u32

  elemental integer(int32) function add32(a, b)
    integer(int32), intent(in) :: a, b
    add32 = u32(int(a, int64) + int(b, int64))
  end function add32

  integer(int32) function pack4(a)
    integer(int8), intent(in) :: a(4)
    pack4 = ior(ior(iand(int(a(1), int32), int(z'ff', int32)), &
                    ishft(iand(int(a(2), int32), int(z'ff', int32)), 8)), &
                ior(ishft(iand(int(a(3), int32), int(z'ff', int32)), 16), &
                    ishft(iand(int(a(4), int32), int(z'ff', int32)), 24)))
  end function pack4
end program main
