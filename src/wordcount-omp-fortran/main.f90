! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  character(len=*), parameter :: tab = "abcdefghijklmnopqrstuvwxyz ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: raw_input = &
    "  But the raven, sitting lonely on the placid bust, spoke only," // new_line('a') // &
    "  That one word, as if his soul in that one word he did outpour." // new_line('a') // &
    "  Nothing further then he uttered - not a feather then he fluttered -" // new_line('a') // &
    "  Till I scarcely more than muttered `Other friends have flown before -" // new_line('a') // &
    "  On the morrow he will leave me, as my hopes have flown before.'" // new_line('a') // &
    "  Then the bird said, `Nevermore.'" // new_line('a')

  character(len=256) :: arg0, arg1
  character(len=1), allocatable :: input(:), random_input(:)
  integer(int64) :: wc_host, wc_device
  integer(int32) :: repeat, i
  integer(int64) :: len_large, n
  integer(int64) :: start_ns, end_ns, elapsed_ns
  integer(int64), save :: rng_state = 1_int64
  real(real64) :: avg_s
  character(len=32) :: avg_text
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat

  print '(A)', 'Text sample:'
  write(*, '(A)') raw_input

  call string_to_char_array(raw_input, input)
  wc_host = word_count_reference(input)
  print '(A,I0,A)', 'Host: Text sample contains ', wc_host, ' words'

  wc_device = word_count(input)
  print '(A,I0,A)', 'Device: Text sample contains ', wc_device, ' words'

  print '(A)', 'Test word count with random inputs'
  call rng_seed(123_int64)
  ok = .true.

  n = 1_int64
  do while (n <= 100000000_int64)
    allocate(random_input(n))
    call fill_random_text(random_input)
    if (word_count_reference(random_input) /= word_count(random_input)) then
      ok = .false.
      deallocate(random_input)
      exit
    end if
    deallocate(random_input)
    n = n * 10_int64
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  len_large = 1024_int64 * 1024_int64 * 256_int64
  allocate(random_input(len_large))
  call fill_random_text(random_input)

  print '(A,I0)', 'Performance evaluation for random texts of character length ', len_large
  start_ns = int(omp_get_wtime() * 1.0d9, int64)
  do i = 0, repeat - 1
    wc_device = word_count(random_input)
  end do
  end_ns = int(omp_get_wtime() * 1.0d9, int64)
  elapsed_ns = end_ns - start_ns
  avg_s = real(elapsed_ns, real64) * 1.0d-9 / real(repeat, real64)
  avg_text = format_default_float(avg_s)
  print '(3A)', 'Average time of word count: ', trim(avg_text), ' (s)'

  deallocate(input)
  deallocate(random_input)

contains

  integer(int64) function rng_next() result(v)
    rng_state = mod(1103515245_int64 * rng_state + 12345_int64, 2147483648_int64)
    v = rng_state
  end function rng_next

  subroutine rng_seed(seed)
    integer(int64), intent(in) :: seed
    rng_state = mod(seed, 2147483648_int64)
    if (rng_state < 0_int64) rng_state = rng_state + 2147483648_int64
    if (rng_state == 0_int64) rng_state = 1_int64
  end subroutine rng_seed

  character(len=32) function format_default_float(value) result(text)
    real(real64), intent(in) :: value
    text = ''
    write(text, '(F0.7)') value
    if (text(1:1) == '.') text = '0' // text(1:len(text)-1)
  end function format_default_float

  logical function is_alpha(c) result(res)
    character(len=1), intent(in) :: c
    integer :: code
    code = iachar(c)
    res = (code >= iachar('A') .and. code <= iachar('z'))
  end function is_alpha

  integer(int64) function word_count_reference(chars) result(wc)
    character(len=1), intent(in) :: chars(:)
    integer(int64) :: i, sz
    sz = int(size(chars), int64)
    if (sz == 0_int64) then
      wc = 0_int64
      return
    end if
    wc = 0_int64
    do i = 1_int64, sz - 1_int64
      if ((.not. is_alpha(chars(i))) .and. is_alpha(chars(i + 1_int64))) wc = wc + 1_int64
    end do
    if (is_alpha(chars(1))) wc = wc + 1_int64
  end function word_count_reference

  integer(int64) function word_count(chars) result(wc)
    character(len=1), intent(in) :: chars(:)
    integer(int64) :: i, sz
    integer :: left_code, right_code
    sz = int(size(chars), int64)
    if (sz == 0_int64) then
      wc = 0_int64
      return
    end if

    wc = 0_int64
    !$omp target data map(to: chars(1:sz))
    !$omp target teams distribute parallel do thread_limit(256) reduction(+:wc)
    do i = 1_int64, sz - 1_int64
      left_code = iachar(chars(i))
      right_code = iachar(chars(i + 1_int64))
      if ((left_code < iachar('A') .or. left_code > iachar('z')) .and. &
          (right_code >= iachar('A') .and. right_code <= iachar('z'))) then
        wc = wc + 1_int64
      end if
    end do
    !$omp end target teams distribute parallel do
    !$omp end target data

    if (is_alpha(chars(1))) wc = wc + 1_int64
  end function word_count

  subroutine fill_random_text(chars)
    character(len=1), intent(out) :: chars(:)
    integer(int64) :: i, idx, tlen
    tlen = int(len(tab), int64)
    do i = 1_int64, int(size(chars), int64)
      idx = mod(rng_next(), tlen) + 1_int64
      chars(i) = tab(idx:idx)
    end do
  end subroutine fill_random_text

  subroutine string_to_char_array(s, chars)
    character(len=*), intent(in) :: s
    character(len=1), allocatable, intent(out) :: chars(:)
    integer :: i, nchar
    nchar = len(s)
    allocate(chars(nchar))
    do i = 1, nchar
      chars(i) = s(i:i)
    end do
  end subroutine string_to_char_array

end program main
