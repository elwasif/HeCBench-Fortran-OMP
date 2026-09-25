! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer, parameter :: p_x = 5, p_y = 1, modulus = 17, curve_a = 2
  character(len=256) :: arg0, arg1, arg2
  integer :: num_pk, repeat, iter
  integer, allocatable :: pk_slow_x(:), pk_slow_y(:), pk_fast_x(:), pk_fast_y(:)
  real(real64) :: start_time, end_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <positive number of keys> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) num_pk
  read(arg2, *) repeat
  if (num_pk <= 0 .or. repeat <= 0) stop 1

  allocate(pk_slow_x(num_pk), pk_slow_y(num_pk), pk_fast_x(num_pk), pk_fast_y(num_pk))

  !$omp target data map(from: pk_slow_x(1:num_pk), pk_slow_y(1:num_pk), &
  !$omp& pk_fast_x(1:num_pk), pk_fast_y(1:num_pk))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call k_slow(18, p_x, p_y, pk_slow_x, pk_slow_y, modulus, curve_a, num_pk)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average time (slow kernel): ', &
    (end_time - start_time) / real(repeat, real64), ' s'

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call k_fast(18, p_x, p_y, pk_fast_x, pk_fast_y, modulus, curve_a, num_pk)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average time (fast kernel): ', &
    (end_time - start_time) / real(repeat, real64), ' s'
  !$omp end target data

  if (all(pk_slow_x == pk_fast_x) .and. all(pk_slow_y == pk_fast_y)) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(pk_slow_x, pk_slow_y, pk_fast_x, pk_fast_y)

contains

  subroutine k_slow(sk, px, py, tx, ty, m, ca, num_pk)
    integer, intent(in) :: sk, px, py, m, ca, num_pk
    integer, intent(out) :: tx(:), ty(:)
    integer :: i

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, num_pk
      call make_pk_slow(sk, px, py, tx(i), ty(i), m, ca)
    end do
    !$omp end target teams distribute parallel do
  end subroutine k_slow

  subroutine k_fast(sk, px, py, tx, ty, m, ca, num_pk)
    integer, intent(in) :: sk, px, py, m, ca, num_pk
    integer, intent(out) :: tx(:), ty(:)
    integer :: i

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, num_pk
      call make_pk_fast(sk, px, py, tx(i), ty(i), m, ca)
    end do
    !$omp end target teams distribute parallel do
  end subroutine k_fast

  integer function ext_euclidian_alg(a_in, b_in, x, y)
    integer, intent(in) :: a_in, b_in
    integer, intent(out) :: x, y
    integer :: x1, y1, a1, b1, s, t, q

    x = 1
    y = 0
    x1 = 0
    y1 = 1
    a1 = a_in
    b1 = b_in
    do while (b1 /= 0)
      q = a1 / b1
      s = x1
      t = x - q * x1
      x = s
      x1 = t
      s = y1
      t = y - q * y1
      y = s
      y1 = t
      s = b1
      t = a1 - q * b1
      a1 = s
      b1 = t
    end do
    ext_euclidian_alg = a1
  end function ext_euclidian_alg

  integer function make_positive(a_in, m)
    integer, intent(in) :: a_in, m
    integer :: a

    a = a_in
    do while (a < 0)
      a = a + m
    end do
    make_positive = mod(a, m)
  end function make_positive

  integer function find_inverse(a, m)
    integer, intent(in) :: a, m
    integer :: t, s, ignored

    ignored = ext_euclidian_alg(a, m, t, s)
    find_inverse = make_positive(t, m)
  end function find_inverse

  subroutine point_addition(m, x1, y1, x2, y2, x3, y3)
    integer, intent(in) :: m, x1, y1, x2, y2
    integer, intent(out) :: x3, y3
    integer :: temp, slope

    temp = make_positive(x2 - x1, m)
    slope = make_positive((y2 - y1) * find_inverse(temp, m), m)
    x3 = make_positive(slope * slope - x1 - x2, m)
    y3 = make_positive(slope * (x1 - x3) - y1, m)
  end subroutine point_addition

  subroutine point_doubling(m, ca, x1, y1, x3, y3)
    integer, intent(in) :: m, ca, x1, y1
    integer, intent(out) :: x3, y3
    integer :: slope

    slope = (3 * x1 * x1 + ca) * find_inverse(2 * y1, m)
    x3 = make_positive(slope * slope - 2 * x1, m)
    y3 = make_positive(slope * (x1 - x3) - y1, m)
  end subroutine point_doubling

  integer function first_set_bit(n)
    integer, intent(in) :: n
    integer :: i

    do i = bit_size(n) - 1, 0, -1
      if (btest(n, i)) then
        first_set_bit = i
        return
      end if
    end do
    first_set_bit = 0
  end function first_set_bit

  subroutine make_pk_fast(sk, px, py, tx, ty, m, ca)
    integer, intent(in) :: sk, px, py, m, ca
    integer, intent(out) :: tx, ty
    integer :: i, cur_x, cur_y, next_x, next_y

    cur_x = px
    cur_y = py
    do i = first_set_bit(sk) - 1, 0, -1
      call point_doubling(m, ca, cur_x, cur_y, next_x, next_y)
      cur_x = next_x
      cur_y = next_y
      if (btest(sk, i)) then
        call point_addition(m, cur_x, cur_y, px, py, next_x, next_y)
        cur_x = next_x
        cur_y = next_y
      end if
    end do
    tx = cur_x
    ty = cur_y
  end subroutine make_pk_fast

  subroutine make_pk_slow(sk, px, py, tx, ty, m, ca)
    integer, intent(in) :: sk, px, py, m, ca
    integer, intent(out) :: tx, ty
    integer :: remaining, cur_x, cur_y, next_x, next_y

    call point_doubling(m, ca, px, py, cur_x, cur_y)
    remaining = sk - 2
    do while (remaining > 0)
      call point_addition(m, cur_x, cur_y, px, py, next_x, next_y)
      cur_x = next_x
      cur_y = next_y
      remaining = remaining - 1
    end do
    tx = cur_x
    ty = cur_y
  end subroutine make_pk_slow

end program main
