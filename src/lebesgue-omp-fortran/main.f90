! SPDX-License-Identifier: CC0-1.0
module lebesgue_mod
  use omp_lib
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793_dp

contains

  subroutine generate_points(kind_id, n, x)
    integer, intent(in) :: kind_id, n
    real(dp), intent(out) :: x(n)
    integer :: i
    real(dp) :: angle

    select case (kind_id)
    case (1)
      do i = 1, n
        angle = pi * real(2 * i - 1, dp) / real(2 * n, dp)
        x(i) = cos(angle)
      end do
    case (2)
      if (n == 1) then
        x(1) = 0.0_dp
      else
        do i = 1, n
          angle = pi * real(n - i, dp) / real(n - 1, dp)
          x(i) = cos(angle)
        end do
      end if
    case (3)
      do i = 1, n
        angle = pi * real(2 * n - 2 * i + 1, dp) / real(2 * n + 1, dp)
        x(i) = cos(angle)
      end do
    case (4)
      do i = 1, n
        angle = pi * real(2 * n - 2 * i + 2, dp) / real(2 * n + 1, dp)
        x(i) = cos(angle)
      end do
    case (5)
      do i = 1, n
        x(i) = real(-n + 1 + 2 * (i - 1), dp) / real(n + 1, dp)
      end do
    case (6)
      if (n == 1) then
        x(1) = 0.0_dp
      else
        do i = 1, n
          x(i) = real(-n + 1 + 2 * (i - 1), dp) / real(n - 1, dp)
        end do
      end if
    case (7)
      do i = 1, n
        x(i) = real(-n + 1 + 2 * (i - 1), dp) / real(n, dp)
      end do
    case (8)
      do i = 1, n
        angle = pi * real(2 * n - 2 * i + 1, dp) / real(2 * n, dp)
        x(i) = cos(angle)
      end do
    case (9)
      do i = 1, n
        angle = pi * real(n - i + 1, dp) / real(n + 1, dp)
        x(i) = cos(angle)
      end do
    end select
  end subroutine generate_points

  subroutine r8vec_linspace_new(n, a, b, x)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x(n)
    integer :: i

    if (n == 1) then
      x(1) = 0.5_dp * (a + b)
    else
      do i = 1, n
        x(i) = (real(n - i, dp) * a + real(i - 1, dp) * b) / real(n - 1, dp)
      end do
    end if
  end subroutine r8vec_linspace_new

  function lebesgue_function(n, x, nfun, xfun) result(lmax)
    integer, intent(in) :: n, nfun
    real(dp), intent(in) :: x(n), xfun(nfun)
    real(dp) :: lmax
    integer :: j, i1, i2
    real(dp) :: t
    real(dp), allocatable :: linterp(:)

    lmax = 0.0_dp
    allocate(linterp(n * nfun))

    !$omp target data map(tofrom:lmax) &
    !$omp& map(to:x(1:n), xfun(1:nfun)) &
    !$omp& map(alloc:linterp(1:n * nfun))
    !$omp target teams distribute parallel do thread_limit(256) reduction(max:lmax) private(i1, i2, t)
    do j = 1, nfun
      t = 0.0_dp
      do i1 = 1, n
        linterp((i1 - 1) * nfun + j) = 1.0_dp
        do i2 = 1, n
          if (i1 /= i2) then
            linterp((i1 - 1) * nfun + j) = linterp((i1 - 1) * nfun + j) * &
              (xfun(j) - x(i2)) / (x(i1) - x(i2))
          end if
        end do
        t = t + abs(linterp((i1 - 1) * nfun + j))
      end do
      lmax = max(lmax, t)
    end do
    !$omp end target teams distribute parallel do
    !$omp end target data

    deallocate(linterp)
  end function lebesgue_function

  function lebesgue_constant(n, x, nfun, xfun) result(lmax)
    integer, intent(in) :: n, nfun
    real(dp), intent(in) :: x(n), xfun(nfun)
    real(dp) :: lmax

    if (n > 1) then
      lmax = lebesgue_function(n, x, nfun, xfun)
    else
      lmax = 1.0_dp
    end if
  end function lebesgue_constant

  logical function verify_result(res, n, x, nfun, xfun)
    integer, intent(in) :: n, nfun
    real(dp), intent(in) :: res, x(n), xfun(nfun)
    integer :: j, i1, i2
    real(dp) :: lmax, t, value

    lmax = 0.0_dp
    do j = 1, nfun
      t = 0.0_dp
      do i1 = 1, n
        value = 1.0_dp
        do i2 = 1, n
          if (i1 /= i2) then
            value = value * (xfun(j) - x(i2)) / (x(i1) - x(i2))
          end if
        end do
        t = t + abs(value)
      end do
      lmax = max(lmax, t)
    end do

    verify_result = abs(res - lmax) <= 1.0d-6
  end function verify_result

  subroutine timestamp()
    integer :: values(8)
    integer :: hour12
    character(len=9), dimension(12), parameter :: months = &
      [character(len=9) :: 'January  ', 'February ', 'March    ', 'April    ', 'May      ', 'June     ', &
                            'July     ', 'August   ', 'September', 'October  ', 'November ', 'December ']
    character(len=2) :: ampm

    call date_and_time(values=values)
    hour12 = modulo(values(5), 12)
    if (hour12 == 0) hour12 = 12
    if (values(5) < 12) then
      ampm = 'AM'
    else
      ampm = 'PM'
    end if

    write(*,'(i2.2,1x,a,1x,i4,1x,i2.2,a1,i2.2,a1,i2.2,1x,a)') &
      values(3), trim(months(values(2))), values(1), hour12, ':', values(6), ':', values(7), ampm
  end subroutine timestamp

  subroutine run_test(test_id, heading, description, nfun)
    integer, intent(in) :: test_id, nfun
    character(len=*), intent(in) :: heading, description
    integer, parameter :: n_max = 11
    integer :: n
    logical :: ok
    real(dp) :: total_time, t_start, t_end
    real(dp), allocatable :: l(:), x(:), xfun(:)

    allocate(xfun(nfun), l(n_max))
    call r8vec_linspace_new(nfun, -1.0_dp, 1.0_dp, xfun)

    print *
    print '(a)', trim(heading)
    print '(a)', '  ' // trim(description)

    total_time = 0.0_dp
    ok = .true.

    do n = 1, n_max
      allocate(x(n))
      call generate_points(test_id, n, x)

      t_start = omp_get_wtime()
      l(n) = lebesgue_constant(n, x, nfun, xfun)
      t_end = omp_get_wtime()
      total_time = total_time + (t_end - t_start)

      ok = ok .and. verify_result(l(n), n, x, nfun, xfun)
      deallocate(x)
    end do

    write(*,'(a,f10.6,a)') '  Total kernel execution time ', total_time, ' (s)'
    if (ok) then
      print '(a)', '  PASS'
    else
      print '(a)', '  FAIL'
    end if

    deallocate(l, xfun)
  end subroutine run_test

end module lebesgue_mod

program main
  use lebesgue_mod
  implicit none

  integer :: argc, ios, nfun, repeat, i
  character(len=64) :: arg

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, arg)
    write(*,'(a,a,a)') 'Usage: ', trim(arg), ' <number of points in an interval> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=ios) nfun
  if (ios /= 0) stop 1

  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) repeat
  if (ios /= 0) stop 1

  print *
  print '(a)', 'LEBESGUE_TEST'

  do i = 1, repeat
    call timestamp()
    call run_test(1, 'LEBESGUE_TEST01:', 'Analyze Chebyshev1 points.', nfun)
    call run_test(2, 'LEBESGUE_TEST02:', 'Analyze Chebyshev2 points.', nfun)
    call run_test(3, 'LEBESGUE_TEST03:', 'Analyze Chebyshev3 points.', nfun)
    call run_test(4, 'LEBESGUE_TEST04:', 'Analyze Chebyshev4 points.', nfun)
    call run_test(5, 'LEBESGUE_TEST05:', 'Analyze Equidistant1 points.', nfun)
    call run_test(6, 'LEBESGUE_TEST06:', 'Analyze Equidistant2 points.', nfun)
    call run_test(7, 'LEBESGUE_TEST07:', 'Analyze Equidistant3 points.', nfun)
    call run_test(8, 'LEBESGUE_TEST08:', 'Analyze Fejer1 points.', nfun)
    call run_test(9, 'LEBESGUE_TEST09:', 'Analyze Fejer2 points.', nfun)
    call timestamp()
  end do
end program main
