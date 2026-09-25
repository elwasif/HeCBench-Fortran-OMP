! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_signed_char
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
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

  character(len=256) :: arg0, arg1, arg2, arg3, arg4, arg5, arg6
  integer :: rows, cols, ncases, ncontrols, nthreads, repeat
  integer(int64) :: data_size
  integer(c_signed_char), allocatable :: data_t(:)
  real(real32), allocatable :: h_results(:), cpu_results(:)
  integer(int64) :: idx
  integer :: k, error_count
  real(real64) :: start_time, end_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 6) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <rows> <cols> <cases> <controls> <threads> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  call get_command_argument(6, arg6)
  read(arg1, *) rows
  read(arg2, *) cols
  read(arg3, *) ncases
  read(arg4, *) ncontrols
  read(arg5, *) nthreads
  read(arg6, *) repeat

  write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'Individuals=', rows, ' SNPs=', cols, &
    ' cases=', ncases, ' controls=', ncontrols, ' nthreads=', nthreads

  data_size = int(rows, int64) * int(cols, int64)
  write(*,'(A,I0)') 'Size of the data = ', data_size

  allocate(data_t(data_size), h_results(cols), cpu_results(cols))
  call c_srand(19937_c_int)
  do idx = 1, data_size
    data_t(idx) = int(ichar('0') + modulo(c_rand(), 3_c_int), c_signed_char)
  end do

  !$omp target data map(to: data_t(1:data_size)) map(from: h_results(1:cols))
  start_time = omp_get_wtime()
  do k = 1, repeat
    call chi2_kernel(rows, cols, ncases, ncontrols, nthreads, data_t, h_results)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time = ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'
  !$omp end target data

  start_time = omp_get_wtime()
  call cpu_kernel(rows, cols, ncases, ncontrols, data_t, cpu_results)
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Host execution time = ', end_time - start_time, ' (s)'

  error_count = 0
  do k = 1, cols
    if (abs(cpu_results(k) - h_results(k)) > 1.0e-4_real32) error_count = error_count + 1
  end do

  if (error_count == 0) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(data_t, h_results, cpu_results)

contains

  subroutine chi2_kernel(rows, cols, ncases, ncontrols, nthreads, data_t, results)
    integer, intent(in) :: rows, cols, ncases, ncontrols, nthreads
    integer(c_signed_char), intent(in) :: data_t(:)
    real(real32), intent(out) :: results(:)
    integer :: col, m, p, y
    integer :: case0, case1, case2, control0, control1, control2
    integer :: tot_cases, tot_controls, total
    real(real32) :: chisquare, expv, c_expected, con_expected, numerator1, numerator2

    !$omp target teams distribute parallel do simd thread_limit(nthreads) &
    !$omp& private(col, m, p, y, case0, case1, case2, control0, control1, control2, &
    !$omp& tot_cases, tot_controls, total, chisquare, expv, c_expected, con_expected, &
    !$omp& numerator1, numerator2)
    do col = 1, cols
      case0 = 1
      case1 = 1
      case2 = 1
      control0 = 1
      control1 = 1
      control2 = 1
      chisquare = 0.0_real32

      do m = 1, ncases
        y = int(data_t(int(m - 1, int64) * int(cols, int64) + int(col, int64)))
        if (y == ichar('0')) then
          case0 = case0 + 1
        else if (y == ichar('1')) then
          case1 = case1 + 1
        else if (y == ichar('2')) then
          case2 = case2 + 1
        end if
      end do

      do m = ncases + 1, ncases + ncontrols
        y = int(data_t(int(m - 1, int64) * int(cols, int64) + int(col, int64)))
        if (y == ichar('0')) then
          control0 = control0 + 1
        else if (y == ichar('1')) then
          control1 = control1 + 1
        else if (y == ichar('2')) then
          control2 = control2 + 1
        end if
      end do

      tot_cases = 1
      tot_controls = 1
      tot_cases = tot_cases + case0 + case1 + case2
      tot_controls = tot_controls + control0 + control1 + control2
      total = tot_cases + tot_controls

      do p = 1, 3
        select case (p)
        case (1)
          expv = real(case0 + control0, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case0, real32) - c_expected
          numerator2 = real(control0, real32) - con_expected
        case (2)
          expv = real(case1 + control1, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case1, real32) - c_expected
          numerator2 = real(control1, real32) - con_expected
        case default
          expv = real(case2 + control2, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case2, real32) - c_expected
          numerator2 = real(control2, real32) - con_expected
        end select
        chisquare = chisquare + numerator1 * numerator1 / c_expected + &
          numerator2 * numerator2 / con_expected
      end do
      results(col) = chisquare
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine chi2_kernel

  subroutine cpu_kernel(rows, cols, ncases, ncontrols, data_t, results)
    integer, intent(in) :: rows, cols, ncases, ncontrols
    integer(c_signed_char), intent(in) :: data_t(:)
    real(real32), intent(out) :: results(:)
    integer :: col, m, p, y
    integer :: case0, case1, case2, control0, control1, control2
    integer :: tot_cases, tot_controls, total
    real(real32) :: chisquare, expv, c_expected, con_expected, numerator1, numerator2

    tot_cases = 1
    tot_controls = 1
    do col = 1, cols
      case0 = 1
      case1 = 1
      case2 = 1
      control0 = 1
      control1 = 1
      control2 = 1
      chisquare = 0.0_real32

      do m = 1, ncases
        y = int(data_t(int(m - 1, int64) * int(cols, int64) + int(col, int64)))
        if (y == ichar('0')) then
          case0 = case0 + 1
        else if (y == ichar('1')) then
          case1 = case1 + 1
        else if (y == ichar('2')) then
          case2 = case2 + 1
        end if
      end do

      do m = ncases + 1, ncases + ncontrols
        y = int(data_t(int(m - 1, int64) * int(cols, int64) + int(col, int64)))
        if (y == ichar('0')) then
          control0 = control0 + 1
        else if (y == ichar('1')) then
          control1 = control1 + 1
        else if (y == ichar('2')) then
          control2 = control2 + 1
        end if
      end do

      tot_cases = tot_cases + case0 + case1 + case2
      tot_controls = tot_controls + control0 + control1 + control2
      total = tot_cases + tot_controls

      do p = 1, 3
        select case (p)
        case (1)
          expv = real(case0 + control0, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case0, real32) - c_expected
          numerator2 = real(control0, real32) - con_expected
        case (2)
          expv = real(case1 + control1, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case1, real32) - c_expected
          numerator2 = real(control1, real32) - con_expected
        case default
          expv = real(case2 + control2, real32)
          c_expected = real(tot_cases, real32) * expv / real(total, real32)
          con_expected = real(tot_controls, real32) * expv / real(total, real32)
          numerator1 = real(case2, real32) - c_expected
          numerator2 = real(control2, real32) - con_expected
        end select
        chisquare = chisquare + numerator1 * numerator1 / c_expected + &
          numerator2 * numerator2 / con_expected
      end do
      results(col) = chisquare
    end do
  end subroutine cpu_kernel

end program main
