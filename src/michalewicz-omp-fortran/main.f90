! SPDX-License-Identifier: CC0-1.0
module michalewicz_mod
  use iso_c_binding, only: c_int, c_ptr, c_size_t, c_loc
  use iso_fortran_env, only: int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  real(real32), parameter :: pi32 = acos(-1.0_real32)

  interface
    subroutine mt19937_reset(seed) bind(C, name='michalewicz_mt19937_reset')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine mt19937_reset

    subroutine mt19937_fill(values, count) bind(C, name='michalewicz_mt19937_fill')
      import :: c_ptr, c_size_t
      type(c_ptr), value :: values
      integer(c_size_t), value :: count
    end subroutine mt19937_fill
  end interface

contains

  subroutine fill_values(values)
    real(real32), target, intent(out) :: values(0:)

    call mt19937_fill(c_loc(values(0)), int(size(values), c_size_t))
  end subroutine fill_values

  pure real(real32) function michalewicz_cpu(values, offset, dim) result(value)
    real(real32), intent(in) :: values(0:)
    integer(int64), intent(in) :: offset
    integer, intent(in) :: dim
    integer :: k
    real(real32) :: x, a, b, c

    value = 0.0_real32
    do k = 0, dim - 1
      x = values(offset + k)
      a = sin(x)
      b = sin((real(k + 1, real32) * x * x) / pi32)
      c = b ** 20
      value = value + a * c
    end do
    value = -value
  end function michalewicz_cpu

  real(real32) function cpu_minimum(values, n, dim) result(min_value)
    real(real32), intent(in) :: values(0:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: dim
    integer(int64) :: j

    min_value = 0.0_real32
    do j = 0, n - 1
      min_value = min(min_value, michalewicz_cpu(values, j * dim, dim))
    end do
  end function cpu_minimum

  subroutine print_error(value, dim)
    real(real32), intent(in) :: value
    integer, intent(in) :: dim
    real(real32) :: true_min
    character(len=32) :: value_text, error_text

    value_text = fixed6(value)
    write(*,'(A,A)') 'Global minima = ', trim(value_text)
    true_min = 0.0_real32
    if (dim == 2) then
      true_min = -1.8013_real32
    else if (dim == 5) then
      true_min = -4.687658_real32
    else if (dim == 10) then
      true_min = -9.66015_real32
    end if
    error_text = fixed6(abs(true_min - value))
    write(*,'(A,A)') 'Error = ', trim(error_text)
  end subroutine print_error

  pure character(len=32) function fixed6(value) result(text)
    real(real32), intent(in) :: value
    character(len=32) :: tmp

    write(tmp, '(F0.6)') value
    tmp = adjustl(tmp)
    if (tmp(1:1) == '.') then
      text = '0' // trim(tmp)
    else if (tmp(1:2) == '-.') then
      text = '-0' // trim(tmp(2:))
    else
      text = trim(tmp)
    end if
  end function fixed6

  subroutine run_dimension(n, repeat, dim)
    integer(int64), intent(in) :: n
    integer, intent(in) :: repeat, dim
    integer(int64) :: total_size, j
    integer :: iter, k
    real(real32), allocatable :: values(:)
    real(real32) :: min_value, reference_value, x, a, b, c, candidate
    real(real64) :: start_time, elapsed_ns

    total_size = n * int(dim, int64)
    allocate(values(0:total_size - 1))
    call fill_values(values)

    min_value = 0.0_real32

    !$omp target data map(to: values(0:total_size - 1)) map(tofrom: min_value)
    start_time = omp_get_wtime()
    do iter = 1, repeat
      !$omp target teams distribute parallel do thread_limit(block_size) reduction(min: min_value) &
      !$omp& private(k, x, a, b, c, candidate)
      do j = 0, n - 1
        candidate = 0.0_real32
        do k = 0, dim - 1
          x = values(j * dim + k)
          a = sin(x)
          b = sin((real(k + 1, real32) * x * x) / pi32)
          c = b ** 20
          candidate = candidate + a * c
        end do
        candidate = -candidate
        min_value = min(min_value, candidate)
      end do
      !$omp end target teams distribute parallel do
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,I0,A,F0.6,A)') 'Average execution time of kernel (dim = ', dim, '): ', &
      elapsed_ns * 1.0e-3_real64 / real(repeat, real64), ' (us)'
    !$omp end target data

    reference_value = cpu_minimum(values, n, dim)
    if (abs(reference_value - min_value) > 1.0e-4_real32) then
      write(*,'(A,I0,A,F0.6,A,F0.6)') 'FAIL dim ', dim, ': CPU reference ', reference_value, &
        ' GPU result ', min_value
      stop 2
    end if

    call print_error(min_value, dim)
    deallocate(values)
  end subroutine run_dimension

end module michalewicz_mod

program main
  use iso_fortran_env, only: int64
  use michalewicz_mod
  implicit none

  integer :: argc, repeat, d
  integer(int64) :: n
  integer, parameter :: dims(3) = [2, 5, 10]
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage: ', trim(prog), ' <number of vectors> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  n = parse_atol(arg)
  call get_command_argument(2, arg)
  repeat = parse_atoi(arg)

  call mt19937_reset(19937_c_int)
  do d = 1, size(dims)
    call run_dimension(n, repeat, dims(d))
  end do
contains

  integer(int64) function parse_atol(text) result(value)
    character(len=*), intent(in) :: text
    integer :: i, sign, digit

    value = 0_int64
    i = 1
    do while (i <= len_trim(text) .and. is_space(text(i:i)))
      i = i + 1
    end do

    sign = 1
    if (i <= len_trim(text)) then
      if (text(i:i) == '-') then
        sign = -1
        i = i + 1
      else if (text(i:i) == '+') then
        i = i + 1
      end if
    end if

    do while (i <= len_trim(text))
      digit = iachar(text(i:i)) - iachar('0')
      if (digit < 0 .or. digit > 9) exit
      value = value * 10_int64 + int(digit, int64)
      i = i + 1
    end do
    value = value * int(sign, int64)
  end function parse_atol

  integer function parse_atoi(text) result(value)
    character(len=*), intent(in) :: text

    value = int(parse_atol(text))
  end function parse_atoi

  logical function is_space(ch) result(space)
    character(len=1), intent(in) :: ch
    integer :: code

    code = iachar(ch)
    space = ch == ' ' .or. (code >= 9 .and. code <= 13)
  end function is_space
end program main
