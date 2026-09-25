! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer :: array_size
  integer :: num_times

  array_size = 33554432
  num_times = 100
  call parse_arguments(array_size, num_times)
  call run_float(array_size, num_times)
  call run_double(array_size, num_times)

contains

  subroutine parse_arguments(array_size, num_times)
    integer, intent(inout) :: array_size, num_times
    integer :: argc, i
    character(len=128) :: arg, value
    argc = command_argument_count()
    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('-s', '--arraysize')
        if (i == argc) stop 1
        call get_command_argument(i + 1, value)
        read(value, *) array_size
        if (array_size <= 0) stop 1
        i = i + 2
      case ('-n', '--numtimes')
        if (i == argc) stop 1
        call get_command_argument(i + 1, value)
        read(value, *) num_times
        if (num_times < 2) then
          print '(A)', 'Number of times must be 2 or more'
          stop 1
        end if
        i = i + 2
      case ('-h', '--help')
        print *
        print '(A)', 'Usage: ./main [OPTIONS]'
        print *
        print '(A)', 'Options:'
        print '(A)', '  -h  --help               Print the message'
        print '(A)', '  -s  --arraysize  SIZE    Use SIZE elements in the array'
        print '(A)', '  -n  --numtimes   NUM     Run the test NUM times (NUM >= 2)'
        print *
        stop
      case default
        write(*, '(A,A,A)') "Unrecognized argument '", trim(arg), "' (try '--help')"
        stop 1
      end select
    end do
    if (mod(array_size, 256) /= 0) then
      print '(A)', 'Array size must be a multiple of 256'
      stop 1
    end if
  end subroutine parse_arguments

  subroutine run_float(array_size, num_times)
    integer, intent(in) :: array_size, num_times
    real(real32), allocatable :: a(:), b(:), c(:)
    real(real64) :: timings(6, num_times)
    integer :: k

    allocate(a(array_size), b(array_size), c(array_size))
    print '(A,I0,A)', 'Running kernels ', num_times, ' times'
    print '(A)', 'Precision: float'
    write(*, '(5A)') 'Array size: ', trim(lead_zero(array_size * 4.0e-6_real64, 1)), ' MB (=', &
      trim(lead_zero(array_size * 4.0e-9_real64, 1)), ' GB)'
    write(*, '(5A)') 'Total size: ', trim(lead_zero(3.0_real64 * array_size * 4.0e-6_real64, 1)), ' MB (=', &
      trim(lead_zero(3.0_real64 * array_size * 4.0e-9_real64, 1)), ' GB)'
    !$omp target data map(alloc: a(1:array_size), b(1:array_size), c(1:array_size))
      call init_float(a, b, c, array_size)
      do k = 1, num_times
        timings(1, k) = time_copy_float(a, c, array_size)
        timings(2, k) = time_mul_float(b, c, array_size)
        timings(3, k) = time_add_float(a, b, c, array_size)
        timings(4, k) = time_triad_float(a, b, c, array_size)
        timings(5, k) = time_dot_float(a, b, array_size)
        timings(6, k) = time_nstream_float(a, b, c, array_size)
      end do
      call print_table(timings, num_times, array_size, 4)
    !$omp end target data
    deallocate(a, b, c)
  end subroutine run_float

  subroutine run_double(array_size, num_times)
    integer, intent(in) :: array_size, num_times
    real(real64), allocatable :: a(:), b(:), c(:)
    real(real64) :: timings(6, num_times)
    integer :: k

    allocate(a(array_size), b(array_size), c(array_size))
    print '(A,I0,A)', 'Running kernels ', num_times, ' times'
    print '(A)', 'Precision: double'
    write(*, '(5A)') 'Array size: ', trim(lead_zero(array_size * 8.0e-6_real64, 1)), ' MB (=', &
      trim(lead_zero(array_size * 8.0e-9_real64, 1)), ' GB)'
    write(*, '(5A)') 'Total size: ', trim(lead_zero(3.0_real64 * array_size * 8.0e-6_real64, 1)), ' MB (=', &
      trim(lead_zero(3.0_real64 * array_size * 8.0e-9_real64, 1)), ' GB)'
    !$omp target data map(alloc: a(1:array_size), b(1:array_size), c(1:array_size))
      call init_double(a, b, c, array_size)
      do k = 1, num_times
        timings(1, k) = time_copy_double(a, c, array_size)
        timings(2, k) = time_mul_double(b, c, array_size)
        timings(3, k) = time_add_double(a, b, c, array_size)
        timings(4, k) = time_triad_double(a, b, c, array_size)
        timings(5, k) = time_dot_double(a, b, array_size)
        timings(6, k) = time_nstream_double(a, b, c, array_size)
      end do
      call print_table(timings, num_times, array_size, 8)
    !$omp end target data
    deallocate(a, b, c)
  end subroutine run_double

  subroutine init_float(a, b, c, n)
    real(real32), intent(inout) :: a(:), b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = 0.1_real32
      b(i) = 0.2_real32
      c(i) = 0.0_real32
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine init_float

  subroutine init_double(a, b, c, n)
    real(real64), intent(inout) :: a(:), b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = 0.1_real64
      b(i) = 0.2_real64
      c(i) = 0.0_real64
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine init_double

  real(real64) function time_copy_float(a, c, n)
    real(real32), intent(in) :: a(:)
    real(real32), intent(inout) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      c(i) = a(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_copy_float = omp_get_wtime() - t
  end function time_copy_float

  real(real64) function time_mul_float(b, c, n)
    real(real32), intent(inout) :: b(:)
    real(real32), intent(in) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      b(i) = 0.4_real32 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_mul_float = omp_get_wtime() - t
  end function time_mul_float

  real(real64) function time_add_float(a, b, c, n)
    real(real32), intent(in) :: a(:), b(:)
    real(real32), intent(inout) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      c(i) = a(i) + b(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_add_float = omp_get_wtime() - t
  end function time_add_float

  real(real64) function time_triad_float(a, b, c, n)
    real(real32), intent(inout) :: a(:)
    real(real32), intent(in) :: b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = b(i) + 0.4_real32 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_triad_float = omp_get_wtime() - t
  end function time_triad_float

  real(real64) function time_dot_float(a, b, n)
    real(real32), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    integer :: i
    real(real32) :: sum
    real(real64) :: t
    sum = 0.0_real32
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd map(tofrom: sum) thread_limit(256) reduction(+:sum)
    do i = 1, n
      sum = sum + a(i) * b(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_dot_float = omp_get_wtime() - t
  end function time_dot_float

  real(real64) function time_nstream_float(a, b, c, n)
    real(real32), intent(inout) :: a(:)
    real(real32), intent(in) :: b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = a(i) + b(i) + 0.4_real32 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_nstream_float = omp_get_wtime() - t
  end function time_nstream_float

  real(real64) function time_copy_double(a, c, n)
    real(real64), intent(in) :: a(:)
    real(real64), intent(inout) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      c(i) = a(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_copy_double = omp_get_wtime() - t
  end function time_copy_double

  real(real64) function time_mul_double(b, c, n)
    real(real64), intent(inout) :: b(:)
    real(real64), intent(in) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      b(i) = 0.4_real64 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_mul_double = omp_get_wtime() - t
  end function time_mul_double

  real(real64) function time_add_double(a, b, c, n)
    real(real64), intent(in) :: a(:), b(:)
    real(real64), intent(inout) :: c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      c(i) = a(i) + b(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_add_double = omp_get_wtime() - t
  end function time_add_double

  real(real64) function time_triad_double(a, b, c, n)
    real(real64), intent(inout) :: a(:)
    real(real64), intent(in) :: b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = b(i) + 0.4_real64 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_triad_double = omp_get_wtime() - t
  end function time_triad_double

  real(real64) function time_dot_double(a, b, n)
    real(real64), intent(in) :: a(:), b(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: sum, t
    sum = 0.0_real64
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd map(tofrom: sum) thread_limit(256) reduction(+:sum)
    do i = 1, n
      sum = sum + a(i) * b(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_dot_double = omp_get_wtime() - t
  end function time_dot_double

  real(real64) function time_nstream_double(a, b, c, n)
    real(real64), intent(inout) :: a(:)
    real(real64), intent(in) :: b(:), c(:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: t
    t = omp_get_wtime()
    !$omp target teams distribute parallel do simd thread_limit(256)
    do i = 1, n
      a(i) = a(i) + b(i) + 0.4_real64 * c(i)
    end do
    !$omp end target teams distribute parallel do simd
    time_nstream_double = omp_get_wtime() - t
  end function time_nstream_double

  subroutine print_table(timings, num_times, array_size, bytes_per_value)
    real(real64), intent(in) :: timings(6, num_times)
    integer, intent(in) :: num_times, array_size, bytes_per_value
    character(len=8), parameter :: labels(6) = ['Copy    ', 'Mul     ', 'Add     ', 'Triad   ', 'Dot     ', 'Nstream ']
    integer, parameter :: factors(6) = [2, 2, 3, 3, 2, 4]
    integer :: i
    real(real64) :: min_t, max_t, avg_t, bandwidth
    character(len=12) :: label_field, bandwidth_field, min_field, max_field, avg_field
    print '(A)', 'Function    MBytes/sec  Min (sec)   Max         Average     '
    do i = 1, 6
      min_t = minval(timings(i, 2:num_times))
      max_t = maxval(timings(i, 2:num_times))
      avg_t = sum(timings(i, 2:num_times)) / real(num_times - 1, real64)
      bandwidth = 1.0e-6_real64 * real(factors(i), real64) * &
        real(bytes_per_value, real64) * real(array_size, real64) / min_t
      call fixed_field(bandwidth, 3, bandwidth_field)
      call fixed_field(min_t, 5, min_field)
      call fixed_field(max_t, 5, max_field)
      call fixed_field(avg_t, 5, avg_field)
      label_field = adjustl(labels(i))
      write(*, '(5A)') label_field, bandwidth_field, min_field, max_field, avg_field
    end do
    print *
  end subroutine print_table

  subroutine fixed_field(value, digits, field)
    real(real64), intent(in) :: value
    integer, intent(in) :: digits
    character(len=12), intent(out) :: field

    field = adjustl(lead_zero(value, digits))
  end subroutine fixed_field

  character(len=32) function lead_zero(value, digits)
    real(real64), intent(in) :: value
    integer, intent(in) :: digits
    character(len=32) :: buffer

    select case (digits)
    case (1)
      write(buffer, '(F0.1)') value
    case (3)
      write(buffer, '(F0.3)') value
    case default
      write(buffer, '(F0.5)') value
    end select
    if (buffer(1:1) == '.') then
      lead_zero = '0' // trim(buffer)
    else if (buffer(1:2) == '-.') then
      lead_zero = '-0' // trim(buffer(2:))
    else
      lead_zero = trim(buffer)
    end if
  end function lead_zero

end program main
