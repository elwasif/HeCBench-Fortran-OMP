! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int, c_float, c_int64_t
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    subroutine silu_make_random_float(array, n) bind(C, name='silu_make_random_float')
      import :: c_float, c_int64_t
      real(c_float), intent(out) :: array(*)
      integer(c_int64_t), value :: n
    end subroutine silu_make_random_float
  end interface

  integer, parameter :: block_sizes(5) = [64, 128, 256, 512, 1024]
  integer :: b, c, h, w, repeat, bs, i
  integer(int64) :: n
  real(real32), allocatable :: x(:), dout(:), out_ref(:), dx_ref(:), d_out(:), d_dx(:)
  real(real64) :: elapsed_ms

  if (command_argument_count() /= 5) then
    print '(A)', 'Usage: ./main <B> <C> <H> <W> <repeat>'
    stop 1
  end if

  b = read_arg(1)
  c = read_arg(2)
  h = read_arg(3)
  w = read_arg(4)
  repeat = read_arg(5)
  if (b <= 0 .or. c <= 0 .or. h <= 0 .or. w <= 0 .or. repeat <= 0) stop 1
  n = int(b, int64) * int(c, int64) * int(h, int64) * int(w, int64)

  allocate(x(n), dout(n), out_ref(n), dx_ref(n), d_out(n), d_dx(n))
  call c_srand(0_c_int)
  call silu_make_random_float(x, n)
  call silu_make_random_float(dout, n)

  call silu_forward_reference(x, out_ref, n)
  call silu_backward_reference(dout, x, dx_ref, n)

  !$omp target data map(to: x(1:n), dout(1:n)) map(alloc: d_out(1:n), d_dx(1:n))
  print '(A)', 'Checking forward pass'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    print '(A,I0)', 'Checking block size ', bs
    call silu_forward(x, d_out, n, bs)
    call validate_result(d_out, out_ref, n)
  end do

  print '(A)', 'Checking forward2 pass'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    print '(A,I0)', 'Checking block size ', bs
    call silu_forward2(x, d_out, n, bs)
    call validate_result(d_out, out_ref, n)
  end do

  print '(A)', 'Checking backward pass'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    print '(A,I0)', 'Checking block size ', bs
    call silu_backward(dout, x, d_dx, n, bs)
    call validate_result(d_dx, dx_ref, n)
  end do

  print '(A)', 'Checking backward2 pass'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    print '(A,I0)', 'Checking block size ', bs
    call silu_backward2(dout, x, d_dx, n, bs)
    call validate_result(d_dx, dx_ref, n)
  end do

  print '(A)'
  print '(A)', 'Forward pass benchmarks:'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    elapsed_ms = benchmark_forward(x, d_out, n, bs, repeat)
    print '(A,I4,A,F6.4,A)', 'block_size ', bs, ' | time ', elapsed_ms, ' ms'
  end do

  print '(A)'
  print '(A)', 'Forward2 pass benchmarks:'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    elapsed_ms = benchmark_forward2(x, d_out, n, bs, repeat)
    print '(A,I4,A,F6.4,A)', 'block_size ', bs, ' | time ', elapsed_ms, ' ms'
  end do

  print '(A)'
  print '(A)', 'Backward pass benchmarks:'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    elapsed_ms = benchmark_backward(dout, x, d_dx, n, bs, repeat)
    print '(A,I4,A,F6.4,A)', 'block_size ', bs, ' | time ', elapsed_ms, ' ms'
  end do

  print '(A)'
  print '(A)', 'Backward2 pass benchmarks:'
  do i = 1, size(block_sizes)
    bs = block_sizes(i)
    elapsed_ms = benchmark_backward2(dout, x, d_dx, n, bs, repeat)
    print '(A,I4,A,F6.4,A)', 'block_size ', bs, ' | time ', elapsed_ms, ' ms'
  end do
  !$omp end target data

  deallocate(x, dout, out_ref, dx_ref, d_out, d_dx)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  real(real32) function silu_value(xv)
    real(real32), intent(in) :: xv
    silu_value = xv / (1.0_real32 + exp(-xv))
  end function silu_value

  real(real32) function silu_grad(xv)
    real(real32), intent(in) :: xv
    real(real32) :: sig
    sig = 1.0_real32 / (1.0_real32 + exp(-xv))
    silu_grad = sig * (1.0_real32 + xv * (1.0_real32 - sig))
  end function silu_grad

  subroutine silu_forward_reference(x, out, n)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer(int64), intent(in) :: n
    integer :: i
    do i = 1, int(n)
      out(i) = silu_value(x(i))
    end do
  end subroutine silu_forward_reference

  subroutine silu_backward_reference(dout, x, dx, n)
    real(real32), intent(in) :: dout(:), x(:)
    real(real32), intent(out) :: dx(:)
    integer(int64), intent(in) :: n
    integer :: i
    do i = 1, int(n)
      dx(i) = dout(i) * silu_grad(x(i))
    end do
  end subroutine silu_backward_reference

  subroutine silu_forward(x, out, n, block_size)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size
    integer :: i, num_teams
    num_teams = get_num_teams(n, block_size)
    !$omp target teams distribute parallel do num_teams(num_teams) thread_limit(block_size)
    do i = 1, int(n)
      out(i) = silu_value(x(i))
    end do
    !$omp end target teams distribute parallel do
  end subroutine silu_forward

  subroutine silu_forward2(x, out, n, block_size)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size
    integer :: i, k, vec_count, tail_start, num_teams, tail_threads
    vec_count = int(n) / 4
    tail_start = 4 * vec_count + 1
    num_teams = get_num_teams(int(vec_count, int64), block_size)
    !$omp target teams distribute parallel do num_teams(num_teams) thread_limit(block_size) private(k)
    do i = 0, vec_count - 1
      do k = 1, 4
        out(4 * i + k) = silu_value(x(4 * i + k))
      end do
    end do
    !$omp end target teams distribute parallel do
    tail_threads = int(n) - 4 * vec_count
    if (tail_threads > 0) then
      !$omp target teams distribute parallel do num_teams(1) thread_limit(tail_threads)
      do i = tail_start, int(n)
        out(i) = silu_value(x(i))
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine silu_forward2

  subroutine silu_backward(dout, x, dx, n, block_size)
    real(real32), intent(in) :: dout(:), x(:)
    real(real32), intent(out) :: dx(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size
    integer :: i, num_teams
    num_teams = get_num_teams(n, block_size)
    !$omp target teams distribute parallel do num_teams(num_teams) thread_limit(block_size)
    do i = 1, int(n)
      dx(i) = dout(i) * silu_grad(x(i))
    end do
    !$omp end target teams distribute parallel do
  end subroutine silu_backward

  subroutine silu_backward2(dout, x, dx, n, block_size)
    real(real32), intent(in) :: dout(:), x(:)
    real(real32), intent(out) :: dx(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size
    integer :: i, k, vec_count, tail_start, num_teams, tail_threads
    vec_count = int(n) / 4
    tail_start = 4 * vec_count + 1
    num_teams = get_num_teams(int(vec_count, int64), block_size)
    !$omp target teams distribute parallel do num_teams(num_teams) thread_limit(block_size) private(k)
    do i = 0, vec_count - 1
      do k = 1, 4
        dx(4 * i + k) = dout(4 * i + k) * silu_grad(x(4 * i + k))
      end do
    end do
    !$omp end target teams distribute parallel do
    tail_threads = int(n) - 4 * vec_count
    if (tail_threads > 0) then
      !$omp target teams distribute parallel do num_teams(1) thread_limit(tail_threads)
      do i = tail_start, int(n)
        dx(i) = dout(i) * silu_grad(x(i))
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine silu_backward2

  integer function get_num_teams(num_elements, num_threads)
    integer(int64), intent(in) :: num_elements
    integer, intent(in) :: num_threads
    get_num_teams = int((num_elements + int(num_threads, int64) - 1_int64) / int(num_threads, int64))
  end function get_num_teams

  subroutine validate_result(device_result, cpu_reference, n)
    real(real32), intent(inout) :: device_result(:)
    real(real32), intent(in) :: cpu_reference(:)
    integer(int64), intent(in) :: n
    integer :: i, faults
    !$omp target update from(device_result(1:n))
    faults = 0
    do i = 1, int(n)
      if (abs(cpu_reference(i) - device_result(i)) > 1.0e-4_real32) then
        faults = faults + 1
        if (faults >= 10) exit
      end if
    end do
    if (faults == 0) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end subroutine validate_result

  real(real64) function benchmark_forward(x, out, n, block_size, repeat)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size, repeat
    integer :: i
    real(real64) :: start_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call silu_forward(x, out, n, block_size)
    end do
    benchmark_forward = ((omp_get_wtime() - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_forward

  real(real64) function benchmark_forward2(x, out, n, block_size, repeat)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: out(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size, repeat
    integer :: i
    real(real64) :: start_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call silu_forward2(x, out, n, block_size)
    end do
    benchmark_forward2 = ((omp_get_wtime() - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_forward2

  real(real64) function benchmark_backward(dout, x, dx, n, block_size, repeat)
    real(real32), intent(in) :: dout(:), x(:)
    real(real32), intent(out) :: dx(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size, repeat
    integer :: i
    real(real64) :: start_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call silu_backward(dout, x, dx, n, block_size)
    end do
    benchmark_backward = ((omp_get_wtime() - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_backward

  real(real64) function benchmark_backward2(dout, x, dx, n, block_size, repeat)
    real(real32), intent(in) :: dout(:), x(:)
    real(real32), intent(out) :: dx(:)
    integer(int64), intent(in) :: n
    integer, intent(in) :: block_size, repeat
    integer :: i
    real(real64) :: start_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call silu_backward2(dout, x, dx, n, block_size)
    end do
    benchmark_backward2 = ((omp_get_wtime() - start_time) * 1.0e3_real64) / real(repeat, real64)
  end function benchmark_backward2

end program main
