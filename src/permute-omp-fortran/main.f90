! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: block_sizes(6) = [32, 64, 128, 256, 512, 1024]
  integer, parameter :: t_size = 1024
  integer, parameter :: c_size = 768
  integer, parameter :: num_heads = 12
  real(real32), parameter :: rand_max_real = 2147483647.0_real32

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: batch_size, repeat_times, block_size, i
  integer(int64) :: s, total_in
  real(real32), allocatable :: inp(:), out(:), q(:), k(:), v(:)
  real(real64) :: start_time, elapsed_ms

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <batch size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) batch_size
  call get_command_argument(2, arg); read(arg, *) repeat_times
  if (batch_size <= 0 .or. repeat_times <= 0) stop 1

  s = int(batch_size, int64) * int(t_size, int64) * int(c_size, int64)
  total_in = 3_int64 * s
  allocate(inp(total_in), out(total_in), q(s), k(s), v(s))

  call c_srand(1_c_int)
  call fill_random(inp)
  call fill_random(out)
  call fill_random(q)
  call fill_random(k)
  call fill_random(v)

  call permute_cpu(inp, q, k, v, batch_size, t_size, c_size, num_heads)

  !$omp target data map(to: inp(1:total_in)) map(tofrom: out(1:total_in))
  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    write(*,'(A,I0,A)') 'Checking block size ', block_size, '.'
    call permute(out, inp, batch_size, t_size, c_size, num_heads, block_size)
    !$omp target update from(out(1:total_in))
    call validate_result(out(1:s), q, 'q', s, 1.0e-6_real32)
    call validate_result(out(s + 1:2 * s), k, 'k', s, 1.0e-6_real32)
    call validate_result(out(2 * s + 1:3 * s), v, 'v', s, 1.0e-6_real32)
  end do

  write(*,'(A)') 'All results match. Starting benchmarks.'
  write(*,*)

  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    start_time = omp_get_wtime()
    call run_repeated(out, inp, batch_size, t_size, c_size, num_heads, block_size, repeat_times)
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat_times, real64)
    write(*,'(A,I4,A,F0.6,A)') 'block_size ', block_size, ' | time ', elapsed_ms, ' ms'
  end do
  !$omp end target data

  deallocate(inp, out, q, k, v)

contains

  subroutine fill_random(values)
    real(real32), intent(out) :: values(:)
    integer :: i

    do i = 1, size(values)
      ! Preserve the C++ original's implicit srand(1)/rand() input stream.
      values(i) = real(c_rand(), real32) / rand_max_real * 2.0_real32 - 1.0_real32
    end do
  end subroutine fill_random

  subroutine permute_cpu(inp, q, k, v, bsz, tsz, csz, nh)
    real(real32), intent(in) :: inp(:)
    real(real32), intent(out) :: q(:), k(:), v(:)
    integer, intent(in) :: bsz, tsz, csz, nh
    integer :: b, n, t, c, idx, head_size

    head_size = csz / nh
    idx = 1
    do b = 0, bsz - 1
      do n = 0, nh - 1
        do t = 0, tsz - 1
          do c = n * head_size, (n + 1) * head_size - 1
            q(idx) = inp(b * tsz * 3 * csz + t * 3 * csz + c + 1)
            k(idx) = inp(b * tsz * 3 * csz + t * 3 * csz + csz + c + 1)
            v(idx) = inp(b * tsz * 3 * csz + t * 3 * csz + 2 * csz + c + 1)
            idx = idx + 1
          end do
        end do
      end do
    end do
  end subroutine permute_cpu

  subroutine permute(out, inp, bsz, tsz, csz, nh, block_size)
    real(real32), intent(out) :: out(:)
    real(real32), intent(in) :: inp(:)
    integer, intent(in) :: bsz, tsz, csz, nh, block_size
    integer :: idx, b, rest, nh_idx, n, d_idx, input_idx, head_size
    integer(int64) :: total_threads, num_blocks, plane

    head_size = csz / nh
    total_threads = int(bsz, int64) * int(tsz, int64) * int(csz, int64)
    num_blocks = (total_threads + int(block_size, int64) - 1_int64) / int(block_size, int64)
    plane = total_threads

    !$omp target teams distribute parallel do private(b, rest, nh_idx, n, d_idx, input_idx) &
    !$omp& num_teams(num_blocks) thread_limit(block_size) num_threads(block_size)
    do idx = 0, int(total_threads) - 1
      b = idx / (csz * tsz)
      rest = modulo(idx, csz * tsz)
      nh_idx = rest / (tsz * head_size)
      rest = modulo(rest, tsz * head_size)
      n = rest / head_size
      d_idx = modulo(rest, head_size)
      input_idx = b * tsz * 3 * csz + n * 3 * csz + nh_idx * head_size + d_idx

      out(idx + 1) = inp(input_idx + 1)
      out(plane + idx + 1) = inp(input_idx + csz + 1)
      out(2 * plane + idx + 1) = inp(input_idx + 2 * csz + 1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine permute

  subroutine run_repeated(out, inp, bsz, tsz, csz, nh, block_size, repeat_times)
    real(real32), intent(out) :: out(:)
    real(real32), intent(in) :: inp(:)
    integer, intent(in) :: bsz, tsz, csz, nh, block_size, repeat_times
    integer :: iter

    do iter = 1, repeat_times
      call permute(out, inp, bsz, tsz, csz, nh, block_size)
    end do
  end subroutine run_repeated

  subroutine validate_result(device_result, cpu_reference, name, num_elements, tolerance)
    real(real32), intent(in) :: device_result(:), cpu_reference(:), tolerance
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: num_elements
    integer(int64) :: i
    integer :: nfaults
    real(real32) :: fp_epsilon, t_eff

    nfaults = 0
    fp_epsilon = epsilon(1.0_real32)
    do i = 1, num_elements
      if (cpu_reference(i) /= cpu_reference(i) .or. abs(cpu_reference(i)) > huge(cpu_reference(i))) cycle
      t_eff = tolerance + abs(cpu_reference(i)) * fp_epsilon
      if (abs(cpu_reference(i) - device_result(i)) > t_eff) then
        write(*,'(A,A,A,I0,A,F0.6,A,F0.6)') 'Mismatch of ', trim(name), ' at ', i - 1, &
          ': CPU_ref: ', cpu_reference(i), ' vs GPU: ', device_result(i)
        nfaults = nfaults + 1
        if (nfaults >= 10) stop 1
      end if
    end do
    if (nfaults > 0) stop 1
  end subroutine validate_result

end program main
