! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_int64_t
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine moe_sum_initialize_input(input, input_size, topk) bind(C, name="moe_sum_initialize_input")
      import :: c_float, c_int, c_int64_t
      real(c_float), intent(out) :: input(*)
      integer(c_int64_t), value :: input_size
      integer(c_int), value :: topk
    end subroutine moe_sum_initialize_input
  end interface

  character(len=256) :: arg0, arg
  integer(int32) :: num_tokens, hidden_size, repeat
  integer(int64) :: output_size, input_size
  real(real32), allocatable :: input(:), output(:), output_vec(:), r_output(:)
  real(real32) :: scalar_bandwidth
  integer :: topk

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of tokens> <hidden size> <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) num_tokens
  call get_command_argument(2, arg)
  read(arg, *) hidden_size
  call get_command_argument(3, arg)
  read(arg, *) repeat

  if (modulo(hidden_size, 4) /= 0) then
    write(*,'(A)') 'Hidden size is a multiple of four'
    stop 1
  end if

  output_size = int(num_tokens, int64) * int(hidden_size, int64)
  allocate(output(output_size), output_vec(output_size), r_output(output_size))

  !$omp target enter data map(alloc: output(1:output_size), output_vec(1:output_size))
  do topk = 2, 4
    input_size = output_size * int(topk, int64)
    allocate(input(input_size))
    !$omp target enter data map(to: input(1:input_size))
    call initialize_input(input, topk)
    call moe_sum_ref(topk, r_output, input, num_tokens, hidden_size)
    !$omp target update to(input(1:input_size))

    call run_moe_sum(topk, input, output, r_output, num_tokens, hidden_size, repeat, scalar_bandwidth)
    call run_moe_sum_vec(topk, input, output, output_vec, num_tokens, hidden_size, repeat, scalar_bandwidth)

    !$omp target exit data map(delete: input(1:input_size))
    deallocate(input)
  end do
  !$omp target exit data map(delete: output(1:output_size), output_vec(1:output_size))

  deallocate(output, output_vec, r_output)

contains

  subroutine initialize_input(input, topk)
    real(real32), intent(out) :: input(:)
    integer, intent(in) :: topk

    call moe_sum_initialize_input(input, int(size(input, kind=int64), c_int64_t), int(topk, c_int))
  end subroutine initialize_input

  subroutine moe_sum_ref(topk, out, input, num_tokens, hidden_size)
    integer, intent(in) :: topk, num_tokens, hidden_size
    real(real32), intent(out) :: out(:)
    real(real32), intent(in) :: input(:)
    integer :: block, idx, k
    integer(int64) :: output_base, input_base
    real(real32) :: x

    do block = 0, num_tokens - 1
      output_base = int(block, int64) * int(hidden_size, int64)
      input_base = output_base * int(topk, int64)
      do idx = 0, hidden_size - 1
        x = 0.0_real32
        do k = 0, topk - 1
          x = x + input(input_base + int(k * hidden_size + idx, int64) + 1)
        end do
        out(output_base + int(idx, int64) + 1) = x
      end do
    end do
  end subroutine moe_sum_ref

  subroutine run_moe_sum(topk, input, output, r_output, num_tokens, hidden_size, repeat, bandwidth)
    integer, intent(in) :: topk, num_tokens, hidden_size, repeat
    real(real32), intent(in) :: input(:), r_output(:)
    real(real32), intent(inout) :: output(:)
    real(real32), intent(out) :: bandwidth
    integer :: iter, block_size
    real(real64) :: start_time, elapsed_ns
    real(real32) :: io_bytes
    logical :: ok

    block_size = min(hidden_size, 1024)
    do iter = 1, 100
      call moe_sum_kernel(topk, output, input, hidden_size, num_tokens, block_size)
    end do

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call moe_sum_kernel(topk, output, input, hidden_size, num_tokens, block_size)
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,I0,A,F0.6,A)') 'Average execution time of kernel (TopK = ', topk, '): ', &
         elapsed_ns * 1.0e-3_real64 / real(repeat, real64), ' (us)'

    !$omp target update from(output(1:int(num_tokens, int64) * int(hidden_size, int64)))
    ok = all(abs(output - r_output) <= 1.0e-4_real32)
    write(*,'(A)') merge('PASS', 'FAIL', ok)
    io_bytes = real(repeat, real32) * &
         real(size(input, kind=int64) + size(output, kind=int64), real32) * 4.0_real32
    bandwidth = io_bytes / real(elapsed_ns, real32)
    write(*,'(A,F0.6,A)') 'Kernel bandwidth: ', bandwidth, ' GB/s '
  end subroutine run_moe_sum

  subroutine run_moe_sum_vec(topk, input, output, output_vec, num_tokens, hidden_size, repeat, scalar_bandwidth)
    integer, intent(in) :: topk, num_tokens, hidden_size, repeat
    real(real32), intent(in) :: scalar_bandwidth
    real(real32), intent(in) :: input(:), output(:)
    real(real32), intent(inout) :: output_vec(:)
    integer :: iter, block_size
    real(real64) :: start_time, elapsed_ns
    real(real32) :: io_bytes, bandwidth, bandwidth_vec, pct
    logical :: ok

    block_size = min(hidden_size / 4, 1024)
    do iter = 1, 100
      call moe_sum_kernel_vec(topk, output_vec, input, hidden_size, num_tokens, block_size)
    end do

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call moe_sum_kernel_vec(topk, output_vec, input, hidden_size, num_tokens, block_size)
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,I0,A,F0.6,A)') 'Average execution time of vec4 kernel (TopK = ', topk, '): ', &
         elapsed_ns * 1.0e-3_real64 / real(repeat, real64), ' (us)'

    !$omp target update from(output_vec(1:int(num_tokens, int64) * int(hidden_size, int64)))
    ok = bitwise_equal(output, output_vec)
    write(*,'(A)') merge('PASS', 'FAIL', ok)
    io_bytes = real(repeat, real32) * &
         real(size(input, kind=int64) + size(output_vec, kind=int64), real32) * 4.0_real32
    bandwidth_vec = io_bytes / real(elapsed_ns, real32)
    bandwidth = scalar_bandwidth
    pct = 100.0_real32 * (bandwidth_vec - bandwidth) / bandwidth
    write(*,'(A,F0.6,A,F0.6,A)') 'Kernel(vec4) bandwidth: ', bandwidth_vec, ' GB/s (', pct, '%)'
  end subroutine run_moe_sum_vec

  logical function bitwise_equal(lhs, rhs) result(equal)
    real(real32), intent(in) :: lhs(:), rhs(:)
    integer(int32), allocatable :: lhs_bits(:), rhs_bits(:)

    allocate(lhs_bits(size(lhs)), rhs_bits(size(rhs)))
    lhs_bits = transfer(lhs, lhs_bits)
    rhs_bits = transfer(rhs, rhs_bits)
    equal = all(lhs_bits == rhs_bits)
    deallocate(lhs_bits, rhs_bits)
  end function bitwise_equal

  subroutine moe_sum_kernel(topk, out, input, d, num_blocks, block_size)
    integer, intent(in) :: topk, d, num_blocks, block_size
    real(real32), intent(inout) :: out(:)
    real(real32), intent(in) :: input(:)
    integer :: block, idx, k
    integer(int64) :: output_base, input_base
    real(real32) :: x

    !$omp target teams distribute parallel do collapse(2) num_threads(block_size) private(output_base, input_base, k, x)
    do block = 0, num_blocks - 1
      do idx = 0, d - 1
        output_base = int(block, int64) * int(d, int64)
        input_base = output_base * int(topk, int64)
        x = 0.0_real32
        do k = 0, topk - 1
          x = x + input(input_base + int(k * d + idx, int64) + 1)
        end do
        out(output_base + int(idx, int64) + 1) = x
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine moe_sum_kernel

  subroutine moe_sum_kernel_vec(topk, out, input, d, num_blocks, block_size)
    integer, intent(in) :: topk, d, num_blocks, block_size
    real(real32), intent(inout) :: out(:)
    real(real32), intent(in) :: input(:)
    integer :: block, idx, k, lane, dv
    integer(int64) :: output_base, input_base
    real(real32) :: acc(4)

    dv = d / 4
    !$omp target teams distribute parallel do collapse(2) num_threads(block_size) private(output_base, input_base, k, lane, acc)
    do block = 0, num_blocks - 1
      do idx = 0, dv - 1
        output_base = int(block, int64) * int(d, int64)
        input_base = output_base * int(topk, int64)
        acc = 0.0_real32
        do k = 0, topk - 1
          do lane = 0, 3
            acc(lane + 1) = acc(lane + 1) + input(input_base + int(k * d + idx * 4 + lane, int64) + 1)
          end do
        end do
        do lane = 0, 3
          out(output_base + int(idx * 4 + lane, int64) + 1) = acc(lane + 1)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine moe_sum_kernel_vec

end program main
