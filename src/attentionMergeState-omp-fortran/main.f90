! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg
  integer :: num_tokens, num_heads, head_size, repeat

  if (command_argument_count() /= 4) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <number of tokens> <number of heads> <head size> '
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) num_tokens
  call get_command_argument(2, arg)
  read(arg, *) num_heads
  call get_command_argument(3, arg)
  read(arg, *) head_size
  call get_command_argument(4, arg)
  read(arg, *) repeat
  if (num_tokens <= 0 .or. num_heads <= 0 .or. head_size <= 0 .or. repeat <= 0) stop 1

  write(*,*)
  write(*,'(A,I0,A,I0,A,I0)') '#tokens ', num_tokens, ', #heads ', num_heads, ', head dimension ', head_size
  write(*,'(A)', advance='no') 'output dtype FP32: '
  call merge_attn_states_launcher(repeat, num_tokens, num_heads, head_size)
  write(*,'(A)', advance='no') 'output dtype FP16: '
  call merge_attn_states_launcher(repeat, num_tokens, num_heads, head_size)
  write(*,'(A)') '----------------------------------------------------'

contains

  subroutine merge_attn_states_launcher(repeat, num_tokens, num_heads, head_size)
    integer, intent(in) :: repeat, num_tokens, num_heads, head_size
    integer(int64) :: output_size, lse_size
    real(real32), allocatable :: prefix_output(:), suffix_output(:), output(:), ref_output(:)
    real(real32), allocatable :: prefix_lse(:), suffix_lse(:), lse(:), ref_lse(:)
    real(real64) :: start_time, end_time
    integer :: i
    logical :: ok

    output_size = int(num_tokens, int64) * int(num_heads, int64) * int(head_size, int64)
    lse_size = int(num_tokens, int64) * int(num_heads, int64)

    allocate(prefix_output(0:output_size - 1), suffix_output(0:output_size - 1), &
             output(0:output_size - 1), ref_output(0:output_size - 1), &
             prefix_lse(0:lse_size - 1), suffix_lse(0:lse_size - 1), &
             lse(0:lse_size - 1), ref_lse(0:lse_size - 1))

    !$omp target data map(alloc: prefix_output(0:output_size - 1), suffix_output(0:output_size - 1), &
    !$omp& output(0:output_size - 1), lse(0:lse_size - 1), &
    !$omp& prefix_lse(0:lse_size - 1), suffix_lse(0:lse_size - 1))
    call uniform_fill_kernel(prefix_output, output_size, -1.0_real32 / sqrt(real(head_size, real32)), &
                             1.0_real32 / sqrt(real(head_size, real32)), 1234_int32)
    call uniform_fill_kernel(suffix_output, output_size, -1.0_real32 / sqrt(real(head_size, real32)), &
                             1.0_real32 / sqrt(real(head_size, real32)), 1234_int32)
    call uniform_fill_kernel(prefix_lse, lse_size, -1.0_real32 / sqrt(real(head_size, real32)), &
                             1.0_real32 / sqrt(real(head_size, real32)), 1234_int32)
    call uniform_fill_kernel(suffix_lse, lse_size, -1.0_real32 / sqrt(real(head_size, real32)), &
                             1.0_real32 / sqrt(real(head_size, real32)), 1234_int32)

    !$omp target update from(prefix_output(0:output_size - 1), suffix_output(0:output_size - 1), &
    !$omp& prefix_lse(0:lse_size - 1), suffix_lse(0:lse_size - 1))

    call merge_reference(ref_output, prefix_output, suffix_output, ref_lse, prefix_lse, suffix_lse, &
                         num_tokens, num_heads, head_size)
    do i = 1, 100
      call merge_kernel(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                        num_tokens, num_heads, head_size)
    end do
    !$omp target update from(output(0:output_size - 1), lse(0:lse_size - 1))

    ok = check_outputs(output, ref_output, lse, ref_lse, output_size, lse_size)

    do i = 1, 100
      call merge_kernel2(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                         num_tokens, num_heads, head_size)
    end do
    !$omp target update from(output(0:output_size - 1), lse(0:lse_size - 1))

    ok = ok .and. check_outputs(output, ref_output, lse, ref_lse, output_size, lse_size)
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    start_time = omp_get_wtime()
    do i = 1, repeat
      call merge_kernel(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                        num_tokens, num_heads, head_size)
    end do
    end_time = omp_get_wtime()
    write(*,'(A,F0.6,A)') 'Average execution time of the kernel: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'

    start_time = omp_get_wtime()
    do i = 1, repeat
      call merge_kernel2(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                         num_tokens, num_heads, head_size)
    end do
    end_time = omp_get_wtime()
    write(*,'(A,F0.6,A)') 'Average execution time of the kernel2: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp end target data

    deallocate(prefix_output, suffix_output, output, ref_output, prefix_lse, suffix_lse, lse, ref_lse)
  end subroutine merge_attn_states_launcher

  subroutine uniform_fill_kernel(data, n, low, high, seed)
    real(real32), intent(out) :: data(0:)
    integer(int64), intent(in) :: n
    real(real32), intent(in) :: low, high
    integer(int32), intent(in) :: seed
    integer(int64) :: idx
    integer(int32) :: rng
    real(real32) :: u, v

    !$omp target teams distribute parallel do private(idx, rng, u, v)
    do idx = 0_int64, n - 1_int64
      rng = ieor(seed, int(idx, int32))
      rng = xorshift32(rng)
      u = real(ishft(rng, -8), real32) * 5.9604644775390625e-8_real32
      v = low + (high - low) * u
      data(idx) = v
    end do
    !$omp end target teams distribute parallel do
  end subroutine uniform_fill_kernel

  integer(int32) function xorshift32(state)
    integer(int32), intent(inout) :: state

    state = ieor(state, ishft(state, 13))
    state = ieor(state, ishft(state, -17))
    state = ieor(state, ishft(state, 5))
    xorshift32 = state
  end function xorshift32

  subroutine merge_kernel(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                          num_tokens, num_heads, head_size)
    real(real32), intent(inout) :: output(0:), lse(0:)
    real(real32), intent(in) :: prefix_output(0:), suffix_output(0:), prefix_lse(0:), suffix_lse(0:)
    integer, intent(in) :: num_tokens, num_heads, head_size
    integer :: t, h, d, lse_idx, base
    real(real32) :: p_lse, s_lse, max_lse, p_exp, s_exp, out_se, p_scale, s_scale

    !$omp target teams distribute collapse(2) num_teams(num_tokens * num_heads) &
    !$omp& private(t, h, d, lse_idx, base, p_lse, s_lse, max_lse, p_exp, s_exp, out_se, p_scale, s_scale)
    do t = 0, num_tokens - 1
      do h = 0, num_heads - 1
        lse_idx = t * num_heads + h
        p_lse = prefix_lse(lse_idx)
        s_lse = suffix_lse(lse_idx)
        max_lse = max(p_lse, s_lse)
        p_exp = exp(p_lse - max_lse)
        s_exp = exp(s_lse - max_lse)
        out_se = p_exp + s_exp
        lse(lse_idx) = log(out_se) + max_lse
        p_scale = p_exp / out_se
        s_scale = s_exp / out_se
        base = t * num_heads * head_size + h * head_size
        !$omp parallel do private(d)
        do d = 0, head_size - 1
          output(base + d) = prefix_output(base + d) * p_scale + suffix_output(base + d) * s_scale
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine merge_kernel

  subroutine merge_kernel2(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                           num_tokens, num_heads, head_size)
    real(real32), intent(inout) :: output(0:), lse(0:)
    real(real32), intent(in) :: prefix_output(0:), suffix_output(0:), prefix_lse(0:), suffix_lse(0:)
    integer, intent(in) :: num_tokens, num_heads, head_size
    integer :: bid, lane, global_idx, token_head_idx, pack_idx, pack_offset
    integer :: pack_size, threads_per_head, token_head_threads, num_threads, grid
    integer :: i, t, h, lse_idx, base
    real(real32) :: p_lse, s_lse, max_lse, p_exp, s_exp, out_se, p_scale, s_scale

    pack_size = 16 / 4
    threads_per_head = head_size / pack_size
    token_head_threads = num_tokens * num_heads * threads_per_head
    num_threads = 128
    grid = (token_head_threads + num_threads - 1) / num_threads

    !$omp target teams distribute parallel do collapse(2) num_teams(grid) num_threads(num_threads) &
    !$omp& private(bid, lane, global_idx, token_head_idx, pack_idx, pack_offset, i, t, h, lse_idx, base, &
    !$omp& p_lse, s_lse, max_lse, p_exp, s_exp, out_se, p_scale, s_scale)
    do bid = 0, grid - 1
      do lane = 0, num_threads - 1
        global_idx = bid * num_threads + lane
        if (global_idx < token_head_threads) then
          token_head_idx = global_idx / threads_per_head
          pack_idx = mod(global_idx, threads_per_head)
          t = token_head_idx / num_heads
          h = mod(token_head_idx, num_heads)
          pack_offset = pack_idx * pack_size
          base = t * num_heads * head_size + h * head_size
          lse_idx = t * num_heads + h

          p_lse = prefix_lse(lse_idx)
          s_lse = suffix_lse(lse_idx)
          max_lse = max(p_lse, s_lse)
          p_exp = exp(p_lse - max_lse)
          s_exp = exp(s_lse - max_lse)
          out_se = p_exp + s_exp
          p_scale = p_exp / out_se
          s_scale = s_exp / out_se

          if (pack_offset < head_size) then
            do i = 0, pack_size - 1
              output(base + pack_offset + i) = prefix_output(base + pack_offset + i) * p_scale + &
                                               suffix_output(base + pack_offset + i) * s_scale
            end do
          end if

          if (pack_idx == 0) lse(lse_idx) = log(out_se) + max_lse
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine merge_kernel2

  subroutine merge_reference(output, prefix_output, suffix_output, lse, prefix_lse, suffix_lse, &
                             num_tokens, num_heads, head_size)
    real(real32), intent(out) :: output(0:), lse(0:)
    real(real32), intent(in) :: prefix_output(0:), suffix_output(0:), prefix_lse(0:), suffix_lse(0:)
    integer, intent(in) :: num_tokens, num_heads, head_size
    integer :: t, h, d, lse_idx, base
    real(real32) :: p_lse, s_lse, max_lse, p_exp, s_exp, out_se, p_scale, s_scale

    do t = 0, num_tokens - 1
      do h = 0, num_heads - 1
        lse_idx = t * num_heads + h
        p_lse = prefix_lse(lse_idx)
        s_lse = suffix_lse(lse_idx)
        max_lse = max(p_lse, s_lse)
        p_exp = exp(p_lse - max_lse)
        s_exp = exp(s_lse - max_lse)
        out_se = p_exp + s_exp
        lse(lse_idx) = log(out_se) + max_lse
        p_scale = p_exp / out_se
        s_scale = s_exp / out_se
        base = t * num_heads * head_size + h * head_size
        do d = 0, head_size - 1
          output(base + d) = prefix_output(base + d) * p_scale + suffix_output(base + d) * s_scale
        end do
      end do
    end do
  end subroutine merge_reference

  logical function check_outputs(output, ref_output, lse, ref_lse, output_size, lse_size)
    real(real32), intent(in) :: output(0:), ref_output(0:), lse(0:), ref_lse(0:)
    integer(int64), intent(in) :: output_size, lse_size
    integer(int64) :: i

    check_outputs = .true.
    do i = 0_int64, output_size - 1_int64
      if (abs(output(i) - ref_output(i)) > 1.0e-3_real32) then
        check_outputs = .false.
        return
      end if
    end do
    do i = 0_int64, lse_size - 1_int64
      if (abs(lse(i) - ref_lse(i)) > 1.0e-3_real32) then
        check_outputs = .false.
        return
      end if
    end do
  end function check_outputs

end program main
