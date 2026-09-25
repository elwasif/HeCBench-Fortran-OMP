! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use omp_lib
  implicit none

  interface
    subroutine moe_initialize_inputs(gating_output, topk_indices, topk_indices_ref, num_tokens, num_experts, topk) &
        bind(C, name="moe_initialize_inputs")
      import :: c_float, c_int
      real(c_float), intent(out) :: gating_output(*)
      integer(c_int), intent(out) :: topk_indices(*), topk_indices_ref(*)
      integer(c_int), value :: num_tokens, num_experts, topk
    end subroutine moe_initialize_inputs
  end interface

  character(len=256) :: arg0, arg
  integer :: num_tokens, num_experts, topk, repeat
  integer(int64) :: output_size, index_size
  real(real32), allocatable :: gating_output(:), softmax_workspace(:)
  real(real32), allocatable :: topk_weights(:), topk_weights_ref(:)
  integer(int32), allocatable :: topk_indices(:), topk_indices_ref(:)
  integer(int32), allocatable :: token_expert_indices(:), token_expert_indices_ref(:)
  integer :: i
  logical :: ok
  real(real64) :: start_time, elapsed_us

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of tokens> <number of experts> <top K> <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg); read(arg, *) num_tokens
  call get_command_argument(2, arg); read(arg, *) num_experts
  call get_command_argument(3, arg); read(arg, *) topk
  call get_command_argument(4, arg); read(arg, *) repeat

  output_size = int(num_tokens, int64) * int(num_experts, int64)
  index_size = int(num_tokens, int64) * int(topk, int64)
  allocate(gating_output(output_size), softmax_workspace(output_size))
  allocate(topk_weights(index_size), topk_weights_ref(index_size))
  allocate(topk_indices(index_size), topk_indices_ref(index_size))
  allocate(token_expert_indices(index_size), token_expert_indices_ref(index_size))

  call initialize_inputs(gating_output, topk_indices, topk_indices_ref)
  call softmax_ref(gating_output, softmax_workspace, num_tokens, num_experts)
  call topk_ref(softmax_workspace, topk_weights_ref, topk_indices_ref, token_expert_indices_ref, num_tokens, num_experts, topk)

  !$omp target data map(to: gating_output(1:output_size), topk_indices(1:index_size)) &
  !$omp& map(alloc: softmax_workspace(1:output_size), topk_weights(1:index_size), token_expert_indices(1:index_size))
  call softmax_kernel(gating_output, softmax_workspace, num_tokens, num_experts)
  call topk_kernel(softmax_workspace, topk_weights, topk_indices, token_expert_indices, num_tokens, num_experts, topk)
  !$omp target update from(topk_weights(1:index_size), topk_indices(1:index_size), token_expert_indices(1:index_size))

  ok = all(topk_indices == topk_indices_ref) .and. all(token_expert_indices == token_expert_indices_ref)
  do i = 1, int(index_size)
    if (abs(topk_weights(i) - topk_weights_ref(i)) > 1.0e-3_real32) ok = .false.
  end do
  write(*,'(A)') merge('PASS', 'FAIL', ok)

  start_time = omp_get_wtime()
  do i = 1, repeat
    call softmax_kernel(gating_output, softmax_workspace, num_tokens, num_experts)
    call topk_kernel(softmax_workspace, topk_weights, topk_indices, token_expert_indices, num_tokens, num_experts, topk)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average execution time of kernels: ', elapsed_us, ' (us)'
  !$omp end target data

  deallocate(gating_output, softmax_workspace, topk_weights, topk_weights_ref)
  deallocate(topk_indices, topk_indices_ref, token_expert_indices, token_expert_indices_ref)

contains

  subroutine initialize_inputs(gating_output, topk_indices, topk_indices_ref)
    real(real32), intent(out) :: gating_output(:)
    integer(int32), intent(out) :: topk_indices(:), topk_indices_ref(:)

    call moe_initialize_inputs(gating_output, topk_indices, topk_indices_ref, &
                               int(num_tokens, c_int), int(num_experts, c_int), int(topk, c_int))
  end subroutine initialize_inputs

  subroutine softmax_ref(input, output, num_tokens, num_experts)
    real(real32), intent(in) :: input(:)
    real(real32), intent(out) :: output(:)
    integer, intent(in) :: num_tokens, num_experts
    integer :: token, expert, base
    real(real32) :: max_val, sum_val

    do token = 0, num_tokens - 1
      base = token * num_experts
      max_val = -huge(1.0_real32)
      do expert = 0, num_experts - 1
        max_val = max(max_val, input(base + expert + 1))
      end do
      sum_val = 0.0_real32
      do expert = 0, num_experts - 1
        sum_val = sum_val + exp(input(base + expert + 1) - max_val)
      end do
      do expert = 0, num_experts - 1
        output(base + expert + 1) = exp(input(base + expert + 1) - max_val) / sum_val
      end do
    end do
  end subroutine softmax_ref

  subroutine topk_ref(input, weights, indices, source_rows, num_tokens, num_experts, topk)
    real(real32), intent(in) :: input(:)
    real(real32), intent(out) :: weights(:)
    integer(int32), intent(inout) :: indices(:)
    integer(int32), intent(out) :: source_rows(:)
    integer, intent(in) :: num_tokens, num_experts, topk
    integer :: token, k_idx, expert, prior_k, best_key, idx, base
    real(real32) :: best_val, val
    logical :: used

    do token = 0, num_tokens - 1
      base = token * num_experts
      do k_idx = 0, topk - 1
        best_key = 0
        best_val = -1.0_real32
        do expert = 0, num_experts - 1
          used = .false.
          do prior_k = 0, k_idx - 1
            if (indices(token * topk + prior_k + 1) == expert) used = .true.
          end do
          if (.not. used) then
            val = input(base + expert + 1)
            if (val > best_val .or. (val == best_val .and. expert < best_key)) then
              best_val = val
              best_key = expert
            end if
          end if
        end do
        idx = token * topk + k_idx + 1
        weights(idx) = best_val
        indices(idx) = int(best_key, int32)
        source_rows(idx) = int(k_idx * num_tokens + token, int32)
      end do
    end do
  end subroutine topk_ref

  subroutine softmax_kernel(input, output, num_tokens, num_experts)
    real(real32), intent(in) :: input(:)
    real(real32), intent(inout) :: output(:)
    integer, intent(in) :: num_tokens, num_experts
    integer :: token, expert, base
    real(real32) :: max_val, sum_val

    !$omp target teams distribute parallel do private(expert, base, max_val, sum_val) num_teams(num_tokens) num_threads(256)
    do token = 0, num_tokens - 1
      base = token * num_experts
      max_val = -huge(1.0_real32)
      !$omp parallel do reduction(max:max_val)
      do expert = 0, num_experts - 1
        max_val = max(max_val, input(base + expert + 1))
      end do
      sum_val = 0.0_real32
      !$omp parallel do reduction(+:sum_val)
      do expert = 0, num_experts - 1
        sum_val = sum_val + exp(input(base + expert + 1) - max_val)
      end do
      !$omp parallel do
      do expert = 0, num_experts - 1
        output(base + expert + 1) = exp(input(base + expert + 1) - max_val) / sum_val
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine softmax_kernel

  subroutine topk_kernel(input, weights, indices, source_rows, num_tokens, num_experts, topk)
    real(real32), intent(in) :: input(:)
    real(real32), intent(inout) :: weights(:)
    integer(int32), intent(inout) :: indices(:), source_rows(:)
    integer, intent(in) :: num_tokens, num_experts, topk
    integer :: token, k_idx, expert, prior_k, best_key, idx, base
    real(real32) :: best_val, val
    logical :: used

    !$omp target teams distribute parallel do private(k_idx, expert, prior_k, best_key, idx, base, best_val, val, used) num_threads(256)
    do token = 0, num_tokens - 1
      base = token * num_experts
      do k_idx = 0, topk - 1
        best_key = 0
        best_val = -1.0_real32
        do expert = 0, num_experts - 1
          used = .false.
          do prior_k = 0, k_idx - 1
            if (indices(token * topk + prior_k + 1) == expert) used = .true.
          end do
          if (.not. used) then
            val = input(base + expert + 1)
            if (val > best_val .or. (val == best_val .and. expert < best_key)) then
              best_val = val
              best_key = expert
            end if
          end if
        end do
        idx = token * topk + k_idx + 1
        weights(idx) = best_val
        indices(idx) = int(best_key, int32)
        source_rows(idx) = int(k_idx * num_tokens + token, int32)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine topk_kernel

end program main
