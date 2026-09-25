! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_size_t
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand

    subroutine fill_logprob_logits(logits, logits_size) bind(C, name='fill_logprob_logits')
      import :: c_float, c_size_t
      real(c_float), intent(out) :: logits(*)
      integer(c_size_t), value :: logits_size
    end subroutine fill_logprob_logits
  end interface

  character(len=256) :: arg0, arg
  integer :: max_length, batch_size, vocab_size, repeat
  integer :: vocab_size_padded, block_size, i
  integer(int64) :: logits_size, log_probs_size, ids_size
  real(real32), allocatable :: logits(:), log_probs(:), log_probs_ref(:)
  real(real32), allocatable :: cum_log_probs(:), cum_log_probs_ref(:)
  integer, allocatable :: ids(:), lengths(:)
  real(real64) :: start_time, end_time
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <max_seq_len> <batch_size> <vocab_size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) max_length
  call get_command_argument(2, arg)
  read(arg, *) batch_size
  call get_command_argument(3, arg)
  read(arg, *) vocab_size
  call get_command_argument(4, arg)
  read(arg, *) repeat

  if (max_length <= 1 .or. batch_size <= 0 .or. vocab_size <= 0 .or. repeat <= 0) stop 1

  vocab_size_padded = ((vocab_size + 31) / 32) * 32
  logits_size = int(batch_size, int64) * int(max_length, int64) * int(vocab_size_padded, int64)
  log_probs_size = int(batch_size, int64) * int(max_length - 1, int64)
  ids_size = int(batch_size, int64) * int(max_length, int64)

  allocate(logits(0:logits_size - 1), log_probs(0:log_probs_size - 1), &
           log_probs_ref(0:log_probs_size - 1), &
           cum_log_probs(0:batch_size - 1), cum_log_probs_ref(0:batch_size - 1), &
           ids(0:ids_size - 1), lengths(0:batch_size - 1))

  call initialize_inputs(logits, ids, lengths, logits_size, ids_size, batch_size, max_length, vocab_size)
  log_probs = 0.0_real32
  log_probs_ref = 0.0_real32
  cum_log_probs = 0.0_real32
  cum_log_probs_ref = 0.0_real32

  if (vocab_size < 1024) then
    block_size = ((vocab_size + 31) / 32) * 32
  else
    block_size = 1024
  end if

  !$omp target data map(to: logits(0:logits_size - 1), ids(0:ids_size - 1), lengths(0:batch_size - 1)) &
  !$omp& map(alloc: log_probs(0:log_probs_size - 1), cum_log_probs(0:batch_size - 1))
  call log_probs_kernel(log_probs, logits, ids, lengths, max_length, batch_size, &
                        vocab_size, vocab_size_padded, block_size)
  call accumulate_log_probs(cum_log_probs, log_probs, lengths, max_length, batch_size, block_size)
  !$omp target update from(log_probs(0:log_probs_size - 1), cum_log_probs(0:batch_size - 1))

  call log_probs_cpu(log_probs_ref, logits, ids, lengths, max_length, batch_size, &
                     vocab_size, vocab_size_padded, cum_log_probs_ref)

  ok = .true.
  do i = 0, int(log_probs_size - 1_int64)
    if (abs(log_probs(i) - log_probs_ref(i)) > 1.0e-3_real32) then
      write(*,'(A,I0,A,F0.6,A,F0.6)') 'log_probs mismatch @', i, ': ', log_probs(i), ' != ', log_probs_ref(i)
      ok = .false.
      exit
    end if
  end do
  do i = 0, batch_size - 1
    if (abs(cum_log_probs(i) - cum_log_probs_ref(i)) > 1.0e-1_real32) then
      write(*,'(A,I0,A,F0.6,A,F0.6)') 'cum_log_probs mismatch @', i, ': ', cum_log_probs(i), ' != ', cum_log_probs_ref(i)
      ok = .false.
    end if
  end do
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  start_time = omp_get_wtime()
  do i = 1, repeat
    call log_probs_kernel(log_probs, logits, ids, lengths, max_length, batch_size, &
                          vocab_size, vocab_size_padded, block_size)
    call accumulate_log_probs(cum_log_probs, log_probs, lengths, max_length, batch_size, block_size)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time of kernels: ', &
    (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
  !$omp end target data

  deallocate(logits, log_probs, log_probs_ref, cum_log_probs, cum_log_probs_ref, ids, lengths)

contains

  subroutine initialize_inputs(logits, ids, lengths, logits_size, ids_size, batch_size, max_length, vocab_size)
    real(real32), intent(out) :: logits(0:)
    integer, intent(out) :: ids(0:), lengths(0:)
    integer(int64), intent(in) :: logits_size, ids_size
    integer, intent(in) :: batch_size, max_length, vocab_size
    integer(int64) :: i
    integer :: b

    call fill_logprob_logits(logits, int(logits_size, c_size_t))
    do b = 0, batch_size - 1
      lengths(b) = max_length
    end do
    call c_srand(123_c_int)
    do i = 0_int64, ids_size - 1_int64
      ids(i) = mod(c_rand(), vocab_size)
    end do
  end subroutine initialize_inputs

  subroutine log_probs_kernel(log_probs, logits, ids, lengths, max_input_length, batch_size, &
                              vocab_size, vocab_size_padded, block_size)
    real(real32), intent(out) :: log_probs(0:)
    real(real32), intent(in) :: logits(0:)
    integer, intent(in) :: ids(0:), lengths(0:)
    integer, intent(in) :: max_input_length, batch_size, vocab_size, vocab_size_padded, block_size
    integer :: bidx, step, j, idx, token_idx
    integer(int64) :: logits_base
    real(real32) :: max_val, sum_exp, val

    !$omp target teams distribute collapse(2) num_teams(batch_size * (max_input_length - 1)) &
    !$omp& thread_limit(block_size) private(bidx, step, j, idx, token_idx, logits_base, max_val, sum_exp, val)
    do bidx = 0, batch_size - 1
      do step = 0, max_input_length - 2
        if (step < lengths(bidx) - 1) then
          logits_base = int(bidx, int64) * int(max_input_length, int64) * int(vocab_size_padded, int64) + &
                        int(step, int64) * int(vocab_size_padded, int64)
          max_val = -huge(1.0_real32)
          !$omp parallel do reduction(max:max_val) num_threads(block_size) private(j, val)
          do j = 0, vocab_size - 1
            val = logits(logits_base + int(j, int64))
            if (val > max_val) max_val = val
          end do
          !$omp end parallel do

          sum_exp = 0.0_real32
          !$omp parallel do reduction(+:sum_exp) num_threads(block_size) private(j)
          do j = 0, vocab_size - 1
            sum_exp = sum_exp + exp(logits(logits_base + int(j, int64)) - max_val)
          end do
          !$omp end parallel do

          idx = step + bidx * (max_input_length - 1)
          token_idx = step + 1 + bidx * max_input_length
          log_probs(idx) = logits(logits_base + int(ids(token_idx), int64)) - max_val - log(sum_exp + 1.0e-9_real32)
        end if
      end do
    end do
    !$omp end target teams distribute
  end subroutine log_probs_kernel

  subroutine accumulate_log_probs(cum_log_probs, log_probs, lengths, max_input_length, batch_size, block_size)
    real(real32), intent(out) :: cum_log_probs(0:)
    real(real32), intent(in) :: log_probs(0:)
    integer, intent(in) :: lengths(0:)
    integer, intent(in) :: max_input_length, batch_size, block_size
    integer :: bidx, step
    real(real32) :: accum

    !$omp target teams distribute num_teams(batch_size) private(bidx, step, accum)
    do bidx = 0, batch_size - 1
      accum = 0.0_real32
      !$omp parallel do reduction(+:accum) num_threads(block_size) private(step)
      do step = 0, lengths(bidx) - 2
        accum = accum + log_probs(step + bidx * (max_input_length - 1))
      end do
      !$omp end parallel do
      cum_log_probs(bidx) = accum
    end do
    !$omp end target teams distribute
  end subroutine accumulate_log_probs

  subroutine log_probs_cpu(log_probs, logits, ids, lengths, max_input_length, batch_size, &
                           vocab_size, vocab_size_padded, cum_log_probs)
    real(real32), intent(out) :: log_probs(0:), cum_log_probs(0:)
    real(real32), intent(in) :: logits(0:)
    integer, intent(in) :: ids(0:), lengths(0:)
    integer, intent(in) :: max_input_length, batch_size, vocab_size, vocab_size_padded
    integer :: bidx, step, j, idx, token_idx
    integer(int64) :: logits_base
    real(real32) :: max_val, sum_exp, accum, val

    do bidx = 0, batch_size - 1
      accum = 0.0_real32
      do step = 0, lengths(bidx) - 2
        logits_base = int(bidx, int64) * int(max_input_length, int64) * int(vocab_size_padded, int64) + &
                      int(step, int64) * int(vocab_size_padded, int64)
        max_val = -huge(1.0_real32)
        do j = 0, vocab_size - 1
          val = logits(logits_base + int(j, int64))
          if (val > max_val) max_val = val
        end do

        sum_exp = 0.0_real32
        do j = 0, vocab_size - 1
          sum_exp = sum_exp + exp(logits(logits_base + int(j, int64)) - max_val)
        end do

        idx = step + bidx * (max_input_length - 1)
        token_idx = step + 1 + bidx * max_input_length
        log_probs(idx) = logits(logits_base + int(ids(token_idx), int64)) - max_val - log(sum_exp + 1.0e-9_real32)
        accum = accum + log_probs(idx)
      end do
      cum_log_probs(bidx) = accum
    end do
  end subroutine log_probs_cpu

end program main
