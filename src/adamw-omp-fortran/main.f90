! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_float, c_int8_t, c_int64_t
  use omp_lib
  implicit none

  integer, parameter :: block_size = 64
  real(real32), parameter :: lr = 1.0e-3_real32
  real(real32), parameter :: weight_decay = 1.0e-2_real32
  real(real32), parameter :: beta1 = 0.9_real32
  real(real32), parameter :: beta2 = 0.999_real32
  real(real32), parameter :: eps = 1.0e-8_real32
  real(real32), parameter :: resid_beta1 = 1.0_real32 - beta1
  real(real32), parameter :: resid_beta2 = 1.0_real32 - beta2
  real(real32), parameter :: weight_decay_update = 1.0_real32 - lr * weight_decay
  real(real32), parameter :: exp_qmap(0:15) = [ &
    -0.8875_real32, -0.6625_real32, -0.4375_real32, -0.2125_real32, &
    -0.0775_real32, -0.0325_real32, -0.0055_real32,  0.0000_real32, &
     0.0055_real32,  0.0325_real32,  0.0775_real32,  0.2125_real32, &
     0.4375_real32,  0.6625_real32,  0.8875_real32,  1.0000_real32]
  real(real32), parameter :: exp_qmidpt(0:14) = [ &
    -0.775_real32, -0.55_real32, -0.325_real32, -0.145_real32, -0.055_real32, &
    -0.019_real32, -0.00275_real32, 0.00275_real32, 0.019_real32, 0.055_real32, &
     0.145_real32, 0.325_real32, 0.55_real32, 0.775_real32, 0.94375_real32]
  real(real32), parameter :: sq_qmap(0:15) = [ &
    0.0625_real32, 0.1250_real32, 0.1875_real32, 0.2500_real32, &
    0.3125_real32, 0.3750_real32, 0.4375_real32, 0.5000_real32, &
    0.5625_real32, 0.6250_real32, 0.6875_real32, 0.7500_real32, &
    0.8125_real32, 0.8750_real32, 0.9375_real32, 1.0000_real32]
  real(real32), parameter :: sq_qmidpt(0:14) = [ &
    0.09375_real32, 0.15625_real32, 0.21875_real32, 0.28125_real32, &
    0.34375_real32, 0.40625_real32, 0.46875_real32, 0.53125_real32, &
    0.59375_real32, 0.65625_real32, 0.71875_real32, 0.78125_real32, &
    0.84375_real32, 0.90625_real32, 0.96875_real32]

  interface
    subroutine adamw_initialize_inputs_c(vector_size, float_size, g, p, p_ref, m_qscale, v_qscale, &
        m_qscale_ref, v_qscale_ref, m, v, m_ref, v_ref) bind(C, name='adamw_initialize_inputs')
      import :: c_float, c_int8_t, c_int64_t
      integer(c_int64_t), value :: vector_size, float_size
      real(c_float) :: g(0:*), p(0:*), p_ref(0:*), m_qscale(0:*), v_qscale(0:*)
      real(c_float) :: m_qscale_ref(0:*), v_qscale_ref(0:*)
      integer(c_int8_t) :: m(0:*), v(0:*), m_ref(0:*), v_ref(0:*)
    end subroutine adamw_initialize_inputs_c
  end interface

  character(len=256) :: arg0, arg
  integer(int64) :: vector_size, float_size
  integer :: time_step
  integer :: step
  real(real32), allocatable :: g(:), p(:), p_ref(:), m_qscale(:), v_qscale(:)
  real(real32), allocatable :: m_qscale_ref(:), v_qscale_ref(:)
  integer(int8), allocatable :: m(:), v(:), m_ref(:), v_ref(:)
  real(real32) :: correction1, correction2_sqrt, step_size
  real(real32) :: absmax_error
  real(real64) :: start_time, end_time, avg_ms

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <vector size> <number of time steps>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) vector_size
  call get_command_argument(2, arg)
  read(arg, *) time_step
  if (vector_size <= 0_int64 .or. time_step <= 0) stop 1

  float_size = 2_int64 * vector_size
  allocate(g(0:float_size - 1), p(0:float_size - 1), p_ref(0:float_size - 1))
  allocate(m_qscale(0:float_size - 1), v_qscale(0:float_size - 1))
  allocate(m_qscale_ref(0:float_size - 1), v_qscale_ref(0:float_size - 1))
  allocate(m(0:vector_size - 1), v(0:vector_size - 1))
  allocate(m_ref(0:vector_size - 1), v_ref(0:vector_size - 1))

  call initialize_inputs(vector_size, float_size, g, p, p_ref, m_qscale, v_qscale, &
    m_qscale_ref, v_qscale_ref, m, v, m_ref, v_ref)

  !$omp target data map(to: g(0:float_size - 1), p(0:float_size - 1), &
  !$omp& m_qscale(0:float_size - 1), v_qscale(0:float_size - 1), &
  !$omp& m(0:vector_size - 1), v(0:vector_size - 1))
  do step = 1, time_step
    call step_constants(step, correction1, correction2_sqrt, step_size)
    call fused_4bit_kernel(p, g, m_qscale, v_qscale, m, v, vector_size, &
      correction2_sqrt, step_size)
    call reference_kernel(p_ref, g, m_qscale_ref, v_qscale_ref, m_ref, v_ref, vector_size, &
      correction2_sqrt, step_size)
  end do
  !$omp target update from(p(0:float_size - 1))

  absmax_error = maxval(abs(p - p_ref))
  if (absmax_error < 10.0_real32) then
    write(*,'(A,F8.6)') 'Absolute maximum error: ', absmax_error
  else
    write(*,'(A,F0.6)') 'Absolute maximum error: ', absmax_error
  end if
  if (absmax_error > 1.0e-3_real32) then
    write(*,'(A)') 'FAIL'
  else
    write(*,'(A)') 'PASS'
  end if

  start_time = omp_get_wtime()
  do step = 1, time_step
    call step_constants(step, correction1, correction2_sqrt, step_size)
    call fused_4bit_kernel(p, g, m_qscale, v_qscale, m, v, vector_size, &
      correction2_sqrt, step_size)
  end do
  end_time = omp_get_wtime()
  avg_ms = (end_time - start_time) * 1.0e3_real64 / real(time_step, real64)
  if (avg_ms < 10.0_real64) then
    write(*,'(A,F8.6,A)') 'Average kernel execution time ', avg_ms, ' (ms)'
  else
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', avg_ms, ' (ms)'
  end if
  !$omp end target data

  deallocate(g, p, p_ref, m_qscale, v_qscale, m_qscale_ref, v_qscale_ref, m, v, m_ref, v_ref)

contains

  subroutine initialize_inputs(vector_size, float_size, g, p, p_ref, m_qscale, v_qscale, &
      m_qscale_ref, v_qscale_ref, m, v, m_ref, v_ref)
    integer(int64), intent(in) :: vector_size, float_size
    real(real32), intent(out) :: g(0:), p(0:), p_ref(0:), m_qscale(0:), v_qscale(0:)
    real(real32), intent(out) :: m_qscale_ref(0:), v_qscale_ref(0:)
    integer(int8), intent(out) :: m(0:), v(0:), m_ref(0:), v_ref(0:)

    call adamw_initialize_inputs_c(vector_size, float_size, g, p, p_ref, m_qscale, v_qscale, &
      m_qscale_ref, v_qscale_ref, m, v, m_ref, v_ref)
  end subroutine initialize_inputs

  pure integer(int8) function to_int8(value)
    integer(int32), intent(in) :: value
    integer(int32) :: wrapped

    wrapped = iand(value, 255_int32)
    if (wrapped > 127_int32) wrapped = wrapped - 256_int32
    to_int8 = int(wrapped, int8)
  end function to_int8

  pure integer(int32) function to_unsigned8(value)
    integer(int8), intent(in) :: value

    to_unsigned8 = int(value, int32)
    if (to_unsigned8 < 0_int32) to_unsigned8 = to_unsigned8 + 256_int32
  end function to_unsigned8

  subroutine step_constants(step, correction1, correction2_sqrt, step_size)
    integer, intent(in) :: step
    real(real32), intent(out) :: correction1, correction2_sqrt, step_size

    correction1 = 1.0_real32 - beta1 ** step
    correction2_sqrt = sqrt(1.0_real32 - beta2 ** step)
    step_size = lr / correction1
  end subroutine step_constants

  subroutine fused_4bit_kernel(p, g, exp_qscale, sq_qscale, exp_state, sq_state, total_size, &
      correction2_sqrt, step_size)
    real(real32), intent(inout) :: p(0:), exp_qscale(0:), sq_qscale(0:)
    real(real32), intent(in) :: g(0:), correction2_sqrt, step_size
    integer(int8), intent(inout) :: exp_state(0:), sq_state(0:)
    integer(int64), intent(in) :: total_size
    integer(int64) :: block_id, global_id, num_blocks
    integer :: thread_id, exp_full, sq_full, exp_left_index, sq_left_index
    integer :: exp_right_index, sq_right_index, q_exp_left, q_sq_left, q_exp_right, q_sq_right
    integer(int32) :: packed_exp, packed_sq
    real(real32) :: absmax_exp, absmax_sq, exp_avg_qscale
    real(real32) :: p_left, p_right, g_left, g_right
    real(real32) :: exp_left, exp_right, sq_left, sq_right
    real(real32) :: local_exp_left(0:block_size - 1), local_exp_right(0:block_size - 1)
    real(real32) :: local_sq_left(0:block_size - 1), local_sq_right(0:block_size - 1)
    integer :: low, high, mid
    real(real32) :: mapped_value

    num_blocks = (total_size + block_size - 1_int64) / block_size
    !$omp target teams distribute num_teams(num_blocks) &
    !$omp& private(block_id, thread_id, global_id, exp_full, sq_full, exp_left_index, sq_left_index, &
    !$omp& exp_right_index, sq_right_index, q_exp_left, q_sq_left, q_exp_right, packed_exp, packed_sq, &
    !$omp& absmax_exp, absmax_sq, exp_avg_qscale, p_left, p_right, g_left, g_right, exp_left, exp_right, &
    !$omp& sq_left, sq_right, local_exp_left, local_exp_right, local_sq_left, local_sq_right, &
    !$omp& low, high, mid, mapped_value)
    do block_id = 0_int64, num_blocks - 1_int64
      absmax_exp = 0.0_real32
      absmax_sq = 0.0_real32
      local_exp_left = 0.0_real32
      local_exp_right = 0.0_real32
      local_sq_left = 0.0_real32
      local_sq_right = 0.0_real32

      !$omp parallel do reduction(max:absmax_sq, absmax_exp) num_threads(block_size) &
      !$omp& private(thread_id, global_id, exp_full, sq_full, exp_left_index, sq_left_index, &
      !$omp& exp_right_index, sq_right_index, exp_avg_qscale, p_left, p_right, g_left, g_right, &
      !$omp& exp_left, exp_right, sq_left, sq_right)
      do thread_id = 0, block_size - 1
        global_id = block_id * block_size + int(thread_id, int64)
        if (global_id < total_size) then
          exp_full = to_unsigned8(exp_state(global_id))
          sq_full = to_unsigned8(sq_state(global_id))
          exp_left_index = iand(exp_full, 15)
          sq_left_index = iand(sq_full, 15)
          exp_right_index = iand(ishft(exp_full, -4), 15)
          sq_right_index = iand(ishft(sq_full, -4), 15)

          p_left = p(2_int64 * global_id) * weight_decay_update
          p_right = p(2_int64 * global_id + 1_int64) * weight_decay_update
          g_left = g(2_int64 * global_id)
          g_right = g(2_int64 * global_id + 1_int64)
          exp_avg_qscale = exp_qscale(block_id)

          exp_left = exp_qmap(exp_left_index) * exp_avg_qscale
          exp_left = beta1 * exp_left + resid_beta1 * g_left
          sq_left = sq_qmap(sq_left_index) * sq_qscale(block_id)
          sq_left = beta2 * sq_left + resid_beta2 * (g_left * g_left)
          p(2_int64 * global_id) = p_left - (step_size * (exp_left / (sqrt(sq_left) / correction2_sqrt + eps)))

          exp_right = exp_qmap(exp_right_index) * exp_avg_qscale
          exp_right = beta1 * exp_right + resid_beta1 * g_right
          sq_right = sq_qmap(sq_right_index) * sq_qscale(block_id)
          sq_right = beta2 * sq_right + resid_beta2 * (g_right * g_right)
          p(2_int64 * global_id + 1_int64) = p_right - (step_size * (exp_right / (sqrt(sq_right) / correction2_sqrt + eps)))

          local_exp_left(thread_id) = exp_left
          local_exp_right(thread_id) = exp_right
          local_sq_left(thread_id) = sq_left
          local_sq_right(thread_id) = sq_right
          absmax_exp = max(absmax_exp, max(exp_left, exp_right))
          absmax_sq = max(absmax_sq, max(sq_left, sq_right))
        end if
      end do
      !$omp end parallel do

      exp_qscale(block_id) = absmax_exp
      sq_qscale(block_id) = absmax_sq

      !$omp parallel do num_threads(block_size) &
      !$omp& private(thread_id, global_id, q_exp_left, q_sq_left, q_exp_right, q_sq_right, &
      !$omp& packed_exp, packed_sq, low, high, mid, mapped_value)
      do thread_id = 0, block_size - 1
        global_id = block_id * block_size + int(thread_id, int64)
        if (global_id < total_size) then
          mapped_value = local_exp_left(thread_id) / absmax_exp
          q_exp_left = q_index_exp(mapped_value, low, high, mid)
          mapped_value = local_sq_left(thread_id) / absmax_sq
          q_sq_left = q_index_sq(mapped_value, low, high, mid)
          mapped_value = local_exp_right(thread_id) / absmax_exp
          q_exp_right = q_index_exp(mapped_value, low, high, mid)
          mapped_value = local_sq_right(thread_id) / absmax_sq
          q_sq_right = q_index_sq(mapped_value, low, high, mid)

          packed_exp = ior(iand(q_exp_left, 15), ishft(iand(q_exp_right, 15), 4))
          packed_sq = ior(iand(q_sq_left, 15), ishft(iand(q_sq_right, 15), 4))
          exp_state(global_id) = to_int8(packed_exp)
          sq_state(global_id) = to_int8(packed_sq)
        end if
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine fused_4bit_kernel

  subroutine reference_kernel(p, g, exp_qscale, sq_qscale, exp_state, sq_state, total_size, &
      correction2_sqrt, step_size)
    real(real32), intent(inout) :: p(0:), exp_qscale(0:), sq_qscale(0:)
    real(real32), intent(in) :: g(0:), correction2_sqrt, step_size
    integer(int8), intent(inout) :: exp_state(0:), sq_state(0:)
    integer(int64), intent(in) :: total_size
    integer(int64) :: block_id, global_id, num_blocks
    integer :: thread_id, exp_full, sq_full, exp_left_index, sq_left_index
    integer :: exp_right_index, sq_right_index, q_exp_left, q_sq_left, q_exp_right, q_sq_right
    integer(int32) :: packed_exp, packed_sq
    real(real32) :: absmax_exp, absmax_sq, exp_avg_qscale
    real(real32) :: p_left, p_right, g_left, g_right
    real(real32) :: exp_left, exp_right, sq_left, sq_right
    real(real32) :: local_exp_left(0:block_size - 1), local_exp_right(0:block_size - 1)
    real(real32) :: local_sq_left(0:block_size - 1), local_sq_right(0:block_size - 1)

    num_blocks = (total_size + block_size - 1_int64) / block_size
    do block_id = 0_int64, num_blocks - 1_int64
      absmax_exp = 0.0_real32
      absmax_sq = 0.0_real32
      local_exp_left = 0.0_real32
      local_exp_right = 0.0_real32
      local_sq_left = 0.0_real32
      local_sq_right = 0.0_real32

      do thread_id = 0, block_size - 1
        global_id = block_id * block_size + int(thread_id, int64)
        if (global_id >= total_size) exit
        exp_full = to_unsigned8(exp_state(global_id))
        sq_full = to_unsigned8(sq_state(global_id))
        exp_left_index = iand(exp_full, 15)
        sq_left_index = iand(sq_full, 15)
        exp_right_index = iand(ishft(exp_full, -4), 15)
        sq_right_index = iand(ishft(sq_full, -4), 15)

        p_left = p(2_int64 * global_id) * weight_decay_update
        p_right = p(2_int64 * global_id + 1_int64) * weight_decay_update
        g_left = g(2_int64 * global_id)
        g_right = g(2_int64 * global_id + 1_int64)
        exp_avg_qscale = exp_qscale(block_id)

        exp_left = exp_qmap(exp_left_index) * exp_avg_qscale
        exp_left = beta1 * exp_left + resid_beta1 * g_left
        sq_left = sq_qmap(sq_left_index) * sq_qscale(block_id)
        sq_left = beta2 * sq_left + resid_beta2 * (g_left * g_left)
        p(2_int64 * global_id) = p_left - (step_size * (exp_left / (sqrt(sq_left) / correction2_sqrt + eps)))

        exp_right = exp_qmap(exp_right_index) * exp_avg_qscale
        exp_right = beta1 * exp_right + resid_beta1 * g_right
        sq_right = sq_qmap(sq_right_index) * sq_qscale(block_id)
        sq_right = beta2 * sq_right + resid_beta2 * (g_right * g_right)
        p(2_int64 * global_id + 1_int64) = p_right - (step_size * (exp_right / (sqrt(sq_right) / correction2_sqrt + eps)))

        local_exp_left(thread_id) = exp_left
        local_exp_right(thread_id) = exp_right
        local_sq_left(thread_id) = sq_left
        local_sq_right(thread_id) = sq_right
        absmax_exp = max(absmax_exp, max(exp_left, exp_right))
        absmax_sq = max(absmax_sq, max(sq_left, sq_right))
      end do

      exp_qscale(block_id) = absmax_exp
      sq_qscale(block_id) = absmax_sq

      do thread_id = 0, block_size - 1
        global_id = block_id * block_size + int(thread_id, int64)
        if (global_id >= total_size) exit
        q_exp_left = q_index_exp(local_exp_left(thread_id) / absmax_exp)
        q_sq_left = q_index_sq(local_sq_left(thread_id) / absmax_sq)
        q_exp_right = q_index_exp(local_exp_right(thread_id) / absmax_exp)
        q_sq_right = q_index_sq(local_sq_right(thread_id) / absmax_sq)
        packed_exp = ior(iand(q_exp_left, 15), ishft(iand(q_exp_right, 15), 4))
        packed_sq = ior(iand(q_sq_left, 15), ishft(iand(q_sq_right, 15), 4))
        exp_state(global_id) = to_int8(packed_exp)
        sq_state(global_id) = to_int8(packed_sq)
      end do
    end do
  end subroutine reference_kernel

  pure integer function q_index_exp(x, low_arg, high_arg, mid_arg)
    real(real32), intent(in) :: x
    integer, intent(out), optional :: low_arg, high_arg, mid_arg
    integer :: low, high, mid

    if (x <= exp_qmap(0)) then
      q_index_exp = 0
      return
    end if
    if (exp_qmap(15) <= x) then
      q_index_exp = 15
      return
    end if
    low = 0
    high = 15
    do while (low < high)
      mid = ishft(low + high, -1)
      if (exp_qmap(mid) <= x) then
        low = mid + 1
      else
        high = mid
      end if
    end do
    if (present(low_arg)) low_arg = low
    if (present(high_arg)) high_arg = high
    if (present(mid_arg)) mid_arg = mid
    if (exp_qmidpt(low - 1) < x) then
      q_index_exp = low
    else
      q_index_exp = low - 1
    end if
  end function q_index_exp

  pure integer function q_index_sq(x, low_arg, high_arg, mid_arg)
    real(real32), intent(in) :: x
    integer, intent(out), optional :: low_arg, high_arg, mid_arg
    integer :: low, high, mid

    if (x <= sq_qmap(0)) then
      q_index_sq = 0
      return
    end if
    if (sq_qmap(15) <= x) then
      q_index_sq = 15
      return
    end if
    low = 0
    high = 15
    do while (low < high)
      mid = ishft(low + high, -1)
      if (sq_qmap(mid) <= x) then
        low = mid + 1
      else
        high = mid
      end if
    end do
    if (present(low_arg)) low_arg = low
    if (present(high_arg)) high_arg = high
    if (present(mid_arg)) mid_arg = mid
    if (sq_qmidpt(low - 1) < x) then
      q_index_sq = low
    else
      q_index_sq = low - 1
    end if
  end function q_index_sq

end program main
