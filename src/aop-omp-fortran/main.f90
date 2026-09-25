! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_double, c_size_t
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: r_w_matrices_smem_slots = 12
  integer, parameter :: max_grid_size = 2048
  integer, parameter :: temp_storage_size = 4 * max_grid_size

  integer :: num_timesteps, num_paths_k, num_paths, num_runs
  real(real64) :: t_expiry, strike, s0, rate, sigma, dt
  logical :: price_put
  real(real64), allocatable :: samples(:), paths(:), svds(:), temp_storage(:)
  integer(int32), allocatable :: all_out_of_the_money(:)
  real(real64) :: h_price, cpu_price, ref_price, start_time, end_time, total_elapsed_ms
  integer :: run

  interface
    subroutine aop_reset_rng() bind(C, name='aop_reset_rng')
    end subroutine aop_reset_rng

    subroutine aop_fill_samples(samples, count) bind(C, name='aop_fill_samples')
      import :: c_double, c_size_t
      real(c_double), intent(out) :: samples(*)
      integer(c_size_t), value :: count
    end subroutine aop_fill_samples
  end interface

  num_timesteps = 100
  num_paths_k = 32
  num_runs = 1
  t_expiry = 1.00_real64
  strike = 4.00_real64
  s0 = 3.60_real64
  rate = 0.06_real64
  sigma = 0.20_real64
  price_put = .true.

  call parse_args(num_timesteps, num_paths_k, num_runs, t_expiry, s0, strike, rate, sigma, price_put)

  write(*,'(A)') '=============='
  write(*,'(A,I0)') 'Num Timesteps         : ', num_timesteps
  write(*,'(A,I0,A)') 'Num Paths             : ', num_paths_k, 'K'
  write(*,'(A,I0)') 'Num Runs              : ', num_runs
  write(*,'(A,F8.6)') 'T                     : ', t_expiry
  write(*,'(A,F8.6)') 'S0                    : ', s0
  write(*,'(A,F8.6)') 'K                     : ', strike
  write(*,'(A,F8.6)') 'r                     : ', rate
  write(*,'(A,F8.6)') 'sigma                 : ', sigma
  if (price_put) then
    write(*,'(A)') 'Option Type           : American Put'
  else
    write(*,'(A)') 'Option Type           : American Call'
  end if

  num_paths = num_paths_k * 1024
  dt = t_expiry / real(num_timesteps, real64)

  allocate(samples(0:num_timesteps * num_paths - 1))
  allocate(paths(0:num_timesteps * num_paths - 1))
  allocate(svds(0:16 * num_timesteps - 1))
  allocate(all_out_of_the_money(0:num_timesteps - 1))
  allocate(temp_storage(0:temp_storage_size - 1))

  h_price = 0.0_real64
  cpu_price = 0.0_real64
  total_elapsed_ms = 0.0_real64
  call aop_reset_rng()

  !$omp target data map(alloc: samples(0:num_timesteps*num_paths-1), paths(0:num_timesteps*num_paths-1), &
  !$omp& svds(0:16*num_timesteps-1), all_out_of_the_money(0:num_timesteps-1), temp_storage(0:temp_storage_size-1))
    do run = 1, num_runs
      call aop_fill_samples(samples, int(num_timesteps * num_paths, c_size_t))

      start_time = omp_get_wtime()
      call do_run(samples, num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, &
        paths, svds, all_out_of_the_money, temp_storage, h_price)
      end_time = omp_get_wtime()
      total_elapsed_ms = total_elapsed_ms + 1000.0_real64 * (end_time - start_time)

      if (run == 1) then
        call do_run_cpu(samples, num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, cpu_price)
        if (abs(cpu_price - h_price) > max(1.0e-7_real64, 1.0e-7_real64 * abs(cpu_price))) then
          write(*,'(A,1X,ES16.8,1X,ES16.8)') 'Fortran validation failed:', h_price, cpu_price
          stop 2
        end if
      end if
    end do
  !$omp end target data

  write(*,'(A)') '=============='
  write(*,'(A,F10.8)') 'GPU Longstaff-Schwartz: ', h_price

  if (price_put) then
    ref_price = binomial_tree(num_timesteps, .true., strike, dt, s0, rate, sigma)
  else
    ref_price = binomial_tree(num_timesteps, .false., strike, dt, s0, rate, sigma)
  end if
  write(*,'(A,F10.8)') 'Binonmial             : ', ref_price

  if (price_put) then
    ref_price = black_scholes_merton_put(t_expiry, strike, s0, rate, sigma)
  else
    ref_price = black_scholes_merton_call(t_expiry, strike, s0, rate, sigma)
  end if
  write(*,'(A,F10.8)') 'European Price        : ', ref_price
  write(*,'(A)') '=============='
  write(*,'(A,F0.3,A)') 'elapsed time for each run         : ', total_elapsed_ms / real(num_runs, real64), 'ms'
  write(*,'(A)') '=============='

  deallocate(samples, paths, svds, all_out_of_the_money, temp_storage)

contains

  subroutine parse_args(num_timesteps, num_paths, num_runs, t_expiry, s0, strike, rate, sigma, price_put)
    integer, intent(inout) :: num_timesteps, num_paths, num_runs
    real(real64), intent(inout) :: t_expiry, s0, strike, rate, sigma
    logical, intent(inout) :: price_put
    integer :: argc, idx
    character(len=256) :: arg, value

    argc = command_argument_count()
    idx = 1
    do while (idx <= argc)
      call get_command_argument(idx, arg)
      select case (trim(arg))
      case ('-timesteps')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) num_timesteps
        idx = idx + 2
      case ('-paths')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) num_paths
        idx = idx + 2
      case ('-runs')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) num_runs
        idx = idx + 2
      case ('-T')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) t_expiry
        idx = idx + 2
      case ('-S0')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) s0
        idx = idx + 2
      case ('-K')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) strike
        idx = idx + 2
      case ('-r')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) rate
        idx = idx + 2
      case ('-sigma')
        call require_value(idx, argc, arg)
        call get_command_argument(idx + 1, value)
        read(value, *) sigma
        idx = idx + 2
      case ('-call')
        price_put = .false.
        idx = idx + 1
      case default
        write(*,'(A,A,A)') 'Unknown option ', trim(arg), '. Aborting!!!'
        stop 1
      end select
    end do

    if (num_timesteps <= 1 .or. num_paths <= 0 .or. num_runs <= 0) stop 1
  end subroutine parse_args

  subroutine require_value(idx, argc, arg)
    integer, intent(in) :: idx, argc
    character(len=*), intent(in) :: arg
    if (idx == argc) then
      write(*,'(A,A)') 'Missing value for ', trim(arg)
      stop 1
    end if
  end subroutine require_value

  real(real64) function payoff_value(is_put, strike, s) result(value)
    logical, intent(in) :: is_put
    real(real64), intent(in) :: strike, s
    if (is_put) then
      value = max(strike - s, 0.0_real64)
    else
      value = max(s - strike, 0.0_real64)
    end if
  end function payoff_value

  integer function is_in_the_money(is_put, strike, s) result(value)
    logical, intent(in) :: is_put
    real(real64), intent(in) :: strike, s
    if (is_put) then
      value = merge(1, 0, s < strike)
    else
      value = merge(1, 0, s > strike)
    end if
  end function is_in_the_money

  subroutine do_run(samples, num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, &
      paths, svds, all_out_of_the_money, temp_storage, h_price)
    real(real64), intent(inout) :: samples(0:)
    integer, intent(in) :: num_timesteps, num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, dt, s0, rate, sigma
    real(real64), intent(inout) :: paths(0:), svds(0:), temp_storage(0:)
    integer(int32), intent(inout) :: all_out_of_the_money(0:)
    real(real64), intent(out) :: h_price
    integer :: timestep
    real(real64) :: exp_min_r_dt

    !$omp target update to(samples(0:num_timesteps*num_paths-1))

    call generate_paths_kernel(num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, samples, paths)
    call reset_flags_kernel(num_timesteps, all_out_of_the_money)
    call prepare_svd_kernel(num_timesteps - 1, num_paths, 4, price_put, strike, paths, all_out_of_the_money, svds)

    exp_min_r_dt = exp(-rate * dt)
    do timestep = num_timesteps - 2, 0, -1
      call compute_beta_kernel(num_paths, price_put, strike, svds(16 * timestep:), paths(timestep * num_paths:), &
        paths((num_timesteps - 1) * num_paths:), all_out_of_the_money(timestep:), temp_storage)
      call update_cashflow_kernel(num_paths, price_put, strike, exp_min_r_dt, temp_storage, &
        paths(timestep * num_paths:), all_out_of_the_money(timestep:), paths((num_timesteps - 1) * num_paths:))
    end do

    call compute_sums_kernel(num_paths, paths((num_timesteps - 1) * num_paths:), exp_min_r_dt, h_price)
  end subroutine do_run

  subroutine generate_paths_kernel(num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, samples, paths)
    integer, intent(in) :: num_timesteps, num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, dt, s0, rate, sigma
    real(real64), intent(in) :: samples(0:)
    real(real64), intent(out) :: paths(0:)
    integer :: path, timestep, offset
    real(real64) :: r_min_half_sigma_sq_dt, sigma_sqrt_dt, s

    r_min_half_sigma_sq_dt = (rate - 0.5_real64 * sigma * sigma) * dt
    sigma_sqrt_dt = sigma * sqrt(dt)
    !$omp target teams distribute parallel do thread_limit(256) private(timestep, offset, s)
    do path = 0, num_paths - 1
      s = s0
      offset = path
      do timestep = 0, num_timesteps - 2
        s = s * exp(r_min_half_sigma_sq_dt + sigma_sqrt_dt * samples(offset))
        paths(offset) = s
        offset = offset + num_paths
      end do
      s = s * exp(r_min_half_sigma_sq_dt + sigma_sqrt_dt * samples(offset))
      paths(offset) = payoff_value(price_put, strike, s)
    end do
    !$omp end target teams distribute parallel do
  end subroutine generate_paths_kernel

  subroutine reset_flags_kernel(num_timesteps, all_out_of_the_money)
    integer, intent(in) :: num_timesteps
    integer(int32), intent(out) :: all_out_of_the_money(0:)
    integer :: i
    !$omp target teams distribute parallel do thread_limit(256)
    do i = 0, num_timesteps - 1
      all_out_of_the_money(i) = 0_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine reset_flags_kernel

  subroutine prepare_svd_kernel(num_teams, num_paths, min_in_the_money, price_put, strike, paths, all_out_of_the_money, svds)
    integer, intent(in) :: num_teams, num_paths, min_in_the_money
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike
    real(real64), intent(in) :: paths(0:)
    integer(int32), intent(inout) :: all_out_of_the_money(0:)
    real(real64), intent(out) :: svds(0:)
    integer, parameter :: num_threads_per_block = 256
    integer :: scan_input(0:num_threads_per_block - 1), scan_output(0:num_threads_per_block)
    integer :: timestep, path, offset, m, found_paths, slot, lid, bid, in_the_money
    integer :: partial_sum, total_sum, lsum, not_enough_paths
    real(real64) :: s, x, x_sq, sums(4), lsums(4), smem_svds(r_w_matrices_smem_slots)

    !$omp target teams num_teams(num_teams) thread_limit(num_threads_per_block) &
    !$omp& private(scan_input, scan_output, smem_svds, lsums, lsum)
    !$omp parallel private(lid, bid, timestep, offset, m, sums, found_paths, path, s, in_the_money, partial_sum, &
    !$omp& total_sum, x, x_sq, not_enough_paths, slot)
      lid = omp_get_thread_num()
      bid = omp_get_team_num()
      timestep = bid
      offset = timestep * num_paths
      sums = 0.0_real64
      m = 0
      found_paths = 0
      if (lid < r_w_matrices_smem_slots) smem_svds(lid + 1) = 0.0_real64
      !$omp barrier

      do path = lid, num_paths - 1, num_threads_per_block
        s = paths(offset + path)
        in_the_money = is_in_the_money(price_put, strike, s)

        scan_input(lid) = in_the_money
        !$omp barrier
        if (lid == 0) then
          scan_output(0) = 0
          do slot = 1, num_threads_per_block
            scan_output(slot) = scan_output(slot - 1) + scan_input(slot - 1)
          end do
        end if
        !$omp barrier
        partial_sum = scan_output(lid)
        total_sum = scan_output(num_threads_per_block)

        if (found_paths < 3) then
          if (in_the_money /= 0 .and. found_paths + partial_sum < 3) smem_svds(found_paths + partial_sum + 1) = s
          !$omp barrier
          found_paths = found_paths + total_sum
        end if

        if (lid == 0) lsum = 0
        !$omp barrier
        !$omp atomic update
        lsum = ior(lsum, in_the_money)
        !$omp barrier
        if (lsum == 0) cycle

        m = m + in_the_money
        x = 0.0_real64
        x_sq = 0.0_real64
        if (in_the_money /= 0) then
          x = s
          x_sq = s * s
        end if
        sums(1) = sums(1) + x
        sums(2) = sums(2) + x_sq
        sums(3) = sums(3) + x_sq * x
        sums(4) = sums(4) + x_sq * x_sq
      end do
      !$omp barrier

      if (lid == 0) lsum = 0
      !$omp barrier
      !$omp atomic update
      lsum = lsum + m
      !$omp barrier

      not_enough_paths = 0
      if (lid == 0 .and. lsum < min_in_the_money) not_enough_paths = 1
      !$omp barrier

      if (not_enough_paths /= 0) then
        if (lid == 0) all_out_of_the_money(bid) = 1_int32
      else
        if (lid == 0) lsums = 0.0_real64
        !$omp barrier
        !$omp atomic update
        lsums(1) = lsums(1) + sums(1)
        !$omp atomic update
        lsums(2) = lsums(2) + sums(2)
        !$omp atomic update
        lsums(3) = lsums(3) + sums(3)
        !$omp barrier

        if (lid == 0) call svd_3x3(lsum, lsums, smem_svds)
        !$omp barrier
        if (lid < r_w_matrices_smem_slots) svds(16 * bid + lid) = smem_svds(lid + 1)
      end if
    !$omp end parallel
    !$omp end target teams
  end subroutine prepare_svd_kernel

  subroutine compute_beta_kernel(num_paths, price_put, strike, svd, paths, cashflows, all_out_of_the_money, beta)
    integer, intent(in) :: num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike
    real(real64), intent(in) :: svd(0:), paths(0:), cashflows(0:)
    integer(int32), intent(in) :: all_out_of_the_money(0:)
    real(real64), intent(inout) :: beta(0:)
    integer :: path, in_money
    real(real64) :: b0, b1, b2, s, q1i, q2i, cashflow
    real(real64) :: r00, r01, r02, r11, r12, r22, w00, w01, w02, w11, w12, w22
    real(real64) :: inv_r00, inv_r11, inv_r22, inv_r01, inv_r02, inv_r12, inv_w00, wi0, wi1, wi2

    b0 = 0.0_real64
    b1 = 0.0_real64
    b2 = 0.0_real64
    !$omp target teams distribute parallel do thread_limit(128) reduction(+:b0,b1,b2) &
    !$omp& private(s, q1i, q2i, cashflow, in_money, r00, r01, r02, r11, r12, r22, w00, w01, w02, w11, w12, w22, &
    !$omp& inv_r00, inv_r11, inv_r22, inv_r01, inv_r02, inv_r12, inv_w00, wi0, wi1, wi2)
    do path = 0, num_paths - 1
      if (all_out_of_the_money(0) == 0_int32) then
        r00 = svd(0); r01 = svd(1); r02 = svd(2)
        r11 = svd(3); r12 = svd(4); r22 = svd(5)
        w00 = svd(6); w01 = svd(7); w02 = svd(8)
        w11 = svd(9); w12 = svd(10); w22 = svd(11)
        inv_r00 = merge(1.0_real64 / r00, 0.0_real64, r00 /= 0.0_real64)
        inv_r11 = merge(1.0_real64 / r11, 0.0_real64, r11 /= 0.0_real64)
        inv_r22 = merge(1.0_real64 / r22, 0.0_real64, r22 /= 0.0_real64)
        inv_r01 = inv_r00 * inv_r11 * r01
        inv_r02 = inv_r00 * inv_r22 * r02
        inv_r12 = inv_r22 * r12
        inv_w00 = w00 * inv_r00
        s = paths(path)
        in_money = is_in_the_money(price_put, strike, s)
        q1i = inv_r11 * s - inv_r01
        q2i = inv_r22 * s * s - inv_r02 - q1i * inv_r12
        wi0 = inv_w00 + w01 * q1i + w02 * q2i
        wi1 =           w11 * q1i + w12 * q2i
        wi2 =                       w22 * q2i
        if (in_money /= 0) then
          cashflow = cashflows(path)
        else
          cashflow = 0.0_real64
        end if
        b0 = b0 + wi0 * cashflow
        b1 = b1 + wi1 * cashflow
        b2 = b2 + wi2 * cashflow
      end if
    end do
    !$omp end target teams distribute parallel do
    beta(0) = b0
    beta(1) = b1
    beta(2) = b2
  end subroutine compute_beta_kernel

  subroutine update_cashflow_kernel(num_paths, price_put, strike, exp_min_r_dt, beta, paths, all_out_of_the_money, cashflows)
    integer, intent(in) :: num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, exp_min_r_dt
    real(real64), intent(in) :: beta(0:), paths(0:)
    integer(int32), intent(in) :: all_out_of_the_money(0:)
    real(real64), intent(inout) :: cashflows(0:)
    integer :: path
    real(real64) :: beta0, beta1, beta2, old_cashflow, s, s2, payoff, estimated_payoff

    beta0 = beta(0)
    beta1 = beta(1)
    beta2 = beta(2)
    !$omp target teams distribute parallel do thread_limit(128) private(old_cashflow, s, s2, payoff, estimated_payoff)
    do path = 0, num_paths - 1
      old_cashflow = exp_min_r_dt * cashflows(path)
      if (all_out_of_the_money(0) /= 0_int32) then
        cashflows(path) = old_cashflow
      else
        s = paths(path)
        s2 = s * s
        payoff = payoff_value(price_put, strike, s)
        estimated_payoff = (beta0 + beta1 * s + beta2 * s2) * exp_min_r_dt
        if (payoff <= 1.0e-8_real64 .or. payoff <= estimated_payoff) payoff = old_cashflow
        cashflows(path) = payoff
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_cashflow_kernel

  subroutine compute_sums_kernel(num_paths, cashflows, exp_min_r_dt, price)
    integer, intent(in) :: num_paths
    real(real64), intent(in) :: cashflows(0:), exp_min_r_dt
    real(real64), intent(out) :: price
    integer :: path
    real(real64) :: sum

    sum = 0.0_real64
    !$omp target teams distribute parallel do thread_limit(128) reduction(+:sum)
    do path = 0, num_paths - 1
      sum = sum + cashflows(path)
    end do
    !$omp end target teams distribute parallel do
    price = exp_min_r_dt * sum / real(num_paths, real64)
  end subroutine compute_sums_kernel

  subroutine do_run_cpu(samples, num_timesteps, num_paths, price_put, strike, dt, s0, rate, sigma, h_price)
    real(real64), intent(in) :: samples(0:)
    integer, intent(in) :: num_timesteps, num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, dt, s0, rate, sigma
    real(real64), intent(out) :: h_price
    real(real64), allocatable :: c_paths(:), c_svds(:), c_beta(:)
    integer(int32), allocatable :: c_all_out(:)
    integer :: path, timestep, offset, m, found_paths, slot
    real(real64) :: drift, vol, s, sums(4), smem_svds(r_w_matrices_smem_slots), exp_min_r_dt

    allocate(c_paths(0:num_timesteps * num_paths - 1), c_svds(0:16 * num_timesteps - 1), c_beta(0:2))
    allocate(c_all_out(0:num_timesteps - 1))
    drift = (rate - 0.5_real64 * sigma * sigma) * dt
    vol = sigma * sqrt(dt)
    do path = 0, num_paths - 1
      s = s0
      offset = path
      do timestep = 0, num_timesteps - 2
        s = s * exp(drift + vol * samples(offset))
        c_paths(offset) = s
        offset = offset + num_paths
      end do
      s = s * exp(drift + vol * samples(offset))
      c_paths(offset) = payoff_value(price_put, strike, s)
    end do

    c_all_out = 0_int32
    do timestep = 0, num_timesteps - 2
      sums = 0.0_real64
      smem_svds = 0.0_real64
      m = 0
      found_paths = 0
      offset = timestep * num_paths
      do path = 0, num_paths - 1
        s = c_paths(offset + path)
        if (is_in_the_money(price_put, strike, s) /= 0) then
          if (found_paths < 3) then
            smem_svds(found_paths + 1) = s
            found_paths = found_paths + 1
          end if
          m = m + 1
          sums(1) = sums(1) + s
          sums(2) = sums(2) + s * s
          sums(3) = sums(3) + s * s * s
        end if
      end do
      if (m < 4) then
        c_all_out(timestep) = 1_int32
      else
        call svd_3x3(m, sums, smem_svds)
        do slot = 0, r_w_matrices_smem_slots - 1
          c_svds(16 * timestep + slot) = smem_svds(slot + 1)
        end do
      end if
    end do

    exp_min_r_dt = exp(-rate * dt)
    do timestep = num_timesteps - 2, 0, -1
      call compute_beta_cpu(num_paths, price_put, strike, c_svds(16 * timestep:), c_paths(timestep * num_paths:), &
        c_paths((num_timesteps - 1) * num_paths:), c_all_out(timestep), c_beta)
      call update_cashflow_cpu(num_paths, price_put, strike, exp_min_r_dt, c_beta, c_paths(timestep * num_paths:), &
        c_all_out(timestep), c_paths((num_timesteps - 1) * num_paths:))
    end do
    h_price = exp_min_r_dt * sum(c_paths((num_timesteps - 1) * num_paths:(num_timesteps * num_paths - 1))) / real(num_paths, real64)
    deallocate(c_paths, c_svds, c_beta, c_all_out)
  end subroutine do_run_cpu

  subroutine compute_beta_cpu(num_paths, price_put, strike, svd, paths, cashflows, all_out, beta)
    integer, intent(in) :: num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, svd(0:), paths(0:), cashflows(0:)
    integer(int32), intent(in) :: all_out
    real(real64), intent(out) :: beta(0:)
    integer :: path, in_money
    real(real64) :: s, q1i, q2i, cashflow
    real(real64) :: r00, r01, r02, r11, r12, r22, w00, w01, w02, w11, w12, w22
    real(real64) :: inv_r00, inv_r11, inv_r22, inv_r01, inv_r02, inv_r12, inv_w00, wi0, wi1, wi2

    beta = 0.0_real64
    if (all_out /= 0_int32) return
    r00 = svd(0); r01 = svd(1); r02 = svd(2)
    r11 = svd(3); r12 = svd(4); r22 = svd(5)
    w00 = svd(6); w01 = svd(7); w02 = svd(8)
    w11 = svd(9); w12 = svd(10); w22 = svd(11)
    inv_r00 = merge(1.0_real64 / r00, 0.0_real64, r00 /= 0.0_real64)
    inv_r11 = merge(1.0_real64 / r11, 0.0_real64, r11 /= 0.0_real64)
    inv_r22 = merge(1.0_real64 / r22, 0.0_real64, r22 /= 0.0_real64)
    inv_r01 = inv_r00 * inv_r11 * r01
    inv_r02 = inv_r00 * inv_r22 * r02
    inv_r12 = inv_r22 * r12
    inv_w00 = w00 * inv_r00
    do path = 0, num_paths - 1
      s = paths(path)
      in_money = is_in_the_money(price_put, strike, s)
      q1i = inv_r11 * s - inv_r01
      q2i = inv_r22 * s * s - inv_r02 - q1i * inv_r12
      wi0 = inv_w00 + w01 * q1i + w02 * q2i
      wi1 =           w11 * q1i + w12 * q2i
      wi2 =                       w22 * q2i
      if (in_money /= 0) then
        cashflow = cashflows(path)
      else
        cashflow = 0.0_real64
      end if
      beta(0) = beta(0) + wi0 * cashflow
      beta(1) = beta(1) + wi1 * cashflow
      beta(2) = beta(2) + wi2 * cashflow
    end do
  end subroutine compute_beta_cpu

  subroutine update_cashflow_cpu(num_paths, price_put, strike, exp_min_r_dt, beta, paths, all_out, cashflows)
    integer, intent(in) :: num_paths
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, exp_min_r_dt, beta(0:), paths(0:)
    integer(int32), intent(in) :: all_out
    real(real64), intent(inout) :: cashflows(0:)
    integer :: path
    real(real64) :: old_cashflow, s, payoff, estimated_payoff

    do path = 0, num_paths - 1
      old_cashflow = exp_min_r_dt * cashflows(path)
      if (all_out /= 0_int32) then
        cashflows(path) = old_cashflow
      else
        s = paths(path)
        payoff = payoff_value(price_put, strike, s)
        estimated_payoff = (beta(0) + beta(1) * s + beta(2) * s * s) * exp_min_r_dt
        if (payoff <= 1.0e-8_real64 .or. payoff <= estimated_payoff) payoff = old_cashflow
        cashflows(path) = payoff
      end if
    end do
  end subroutine update_cashflow_cpu

  subroutine assemble_r(m, sums, smem_svds)
    integer, intent(in) :: m
    real(real64), intent(in) :: sums(4)
    real(real64), intent(inout) :: smem_svds(r_w_matrices_smem_slots)
    real(real64) :: x0, x1, x2, x0_sq, sum1, sum2, sum3, sum4, m_as_dbl
    real(real64) :: sigma_l, mu, v0, v0_sq, beta, inv_v0, one_min_beta, beta_div_v0
    real(real64) :: beta_div_v0_sq, c1, c2, x1_sq, c3, c4, c5, x2_sq

    x0 = smem_svds(1)
    x1 = smem_svds(2)
    x2 = smem_svds(3)
    x0_sq = x0 * x0
    sum1 = sums(1) - x0
    sum2 = sums(2) - x0_sq
    sum3 = sums(3) - x0_sq * x0
    sum4 = sums(4) - x0_sq * x0_sq
    m_as_dbl = real(m, real64)
    sigma_l = m_as_dbl - 1.0_real64
    mu = sqrt(m_as_dbl)
    v0 = -sigma_l / (1.0_real64 + mu)
    v0_sq = v0 * v0
    beta = 2.0_real64 * v0_sq / (sigma_l + v0_sq)
    inv_v0 = 1.0_real64 / v0
    one_min_beta = 1.0_real64 - beta
    beta_div_v0 = beta * inv_v0
    smem_svds(1) = mu
    smem_svds(2) = one_min_beta * x0 - beta_div_v0 * sum1
    smem_svds(3) = one_min_beta * x0_sq - beta_div_v0 * sum2
    beta_div_v0_sq = beta_div_v0 * inv_v0
    c1 = beta_div_v0_sq * sum1 + beta_div_v0 * x0
    c2 = beta_div_v0_sq * sum2 + beta_div_v0 * x0_sq

    x1_sq = x1 * x1
    sum1 = sum1 - x1
    sum2 = sum2 - x1_sq
    sum3 = sum3 - x1_sq * x1
    sum4 = sum4 - x1_sq * x1_sq
    x0 = x1 - c1
    x0_sq = x0 * x0
    sigma_l = sum2 - 2.0_real64 * c1 * sum1 + (m_as_dbl - 2.0_real64) * c1 * c1
    if (abs(sigma_l) < 1.0e-16_real64) then
      beta = 0.0_real64
      v0 = 1.0_real64
    else
      mu = sqrt(x0_sq + sigma_l)
      if (x0 <= 0.0_real64) then
        v0 = x0 - mu
      else
        v0 = -sigma_l / (x0 + mu)
      end if
      v0_sq = v0 * v0
      beta = 2.0_real64 * v0_sq / (sigma_l + v0_sq)
    end if
    inv_v0 = 1.0_real64 / v0
    beta_div_v0 = beta * inv_v0
    c3 = (sum3 - c1 * sum2 - c2 * sum1 + (m_as_dbl - 2.0_real64) * c1 * c2) * beta_div_v0
    c4 = (x1_sq - c2) * beta_div_v0 + c3 * inv_v0
    c5 = c1 * c4 - c2
    one_min_beta = 1.0_real64 - beta
    smem_svds(4) = one_min_beta * x0 - beta_div_v0 * sigma_l
    smem_svds(5) = one_min_beta * (x1_sq - c2) - c3

    x2_sq = x2 * x2
    sum1 = sum1 - x2
    sum2 = sum2 - x2_sq
    sum3 = sum3 - x2_sq * x2
    sum4 = sum4 - x2_sq * x2_sq
    x0 = x2_sq - c4 * x2 + c5
    sigma_l = sum4 - 2.0_real64 * c4 * sum3 + (c4 * c4 + 2.0_real64 * c5) * sum2 - &
      2.0_real64 * c4 * c5 * sum1 + (m_as_dbl - 3.0_real64) * c5 * c5
    if (abs(sigma_l) < 1.0e-12_real64) then
      beta = 0.0_real64
      v0 = 1.0_real64
    else
      mu = sqrt(x0 * x0 + sigma_l)
      if (x0 <= 0.0_real64) then
        v0 = x0 - mu
      else
        v0 = -sigma_l / (x0 + mu)
      end if
      v0_sq = v0 * v0
      beta = 2.0_real64 * v0_sq / (sigma_l + v0_sq)
    end if
    smem_svds(6) = (1.0_real64 - beta) * x0 - (beta / v0) * sigma_l
  end subroutine assemble_r

  real(real64) function off_diag_norm(a01, a02, a12) result(value)
    real(real64), intent(in) :: a01, a02, a12
    value = sqrt(2.0_real64 * (a01 * a01 + a02 * a02 + a12 * a12))
  end function off_diag_norm

  subroutine swap_real(x, y)
    real(real64), intent(inout) :: x, y
    real(real64) :: t
    t = x
    x = y
    y = t
  end subroutine swap_real

  subroutine svd_3x3(m, sums, smem_svds)
    integer, intent(in) :: m
    real(real64), intent(in) :: sums(4)
    real(real64), intent(inout) :: smem_svds(r_w_matrices_smem_slots)
    integer :: iter
    real(real64) :: r00, r01, r02, r11, r12, r22
    real(real64) :: a00, a01, a02, a11, a12, a22
    real(real64) :: v00, v01, v02, v10, v11, v12, v20, v21, v22
    real(real64) :: c, s, b00, b01, b02, b10, b11, b12, b20, b21, b22
    real(real64) :: tau, sgn, t, inv_s0, inv_s1, inv_s2
    real(real64) :: u00, u01, u02, u10, u11, u12, u20, u21, u22
    real(real64) :: bb00, bb01, bb02, bb11, bb12, bb22

    call assemble_r(m, sums, smem_svds)
    r00 = smem_svds(1); r01 = smem_svds(2); r02 = smem_svds(3)
    r11 = smem_svds(4); r12 = smem_svds(5); r22 = smem_svds(6)
    a00 = r00 * r00
    a01 = r00 * r01
    a02 = r00 * r02
    a11 = r01 * r01 + r11 * r11
    a12 = r01 * r02 + r11 * r12
    a22 = r02 * r02 + r12 * r12 + r22 * r22
    v00 = 1.0_real64; v01 = 0.0_real64; v02 = 0.0_real64
    v10 = 0.0_real64; v11 = 1.0_real64; v12 = 0.0_real64
    v20 = 0.0_real64; v21 = 0.0_real64; v22 = 1.0_real64

    iter = 0
    do while (off_diag_norm(a01, a02, a12) >= 1.0e-12_real64 .and. iter < 16)
      iter = iter + 1
      c = 1.0_real64; s = 0.0_real64
      if (a01 /= 0.0_real64) then
        tau = (a11 - a00) / (2.0_real64 * a01)
        sgn = merge(-1.0_real64, 1.0_real64, tau < 0.0_real64)
        t = sgn / (sgn * tau + sqrt(1.0_real64 + tau * tau))
        c = 1.0_real64 / sqrt(1.0_real64 + t * t)
        s = t * c
      end if
      b00 = c * a00 - s * a01
      b01 = s * a00 + c * a01
      b10 = c * a01 - s * a11
      b11 = s * a01 + c * a11
      b02 = a02
      a00 = c * b00 - s * b10
      a01 = c * b01 - s * b11
      a11 = s * b01 + c * b11
      a02 = c * b02 - s * a12
      a12 = s * b02 + c * a12
      b00 = c * v00 - s * v01; v01 = s * v00 + c * v01; v00 = b00
      b10 = c * v10 - s * v11; v11 = s * v10 + c * v11; v10 = b10
      b20 = c * v20 - s * v21; v21 = s * v20 + c * v21; v20 = b20

      c = 1.0_real64; s = 0.0_real64
      if (a02 /= 0.0_real64) then
        tau = (a22 - a00) / (2.0_real64 * a02)
        sgn = merge(-1.0_real64, 1.0_real64, tau < 0.0_real64)
        t = sgn / (sgn * tau + sqrt(1.0_real64 + tau * tau))
        c = 1.0_real64 / sqrt(1.0_real64 + t * t)
        s = t * c
      end if
      b00 = c * a00 - s * a02
      b01 = c * a01 - s * a12
      b02 = s * a00 + c * a02
      b20 = c * a02 - s * a22
      b22 = s * a02 + c * a22
      a00 = c * b00 - s * b20
      a12 = s * a01 + c * a12
      a02 = c * b02 - s * b22
      a22 = s * b02 + c * b22
      a01 = b01
      b00 = c * v00 - s * v02; v02 = s * v00 + c * v02; v00 = b00
      b10 = c * v10 - s * v12; v12 = s * v10 + c * v12; v10 = b10
      b20 = c * v20 - s * v22; v22 = s * v20 + c * v22; v20 = b20

      c = 1.0_real64; s = 0.0_real64
      if (a12 /= 0.0_real64) then
        tau = (a22 - a11) / (2.0_real64 * a12)
        sgn = merge(-1.0_real64, 1.0_real64, tau < 0.0_real64)
        t = sgn / (sgn * tau + sqrt(1.0_real64 + tau * tau))
        c = 1.0_real64 / sqrt(1.0_real64 + t * t)
        s = t * c
      end if
      b02 = s * a01 + c * a02
      b11 = c * a11 - s * a12
      b12 = s * a11 + c * a12
      b21 = c * a12 - s * a22
      b22 = s * a12 + c * a22
      a01 = c * a01 - s * a02
      a02 = b02
      a11 = c * b11 - s * b21
      a12 = c * b12 - s * b22
      a22 = s * b12 + c * b22
      b01 = c * v01 - s * v02; v02 = s * v01 + c * v02; v01 = b01
      b11 = c * v11 - s * v12; v12 = s * v11 + c * v12; v11 = b11
      b21 = c * v21 - s * v22; v22 = s * v21 + c * v22; v21 = b21
    end do

    if (a00 < a11) then
      call swap_real(a00, a11); call swap_real(v00, v01); call swap_real(v10, v11); call swap_real(v20, v21)
    end if
    if (a00 < a22) then
      call swap_real(a00, a22); call swap_real(v00, v02); call swap_real(v10, v12); call swap_real(v20, v22)
    end if
    if (a11 < a22) then
      call swap_real(a11, a22); call swap_real(v01, v02); call swap_real(v11, v12); call swap_real(v21, v22)
    end if

    inv_s0 = merge(1.0_real64 / a00, 0.0_real64, abs(a00) >= 1.0e-12_real64)
    inv_s1 = merge(1.0_real64 / a11, 0.0_real64, abs(a11) >= 1.0e-12_real64)
    inv_s2 = merge(1.0_real64 / a22, 0.0_real64, abs(a22) >= 1.0e-12_real64)
    u00 = v00 * inv_s0; u01 = v01 * inv_s1; u02 = v02 * inv_s2
    u10 = v10 * inv_s0; u11 = v11 * inv_s1; u12 = v12 * inv_s2
    u20 = v20 * inv_s0; u21 = v21 * inv_s1; u22 = v22 * inv_s2
    bb00 = u00 * v00 + u01 * v01 + u02 * v02
    bb01 = u00 * v10 + u01 * v11 + u02 * v12
    bb02 = u00 * v20 + u01 * v21 + u02 * v22
    bb11 = u10 * v10 + u11 * v11 + u12 * v12
    bb12 = u10 * v20 + u11 * v21 + u12 * v22
    bb22 = u20 * v20 + u21 * v21 + u22 * v22
    smem_svds(7) = bb00 * r00 + bb01 * r01 + bb02 * r02
    smem_svds(8) =              bb01 * r11 + bb02 * r12
    smem_svds(9) =                           bb02 * r22
    smem_svds(10) =             bb11 * r11 + bb12 * r12
    smem_svds(11) =                          bb12 * r22
    smem_svds(12) =                          bb22 * r22
  end subroutine svd_3x3

  real(real64) function binomial_tree(num_timesteps, price_put, strike, dt, s0, rate, sigma) result(value)
    integer, intent(in) :: num_timesteps
    logical, intent(in) :: price_put
    real(real64), intent(in) :: strike, dt, s0, rate, sigma
    real(real64), allocatable :: tree(:)
    real(real64) :: u, d, a, p, k, expected, earlyex
    integer :: t, i

    allocate(tree(0:num_timesteps))
    u = exp(sigma * sqrt(dt))
    d = exp(-sigma * sqrt(dt))
    a = exp(rate * dt)
    p = (a - d) / (u - d)
    k = d ** num_timesteps
    do t = 0, num_timesteps
      tree(t) = payoff_value(price_put, strike, s0 * k)
      k = k * u * u
    end do
    do t = num_timesteps - 1, 0, -1
      k = d ** t
      do i = 0, t
        expected = exp(-rate * dt) * (p * tree(i + 1) + (1.0_real64 - p) * tree(i))
        earlyex = payoff_value(price_put, strike, s0 * k)
        tree(i) = max(earlyex, expected)
        k = k * u * u
      end do
    end do
    value = tree(0)
    deallocate(tree)
  end function binomial_tree

  real(real64) function normcdf(x) result(value)
    real(real64), intent(in) :: x
    value = (1.0_real64 + erf(x / sqrt(2.0_real64))) / 2.0_real64
  end function normcdf

  real(real64) function black_scholes_merton_put(t_expiry, strike, s0, rate, sigma) result(value)
    real(real64), intent(in) :: t_expiry, strike, s0, rate, sigma
    real(real64) :: d1, d2
    d1 = (log(s0 / strike) + (rate + 0.5_real64 * sigma * sigma) * t_expiry) / (sigma * sqrt(t_expiry))
    d2 = d1 - sigma * sqrt(t_expiry)
    value = strike * exp(-rate * t_expiry) * normcdf(-d2) - s0 * normcdf(-d1)
  end function black_scholes_merton_put

  real(real64) function black_scholes_merton_call(t_expiry, strike, s0, rate, sigma) result(value)
    real(real64), intent(in) :: t_expiry, strike, s0, rate, sigma
    real(real64) :: d1, d2
    d1 = (log(s0 / strike) + (rate + 0.5_real64 * sigma * sigma) * t_expiry) / (sigma * sqrt(t_expiry))
    d2 = d1 - sigma * sqrt(t_expiry)
    value = s0 * normcdf(d1) - strike * exp(-rate * t_expiry) * normcdf(d2)
  end function black_scholes_merton_call

end program main
