! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: rd = 8
  integer, parameter :: rd2 = rd * rd
  real(real64), parameter :: tolerance = 1.0e-3_real64
  real(real64), parameter :: c_rand_max = 2147483647.0_real64

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

  character(len=256) :: arg
  integer :: nseries, nobs, fc_steps, repeat
  integer :: rd2_word, rd_word, nobs_word, ns_word, fc_word
  real(real64), allocatable :: rqr(:), tmat(:), pmat(:), zvec(:), alpha(:), ys(:), mu(:)
  real(real64), allocatable :: vs(:), fs(:), sum_logfs(:), fc(:), f_fc(:), f_fc_ref(:)
  integer :: n_diff, iter
  real(real64) :: start_time, end_time
  logical :: ok

  if (command_argument_count() /= 4) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <#series> <#observations> <forcast steps> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nseries
  call get_command_argument(2, arg)
  read(arg, *) nobs
  call get_command_argument(3, arg)
  read(arg, *) fc_steps
  call get_command_argument(4, arg)
  read(arg, *) repeat

  if (nseries <= 0 .or. nobs <= 0 .or. fc_steps < 0 .or. repeat <= 0) stop 1

  rd2_word = nseries * rd2
  rd_word = nseries * rd
  nobs_word = nseries * nobs
  ns_word = nseries
  fc_word = fc_steps * nseries

  allocate(rqr(0:rd2_word - 1), tmat(0:rd2_word - 1), pmat(0:rd2_word - 1))
  allocate(zvec(0:rd_word - 1), alpha(0:rd_word - 1))
  allocate(ys(0:nobs_word - 1), mu(0:ns_word - 1))
  allocate(vs(0:nobs_word - 1), fs(0:nobs_word - 1), sum_logfs(0:ns_word - 1))
  allocate(fc(0:fc_word - 1), f_fc(0:fc_word - 1), f_fc_ref(0:fc_word - 1))

  call initialize_inputs(nseries, nobs, rqr, tmat, pmat, zvec, alpha, ys, mu)

  !$omp target data map(to: rqr(0:rd2_word - 1), tmat(0:rd2_word - 1), pmat(0:rd2_word - 1), &
  !$omp& zvec(0:rd_word - 1), alpha(0:rd_word - 1), ys(0:nobs_word - 1), mu(0:ns_word - 1)) &
  !$omp& map(alloc: vs(0:nobs_word - 1), fs(0:nobs_word - 1), sum_logfs(0:ns_word - 1), &
  !$omp& fc(0:fc_word - 1), f_fc(0:fc_word - 1))
  do n_diff = 0, rd - 1
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call kalman_device(ys, nobs, tmat, zvec, rqr, pmat, alpha, mu, nseries, vs, fs, &
        sum_logfs, n_diff, fc_steps, fc, f_fc)
    end do
    end_time = omp_get_wtime()
    write(*,'(A,I0,A,F0.6,A)') 'Average kernel execution time (n_diff = ', n_diff, '): ', &
      (end_time - start_time) / real(repeat, real64), ' (s)'
    if (fc_word > 0) then
      !$omp target update from(f_fc(0:fc_word - 1))
    end if
    call kalman_host(ys, nobs, tmat, zvec, rqr, pmat, alpha, mu, nseries, vs, fs, &
      sum_logfs, n_diff, fc_steps, fc, f_fc_ref)
    ok = compare_arrays(f_fc, f_fc_ref, fc_word)
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if
  end do
  !$omp end target data

  deallocate(f_fc_ref, f_fc, fc, sum_logfs, fs, vs, mu, ys, alpha, zvec, pmat, tmat, rqr)

contains

  subroutine initialize_inputs(nseries, nobs, rqr, tmat, pmat, zvec, alpha, ys, mu)
    integer, intent(in) :: nseries, nobs
    real(real64), intent(out) :: rqr(0:), tmat(0:), pmat(0:), zvec(0:), alpha(0:), ys(0:), mu(0:)
    integer :: i

    call c_srand(123_c_int)
    do i = 0, rd2 * nseries - 1
      rqr(i) = real(c_rand(), real64) / c_rand_max
    end do
    do i = 0, rd2 * nseries - 1
      tmat(i) = 1.0_real64
    end do
    do i = 0, rd2 * nseries - 1
      pmat(i) = real(c_rand(), real64) / c_rand_max
    end do
    do i = 0, rd * nseries - 1
      zvec(i) = real(c_rand(), real64) / c_rand_max
    end do
    do i = 0, rd * nseries - 1
      alpha(i) = real(c_rand(), real64) / c_rand_max
    end do
    do i = 0, nobs * nseries - 1
      ys(i) = real(c_rand(), real64) / c_rand_max
    end do
    do i = 0, nseries - 1
      mu(i) = real(c_rand(), real64) / c_rand_max
    end do
  end subroutine initialize_inputs

  subroutine kalman_device(ys, nobs, tmat, zvec, rqr, pmat, alpha, mu, batch_size, vs, fs, &
      sum_logfs, n_diff, fc_steps, fc, f_fc)
    real(real64), intent(in) :: ys(0:), tmat(0:), zvec(0:), rqr(0:), pmat(0:), alpha(0:), mu(0:)
    real(real64), intent(out) :: vs(0:), fs(0:), sum_logfs(0:), fc(0:), f_fc(0:)
    integer, intent(in) :: nobs, batch_size, n_diff, fc_steps
    integer :: bid, it, i, j, k, offset_rd, offset_rd2
    real(real64) :: l_rqr(rd2), l_t(rd2), l_z(rd), l_p(rd2), l_alpha(rd), l_k(rd)
    real(real64) :: l_tmp(rd2), l_tp(rd2), vs_it, fs_it, inv_fs, sum_value, pred, local_mu

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(bid, it, i, j, k, offset_rd, offset_rd2, l_rqr, l_t, l_z, l_p, l_alpha, &
    !$omp& l_k, l_tmp, l_tp, vs_it, fs_it, inv_fs, sum_value, pred, local_mu)
    do bid = 0, batch_size - 1
      offset_rd = bid * rd
      offset_rd2 = bid * rd2

      do i = 1, rd2
        l_rqr(i) = rqr(offset_rd2 + i - 1)
        l_t(i) = tmat(offset_rd2 + i - 1)
        l_p(i) = pmat(offset_rd2 + i - 1)
        l_tmp(i) = 0.0_real64
        l_tp(i) = 0.0_real64
      end do
      do i = 1, rd
        if (n_diff > 0) then
          l_z(i) = zvec(offset_rd + i - 1)
        else
          l_z(i) = 0.0_real64
        end if
        l_alpha(i) = alpha(offset_rd + i - 1)
        l_k(i) = 0.0_real64
      end do

      sum_value = 0.0_real64
      local_mu = mu(bid)

      do it = 0, nobs - 1
        vs_it = ys(bid * nobs + it)
        if (n_diff == 0) then
          vs_it = vs_it - l_alpha(1)
        else
          do i = 1, rd
            vs_it = vs_it - l_alpha(i) * l_z(i)
          end do
        end if
        vs(bid * nobs + it) = vs_it

        if (n_diff == 0) then
          fs_it = l_p(1)
        else
          fs_it = 0.0_real64
          do i = 1, rd
            do j = 1, rd
              fs_it = fs_it + l_p(j + (i - 1) * rd) * l_z(i) * l_z(j)
            end do
          end do
        end if
        fs(bid * nobs + it) = fs_it
        if (it >= n_diff) sum_value = sum_value + log(fs_it)

        do i = 1, rd
          do j = 1, rd
            l_tp(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_tp(i + (j - 1) * rd) = l_tp(i + (j - 1) * rd) + &
                l_t(i + (k - 1) * rd) * l_p(k + (j - 1) * rd)
            end do
          end do
        end do

        inv_fs = 1.0_real64 / fs_it
        if (n_diff == 0) then
          do i = 1, rd
            l_k(i) = inv_fs * l_tp(i)
          end do
        else
          do i = 1, rd
            l_k(i) = 0.0_real64
            do j = 1, rd
              l_k(i) = l_k(i) + l_tp(i + (j - 1) * rd) * l_z(j)
            end do
            l_k(i) = inv_fs * l_k(i)
          end do
        end if

        do i = 1, rd
          l_tmp(i) = 0.0_real64
          do j = 1, rd
            l_tmp(i) = l_tmp(i) + l_t(i + (j - 1) * rd) * l_alpha(j)
          end do
        end do
        do i = 1, rd
          l_alpha(i) = l_tmp(i) + l_k(i) * vs_it
        end do
        l_alpha(n_diff + 1) = l_alpha(n_diff + 1) + local_mu

        do i = 1, rd2
          l_tmp(i) = l_t(i)
        end do
        if (n_diff == 0) then
          do i = 1, rd
            l_tmp(i) = l_tmp(i) - l_k(i)
          end do
        else
          do i = 1, rd
            do j = 1, rd
              l_tmp(j + (i - 1) * rd) = l_tmp(j + (i - 1) * rd) - l_k(i) * l_z(j)
            end do
          end do
        end if

        do i = 1, rd
          do j = 1, rd
            l_p(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + &
                l_tp(i + (k - 1) * rd) * l_tmp(j + (k - 1) * rd)
            end do
            l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + l_rqr(i + (j - 1) * rd)
          end do
        end do
      end do
      sum_logfs(bid) = sum_value

      do it = 0, fc_steps - 1
        if (n_diff == 0) then
          pred = l_alpha(1)
        else
          pred = 0.0_real64
          do i = 1, rd
            pred = pred + l_alpha(i) * l_z(i)
          end do
        end if
        fc(bid * fc_steps + it) = pred

        do i = 1, rd
          l_tmp(i) = 0.0_real64
          do j = 1, rd
            l_tmp(i) = l_tmp(i) + l_t(i + (j - 1) * rd) * l_alpha(j)
          end do
        end do
        do i = 1, rd
          l_alpha(i) = l_tmp(i)
        end do
        l_alpha(n_diff + 1) = l_alpha(n_diff + 1) + local_mu

        if (n_diff == 0) then
          f_fc(bid * fc_steps + it) = l_p(1)
        else
          fs_it = 0.0_real64
          do i = 1, rd
            do j = 1, rd
              fs_it = fs_it + l_p(j + (i - 1) * rd) * l_z(i) * l_z(j)
            end do
          end do
          f_fc(bid * fc_steps + it) = fs_it
        end if

        do i = 1, rd
          do j = 1, rd
            l_tp(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_tp(i + (j - 1) * rd) = l_tp(i + (j - 1) * rd) + &
                l_t(i + (k - 1) * rd) * l_p(k + (j - 1) * rd)
            end do
          end do
        end do
        do i = 1, rd
          do j = 1, rd
            l_p(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + &
                l_tp(i + (k - 1) * rd) * l_t(j + (k - 1) * rd)
            end do
            l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + l_rqr(i + (j - 1) * rd)
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kalman_device

  subroutine kalman_host(ys, nobs, tmat, zvec, rqr, pmat, alpha, mu, batch_size, vs, fs, &
      sum_logfs, n_diff, fc_steps, fc, f_fc)
    real(real64), intent(in) :: ys(0:), tmat(0:), zvec(0:), rqr(0:), pmat(0:), alpha(0:), mu(0:)
    real(real64), intent(out) :: vs(0:), fs(0:), sum_logfs(0:), fc(0:), f_fc(0:)
    integer, intent(in) :: nobs, batch_size, n_diff, fc_steps
    integer :: bid, it, i, j, k, offset_rd, offset_rd2
    real(real64) :: l_rqr(rd2), l_t(rd2), l_z(rd), l_p(rd2), l_alpha(rd), l_k(rd)
    real(real64) :: l_tmp(rd2), l_tp(rd2), vs_it, fs_it, inv_fs, sum_value, pred, local_mu

    do bid = 0, batch_size - 1
      offset_rd = bid * rd
      offset_rd2 = bid * rd2

      do i = 1, rd2
        l_rqr(i) = rqr(offset_rd2 + i - 1)
        l_t(i) = tmat(offset_rd2 + i - 1)
        l_p(i) = pmat(offset_rd2 + i - 1)
        l_tmp(i) = 0.0_real64
        l_tp(i) = 0.0_real64
      end do
      do i = 1, rd
        if (n_diff > 0) then
          l_z(i) = zvec(offset_rd + i - 1)
        else
          l_z(i) = 0.0_real64
        end if
        l_alpha(i) = alpha(offset_rd + i - 1)
        l_k(i) = 0.0_real64
      end do

      sum_value = 0.0_real64
      local_mu = mu(bid)

      do it = 0, nobs - 1
        vs_it = ys(bid * nobs + it)
        if (n_diff == 0) then
          vs_it = vs_it - l_alpha(1)
        else
          do i = 1, rd
            vs_it = vs_it - l_alpha(i) * l_z(i)
          end do
        end if
        vs(bid * nobs + it) = vs_it

        if (n_diff == 0) then
          fs_it = l_p(1)
        else
          fs_it = 0.0_real64
          do i = 1, rd
            do j = 1, rd
              fs_it = fs_it + l_p(j + (i - 1) * rd) * l_z(i) * l_z(j)
            end do
          end do
        end if
        fs(bid * nobs + it) = fs_it
        if (it >= n_diff) sum_value = sum_value + log(fs_it)

        do i = 1, rd
          do j = 1, rd
            l_tp(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_tp(i + (j - 1) * rd) = l_tp(i + (j - 1) * rd) + &
                l_t(i + (k - 1) * rd) * l_p(k + (j - 1) * rd)
            end do
          end do
        end do

        inv_fs = 1.0_real64 / fs_it
        if (n_diff == 0) then
          do i = 1, rd
            l_k(i) = inv_fs * l_tp(i)
          end do
        else
          do i = 1, rd
            l_k(i) = 0.0_real64
            do j = 1, rd
              l_k(i) = l_k(i) + l_tp(i + (j - 1) * rd) * l_z(j)
            end do
            l_k(i) = inv_fs * l_k(i)
          end do
        end if

        do i = 1, rd
          l_tmp(i) = 0.0_real64
          do j = 1, rd
            l_tmp(i) = l_tmp(i) + l_t(i + (j - 1) * rd) * l_alpha(j)
          end do
        end do
        do i = 1, rd
          l_alpha(i) = l_tmp(i) + l_k(i) * vs_it
        end do
        l_alpha(n_diff + 1) = l_alpha(n_diff + 1) + local_mu

        do i = 1, rd2
          l_tmp(i) = l_t(i)
        end do
        if (n_diff == 0) then
          do i = 1, rd
            l_tmp(i) = l_tmp(i) - l_k(i)
          end do
        else
          do i = 1, rd
            do j = 1, rd
              l_tmp(j + (i - 1) * rd) = l_tmp(j + (i - 1) * rd) - l_k(i) * l_z(j)
            end do
          end do
        end if

        do i = 1, rd
          do j = 1, rd
            l_p(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + &
                l_tp(i + (k - 1) * rd) * l_tmp(j + (k - 1) * rd)
            end do
            l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + l_rqr(i + (j - 1) * rd)
          end do
        end do
      end do
      sum_logfs(bid) = sum_value

      do it = 0, fc_steps - 1
        if (n_diff == 0) then
          pred = l_alpha(1)
        else
          pred = 0.0_real64
          do i = 1, rd
            pred = pred + l_alpha(i) * l_z(i)
          end do
        end if
        fc(bid * fc_steps + it) = pred

        do i = 1, rd
          l_tmp(i) = 0.0_real64
          do j = 1, rd
            l_tmp(i) = l_tmp(i) + l_t(i + (j - 1) * rd) * l_alpha(j)
          end do
        end do
        do i = 1, rd
          l_alpha(i) = l_tmp(i)
        end do
        l_alpha(n_diff + 1) = l_alpha(n_diff + 1) + local_mu

        if (n_diff == 0) then
          f_fc(bid * fc_steps + it) = l_p(1)
        else
          fs_it = 0.0_real64
          do i = 1, rd
            do j = 1, rd
              fs_it = fs_it + l_p(j + (i - 1) * rd) * l_z(i) * l_z(j)
            end do
          end do
          f_fc(bid * fc_steps + it) = fs_it
        end if

        do i = 1, rd
          do j = 1, rd
            l_tp(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_tp(i + (j - 1) * rd) = l_tp(i + (j - 1) * rd) + &
                l_t(i + (k - 1) * rd) * l_p(k + (j - 1) * rd)
            end do
          end do
        end do
        do i = 1, rd
          do j = 1, rd
            l_p(i + (j - 1) * rd) = 0.0_real64
            do k = 1, rd
              l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + &
                l_tp(i + (k - 1) * rd) * l_t(j + (k - 1) * rd)
            end do
            l_p(i + (j - 1) * rd) = l_p(i + (j - 1) * rd) + l_rqr(i + (j - 1) * rd)
          end do
        end do
      end do
    end do
  end subroutine kalman_host

  logical function compare_arrays(lhs, rhs, n) result(ok)
    real(real64), intent(in) :: lhs(0:), rhs(0:)
    integer, intent(in) :: n
    integer :: i

    ok = .true.
    do i = 0, n - 1
      if (abs(lhs(i) - rhs(i)) > tolerance) then
        ok = .false.
        exit
      end if
    end do
  end function compare_arrays

end program main
