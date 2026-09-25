! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer, parameter :: max_wg_size = 256
  real(real32), parameter :: score_max = huge(1.0_real32)
  integer :: num_pac, num_snp, iteration, block_snp
  integer :: pp_zeros, pp_ones, phen_ones, mask_zeros, mask_ones
  integer(int32), allocatable :: snp_data(:), ph_data(:)
  integer(int32), allocatable :: snp_trans(:)
  integer(int32), allocatable :: bin_zeros(:), bin_ones(:), bin_zeros_trans(:), bin_ones_trans(:)
  real(real32), allocatable :: scores(:), scores_ref(:)
  integer :: p1, p2
  logical :: ok
  real(real64) :: elapsed

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <num_pac> <num_snp> <iteration>'
    stop 1
  end if

  num_pac = read_arg(1)
  num_snp = read_arg(2)
  iteration = read_arg(3)
  if (num_pac <= 0 .or. num_snp <= 1 .or. iteration <= 0) then
    print '(A)', 'invalid arguments'
    stop 1
  end if

  block_snp = 64
  allocate(snp_data(num_pac * num_snp), ph_data(num_pac), snp_trans(num_pac * num_snp))
  call fill_inputs(snp_data, ph_data, num_pac, num_snp)
  call transpose_snp(snp_data, snp_trans, num_pac, num_snp)

  phen_ones = count(ph_data == 1_int32)
  pp_zeros = (num_pac - phen_ones + 31) / 32
  pp_ones = (phen_ones + 31) / 32
  if (pp_zeros <= 0 .or. pp_ones <= 0) then
    print '(A)', 'invalid phenotype split'
    stop 1
  end if

  allocate(bin_zeros(num_snp * pp_zeros * 2), bin_ones(num_snp * pp_ones * 2))
  allocate(bin_zeros_trans(num_snp * pp_zeros * 2), bin_ones_trans(num_snp * pp_ones * 2))
  bin_zeros = 0_int32
  bin_ones = 0_int32
  call pack_snp(snp_trans, ph_data, bin_zeros, bin_ones, num_pac, num_snp, pp_zeros, pp_ones)
  mask_zeros = make_tail_mask(num_pac - phen_ones, pp_zeros)
  mask_ones = make_tail_mask(phen_ones, pp_ones)
  call transpose_binary(bin_zeros, bin_zeros_trans, num_snp, pp_zeros)
  call transpose_binary(bin_ones, bin_ones_trans, num_snp, pp_ones)

  allocate(scores(num_snp * num_snp), scores_ref(num_snp * num_snp))
  scores = score_max
  scores_ref = score_max

  elapsed = compute_scores_device(bin_zeros_trans, bin_ones_trans, scores, num_snp, pp_zeros, pp_ones, &
                                  mask_zeros, mask_ones, block_snp, iteration)
  write(*, '(A,F0.6,A)') 'Average kernel execution time: ', elapsed / real(iteration, real64), ' (s)'

  p1 = min_score(scores, num_snp, num_snp)
  call compute_scores_reference(bin_zeros_trans, bin_ones_trans, scores_ref, num_snp, pp_zeros, pp_ones, mask_zeros, mask_ones)
  p2 = min_score(scores_ref, num_snp, num_snp)
  ok = (p1 == p2) .and. (abs(scores(p1) - scores_ref(p2)) < 1.0e-3_real32)
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(snp_data, ph_data, snp_trans, bin_zeros, bin_ones, bin_zeros_trans, bin_ones_trans, scores, scores_ref)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_inputs(snp_data, ph_data, num_pac, num_snp)
    integer(int32), intent(out) :: snp_data(:), ph_data(:)
    integer, intent(in) :: num_pac, num_snp
    integer :: i, j
    call c_srand(100_c_int)
    do i = 0, num_pac - 1
      do j = 0, num_snp - 1
        snp_data(i * num_snp + j + 1) = int(mod(c_rand(), 3_c_int), int32)
      end do
    end do
    do i = 0, num_pac - 1
      ph_data(i + 1) = int(mod(c_rand(), 2_c_int), int32)
    end do
  end subroutine fill_inputs

  subroutine transpose_snp(snp_data, snp_trans, num_pac, num_snp)
    integer(int32), intent(in) :: snp_data(:)
    integer(int32), intent(out) :: snp_trans(:)
    integer, intent(in) :: num_pac, num_snp
    integer :: i, j
    do i = 0, num_pac - 1
      do j = 0, num_snp - 1
        snp_trans(j * num_pac + i + 1) = snp_data(i * num_snp + j + 1)
      end do
    end do
  end subroutine transpose_snp

  subroutine pack_snp(snp_trans, ph_data, bin_zeros, bin_ones, num_pac, num_snp, pp_zeros, pp_ones)
    integer(int32), intent(in) :: snp_trans(:), ph_data(:)
    integer(int32), intent(inout) :: bin_zeros(:), bin_ones(:)
    integer, intent(in) :: num_pac, num_snp, pp_zeros, pp_ones
    integer :: i, j, x_zeros, x_ones, n_zeros, n_ones, temp, base
    do i = 0, num_snp - 1
      x_zeros = -1
      x_ones = -1
      n_zeros = 0
      n_ones = 0
      do j = 0, num_pac - 1
        temp = int(snp_trans(i * num_pac + j + 1))
        if (ph_data(j + 1) == 1_int32) then
          if (mod(n_ones, 32) == 0) x_ones = x_ones + 1
          base = i * pp_ones * 2 + x_ones * 2
          bin_ones(base + 1) = shiftl(bin_ones(base + 1), 1)
          bin_ones(base + 2) = shiftl(bin_ones(base + 2), 1)
          if (temp == 0 .or. temp == 1) bin_ones(base + temp + 1) = ior(bin_ones(base + temp + 1), 1_int32)
          n_ones = n_ones + 1
        else
          if (mod(n_zeros, 32) == 0) x_zeros = x_zeros + 1
          base = i * pp_zeros * 2 + x_zeros * 2
          bin_zeros(base + 1) = shiftl(bin_zeros(base + 1), 1)
          bin_zeros(base + 2) = shiftl(bin_zeros(base + 2), 1)
          if (temp == 0 .or. temp == 1) bin_zeros(base + temp + 1) = ior(bin_zeros(base + temp + 1), 1_int32)
          n_zeros = n_zeros + 1
        end if
      end do
    end do
  end subroutine pack_snp

  integer function make_tail_mask(count_items, packets)
    integer, intent(in) :: count_items, packets
    integer :: x
    make_tail_mask = not(0_int32)
    do x = count_items, packets * 32 - 1
      make_tail_mask = shiftr(make_tail_mask, 1)
    end do
  end function make_tail_mask

  subroutine transpose_binary(src, dst, num_snp, packets)
    integer(int32), intent(in) :: src(:)
    integer(int32), intent(out) :: dst(:)
    integer, intent(in) :: num_snp, packets
    integer :: i, j
    do i = 0, num_snp - 1
      do j = 0, packets - 1
        dst((j * num_snp + i) * 2 + 1) = src((i * packets + j) * 2 + 1)
        dst((j * num_snp + i) * 2 + 2) = src((i * packets + j) * 2 + 2)
      end do
    end do
  end subroutine transpose_binary

  real(real64) function compute_scores_device(data_zeros, data_ones, scores, num_snp, pp_zeros, pp_ones, &
                                              mask_zeros, mask_ones, block_snp, iteration)
    integer(int32), intent(in) :: data_zeros(:), data_ones(:)
    real(real32), intent(inout) :: scores(:)
    integer, intent(in) :: num_snp, pp_zeros, pp_ones, mask_zeros, mask_ones, block_snp, iteration
    integer :: iter, num_snp_m, i, j, p, k, tid, base_i, base_j, step, n0, n1, n2
    integer(int32) :: ft(18), t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2, active_mask
    real(real32) :: score, g0, g1, g2
    real(real64) :: start_time, end_time

    num_snp_m = num_snp
    do while (mod(num_snp_m, block_snp) /= 0)
      num_snp_m = num_snp_m + 1
    end do

    !$omp target data map(to: data_zeros(1:num_snp*pp_zeros*2), data_ones(1:num_snp*pp_ones*2)) map(tofrom: scores(1:num_snp*num_snp))
      start_time = omp_get_wtime()
      do iter = 1, iteration
        !$omp target teams distribute parallel do collapse(2) thread_limit(block_snp) private(i, j, p, k, tid, ft, score, &
        !$omp& t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2, active_mask, base_i, base_j, step, n0, n1, n2, g0, g1, g2)
        do i = 0, num_snp_m - 1
          do j = 0, num_snp_m - 1
            score = score_max
            if (j > i .and. i < num_snp .and. j < num_snp) then
              ft = 0_int32
              base_i = i * 2
              base_j = j * 2
              step = 2 * num_snp
              do p = 0, 2 * pp_zeros * num_snp - step - 1, step
                active_mask = not(0_int32)
                di2 = iand(not(ior(data_zeros(base_i + p + 1), data_zeros(base_i + p + 2))), active_mask)
                dj2 = iand(not(ior(data_zeros(base_j + p + 1), data_zeros(base_j + p + 2))), active_mask)
                t00 = iand(data_zeros(base_i + p + 1), data_zeros(base_j + p + 1))
                t01 = iand(data_zeros(base_i + p + 1), data_zeros(base_j + p + 2))
                t02 = iand(data_zeros(base_i + p + 1), dj2)
                t10 = iand(data_zeros(base_i + p + 2), data_zeros(base_j + p + 1))
                t11 = iand(data_zeros(base_i + p + 2), data_zeros(base_j + p + 2))
                t12 = iand(data_zeros(base_i + p + 2), dj2)
                t20 = iand(di2, data_zeros(base_j + p + 1))
                t21 = iand(di2, data_zeros(base_j + p + 2))
                t22 = iand(di2, dj2)
                ft(1) = ft(1) + popcnt(t00)
                ft(2) = ft(2) + popcnt(t01)
                ft(3) = ft(3) + popcnt(t02)
                ft(4) = ft(4) + popcnt(t10)
                ft(5) = ft(5) + popcnt(t11)
                ft(6) = ft(6) + popcnt(t12)
                ft(7) = ft(7) + popcnt(t20)
                ft(8) = ft(8) + popcnt(t21)
                ft(9) = ft(9) + popcnt(t22)
              end do
              p = 2 * pp_zeros * num_snp - step
              active_mask = int(mask_zeros, int32)
              di2 = iand(not(ior(data_zeros(base_i + p + 1), data_zeros(base_i + p + 2))), active_mask)
              dj2 = iand(not(ior(data_zeros(base_j + p + 1), data_zeros(base_j + p + 2))), active_mask)
              t00 = iand(data_zeros(base_i + p + 1), data_zeros(base_j + p + 1))
              t01 = iand(data_zeros(base_i + p + 1), data_zeros(base_j + p + 2))
              t02 = iand(data_zeros(base_i + p + 1), dj2)
              t10 = iand(data_zeros(base_i + p + 2), data_zeros(base_j + p + 1))
              t11 = iand(data_zeros(base_i + p + 2), data_zeros(base_j + p + 2))
              t12 = iand(data_zeros(base_i + p + 2), dj2)
              t20 = iand(di2, data_zeros(base_j + p + 1))
              t21 = iand(di2, data_zeros(base_j + p + 2))
              t22 = iand(di2, dj2)
              ft(1) = ft(1) + popcnt(t00)
              ft(2) = ft(2) + popcnt(t01)
              ft(3) = ft(3) + popcnt(t02)
              ft(4) = ft(4) + popcnt(t10)
              ft(5) = ft(5) + popcnt(t11)
              ft(6) = ft(6) + popcnt(t12)
              ft(7) = ft(7) + popcnt(t20)
              ft(8) = ft(8) + popcnt(t21)
              ft(9) = ft(9) + popcnt(t22)

              do p = 0, 2 * pp_ones * num_snp - step - 1, step
                active_mask = not(0_int32)
                di2 = iand(not(ior(data_ones(base_i + p + 1), data_ones(base_i + p + 2))), active_mask)
                dj2 = iand(not(ior(data_ones(base_j + p + 1), data_ones(base_j + p + 2))), active_mask)
                t00 = iand(data_ones(base_i + p + 1), data_ones(base_j + p + 1))
                t01 = iand(data_ones(base_i + p + 1), data_ones(base_j + p + 2))
                t02 = iand(data_ones(base_i + p + 1), dj2)
                t10 = iand(data_ones(base_i + p + 2), data_ones(base_j + p + 1))
                t11 = iand(data_ones(base_i + p + 2), data_ones(base_j + p + 2))
                t12 = iand(data_ones(base_i + p + 2), dj2)
                t20 = iand(di2, data_ones(base_j + p + 1))
                t21 = iand(di2, data_ones(base_j + p + 2))
                t22 = iand(di2, dj2)
                ft(10) = ft(10) + popcnt(t00)
                ft(11) = ft(11) + popcnt(t01)
                ft(12) = ft(12) + popcnt(t02)
                ft(13) = ft(13) + popcnt(t10)
                ft(14) = ft(14) + popcnt(t11)
                ft(15) = ft(15) + popcnt(t12)
                ft(16) = ft(16) + popcnt(t20)
                ft(17) = ft(17) + popcnt(t21)
                ft(18) = ft(18) + popcnt(t22)
              end do
              p = 2 * pp_ones * num_snp - step
              active_mask = int(mask_ones, int32)
              di2 = iand(not(ior(data_ones(base_i + p + 1), data_ones(base_i + p + 2))), active_mask)
              dj2 = iand(not(ior(data_ones(base_j + p + 1), data_ones(base_j + p + 2))), active_mask)
              t00 = iand(data_ones(base_i + p + 1), data_ones(base_j + p + 1))
              t01 = iand(data_ones(base_i + p + 1), data_ones(base_j + p + 2))
              t02 = iand(data_ones(base_i + p + 1), dj2)
              t10 = iand(data_ones(base_i + p + 2), data_ones(base_j + p + 1))
              t11 = iand(data_ones(base_i + p + 2), data_ones(base_j + p + 2))
              t12 = iand(data_ones(base_i + p + 2), dj2)
              t20 = iand(di2, data_ones(base_j + p + 1))
              t21 = iand(di2, data_ones(base_j + p + 2))
              t22 = iand(di2, dj2)
              ft(10) = ft(10) + popcnt(t00)
              ft(11) = ft(11) + popcnt(t01)
              ft(12) = ft(12) + popcnt(t02)
              ft(13) = ft(13) + popcnt(t10)
              ft(14) = ft(14) + popcnt(t11)
              ft(15) = ft(15) + popcnt(t12)
              ft(16) = ft(16) + popcnt(t20)
              ft(17) = ft(17) + popcnt(t21)
              ft(18) = ft(18) + popcnt(t22)

              score = 0.0_real32
              do k = 0, 8
                n0 = ft(k + 1) + ft(10 + k) + 1
                n1 = ft(k + 1)
                n2 = ft(10 + k)
                if (n0 == 0) then
                  g0 = 0.0_real32
                else
                  g0 = (real(n0, real32) + 0.5_real32) * log(real(n0, real32)) - (real(n0, real32) - 1.0_real32)
                end if
                if (n1 == 0) then
                  g1 = 0.0_real32
                else
                  g1 = (real(n1, real32) + 0.5_real32) * log(real(n1, real32)) - (real(n1, real32) - 1.0_real32)
                end if
                if (n2 == 0) then
                  g2 = 0.0_real32
                else
                  g2 = (real(n2, real32) + 0.5_real32) * log(real(n2, real32)) - (real(n2, real32) - 1.0_real32)
                end if
                score = score + g0 - g1 - g2
              end do
              score = abs(score)
              if (score == 0.0_real32) score = score_max
              tid = i * num_snp + j + 1
              scores(tid) = score
            end if
          end do
        end do
        !$omp end target teams distribute parallel do
      end do
      end_time = omp_get_wtime()
    !$omp end target data
    compute_scores_device = end_time - start_time
  end function compute_scores_device

  subroutine compute_scores_reference(data_zeros, data_ones, scores, num_snp, pp_zeros, pp_ones, mask_zeros, mask_ones)
    integer(int32), intent(in) :: data_zeros(:), data_ones(:)
    real(real32), intent(inout) :: scores(:)
    integer, intent(in) :: num_snp, pp_zeros, pp_ones, mask_zeros, mask_ones
    integer :: i, j, k
    integer(int32) :: ft(18)
    real(real32) :: score
    do i = 0, num_snp - 1
      do j = 0, num_snp - 1
        if (j > i) then
          ft = 0_int32
          call accumulate_pair(data_zeros, ft, i, j, num_snp, pp_zeros, mask_zeros, 0)
          call accumulate_pair(data_ones, ft, i, j, num_snp, pp_ones, mask_ones, 9)
          score = 0.0_real32
          do k = 0, 8
            score = score + gammafunction(ft(k + 1) + ft(10 + k) + 1_int32) - &
              gammafunction(ft(k + 1)) - gammafunction(ft(10 + k))
          end do
          score = abs(score)
          if (score == 0.0_real32) score = score_max
          scores(i * num_snp + j + 1) = score
        end if
      end do
    end do
  end subroutine compute_scores_reference

  subroutine accumulate_pair(data, ft, i_snp, j_snp, num_snp, packets, mask, ft_offset)
    integer(int32), intent(in) :: data(:), mask
    integer(int32), intent(inout) :: ft(18)
    integer, intent(in) :: i_snp, j_snp, num_snp, packets, ft_offset
    integer :: p, base_i, base_j, step
    integer(int32) :: t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2

    base_i = i_snp * 2
    base_j = j_snp * 2
    step = 2 * num_snp
    do p = 0, 2 * packets * num_snp - step - 1, step
      call pair_counts(data, base_i + p, base_j + p, not(0_int32), t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2)
      call add_counts(ft, ft_offset, t00, t01, t02, t10, t11, t12, t20, t21, t22)
    end do
    p = 2 * packets * num_snp - step
    call pair_counts(data, base_i + p, base_j + p, int(mask, int32), t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2)
    call add_counts(ft, ft_offset, t00, t01, t02, t10, t11, t12, t20, t21, t22)
  end subroutine accumulate_pair

  subroutine pair_counts(data, idx_i, idx_j, mask, t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2)
    integer(int32), intent(in) :: data(:), mask
    integer, intent(in) :: idx_i, idx_j
    integer(int32), intent(out) :: t00, t01, t02, t10, t11, t12, t20, t21, t22, di2, dj2
    di2 = iand(not(ior(data(idx_i + 1), data(idx_i + 2))), mask)
    dj2 = iand(not(ior(data(idx_j + 1), data(idx_j + 2))), mask)
    t00 = iand(data(idx_i + 1), data(idx_j + 1))
    t01 = iand(data(idx_i + 1), data(idx_j + 2))
    t02 = iand(data(idx_i + 1), dj2)
    t10 = iand(data(idx_i + 2), data(idx_j + 1))
    t11 = iand(data(idx_i + 2), data(idx_j + 2))
    t12 = iand(data(idx_i + 2), dj2)
    t20 = iand(di2, data(idx_j + 1))
    t21 = iand(di2, data(idx_j + 2))
    t22 = iand(di2, dj2)
  end subroutine pair_counts

  subroutine add_counts(ft, ft_offset, t00, t01, t02, t10, t11, t12, t20, t21, t22)
    integer(int32), intent(inout) :: ft(18)
    integer, intent(in) :: ft_offset
    integer(int32), intent(in) :: t00, t01, t02, t10, t11, t12, t20, t21, t22
    ft(ft_offset + 1) = ft(ft_offset + 1) + popcnt(t00)
    ft(ft_offset + 2) = ft(ft_offset + 2) + popcnt(t01)
    ft(ft_offset + 3) = ft(ft_offset + 3) + popcnt(t02)
    ft(ft_offset + 4) = ft(ft_offset + 4) + popcnt(t10)
    ft(ft_offset + 5) = ft(ft_offset + 5) + popcnt(t11)
    ft(ft_offset + 6) = ft(ft_offset + 6) + popcnt(t12)
    ft(ft_offset + 7) = ft(ft_offset + 7) + popcnt(t20)
    ft(ft_offset + 8) = ft(ft_offset + 8) + popcnt(t21)
    ft(ft_offset + 9) = ft(ft_offset + 9) + popcnt(t22)
  end subroutine add_counts

  real(real32) function gammafunction(n)
    integer(int32), intent(in) :: n
    if (n == 0_int32) then
      gammafunction = 0.0_real32
    else
      gammafunction = (real(n, real32) + 0.5_real32) * log(real(n, real32)) - (real(n, real32) - 1.0_real32)
    end if
  end function gammafunction

  integer function min_score(scores, nrows, ncols)
    real(real32), intent(in) :: scores(:)
    integer, intent(in) :: nrows, ncols
    real(real32) :: score
    integer :: idx
    score = scores(1)
    min_score = 1
    do idx = 2, nrows * ncols
      if (score > scores(idx)) then
        score = scores(idx)
        min_score = idx
      end if
    end do
  end function min_score

end program main
