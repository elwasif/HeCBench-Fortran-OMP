! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_double, c_int64_t
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: ncases = 3
  integer, parameter :: exact_cases(ncases) = [1000, 0, 1000]
  integer, parameter :: sampled_cases(ncases) = [0, 1000, 1000]
  integer, parameter :: ncols_cases(ncases) = [2000, 2000, 2000]
  integer, parameter :: background_cases(ncases) = [10, 10, 10]
  integer, parameter :: max_samples_cases(ncases) = [11, 11, 11]
  integer(c_int64_t), parameter :: seed_cases(ncases) = [1234_c_int64_t, 1234_c_int64_t, 1234_c_int64_t]
  integer :: repeat, case_id, r
  real(real64) :: total_time

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <repeat>'
    stop 1
  end if
  repeat = read_arg(1)

  do case_id = 1, ncases
    total_time = 0.0_real64
    do r = 1, repeat
      call run_case(exact_cases(case_id), sampled_cases(case_id), ncols_cases(case_id), &
                    background_cases(case_id), max_samples_cases(case_id), seed_cases(case_id), total_time)
    end do
    write(*, '(A,F0.6,A)') 'Average execution time of kernels: ', (total_time * 1.0e6_real64) / real(repeat, real64), ' (us)'
  end do

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  function LCG_random_double(seed) result(random_value)
    integer(c_int64_t), intent(inout) :: seed
    real(c_double) :: random_value
    integer(c_int64_t), parameter :: a = 2806196910506780709_c_int64_t

    seed = a * seed + 1_c_int64_t
    seed = iand(seed, not(ishft(1_c_int64_t, 63)))
    random_value = real(seed, c_double) / real(ishft(1_c_int64_t, 62), c_double) / 2.0_c_double
  end function LCG_random_double

  subroutine run_case(nrows_exact, nrows_sampled, ncols, nrows_background, max_samples, seed, total_time)
    integer, intent(in) :: nrows_exact, nrows_sampled, ncols, nrows_background, max_samples
    integer(c_int64_t), intent(in) :: seed
    real(real64), intent(inout) :: total_time
    integer :: nrows_x, i, j
    real(real32) :: sent_value
    real(real32), allocatable :: background(:), observation(:), x(:), dataset(:)
    integer, allocatable :: nsamples(:)
    real(real64) :: start_time, end_time

    nrows_x = nrows_exact + nrows_sampled
    allocate(background(nrows_background * ncols), observation(ncols), x(nrows_x * ncols))
    allocate(dataset(nrows_x * nrows_background * ncols))
    allocate(nsamples(max(1, nrows_sampled / 2)))

    sent_value = real(nrows_x * nrows_background * ncols * 100, real32)
    observation = sent_value
    do i = 0, nrows_background - 1
      do j = 0, ncols - 1
        background(i * ncols + j + 1) = real(i * 2 + 1, real32)
      end do
    end do

    x = 0.0_real32
    do i = 0, nrows_exact - 1
      x(i * ncols + i + 1) = 1.0_real32
      x(i * ncols + i + 2) = 1.0_real32
    end do

    do i = 0, nrows_sampled / 2 - 1
      nsamples(i + 1) = max_samples - mod(i, 2)
    end do
    dataset = 0.0_real32

    !$omp target data map(to: background(1:nrows_background*ncols), observation(1:ncols), nsamples(1:max(1,nrows_sampled/2))) &
    !$omp& map(tofrom: x(1:nrows_x*ncols)) map(from: dataset(1:nrows_x*nrows_background*ncols))
      start_time = omp_get_wtime()
      if (nrows_exact > 0) then
        call fill_exact(background, observation, x, dataset, nrows_exact, ncols, nrows_background)
      end if
      if (nrows_sampled > 0) then
        call fill_sampled(background, observation, nsamples, x, dataset, nrows_exact, nrows_sampled, ncols, nrows_background, seed)
      end if
      end_time = omp_get_wtime()
    !$omp end target data
    total_time = total_time + (end_time - start_time)

    call validate_case(x, dataset, nsamples, sent_value, nrows_exact, nrows_sampled, ncols, nrows_background)
    deallocate(background, observation, x, dataset, nsamples)
  end subroutine run_case

  subroutine fill_exact(background, observation, x, dataset, nrows_exact, ncols, nrows_background)
    real(real32), intent(in) :: background(:), observation(:), x(:)
    real(real32), intent(out) :: dataset(:)
    integer, intent(in) :: nrows_exact, ncols, nrows_background
    integer :: gid, col, row, row_idx, curr_x, nthreads

    nthreads = min(256, ncols)
    !$omp target teams num_teams(nrows_exact) thread_limit(nthreads)
    !$omp parallel private(gid, col, row, row_idx, curr_x)
      gid = omp_get_team_num()
      col = omp_get_thread_num()
      row = gid * ncols
      do while (col < ncols)
        curr_x = int(x(row + col + 1))
        do row_idx = gid * nrows_background, gid * nrows_background + nrows_background - 1
          if (curr_x == 0) then
            dataset(row_idx * ncols + col + 1) = background(mod(row_idx, nrows_background) * ncols + col + 1)
          else
            dataset(row_idx * ncols + col + 1) = observation(col + 1)
          end if
        end do
        col = col + omp_get_num_threads()
      end do
    !$omp end parallel
    !$omp end target teams
  end subroutine fill_exact

  subroutine fill_sampled(background, observation, nsamples, x, dataset, nrows_exact, nrows_sampled, ncols, nrows_background, seed)
    real(real32), intent(in) :: background(:), observation(:)
    integer, intent(in) :: nsamples(:)
    real(real32), intent(inout) :: x(:)
    real(real32), intent(inout) :: dataset(:)
    integer, intent(in) :: nrows_exact, nrows_sampled, ncols, nrows_background
    integer(c_int64_t), intent(in) :: seed
    integer :: bid, tid, k_blk, rand_idx, col_idx, curr_x, bg_row_idx, nthreads
    integer(c_int64_t) :: seed_value
    real(real32) :: old_x

    nthreads = min(256, ncols)
    !$omp target teams num_teams(nrows_sampled / 2) thread_limit(nthreads)
    !$omp parallel private(bid, tid, k_blk, rand_idx, col_idx, curr_x, bg_row_idx, old_x, seed_value)
      bid = omp_get_team_num()
      tid = omp_get_thread_num()
      seed_value = seed
      k_blk = nsamples(bid + 1)
      if (tid < k_blk) then
        rand_idx = int(LCG_random_double(seed_value) * real(ncols, c_double))
        do
          !$omp atomic capture
          old_x = x((nrows_exact + 2 * bid) * ncols + rand_idx + 1)
          x((nrows_exact + 2 * bid) * ncols + rand_idx + 1) = 1.0_real32
          !$omp end atomic
          if (old_x == 0.0_real32) exit
          rand_idx = int(LCG_random_double(seed_value) * real(ncols, c_double))
        end do
      end if
      !$omp barrier

      col_idx = tid
      do while (col_idx < ncols)
        curr_x = int(x((nrows_exact + 2 * bid) * ncols + col_idx + 1))
        x((nrows_exact + 2 * bid + 1) * ncols + col_idx + 1) = real(1 - curr_x, real32)
        do bg_row_idx = 2 * bid * nrows_background, 2 * bid * nrows_background + nrows_background - 1
          if (curr_x == 0) then
            dataset((nrows_exact * nrows_background + bg_row_idx) * ncols + col_idx + 1) = &
              background(mod(bg_row_idx, nrows_background) * ncols + col_idx + 1)
          else
            dataset((nrows_exact * nrows_background + bg_row_idx) * ncols + col_idx + 1) = observation(col_idx + 1)
          end if
        end do

        do bg_row_idx = (2 * bid + 1) * nrows_background, (2 * bid + 2) * nrows_background - 1
          if (curr_x == 0) then
            dataset((nrows_exact * nrows_background + bg_row_idx) * ncols + col_idx + 1) = observation(col_idx + 1)
          else
            dataset((nrows_exact * nrows_background + bg_row_idx) * ncols + col_idx + 1) = &
              background(mod(bg_row_idx, nrows_background) * ncols + col_idx + 1)
          end if
        end do
        col_idx = col_idx + omp_get_num_threads()
      end do
    !$omp end parallel
    !$omp end target teams
  end subroutine fill_sampled

  subroutine validate_case(x, dataset, nsamples, sent_value, nrows_exact, nrows_sampled, ncols, nrows_background)
    real(real32), intent(in) :: x(:), dataset(:), sent_value
    integer, intent(in) :: nsamples(:), nrows_exact, nrows_sampled, ncols, nrows_background
    logical :: test_sampled_x, test_scatter_exact, test_scatter_sampled
    integer :: i, j, k, counter, sample_idx, compliment_ctr

    test_sampled_x = .true.
    sample_idx = 1
    do i = nrows_exact * ncols + 1, (nrows_exact + nrows_sampled) * ncols, 2 * ncols
      if (nrows_sampled <= 0) exit
      counter = 0
      do k = i, i + ncols - 1
        if (x(k) == 1.0_real32) counter = counter + 1
      end do
      test_sampled_x = test_sampled_x .and. (counter == nsamples(sample_idx))
      counter = 0
      do k = i + ncols, i + 2 * ncols - 1
        if (x(k) == 1.0_real32) counter = counter + 1
      end do
      test_sampled_x = test_sampled_x .and. (counter == ncols - nsamples(sample_idx))
      sample_idx = sample_idx + 1
    end do

    test_scatter_exact = .true.
    do i = 0, nrows_exact - 1
      do j = i * nrows_background * ncols + 1, (i + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_exact = test_scatter_exact .and. (counter == 2)
        if (.not. test_scatter_exact) then
          write(*, '(A,I0,A)') 'test_scatter_exact counter failed with: ', counter, ', expected value was 2.'
          exit
        end if
      end do
      if (.not. test_scatter_exact) exit
    end do

    test_scatter_sampled = .true.
    compliment_ctr = 0
    do i = nrows_exact, nrows_exact + nrows_sampled / 2 - 1
      do j = (i + compliment_ctr) * nrows_background * ncols + 1, (i + compliment_ctr + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_sampled = test_scatter_sampled .and. (counter == nsamples(i - nrows_exact + 1))
        if (.not. test_scatter_sampled) then
          write(*, '(A,I0,A,I0,A)') 'test_scatter_sampled counter failed with: ', counter, &
            ', expected value was ', nsamples(i - nrows_exact + 1), '.'
          exit
        end if
      end do

      compliment_ctr = compliment_ctr + 1
      do j = (i + compliment_ctr) * nrows_background * ncols + 1, (i + compliment_ctr + 1) * nrows_background * ncols, ncols
        counter = 0
        do k = j, j + ncols - 1
          if (dataset(k) == sent_value) counter = counter + 1
        end do
        test_scatter_sampled = test_scatter_sampled .and. (counter == ncols - nsamples(i - nrows_exact + 1))
        if (.not. test_scatter_sampled) then
          write(*, '(A,I0,A,I0,A)') 'test_scatter_sampled counter failed with: ', counter, &
            ', expected value was ', ncols - nsamples(i - nrows_exact + 1), '.'
          exit
        end if
      end do
      if (.not. test_scatter_sampled) exit
    end do

    if (.not. test_sampled_x) print '(A)', 'test_sampled_X failed'
  end subroutine validate_case

end program main
