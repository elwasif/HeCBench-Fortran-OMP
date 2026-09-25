! SPDX-License-Identifier: CC0-1.0
program stddev_main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer :: d, n, repeat, i
  character(len=64) :: arg
  real(real32), allocatable :: data(:), std(:), std_ref(:)
  real(real64) :: start_time, elapsed_s
  logical :: ok
  integer(c_int), parameter :: rand_max = 2147483647_c_int

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name="rand")
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <D> <N> <repeat>")')
    write(*,'("D: number of columns of data (must be a multiple of 32)")')
    write(*,'("N: number of rows of data (at least one row)")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) d
  call get_command_argument(2, arg)
  read(arg, *) n
  call get_command_argument(3, arg)
  read(arg, *) repeat

  allocate(data(0:d*n-1), std(0:d-1), std_ref(0:d-1))
  call c_srand(123_c_int)
  do i = 0, d * n - 1
    data(i) = real(c_rand(), real32) / real(rand_max, real32)
  end do

  !$omp target data map(to: data) map(from: std)
  call stddev_kernel(std, data, d, n, .true.)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call stddev_kernel(std, data, d, n, .true.)
  end do
  elapsed_s = (omp_get_wtime() - start_time) / real(repeat, real64)
  !$omp end target data

  write(*,'("Average execution time of stddev kernels: ",F8.6," (s)")') elapsed_s

  call stddev_ref_kernel(std_ref, data, d, n, .true.)
  ok = maxval(abs(std - std_ref)) <= 1.0e-3_real32
  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine stddev_kernel(std, data, d, n, sample)
    real(real32), intent(out) :: std(0:)
    real(real32), intent(in) :: data(0:)
    integer, intent(in) :: d, n
    logical, intent(in) :: sample
    integer, parameter :: tpb = 256
    integer, parameter :: rows_per_thread = 4
    integer, parameter :: cols_per_blk = 32
    integer, parameter :: rows_per_blk = (tpb / cols_per_blk) * rows_per_thread
    integer :: team_x, team_y, teams, rows_per_blk_per_iter
    integer :: tx, bx, by, grid_dim_x
    integer :: this_col_id, this_row_id, col_id, row_id, stride
    integer :: col, i, sample_size
    real(real32) :: thread_data, val
    real(real32) :: sstd(0:cols_per_blk-1)

    team_x = (n + rows_per_blk - 1) / rows_per_blk
    team_y = (d + cols_per_blk - 1) / cols_per_blk
    teams = team_x * team_y
    rows_per_blk_per_iter = tpb / cols_per_blk
    sample_size = merge(n - 1, n, sample)

    !$omp target teams distribute parallel do thread_limit(256)
    do col = 0, d - 1
      std(col) = 0.0_real32
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams num_teams(teams) thread_limit(tpb) private(sstd)
    !$omp parallel private(tx, bx, by, grid_dim_x, this_col_id, this_row_id, col_id, row_id, stride, i, val, thread_data)
    tx = omp_get_thread_num()
    bx = mod(omp_get_team_num(), team_x)
    by = omp_get_team_num() / team_x
    grid_dim_x = team_x

    this_col_id = mod(tx, cols_per_blk)
    this_row_id = tx / cols_per_blk
    col_id = this_col_id + by * cols_per_blk
    row_id = this_row_id + bx * rows_per_blk_per_iter
    thread_data = 0.0_real32
    stride = rows_per_blk_per_iter * grid_dim_x

    do i = row_id, n - 1, stride
      if (col_id < d) then
        val = data(i * d + col_id)
      else
        val = 0.0_real32
      end if
      thread_data = thread_data + val * val
    end do

    if (tx < cols_per_blk) sstd(tx) = 0.0_real32
    !$omp barrier

    !$omp atomic update
    sstd(this_col_id) = sstd(this_col_id) + thread_data
    !$omp barrier

    if (tx < cols_per_blk .and. col_id < d) then
      !$omp atomic update
      std(col_id) = std(col_id) + sstd(this_col_id)
    end if
    !$omp end parallel
    !$omp end target teams

    !$omp target teams distribute parallel do thread_limit(tpb)
    do col = 0, d - 1
      std(col) = sqrt(std(col) / real(sample_size, real32))
    end do
    !$omp end target teams distribute parallel do
  end subroutine stddev_kernel

  subroutine stddev_ref_kernel(std, data, d, n, sample)
    real(real32), intent(out) :: std(0:)
    real(real32), intent(in) :: data(0:)
    integer, intent(in) :: d, n
    logical, intent(in) :: sample
    integer :: col, row, sample_size
    real(real32) :: sumsq

    sample_size = merge(n - 1, n, sample)
    do col = 0, d - 1
      sumsq = 0.0_real32
      do row = 0, n - 1
        sumsq = sumsq + data(row * d + col) * data(row * d + col)
      end do
      std(col) = sqrt(sumsq / real(sample_size, real32))
    end do
  end subroutine stddev_ref_kernel

end program stddev_main
