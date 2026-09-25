! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: vector_size = 8 * 1024 * 1024
  integer, parameter :: granularity = 8
  integer, parameter :: fusion_degree = 4
  integer, parameter :: block_dim = 256
  real(real32), parameter :: seed = 0.1_real32

  character(len=128) :: arg1, arg2
  integer :: compute_iterations, repeat
  integer :: datasize_mb
  logical :: ok

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage: ./main <compute iterations> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) compute_iterations
  read(arg2, *) repeat

  datasize_mb = vector_size * 4 / (1024 * 1024)
  print '(A,I0,A)', 'Buffer size: ', datasize_mb, 'MB'
  call mixbench_gpu(vector_size, compute_iterations, repeat, ok)
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

contains

  subroutine mixbench_gpu(size, compute_iterations, repeat, ok)
    integer, intent(in) :: size, compute_iterations, repeat
    logical, intent(out) :: ok
    real(real32), allocatable :: cd(:)
    integer :: reduced_grid_size, grid_dim, i
    real(real64) :: start_time, end_time, elapsed_s

    print '(A)', 'Trade-off type:compute with global memory (block strided)'
    allocate(cd(size))
    cd = 0.0_real32

    reduced_grid_size = size / granularity / 128
    grid_dim = reduced_grid_size / block_dim

    !$omp target data map(tofrom: cd(1:size))
    do i = 1, repeat
      call benchmark_func(cd, grid_dim, block_dim, compute_iterations)
    end do

    start_time = omp_get_wtime()
    do i = 1, repeat
      call benchmark_func(cd, grid_dim, block_dim, compute_iterations)
    end do
    end_time = omp_get_wtime()
    elapsed_s = end_time - start_time
    print '(A,F0.6,A)', 'Total kernel execution time: ', elapsed_s, ' (s)'
    !$omp end target data

    ok = .true.
    do i = 1, size
      if (cd(i) /= 0.0_real32) then
        if (abs(cd(i) - 0.050807_real32) > 1.0e-6_real32) then
          ok = .false.
          print '(A,I0,A,F0.6)', 'Verification failed at index ', i - 1, ': ', cd(i)
          exit
        end if
      end if
    end do

    deallocate(cd)
  end subroutine mixbench_gpu

  subroutine benchmark_func(cd, grid_dim, block_dim, compute_iterations)
    real(real32), intent(inout) :: cd(:)
    integer, intent(in) :: grid_dim, block_dim, compute_iterations
    integer :: team, thread, k, j, iter, idx, stride, big_stride
    real(real32) :: tmps(granularity), sum_value

    !$omp target teams num_teams(grid_dim) thread_limit(block_dim) private(team, thread, k, j, iter, idx, stride, big_stride, tmps, sum_value)
    !$omp parallel private(team, thread, k, j, iter, idx, stride, big_stride, tmps, sum_value)
    team = omp_get_team_num()
    thread = omp_get_thread_num()
    stride = block_dim
    idx = team * block_dim * granularity + thread + 1
    big_stride = grid_dim * block_dim * granularity

    do k = 0, fusion_degree - 1
      do j = 0, granularity - 1
        tmps(j + 1) = cd(idx + j * stride + k * big_stride)
        do iter = 1, compute_iterations
          tmps(j + 1) = tmps(j + 1) * tmps(j + 1) + seed
        end do
      end do

      sum_value = 0.0_real32
      do j = 1, granularity, 2
        sum_value = sum_value + tmps(j) * tmps(j + 1)
      end do

      do j = 0, granularity - 1
        cd(idx + k * big_stride) = sum_value
      end do
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine benchmark_func

end program main
