! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer, parameter :: row_size = 17
  real(real64), parameter :: lower_limit = 0.0_real64
  real(real64), parameter :: upper_limit = 15.0_real64
  real(real64), parameter :: eps = 1.0e-7_real64

  character(len=128) :: arg1, arg2, arg3
  integer :: nwg, wgs, repeat, iter, k
  real(real64), allocatable :: result(:)
  real(real64) :: start_time, end_time, elapsed_s, d_sum, ref_sum

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <number of work-groups> <work-group size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) nwg
  read(arg2, *) wgs
  read(arg3, *) repeat

  allocate(result(nwg))
  result = 0.0_real64
  d_sum = 0.0_real64

  !$omp target data map(from: result(1:nwg))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call compute_segments(result, nwg, wgs)
    !$omp target update from(result(1:nwg))
    d_sum = 0.0_real64
    do k = 1, nwg
      d_sum = d_sum + result(k)
    end do
  end do
  end_time = omp_get_wtime()
  elapsed_s = (end_time - start_time) / real(repeat, real64)
  print '(A,F0.6,A)', 'Average kernel execution time: ', elapsed_s, ' (s)'
  !$omp end target data

  ref_sum = romberg_reference(lower_limit, upper_limit, row_size, eps)
  if (abs(d_sum - ref_sum) > eps) then
    print '(A)', 'FAIL'
  else
    print '(A)', 'PASS'
  end if

  deallocate(result)

contains

  subroutine compute_segments(result, nwg, wgs)
    real(real64), intent(inout) :: result(:)
    integer, intent(in) :: nwg, wgs
    real(real64) :: a, b

    a = lower_limit
    b = upper_limit
    !$omp target teams num_teams(nwg) thread_limit(wgs) firstprivate(a, b, nwg)
    block
      real(real64) :: smem(row_size * 64)

      !$omp parallel
      block
        integer :: threadIdx_x, blockIdx_x, gridDim_x, blockDim_x
        integer :: i, k, col, row, max_eval
        real(real64) :: diff, step, sum, local_col(row_size), a_seg, b_seg

        threadIdx_x = omp_get_thread_num()
        blockIdx_x = omp_get_team_num()
        gridDim_x = omp_get_num_teams()
        blockDim_x = omp_get_num_threads()
        diff = (b - a) / real(gridDim_x, real64)
        max_eval = ishft(1, row_size - 1)
        b_seg = a + real(blockIdx_x + 1, real64) * diff
        a_seg = a + real(blockIdx_x, real64) * diff

        step = (b_seg - a_seg) / real(max_eval, real64)

        local_col = 0.0_real64
        if (threadIdx_x == 0) then
          k = blockDim_x
          local_col(1) = integrand(a_seg) + integrand(b_seg)
        else
          k = threadIdx_x
        end if

        do while (k < max_eval)
          local_col(row_size - getFirstSetBitPos(k) + 1) = local_col(row_size - getFirstSetBitPos(k) + 1) + &
              2.0_real64 * integrand(a_seg + step * real(k, real64))
          k = k + blockDim_x
        end do

        do i = 1, row_size
          smem(row_size * threadIdx_x + i) = local_col(i)
        end do
        !$omp barrier

        if (threadIdx_x < row_size) then
          sum = 0.0_real64
          do i = threadIdx_x, blockDim_x * row_size - 1, row_size
            sum = sum + smem(i + 1)
          end do
          smem(threadIdx_x + 1) = sum
        end if
        !$omp barrier

        if (threadIdx_x == 0) then
          local_col(1) = smem(1)

          do k = 2, row_size
            local_col(k) = local_col(k - 1) + smem(k)
          end do

          do k = 1, row_size
            local_col(k) = local_col(k) * (b_seg - a_seg) / real(ishft(1, k), real64)
          end do

          do col = 1, row_size - 1
            do row = row_size, col + 1, -1
              local_col(row) = local_col(row) + (local_col(row) - local_col(row - 1)) / &
                  real(ishft(1, 2 * col - 1) - 1, real64)
            end do
          end do

          result(blockIdx_x + 1) = local_col(row_size)
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine compute_segments

  integer function getFirstSetBitPos(n)
    integer, intent(in) :: n
    integer :: value

    value = iand(n, -n)
    getFirstSetBitPos = 1
    do while (value > 1)
      value = ishft(value, -1)
      getFirstSetBitPos = getFirstSetBitPos + 1
    end do
  end function getFirstSetBitPos

  real(real64) function integrand(x)
    real(real64), intent(in) :: x

    integrand = exp(x) * sin(x)
  end function integrand

  real(real64) function romberg_reference(a, b, max_steps, acc)
    real(real64), intent(in) :: a, b, acc
    integer, intent(in) :: max_steps
    real(real64) :: r_prev(row_size), r_curr(row_size)
    real(real64) :: h, c, n_k
    integer :: i, j, ep

    r_prev = 0.0_real64
    r_curr = 0.0_real64
    h = b - a
    r_prev(1) = (integrand(a) + integrand(b)) * h * 0.5_real64

    do i = 2, max_steps
      h = h * 0.5_real64
      c = 0.0_real64
      ep = ishft(1, i - 2)
      do j = 1, ep
        c = c + integrand(a + real(2 * j - 1, real64) * h)
      end do
      r_curr(1) = h * c + 0.5_real64 * r_prev(1)

      do j = 2, i
        n_k = real(ishft(1, 2 * (j - 1)), real64)
        r_curr(j) = (n_k * r_curr(j - 1) - r_prev(j - 1)) / (n_k - 1.0_real64)
      end do

      if (i > 3 .and. abs(r_prev(i - 1) - r_curr(i)) < acc) then
        romberg_reference = r_curr(i - 1)
        return
      end if

      r_prev = r_curr
    end do

    romberg_reference = r_prev(max_steps)
  end function romberg_reference

end program main
