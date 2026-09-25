! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: size, repeat
  integer, parameter :: block_size_x = 32
  integer, parameter :: block_size_y = 8
  integer, parameter :: tile_size_x = 4
  integer, parameter :: tile_size_y = 4
  integer, parameter :: reduce_block_size = 256

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage ./', trim(arg0), ' <size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) size
  read(arg2, *) repeat

  print '(A)', 'Test single precision'
  call test_real32(size, repeat)

  print '(A)', 'Test double precision'
  call test_real64(size, repeat)

contains

  subroutine test_real32(size, repeat)
    integer, intent(in) :: size, repeat
    integer :: i, nblocks, nteams_x, nteams_y
    real(real32), allocatable :: a(:), b(:), scale_a(:), scale_b(:), cost(:)
    real(real32) :: output, host_output
    real(real64) :: start_time, end_time, avg_time, analytical

    nteams_x = size / (block_size_x * tile_size_x)
    nteams_y = size / (block_size_y * tile_size_y)
    nblocks = nteams_x * nteams_y

    allocate(a(size * 2), b(size * 2), scale_a(size), scale_b(size), cost(nblocks))
    a = 1.0_real32
    b = 0.0_real32
    scale_a = 1.0_real32
    scale_b = 1.0_real32

    !$omp target data map(to: a(1:size*2), b(1:size*2), scale_a(1:size), scale_b(1:size)) &
    !$omp& map(alloc: cost(1:nblocks))
    start_time = omp_get_wtime()
    do i = 1, repeat
      call distance_real32(a, b, size, size, scale_a, scale_b, cost, nteams_x, nblocks)
      output = reduce_cross_term_real32(cost, size, size, nblocks)
    end do
    end_time = omp_get_wtime()
    avg_time = (end_time - start_time) / real(repeat, real64)
    print '(A,F8.6,A)', 'Average kernel execution time ', avg_time, ' (s)'
    print '(A,F0.6)', '    device result: ', real(output, real64)

    host_output = host_cost_real32(a, b, scale_a, scale_b, size)
    print '(A,F0.6)', '      host result: ', real(host_output, real64)
    analytical = real(size, real64) * real(size, real64) * exp(-1.0_real64)
    print '(A,F0.6)', 'analytical result: ', analytical
    print '(A)'
    !$omp end target data

    deallocate(a, b, scale_a, scale_b, cost)
  end subroutine test_real32

  subroutine test_real64(size, repeat)
    integer, intent(in) :: size, repeat
    integer :: i, nblocks, nteams_x, nteams_y
    real(real64), allocatable :: a(:), b(:), scale_a(:), scale_b(:), cost(:)
    real(real64) :: output, host_output, start_time, end_time, avg_time, analytical

    nteams_x = size / (block_size_x * tile_size_x)
    nteams_y = size / (block_size_y * tile_size_y)
    nblocks = nteams_x * nteams_y

    allocate(a(size * 2), b(size * 2), scale_a(size), scale_b(size), cost(nblocks))
    a = 1.0_real64
    b = 0.0_real64
    scale_a = 1.0_real64
    scale_b = 1.0_real64

    !$omp target data map(to: a(1:size*2), b(1:size*2), scale_a(1:size), scale_b(1:size)) &
    !$omp& map(alloc: cost(1:nblocks))
    start_time = omp_get_wtime()
    do i = 1, repeat
      call distance_real64(a, b, size, size, scale_a, scale_b, cost, nteams_x, nblocks)
      output = reduce_cross_term_real64(cost, size, size, nblocks)
    end do
    end_time = omp_get_wtime()
    avg_time = (end_time - start_time) / real(repeat, real64)
    print '(A,F8.6,A)', 'Average kernel execution time ', avg_time, ' (s)'
    print '(A,F0.6)', '    device result: ', output

    host_output = host_cost_real64(a, b, scale_a, scale_b, size)
    print '(A,F0.6)', '      host result: ', host_output
    analytical = real(size, real64) * real(size, real64) * exp(-1.0_real64)
    print '(A,F0.6)', 'analytical result: ', analytical
    print '(A)'
    !$omp end target data

    deallocate(a, b, scale_a, scale_b, cost)
  end subroutine test_real64

  subroutine distance_real32(a, b, m, n, scale_a, scale_b, cross_term, nteams_x, nblocks)
    real(real32), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: m, n, nteams_x, nblocks
    real(real32), intent(out) :: cross_term(:)
    real(real32) :: sh_a(2, block_size_x * tile_size_x)
    real(real32) :: sh_b(2, block_size_y * tile_size_y)
    real(real32) :: sh_scale_a(block_size_x * tile_size_x)
    real(real32) :: sh_scale_b(block_size_y * tile_size_y)
    real(real32) :: sum, s_cross_term, dist_ij
    integer :: tx, ty, bx, by, ii, jj, ti, tj, d, team

    !$omp target teams num_teams(nblocks) thread_limit(block_size_x * block_size_y) &
    !$omp& private(sh_a, sh_b, sh_scale_a, sh_scale_b, sum)
    sum = 0.0_real32
    !$omp parallel private(tx, ty, bx, by, ii, jj, ti, tj, d, s_cross_term, dist_ij, team)
    tx = mod(omp_get_thread_num(), block_size_x)
    ty = omp_get_thread_num() / block_size_x
    team = omp_get_team_num()
    bx = mod(team, nteams_x)
    by = team / nteams_x
    ii = tx + bx * block_size_x * tile_size_x
    jj = ty + by * block_size_y * tile_size_y

    do d = 1, 2
      do ti = 0, tile_size_x - 1
        sh_a(d, tx + ti * block_size_x + 1) = a(ii + ti * block_size_x + 1 + (d - 1) * m)
      end do
      if (tx == 0) then
        do tj = 0, tile_size_y - 1
          sh_b(d, ty + tj * block_size_y + 1) = b(jj + tj * block_size_y + 1 + (d - 1) * n)
        end do
      end if
    end do
    do ti = 0, tile_size_x - 1
      sh_scale_a(tx + ti * block_size_x + 1) = scale_a(ii + ti * block_size_x + 1)
    end do
    if (tx == 0) then
      do tj = 0, tile_size_y - 1
        sh_scale_b(ty + tj * block_size_y + 1) = scale_b(jj + tj * block_size_y + 1)
      end do
    end if

    s_cross_term = 0.0_real32
    do ti = 0, tile_size_x - 1
      do tj = 0, tile_size_y - 1
        if (ii + ti * block_size_x < m .and. jj + tj * block_size_y < n) then
          dist_ij = 0.0_real32
          do d = 1, 2
            dist_ij = dist_ij + &
              (sh_a(d, tx + ti * block_size_x + 1) - sh_b(d, ty + tj * block_size_y + 1)) * &
              (sh_a(d, tx + ti * block_size_x + 1) - sh_b(d, ty + tj * block_size_y + 1))
          end do
          s_cross_term = s_cross_term + &
            exp(-dist_ij / (sh_scale_a(tx + ti * block_size_x + 1) + sh_scale_b(ty + tj * block_size_y + 1)))
        end if
      end do
    end do

    !$omp atomic update
    sum = sum + s_cross_term
    !$omp barrier
    if (tx == 0 .and. ty == 0) cross_term(by * nteams_x + bx + 1) = sum
    !$omp end parallel
    !$omp end target teams
  end subroutine distance_real32

  subroutine distance_real64(a, b, m, n, scale_a, scale_b, cross_term, nteams_x, nblocks)
    real(real64), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: m, n, nteams_x, nblocks
    real(real64), intent(out) :: cross_term(:)
    real(real64) :: sh_a(2, block_size_x * tile_size_x)
    real(real64) :: sh_b(2, block_size_y * tile_size_y)
    real(real64) :: sh_scale_a(block_size_x * tile_size_x)
    real(real64) :: sh_scale_b(block_size_y * tile_size_y)
    real(real64) :: sum, s_cross_term, dist_ij
    integer :: tx, ty, bx, by, ii, jj, ti, tj, d, team

    !$omp target teams num_teams(nblocks) thread_limit(block_size_x * block_size_y) &
    !$omp& private(sh_a, sh_b, sh_scale_a, sh_scale_b, sum)
    sum = 0.0_real64
    !$omp parallel private(tx, ty, bx, by, ii, jj, ti, tj, d, s_cross_term, dist_ij, team)
    tx = mod(omp_get_thread_num(), block_size_x)
    ty = omp_get_thread_num() / block_size_x
    team = omp_get_team_num()
    bx = mod(team, nteams_x)
    by = team / nteams_x
    ii = tx + bx * block_size_x * tile_size_x
    jj = ty + by * block_size_y * tile_size_y

    do d = 1, 2
      do ti = 0, tile_size_x - 1
        sh_a(d, tx + ti * block_size_x + 1) = a(ii + ti * block_size_x + 1 + (d - 1) * m)
      end do
      if (tx == 0) then
        do tj = 0, tile_size_y - 1
          sh_b(d, ty + tj * block_size_y + 1) = b(jj + tj * block_size_y + 1 + (d - 1) * n)
        end do
      end if
    end do
    do ti = 0, tile_size_x - 1
      sh_scale_a(tx + ti * block_size_x + 1) = scale_a(ii + ti * block_size_x + 1)
    end do
    if (tx == 0) then
      do tj = 0, tile_size_y - 1
        sh_scale_b(ty + tj * block_size_y + 1) = scale_b(jj + tj * block_size_y + 1)
      end do
    end if

    s_cross_term = 0.0_real64
    do ti = 0, tile_size_x - 1
      do tj = 0, tile_size_y - 1
        if (ii + ti * block_size_x < m .and. jj + tj * block_size_y < n) then
          dist_ij = 0.0_real64
          do d = 1, 2
            dist_ij = dist_ij + &
              (sh_a(d, tx + ti * block_size_x + 1) - sh_b(d, ty + tj * block_size_y + 1)) * &
              (sh_a(d, tx + ti * block_size_x + 1) - sh_b(d, ty + tj * block_size_y + 1))
          end do
          s_cross_term = s_cross_term + &
            exp(-dist_ij / (sh_scale_a(tx + ti * block_size_x + 1) + sh_scale_b(ty + tj * block_size_y + 1)))
        end if
      end do
    end do

    !$omp atomic update
    sum = sum + s_cross_term
    !$omp barrier
    if (tx == 0 .and. ty == 0) cross_term(by * nteams_x + bx + 1) = sum
    !$omp end parallel
    !$omp end target teams
  end subroutine distance_real64

  function reduce_cross_term_real32(cross_term, m, n, nblocks) result(total)
    real(real32), intent(in) :: cross_term(:)
    integer, intent(in) :: m, n, nblocks
    real(real32) :: total
    integer :: i

    total = 0.0_real32
    !$omp target teams distribute parallel do reduction(+:total) thread_limit(reduce_block_size) map(tofrom: total)
    do i = 1, nblocks
      total = total + cross_term(i)
    end do
    !$omp end target teams distribute parallel do
  end function reduce_cross_term_real32

  function reduce_cross_term_real64(cross_term, m, n, nblocks) result(total)
    real(real64), intent(in) :: cross_term(:)
    integer, intent(in) :: m, n, nblocks
    real(real64) :: total
    integer :: i

    total = 0.0_real64
    !$omp target teams distribute parallel do reduction(+:total) thread_limit(reduce_block_size) map(tofrom: total)
    do i = 1, nblocks
      total = total + cross_term(i)
    end do
    !$omp end target teams distribute parallel do
  end function reduce_cross_term_real64

  function host_cost_real32(a, b, scale_a, scale_b, size) result(total32)
    real(real32), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real32) :: total32
    real(real64) :: total64
    integer :: i, j
    real(real32) :: dist

    total64 = 0.0_real64
    do i = 1, size
      do j = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total64 = total64 + real(exp(-dist / (scale_a(i) + scale_b(j))), real64)
      end do
    end do
    total32 = real(total64, real32)
  end function host_cost_real32

  function host_cost_real64(a, b, scale_a, scale_b, size) result(total)
    real(real64), intent(in) :: a(:), b(:), scale_a(:), scale_b(:)
    integer, intent(in) :: size
    real(real64) :: total
    integer :: i, j
    real(real64) :: dist

    total = 0.0_real64
    do i = 1, size
      do j = 1, size
        dist = (a(i) - b(j)) * (a(i) - b(j)) + &
               (a(i + size) - b(j + size)) * (a(i + size) - b(j + size))
        total = total + exp(-dist / (scale_a(i) + scale_b(j)))
      end do
    end do
  end function host_cost_real64

end program main
