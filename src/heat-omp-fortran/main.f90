! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer :: n, nsteps, t
  real(real64) :: start_total, stop_total, alpha, length, dx, dt, r, r2
  real(real64) :: tic, toc, norm, bandwidth, pi
  real(real64), allocatable :: u(:), u_tmp(:)

  start_total = omp_get_wtime()
  n = 1000
  nsteps = 10
  if (command_argument_count() == 2) then
    n = read_arg(1)
    nsteps = read_arg(2)
  end if
  if (n < 0 .or. nsteps < 0) stop 1

  alpha = 0.1_real64
  length = 1000.0_real64
  dx = length / real(n + 1, real64)
  dt = 0.5_real64 / real(nsteps, real64)
  r = alpha * dt / (dx * dx)
  r2 = 1.0_real64 - 4.0_real64 * r
  pi = acos(-1.0_real64)

  print '(A)'
  print '(A)', ' MMS heat equation'
  print '(A)'
  print '(A)', '--------------------'
  print '(A)', 'Problem input'
  print '(A)'
  print '(A,I0,A,I0)', ' Grid size: ', n, ' x ', n
  print '(A,ES12.6E2)', ' Cell width: ', dx
  print '(A,F11.6,A,F11.6)', ' Grid length: ', length, ' x ', length
  print '(A)'
  print '(A,ES12.6E2)', ' Alpha: ', alpha
  print '(A)'
  print '(A,I0)', ' Steps: ', nsteps
  print '(A,ES12.6E2)', ' Total time: ', dt * real(nsteps, real64)
  print '(A,ES12.6E2)', ' Time step: ', dt
  print '(A)', '--------------------'
  print '(A)', 'Stability'
  print '(A)'
  print '(A,F8.6)', ' r value: ', r
  if (r > 0.5_real64) print '(A)', ' Warning: unstable'
  print '(A)', '--------------------'

  allocate(u(n * n), u_tmp(n * n))

  !$omp target data map(tofrom: u(1:n*n), u_tmp(1:n*n))
  call initialize_grid(u, u_tmp, n, dx, length, pi)
  tic = omp_get_wtime()
  do t = 1, nsteps
    if (mod(t, 2) == 1) then
      call solve_step(u, u_tmp, n, r, r2)
    else
      call solve_step(u_tmp, u, n, r, r2)
    end if
  end do
  toc = omp_get_wtime()
  !$omp end target data

  if (mod(nsteps, 2) == 0) then
    norm = l2norm(n, u, nsteps, dt, alpha, dx, length, pi)
  else
    norm = l2norm(n, u_tmp, nsteps, dt, alpha, dx, length, pi)
  end if

  stop_total = omp_get_wtime()
  bandwidth = 1.0e-9_real64 * 2.0_real64 * real(n, real64) * real(n, real64) * real(nsteps, real64) * 8.0_real64 / (toc - tic)

  print '(A)', 'Results'
  print '(A)'
  print '(A,ES12.6E2)', 'Error (L2norm): ', norm
  print '(A,F8.6)', 'Solve time (s): ', toc - tic
  print '(A,F8.6)', 'Total time (s): ', stop_total - start_total
  print '(A,F8.6)', 'Bandwidth (GB/s): ', bandwidth
  print '(A)', '--------------------'

  deallocate(u, u_tmp)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_grid(u, u_tmp, n, dx, length, pi)
    real(real64), intent(out) :: u(:), u_tmp(:)
    integer, intent(in) :: n
    real(real64), intent(in) :: dx, length, pi
    integer :: i, j
    real(real64) :: x, y
    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(x, y)
    do j = 0, n - 1
      do i = 0, n - 1
        y = real(j + 1, real64) * dx
        x = real(i + 1, real64) * dx
        u(i + j * n + 1) = sin(pi * x / length) * sin(pi * y / length)
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do collapse(2) thread_limit(256)
    do j = 0, n - 1
      do i = 0, n - 1
        u_tmp(i + j * n + 1) = 0.0_real64
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine initialize_grid

  subroutine solve_step(u, u_tmp, n, r, r2)
    real(real64), intent(in) :: u(:)
    real(real64), intent(out) :: u_tmp(:)
    integer, intent(in) :: n
    real(real64), intent(in) :: r, r2
    integer :: i, j, idx
    real(real64) :: east, west, north, south
    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(idx, east, west, north, south)
    do j = 0, n - 1
      do i = 0, n - 1
        idx = i + j * n + 1
        if (i < n - 1) then
          east = u(idx + 1)
        else
          east = 0.0_real64
        end if
        if (i > 0) then
          west = u(idx - 1)
        else
          west = 0.0_real64
        end if
        if (j < n - 1) then
          north = u(idx + n)
        else
          north = 0.0_real64
        end if
        if (j > 0) then
          south = u(idx - n)
        else
          south = 0.0_real64
        end if
        u_tmp(idx) = r2 * u(idx) + r * (east + west + north + south)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine solve_step

  real(real64) function solution(t, x, y, alpha, length, pi)
    real(real64), intent(in) :: t, x, y, alpha, length, pi
    solution = exp(-2.0_real64 * alpha * pi * pi * t / (length * length)) * sin(pi * x / length) * sin(pi * y / length)
  end function solution

  real(real64) function l2norm(n, u, nsteps, dt, alpha, dx, length, pi)
    integer, intent(in) :: n, nsteps
    real(real64), intent(in) :: u(:), dt, alpha, dx, length, pi
    integer :: i, j
    real(real64) :: time, x, y, answer
    time = dt * real(nsteps, real64)
    l2norm = 0.0_real64
    y = dx
    do j = 0, n - 1
      x = dx
      do i = 0, n - 1
        answer = solution(time, x, y, alpha, length, pi)
        l2norm = l2norm + (u(i + j * n + 1) - answer) * (u(i + j * n + 1) - answer)
        x = x + dx
      end do
      y = y + dx
    end do
    l2norm = sqrt(l2norm)
  end function l2norm

end program main
