! SPDX-License-Identifier: CC0-1.0
program burger
  use iso_fortran_env, only: real64
  use omp_lib
  implicit none

  integer :: argc, x_points, y_points, num_itrs
  character(len=64) :: arg
  real(real64), allocatable :: x(:), y(:), u(:), v(:), u_new(:), v_new(:), d_u(:), d_v(:)
  real(real64) :: x_len, y_len, del_x, del_y, nu, sigma, del_t
  integer :: i, j, itr, n
  real(real64) :: start_time, end_time
  character(len=32) :: time_text
  logical :: ok

  argc = command_argument_count()
  if (argc /= 3) then
    call get_command_argument(0, arg)
    print '(A,A,A)', 'Usage: ', trim(arg), ' <dim_x> <dim_y> <nt>'
    print '(A)', 'dim_x: number of grid points in the x axis'
    print '(A)', 'dim_y: number of grid points in the y axis'
    print '(A)', 'nt: number of time steps'
    stop 255
  end if

  call get_command_argument(1, arg)
  read(arg, *) x_points
  call get_command_argument(2, arg)
  read(arg, *) y_points
  call get_command_argument(3, arg)
  read(arg, *) num_itrs

  x_len = 2.0_real64
  y_len = 2.0_real64
  del_x = x_len / real(x_points - 1, real64)
  del_y = y_len / real(y_points - 1, real64)
  nu = 0.01_real64
  sigma = 0.0009_real64
  del_t = sigma * del_x * del_y / nu
  n = x_points * y_points

  allocate(x(x_points), y(y_points), u(n), v(n), u_new(n), v_new(n), d_u(n), d_v(n))

  print '(A)', "2D Burger's equation"
  print '(A,I0,A,I0)', 'Grid dimension: x = ', x_points, ' y = ', y_points
  print '(A,I0)', 'Number of time steps: ', num_itrs

  do i = 1, x_points
    x(i) = real(i - 1, real64) * del_x
  end do
  do i = 1, y_points
    y(i) = real(i - 1, real64) * del_y
  end do
  call initialize_fields(x_points, y_points, x, y, u, v, u_new, v_new)

  !$omp target data map(to: u_new(1:n), v_new(1:n)) map(tofrom: u(1:n), v(1:n))
    start_time = omp_get_wtime()

    do itr = 1, num_itrs
      call device_step(x_points, y_points, del_x, del_y, del_t, nu, u, v, u_new, v_new)
      call device_boundaries(x_points, y_points, u_new, v_new)
      call device_copy(x_points, y_points, u, v, u_new, v_new)
    end do

    end_time = omp_get_wtime()
    write(time_text, '(F0.6)') end_time - start_time
    if (time_text(1:1) == '.') time_text = '0' // trim(time_text)
    print '(A,A,A)', 'Total kernel execution time ', trim(time_text), ' (s)'
  !$omp end target data

  d_u = u
  d_v = v

  print '(A)', 'Serial computing for verification...'
  call initialize_fields(x_points, y_points, x, y, u, v, u_new, v_new)
  do itr = 1, num_itrs
    call host_step(x_points, y_points, del_x, del_y, del_t, nu, u, v, u_new, v_new)
    call host_boundaries(x_points, y_points, u_new, v_new)
    u = u_new
    v = v_new
  end do

  ok = .true.
  do i = 1, n
    if (abs(d_u(i) - u(i)) > 1.0e-6_real64 .or. abs(d_v(i) - v(i)) > 1.0e-6_real64) ok = .false.
  end do
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(x, y, u, v, u_new, v_new, d_u, d_v)

contains

  integer function idx(row, col, x_points) result(pos)
    integer, intent(in) :: row, col, x_points
    pos = row * x_points + col + 1
  end function idx

  subroutine initialize_fields(x_points, y_points, x, y, u, v, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(out) :: u(:), v(:), u_new(:), v_new(:)
    integer :: i, j, p

    do i = 0, y_points - 1
      do j = 0, x_points - 1
        p = idx(i, j, x_points)
        u(p) = 1.0_real64
        v(p) = 1.0_real64
        u_new(p) = 1.0_real64
        v_new(p) = 1.0_real64
        if (x(j + 1) > 0.5_real64 .and. x(j + 1) < 1.0_real64 .and. &
            y(i + 1) > 0.5_real64 .and. y(i + 1) < 1.0_real64) then
          u(p) = 2.0_real64
          v(p) = 2.0_real64
          u_new(p) = 2.0_real64
          v_new(p) = 2.0_real64
        end if
      end do
    end do
  end subroutine initialize_fields

  subroutine device_step(x_points, y_points, del_x, del_y, del_t, nu, u, v, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(in) :: del_x, del_y, del_t, nu
    real(real64), intent(inout) :: u(:), v(:), u_new(:), v_new(:)
    integer :: i, j, p

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(i,j,p) nowait
    do i = 1, y_points - 2
      do j = 1, x_points - 2
        p = idx(i, j, x_points)
        u_new(p) = u(p) + (nu * del_t / (del_x * del_x)) * &
            (u(idx(i, j + 1, x_points)) + u(idx(i, j - 1, x_points)) - 2.0_real64 * u(p)) + &
            (nu * del_t / (del_y * del_y)) * &
            (u(idx(i + 1, j, x_points)) + u(idx(i - 1, j, x_points)) - 2.0_real64 * u(p)) - &
            (del_t / del_x) * u(p) * (u(p) - u(idx(i, j - 1, x_points))) - &
            (del_t / del_y) * v(p) * (u(p) - u(idx(i - 1, j, x_points)))

        v_new(p) = v(p) + (nu * del_t / (del_x * del_x)) * &
            (v(idx(i, j + 1, x_points)) + v(idx(i, j - 1, x_points)) - 2.0_real64 * v(p)) + &
            (nu * del_t / (del_y * del_y)) * &
            (v(idx(i + 1, j, x_points)) + v(idx(i - 1, j, x_points)) - 2.0_real64 * v(p)) - &
            (del_t / del_x) * u(p) * (v(p) - v(idx(i, j - 1, x_points))) - &
            (del_t / del_y) * v(p) * (v(p) - v(idx(i - 1, j, x_points)))
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine device_step

  subroutine device_boundaries(x_points, y_points, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(inout) :: u_new(:), v_new(:)
    integer :: i, j

    !$omp target teams distribute parallel do thread_limit(256) private(i) nowait
    do i = 0, x_points - 1
      u_new(idx(0, i, x_points)) = 1.0_real64
      v_new(idx(0, i, x_points)) = 1.0_real64
      u_new(idx(y_points - 1, i, x_points)) = 1.0_real64
      v_new(idx(y_points - 1, i, x_points)) = 1.0_real64
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do thread_limit(256) private(j) nowait
    do j = 0, y_points - 1
      u_new(idx(j, 0, x_points)) = 1.0_real64
      v_new(idx(j, 0, x_points)) = 1.0_real64
      u_new(idx(j, x_points - 1, x_points)) = 1.0_real64
      v_new(idx(j, x_points - 1, x_points)) = 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine device_boundaries

  subroutine device_copy(x_points, y_points, u, v, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(inout) :: u(:), v(:), u_new(:), v_new(:)
    integer :: i, j, p

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(i,j,p)
    do i = 0, y_points - 1
      do j = 0, x_points - 1
        p = idx(i, j, x_points)
        u(p) = u_new(p)
        v(p) = v_new(p)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine device_copy

  subroutine host_step(x_points, y_points, del_x, del_y, del_t, nu, u, v, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(in) :: del_x, del_y, del_t, nu
    real(real64), intent(inout) :: u(:), v(:), u_new(:), v_new(:)
    integer :: i, j, p

    do i = 1, y_points - 2
      do j = 1, x_points - 2
        p = idx(i, j, x_points)
        u_new(p) = u(p) + (nu * del_t / (del_x * del_x)) * &
            (u(idx(i, j + 1, x_points)) + u(idx(i, j - 1, x_points)) - 2.0_real64 * u(p)) + &
            (nu * del_t / (del_y * del_y)) * &
            (u(idx(i + 1, j, x_points)) + u(idx(i - 1, j, x_points)) - 2.0_real64 * u(p)) - &
            (del_t / del_x) * u(p) * (u(p) - u(idx(i, j - 1, x_points))) - &
            (del_t / del_y) * v(p) * (u(p) - u(idx(i - 1, j, x_points)))

        v_new(p) = v(p) + (nu * del_t / (del_x * del_x)) * &
            (v(idx(i, j + 1, x_points)) + v(idx(i, j - 1, x_points)) - 2.0_real64 * v(p)) + &
            (nu * del_t / (del_y * del_y)) * &
            (v(idx(i + 1, j, x_points)) + v(idx(i - 1, j, x_points)) - 2.0_real64 * v(p)) - &
            (del_t / del_x) * u(p) * (v(p) - v(idx(i, j - 1, x_points))) - &
            (del_t / del_y) * v(p) * (v(p) - v(idx(i - 1, j, x_points)))
      end do
    end do
  end subroutine host_step

  subroutine host_boundaries(x_points, y_points, u_new, v_new)
    integer, intent(in) :: x_points, y_points
    real(real64), intent(inout) :: u_new(:), v_new(:)
    integer :: i, j

    do i = 0, x_points - 1
      u_new(idx(0, i, x_points)) = 1.0_real64
      v_new(idx(0, i, x_points)) = 1.0_real64
      u_new(idx(y_points - 1, i, x_points)) = 1.0_real64
      v_new(idx(y_points - 1, i, x_points)) = 1.0_real64
    end do
    do j = 0, y_points - 1
      u_new(idx(j, 0, x_points)) = 1.0_real64
      v_new(idx(j, 0, x_points)) = 1.0_real64
      u_new(idx(j, x_points - 1, x_points)) = 1.0_real64
      v_new(idx(j, x_points - 1, x_points)) = 1.0_real64
    end do
  end subroutine host_boundaries

end program burger
