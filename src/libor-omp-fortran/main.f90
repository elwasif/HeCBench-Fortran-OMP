! SPDX-License-Identifier: CC0-1.0
program libor_omp_fortran
  use iso_fortran_env, only: real32
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: sp = real32
  integer, parameter :: block_size = 64
  integer, parameter :: grid_size = 1500
  integer, parameter :: nn = 80
  integer, parameter :: nmat = 40
  integer, parameter :: l2_size = 3280
  integer, parameter :: nopt = 15
  integer, parameter :: npath = 96000

  integer :: repeat, arg_status, i
  integer :: maturities(nopt)
  real(sp) :: lambda(nn), swaprates(nopt), delta
  real(sp), allocatable :: h_v(:), h_lb(:)
  real(8) :: start_time, end_time, v_avg, lb_avg
  character(len=64) :: arg
  logical :: ok

  if (command_argument_count() /= 1) then
    call get_command_argument(0, arg)
    write(*,'("Usage: ",A," <repeat>")') trim(arg)
    stop 1
  end if

  call get_command_argument(1, arg, status=arg_status)
  read(arg, *, iostat=arg_status) repeat
  if (arg_status /= 0 .or. repeat <= 0) then
    call get_command_argument(0, arg)
    write(*,'("Usage: ",A," <repeat>")') trim(arg)
    stop 1
  end if

  delta = 0.25_sp
  lambda = 0.2_sp
  maturities = [4, 4, 4, 8, 8, 8, 20, 20, 20, 28, 28, 28, 40, 40, 40]
  swaprates = [0.045_sp, 0.05_sp, 0.055_sp, 0.045_sp, 0.05_sp, 0.055_sp, &
               0.045_sp, 0.05_sp, 0.055_sp, 0.045_sp, 0.05_sp, 0.055_sp, &
               0.045_sp, 0.05_sp, 0.055_sp]
  ok = .true.

  allocate(h_v(npath), h_lb(npath))

  !$omp target data map(to: maturities(1:nopt), swaprates(1:nopt), lambda(1:nn)) &
  !$omp& map(alloc: h_v(1:npath), h_lb(1:npath))
    start_time = omp_get_wtime()
    do i = 1, repeat
      call run_value_kernel(h_v, lambda, maturities, swaprates, delta)
    end do
    end_time = omp_get_wtime()
    write(*,'("Average kernel execution time : ",F8.6," (s)")') (end_time - start_time) / real(repeat, 8)

    !$omp target update from(h_v(1:npath))
    v_avg = sum(real(h_v, 8)) / real(npath, 8)
    if (abs(v_avg - 224.323d0) > 1.0d-3) then
      ok = .false.
      write(*,'("Expected: 224.323 Actual ",F15.3)') v_avg
    end if

    start_time = omp_get_wtime()
    do i = 1, repeat
      call run_greeks_kernel(h_v, h_lb, lambda, maturities, swaprates, delta)
    end do
    end_time = omp_get_wtime()
    write(*,'("Average kernel execution time : ",F8.6," (s)")') (end_time - start_time) / real(repeat, 8)

    !$omp target update from(h_lb(1:npath), h_v(1:npath))
  !$omp end target data

  v_avg = sum(real(h_v, 8)) / real(npath, 8)
  lb_avg = sum(real(h_lb, 8)) / real(npath, 8)

  if (abs(v_avg - 224.323d0) > 1.0d-3) then
    ok = .false.
    write(*,'("Expected: 224.323 Actual ",F15.3)') v_avg
  end if
  if (abs(lb_avg - 21.348d0) > 1.0d-3) then
    ok = .false.
    write(*,'("Expected:  21.348 Actual ",F15.3)') lb_avg
  end if

  deallocate(h_v, h_lb)
  if (.not. ok) stop 1

contains

  subroutine run_value_kernel(h_v, lambda, maturities, swaprates, delta)
    real(sp), intent(inout) :: h_v(npath)
    real(sp), intent(in) :: lambda(nn), swaprates(nopt), delta
    integer, intent(in) :: maturities(nopt)
    integer :: tid

    !$omp target teams distribute parallel do num_teams(grid_size) thread_limit(block_size)
    do tid = 1, grid_size * block_size
      block
        integer :: path, j, n, m, opt, thread_n
        real(sp) :: l(nn), z(nn), b(nmat), s(nmat)
        real(sp) :: sqez, lam, con1, accum, vrat, value, discount, swapval

        thread_n = grid_size * block_size
        do path = tid, npath, thread_n
          do j = 1, nn
            z(j) = 0.3_sp
            l(j) = 0.05_sp
          end do

          do n = 1, nmat
            sqez = sqrt(delta) * z(n)
            accum = 0.0_sp
            do j = n + 1, nn
              lam = lambda(j - n)
              con1 = delta * lam
              accum = accum + con1 * l(j) / (1.0_sp + delta * l(j))
              vrat = exp(con1 * accum + lam * (sqez - 0.5_sp * con1))
              l(j) = l(j) * vrat
            end do
          end do

          discount = 1.0_sp
          accum = 0.0_sp
          do n = nmat + 1, nn
            discount = discount / (1.0_sp + delta * l(n))
            accum = accum + delta * discount
            b(n - nmat) = discount
            s(n - nmat) = accum
          end do

          value = 0.0_sp
          do opt = 1, nopt
            m = maturities(opt)
            swapval = b(m) + swaprates(opt) * s(m) - 1.0_sp
            if (swapval < 0.0_sp) value = value - 100.0_sp * swapval
          end do

          discount = 1.0_sp
          do n = 1, nmat
            discount = discount / (1.0_sp + delta * l(n))
          end do
          h_v(path) = discount * value
        end do
      end block
    end do
    !$omp end target teams distribute parallel do
  end subroutine run_value_kernel

  subroutine run_greeks_kernel(h_v, h_lb, lambda, maturities, swaprates, delta)
    real(sp), intent(inout) :: h_v(npath), h_lb(npath)
    real(sp), intent(in) :: lambda(nn), swaprates(nopt), delta
    integer, intent(in) :: maturities(nopt)
    integer :: tid

    !$omp target teams distribute parallel do num_teams(grid_size) thread_limit(block_size)
    do tid = 1, grid_size * block_size
      block
        integer :: path, j, n, m, opt, thread_n
        real(sp) :: l(nn), l_b(nn), l2(l2_size), z(nn)
        real(sp) :: b(nmat), s(nmat), b_b(nmat), s_b(nmat)
        real(sp) :: sqez, lam, con1, accum, vrat, value, discount, swapval
        real(sp) :: faci, v1

        thread_n = grid_size * block_size
        do path = tid, npath, thread_n
          do j = 1, nn
            z(j) = 0.3_sp
            l(j) = 0.05_sp
            l2(j) = l(j)
            l_b(j) = 0.0_sp
          end do

          do n = 1, nmat
            sqez = sqrt(delta) * z(n)
            accum = 0.0_sp
            do j = n + 1, nn
              lam = lambda(j - n)
              con1 = delta * lam
              accum = accum + con1 * l(j) / (1.0_sp + delta * l(j))
              vrat = exp(con1 * accum + lam * (sqez - 0.5_sp * con1))
              l(j) = l(j) * vrat
              l2(j + n * nn) = l(j)
            end do
          end do

          discount = 1.0_sp
          accum = 0.0_sp
          do m = 1, nn - nmat
            n = m + nmat
            discount = discount / (1.0_sp + delta * l(n))
            accum = accum + delta * discount
            b(m) = discount
            s(m) = accum
          end do

          value = 0.0_sp
          b_b = 0.0_sp
          s_b = 0.0_sp
          do opt = 1, nopt
            m = maturities(opt)
            swapval = b(m) + swaprates(opt) * s(m) - 1.0_sp
            if (swapval < 0.0_sp) then
              value = value - 100.0_sp * swapval
              s_b(m) = s_b(m) - 100.0_sp * swaprates(opt)
              b_b(m) = b_b(m) - 100.0_sp
            end if
          end do

          do m = nn - nmat, 1, -1
            n = m + nmat
            b_b(m) = b_b(m) + delta * s_b(m)
            l_b(n) = -b_b(m) * b(m) * (delta / (1.0_sp + delta * l(n)))
            if (m > 1) then
              s_b(m - 1) = s_b(m - 1) + s_b(m)
              b_b(m - 1) = b_b(m - 1) + b_b(m) / (1.0_sp + delta * l(n))
            end if
          end do

          discount = 1.0_sp
          do n = 1, nmat
            discount = discount / (1.0_sp + delta * l(n))
          end do
          value = discount * value

          do n = 1, nmat
            l_b(n) = -value * delta / (1.0_sp + delta * l(n))
          end do
          do n = nmat + 1, nn
            l_b(n) = discount * l_b(n)
          end do

          do n = nmat, 1, -1
            v1 = 0.0_sp
            do j = nn, n + 1, -1
              v1 = v1 + lambda(j - n) * l2(j + n * nn) * l_b(j)
              faci = delta / (1.0_sp + delta * l2(j + (n - 1) * nn))
              l_b(j) = l_b(j) * (l2(j + n * nn) / l2(j + (n - 1) * nn)) + &
                       v1 * lambda(j - n) * faci * faci
            end do
          end do

          h_v(path) = value
          h_lb(path) = l_b(nn)
        end do
      end block
    end do
    !$omp end target teams distribute parallel do
  end subroutine run_greeks_kernel

end program libor_omp_fortran
