! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  integer, parameter :: width = 256
  integer, parameter :: height = 256
  integer, parameter :: n = width * height
  real(real64), parameter :: omega = 1.2_real64
  real(real64), parameter :: epsilon = 1.0e-3_real64

  integer :: iterations
  integer, allocatable :: cell_type(:)
  real(real64), allocatable :: if0(:), of0(:), if1234(:), of1234(:), if5678(:), of5678(:)
  real(real64), allocatable :: ref0(:), ref1234(:), ref5678(:)
  real(real64) :: e(2, 9), w(9)
  logical :: ok

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

  if (command_argument_count() /= 1) then
    write(*,'(A)') 'Usage %s <iterations>'
    stop 1
  end if

  iterations = read_int_arg(1)

  allocate(cell_type(n))
  allocate(if0(n), of0(n), if1234(4*n), of1234(4*n), if5678(4*n), of5678(4*n))
  allocate(ref0(n), ref1234(4*n), ref5678(4*n))

  call init_constants(e, w)
  call init_problem(e, w, cell_type, if0, if1234, if5678)

  call reference_lbm(iterations, e, w, cell_type, if0, if1234, if5678, ref0, ref1234, ref5678)
  call fluid_sim(iterations, w, cell_type, if0, if1234, if5678, of0, of1234, of5678)

  ok = verify_results(of0, of1234, of5678, ref0, ref1234, ref5678)
  write(*,'(A)') merge('PASS', 'FAIL', ok)

  deallocate(cell_type)
  deallocate(if0, of0, if1234, of1234, if5678, of5678)
  deallocate(ref0, ref1234, ref5678)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  subroutine init_constants(e, w)
    real(real64), intent(out) :: e(2, 9), w(9)
    e(:, 1) = [0.0_real64, 0.0_real64]
    e(:, 2) = [1.0_real64, 0.0_real64]
    e(:, 3) = [0.0_real64, 1.0_real64]
    e(:, 4) = [-1.0_real64, 0.0_real64]
    e(:, 5) = [0.0_real64, -1.0_real64]
    e(:, 6) = [1.0_real64, 1.0_real64]
    e(:, 7) = [-1.0_real64, 1.0_real64]
    e(:, 8) = [-1.0_real64, -1.0_real64]
    e(:, 9) = [1.0_real64, -1.0_real64]
    w = [4.0_real64 / 9.0_real64, &
         1.0_real64 / 9.0_real64, 1.0_real64 / 9.0_real64, &
         1.0_real64 / 9.0_real64, 1.0_real64 / 9.0_real64, &
         1.0_real64 / 36.0_real64, 1.0_real64 / 36.0_real64, &
         1.0_real64 / 36.0_real64, 1.0_real64 / 36.0_real64]
  end subroutine init_constants

  subroutine init_problem(e, w, cell_type, if0, if1234, if5678)
    real(real64), intent(in) :: e(2, 9), w(9)
    integer, intent(out) :: cell_type(:)
    real(real64), intent(out) :: if0(:), if1234(:), if5678(:)
    integer :: x, y, pos, base, den
    real(real64) :: u0(2)

    call c_srand(123_c_int)
    u0 = [0.01_real64, 0.01_real64]
    do y = 1, height
      do x = 1, width
        pos = x + (y - 1) * width
        base = 4 * (pos - 1)
        den = mod(c_rand(), 10_c_int) + 1
        if0(pos) = compute_feq(real(den, real64), w(1), e(:, 1), u0)
        if1234(base + 1) = compute_feq(real(den, real64), w(2), e(:, 2), u0)
        if1234(base + 2) = compute_feq(real(den, real64), w(3), e(:, 3), u0)
        if1234(base + 3) = compute_feq(real(den, real64), w(4), e(:, 4), u0)
        if1234(base + 4) = compute_feq(real(den, real64), w(5), e(:, 5), u0)
        if5678(base + 1) = compute_feq(real(den, real64), w(6), e(:, 6), u0)
        if5678(base + 2) = compute_feq(real(den, real64), w(7), e(:, 7), u0)
        if5678(base + 3) = compute_feq(real(den, real64), w(8), e(:, 8), u0)
        if5678(base + 4) = compute_feq(real(den, real64), w(9), e(:, 9), u0)
        if (x == 1 .or. x == width .or. y == 1 .or. y == height) then
          cell_type(pos) = 1
        else
          cell_type(pos) = 0
        end if
      end do
    end do
  end subroutine init_problem

  real(real64) function compute_feq(rho, weight, dir, velocity) result(value)
    real(real64), intent(in) :: rho, weight, dir(2), velocity(2)
    real(real64) :: u2, eu
    u2 = velocity(1) * velocity(1) + velocity(2) * velocity(2)
    eu = dir(1) * velocity(1) + dir(2) * velocity(2)
    value = rho * weight * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
  end function compute_feq

  subroutine fluid_sim(iterations, w, cell_type, if0, if1234, if5678, of0, of1234, of5678)
    integer, intent(in) :: iterations, cell_type(:)
    real(real64), intent(in) :: w(:)
    real(real64), intent(inout) :: if0(:), if1234(:), if5678(:)
    real(real64), intent(out) :: of0(:), of1234(:), of5678(:)
    integer :: iter
    real(real64) :: start_time, end_time, elapsed

    of0 = if0
    of1234 = if1234
    of5678 = if5678

    !$omp target data map(to: cell_type(1:n), w(1:9), if0(1:n), of0(1:n), &
    !$omp& if1234(1:4*n), of1234(1:4*n), if5678(1:4*n), of5678(1:4*n))
    start_time = omp_get_wtime()
    do iter = 1, iterations
      if (mod(iter, 2) == 1) then
        call lbm_device(if0, of0, if1234, of1234, if5678, of5678, cell_type, w)
      else
        call lbm_device(of0, if0, of1234, if1234, of5678, if5678, cell_type, w)
      end if
    end do
    end_time = omp_get_wtime()

    elapsed = (end_time - start_time) / real(iterations, real64)
    write(*,'(A,F8.6,A)') 'Average kernel execution time ', elapsed, ' (s)'

    if (mod(iterations, 2) == 0) then
      !$omp target update from(if0(1:n), if1234(1:4*n), if5678(1:4*n))
    else
      !$omp target update from(of0(1:n), of1234(1:4*n), of5678(1:4*n))
    end if
    !$omp end target data

    if (mod(iterations, 2) == 0) then
      of0 = if0
      of1234 = if1234
      of5678 = if5678
    end if
  end subroutine fluid_sim

  subroutine lbm_device(in0, out0, in1234, out1234, in5678, out5678, cell_type, w)
    real(real64), intent(in) :: in0(:), in1234(:), in5678(:), w(:)
    real(real64), intent(inout) :: out0(:), out1234(:), out5678(:)
    integer, intent(in) :: cell_type(:)
    integer :: x, y, pos, base, dst, dst_base
    real(real64) :: f0, f1, f2, f3, f4, f5, f6, f7, f8
    real(real64) :: e0, e1, e2, e3, e4, e5, e6, e7, e8
    real(real64) :: rho, ux, uy, u2, eu

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) &
    !$omp& private(x, y, pos, base, dst, dst_base, f0, f1, f2, f3, f4, f5, f6, f7, f8, &
    !$omp& e0, e1, e2, e3, e4, e5, e6, e7, e8, rho, ux, uy, u2, eu)
    do y = 1, height
      do x = 1, width
        pos = x + (y - 1) * width
        base = 4 * (pos - 1)
        f0 = in0(pos)
        f1 = in1234(base + 1)
        f2 = in1234(base + 2)
        f3 = in1234(base + 3)
        f4 = in1234(base + 4)
        f5 = in5678(base + 1)
        f6 = in5678(base + 2)
        f7 = in5678(base + 3)
        f8 = in5678(base + 4)

        if (cell_type(pos) == 1) then
          e0 = f0
          e1 = f3
          e2 = f4
          e3 = f1
          e4 = f2
          e5 = f7
          e6 = f8
          e7 = f5
          e8 = f6
        else
          rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
          ux = (f1 - f3 + f5 - f6 - f7 + f8) / rho
          uy = (f2 - f4 + f5 + f6 - f7 - f8) / rho
          u2 = ux * ux + uy * uy

          eu = 0.0_real64
          e0 = rho * w(1) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = ux
          e1 = rho * w(2) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = uy
          e2 = rho * w(3) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = -ux
          e3 = rho * w(4) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = -uy
          e4 = rho * w(5) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = ux + uy
          e5 = rho * w(6) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = -ux + uy
          e6 = rho * w(7) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = -ux - uy
          e7 = rho * w(8) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)
          eu = ux - uy
          e8 = rho * w(9) * (1.0_real64 + 3.0_real64 * eu + 4.5_real64 * eu * eu - 1.5_real64 * u2)

          e0 = (1.0_real64 - omega) * f0 + omega * e0
          e1 = (1.0_real64 - omega) * f1 + omega * e1
          e2 = (1.0_real64 - omega) * f2 + omega * e2
          e3 = (1.0_real64 - omega) * f3 + omega * e3
          e4 = (1.0_real64 - omega) * f4 + omega * e4
          e5 = (1.0_real64 - omega) * f5 + omega * e5
          e6 = (1.0_real64 - omega) * f6 + omega * e6
          e7 = (1.0_real64 - omega) * f7 + omega * e7
          e8 = (1.0_real64 - omega) * f8 + omega * e8
        end if

        if (x > 1 .and. x < width .and. y > 1 .and. y < height) then
          out0(pos) = e0
          dst = (x + 1) + (y - 1) * width
          dst_base = 4 * (dst - 1)
          out1234(dst_base + 1) = e1
          dst = x + y * width
          dst_base = 4 * (dst - 1)
          out1234(dst_base + 2) = e2
          dst = (x - 1) + (y - 1) * width
          dst_base = 4 * (dst - 1)
          out1234(dst_base + 3) = e3
          dst = x + (y - 2) * width
          dst_base = 4 * (dst - 1)
          out1234(dst_base + 4) = e4
          dst = (x + 1) + y * width
          dst_base = 4 * (dst - 1)
          out5678(dst_base + 1) = e5
          dst = (x - 1) + y * width
          dst_base = 4 * (dst - 1)
          out5678(dst_base + 2) = e6
          dst = (x - 1) + (y - 2) * width
          dst_base = 4 * (dst - 1)
          out5678(dst_base + 3) = e7
          dst = (x + 1) + (y - 2) * width
          dst_base = 4 * (dst - 1)
          out5678(dst_base + 4) = e8
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine lbm_device

  subroutine reference_lbm(iterations, e, w, cell_type, if0_src, if1234_src, if5678_src, out0, out1234, out5678)
    integer, intent(in) :: iterations, cell_type(:)
    real(real64), intent(in) :: e(2, 9), w(9), if0_src(:), if1234_src(:), if5678_src(:)
    real(real64), intent(out) :: out0(:), out1234(:), out5678(:)
    real(real64), allocatable :: in0(:), in1234(:), in5678(:), ef0(:), ef1234(:), ef5678(:)
    integer :: iter, x, y, pos, base, k, dst, dst_base
    real(real64) :: f0, rho, ux, uy, vel(2)

    allocate(in0(n), in1234(4*n), in5678(4*n), ef0(n), ef1234(4*n), ef5678(4*n))
    in0 = if0_src
    in1234 = if1234_src
    in5678 = if5678_src
    out0 = if0_src
    out1234 = if1234_src
    out5678 = if5678_src

    do iter = 1, iterations
      do y = 1, height
        do x = 1, width
          pos = x + (y - 1) * width
          base = 4 * (pos - 1)
          if (cell_type(pos) == 1) then
            ef0(pos) = in0(pos)
            ef1234(base + 1) = in1234(base + 3)
            ef1234(base + 2) = in1234(base + 4)
            ef1234(base + 3) = in1234(base + 1)
            ef1234(base + 4) = in1234(base + 2)
            ef5678(base + 1) = in5678(base + 3)
            ef5678(base + 2) = in5678(base + 4)
            ef5678(base + 3) = in5678(base + 1)
            ef5678(base + 4) = in5678(base + 2)
          else
            f0 = in0(pos)
            rho = f0 + sum(in1234(base + 1:base + 4)) + sum(in5678(base + 1:base + 4))
            ux = in1234(base + 1) * e(1, 2) + in1234(base + 2) * e(1, 3) + &
                 in1234(base + 3) * e(1, 4) + in1234(base + 4) * e(1, 5) + &
                 in5678(base + 1) * e(1, 6) + in5678(base + 2) * e(1, 7) + &
                 in5678(base + 3) * e(1, 8) + in5678(base + 4) * e(1, 9)
            uy = in1234(base + 1) * e(2, 2) + in1234(base + 2) * e(2, 3) + &
                 in1234(base + 3) * e(2, 4) + in1234(base + 4) * e(2, 5) + &
                 in5678(base + 1) * e(2, 6) + in5678(base + 2) * e(2, 7) + &
                 in5678(base + 3) * e(2, 8) + in5678(base + 4) * e(2, 9)
            vel = [ux / rho, uy / rho]
            ef0(pos) = (1.0_real64 - omega) * in0(pos) + omega * compute_feq(rho, w(1), e(:, 1), vel)
            do k = 1, 4
              ef1234(base + k) = (1.0_real64 - omega) * in1234(base + k) + &
                  omega * compute_feq(rho, w(k + 1), e(:, k + 1), vel)
              ef5678(base + k) = (1.0_real64 - omega) * in5678(base + k) + &
                  omega * compute_feq(rho, w(k + 5), e(:, k + 5), vel)
            end do
          end if
        end do
      end do

      do y = 2, height - 1
        do x = 2, width - 1
          pos = x + (y - 1) * width
          base = 4 * (pos - 1)
          do k = 1, 9
            dst = (x + int(e(1, k))) + (y + int(e(2, k)) - 1) * width
            dst_base = 4 * (dst - 1)
            select case (k)
            case (1)
              out0(dst) = ef0(pos)
            case (2:5)
              out1234(dst_base + k - 1) = ef1234(base + k - 1)
            case (6:9)
              out5678(dst_base + k - 5) = ef5678(base + k - 5)
            end select
          end do
        end do
      end do

      in0 = out0
      in1234 = out1234
      in5678 = out5678
    end do

    deallocate(in0, in1234, in5678, ef0, ef1234, ef5678)
  end subroutine reference_lbm

  logical function verify_results(out0, out1234, out5678, ref0, ref1234, ref5678) result(ok)
    real(real64), intent(in) :: out0(:), out1234(:), out5678(:), ref0(:), ref1234(:), ref5678(:)
    integer :: idx
    ok = .true.
    do idx = 1, n
      if (out0(idx) - ref0(idx) > epsilon) ok = .false.
    end do
    do idx = 1, 4 * n
      if (out1234(idx) - ref1234(idx) > epsilon) ok = .false.
      if (out5678(idx) - ref5678(idx) > epsilon) ok = .false.
    end do
  end function verify_results

end program main
