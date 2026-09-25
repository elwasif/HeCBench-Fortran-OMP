! SPDX-License-Identifier: CC0-1.0
module cooling_mod
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
contains
  real(real64) function primordial_cool(n, temp, heat_flag) result(cool)
    real(real64), intent(in) :: n, temp
    integer, intent(in) :: heat_flag
    real(real64) :: n_h, y_mass, y, g_ff
    real(real64) :: n_h0, n_hp, n_he0, n_hep, n_hepp, n_e, n_e_old
    real(real64) :: alpha_hp, alpha_hep, alpha_d, alpha_hepp
    real(real64) :: gamma_eh0, gamma_ehe0, gamma_ehep
    real(real64) :: le_h0, le_hep, li_h0, li_he0, li_hep
    real(real64) :: lr_hp, lr_hep, lr_hepp, ld_hep, l_ff
    real(real64) :: gamma_lh0, gamma_lhe0, gamma_lhep, e_h0, e_he0, e_hep, h_rate
    real(real64) :: diff, tol
    integer :: iter, n_iter

    y_mass = 0.24_real64
    y = y_mass / (4.0_real64 - 4.0_real64 * y_mass)
    n_h = n

    alpha_hp = 8.4e-11_real64 * (1.0_real64 / sqrt(temp)) * (temp / 1.0e3_real64)**(-0.2_real64) * &
        (1.0_real64 / (1.0_real64 + (temp / 1.0e6_real64)**0.7_real64))
    alpha_hep = 1.5e-10_real64 * temp**(-0.6353_real64)
    alpha_d = 1.9e-3_real64 * temp**(-1.5_real64) * exp(-470000.0_real64 / temp) * &
        (1.0_real64 + 0.3_real64 * exp(-94000.0_real64 / temp))
    alpha_hepp = 3.36e-10_real64 * (1.0_real64 / sqrt(temp)) * (temp / 1.0e3_real64)**(-0.2_real64) * &
        (1.0_real64 / (1.0_real64 + (temp / 1.0e6_real64)**0.7_real64))
    gamma_eh0 = 5.85e-11_real64 * sqrt(temp) * exp(-157809.1_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64)))
    gamma_ehe0 = 2.38e-11_real64 * sqrt(temp) * exp(-285335.4_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64)))
    gamma_ehep = 5.68e-12_real64 * sqrt(temp) * exp(-631515.0_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64)))
    gamma_lh0 = 3.19851e-13_real64
    gamma_lhe0 = 3.13029e-13_real64
    gamma_lhep = 2.00541e-14_real64
    e_h0 = 2.4796e-24_real64
    e_he0 = 6.86167e-24_real64
    e_hep = 6.21868e-25_real64

    n_e = n_h
    n_iter = 20
    tol = 1.0e-6_real64
    if (heat_flag /= 0) then
      do iter = 1, n_iter
        n_e_old = n_e
        n_h0 = n_h * alpha_hp / (alpha_hp + gamma_eh0 + gamma_lh0 / n_e)
        n_hp = n_h - n_h0
        n_hep = y * n_h / (1.0_real64 + (alpha_hep + alpha_d) / (gamma_ehe0 + gamma_lhe0 / n_e) + &
            (gamma_ehep + gamma_lhep / n_e) / alpha_hepp)
        n_he0 = n_hep * (alpha_hep + alpha_d) / (gamma_ehe0 + gamma_lhe0 / n_e)
        n_hepp = n_hep * (gamma_ehep + gamma_lhep / n_e) / alpha_hepp
        n_e = n_hp + n_hep + 2.0_real64 * n_hepp
        diff = abs(n_e_old - n_e)
        if (diff < tol) exit
      end do
    else
      n_h0 = n_h * alpha_hp / (alpha_hp + gamma_eh0)
      n_hp = n_h - n_h0
      n_hep = y * n_h / (1.0_real64 + (alpha_hep + alpha_d) / gamma_ehe0 + gamma_ehep / alpha_hepp)
      n_he0 = n_hep * (alpha_hep + alpha_d) / gamma_ehe0
      n_hepp = n_hep * gamma_ehep / alpha_hepp
      n_e = n_hp + n_hep + 2.0_real64 * n_hepp
    end if

    le_h0 = 7.50e-19_real64 * exp(-118348.0_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64))) * n_e * n_h0
    le_hep = 5.54e-17_real64 * temp**(-0.397_real64) * exp(-473638.0_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64))) * n_e * n_hep
    li_h0 = 1.27e-21_real64 * sqrt(temp) * exp(-157809.1_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64))) * n_e * n_h0
    li_he0 = 9.38e-22_real64 * sqrt(temp) * exp(-285335.4_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64))) * n_e * n_he0
    li_hep = 4.95e-22_real64 * sqrt(temp) * exp(-631515.0_real64 / temp) * &
        (1.0_real64 / (1.0_real64 + sqrt(temp / 1.0e5_real64))) * n_e * n_hep
    lr_hp = 8.70e-27_real64 * sqrt(temp) * (temp / 1.0e3_real64)**(-0.2_real64) * &
        (1.0_real64 / (1.0_real64 + (temp / 1.0e6_real64)**0.7_real64)) * n_e * n_hp
    lr_hep = 1.55e-26_real64 * temp**0.3647_real64 * n_e * n_hep
    lr_hepp = 3.48e-26_real64 * sqrt(temp) * (temp / 1.0e3_real64)**(-0.2_real64) * &
        (1.0_real64 / (1.0_real64 + (temp / 1.0e6_real64)**0.7_real64)) * n_e * n_hepp
    ld_hep = 1.24e-13_real64 * temp**(-1.5_real64) * exp(-470000.0_real64 / temp) * &
        (1.0_real64 + 0.3_real64 * exp(-94000.0_real64 / temp)) * n_e * n_hep
    g_ff = 1.1_real64 + 0.34_real64 * exp(-((5.5_real64 - log(temp)) * (5.5_real64 - log(temp))) / 3.0_real64)
    l_ff = 1.42e-27_real64 * g_ff * sqrt(temp) * (n_hp + n_hep + 4.0_real64 * n_hepp) * n_e

    cool = le_h0 + le_hep + li_h0 + li_he0 + li_hep + lr_hp + lr_hep + lr_hepp + ld_hep + l_ff
    h_rate = 0.0_real64
    if (heat_flag /= 0) h_rate = n_h0 * e_h0 + n_he0 * e_he0 + n_hep * e_hep
    cool = cool - h_rate
  end function primordial_cool
end module cooling_mod

program cooling_main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  use cooling_mod
  implicit none

  integer :: num, repeat, i
  real(real64), parameter :: density = 0.0899_real64
  real(real64), allocatable :: temp(:), host_result(:), device_result(:)
  real(real64) :: start_time, elapsed_ms
  character(len=64) :: arg
  logical :: ok

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <number of points> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) num
  call get_command_argument(2, arg)
  read(arg, *) repeat

  allocate(temp(0:num-1), host_result(0:num-1), device_result(0:num-1))
  do i = 0, num - 1
    temp(i) = -275.0_real64 + real(i, real64) * 275.0_real64 * 2.0_real64 / real(num, real64)
  end do

  !$omp target data map(to: temp) map(from: device_result)
  do i = 1, repeat
    call cool_kernel(num, density, temp, device_result, 1)
  end do

  start_time = omp_get_wtime()
  do i = 1, repeat
    call cool_kernel(num, density, temp, device_result, 1)
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
  !$omp end target data

  write(*,'("Average kernel execution time ",F0.6," (ms)")') elapsed_ms

  call reference_kernel(num, density, temp, host_result, 1)
  ok = check_values(device_result, host_result, num)
  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine cool_kernel(num, density, temp, result, heat_flag)
    integer, intent(in) :: num, heat_flag
    real(real64), intent(in) :: density, temp(0:)
    real(real64), intent(out) :: result(0:)
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256)
    do idx = 0, num - 1
      result(idx) = primordial_cool(density, temp(idx), heat_flag)
    end do
    !$omp end target teams distribute parallel do
  end subroutine cool_kernel

  subroutine reference_kernel(num, density, temp, result, heat_flag)
    integer, intent(in) :: num, heat_flag
    real(real64), intent(in) :: density, temp(0:)
    real(real64), intent(out) :: result(0:)
    integer :: idx

    do idx = 0, num - 1
      result(idx) = primordial_cool(density, temp(idx), heat_flag)
    end do
  end subroutine reference_kernel

  logical function check_values(device_values, reference_values, num)
    real(real64), intent(in) :: device_values(0:), reference_values(0:)
    integer, intent(in) :: num
    integer :: idx

    check_values = .true.
    do idx = 0, num - 1
      if (abs(device_values(idx) - reference_values(idx)) > 1.0e-3_real64) then
        check_values = .false.
        exit
      end if
    end do
  end function check_values

end program cooling_main
