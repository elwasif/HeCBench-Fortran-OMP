! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int8, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name='rand')
      import :: c_int
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: neurons, iterations
  real(real32), allocatable :: ge(:), gi(:), h(:), m(:), n(:), v(:), lastspike(:)
  real(real32), allocatable :: ref_ge(:), ref_gi(:), ref_h(:), ref_m(:), ref_n(:), ref_v(:), ref_lastspike(:)
  integer(int8), allocatable :: not_refract(:), ref_not_refract(:)
  real(real32) :: dt, current_t
  real(real64) :: rsme

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <neurons> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) neurons
  call get_command_argument(2, arg); read(arg, *) iterations
  if (neurons <= 0 .or. iterations <= 0) stop 1

  allocate(ge(neurons), gi(neurons), h(neurons), m(neurons), n(neurons), v(neurons), lastspike(neurons))
  allocate(ref_ge(neurons), ref_gi(neurons), ref_h(neurons), ref_m(neurons), ref_n(neurons), ref_v(neurons))
  allocate(ref_lastspike(neurons), not_refract(neurons), ref_not_refract(neurons))

  write(*,'(A)', advance='no') 'initializing ... '
  call c_srand(2_c_int)
  call initialize_inputs(ge, gi, h, m, n, v, lastspike, not_refract)
  ref_ge = ge
  ref_gi = gi
  ref_h = h
  ref_m = m
  ref_n = n
  ref_v = v
  ref_lastspike = lastspike
  ref_not_refract = not_refract
  dt = 0.0001_real32
  current_t = 0.01_real32
  write(*,'(A)') 'done.'

  call neurongroup_stateupdater_host(ref_ge, ref_gi, ref_h, ref_m, ref_n, ref_v, ref_lastspike, &
                                     dt, current_t, ref_not_refract, neurons, iterations)
  call neurongroup_stateupdater(ge, gi, h, m, n, v, lastspike, dt, current_t, not_refract, neurons, iterations)

  rsme = compute_rsme(ge, gi, h, m, n, v, not_refract, ref_ge, ref_gi, ref_h, ref_m, ref_n, ref_v, &
                      ref_not_refract, neurons)
  write(*,'(A,F8.6)') 'RSME = ', rsme

  deallocate(ge, gi, h, m, n, v, lastspike, ref_ge, ref_gi, ref_h, ref_m, ref_n, ref_v, ref_lastspike)
  deallocate(not_refract, ref_not_refract)

contains

  subroutine initialize_inputs(ge, gi, h, m, n, v, lastspike, not_refract)
    real(real32), intent(out) :: ge(:), gi(:), h(:), m(:), n(:), v(:), lastspike(:)
    integer(int8), intent(out) :: not_refract(:)
    integer :: i

    do i = 2, size(ge)
      ge(i) = 0.15_real32 + merge(0.1_real32, -0.1_real32, next_rand_mod(2) == 0)
      gi(i) = 0.25_real32 + merge(0.2_real32, -0.2_real32, next_rand_mod(2) == 0)
      h(i) = 0.35_real32 + merge(0.3_real32, -0.3_real32, next_rand_mod(2) == 0)
      m(i) = 0.45_real32 + merge(0.4_real32, -0.4_real32, next_rand_mod(2) == 0)
      n(i) = 0.55_real32 + merge(0.5_real32, -0.5_real32, next_rand_mod(2) == 0)
      v(i) = 0.65_real32 + merge(0.6_real32, -0.6_real32, next_rand_mod(2) == 0)
      lastspike(i) = 1.0_real32 / real(next_rand_mod(1000) + 1, real32)
      not_refract(i) = 0_int8
    end do
  end subroutine initialize_inputs

  integer function next_rand_mod(divisor)
    integer, intent(in) :: divisor
    integer(c_int) :: value

    value = c_rand()
    next_rand_mod = modulo(value, divisor)
  end function next_rand_mod

  integer function timestep(time_value, step)
    real(real32), intent(in) :: time_value, step

    timestep = int((time_value + 1.0e-3_real32 * step) / step)
  end function timestep

  subroutine constants(dt, lio_1, lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, &
                       lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, &
                       lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28, &
                       lio_29, lio_30, lio_31, lio_32, lio_33)
    real(real32), intent(in) :: dt
    integer, intent(out) :: lio_1
    real(real32), intent(out) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32), intent(out) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32), intent(out) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32), intent(out) :: lio_29, lio_30, lio_31, lio_32, lio_33

    lio_1 = timestep(0.003_real32, dt)
    lio_2 = 9.939082_real32
    lio_3 = -55.555556_real32
    lio_4 = 0.00001_real32
    lio_5 = -200.0_real32
    lio_6 = -0.02016_real32
    lio_7 = -0.000001_real32
    lio_8 = 0.0_real32
    lio_9 = -250.0_real32
    lio_10 = 0.00416_real32
    lio_11 = 0.02016_real32
    lio_12 = -0.01764_real32
    lio_13 = -0.000001_real32
    lio_14 = 0.000099_real32
    lio_15 = 200.0_real32
    lio_16 = 0.0112_real32
    lio_17 = -0.002016_real32
    lio_18 = 0.0_real32
    lio_19 = 0.00048_real32
    lio_20 = 0.002016_real32
    lio_21 = 132.901474_real32
    lio_22 = -25.0_real32
    lio_23 = exp(-2000.0_real32 * dt)
    lio_24 = -3.0_real32
    lio_25 = -2700.0_real32
    lio_26 = 5000.0_real32
    lio_27 = 0.0_real32
    lio_28 = -400000000.0_real32
    lio_29 = -50.0_real32
    lio_30 = -30000.0_real32
    lio_31 = 100000.0_real32
    lio_32 = 5000000000.0_real32
    lio_33 = exp(-100.0_real32 * dt)
  end subroutine constants

  subroutine neurongroup_stateupdater(ge, gi, h, m, n, v, lastspike, dt, current_t, not_refract, neurons, iterations)
    real(real32), intent(inout) :: ge(:), gi(:), h(:), m(:), n(:), v(:)
    real(real32), intent(in) :: lastspike(:), dt, current_t
    integer(int8), intent(out) :: not_refract(:)
    integer, intent(in) :: neurons, iterations
    integer :: iter, lio_1
    real(real32) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32) :: lio_29, lio_30, lio_31, lio_32, lio_33
    real(real64) :: start_time, elapsed_us

    call constants(dt, lio_1, lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, &
                   lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, &
                   lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28, &
                   lio_29, lio_30, lio_31, lio_32, lio_33)

    !$omp target data map(tofrom: h(1:neurons), m(1:neurons), n(1:neurons), ge(1:neurons), v(1:neurons), gi(1:neurons)) &
    !$omp& map(to: lastspike(1:neurons)) map(from: not_refract(1:neurons))
    start_time = omp_get_wtime()
    do iter = 1, iterations
      call cobahh_step(h, m, n, ge, v, gi, lastspike, not_refract, neurons, dt, current_t, lio_1, lio_2, lio_3, &
                       lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, lio_13, lio_14, &
                       lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, lio_24, &
                       lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(iterations, real64)
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed_us, ' (us)'
    !$omp end target data
  end subroutine neurongroup_stateupdater

  subroutine neurongroup_stateupdater_host(ge, gi, h, m, n, v, lastspike, dt, current_t, not_refract, neurons, iterations)
    real(real32), intent(inout) :: ge(:), gi(:), h(:), m(:), n(:), v(:)
    real(real32), intent(in) :: lastspike(:), dt, current_t
    integer(int8), intent(out) :: not_refract(:)
    integer, intent(in) :: neurons, iterations
    integer :: iter, lio_1
    real(real32) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32) :: lio_29, lio_30, lio_31, lio_32, lio_33

    call constants(dt, lio_1, lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, &
                   lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, &
                   lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28, &
                   lio_29, lio_30, lio_31, lio_32, lio_33)
    do iter = 1, iterations
      call cobahh_step_host(h, m, n, ge, v, gi, lastspike, not_refract, neurons, dt, current_t, lio_1, lio_2, lio_3, &
                            lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, lio_13, &
                            lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, &
                            lio_24, lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    end do
  end subroutine neurongroup_stateupdater_host

  subroutine cobahh_step(h, m, n, ge, v, gi, lastspike, not_refract, neurons, dt, current_t, lio_1, lio_2, lio_3, &
                         lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, lio_13, lio_14, &
                         lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, lio_24, &
                         lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    real(real32), intent(inout) :: h(:), m(:), n(:), ge(:), v(:), gi(:)
    real(real32), intent(in) :: lastspike(:), dt, current_t
    integer(int8), intent(out) :: not_refract(:)
    integer, intent(in) :: neurons, lio_1
    real(real32), intent(in) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32), intent(in) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32), intent(in) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32), intent(in) :: lio_29, lio_30, lio_31, lio_32, lio_33
    integer :: idx
    real(real32) :: hh, mm, nn, gee, vv, gii, last, ba_h, new_h, ba_m, new_m, ba_n, new_n, new_ge
    real(real32) :: ba_v, new_v, new_gi

    !$omp target teams distribute parallel do thread_limit(256) private(hh, mm, nn, gee, vv, gii, last, ba_h, new_h, ba_m, &
    !$omp& new_m, ba_n, new_n, new_ge, ba_v, new_v, new_gi)
    do idx = 1, neurons
      call update_one(h(idx), m(idx), n(idx), ge(idx), v(idx), gi(idx), lastspike(idx), not_refract(idx), dt, current_t, &
                      lio_1, lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, &
                      lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, &
                      lio_24, lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    end do
    !$omp end target teams distribute parallel do
  end subroutine cobahh_step

  subroutine cobahh_step_host(h, m, n, ge, v, gi, lastspike, not_refract, neurons, dt, current_t, lio_1, lio_2, lio_3, &
                              lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, lio_13, lio_14, &
                              lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, lio_24, &
                              lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    real(real32), intent(inout) :: h(:), m(:), n(:), ge(:), v(:), gi(:)
    real(real32), intent(in) :: lastspike(:), dt, current_t
    integer(int8), intent(out) :: not_refract(:)
    integer, intent(in) :: neurons, lio_1
    real(real32), intent(in) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32), intent(in) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32), intent(in) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32), intent(in) :: lio_29, lio_30, lio_31, lio_32, lio_33
    integer :: idx

    do idx = 1, neurons
      call update_one(h(idx), m(idx), n(idx), ge(idx), v(idx), gi(idx), lastspike(idx), not_refract(idx), dt, current_t, &
                      lio_1, lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, &
                      lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, &
                      lio_24, lio_25, lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    end do
  end subroutine cobahh_step_host

  subroutine update_one(h, m, n, ge, v, gi, lastspike, not_refract, dt, current_t, lio_1, lio_2, lio_3, lio_4, &
                        lio_5, lio_6, lio_7, lio_8, lio_9, lio_10, lio_11, lio_12, lio_13, lio_14, lio_15, &
                        lio_16, lio_17, lio_18, lio_19, lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, &
                        lio_26, lio_27, lio_28, lio_29, lio_30, lio_31, lio_32, lio_33)
    real(real32), intent(inout) :: h, m, n, ge, v, gi
    real(real32), intent(in) :: lastspike, dt, current_t
    integer(int8), intent(out) :: not_refract
    integer, intent(in) :: lio_1
    real(real32), intent(in) :: lio_2, lio_3, lio_4, lio_5, lio_6, lio_7, lio_8, lio_9, lio_10
    real(real32), intent(in) :: lio_11, lio_12, lio_13, lio_14, lio_15, lio_16, lio_17, lio_18, lio_19
    real(real32), intent(in) :: lio_20, lio_21, lio_22, lio_23, lio_24, lio_25, lio_26, lio_27, lio_28
    real(real32), intent(in) :: lio_29, lio_30, lio_31, lio_32, lio_33
    real(real32) :: ba_h, new_h, ba_m, new_m, ba_n, new_n, new_ge, ba_v, new_v, new_gi

    not_refract = merge(1_int8, 0_int8, timestep(current_t - lastspike, dt) >= lio_1)
    ba_h = (lio_2 * exp(lio_3 * v)) / (((-4.0_real32) / (0.001_real32 + (lio_4 * exp(lio_5 * v)))) - &
           (lio_2 * exp(lio_3 * v)))
    new_h = (-ba_h) + ((ba_h + h) * exp(dt * (((-4.0_real32) / (0.001_real32 + (lio_4 * exp(lio_5 * v)))) - &
            (lio_2 * exp(lio_3 * v)))))
    ba_m = (((lio_6 / (lio_7 + (lio_8 * exp(lio_9 * v)))) + (lio_10 / (lio_7 + (lio_8 * exp(lio_9 * v))))) - &
           ((0.32_real32 * v) / (lio_7 + (lio_8 * exp(lio_9 * v))))) / (((((lio_11 / (lio_7 + &
           (lio_8 * exp(lio_9 * v)))) + (lio_12 / (lio_13 + (lio_14 * exp(lio_15 * v))))) + &
           (lio_16 / (lio_13 + (lio_14 * exp(lio_15 * v))))) + ((0.32_real32 * v) / (lio_7 + &
           (lio_8 * exp(lio_9 * v))))) - ((lio_10 / (lio_7 + (lio_8 * exp(lio_9 * v)))) + &
           ((0.28_real32 * v) / (lio_13 + (lio_14 * exp(lio_15 * v))))))
    new_m = (-ba_m) + ((ba_m + m) * exp(dt * (((((lio_11 / (lio_7 + (lio_8 * exp(lio_9 * v)))) + &
            (lio_12 / (lio_13 + (lio_14 * exp(lio_15 * v))))) + (lio_16 / (lio_13 + &
            (lio_14 * exp(lio_15 * v))))) + ((0.32_real32 * v) / (lio_7 + (lio_8 * exp(lio_9 * v))))) - &
            ((lio_10 / (lio_7 + (lio_8 * exp(lio_9 * v)))) + ((0.28_real32 * v) / (lio_13 + &
            (lio_14 * exp(lio_15 * v))))))))
    ba_n = (((lio_17 / (lio_7 + (lio_18 * exp(lio_5 * v)))) + (lio_19 / (lio_7 + (lio_18 * exp(lio_5 * v))))) - &
           ((0.032_real32 * v) / (lio_7 + (lio_18 * exp(lio_5 * v))))) / (((lio_20 / (lio_7 + &
           (lio_18 * exp(lio_5 * v)))) + ((0.032_real32 * v) / (lio_7 + (lio_18 * exp(lio_5 * v))))) - &
           ((lio_19 / (lio_7 + (lio_18 * exp(lio_5 * v)))) + (lio_21 * exp(lio_22 * v))))
    new_n = (-ba_n) + ((ba_n + n) * exp(dt * (((lio_20 / (lio_7 + (lio_18 * exp(lio_5 * v)))) + &
            ((0.032_real32 * v) / (lio_7 + (lio_18 * exp(lio_5 * v))))) - ((lio_19 / (lio_7 + &
            (lio_18 * exp(lio_5 * v)))) + (lio_21 * exp(lio_22 * v))))))
    new_ge = lio_23 * ge
    ba_v = (lio_24 + ((((lio_25 * (n * n * n * n)) + (lio_26 * (h * (m * m * m)))) + (lio_27 * ge)) + &
           (lio_28 * gi))) / ((lio_29 + (lio_30 * (n * n * n * n))) - (((lio_31 * (h * (m * m * m))) + &
           (lio_32 * ge)) + (lio_32 * gi)))
    new_v = (-ba_v) + ((ba_v + v) * exp(dt * ((lio_29 + (lio_30 * (n * n * n * n))) - (((lio_31 * &
            (h * (m * m * m))) + (lio_32 * ge)) + (lio_32 * gi)))))
    new_gi = lio_33 * gi

    h = new_h
    m = new_m
    n = new_n
    ge = new_ge
    v = new_v
    gi = new_gi
  end subroutine update_one

  real(real64) function compute_rsme(ge, gi, h, m, n, v, not_refract, ref_ge, ref_gi, ref_h, ref_m, ref_n, ref_v, &
                                     ref_not_refract, neurons)
    real(real32), intent(in) :: ge(:), gi(:), h(:), m(:), n(:), v(:)
    real(real32), intent(in) :: ref_ge(:), ref_gi(:), ref_h(:), ref_m(:), ref_n(:), ref_v(:)
    integer(int8), intent(in) :: not_refract(:), ref_not_refract(:)
    integer, intent(in) :: neurons
    integer :: i
    real(real64) :: accum

    accum = 0.0_real64
    do i = 1, neurons
      accum = accum + real(ge(i) - ref_ge(i), real64) ** 2
      accum = accum + real(gi(i) - ref_gi(i), real64) ** 2
      accum = accum + real(h(i) - ref_h(i), real64) ** 2
      accum = accum + real(m(i) - ref_m(i), real64) ** 2
      accum = accum + real(n(i) - ref_n(i), real64) ** 2
      accum = accum + real(v(i) - ref_v(i), real64) ** 2
      accum = accum + real(int(not_refract(i)) - int(ref_not_refract(i)), real64) ** 2
    end do
    compute_rsme = sqrt(accum / real(neurons, real64))
  end function compute_rsme

end program main
