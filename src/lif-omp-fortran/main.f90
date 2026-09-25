! SPDX-License-Identifier: CC0-1.0
program lif_main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer :: neurons_per_item, num_items, num_steps
  integer :: num_neurons, i, step
  integer(c_int), parameter :: RAND_MAX_C = 2147483647_c_int
  real(real32), parameter :: dt = 0.1_real32
  real(real32), parameter :: tau_rc = 10.0_real32
  real(real32), parameter :: tau_ref = 2.0_real32
  real(real32), allocatable :: encode_result(:), voltage(:), reftime(:), spikes(:)
  real(real32), allocatable :: voltage_host(:), reftime_host(:), spikes_host(:)
  real(real32), allocatable :: bias(:), gain(:)
  real(real32), allocatable :: spike_reftime(:), spike_reftime_host(:)
  integer :: num_spikes, num_spikes_host, compare_count
  real(real64) :: start_time, elapsed_us
  character(len=64) :: arg
  logical :: ok

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name="rand")
      import :: c_int
    end function c_rand
  end interface

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <neurons per item> <num_items> <num_steps>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) neurons_per_item
  call get_command_argument(2, arg)
  read(arg, *) num_items
  call get_command_argument(3, arg)
  read(arg, *) num_steps

  num_neurons = neurons_per_item * num_items
  allocate(encode_result(0:num_items-1), bias(0:neurons_per_item-1), gain(0:neurons_per_item-1))
  allocate(voltage(0:num_neurons-1), reftime(0:num_neurons-1), spikes(0:num_neurons-1))
  allocate(voltage_host(0:num_neurons-1), reftime_host(0:num_neurons-1), spikes_host(0:num_neurons-1))
  allocate(spike_reftime(0:num_neurons-1), spike_reftime_host(0:num_neurons-1))

  call c_srand(123_c_int)
  do i = 0, num_items - 1
    encode_result(i) = rand_unit()
  end do
  do i = 0, num_neurons - 1
    voltage(i) = 1.0_real32 + rand_unit()
    voltage_host(i) = voltage(i)
    reftime(i) = real(mod(c_rand(), 5_c_int), real32) / 10.0_real32
    reftime_host(i) = reftime(i)
  end do
  do i = 0, neurons_per_item - 1
    bias(i) = rand_unit()
    gain(i) = rand_unit() + 0.5_real32
  end do

  !$omp target data map(to: encode_result, bias, gain) map(from: spikes) map(tofrom: voltage, reftime)
  start_time = omp_get_wtime()
  do step = 1, num_steps
    call lif_step(num_neurons, neurons_per_item, dt, encode_result, voltage, reftime, tau_rc, tau_ref, bias, gain, spikes)
  end do
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(num_steps, real64)
  !$omp end target data
  write(*,'("Average kernel execution time: ",F0.6," (us)")') elapsed_us

  do step = 1, num_steps
    call lif_reference(num_neurons, neurons_per_item, dt, encode_result, voltage_host, reftime_host, &
        tau_rc, tau_ref, bias, gain, spikes_host)
  end do

  num_spikes = 0
  num_spikes_host = 0
  do i = 0, num_neurons - 1
    if (spikes(i) == 1.0_real32 / dt) then
      spike_reftime(num_spikes) = reftime(i)
      num_spikes = num_spikes + 1
    end if
    if (spikes_host(i) == 1.0_real32 / dt) then
      spike_reftime_host(num_spikes_host) = reftime_host(i)
      num_spikes_host = num_spikes_host + 1
    end if
  end do

  write(*,'("Number of spikes on host and device: ",I0," ",I0)') num_spikes_host, num_spikes
  compare_count = min(num_spikes, num_spikes_host)
  ok = .true.
  do i = 0, compare_count - 1
    if (abs(spike_reftime(i) - spike_reftime_host(i)) > 0.1_real32) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine lif_step(num_neurons, neurons_per_item, dt, encode_result, voltage_array, reftime_array, &
      tau_rc, tau_ref, bias, gain, spikes)
    integer, intent(in) :: num_neurons, neurons_per_item
    real(real32), intent(in) :: dt, tau_rc, tau_ref
    real(real32), intent(in) :: encode_result(0:), bias(0:), gain(0:)
    real(real32), intent(inout) :: voltage_array(0:), reftime_array(0:)
    real(real32), intent(out) :: spikes(0:)
    integer :: idx, neuron_index, item_index
    real(real32) :: voltage, ref_time, current, dv, spike, mult

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(neuron_index, item_index, voltage, ref_time, current, dv, spike, mult)
    do idx = 0, num_neurons - 1
      neuron_index = mod(idx, neurons_per_item)
      item_index = idx / neurons_per_item
      voltage = voltage_array(idx)
      ref_time = reftime_array(idx)
      current = bias(neuron_index) + gain(neuron_index) * encode_result(item_index)
      dv = -(exp(-dt / tau_rc) - 1.0_real32) * (current - voltage)
      voltage = max(voltage + dv, 0.0_real32)
      ref_time = ref_time - dt
      mult = ref_time
      mult = mult * (-1.0_real32 / dt)
      mult = mult + 1.0_real32
      mult = min(mult, 1.0_real32)
      mult = max(mult, 0.0_real32)
      voltage = voltage * mult
      if (voltage > 1.0_real32) then
        spike = 1.0_real32 / dt
        ref_time = tau_ref + dt * (1.0_real32 - (voltage - 1.0_real32) / dv)
        voltage = 0.0_real32
      else
        spike = 0.0_real32
      end if
      reftime_array(idx) = ref_time
      voltage_array(idx) = voltage
      spikes(idx) = spike
    end do
    !$omp end target teams distribute parallel do
  end subroutine lif_step

  subroutine lif_reference(num_neurons, neurons_per_item, dt, encode_result, voltage_array, reftime_array, &
      tau_rc, tau_ref, bias, gain, spikes)
    integer, intent(in) :: num_neurons, neurons_per_item
    real(real32), intent(in) :: dt, tau_rc, tau_ref
    real(real32), intent(in) :: encode_result(0:), bias(0:), gain(0:)
    real(real32), intent(inout) :: voltage_array(0:), reftime_array(0:)
    real(real32), intent(out) :: spikes(0:)
    integer :: idx, neuron_index, item_index
    real(real32) :: voltage, ref_time, current, dv, spike, mult

    do idx = 0, num_neurons - 1
      neuron_index = mod(idx, neurons_per_item)
      item_index = idx / neurons_per_item
      voltage = voltage_array(idx)
      ref_time = reftime_array(idx)
      current = bias(neuron_index) + gain(neuron_index) * encode_result(item_index)
      dv = -(exp(-dt / tau_rc) - 1.0_real32) * (current - voltage)
      voltage = max(voltage + dv, 0.0_real32)
      ref_time = ref_time - dt
      mult = max(min(ref_time * (-1.0_real32 / dt) + 1.0_real32, 1.0_real32), 0.0_real32)
      voltage = voltage * mult
      if (voltage > 1.0_real32) then
        spike = 1.0_real32 / dt
        ref_time = tau_ref + dt * (1.0_real32 - (voltage - 1.0_real32) / dv)
        voltage = 0.0_real32
      else
        spike = 0.0_real32
      end if
      reftime_array(idx) = ref_time
      voltage_array(idx) = voltage
      spikes(idx) = spike
    end do
  end subroutine lif_reference

  real(real32) function rand_unit()
    rand_unit = real(c_rand(), real32) / real(RAND_MAX_C, real32)
  end function rand_unit

end program lif_main
