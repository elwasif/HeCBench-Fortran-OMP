! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  interface
    function c_log2f(x) bind(C, name='log2f') result(y)
      import :: c_float
      real(c_float), value :: x
      real(c_float) :: y
    end function c_log2f
  end interface

  character(len=512) :: config_path
  character(len=64) :: placeholder
  integer :: unit, repeat, precision_count
  integer(int64) :: ceiling_val, current
  integer :: num_inputs, increment, i, pidx, iter
  integer, allocatable :: precisions(:)
  real(real32), allocatable :: inputs(:), outputs(:, :), ref_vals(:)
  real(real64) :: start_time, end_time, elapsed_us
  real(real32) :: rmse

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <config filename>'
    stop 1
  end if

  call get_command_argument(1, config_path)
  open(newunit=unit, file=trim(config_path), status='old', action='read')
  read(unit, *) placeholder, ceiling_val
  read(unit, *) placeholder, repeat
  read(unit, *) placeholder, precision_count
  allocate(precisions(precision_count))
  read(unit, *) placeholder, (precisions(i), i = 1, precision_count)
  close(unit)

  increment = 1
  num_inputs = int((ceiling_val + int(increment, int64) - 1_int64) / int(increment, int64))
  allocate(inputs(num_inputs), outputs(num_inputs, precision_count), ref_vals(num_inputs))

  current = 1_int64
  do i = 1, num_inputs
    inputs(i) = real(current, real32)
    ref_vals(i) = c_log2f(inputs(i))
    current = current + int(increment, int64)
  end do
  outputs = 0.0_real32

  print '(A,I0)', 'Number of precision counts : ', precision_count
  print '(A,I0)', ' Number of inputs to evaluate for each precision: ', num_inputs
  print '(A,I0)', ' Number of runs for each precision : ', repeat

  !$omp target data map(to: inputs(1:num_inputs), precisions(1:precision_count)) map(from: outputs(1:num_inputs,1:precision_count))
  do pidx = 1, precision_count
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call compute_log_for_precision(inputs, outputs, num_inputs, pidx, precisions(pidx))
    end do
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A)'
    print '(A,I0,A)', 'Iterative approximation with ', precisions(pidx), ' bits of precision'
    print '(A,G0.6,A)', 'Average kernel execution time ', elapsed_us, ' (us)'
  end do
  !$omp end target data

  print '(A)', '-------------- SUMMARY (Device results): --------------'
  print '(A)'
  do pidx = 1, precision_count
    print '(A,I0,A)', '----- Iterative approximation with ', precisions(pidx), ' bits of precision -----'
    rmse = compute_rmse(outputs(:, pidx), ref_vals, num_inputs)
    call print_cxx_default_float('RMSE : ', rmse)
  end do

  deallocate(inputs, outputs, ref_vals, precisions)

contains

  subroutine compute_log_for_precision(inputs, outputs, num_inputs, pidx, precision)
    real(real32), intent(in) :: inputs(:)
    real(real32), intent(inout) :: outputs(:, :)
    integer, intent(in) :: num_inputs, pidx, precision
    integer :: j

    !$omp target teams distribute parallel do thread_limit(256) private(j) firstprivate(num_inputs, pidx, precision)
    do j = 1, num_inputs
      outputs(j, pidx) = binary_log(inputs(j), precision)
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_log_for_precision

  real(real32) function binary_log(input, precision)
    real(real32), intent(in) :: input
    integer, intent(in) :: precision
    integer(int32) :: bits, exponent, m, sum_m, test
    integer(int64) :: denom, prev_denom
    real(real32) :: y, result
    logical :: max_condition_met

    bits = transfer(input, bits)
    exponent = int(ishft(iand(bits, int(Z'7F800000', int32)), -23), int32) - 127_int32
    m = 0_int32
    sum_m = 0_int32
    result = 0.0_real32
    test = ishft(1_int32, exponent)
    y = input / real(test, real32)
    max_condition_met = .false.
    denom = 0_int64
    prev_denom = 0_int64

    do while (((sum_m < precision + 1) .and. (y /= 1.0_real32)) .or. max_condition_met)
      m = 0_int32
      do while ((y < 2.0_real32) .and. (sum_m + m < precision + 1))
        y = y * y
        m = m + 1_int32
      end do

      sum_m = sum_m + m
      prev_denom = denom
      denom = ishft(1_int64, int(sum_m))

      if (sum_m >= precision) exit
      if (prev_denom > denom) then
        max_condition_met = .true.
        exit
      end if

      result = result + 1.0_real32 / real(denom, real32)
      y = y / 2.0_real32
    end do

    binary_log = real(exponent, real32) + result
  end function binary_log

  real(real32) function compute_rmse(values, refs, n)
    real(real32), intent(in) :: values(:), refs(:)
    integer, intent(in) :: n
    integer :: j
    real(real32) :: sum_sq, delta

    sum_sq = 0.0_real32
    do j = 1, n
      delta = values(j) - refs(j)
      sum_sq = sum_sq + delta * delta
    end do
    sum_sq = sum_sq / real(n, real32)
    compute_rmse = sqrt(sum_sq)
  end function compute_rmse

  subroutine print_cxx_default_float(prefix, value)
    character(len=*), intent(in) :: prefix
    real(real32), intent(in) :: value
    character(len=32) :: fmt
    character(len=96) :: number
    integer :: exponent, decimals
    real(real64) :: magnitude

    if (value == 0.0_real32) then
      print '(A,A)', prefix, '0'
      return
    end if

    magnitude = abs(real(value, real64))
    exponent = floor(log10(magnitude))
    decimals = max(0, 6 - exponent - 1)
    write(fmt, '(A,I0,A)') '(F64.', decimals, ')'
    write(number, fmt) value
    print '(A,A)', prefix, trim(adjustl(number))
  end subroutine print_cxx_default_float

end program main
