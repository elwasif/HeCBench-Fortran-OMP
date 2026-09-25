! SPDX-License-Identifier: CC0-1.0
program rfs_main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

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

  integer(c_int), parameter :: c_rand_max = 2147483647_c_int
  integer :: n_arrays, n_elems, n, i
  integer(c_int) :: rand_value
  real(real32), allocatable :: arrays(:), max_val(:), result(:), factor(:), result_ref(:)
  real(real64) :: start_time, elapsed
  character(len=64) :: arg
  logical :: ok

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <number of arrays> <length of each array>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) n_arrays
  call get_command_argument(2, arg)
  read(arg, *) n_elems

  allocate(arrays(0:n_arrays*n_elems-1), max_val(0:n_arrays-1), result(0:n_arrays-1))
  allocate(factor(0:n_arrays-1), result_ref(0:n_arrays-1))

  call c_srand(123_c_int)
  do n = 0, n_arrays - 1
    max_val(n) = 0.0_real32
    do i = 0, n_elems - 1
      rand_value = c_rand()
      arrays(n * n_elems + i) = real(rand_value, real32) / real(c_rand_max, real32)
      rand_value = c_rand()
      if (mod(rand_value, 2_c_int) /= 0) arrays(n * n_elems + i) = -arrays(n * n_elems + i)
      max_val(n) = max(abs(arrays(n * n_elems + i)), max_val(n))
    end do
    factor(n) = create_rounding_factor(max_val(n), n_elems)
  end do

  do n = 0, n_arrays - 1
    result_ref(n) = 0.0_real32
    do i = 0, n_elems - 1
      result_ref(n) = result_ref(n) + truncate_value(factor(n), arrays(n * n_elems + i))
    end do
  end do

  !$omp target data map(to: arrays, max_val) map(alloc: result)
  !$omp target teams distribute parallel do thread_limit(256)
  do i = 0, n_arrays - 1
    result(i) = 0.0_real32
  end do
  !$omp end target teams distribute parallel do

  start_time = omp_get_wtime()
  do n = 0, n_arrays - 1
    call sum_array(factor(n), n_elems, arrays, n * n_elems, result, n)
  end do
  elapsed = (omp_get_wtime() - start_time) / real(n_arrays, real64)
  write(*,'("Average kernel execution time (sumArray): ",F0.6," (s)")') elapsed

  !$omp target update from(result)
  ok = exact_match(result_ref, result, n_arrays)
  write(*,'(A)') merge("PASS", "FAIL", ok)

  start_time = omp_get_wtime()
  call sum_arrays(n_arrays, n_elems, arrays, result, max_val)
  elapsed = omp_get_wtime() - start_time
  write(*,'("Kernel execution time (sumArrays): ",F0.6," (s)")') elapsed

  !$omp target update from(result)
  !$omp end target data
  ok = exact_match(result_ref, result, n_arrays)
  write(*,'(A)') merge("PASS", "FAIL", ok)

contains

  subroutine sum_array(factor, length, values, offset, result, result_index)
    real(real32), intent(in) :: factor
    integer, intent(in) :: length, offset, result_index
    real(real32), intent(in) :: values(0:)
    real(real32), intent(inout) :: result(0:)
    integer :: idx
    real(real32) :: q

    !$omp target teams distribute parallel do num_teams(256) thread_limit(256) private(q)
    do idx = 0, length - 1
      q = (factor + values(offset + idx)) - factor
      !$omp atomic update
      result(result_index) = result(result_index) + q
    end do
    !$omp end target teams distribute parallel do
  end subroutine sum_array

  subroutine sum_arrays(n_arrays, length, values, result, max_val)
    integer, intent(in) :: n_arrays, length
    real(real32), intent(in) :: values(0:), max_val(0:)
    real(real32), intent(out) :: result(0:)
    integer :: array_idx, idx, exp_val
    real(real32) :: delta, factor, accum, q

    !$omp target teams distribute parallel do num_teams(256) thread_limit(256) &
    !$omp& private(idx, exp_val, delta, factor, accum, q)
    do array_idx = 0, n_arrays - 1
      delta = (max_val(array_idx) * real(length, real32)) / &
          (1.0_real32 - 2.0_real32 * real(length, real32) * epsilon(1.0_real32))
      exp_val = exponent(delta)
      factor = scale(1.0_real32, exp_val)
      accum = 0.0_real32
      do idx = length - 1, 0, -1
        q = (factor + values(array_idx * length + idx)) - factor
        accum = accum + q
      end do
      result(array_idx) = accum
    end do
    !$omp end target teams distribute parallel do
  end subroutine sum_arrays

  real(real32) function create_rounding_factor(max_value, length)
    real(real32), intent(in) :: max_value
    integer, intent(in) :: length
    real(real32) :: delta

    delta = (max_value * real(length, real32)) / &
        (1.0_real32 - 2.0_real32 * real(length, real32) * epsilon(1.0_real32))
    create_rounding_factor = scale(1.0_real32, exponent(delta))
  end function create_rounding_factor

  real(real32) function truncate_value(factor, value)
    real(real32), intent(in) :: factor, value
    truncate_value = (factor + value) - factor
  end function truncate_value

  logical function exact_match(reference_values, device_values, size)
    real(real32), intent(in) :: reference_values(0:), device_values(0:)
    integer, intent(in) :: size
    integer :: idx

    exact_match = .true.
    do idx = 0, size - 1
      if (transfer(reference_values(idx), 0_int32) /= transfer(device_values(idx), 0_int32)) then
        exact_match = .false.
        exit
      end if
    end do
  end function exact_match

end program rfs_main
