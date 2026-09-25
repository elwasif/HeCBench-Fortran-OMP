! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  character(len=128) :: arg1, arg2, arg3
  integer :: n, d, repeat, i
  real(real32), allocatable :: key(:), value(:), query(:), host_out(:), device_out(:)
  logical :: ok

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <rows> <columns> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) n
  read(arg2, *) d
  read(arg3, *) repeat
  if (n <= 0 .or. d <= 0 .or. repeat <= 0) stop 1

  allocate(key(n * d), value(n * d), query(d), host_out(d), device_out(d))
  call initialize_inputs(key, value, query, n, d)

  call attention_host(key, value, query, host_out, n, d)
  call attention_device(key, value, query, device_out, n, d, repeat)

  ok = .true.
  do i = 1, d
    if (abs(host_out(i) - device_out(i)) > 1.0e-3_real32) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(key, value, query, host_out, device_out)

contains

  subroutine initialize_inputs(key, value, query, n, d)
    real(real32), intent(out) :: key(:), value(:), query(:)
    integer, intent(in) :: n, d
    integer :: i, j, idx

    do i = 1, n
      do j = 1, d
        idx = (i - 1) * d + j
        key(idx) = real(mod(37 * idx + 11, 2001), real32) * 1.0e-5_real32 - 0.01_real32
        value(idx) = real(mod(53 * idx + 7, 2001), real32) * 1.0e-5_real32 - 0.01_real32
      end do
    end do
    do j = 1, d
      query(j) = real(mod(97 * j + 3, 2001), real32) * 1.0e-5_real32 - 0.01_real32
    end do
  end subroutine initialize_inputs

  subroutine attention_host(key, value, query, output, n, d)
    real(real32), intent(in) :: key(:), value(:), query(:)
    real(real32), intent(out) :: output(:)
    integer, intent(in) :: n, d
    real(real32), allocatable :: dot_product(:), score(:)
    real(real32) :: sum_value
    integer :: i, j, idx

    allocate(dot_product(n), score(n))
    do i = 1, n
      sum_value = 0.0_real32
      do j = 1, d
        idx = (i - 1) * d + j
        sum_value = sum_value + key(idx) * query(j)
      end do
      dot_product(i) = sum_value
    end do

    sum_value = 0.0_real32
    do i = 1, n
      sum_value = sum_value + exp(dot_product(i))
    end do

    do i = 1, n
      score(i) = exp(dot_product(i)) / sum_value
    end do

    do j = 1, d
      sum_value = 0.0_real32
      do i = 1, n
        idx = (i - 1) * d + j
        sum_value = sum_value + score(i) * value(idx)
      end do
      output(j) = sum_value
    end do

    deallocate(dot_product, score)
  end subroutine attention_host

  subroutine attention_device(key, value, query, output, n, d, repeat)
    real(real32), intent(in) :: key(:), value(:), query(:)
    real(real32), intent(out) :: output(:)
    integer, intent(in) :: n, d, repeat
    real(real32), allocatable :: dot_product(:), score(:)
    real(real32) :: exp_sum
    real(real64) :: start_time, end_time, elapsed_ms
    integer :: iter, i, j, idx
    real(real32) :: local_sum

    allocate(dot_product(n), score(n))
    dot_product = 0.0_real32
    score = 0.0_real32
    output = 0.0_real32

    !$omp target data map(to: key(1:n*d), value(1:n*d), query(1:d)) &
    !$omp& map(alloc: dot_product(1:n), score(1:n), exp_sum) map(from: output(1:d))
      start_time = omp_get_wtime()
      do iter = 1, repeat
        exp_sum = 0.0_real32
        !$omp target update to(exp_sum)

        !$omp target teams distribute parallel do thread_limit(256) private(i, j, idx, local_sum)
        do i = 1, n
          local_sum = 0.0_real32
          do j = 1, d
            idx = (i - 1) * d + j
            local_sum = local_sum + key(idx) * query(j)
          end do
          dot_product(i) = local_sum
          !$omp atomic update
          exp_sum = exp_sum + exp(local_sum)
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do thread_limit(256) private(i)
        do i = 1, n
          score(i) = exp(dot_product(i)) / exp_sum
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do thread_limit(256) private(i, j, idx, local_sum)
        do j = 1, d
          local_sum = 0.0_real32
          do i = 1, n
            idx = (i - 1) * d + j
            local_sum = local_sum + score(i) * value(idx)
          end do
          output(j) = local_sum
        end do
        !$omp end target teams distribute parallel do
      end do
      end_time = omp_get_wtime()
    !$omp end target data

    elapsed_ms = (end_time - start_time) * 1.0e3_real64 / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of kernels ', elapsed_ms, ' (ms)'
    deallocate(dot_product, score)
  end subroutine attention_device

end program main
