! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name="rand")
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer, parameter :: map_size = 1024
  integer, parameter :: out_map_size = map_size - 2
  integer, parameter :: local_x = 32
  integer, parameter :: local_y = 8
  integer, parameter :: tile_n = (out_map_size + 1) / 2
  integer, parameter :: global_x = ((tile_n + local_x - 1) / local_x) * local_x
  integer, parameter :: global_y = ((tile_n + local_y - 1) / local_y) * local_y
  integer, parameter :: thread_size = local_x * local_y
  integer, parameter :: input_size = map_size * map_size
  integer, parameter :: output_size = out_map_size * out_map_size
  real(real32), parameter :: percent_diff_error_threshold = 1.05_real32
  real(real32), parameter :: small_float_val = 0.00000001_real32

  real(real32), allocatable :: a(:), b(:), b_host(:), c(:)
  real(real64) :: start_time, end_time, co_start, co_time
  real(real64) :: total_time, ratio
  integer :: cpu_offset, cpu_global_x, gpu_global_x, offset_x, update_count
  logical :: pass, cpu_run, gpu_run

  start_time = omp_get_wtime()
  allocate(a(input_size), b(output_size), b_host(output_size), c(16))

  call initialize_input(a)
  call filter_transformation(c)

  pass = .true.
  co_time = 0.0_real64

  !$omp target data map(to: a(1:input_size), c(1:16)) map(alloc: b(1:output_size))
  do cpu_offset = 0, 100
    cpu_global_x = (cpu_offset * (global_x / local_x) / 100) * local_x
    gpu_global_x = global_x - cpu_global_x
    offset_x = cpu_global_x
    cpu_run = cpu_global_x > 0
    gpu_run = gpu_global_x > 0

    co_start = omp_get_wtime()

    if (gpu_run) then
      call winograd_device(a, b, c, gpu_global_x, global_y, offset_x)
    end if

    if (cpu_run) then
      call winograd_cpu_region(a, b, c, cpu_global_x)
      if (gpu_run) then
        update_count = min(offset_x * 2 * out_map_size, output_size)
        !$omp target update to(b(1:update_count))
      else
        !$omp target update to(b(1:output_size))
      end if
    end if

    !$omp target update from(b(1:output_size))
    co_time = co_time + omp_get_wtime() - co_start

    call winograd_reference(a, b_host, c)
    pass = pass .and. compare_results(b_host, b)
  end do
  !$omp end target data

  write(*,'(A)') merge('PASS', 'FAIL', pass)

  end_time = omp_get_wtime()
  total_time = end_time - start_time
  ratio = 100.0_real64 * co_time / total_time
  write(*,'(A,F0.6,A)') 'Co-execution time: ', co_time, ' s'
  write(*,'(A,F0.6,A)') 'Total time: ', total_time, ' s'
  write(*,'(A,F0.2,A)') 'Ratio of co-execution time to total time: ', ratio, '%'

  deallocate(a, b, b_host, c)

contains

  subroutine initialize_input(a)
    real(real32), intent(out) :: a(:)
    integer :: idx
    call c_srand(1_c_int)
    do idx = 1, size(a)
      a(idx) = real(c_rand(), real32) / 2147483647.0_real32
    end do
  end subroutine initialize_input

  subroutine filter_transformation(transformed_filter)
    real(real32), intent(out) :: transformed_filter(:)
    real(real32) :: filter(3, 3), tmp_filter(4, 3)
    integer :: i, j

    filter(1, 1) = 0.2_real32
    filter(2, 1) = 0.5_real32
    filter(3, 1) = -0.8_real32
    filter(1, 2) = -0.3_real32
    filter(2, 2) = 0.6_real32
    filter(3, 2) = -0.9_real32
    filter(1, 3) = 0.4_real32
    filter(2, 3) = 0.7_real32
    filter(3, 3) = 0.10_real32

    do j = 1, 3
      tmp_filter(1, j) = filter(1, j)
      tmp_filter(2, j) = 0.5_real32 * filter(1, j) + 0.5_real32 * filter(2, j) + 0.5_real32 * filter(3, j)
      tmp_filter(3, j) = 0.5_real32 * filter(1, j) - 0.5_real32 * filter(2, j) + 0.5_real32 * filter(3, j)
      tmp_filter(4, j) = filter(3, j)
    end do

    do i = 1, 4
      transformed_filter((i - 1) * 4 + 1) = tmp_filter(i, 1)
      transformed_filter((i - 1) * 4 + 2) = 0.5_real32 * tmp_filter(i, 1) + &
          0.5_real32 * tmp_filter(i, 2) + 0.5_real32 * tmp_filter(i, 3)
      transformed_filter((i - 1) * 4 + 3) = 0.5_real32 * tmp_filter(i, 1) - &
          0.5_real32 * tmp_filter(i, 2) + 0.5_real32 * tmp_filter(i, 3)
      transformed_filter((i - 1) * 4 + 4) = tmp_filter(i, 3)
    end do
  end subroutine filter_transformation

  subroutine winograd_device(input, output, transformed_filter, tile_i_size, tile_j_size, offset_i)
    real(real32), intent(in) :: input(:), transformed_filter(:)
    real(real32), intent(inout) :: output(:)
    integer, intent(in) :: tile_i_size, tile_j_size, offset_i
    integer :: tile_i, tile_j, i, j, x, y, out_idx
    real(real32) :: input_tile(4, 4), tmp_tile(4, 4), transformed_tile(4, 4)
    real(real32) :: multiplied_tile(4, 4), tmp_tile_1(2, 4), final_tile(2, 2)

    !$omp target teams distribute parallel do collapse(2) thread_limit(thread_size) &
    !$omp& private(tile_i, tile_j, i, j, x, y, out_idx, input_tile, tmp_tile, transformed_tile, &
    !$omp& multiplied_tile, tmp_tile_1, final_tile)
    do tile_j = 0, tile_j_size - 1
      do tile_i = 0, tile_i_size - 1
        do i = 1, 4
          do j = 1, 4
            x = 2 * (tile_i + offset_i) + (i - 1)
            y = 2 * tile_j + (j - 1)
            if (x >= map_size .or. y >= map_size) then
              input_tile(i, j) = 0.0_real32
            else
              input_tile(i, j) = input(x * map_size + y + 1)
            end if
          end do
        end do

        do j = 1, 4
          tmp_tile(1, j) = input_tile(1, j) - input_tile(3, j)
          tmp_tile(2, j) = input_tile(2, j) + input_tile(3, j)
          tmp_tile(3, j) = -input_tile(2, j) + input_tile(3, j)
          tmp_tile(4, j) = input_tile(2, j) - input_tile(4, j)
        end do
        do i = 1, 4
          transformed_tile(i, 1) = tmp_tile(i, 1) - tmp_tile(i, 3)
          transformed_tile(i, 2) = tmp_tile(i, 2) + tmp_tile(i, 3)
          transformed_tile(i, 3) = -tmp_tile(i, 2) + tmp_tile(i, 3)
          transformed_tile(i, 4) = tmp_tile(i, 2) - tmp_tile(i, 4)
        end do

        do i = 1, 4
          do j = 1, 4
            multiplied_tile(i, j) = transformed_tile(i, j) * transformed_filter((i - 1) * 4 + j)
          end do
        end do

        do j = 1, 4
          tmp_tile_1(1, j) = multiplied_tile(1, j) + multiplied_tile(2, j) + multiplied_tile(3, j)
          tmp_tile_1(2, j) = multiplied_tile(2, j) - multiplied_tile(3, j) - multiplied_tile(4, j)
        end do
        do i = 1, 2
          final_tile(i, 1) = tmp_tile_1(i, 1) + tmp_tile_1(i, 2) + tmp_tile_1(i, 3)
          final_tile(i, 2) = tmp_tile_1(i, 2) - tmp_tile_1(i, 3) - tmp_tile_1(i, 4)
        end do

        do i = 1, 2
          do j = 1, 2
            x = 2 * (tile_i + offset_i) + (i - 1)
            y = 2 * tile_j + (j - 1)
            if (x < out_map_size .and. y < out_map_size) then
              out_idx = x * out_map_size + y + 1
              output(out_idx) = final_tile(i, j)
            end if
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine winograd_device

  subroutine winograd_cpu_region(input, output, transformed_filter, cpu_global_x)
    real(real32), intent(in) :: input(:), transformed_filter(:)
    real(real32), intent(inout) :: output(:)
    integer, intent(in) :: cpu_global_x
    integer :: tile_i, tile_j, i, j, x, y, out_idx
    real(real32) :: input_tile(4, 4), tmp_tile(4, 4), transformed_tile(4, 4)
    real(real32) :: multiplied_tile(4, 4), tmp_tile_1(2, 4), final_tile(2, 2)

    !$omp parallel do private(tile_i, tile_j, i, j, x, y, out_idx, input_tile, tmp_tile, transformed_tile, &
    !$omp& multiplied_tile, tmp_tile_1, final_tile)
    do tile_i = 0, cpu_global_x - 1
      do tile_j = 0, tile_n - 1
        call compute_tile_host(input, output, transformed_filter, tile_i, tile_j, 0)
      end do
    end do
    !$omp end parallel do
  end subroutine winograd_cpu_region

  subroutine winograd_reference(input, output, transformed_filter)
    real(real32), intent(in) :: input(:), transformed_filter(:)
    real(real32), intent(out) :: output(:)
    integer :: tile_i, tile_j
    do tile_i = 0, tile_n - 1
      do tile_j = 0, tile_n - 1
        call compute_tile_host(input, output, transformed_filter, tile_i, tile_j, 0)
      end do
    end do
  end subroutine winograd_reference

  subroutine compute_tile_host(input, output, transformed_filter, tile_i, tile_j, offset_i)
    real(real32), intent(in) :: input(:), transformed_filter(:)
    real(real32), intent(inout) :: output(:)
    integer, intent(in) :: tile_i, tile_j, offset_i
    integer :: i, j, x, y, out_idx
    real(real32) :: input_tile(4, 4), tmp_tile(4, 4), transformed_tile(4, 4)
    real(real32) :: multiplied_tile(4, 4), tmp_tile_1(2, 4), final_tile(2, 2)

    do i = 1, 4
      do j = 1, 4
        x = 2 * (tile_i + offset_i) + (i - 1)
        y = 2 * tile_j + (j - 1)
        if (x >= map_size .or. y >= map_size) then
          input_tile(i, j) = 0.0_real32
        else
          input_tile(i, j) = input(x * map_size + y + 1)
        end if
      end do
    end do

    do j = 1, 4
      tmp_tile(1, j) = input_tile(1, j) - input_tile(3, j)
      tmp_tile(2, j) = input_tile(2, j) + input_tile(3, j)
      tmp_tile(3, j) = -input_tile(2, j) + input_tile(3, j)
      tmp_tile(4, j) = input_tile(2, j) - input_tile(4, j)
    end do
    do i = 1, 4
      transformed_tile(i, 1) = tmp_tile(i, 1) - tmp_tile(i, 3)
      transformed_tile(i, 2) = tmp_tile(i, 2) + tmp_tile(i, 3)
      transformed_tile(i, 3) = -tmp_tile(i, 2) + tmp_tile(i, 3)
      transformed_tile(i, 4) = tmp_tile(i, 2) - tmp_tile(i, 4)
    end do

    do i = 1, 4
      do j = 1, 4
        multiplied_tile(i, j) = transformed_tile(i, j) * transformed_filter((i - 1) * 4 + j)
      end do
    end do

    do j = 1, 4
      tmp_tile_1(1, j) = multiplied_tile(1, j) + multiplied_tile(2, j) + multiplied_tile(3, j)
      tmp_tile_1(2, j) = multiplied_tile(2, j) - multiplied_tile(3, j) - multiplied_tile(4, j)
    end do
    do i = 1, 2
      final_tile(i, 1) = tmp_tile_1(i, 1) + tmp_tile_1(i, 2) + tmp_tile_1(i, 3)
      final_tile(i, 2) = tmp_tile_1(i, 2) - tmp_tile_1(i, 3) - tmp_tile_1(i, 4)
    end do

    do i = 1, 2
      do j = 1, 2
        x = 2 * (tile_i + offset_i) + (i - 1)
        y = 2 * tile_j + (j - 1)
        if (x < out_map_size .and. y < out_map_size) then
          out_idx = x * out_map_size + y + 1
          output(out_idx) = final_tile(i, j)
        end if
      end do
    end do
  end subroutine compute_tile_host

  logical function compare_results(expected, actual) result(ok)
    real(real32), intent(in) :: expected(:), actual(:)
    integer :: idx
    ok = .true.
    do idx = 1, output_size
      if (percent_diff(expected(idx), actual(idx)) > percent_diff_error_threshold) then
        ok = .false.
      end if
    end do
  end function compare_results

  real(real32) function percent_diff(val1, val2) result(value)
    real(real32), intent(in) :: val1, val2
    if (abs(val1) < 0.01_real32 .and. abs(val2) < 0.01_real32) then
      value = 0.0_real32
    else
      value = 100.0_real32 * abs(abs(val1 - val2) / abs(val1 + small_float_val))
    end if
  end function percent_diff

end program main
