! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  interface
    subroutine s8n_fill_input_cpp(values, input_size) bind(C, name="s8n_fill_input_cpp")
      import :: c_int
      integer(c_int) :: values(*)
      integer(c_int), value :: input_size
    end subroutine s8n_fill_input_cpp
  end interface

  integer :: b, n, repeat, input_size, output_size, radius
  integer :: i, error_count
  integer, allocatable :: h_xyz(:), h_out(:), h_out2(:), h_out4(:)
  integer, allocatable :: r_out(:), r_out2(:), r_out4(:)
  real(real64) :: start_time, end_time, elapsed_us

  if (command_argument_count() /= 3) then
    print '(A,A)', 'Usage: ', './main <number of batches> <number of points> <repeat>'
    stop 1
  end if

  b = read_arg(1)
  n = read_arg(2)
  repeat = read_arg(3)
  input_size = b * n * 3
  output_size = b * n * 8
  radius = 512
  error_count = 0

  allocate(h_xyz(input_size))
  allocate(h_out(output_size), r_out(output_size))
  allocate(h_out2(output_size * 2), r_out2(output_size * 2))
  allocate(h_out4(output_size * 4), r_out4(output_size * 4))

  call fill_input(h_xyz)

  !$omp target data map(to: h_xyz(1:input_size)) &
  !$omp& map(alloc: h_out(1:output_size), h_out2(1:output_size*2), h_out4(1:output_size*4))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call k_cube_select(b, n, radius, h_xyz, h_out)
  end do
  end_time = omp_get_wtime()
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time of select kernel: ', elapsed_us, ' (us)'

  !$omp target update from(h_out(1:output_size))
  call cube_select(b, n, radius, h_xyz, r_out)
  if (any(h_out /= r_out)) error_count = error_count + 1

  start_time = omp_get_wtime()
  do i = 1, repeat
    call k_cube_select_two(b, n, radius, h_xyz, h_out2)
  end do
  end_time = omp_get_wtime()
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time of select2 kernel: ', elapsed_us, ' (us)'

  !$omp target update from(h_out2(1:output_size*2))
  call cube_select_two(b, n, radius, h_xyz, r_out2)
  if (any(h_out2 /= r_out2)) error_count = error_count + 1

  start_time = omp_get_wtime()
  do i = 1, repeat
    call k_cube_select_four(b, n, radius, h_xyz, h_out4)
  end do
  end_time = omp_get_wtime()
  elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
  print '(A,F0.6,A)', 'Average execution time of select4 kernel: ', elapsed_us, ' (us)'

  !$omp target update from(h_out4(1:output_size*4))
  call cube_select_four(b, n, radius, h_xyz, r_out4)
  if (any(h_out4 /= r_out4)) error_count = error_count + 1
  !$omp end target data

  if (error_count == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(h_xyz, h_out, h_out2, h_out4, r_out, r_out2, r_out4)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_input(values)
    integer, intent(out) :: values(:)
    call s8n_fill_input_cpp(values, int(size(values), c_int))
  end subroutine fill_input

  subroutine k_cube_select(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    integer :: batch_idx, point_idx, other_idx, slot, base_xyz, base_out
    integer :: x, y, z, tx, ty, tz, dist, octant
    integer :: temp_dist(8)
    !$omp target teams distribute num_teams(b) private(base_xyz, base_out)
    do batch_idx = 0, b - 1
      base_xyz = batch_idx * n * 3
      base_out = batch_idx * n * 8
      !$omp parallel do num_threads(512) private(point_idx, other_idx, slot, x, y, z, tx, ty, tz, dist, octant, temp_dist)
      do point_idx = 0, n - 1
        x = xyz_all(base_xyz + point_idx * 3 + 1)
        y = xyz_all(base_xyz + point_idx * 3 + 2)
        z = xyz_all(base_xyz + point_idx * 3 + 3)
        do slot = 1, 8
          temp_dist(slot) = radius
          out_all(base_out + point_idx * 8 + slot) = point_idx
        end do
        do other_idx = 0, n - 1
          if (point_idx /= other_idx) cycle
          tx = xyz_all(base_xyz + other_idx * 3 + 1)
          ty = xyz_all(base_xyz + other_idx * 3 + 2)
          tz = xyz_all(base_xyz + other_idx * 3 + 3)
          dist = squared_distance(x, y, z, tx, ty, tz)
          if (dist > radius) cycle
          octant = merge(4, 0, tx > x) + merge(2, 0, ty > y) + merge(1, 0, tz > z) + 1
          if (dist < temp_dist(octant)) then
            out_all(base_out + point_idx * 8 + octant) = other_idx
            temp_dist(octant) = dist
          end if
        end do
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine k_cube_select

  subroutine k_cube_select_two(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    call select_multi(b, n, radius, xyz_all, out_all, 2, .true.)
  end subroutine k_cube_select_two

  subroutine k_cube_select_four(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    call select_multi(b, n, radius, xyz_all, out_all, 4, .true.)
  end subroutine k_cube_select_four

  subroutine cube_select(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    call select_single_host(b, n, radius, xyz_all, out_all)
  end subroutine cube_select

  subroutine cube_select_two(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    call select_multi(b, n, radius, xyz_all, out_all, 2, .false.)
  end subroutine cube_select_two

  subroutine cube_select_four(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    call select_multi(b, n, radius, xyz_all, out_all, 4, .false.)
  end subroutine cube_select_four

  subroutine select_single_host(b, n, radius, xyz_all, out_all)
    integer, intent(in) :: b, n, radius
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    integer :: batch_idx, point_idx, other_idx, slot, base_xyz, base_out
    integer :: x, y, z, tx, ty, tz, dist, octant
    integer :: temp_dist(8)
    do batch_idx = 0, b - 1
      base_xyz = batch_idx * n * 3
      base_out = batch_idx * n * 8
      do point_idx = 0, n - 1
        x = xyz_all(base_xyz + point_idx * 3 + 1)
        y = xyz_all(base_xyz + point_idx * 3 + 2)
        z = xyz_all(base_xyz + point_idx * 3 + 3)
        do slot = 1, 8
          temp_dist(slot) = radius
          out_all(base_out + point_idx * 8 + slot) = point_idx
        end do
        do other_idx = 0, n - 1
          if (point_idx == other_idx) cycle
          tx = xyz_all(base_xyz + other_idx * 3 + 1)
          ty = xyz_all(base_xyz + other_idx * 3 + 2)
          tz = xyz_all(base_xyz + other_idx * 3 + 3)
          dist = squared_distance(x, y, z, tx, ty, tz)
          if (dist > radius) cycle
          octant = merge(4, 0, tx > x) + merge(2, 0, ty > y) + merge(1, 0, tz > z) + 1
          if (dist < temp_dist(octant)) then
            out_all(base_out + point_idx * 8 + octant) = other_idx
            temp_dist(octant) = dist
          end if
        end do
      end do
    end do
  end subroutine select_single_host

  subroutine select_multi(b, n, radius, xyz_all, out_all, per_octant, on_device)
    integer, intent(in) :: b, n, radius, per_octant
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    logical, intent(in) :: on_device
    if (on_device) then
      call select_multi_device(b, n, radius, xyz_all, out_all, per_octant)
    else
      call select_multi_host(b, n, radius, xyz_all, out_all, per_octant)
    end if
  end subroutine select_multi

  subroutine select_multi_device(b, n, radius, xyz_all, out_all, per_octant)
    integer, intent(in) :: b, n, radius, per_octant
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    integer :: batch_idx, point_idx, other_idx, k, kk, base_xyz, base_out, stride
    integer :: x, y, z, tx, ty, tz, dist, base_slot
    integer :: temp_dist(32)
    logical :: flag
    stride = 8 * per_octant
    !$omp target teams distribute num_teams(b) private(base_xyz, base_out)
    do batch_idx = 0, b - 1
      base_xyz = batch_idx * n * 3
      base_out = batch_idx * n * stride
      !$omp parallel do num_threads(512) private(point_idx, other_idx, k, kk, x, y, z, tx, ty, tz, dist, base_slot, temp_dist, flag)
      do point_idx = 0, n - 1
        call init_slots(out_all, temp_dist, base_out, point_idx, stride, radius)
        x = xyz_all(base_xyz + point_idx * 3 + 1)
        y = xyz_all(base_xyz + point_idx * 3 + 2)
        z = xyz_all(base_xyz + point_idx * 3 + 3)
        do other_idx = 0, n - 1
          if (point_idx == other_idx) cycle
          tx = xyz_all(base_xyz + other_idx * 3 + 1)
          ty = xyz_all(base_xyz + other_idx * 3 + 2)
          tz = xyz_all(base_xyz + other_idx * 3 + 3)
          dist = squared_distance(x, y, z, tx, ty, tz)
          if (dist > radius) cycle
          base_slot = merge(4, 0, tx > x) + merge(2, 0, ty > y) + merge(1, 0, tz > z)
          base_slot = base_slot * per_octant + 1
          flag = .false.
          do k = 0, per_octant - 1
            if (dist < temp_dist(base_slot + k)) flag = .true.
            if (flag) then
              do kk = per_octant - 1, k + 1, -1
                out_all(base_out + point_idx * stride + base_slot + kk) = out_all(base_out + point_idx * stride + base_slot + kk - 1)
                temp_dist(base_slot + kk) = temp_dist(base_slot + kk - 1)
              end do
              out_all(base_out + point_idx * stride + base_slot + k) = other_idx
              temp_dist(base_slot + k) = dist
              exit
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine select_multi_device

  subroutine select_multi_host(b, n, radius, xyz_all, out_all, per_octant)
    integer, intent(in) :: b, n, radius, per_octant
    integer, intent(in) :: xyz_all(:)
    integer, intent(out) :: out_all(:)
    integer :: batch_idx, point_idx, other_idx, k, kk, base_xyz, base_out, stride
    integer :: x, y, z, tx, ty, tz, dist, base_slot
    integer :: temp_dist(32)
    logical :: flag
    stride = 8 * per_octant
    do batch_idx = 0, b - 1
      base_xyz = batch_idx * n * 3
      base_out = batch_idx * n * stride
      do point_idx = 0, n - 1
        call init_slots(out_all, temp_dist, base_out, point_idx, stride, radius)
        x = xyz_all(base_xyz + point_idx * 3 + 1)
        y = xyz_all(base_xyz + point_idx * 3 + 2)
        z = xyz_all(base_xyz + point_idx * 3 + 3)
        do other_idx = 0, n - 1
          if (point_idx == other_idx) cycle
          tx = xyz_all(base_xyz + other_idx * 3 + 1)
          ty = xyz_all(base_xyz + other_idx * 3 + 2)
          tz = xyz_all(base_xyz + other_idx * 3 + 3)
          dist = squared_distance(x, y, z, tx, ty, tz)
          if (dist > radius) cycle
          base_slot = merge(4, 0, tx > x) + merge(2, 0, ty > y) + merge(1, 0, tz > z)
          base_slot = base_slot * per_octant + 1
          flag = .false.
          do k = 0, per_octant - 1
            if (dist < temp_dist(base_slot + k)) flag = .true.
            if (flag) then
              do kk = per_octant - 1, k + 1, -1
                out_all(base_out + point_idx * stride + base_slot + kk) = out_all(base_out + point_idx * stride + base_slot + kk - 1)
                temp_dist(base_slot + kk) = temp_dist(base_slot + kk - 1)
              end do
              out_all(base_out + point_idx * stride + base_slot + k) = other_idx
              temp_dist(base_slot + k) = dist
              exit
            end if
          end do
        end do
      end do
    end do
  end subroutine select_multi_host

  subroutine init_slots(out_all, temp_dist, base_out, point_idx, stride, radius)
    integer, intent(inout) :: out_all(:)
    integer, intent(out) :: temp_dist(:)
    integer, intent(in) :: base_out, point_idx, stride, radius
    integer :: slot
    do slot = 1, stride
      temp_dist(slot) = radius
      out_all(base_out + point_idx * stride + slot) = point_idx
    end do
  end subroutine init_slots

  integer function squared_distance(x, y, z, tx, ty, tz)
    integer, intent(in) :: x, y, z, tx, ty, tz
    squared_distance = (x - tx) * (x - tx) + (y - ty) * (y - ty) + (z - tz) * (z - tz)
  end function squared_distance

end program main
