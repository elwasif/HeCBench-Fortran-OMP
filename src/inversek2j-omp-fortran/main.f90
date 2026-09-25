! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: max_loop = 25
  integer, parameter :: num_joints = 3
  integer, parameter :: num_joints_p1 = num_joints + 1
  integer, parameter :: block_size = 128
  real(real32), parameter :: pi = 3.14159265358979_real32
  real(real32), parameter :: tolerance = 1.0e-3_real32

  character(len=512) :: input_path
  integer :: data_size, iteration, error_count
  real(real32), allocatable :: x_target(:), y_target(:), angle_device(:), angle_cpu(:)

  if (command_argument_count() /= 2) then
    write(*, '(A)', advance='no') 'Usage: ./invkin <input file coefficients> <iterations>'
    print '(A)'
    stop 1
  end if

  call get_command_argument(1, input_path)
  iteration = read_arg(2)

  call read_coordinates(trim(input_path), x_target, y_target, data_size)
  allocate(angle_device(data_size * num_joints), angle_cpu(data_size * num_joints))
  angle_device = 0.0_real32
  angle_cpu = 0.0_real32

  print '(A,I0)', '# Data Size = ', data_size
  print '(A)', '# Coordinates are read from file...'

  call invkin_device(x_target, y_target, angle_device, data_size, iteration)

  call invkin_cpu(x_target, y_target, angle_cpu, data_size)
  error_count = count_angle_errors(angle_device, angle_cpu, data_size)
  if (error_count == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(x_target, y_target, angle_device, angle_cpu)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine read_coordinates(path, x_target, y_target, data_size)
    character(len=*), intent(in) :: path
    real(real32), allocatable, intent(out) :: x_target(:), y_target(:)
    integer, intent(out) :: data_size
    integer :: unit, i, ios
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open coordinate input'
    read(unit, *) data_size
    allocate(x_target(data_size), y_target(data_size))
    do i = 1, data_size
      read(unit, *) x_target(i), y_target(i)
    end do
    close(unit)
  end subroutine read_coordinates

  subroutine invkin_device(x_target, y_target, angles, data_size, iteration)
    real(real32), intent(in) :: x_target(:), y_target(:)
    real(real32), intent(out) :: angles(:)
    integer, intent(in) :: data_size, iteration
    integer :: n, idx, i, iter, curr_loop
    real(real32) :: angle_out(num_joints), x_data(num_joints_p1), y_data(num_joints_p1)
    real(real32) :: curr_x, curr_y, angle, direction, a_dot_b
    real(real32) :: pe_x, pe_y, pc_x, pc_y, diff_pe_pc_x, diff_pe_pc_y
    real(real32) :: diff_tgt_pc_x, diff_tgt_pc_y, len_diff_pe_pc, len_diff_tgt_pc
    real(real32) :: a_x, a_y, b_x, b_y
    real(real64) :: start_time, end_time, avg_us

    !$omp target data map(to: x_target(1:data_size), y_target(1:data_size)) map(from: angles(1:data_size*num_joints))
    start_time = omp_get_wtime()
    do n = 1, iteration
      !$omp target teams distribute parallel do simd thread_limit(block_size) &
      !$omp& private(i, iter, curr_loop, angle_out, x_data, y_data, curr_x, curr_y, angle, direction, a_dot_b) &
      !$omp& private(pe_x, pe_y, pc_x, pc_y, diff_pe_pc_x, diff_pe_pc_y, diff_tgt_pc_x, diff_tgt_pc_y) &
      !$omp& private(len_diff_pe_pc, len_diff_tgt_pc, a_x, a_y, b_x, b_y)
      do idx = 1, data_size
        curr_x = x_target(idx)
        curr_y = y_target(idx)
        do i = 1, num_joints
          angle_out(i) = 0.0_real32
        end do
        do i = 1, num_joints_p1
          x_data(i) = real(i - 1, real32)
          y_data(i) = 0.0_real32
        end do

        do curr_loop = 1, max_loop
          do iter = num_joints, 1, -1
            pe_x = x_data(num_joints_p1)
            pe_y = y_data(num_joints_p1)
            pc_x = x_data(iter)
            pc_y = y_data(iter)
            diff_pe_pc_x = pe_x - pc_x
            diff_pe_pc_y = pe_y - pc_y
            diff_tgt_pc_x = curr_x - pc_x
            diff_tgt_pc_y = curr_y - pc_y
            len_diff_pe_pc = sqrt(diff_pe_pc_x * diff_pe_pc_x + diff_pe_pc_y * diff_pe_pc_y)
            len_diff_tgt_pc = sqrt(diff_tgt_pc_x * diff_tgt_pc_x + diff_tgt_pc_y * diff_tgt_pc_y)
            a_x = diff_pe_pc_x / len_diff_pe_pc
            a_y = diff_pe_pc_y / len_diff_pe_pc
            b_x = diff_tgt_pc_x / len_diff_tgt_pc
            b_y = diff_tgt_pc_y / len_diff_tgt_pc
            a_dot_b = max(-1.0_real32, min(1.0_real32, a_x * b_x + a_y * b_y))
            angle = acos(a_dot_b) * (180.0_real32 / pi)
            direction = a_x * b_y - a_y * b_x
            if (direction < 0.0_real32) angle = -angle
            if (angle > 30.0_real32) then
              angle = 30.0_real32
            else if (angle < -30.0_real32) then
              angle = -30.0_real32
            end if
            angle_out(iter) = angle
            do i = 1, num_joints - 1
              angle_out(i + 1) = angle_out(i + 1) + angle_out(i)
            end do
          end do
        end do

        angles((idx - 1) * num_joints + 1) = angle_out(1)
        angles((idx - 1) * num_joints + 2) = angle_out(2)
        angles((idx - 1) * num_joints + 3) = angle_out(3)
      end do
      !$omp end target teams distribute parallel do simd
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0e6_real64) / real(iteration, real64)
    write(*, '(A,F0.6,A)') 'Average kernel execution time ', avg_us, ' (us)'
    !$omp end target data
  end subroutine invkin_device

  subroutine invkin_cpu(x_target, y_target, angles, data_size)
    real(real32), intent(in) :: x_target(:), y_target(:)
    real(real32), intent(out) :: angles(:)
    integer, intent(in) :: data_size
    integer :: idx, i, iter, curr_loop
    real(real32) :: angle_out(num_joints), x_data(num_joints_p1), y_data(num_joints_p1)
    real(real32) :: curr_x, curr_y, angle, direction, a_dot_b
    real(real32) :: pe_x, pe_y, pc_x, pc_y, diff_pe_pc_x, diff_pe_pc_y
    real(real32) :: diff_tgt_pc_x, diff_tgt_pc_y, len_diff_pe_pc, len_diff_tgt_pc
    real(real32) :: a_x, a_y, b_x, b_y

    do idx = 1, data_size
      curr_x = x_target(idx)
      curr_y = y_target(idx)
      angle_out = 0.0_real32
      do i = 1, num_joints_p1
        x_data(i) = real(i - 1, real32)
        y_data(i) = 0.0_real32
      end do
      do curr_loop = 1, max_loop
        do iter = num_joints, 1, -1
          pe_x = x_data(num_joints_p1)
          pe_y = y_data(num_joints_p1)
          pc_x = x_data(iter)
          pc_y = y_data(iter)
          diff_pe_pc_x = pe_x - pc_x
          diff_pe_pc_y = pe_y - pc_y
          diff_tgt_pc_x = curr_x - pc_x
          diff_tgt_pc_y = curr_y - pc_y
          len_diff_pe_pc = sqrt(diff_pe_pc_x * diff_pe_pc_x + diff_pe_pc_y * diff_pe_pc_y)
          len_diff_tgt_pc = sqrt(diff_tgt_pc_x * diff_tgt_pc_x + diff_tgt_pc_y * diff_tgt_pc_y)
          a_x = diff_pe_pc_x / len_diff_pe_pc
          a_y = diff_pe_pc_y / len_diff_pe_pc
          b_x = diff_tgt_pc_x / len_diff_tgt_pc
          b_y = diff_tgt_pc_y / len_diff_tgt_pc
          a_dot_b = max(-1.0_real32, min(1.0_real32, a_x * b_x + a_y * b_y))
          angle = acos(a_dot_b) * (180.0_real32 / pi)
          direction = a_x * b_y - a_y * b_x
          if (direction < 0.0_real32) angle = -angle
          if (angle > 30.0_real32) then
            angle = 30.0_real32
          else if (angle < -30.0_real32) then
            angle = -30.0_real32
          end if
          angle_out(iter) = angle
          do i = 1, num_joints - 1
            angle_out(i + 1) = angle_out(i + 1) + angle_out(i)
          end do
        end do
      end do
      angles((idx - 1) * num_joints + 1) = angle_out(1)
      angles((idx - 1) * num_joints + 2) = angle_out(2)
      angles((idx - 1) * num_joints + 3) = angle_out(3)
    end do
  end subroutine invkin_cpu

  integer function count_angle_errors(angle_device, angle_cpu, data_size)
    real(real32), intent(in) :: angle_device(:), angle_cpu(:)
    integer, intent(in) :: data_size
    integer :: idx, joint
    count_angle_errors = 0
    do idx = 1, data_size
      do joint = 1, num_joints
        if (abs(angle_device((idx - 1) * num_joints + joint) - angle_cpu((idx - 1) * num_joints + joint)) > tolerance) then
          count_angle_errors = count_angle_errors + 1
          exit
        end if
      end do
    end do
  end function count_angle_errors

end program main
