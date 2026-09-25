! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_long, c_signed_char
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: instances = 224, attributes = 4096, threads = 128
  integer :: iterations, i, attr, instance_id, status, value
  integer(int32), allocatable :: data(:), cpu_distance(:), gpu_distance(:)
  integer(c_signed_char), allocatable :: data_char(:)
  real(real64) :: start_time, end_time, elapsed_us

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_random() bind(C, name='random') result(value)
      import :: c_long
      integer(c_long) :: value
    end function c_random
  end interface

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <iterations>'
    stop 1
  end if
  iterations = read_arg(1)

  allocate(data(instances * attributes), data_char(instances * attributes), &
           cpu_distance(instances * instances), gpu_distance(instances * instances))
  call c_srand(2_c_int)
  do attr = 1, attributes
    do instance_id = 1, instances
      value = int(mod(c_random(), 3_c_long))
      data(attr + attributes * (instance_id - 1)) = int(value, int32)
      data_char(attr + attributes * (instance_id - 1)) = int(value, c_signed_char)
    end do
  end do

  cpu_distance = 0_int32
  start_time = omp_get_wtime()
  call cpu_reference(data, cpu_distance)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'CPU time: ', (end_time - start_time) * 1.0e6_real64, ' (us)'

  status = 0
  !$omp target data map(to: data_char(1:instances*attributes)) map(alloc: gpu_distance(1:instances*instances))
  elapsed_us = 0.0_real64
  do i = 1, iterations
    gpu_distance = 0_int32
    !$omp target update to(gpu_distance(1:instances*instances))
    start_time = omp_get_wtime()
    call distance_register_device(data_char, gpu_distance)
    end_time = omp_get_wtime()
    elapsed_us = elapsed_us + (end_time - start_time) * 1.0e6_real64
  end do
  !$omp target update from(gpu_distance(1:instances*instances))
  print '(A,F0.6,A)', 'Average kernel execution time (w/o shared memory): ', elapsed_us / real(iterations, real64), ' (us)'
  status = print_status(cpu_distance, gpu_distance)

  elapsed_us = 0.0_real64
  do i = 1, iterations
    gpu_distance = 0_int32
    !$omp target update to(gpu_distance(1:instances*instances))
    start_time = omp_get_wtime()
    call distance_shared_device(data_char, gpu_distance)
    end_time = omp_get_wtime()
    elapsed_us = elapsed_us + (end_time - start_time) * 1.0e6_real64
  end do
  !$omp target update from(gpu_distance(1:instances*instances))
  print '(A,F0.6,A)', 'Average kernel execution time (w/ shared memory): ', elapsed_us / real(iterations, real64), ' (us)'
  if (print_status(cpu_distance, gpu_distance) /= 0) status = 1
  !$omp end target data

  deallocate(data_char, data, cpu_distance, gpu_distance)
  stop status

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine cpu_reference(data, distance)
    integer(int32), intent(in) :: data(:)
    integer(int32), intent(inout) :: distance(:)
    integer :: i, j, k, count
    !$omp parallel do collapse(2) private(k, count)
    do j = 1, instances
      do i = 1, instances
        count = 0
        do k = 1, attributes
          if (data((i - 1) * attributes + k) /= data((j - 1) * attributes + k)) count = count + 1
        end do
        distance(i + instances * (j - 1)) = count
      end do
    end do
    !$omp end parallel do
  end subroutine cpu_reference

  subroutine distance_register_device(data_char, distance)
    integer(c_signed_char), intent(in) :: data_char(:)
    integer(int32), intent(inout) :: distance(:)
    integer :: idx, gx, gy, k, count
    integer(c_signed_char) :: jx, jy, jz, jw, kx, ky, kz, kw
    !$omp target teams num_teams(instances*instances) thread_limit(threads) private(idx, gx, gy, k, count, jx, jy, jz, jw, kx, ky, kz, kw)
    !$omp parallel private(idx, gx, gy, k, count, jx, jy, jz, jw, kx, ky, kz, kw)
    idx = omp_get_thread_num()
    gx = mod(omp_get_team_num(), instances)
    gy = omp_get_team_num() / instances

    do k = 4 * idx, attributes - 1, threads * 4
      jx = data_char(k + attributes * gx + 1)
      jy = data_char(k + attributes * gx + 2)
      jz = data_char(k + attributes * gx + 3)
      jw = data_char(k + attributes * gx + 4)
      kx = data_char(k + attributes * gy + 1)
      ky = data_char(k + attributes * gy + 2)
      kz = data_char(k + attributes * gy + 3)
      kw = data_char(k + attributes * gy + 4)

      count = 0
      if (jx /= kx) count = count + 1
      if (jy /= ky) count = count + 1
      if (jz /= kz) count = count + 1
      if (jw /= kw) count = count + 1

      !$omp atomic update
      distance(instances * gx + gy + 1) = distance(instances * gx + gy + 1) + int(count, int32)
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine distance_register_device

  subroutine distance_shared_device(data_char, distance)
    integer(c_signed_char), intent(in) :: data_char(:)
    integer(int32), intent(inout) :: distance(:)
    integer(int32) :: dist(threads)
    integer :: idx, gx, gy, k, count, stride
    integer(c_signed_char) :: jx, jy, jz, jw, kx, ky, kz, kw
    !$omp target teams num_teams(instances*instances) thread_limit(threads) private(dist, idx, gx, gy, k, count, stride, jx, jy, jz, jw, kx, ky, kz, kw)
    !$omp parallel private(idx, gx, gy, k, count, stride, jx, jy, jz, jw, kx, ky, kz, kw) shared(dist)
    idx = omp_get_thread_num()
    gx = mod(omp_get_team_num(), instances)
    gy = omp_get_team_num() / instances

    dist(idx + 1) = 0_int32
    !$omp barrier

    do k = 4 * idx, attributes - 1, threads * 4
      jx = data_char(k + attributes * gx + 1)
      jy = data_char(k + attributes * gx + 2)
      jz = data_char(k + attributes * gx + 3)
      jw = data_char(k + attributes * gx + 4)
      kx = data_char(k + attributes * gy + 1)
      ky = data_char(k + attributes * gy + 2)
      kz = data_char(k + attributes * gy + 3)
      kw = data_char(k + attributes * gy + 4)

      count = 0
      if (jx /= kx) count = count + 1
      if (jy /= ky) count = count + 1
      if (jz /= kz) count = count + 1
      if (jw /= kw) count = count + 1
      dist(idx + 1) = dist(idx + 1) + int(count, int32)
    end do

    !$omp barrier
    stride = threads / 2
    do while (stride > 0)
      if (idx < stride) dist(idx + 1) = dist(idx + 1) + dist(idx + stride + 1)
      !$omp barrier
      stride = stride / 2
    end do

    if (idx == 0) distance(instances * gy + gx + 1) = dist(1)
    !$omp end parallel
    !$omp end target teams
  end subroutine distance_shared_device

  integer function print_status(expected, actual)
    integer(int32), intent(in) :: expected(:), actual(:)
    integer :: i
    logical :: ok
    ok = .true.
    do i = 1, instances * instances
      if (expected(i) /= actual(i)) then
        ok = .false.
        exit
      end if
    end do
    if (ok) then
      print '(A)', 'PASS'
      print_status = 0
    else
      print '(A)', 'FAIL'
      print_status = 1
    end if
  end function print_status

end program main
