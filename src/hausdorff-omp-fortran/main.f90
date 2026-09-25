! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
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

  real(real32), parameter :: rand_max = 2147483647.0_real32
  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: num_a, num_b, repeat, i
  real(real32), allocatable :: h_Apoints(:, :), h_Bpoints(:, :)
  real(real32) :: distance(2), ref_distance, test_distance
  real(real64) :: start_time, end_time, total_time_ms
  logical :: error

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <number of points in space A> <number of points in space B> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) num_a
  read(arg2, *) num_b
  read(arg3, *) repeat

  allocate(h_Apoints(2, num_a), h_Bpoints(2, num_b))

  call c_srand(123_c_int)
  do i = 1, num_a
    h_Apoints(1, i) = real(c_rand(), real32) / rand_max
    h_Apoints(2, i) = real(c_rand(), real32) / rand_max
  end do

  do i = 1, num_b
    h_Bpoints(1, i) = real(c_rand(), real32) / rand_max
    h_Bpoints(2, i) = real(c_rand(), real32) / rand_max
  end do

  distance = -1.0_real32
  total_time_ms = 0.0_real64

  !$omp target data map(to: h_Apoints(1:2,1:num_a), h_Bpoints(1:2,1:num_b)) &
  !$omp& map(tofrom: distance(1:2))
  do i = 1, repeat
    distance = -1.0_real32
    !$omp target update to(distance(1:2))

    start_time = omp_get_wtime()
    call compute_distance(h_Apoints, h_Bpoints, distance(1), num_a, num_b)
    call compute_distance(h_Bpoints, h_Apoints, distance(2), num_b, num_a)
    end_time = omp_get_wtime()

    total_time_ms = total_time_ms + (end_time - start_time) * 1.0e3_real64
  end do
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average execution time of kernels: ', total_time_ms / real(repeat, real64), ' (ms)'

  write(*,'(A)') 'Verifying the result may take a while..'
  ref_distance = hausdorff_distance(h_Apoints, h_Bpoints, num_a, num_b)
  test_distance = max(distance(1), distance(2))

  error = abs(test_distance - ref_distance) > 1.0e-3_real32
  if (error) then
    write(*,'(A)') 'FAIL'
  else
    write(*,'(A)') 'PASS'
  end if

  deallocate(h_Apoints, h_Bpoints)

contains

  pure real(real32) function hd(ap_x, ap_y, bp_x, bp_y) result(value)
    real(real32), intent(in) :: ap_x, ap_y, bp_x, bp_y
    value = (ap_x - bp_x) * (ap_x - bp_x) + (ap_y - bp_y) * (ap_y - bp_y)
  end function hd

  subroutine compute_distance(Apoints, Bpoints, distance, num_a, num_b)
    real(real32), intent(in) :: Apoints(:, :), Bpoints(:, :)
    real(real32), intent(inout) :: distance
    integer, intent(in) :: num_a, num_b
    integer :: i, j
    real(real32) :: d, t, p_x, p_y

    !$omp target teams distribute parallel do reduction(max:distance) thread_limit(256) &
    !$omp& private(i, j, d, t)
    do i = 1, num_a
      d = huge(1.0_real32)
      p_x = Apoints(1, i)
      p_y = Apoints(2, i)
      do j = 1, num_b
        t = hd(p_x, p_y, Bpoints(1, j), Bpoints(2, j))
        if (t < d) d = t
      end do
      if (d > distance) distance = d
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_distance

  real(real32) function directed_distance(Apoints, Bpoints, num_a, num_b) result(dis)
    real(real32), intent(in) :: Apoints(:, :), Bpoints(:, :)
    integer, intent(in) :: num_a, num_b
    integer :: i, j
    real(real32) :: d, t, p_x, p_y

    dis = -1.0_real32
    do i = 1, num_a
      d = huge(1.0_real32)
      p_x = Apoints(1, i)
      p_y = Apoints(2, i)
      do j = 1, num_b
        t = hd(p_x, p_y, Bpoints(1, j), Bpoints(2, j))
        if (t < d) d = t
      end do
      if (d > dis) dis = d
    end do
  end function directed_distance

  real(real32) function hausdorff_distance(Apoints, Bpoints, num_a, num_b) result(dis)
    real(real32), intent(in) :: Apoints(:, :), Bpoints(:, :)
    integer, intent(in) :: num_a, num_b
    real(real32) :: h_ab, h_ba

    h_ab = directed_distance(Apoints, Bpoints, num_a, num_b)
    h_ba = directed_distance(Bpoints, Apoints, num_b, num_a)
    dis = max(h_ab, h_ba)
  end function hausdorff_distance

end program main
