! SPDX-License-Identifier: CC0-1.0
module minmax_types
  use, intrinsic :: iso_fortran_env, only : int64, real32
  implicit none

  type :: vec_2d
    real(real32) :: x
    real(real32) :: y
  end type vec_2d

  type :: min_pair
    real(real32) :: val
    integer(int64) :: idx
  end type min_pair

  type :: max_pair
    real(real32) :: val
    integer(int64) :: idx
  end type max_pair

  type :: minmax_pair
    real(real32) :: min_val
    real(real32) :: max_val
    integer(int64) :: min_idx
    integer(int64) :: max_idx
  end type minmax_pair

contains

  pure function norm2(point) result(value)
    type(vec_2d), intent(in) :: point
    real(real32) :: value

    value = point%x * point%x + point%y * point%y
  end function norm2

end module minmax_types

program minmax
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use minmax_types
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

  integer :: repeat, n
  real(real32) :: box_size
  character(len=64) :: arg
  type(vec_2d), allocatable :: points(:)
  type(vec_2d) :: min_point(2), max_point(2), r_min_point, r_max_point
  real(real64) :: start_time, elapsed_us
  logical :: ok

  if (command_argument_count() /= 2) then
    write(*,'("Usage: ./main <bounding-box size> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) box_size
  call get_command_argument(2, arg)
  read(arg, *) repeat

  call generate_points(box_size, points, n)
  write(*,'("Total number of points: ",I0)') n

  !$omp target data map(to: points)
  start_time = omp_get_wtime()
  call run_separate(points, n, repeat, min_point(1), max_point(1))
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average execution time of omp:min() + omp:max(): ",F0.6," (us)")') elapsed_us

  start_time = omp_get_wtime()
  call run_combined(points, n, repeat, min_point(2), max_point(2))
  elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
  write(*,'("Average execution time of omp:minmax(): ",F0.6," (us)")') elapsed_us
  !$omp end target data

  call cpu_reference(points, n, r_min_point, r_max_point)
  ok = points_equal(min_point(1), r_min_point) .and. &
       points_equal(max_point(1), r_max_point) .and. &
       points_equal(min_point(2), r_min_point) .and. &
       points_equal(max_point(2), r_max_point)

  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine run_separate(points, count, repeat, min_out, max_out)
    type(vec_2d), intent(in) :: points(:)
    integer, intent(in) :: count, repeat
    type(vec_2d), intent(out) :: min_out, max_out
    integer :: r, i
    real(real32) :: min_val, max_val, val
    type(min_pair) :: min_result
    type(max_pair) :: max_result

    do r = 1, repeat
      min_val = huge(1.0_real32)
      min_result = min_pair(huge(1.0_real32), -1_int64)
      !$omp target teams distribute parallel do reduction(min:min_val)
      do i = 1, count
        val = norm2(points(i))
        min_val = min(min_val, val)
      end do
      !$omp end target teams distribute parallel do
      min_result%val = min_val
      min_result%idx = find_index(points, count, min_result%val)

      max_val = -huge(1.0_real32)
      max_result = max_pair(-huge(1.0_real32), -1_int64)
      !$omp target teams distribute parallel do reduction(max:max_val)
      do i = 1, count
        val = norm2(points(i))
        max_val = max(max_val, val)
      end do
      !$omp end target teams distribute parallel do
      max_result%val = max_val
      max_result%idx = find_index(points, count, max_result%val)
    end do

    min_out = points(min_result%idx)
    max_out = points(max_result%idx)
  end subroutine run_separate

  subroutine run_combined(points, count, repeat, min_out, max_out)
    type(vec_2d), intent(in) :: points(:)
    integer, intent(in) :: count, repeat
    type(vec_2d), intent(out) :: min_out, max_out
    integer :: r, i
    real(real32) :: min_val, max_val, val
    type(minmax_pair) :: result

    do r = 1, repeat
      min_val = huge(1.0_real32)
      max_val = -huge(1.0_real32)
      result = minmax_pair(huge(1.0_real32), -huge(1.0_real32), -1_int64, -1_int64)
      !$omp target teams distribute parallel do reduction(min:min_val) reduction(max:max_val)
      do i = 1, count
        val = norm2(points(i))
        min_val = min(min_val, val)
        max_val = max(max_val, val)
      end do
      !$omp end target teams distribute parallel do
      result%min_val = min_val
      result%max_val = max_val
      result%min_idx = find_index(points, count, result%min_val)
      result%max_idx = find_index(points, count, result%max_val)
    end do

    min_out = points(result%min_idx)
    max_out = points(result%max_idx)
  end subroutine run_combined

  subroutine cpu_reference(points, count, min_out, max_out)
    type(vec_2d), intent(in) :: points(:)
    integer, intent(in) :: count
    type(vec_2d), intent(out) :: min_out, max_out
    integer :: i
    integer :: ref_min, ref_max

    ref_min = 1
    ref_max = 1
    do i = 2, count
      if (norm2(points(i)) < norm2(points(ref_min))) ref_min = i
      if (norm2(points(i)) > norm2(points(ref_max))) ref_max = i
    end do
    min_out = points(ref_min)
    max_out = points(ref_max)
  end subroutine cpu_reference

  subroutine generate_points(size, points, total_points)
    real(real32), intent(in) :: size
    type(vec_2d), allocatable, intent(out) :: points(:)
    integer, intent(out) :: total_points
    integer, allocatable :: rect_counts(:)
    real(real32), allocatable :: tlx(:), tly(:), brx(:), bry(:)
    real(real32) :: cur_tlx, cur_tly, cur_brx, cur_bry, area_x, area_y, phi
    integer :: rect_capacity, rects, nrect_points, offset, i, j
    real(real32), parameter :: rand_max = 2147483647.0_real32

    phi = (1.0_real32 + sqrt(5.0_real32)) * 0.5_real32
    cur_tlx = 0.0_real32
    cur_tly = 0.0_real32
    cur_brx = size
    cur_bry = size
    area_x = cur_brx - cur_tlx
    area_y = cur_bry - cur_tly

    rect_capacity = 64
    allocate(rect_counts(rect_capacity), tlx(rect_capacity), tly(rect_capacity), brx(rect_capacity), bry(rect_capacity))
    rects = 0
    total_points = 0

    do while (area_x > 1.0_real32 .and. area_y > 1.0_real32)
      select case (mod(rects, 4))
      case (0)
        cur_brx = cur_tlx - (cur_tlx - cur_brx) / phi
      case (1)
        cur_bry = cur_tly - (cur_tly - cur_bry) / phi
      case (2)
        cur_tlx = cur_tlx + (cur_brx - cur_tlx) / phi
      case default
        cur_tly = cur_tly + (cur_bry - cur_tly) / phi
      end select

      area_x = cur_brx - cur_tlx
      area_y = cur_bry - cur_tly
      nrect_points = int(sqrt(area_x * area_y * 1000000.0_real32))

      rects = rects + 1
      if (rects > rect_capacity) stop 1
      rect_counts(rects) = nrect_points
      tlx(rects) = cur_tlx
      tly(rects) = cur_tly
      brx(rects) = cur_brx
      bry(rects) = cur_bry
      total_points = total_points + nrect_points
    end do

    allocate(points(total_points))
    call c_srand(123_c_int)
    offset = 0
    do i = 1, rects
      do j = 1, rect_counts(i)
        offset = offset + 1
        points(offset)%x = tlx(i) + (brx(i) - tlx(i)) * &
          (real(c_rand(), real32) / rand_max)
        points(offset)%y = tly(i) + (bry(i) - tly(i)) * &
          (real(c_rand(), real32) / rand_max)
      end do
    end do
  end subroutine generate_points

  pure function points_equal(lhs, rhs) result(equal)
    type(vec_2d), intent(in) :: lhs, rhs
    logical :: equal

    equal = lhs%x == rhs%x .and. lhs%y == rhs%y
  end function points_equal

  pure function find_index(points, count, value) result(idx)
    type(vec_2d), intent(in) :: points(:)
    integer, intent(in) :: count
    real(real32), intent(in) :: value
    integer(int64) :: idx
    integer :: i
    real(real32) :: best_diff, diff

    idx = 1_int64
    best_diff = abs(norm2(points(1)) - value)
    do i = 1, count
      diff = abs(norm2(points(i)) - value)
      if (diff < best_diff) then
        idx = int(i, int64)
        best_diff = diff
      end if
    end do
  end function find_index

end program minmax
