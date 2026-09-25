! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  real(c_double), parameter :: rand_max_plus_one = 2147483648.0_c_double
  character(len=256) :: arg0, arg
  integer :: length, iterations
  real(real32), parameter :: tolerance = 0.001_real32
  real(real32), allocatable :: diagonal(:), off_diagonal(:), eig0(:), eig1(:), ref0(:), ref1(:)
  integer, allocatable :: counts(:)
  integer :: active, ref_active, i
  real(real32) :: lower_limit, upper_limit
  real(real64) :: start_time, end_time
  logical :: ok

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <length of the diagonal of the square matrix> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) length
  call get_command_argument(2, arg)
  read(arg, *) iterations
  if (iterations <= 0) stop 1
  if (is_not_power_of_two(length)) length = round_to_power_of_two(length)
  if (length < 256) length = 256

  allocate(diagonal(0:length - 1), off_diagonal(0:length - 2), counts(0:length - 1))
  allocate(eig0(0:2 * length - 1), eig1(0:2 * length - 1), ref0(0:2 * length - 1), ref1(0:2 * length - 1))

  call fill_random(diagonal, length, 123)
  call fill_random(off_diagonal, length - 1, 133)
  call compute_gerschgorin(lower_limit, upper_limit, diagonal, off_diagonal, length)
  call initialize_intervals(eig0, lower_limit, upper_limit, length)
  eig1 = upper_limit

  !$omp target data map(to: diagonal(0:length - 1), off_diagonal(0:length - 2)) &
  !$omp& map(alloc: counts(0:length - 1), eig0(0:2 * length - 1), eig1(0:2 * length - 1))
  if (iterations /= 1) then
    do i = 1, 2
      call run_kernels(diagonal, off_diagonal, counts, eig0, eig1, length, tolerance, active)
    end do
  end if

  write(*,'(A,I0,A)') 'Executing kernel for ', iterations, ' iterations'
  write(*,'(A)') '-------------------------------------------'
  start_time = omp_get_wtime()
  do i = 1, iterations
    call run_kernels(diagonal, off_diagonal, counts, eig0, eig1, length, tolerance, active)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time ', &
    (end_time - start_time) * 1.0e6_real64 / real(iterations, real64), ' (us)'
  !$omp end target data

  call compute_gerschgorin(lower_limit, upper_limit, diagonal, off_diagonal, length)
  call initialize_intervals(ref0, lower_limit, upper_limit, length)
  ref1 = upper_limit
  ref_active = 0
  do while (is_complete_ref(select_interval(ref_active, ref0, ref1), length, tolerance))
    if (ref_active == 0) then
      call eigen_cpu_reference(diagonal, off_diagonal, length, ref0, ref1, tolerance)
    else
      call eigen_cpu_reference(diagonal, off_diagonal, length, ref1, ref0, tolerance)
    end if
    ref_active = 1 - ref_active
  end do

  if (active == 0 .and. ref_active == 0) then
    ok = compare_intervals(eig0, ref0, 2 * length)
  else if (active == 0) then
    ok = compare_intervals(eig0, ref1, 2 * length)
  else if (ref_active == 0) then
    ok = compare_intervals(eig1, ref0, 2 * length)
  else
    ok = compare_intervals(eig1, ref1, 2 * length)
  end if

  if (ok) then
    write(*,'(A/)') 'PASS'
  else
    write(*,'(A/)') 'FAIL'
  end if

  deallocate(diagonal, off_diagonal, counts, eig0, eig1, ref0, ref1)

contains

  subroutine run_kernels(diagonal, off_diagonal, counts, eig0, eig1, length, tolerance, active)
    real(real32), intent(in) :: diagonal(0:), off_diagonal(0:), tolerance
    integer, intent(inout) :: counts(0:)
    real(real32), intent(inout) :: eig0(0:), eig1(0:)
    integer, intent(in) :: length
    integer, intent(out) :: active

    !$omp target update to(eig0(0:2 * length - 1), eig1(0:2 * length - 1))
    active = 0
    do while (is_complete_ref(select_interval(active, eig0, eig1), length, tolerance))
      if (active == 0) then
        call count_kernel(counts, eig0, diagonal, off_diagonal, length)
        call recalc_kernel(eig1, eig0, counts, diagonal, off_diagonal, length, tolerance)
        active = 1
        !$omp target update from(eig1(0:2 * length - 1))
      else
        call count_kernel(counts, eig1, diagonal, off_diagonal, length)
        call recalc_kernel(eig0, eig1, counts, diagonal, off_diagonal, length, tolerance)
        active = 0
        !$omp target update from(eig0(0:2 * length - 1))
      end if
    end do
  end subroutine run_kernels

  subroutine count_kernel(counts, intervals, diagonal, off_diagonal, length)
    integer, intent(inout) :: counts(0:)
    real(real32), intent(in) :: intervals(0:), diagonal(0:), off_diagonal(0:)
    integer, intent(in) :: length
    integer :: gid, i, lower_count, upper_count
    real(real32) :: lower_limit, upper_limit, prev_diff, diff

    !$omp target teams distribute parallel do thread_limit(256) private(gid, i, lower_count, upper_count, lower_limit, upper_limit, prev_diff, diff)
    do gid = 0, length - 1
      lower_limit = intervals(2 * gid)
      upper_limit = intervals(2 * gid + 1)
      lower_count = 0
      prev_diff = diagonal(0) - lower_limit
      if (prev_diff < 0.0_real32) lower_count = lower_count + 1
      do i = 1, length - 1
        diff = (diagonal(i) - lower_limit) - ((off_diagonal(i - 1) * off_diagonal(i - 1)) / prev_diff)
        if (diff < 0.0_real32) lower_count = lower_count + 1
        prev_diff = diff
      end do
      upper_count = 0
      prev_diff = diagonal(0) - upper_limit
      if (prev_diff < 0.0_real32) upper_count = upper_count + 1
      do i = 1, length - 1
        diff = (diagonal(i) - upper_limit) - ((off_diagonal(i - 1) * off_diagonal(i - 1)) / prev_diff)
        if (diff < 0.0_real32) upper_count = upper_count + 1
        prev_diff = diff
      end do
      counts(gid) = upper_count - lower_count
    end do
    !$omp end target teams distribute parallel do
  end subroutine count_kernel

  subroutine recalc_kernel(new_intervals, intervals, counts, diagonal, off_diagonal, length, tolerance)
    real(real32), intent(inout) :: new_intervals(0:)
    real(real32), intent(in) :: intervals(0:), diagonal(0:), off_diagonal(0:), tolerance
    integer, intent(in) :: counts(0:), length
    integer :: gid, lower_id, upper_id, current_index, index, l_id, u_id, i
    integer :: n, base_count
    real(real32) :: mid_value, division_width, lower_bound, upper_bound, prev_diff, diff

    !$omp target teams distribute parallel do thread_limit(256) private(gid, lower_id, upper_id, current_index, index, l_id, u_id, i, n, base_count, mid_value, division_width, lower_bound, upper_bound, prev_diff, diff)
    do gid = 0, length - 1
      lower_id = 2 * gid
      upper_id = lower_id + 1
      current_index = gid
      index = 0
      do while (index < length - 1 .and. current_index >= counts(index))
        current_index = current_index - counts(index)
        index = index + 1
      end do
      l_id = 2 * index
      u_id = l_id + 1
      if (counts(index) == 1) then
        lower_bound = intervals(l_id)
        upper_bound = intervals(u_id)
        mid_value = (upper_bound + lower_bound) / 2.0_real32
        n = 0
        prev_diff = diagonal(0) - mid_value
        if (prev_diff < 0.0_real32) n = n + 1
        do i = 1, length - 1
          diff = (diagonal(i) - mid_value) - ((off_diagonal(i - 1) * off_diagonal(i - 1)) / prev_diff)
          if (diff < 0.0_real32) n = n + 1
          prev_diff = diff
        end do
        base_count = 0
        prev_diff = diagonal(0) - lower_bound
        if (prev_diff < 0.0_real32) base_count = base_count + 1
        do i = 1, length - 1
          diff = (diagonal(i) - lower_bound) - ((off_diagonal(i - 1) * off_diagonal(i - 1)) / prev_diff)
          if (diff < 0.0_real32) base_count = base_count + 1
          prev_diff = diff
        end do
        n = n - base_count
        if (upper_bound - lower_bound < tolerance) then
          new_intervals(lower_id) = lower_bound
          new_intervals(upper_id) = upper_bound
        else if (n == 0) then
          new_intervals(lower_id) = mid_value
          new_intervals(upper_id) = upper_bound
        else
          new_intervals(lower_id) = lower_bound
          new_intervals(upper_id) = mid_value
        end if
      else
        division_width = (intervals(u_id) - intervals(l_id)) / real(counts(index), real32)
        new_intervals(lower_id) = intervals(l_id) + division_width * real(current_index, real32)
        new_intervals(upper_id) = new_intervals(lower_id) + division_width
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine recalc_kernel

  subroutine eigen_cpu_reference(diagonal, off_diagonal, length, intervals, new_intervals, tolerance)
    real(real32), intent(in) :: diagonal(0:), off_diagonal(0:), intervals(0:), tolerance
    real(real32), intent(out) :: new_intervals(0:)
    integer, intent(in) :: length
    integer :: i, j, offset, new_lid, new_uid
    integer :: less_lower, less_upper, num_sub
    real(real32) :: avg_width, lower_bound, upper_bound, mid

    new_intervals = intervals
    offset = 0
    do i = 0, length - 1
      less_lower = sturm_count(diagonal, off_diagonal, length, intervals(2 * i))
      less_upper = sturm_count(diagonal, off_diagonal, length, intervals(2 * i + 1))
      num_sub = less_upper - less_lower
      if (num_sub > 1) then
        avg_width = (intervals(2 * i + 1) - intervals(2 * i)) / real(num_sub, real32)
        do j = 0, num_sub - 1
          new_lid = 2 * (offset + j)
          new_uid = new_lid + 1
          new_intervals(new_lid) = intervals(2 * i) + real(j, real32) * avg_width
          new_intervals(new_uid) = new_intervals(new_lid) + avg_width
        end do
      else if (num_sub == 1) then
        lower_bound = intervals(2 * i)
        upper_bound = intervals(2 * i + 1)
        mid = (lower_bound + upper_bound) / 2.0_real32
        new_lid = 2 * offset
        new_uid = new_lid + 1
        if (upper_bound - lower_bound < tolerance) then
          new_intervals(new_lid) = lower_bound
          new_intervals(new_uid) = upper_bound
        else if (sturm_count(diagonal, off_diagonal, length, mid) == less_upper) then
          new_intervals(new_lid) = lower_bound
          new_intervals(new_uid) = mid
        else
          new_intervals(new_lid) = mid
          new_intervals(new_uid) = upper_bound
        end if
      end if
      offset = offset + num_sub
    end do
  end subroutine eigen_cpu_reference

  integer function sturm_count(diagonal, off_diagonal, length, x)
    real(real32), intent(in) :: diagonal(0:), off_diagonal(0:), x
    integer, intent(in) :: length
    integer :: i
    real(real32) :: prev_diff, diff

    sturm_count = 0
    prev_diff = diagonal(0) - x
    if (prev_diff < 0.0_real32) sturm_count = sturm_count + 1
    do i = 1, length - 1
      diff = (diagonal(i) - x) - ((off_diagonal(i - 1) * off_diagonal(i - 1)) / prev_diff)
      if (diff < 0.0_real32) sturm_count = sturm_count + 1
      prev_diff = diff
    end do
  end function sturm_count

  subroutine compute_gerschgorin(lower_limit, upper_limit, diagonal, off_diagonal, length)
    real(real32), intent(out) :: lower_limit, upper_limit
    real(real32), intent(in) :: diagonal(0:), off_diagonal(0:)
    integer, intent(in) :: length
    integer :: i
    real(real32) :: r

    lower_limit = diagonal(0) - abs(off_diagonal(0))
    upper_limit = diagonal(0) + abs(off_diagonal(0))
    do i = 1, length - 2
      r = abs(off_diagonal(i - 1)) + abs(off_diagonal(i))
      lower_limit = min(lower_limit, diagonal(i) - r)
      upper_limit = max(upper_limit, diagonal(i) + r)
    end do
    lower_limit = min(lower_limit, diagonal(length - 1) - abs(off_diagonal(length - 2)))
    upper_limit = max(upper_limit, diagonal(length - 1) + abs(off_diagonal(length - 2)))
  end subroutine compute_gerschgorin

  subroutine initialize_intervals(intervals, lower_limit, upper_limit, length)
    real(real32), intent(out) :: intervals(0:)
    real(real32), intent(in) :: lower_limit, upper_limit
    integer, intent(in) :: length

    intervals = upper_limit
    intervals(0) = lower_limit
    intervals(1) = upper_limit
  end subroutine initialize_intervals

  logical function is_complete_ref(intervals, length, tolerance)
    real(real32), intent(in) :: intervals(0:), tolerance
    integer, intent(in) :: length
    integer :: i

    is_complete_ref = .false.
    do i = 0, length - 1
      if (intervals(2 * i + 1) - intervals(2 * i) >= tolerance) then
        is_complete_ref = .true.
        return
      end if
    end do
  end function is_complete_ref

  logical function compare_intervals(reference, data, n)
    real(real32), intent(in) :: reference(0:), data(0:)
    integer, intent(in) :: n
    integer :: i
    real(real32) :: error, ref_norm, diff

    error = 0.0_real32
    ref_norm = 0.0_real32
    do i = 1, n - 1
      diff = reference(i) - data(i)
      error = error + diff * diff
      ref_norm = ref_norm + reference(i) * reference(i)
    end do
    if (abs(ref_norm) < 1.0e-7_real32) then
      compare_intervals = .false.
    else
      compare_intervals = sqrt(error) / sqrt(ref_norm) < 1.0e-6_real32
    end if
  end function compare_intervals

  function select_interval(active, eig0, eig1) result(intervals)
    integer, intent(in) :: active
    real(real32), target, intent(in) :: eig0(0:), eig1(0:)
    real(real32), pointer :: intervals(:)

    if (active == 0) then
      intervals => eig0
    else
      intervals => eig1
    end if
  end function select_interval

  subroutine fill_random(values, n, seed)
    real(real32), intent(out) :: values(0:)
    integer, intent(in) :: n, seed
    integer :: i

    call c_srand(int(seed, c_int))
    do i = 0, n - 1
      values(i) = real(256.0_c_double * real(c_rand(), c_double) / rand_max_plus_one, real32)
    end do
  end subroutine fill_random

  logical function is_not_power_of_two(value)
    integer, intent(in) :: value
    is_not_power_of_two = value <= 0 .or. iand(value, -value) /= value
  end function is_not_power_of_two

  integer function round_to_power_of_two(value)
    integer, intent(in) :: value
    round_to_power_of_two = 1
    do while (round_to_power_of_two < value)
      round_to_power_of_two = round_to_power_of_two * 2
    end do
  end function round_to_power_of_two

end program main
