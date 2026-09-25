! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_char, c_int, c_null_char
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

    function c_atoi(str) bind(C, name="atoi") result(value)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: str(*)
      integer(c_int) :: value
    end function c_atoi
  end interface

  integer, parameter :: blocks = 256
  real(real32), parameter :: rand_max = 2147483647.0_real32
  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: m, n, repeat, out_len, grids
  real(real32), allocatable :: subject(:), lower_bound(:), upper_bound(:)
  real(real32), allocatable :: lb(:), lb_h(:), avgs(:), stds(:)
  integer :: i
  logical :: ok
  real(real64) :: start_time, end_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(3A)') 'Usage: ./', trim(arg0), ' <query length> <subject length> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  m = c_atoi(trim(arg1) // c_null_char)
  n = c_atoi(trim(arg2) // c_null_char)
  repeat = c_atoi(trim(arg3) // c_null_char)

  out_len = n - m + 1
  grids = (out_len + blocks - 1) / blocks

  write(*,'(A,I0)') 'Query length = ', m
  write(*,'(A,I0)') 'Subject length = ', n

  allocate(subject(n), lower_bound(n), upper_bound(n))
  allocate(lb(out_len), lb_h(out_len), avgs(out_len), stds(out_len))

  lower_bound = 0.0_real32
  upper_bound = 0.0_real32
  call c_srand(123_c_int)
  do i = 1, n
    subject(i) = real(c_rand(), real32) / rand_max
  end do
  do i = 1, out_len
    avgs(i) = real(c_rand(), real32) / rand_max
  end do
  do i = 1, out_len
    stds(i) = real(c_rand(), real32) / rand_max
  end do
  do i = 1, m
    upper_bound(i) = real(c_rand(), real32) / rand_max
  end do
  do i = 1, m
    lower_bound(i) = real(c_rand(), real32) / rand_max
  end do

  !$omp target data map(to: subject(1:n), avgs(1:out_len), stds(1:out_len), &
  !$omp& lower_bound(1:n), upper_bound(1:n)) map(from: lb(1:out_len))
  start_time = omp_get_wtime()

  do i = 1, repeat
    call keogh_kernel(subject, avgs, stds, lower_bound, upper_bound, lb, m, n, out_len, grids)
  end do

  end_time = omp_get_wtime()
  elapsed = (end_time - start_time) / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed, ' (s)'
  !$omp end target data

  call reference_keogh(subject, avgs, stds, lb_h, lower_bound, upper_bound, m, n)

  ok = .true.
  do i = 1, out_len
    if (abs(lb(i) - lb_h(i)) > 1.0e-2_real32) then
      write(*,'(I0,1X,F0.6,1X,F0.6)') i - 1, lb(i), lb_h(i)
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(lb, lb_h, avgs, stds, subject, lower_bound, upper_bound)

contains

  subroutine keogh_kernel(subject, avgs, stds, lower_bound, upper_bound, lb, m, n, out_len, grids)
    real(real32), intent(in) :: subject(:), avgs(:), stds(:), lower_bound(:), upper_bound(:)
    real(real32), intent(out) :: lb(:)
    integer, intent(in) :: m, n, out_len, grids
    integer :: idx, j
    real(real32) :: residues, avg, std, value, lower, upper

    !$omp target teams distribute num_teams(grids) thread_limit(blocks) &
    !$omp& private(idx, j, residues, avg, std, value, lower, upper)
    do idx = 1, out_len
      residues = 0.0_real32
      avg = avgs(idx)
      std = stds(idx)

      !$omp parallel do reduction(+:residues) private(j, value, lower, upper)
      do j = 1, m
        value = (subject(idx + j - 1) - avg) / std
        lower = value - lower_bound(j)
        upper = value - upper_bound(j)
        if (upper > 0.0_real32) residues = residues + upper * upper
        if (lower < 0.0_real32) residues = residues + lower * lower
      end do
      !$omp end parallel do

      lb(idx) = residues
    end do
    !$omp end target teams distribute
  end subroutine keogh_kernel

  subroutine reference_keogh(subject, avgs, stds, lb_keogh, zlower, zupper, m, n)
    real(real32), intent(in) :: subject(:), avgs(:), stds(:), zlower(:), zupper(:)
    real(real32), intent(out) :: lb_keogh(:)
    integer, intent(in) :: m, n
    integer :: idx, j
    real(real32) :: residues, avg, std, value, lower, upper

    do idx = 1, n - m + 1
      residues = 0.0_real32
      avg = avgs(idx)
      std = stds(idx)
      do j = 1, m
        value = (subject(idx + j - 1) - avg) / std
        lower = value - zlower(j)
        upper = value - zupper(j)
        if (upper > 0.0_real32) residues = residues + upper * upper
        if (lower < 0.0_real32) residues = residues + lower * lower
      end do
      lb_keogh(idx) = residues
    end do
  end subroutine reference_keogh

end program main
