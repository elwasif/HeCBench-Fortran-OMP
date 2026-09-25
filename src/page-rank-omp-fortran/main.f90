! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : error_unit, int32, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  real(real32), parameter :: d_factor = 0.85_real32
  integer, parameter :: default_iter = 1000
  real(real32), parameter :: default_threshold = 1.0e-16_real32

  integer :: n, iter, divisor, t, block_size
  real(real32) :: thresh, max_diff, max_diff_ref
  real(real64) :: ktime, start_time, end_time
  integer(int32), allocatable :: pages(:), noutlinks(:)
  real(real32), allocatable :: maps(:), ranks_gpu(:), ranks_ref(:), diffs(:)
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

  n = 1000
  iter = default_iter
  thresh = default_threshold
  divisor = 2
  call parse_args(n, iter, thresh, divisor)
  if (n <= 1 .or. iter <= 0 .or. divisor <= 0) stop 1

  allocate(pages(n * n), noutlinks(n), maps(n * n), ranks_gpu(n), ranks_ref(n), diffs(n))
  call random_pages(n, pages, noutlinks, divisor)
  ranks_gpu = 1.0_real32 / real(n, real32)
  ranks_ref = ranks_gpu
  diffs = 0.0_real32
  max_diff = 99.0_real32
  max_diff_ref = 99.0_real32
  ktime = 0.0_real64
  block_size = min(n, 256)

  !$omp target data map(to: pages(1:n*n), noutlinks(1:n), ranks_gpu(1:n), diffs(1:n)) map(alloc: maps(1:n*n))
  do t = 1, iter
    if (max_diff < thresh) exit
    start_time = omp_get_wtime()
    call map_device(pages, ranks_gpu, maps, noutlinks, n, block_size)
    call reduce_device(ranks_gpu, maps, diffs, n, block_size)
    end_time = omp_get_wtime()
    ktime = ktime + (end_time - start_time)
    !$omp target update from(diffs(1:n))
    max_diff = maximum_dif(diffs, n)
  end do
  !$omp end target data

  write(error_unit, '(A,F0.6,A,I0)') 'Max difference ', max_diff, ' is reached at iteration ', t

  diffs = 0.0_real32
  do t = 1, iter
    if (max_diff_ref < thresh) exit
    call map_reference(pages, ranks_ref, maps, noutlinks, n)
    call reduce_reference(ranks_ref, maps, diffs, n)
    max_diff_ref = maximum_dif(diffs, n)
  end do

  print '(A,I0,A,I0,A,F8.6,A,F0.6,A)', '"Options": "-n ', n, ' -i ', iter, ' -t ', thresh, &
    '". Total kernel execution time: ', ktime, ' (s)'

  ok = abs(max_diff - max_diff_ref) < 1.0e-3_real32
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(pages, noutlinks, maps, ranks_gpu, ranks_ref, diffs)

contains

  subroutine parse_args(n, iter, thresh, divisor)
    integer, intent(inout) :: n, iter, divisor
    real(real32), intent(inout) :: thresh
    integer :: arg_count, pos
    character(len=256) :: key, value

    arg_count = command_argument_count()
    pos = 1
    do while (pos <= arg_count)
      call get_command_argument(pos, key)
      if (pos == arg_count) exit
      call get_command_argument(pos + 1, value)
      select case (trim(key))
      case ('-n')
        read(value, *) n
      case ('-i')
        read(value, *) iter
      case ('-t')
        read(value, *) thresh
      case ('-q')
        read(value, *) divisor
      case default
        call usage()
        stop 1
      end select
      pos = pos + 2
    end do
  end subroutine parse_args

  subroutine usage()
    print '(A)', 'Usage: ./main [-n number of pages] [-i max iterations] [-t threshold] [-q divisor for zero density]'
  end subroutine usage

  subroutine random_pages(n, pages, noutlinks, divisor)
    integer, intent(in) :: n, divisor
    integer(int32), intent(out) :: pages(:), noutlinks(:)
    integer :: i, j, k

    call c_srand(1_c_int)
    do i = 1, n
      noutlinks(i) = 0_int32
      do j = 1, n
        if (i /= j .and. mod(abs(c_rand()), divisor) == 0) then
          pages((i - 1) * n + j) = 1_int32
          noutlinks(i) = noutlinks(i) + 1_int32
        end if
      end do
      if (noutlinks(i) == 0_int32) then
        do
          k = mod(abs(c_rand()), n) + 1
          if (k /= i) exit
        end do
        pages((i - 1) * n + k) = 1_int32
        noutlinks(i) = 1_int32
      end if
    end do
  end subroutine random_pages

  real(real32) function maximum_dif(diffs, n)
    real(real32), intent(in) :: diffs(:)
    integer, intent(in) :: n
    integer :: i

    maximum_dif = 0.0_real32
    do i = 1, n
      maximum_dif = max(maximum_dif, diffs(i))
    end do
  end function maximum_dif

  subroutine map_device(pages, page_ranks, maps, noutlinks, n, block_size)
    integer(int32), intent(in) :: pages(:), noutlinks(:)
    real(real32), intent(in) :: page_ranks(:)
    real(real32), intent(out) :: maps(:)
    integer, intent(in) :: n, block_size
    integer :: i, j
    real(real32) :: outbound_rank

    !$omp target teams distribute parallel do thread_limit(block_size) private(j, outbound_rank)
    do i = 1, n
      outbound_rank = page_ranks(i) / real(noutlinks(i), real32)
      do j = 1, n
        maps((i - 1) * n + j) = real(pages((i - 1) * n + j), real32) * outbound_rank
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine map_device

  subroutine reduce_device(page_ranks, maps, diffs, n, block_size)
    real(real32), intent(inout) :: page_ranks(:), diffs(:)
    real(real32), intent(in) :: maps(:)
    integer, intent(in) :: n, block_size
    integer :: i, j
    real(real32) :: old_rank, new_rank

    !$omp target teams distribute parallel do thread_limit(block_size) private(i, old_rank, new_rank)
    do j = 1, n
      old_rank = page_ranks(j)
      new_rank = 0.0_real32
      do i = 1, n
        new_rank = new_rank + maps((i - 1) * n + j)
      end do
      new_rank = ((1.0_real32 - d_factor) / real(n, real32)) + (d_factor * new_rank)
      diffs(j) = max(abs(new_rank - old_rank), diffs(j))
      page_ranks(j) = new_rank
    end do
    !$omp end target teams distribute parallel do
  end subroutine reduce_device

  subroutine map_reference(pages, page_ranks, maps, noutlinks, n)
    integer(int32), intent(in) :: pages(:), noutlinks(:)
    real(real32), intent(in) :: page_ranks(:)
    real(real32), intent(out) :: maps(:)
    integer, intent(in) :: n
    integer :: i, j
    real(real32) :: outbound_rank

    do i = 1, n
      outbound_rank = page_ranks(i) / real(noutlinks(i), real32)
      do j = 1, n
        maps((i - 1) * n + j) = real(pages((i - 1) * n + j), real32) * outbound_rank
      end do
    end do
  end subroutine map_reference

  subroutine reduce_reference(page_ranks, maps, diffs, n)
    real(real32), intent(inout) :: page_ranks(:), diffs(:)
    real(real32), intent(in) :: maps(:)
    integer, intent(in) :: n
    integer :: i, j
    real(real32) :: old_rank, new_rank

    do j = 1, n
      old_rank = page_ranks(j)
      new_rank = 0.0_real32
      do i = 1, n
        new_rank = new_rank + maps((i - 1) * n + j)
      end do
      new_rank = ((1.0_real32 - d_factor) / real(n, real32)) + (d_factor * new_rank)
      diffs(j) = max(abs(new_rank - old_rank), diffs(j))
      page_ranks(j) = new_rank
    end do
  end subroutine reduce_reference

end program main
