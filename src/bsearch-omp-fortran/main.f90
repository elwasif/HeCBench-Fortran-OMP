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

  integer :: num_elem, repeat, a_size, z_size, n
  integer :: i
  real(real32), allocatable :: a(:), z(:)
  integer, allocatable :: r(:)

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage ./main <number of elements> <repeat>'
    stop 1
  end if

  num_elem = read_arg(1)
  repeat = read_arg(2)
  a_size = num_elem
  z_size = 2 * a_size
  n = a_size - 1

  allocate(a(a_size), z(z_size), r(z_size))
  do i = 1, a_size
    a(i) = real(i - 1, real32)
  end do
  call fill_queries(z, n)
  r = 0

  !$omp target data map(to: a(1:a_size), z(1:z_size)) map(from: r(1:z_size))
  call bs1(a, z, r, z_size, n, repeat)

  call bs2(a, z, r, z_size, n, repeat)

  call bs3(a, z, r, z_size, n, repeat)

  call bs4(a, z, r, z_size, n, repeat)
  !$omp end target data

  deallocate(a, z, r)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_queries(z, n)
    real(real32), intent(out) :: z(:)
    integer, intent(in) :: n
    integer :: idx
    call c_srand(2_c_int)
    do idx = 1, size(z)
      z(idx) = real(mod(c_rand(), int(n, c_int)), real32)
    end do
  end subroutine fill_queries

  subroutine bs1(a, z, r, z_size, n, repeat)
    real(real32), intent(in) :: a(:), z(:)
    integer, intent(out) :: r(:)
    integer, intent(in) :: z_size, n, repeat
    integer :: rep, idx, low, high, mid
    real(real32) :: value
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do rep = 1, repeat
      !$omp target teams distribute parallel do thread_limit(256) private(value, low, high, mid)
      do idx = 1, z_size
        value = z(idx)
        low = 0
        high = n
        do while (high - low > 1)
          mid = low + (high - low) / 2
          if (value < a(mid + 1)) then
            high = mid
          else
            low = mid
          end if
        end do
        r(idx) = low
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    call print_time('bs1', start_time, end_time, repeat)
  end subroutine bs1

  subroutine bs2(a, z, r, z_size, n, repeat)
    real(real32), intent(in) :: a(:), z(:)
    integer, intent(out) :: r(:)
    integer, intent(in) :: z_size, n, repeat
    integer :: rep, idx, nbits, k, search_idx, candidate
    real(real32) :: value
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do rep = 1, repeat
      !$omp target teams distribute parallel do thread_limit(256) private(nbits, k, search_idx, candidate, value)
      do idx = 1, z_size
        nbits = 0
        do while (ishft(n, -nbits) /= 0)
          nbits = nbits + 1
        end do
        k = ishft(1, nbits - 1)
        value = z(idx)
        if (a(k + 1) <= value) then
          search_idx = k
        else
          search_idx = 0
        end if
        k = ishft(k, -1)
        do while (k /= 0)
          candidate = ior(search_idx, k)
          if (candidate < n .and. value >= a(candidate + 1)) search_idx = candidate
          k = ishft(k, -1)
        end do
        r(idx) = search_idx
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    call print_time('bs2', start_time, end_time, repeat)
  end subroutine bs2

  subroutine bs3(a, z, r, z_size, n, repeat)
    real(real32), intent(in) :: a(:), z(:)
    integer, intent(out) :: r(:)
    integer, intent(in) :: z_size, n, repeat
    integer :: rep, idx, nbits, k, search_idx, candidate, bounded
    real(real32) :: value
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do rep = 1, repeat
      !$omp target teams distribute parallel do thread_limit(256) private(nbits, k, search_idx, candidate, bounded, value)
      do idx = 1, z_size
        nbits = 0
        do while (ishft(n, -nbits) /= 0)
          nbits = nbits + 1
        end do
        k = ishft(1, nbits - 1)
        value = z(idx)
        if (a(k + 1) <= value) then
          search_idx = k
        else
          search_idx = 0
        end if
        k = ishft(k, -1)
        do while (k /= 0)
          candidate = ior(search_idx, k)
          bounded = min(candidate, n)
          if (value >= a(bounded + 1)) search_idx = candidate
          k = ishft(k, -1)
        end do
        r(idx) = search_idx
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    call print_time('bs3', start_time, end_time, repeat)
  end subroutine bs3

  subroutine bs4(a, z, r, z_size, n, repeat)
    real(real32), intent(in) :: a(:), z(:)
    integer, intent(out) :: r(:)
    integer, intent(in) :: z_size, n, repeat
    integer :: rep, nbits, k, lid, gid, p, search_idx, candidate, bounded
    real(real32) :: value
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do rep = 1, repeat
      !$omp target teams num_teams(z_size / 256) thread_limit(256) private(k)
      !$omp parallel private(nbits, lid, gid, p, search_idx, candidate, bounded, value)
      lid = omp_get_thread_num()
      gid = omp_get_team_num() * omp_get_num_threads() + lid
      if (lid == 0) then
        nbits = 0
        do while (ishft(n, -nbits) /= 0)
          nbits = nbits + 1
        end do
        k = ishft(1, nbits - 1)
      end if
      !$omp barrier

      p = k
      value = z(gid + 1)
      if (a(p + 1) <= value) then
        search_idx = p
      else
        search_idx = 0
      end if
      p = ishft(p, -1)
      do while (p /= 0)
        candidate = ior(search_idx, p)
        bounded = min(candidate, n)
        if (value >= a(bounded + 1)) search_idx = candidate
        p = ishft(p, -1)
      end do
      r(gid + 1) = search_idx
      !$omp end parallel
      !$omp end target teams
    end do
    end_time = omp_get_wtime()
    call print_time('bs4', start_time, end_time, repeat)
  end subroutine bs4

  subroutine verify(a, z, r, a_size, z_size, label)
    real(real32), intent(in) :: a(:), z(:)
    integer, intent(inout) :: r(:)
    integer, intent(in) :: a_size, z_size
    character(len=*), intent(in) :: label
    integer :: idx
    do idx = 1, z_size
      if (.not. (r(idx) + 1 < a_size .and. a(r(idx) + 1) <= z(idx) .and. z(idx) < a(r(idx) + 2))) then
        print '(A,A,A,I0,A,I0)', trim(label), ': incorrect result: index = ', '', idx - 1, ' r[index] = ', r(idx)
        exit
      end if
      r(idx) = -1
    end do
  end subroutine verify

  subroutine print_time(label, start_time, end_time, repeat)
    character(len=*), intent(in) :: label
    real(real64), intent(in) :: start_time, end_time
    integer, intent(in) :: repeat
    print '(A,A,A,ES12.6E2,A)', 'Average device execution time (', trim(label), ') ', &
        (end_time - start_time) / real(repeat, real64), ' (s)'
  end subroutine print_time

end program main
