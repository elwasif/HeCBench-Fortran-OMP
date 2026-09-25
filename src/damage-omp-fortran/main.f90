! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: iso_c_binding, only : c_double, c_int64_t
  use omp_lib
  implicit none

  interface
    function damage_lcg_random_double(seed) bind(C, name='damage_lcg_random_double') result(value)
      import :: c_double, c_int64_t
      integer(c_int64_t), intent(inout) :: seed
      real(c_double) :: value
    end function damage_lcg_random_double
  end interface

  integer, parameter :: block_size = 256
  character(len=256) :: arg0, arg
  integer :: n, repeat, m, i, j, s
  integer(c_int64_t) :: seed
  integer, allocatable :: nlist(:), family(:), n_neigh(:)
  real(real64), allocatable :: damage(:)
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of points> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) n
  call get_command_argument(2, arg); read(arg, *) repeat
  if (n <= 0 .or. repeat <= 0) stop 1

  m = (n + block_size - 1) / block_size
  allocate(nlist(n), family(m), n_neigh(m), damage(m))

  seed = 123_c_int64_t
  do i = 1, n
    if (damage_lcg_random_double(seed) > 0.5_real64) then
      nlist(i) = 1
    else
      nlist(i) = -1
    end if
  end do

  do i = 1, m
    s = 0
    do j = (i - 1) * block_size + 1, min(i * block_size, n)
      if (nlist(j) /= -1) s = s + 1
    end do
    family(i) = int(real(s + 1, real64) + real(s, real64) * damage_lcg_random_double(seed))
  end do

  !$omp target data map(to: nlist(1:n), family(1:m)) map(from: n_neigh(1:m), damage(1:m))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call damage_of_node(n, m, nlist, family, n_neigh, damage)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data
  write(*,'(A,F8.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  call validate(n, m, nlist, family, n_neigh, damage)

  !$omp target data map(to: nlist(1:n), family(1:m)) map(from: n_neigh(1:m), damage(1:m))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call damage_of_node_optimized(m, n, nlist, family, n_neigh, damage)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data
  write(*,'(A,F8.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  call validate(n, m, nlist, family, n_neigh, damage)

  deallocate(nlist, family, n_neigh, damage)

contains

  subroutine damage_of_node(n, m, nlist, family, n_neigh, damage)
    integer, intent(in) :: n, m
    integer, intent(in) :: nlist(n), family(m)
    integer, intent(out) :: n_neigh(m)
    real(real64), intent(out) :: damage(m)
    integer :: local_cache(block_size)
    integer :: local_id, local_size, nid, global_id, step, neighbours

    !$omp target teams num_teams((n + block_size - 1) / block_size) thread_limit(block_size) &
    !$omp& private(local_cache)
    !$omp parallel shared(local_cache) private(local_id, local_size, nid, global_id, step, neighbours)
    local_id = omp_get_thread_num()
    local_size = block_size
    nid = omp_get_team_num()
    global_id = nid * local_size + local_id + 1

    if (global_id <= n) then
      if (nlist(global_id) /= -1) then
        local_cache(local_id + 1) = 1
      else
        local_cache(local_id + 1) = 0
      end if
    else
      local_cache(local_id + 1) = 0
    end if

    !$omp barrier

    step = local_size / 2
    do while (step > 0)
      if (local_id < step) then
        local_cache(local_id + 1) = local_cache(local_id + 1) + local_cache(local_id + step + 1)
      end if
      !$omp barrier
      step = step / 2
    end do

    if (local_id == 0 .and. nid < m) then
      neighbours = local_cache(1)
      n_neigh(nid + 1) = neighbours
      damage(nid + 1) = 1.0_real64 - real(neighbours, real64) / real(family(nid + 1), real64)
    end if
    !$omp end parallel
    !$omp end target teams
  end subroutine damage_of_node

  subroutine damage_of_node_optimized(m, n, nlist, family, n_neigh, damage)
    integer, intent(in) :: m, n
    integer, intent(in) :: nlist(n), family(m)
    integer, intent(out) :: n_neigh(m)
    real(real64), intent(out) :: damage(m)
    integer :: nid, idx, lower, upper, sum

    !$omp target teams distribute num_teams(m) private(lower, upper, sum, idx)
    do nid = 1, m
      lower = (nid - 1) * block_size + 1
      upper = min(nid * block_size, n)
      sum = 0

      !$omp parallel do reduction(+:sum) num_threads(block_size)
      do idx = lower, upper
        if (nlist(idx) /= -1) sum = sum + 1
      end do
      !$omp end parallel do

      n_neigh(nid) = sum
      damage(nid) = 1.0_real64 - real(sum, real64) / real(family(nid), real64)
    end do
    !$omp end target teams distribute
  end subroutine damage_of_node_optimized

  subroutine validate(n, m, nlist, family, n_neigh, damage)
    integer, intent(in) :: n, m
    integer, intent(in) :: nlist(n), family(m), n_neigh(m)
    real(real64), intent(in) :: damage(m)
    integer :: gid, j, lower, upper, sum
    real(real64) :: damage_ref
    logical :: ok

    ok = .true.
    do gid = 1, m
      lower = (gid - 1) * block_size + 1
      upper = min(gid * block_size, n)
      sum = 0
      do j = lower, upper
        if (nlist(j) /= -1) sum = sum + 1
      end do
      damage_ref = 1.0_real64 - real(sum, real64) / real(family(gid), real64)
      if (n_neigh(gid) /= sum .or. abs(damage(gid) - damage_ref) > 1.0e-6_real64) then
        ok = .false.
        exit
      end if
    end do

    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if
  end subroutine validate

end program main
