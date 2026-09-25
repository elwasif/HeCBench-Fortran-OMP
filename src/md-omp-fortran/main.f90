! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: max_neighbors = 128
  integer, parameter :: domain_edge = 20
  integer, parameter :: thread_limit_value = 256
  real(real32), parameter :: cutsq = 13.5_real32
  real(real32), parameter :: lj1 = 1.5_real32
  real(real32), parameter :: lj2 = 2.0_real32
  integer, parameter :: prob_sizes(0:3) = [12288, 24576, 36864, 73728]

  character(len=256) :: arg0, arg
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

  integer :: size_class, iteration, n_atom, i, total_pairs
  integer, allocatable :: neighbor_list(:)
  real(real32), allocatable :: pos_x(:), pos_y(:), pos_z(:)
  real(real32), allocatable :: force_x(:), force_y(:), force_z(:)
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)', advance='no') 'usage: ', trim(arg0), ' <class size> <iteration>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) size_class
  call get_command_argument(2, arg); read(arg, *) iteration
  if (size_class < 0 .or. size_class > 3 .or. iteration < 0) stop 1

  n_atom = prob_sizes(size_class)
  allocate(pos_x(n_atom), pos_y(n_atom), pos_z(n_atom))
  allocate(force_x(n_atom), force_y(n_atom), force_z(n_atom))
  allocate(neighbor_list(max_neighbors * n_atom))

  write(*,'(A)') 'Initializing test problem (this can take several minutes for large problems).'

  call c_srand(123_c_int)
  do i = 1, n_atom
    pos_x(i) = real(modulo(c_rand(), domain_edge), real32)
    pos_y(i) = real(modulo(c_rand(), domain_edge), real32)
    pos_z(i) = real(modulo(c_rand(), domain_edge), real32)
  end do

  write(*,'(A)') 'Finished.'
  total_pairs = build_neighbor_list(n_atom, pos_x, pos_y, pos_z, neighbor_list)
  write(*,'(I0,A,I0,A,F0.6,A)') total_pairs, ' of ', n_atom * max_neighbors, &
    ' pairs within cutoff distance = ', 100.0_real64 * real(total_pairs, real64) / &
    real(n_atom * max_neighbors, real64), ' %'

  !$omp target data map(to: pos_x(1:n_atom), pos_y(1:n_atom), pos_z(1:n_atom), &
  !$omp& neighbor_list(1:max_neighbors*n_atom)) map(from: force_x(1:n_atom), &
  !$omp& force_y(1:n_atom), force_z(1:n_atom))
  call md_kernel(pos_x, pos_y, pos_z, force_x, force_y, force_z, neighbor_list, n_atom)
  !$omp target update from(force_x(1:n_atom), force_y(1:n_atom), force_z(1:n_atom))

  write(*,'(A)') 'Performing Correctness Check (may take several minutes)'
  call check_results(force_x, force_y, force_z, pos_x, pos_y, pos_z, neighbor_list, n_atom)

  start_time = omp_get_wtime()
  do i = 1, iteration
    call md_kernel(pos_x, pos_y, pos_z, force_x, force_y, force_z, neighbor_list, n_atom)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  if (iteration > 0) then
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(iteration, real64), ' (s)'
  else
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', 0.0_real64, ' (s)'
  end if

  deallocate(pos_x, pos_y, pos_z, force_x, force_y, force_z, neighbor_list)

contains

  subroutine md_kernel(pos_x, pos_y, pos_z, force_x, force_y, force_z, neighbor_list, n_atom)
    integer, intent(in) :: n_atom
    integer, intent(in) :: neighbor_list(max_neighbors * n_atom)
    real(real32), intent(in) :: pos_x(n_atom), pos_y(n_atom), pos_z(n_atom)
    real(real32), intent(out) :: force_x(n_atom), force_y(n_atom), force_z(n_atom)
    integer :: idx, j, jidx
    real(real32) :: ipos_x, ipos_y, ipos_z, fx, fy, fz
    real(real32) :: delx, dely, delz, r2inv, r6inv, force_c

    !$omp target teams distribute parallel do private(j, jidx, ipos_x, ipos_y, ipos_z, &
    !$omp& fx, fy, fz, delx, dely, delz, r2inv, r6inv, force_c) thread_limit(thread_limit_value)
    do idx = 1, n_atom
      ipos_x = pos_x(idx)
      ipos_y = pos_y(idx)
      ipos_z = pos_z(idx)
      fx = 0.0_real32
      fy = 0.0_real32
      fz = 0.0_real32

      do j = 1, max_neighbors
        jidx = neighbor_list((j - 1) * n_atom + idx) + 1
        delx = ipos_x - pos_x(jidx)
        dely = ipos_y - pos_y(jidx)
        delz = ipos_z - pos_z(jidx)
        r2inv = delx * delx + dely * dely + delz * delz
        if (r2inv > 0.0_real32 .and. r2inv < cutsq) then
          r2inv = 1.0_real32 / r2inv
          r6inv = r2inv * r2inv * r2inv
          force_c = r2inv * r6inv * (lj1 * r6inv - lj2)
          fx = fx + delx * force_c
          fy = fy + dely * force_c
          fz = fz + delz * force_c
        end if
      end do

      force_x(idx) = fx
      force_y(idx) = fy
      force_z(idx) = fz
    end do
    !$omp end target teams distribute parallel do
  end subroutine md_kernel

  integer function build_neighbor_list(n_atom, pos_x, pos_y, pos_z, neighbor_list) result(total_pairs)
    integer, intent(in) :: n_atom
    real(real32), intent(in) :: pos_x(n_atom), pos_y(n_atom), pos_z(n_atom)
    integer, intent(out) :: neighbor_list(max_neighbors * n_atom)
    integer :: i, j, count
    real(real32) :: curr_dist(max_neighbors), dist_ij
    integer :: curr_list(max_neighbors)

    total_pairs = 0
    do i = 1, n_atom
      curr_dist = huge(1.0_real32)
      curr_list = 0
      do j = 1, n_atom
        if (i == j) cycle
        dist_ij = distance2(pos_x, pos_y, pos_z, i, j)
        call insert_in_order(curr_dist, curr_list, j - 1, dist_ij)
      end do
      count = populate_neighbor_list(curr_dist, curr_list, i - 1, n_atom, neighbor_list)
      total_pairs = total_pairs + count
    end do
  end function build_neighbor_list

  real(real32) function distance2(pos_x, pos_y, pos_z, i, j)
    real(real32), intent(in) :: pos_x(:), pos_y(:), pos_z(:)
    integer, intent(in) :: i, j
    real(real32) :: delx, dely, delz

    delx = pos_x(i) - pos_x(j)
    dely = pos_y(i) - pos_y(j)
    delz = pos_z(i) - pos_z(j)
    distance2 = delx * delx + dely * dely + delz * delz
  end function distance2

  subroutine insert_in_order(curr_dist, curr_list, atom_index, dist_ij)
    real(real32), intent(inout) :: curr_dist(max_neighbors)
    integer, intent(inout) :: curr_list(max_neighbors)
    integer, intent(in) :: atom_index
    real(real32), intent(in) :: dist_ij
    integer :: k

    if (dist_ij > curr_dist(max_neighbors)) return

    do k = 1, max_neighbors
      if (dist_ij < curr_dist(k)) then
        if (k < max_neighbors) then
          curr_dist(k + 1:max_neighbors) = curr_dist(k:max_neighbors - 1)
          curr_list(k + 1:max_neighbors) = curr_list(k:max_neighbors - 1)
        end if
        curr_dist(k) = dist_ij
        curr_list(k) = atom_index
        return
      end if
    end do
  end subroutine insert_in_order

  integer function populate_neighbor_list(curr_dist, curr_list, atom_index, n_atom, neighbor_list) result(valid_pairs)
    real(real32), intent(in) :: curr_dist(max_neighbors)
    integer, intent(in) :: curr_list(max_neighbors)
    integer, intent(in) :: atom_index, n_atom
    integer, intent(out) :: neighbor_list(max_neighbors * n_atom)
    integer :: idx

    valid_pairs = 0
    do idx = 1, max_neighbors
      neighbor_list((idx - 1) * n_atom + atom_index + 1) = curr_list(idx)
      if (curr_dist(idx) < cutsq) valid_pairs = valid_pairs + 1
    end do
  end function populate_neighbor_list

  subroutine check_results(force_x, force_y, force_z, pos_x, pos_y, pos_z, neighbor_list, n_atom)
    integer, intent(in) :: n_atom
    integer, intent(in) :: neighbor_list(max_neighbors * n_atom)
    real(real32), intent(in) :: force_x(n_atom), force_y(n_atom), force_z(n_atom)
    real(real32), intent(in) :: pos_x(n_atom), pos_y(n_atom), pos_z(n_atom)
    integer :: i, j, jidx
    real(real32) :: fx, fy, fz, delx, dely, delz, r2inv, r6inv, force_c
    real(real32) :: max_error

    max_error = 0.0_real32
    do i = 1, n_atom
      fx = 0.0_real32
      fy = 0.0_real32
      fz = 0.0_real32
      do j = 1, max_neighbors
        jidx = neighbor_list((j - 1) * n_atom + i) + 1
        delx = pos_x(i) - pos_x(jidx)
        dely = pos_y(i) - pos_y(jidx)
        delz = pos_z(i) - pos_z(jidx)
        r2inv = delx * delx + dely * dely + delz * delz
        if (r2inv > 0.0_real32 .and. r2inv < cutsq) then
          r2inv = 1.0_real32 / r2inv
          r6inv = r2inv * r2inv * r2inv
          force_c = r2inv * r6inv * (lj1 * r6inv - lj2)
          fx = fx + delx * force_c
          fy = fy + dely * force_c
          fz = fz + delz * force_c
        end if
      end do
      max_error = max(max_error, abs(fx - force_x(i)))
      max_error = max(max_error, abs(fy - force_y(i)))
      max_error = max(max_error, abs(fz - force_z(i)))
    end do
    write(*,'(A,ES12.5)') 'Max error between host and device: ', max_error
  end subroutine check_results

end program main
