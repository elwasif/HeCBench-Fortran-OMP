! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  real(real32), parameter :: a1 = 1.5_real32
  real(real32), parameter :: a2 = 2.0_real32
  real(real32), parameter :: a3 = 2.5_real32
  real(real32), parameter :: a4 = 3.0_real32
  real(real32), parameter :: a5 = 3.5_real32
  real(real32), parameter :: r_min = 0.0_real32
  real(real32), parameter :: r_max = 2.0_real32
  real(real32), parameter :: pi_value = 3.1415926_real32
  real(real32), parameter :: eps = 1.0_real32
  real(real32), parameter :: rand_max = 2147483647.0_real32

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

  character(len=256) :: arg0, arg
  integer :: numk, check, iterations, dnum, inum, i
  real(real32), parameter :: width = 2000.0_real32, height = 2000.0_real32
  real(real32), parameter :: area = width * height
  real(real32), allocatable :: dx(:), dy(:), dz(:), avg_dist(:), ix(:), iy(:), iz(:), h_iz(:)
  real(real64) :: start_time, elapsed
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <pts> <check> <iterations>'
    write(*,'(A)') 'pts: number of points (unit: 1K)'
    write(*,'(A)') 'check: enable verification when the value is 1'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) numk
  call get_command_argument(2, arg); read(arg, *) check
  call get_command_argument(3, arg); read(arg, *) iterations
  if (numk <= 0 .or. iterations <= 0) stop 1

  dnum = numk * 1024
  inum = dnum
  allocate(dx(dnum), dy(dnum), dz(dnum), avg_dist(dnum), ix(inum), iy(inum), iz(inum), h_iz(inum))

  call c_srand(123_c_int)
  do i = 1, dnum
    dx(i) = rand_float() * 1000.0_real32
    dy(i) = rand_float() * 1000.0_real32
    dz(i) = rand_float() * 1000.0_real32
  end do
  do i = 1, inum
    ix(i) = rand_float() * 1000.0_real32
    iy(i) = rand_float() * 1000.0_real32
    iz(i) = 0.0_real32
  end do
  do i = 1, dnum
    avg_dist(i) = rand_float() * 3.0_real32
  end do
  h_iz = 0.0_real32

  write(*,'(A,I0,A)') 'Size = : ', numk, ' K '
  write(*,'(A,I0)') 'dnum = : ', dnum
  write(*,'(A,I0)') 'inum = : ', inum

  if (check /= 0) then
    write(*,'(A)') 'Verification enabled'
    call reference_kernel(dx, dy, dz, dnum, ix, iy, h_iz, inum, area, avg_dist)
  else
    write(*,'(A)') 'Verification disabled'
  end if

  !$omp target data map(to: dx(1:dnum), dy(1:dnum), dz(1:dnum), ix(1:inum), iy(1:inum), &
  !$omp& avg_dist(1:dnum)) map(alloc: iz(1:inum))
  call aidw_kernel(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
  !$omp target update from(iz(1:inum))

  if (check /= 0) then
    ok = verify(iz, h_iz, inum, eps)
    write(*,'(A)') merge('PASS', 'FAIL', ok)
  end if

  call aidw_kernel_tiled(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
  !$omp target update from(iz(1:inum))
  if (check /= 0) then
    ok = verify(iz, h_iz, inum, eps)
    write(*,'(A)') merge('PASS', 'FAIL', ok)
  end if

  start_time = omp_get_wtime()
  do i = 1, iterations
    call aidw_kernel(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
  end do
  elapsed = omp_get_wtime() - start_time
  write(*,'(A,F0.6,A)') 'Average execution time of AIDW_Kernel       ', &
    elapsed / real(iterations, real64), ' (s)'

  start_time = omp_get_wtime()
  do i = 1, iterations
    call aidw_kernel_tiled(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
  end do
  elapsed = omp_get_wtime() - start_time
  write(*,'(A,F0.6,A)') 'Average execution time of AIDW_Kernel_Tiled ', &
    elapsed / real(iterations, real64), ' (s)'
  !$omp end target data

  deallocate(dx, dy, dz, avg_dist, ix, iy, iz, h_iz)

contains

  real(real32) function rand_float() result(value)
    value = real(c_rand(), real32) / rand_max
  end function rand_float

  subroutine aidw_kernel(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
    real(real32), intent(in) :: dx(:), dy(:), dz(:), ix(:), iy(:), avg_dist(:), area
    real(real32), intent(inout) :: iz(:)
    integer, intent(in) :: dnum, inum
    integer :: tid

    !$omp target teams distribute parallel do thread_limit(block_size) &
    !$omp& map(to: dx(1:dnum), dy(1:dnum), dz(1:dnum), ix(1:inum), iy(1:inum), &
    !$omp& avg_dist(1:inum)) map(tofrom: iz(1:inum))
    do tid = 1, inum
      iz(tid) = interpolate_point(tid, dx, dy, dz, dnum, ix, iy, area, avg_dist)
    end do
    !$omp end target teams distribute parallel do
  end subroutine aidw_kernel

  subroutine aidw_kernel_tiled(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
    real(real32), intent(in) :: dx(:), dy(:), dz(:), ix(:), iy(:), avg_dist(:), area
    real(real32), intent(inout) :: iz(:)
    integer, intent(in) :: dnum, inum
    integer :: lid, tid, part, m, e, num_threads
    real(real32) :: sdx(block_size), sdy(block_size), sdz(block_size)
    real(real32) :: dist, t, alpha, sum_up, sum_dn, six_s, siy_s
    real(real32) :: r_obs, r_exp, r_s0, u_r, six_t, siy_t

    !$omp target teams num_teams((inum + block_size - 1) / block_size) thread_limit(block_size) &
    !$omp& map(to: dx(1:dnum), dy(1:dnum), dz(1:dnum), ix(1:inum), iy(1:inum), &
    !$omp& avg_dist(1:inum)) map(tofrom: iz(1:inum))
    !$omp parallel private(lid, tid, part, m, e, num_threads, dist, t, alpha, sum_up, sum_dn, &
    !$omp& six_s, siy_s, r_obs, r_exp, r_s0, u_r, six_t, siy_t)
    lid = omp_get_thread_num()
    tid = omp_get_team_num() * block_size + lid + 1
    if (tid <= inum) then
      dist = 0.0_real32
      t = 0.0_real32
      alpha = 0.0_real32
      part = (dnum - 1) / block_size
      sum_up = 0.0_real32
      sum_dn = 0.0_real32

      r_obs = avg_dist(tid)
      r_exp = 1.0_real32 / (2.0_real32 * sqrt(real(dnum, real32) / area))
      r_s0 = r_obs / r_exp
      u_r = 0.0_real32
      if (r_s0 >= r_min) u_r = 0.5_real32 - 0.5_real32 * cos(pi_value / r_max * (r_s0 - r_min))
      if (r_s0 >= r_max) u_r = 1.0_real32

      if (u_r >= 0.0_real32 .and. u_r <= 0.1_real32) alpha = a1
      if (u_r > 0.1_real32 .and. u_r <= 0.3_real32) alpha = a1 * (1.0_real32 - 5.0_real32 * (u_r - 0.1_real32)) + a2 * 5.0_real32 * (u_r - 0.1_real32)
      if (u_r > 0.3_real32 .and. u_r <= 0.5_real32) alpha = a3 * 5.0_real32 * (u_r - 0.3_real32) + a1 * (1.0_real32 - 5.0_real32 * (u_r - 0.3_real32))
      if (u_r > 0.5_real32 .and. u_r <= 0.7_real32) alpha = a3 * (1.0_real32 - 5.0_real32 * (u_r - 0.5_real32)) + a4 * 5.0_real32 * (u_r - 0.5_real32)
      if (u_r > 0.7_real32 .and. u_r <= 0.9_real32) alpha = a5 * 5.0_real32 * (u_r - 0.7_real32) + a4 * (1.0_real32 - 5.0_real32 * (u_r - 0.7_real32))
      if (u_r > 0.9_real32 .and. u_r <= 1.0_real32) alpha = a5
      alpha = alpha * 0.5_real32

      six_t = ix(tid)
      siy_t = iy(tid)
      do m = 0, part
        num_threads = min(block_size, dnum - block_size * m)
        if (lid < num_threads) then
          sdx(lid + 1) = dx(lid + block_size * m + 1)
          sdy(lid + 1) = dy(lid + block_size * m + 1)
          sdz(lid + 1) = dz(lid + block_size * m + 1)
        end if
        !$omp barrier
        do e = 1, block_size
          six_s = six_t - sdx(e)
          siy_s = siy_t - sdy(e)
          dist = six_s * six_s + siy_s * siy_s
          t = 1.0_real32 / (dist ** alpha)
          sum_dn = sum_dn + t
          sum_up = sum_up + t * sdz(e)
        end do
        !$omp barrier
      end do
      iz(tid) = sum_up / sum_dn
    end if
    !$omp end parallel
    !$omp end target teams
  end subroutine aidw_kernel_tiled

  subroutine reference_kernel(dx, dy, dz, dnum, ix, iy, iz, inum, area, avg_dist)
    real(real32), intent(in) :: dx(:), dy(:), dz(:), ix(:), iy(:), avg_dist(:), area
    real(real32), intent(out) :: iz(:)
    integer, intent(in) :: dnum, inum
    integer :: tid

    !$omp parallel do
    do tid = 1, inum
      iz(tid) = interpolate_point(tid, dx, dy, dz, dnum, ix, iy, area, avg_dist)
    end do
    !$omp end parallel do
  end subroutine reference_kernel

  real(real32) function interpolate_point(tid, dx, dy, dz, dnum, ix, iy, area, avg_dist) result(value)
    integer, intent(in) :: tid, dnum
    real(real32), intent(in) :: dx(:), dy(:), dz(:), ix(:), iy(:), avg_dist(:), area
    integer :: j
    real(real32) :: sum, dist, t, z, alpha, r_obs, r_exp, r_s0, u_r

    sum = 0.0_real32
    z = 0.0_real32
    alpha = 0.0_real32
    r_obs = avg_dist(tid)
    r_exp = 1.0_real32 / (2.0_real32 * sqrt(real(dnum, real32) / area))
    r_s0 = r_obs / r_exp
    u_r = 0.0_real32
    if (r_s0 >= r_min) u_r = 0.5_real32 - 0.5_real32 * cos(pi_value / r_max * (r_s0 - r_min))
    if (r_s0 >= r_max) u_r = 1.0_real32

    if (u_r >= 0.0_real32 .and. u_r <= 0.1_real32) alpha = a1
    if (u_r > 0.1_real32 .and. u_r <= 0.3_real32) &
      alpha = a1 * (1.0_real32 - 5.0_real32 * (u_r - 0.1_real32)) + a2 * 5.0_real32 * (u_r - 0.1_real32)
    if (u_r > 0.3_real32 .and. u_r <= 0.5_real32) &
      alpha = a3 * 5.0_real32 * (u_r - 0.3_real32) + a1 * (1.0_real32 - 5.0_real32 * (u_r - 0.3_real32))
    if (u_r > 0.5_real32 .and. u_r <= 0.7_real32) &
      alpha = a3 * (1.0_real32 - 5.0_real32 * (u_r - 0.5_real32)) + a4 * 5.0_real32 * (u_r - 0.5_real32)
    if (u_r > 0.7_real32 .and. u_r <= 0.9_real32) &
      alpha = a5 * 5.0_real32 * (u_r - 0.7_real32) + a4 * (1.0_real32 - 5.0_real32 * (u_r - 0.7_real32))
    if (u_r > 0.9_real32 .and. u_r <= 1.0_real32) alpha = a5
    alpha = alpha * 0.5_real32

    do j = 1, dnum
      dist = (ix(tid) - dx(j)) * (ix(tid) - dx(j)) + (iy(tid) - dy(j)) * (iy(tid) - dy(j))
      t = 1.0_real32 / (dist ** alpha)
      sum = sum + t
      z = z + dz(j) * t
    end do
    value = z / sum
  end function interpolate_point

  logical function verify(gold, test, len, tolerance) result(ok)
    real(real32), intent(in) :: gold(:), test(:), tolerance
    integer, intent(in) :: len
    integer :: i

    ok = .true.
    do i = 1, len
      if (abs(gold(i) - test(i)) > tolerance) then
        write(*,'(I0,1X,F0.6,1X,F0.6)') i - 1, gold(i), test(i)
        ok = .false.
        return
      end if
    end do
  end function verify

end program main
