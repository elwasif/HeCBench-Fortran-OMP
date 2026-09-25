! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: dim = 30_int32
  real(real32), parameter :: start_range_min = -5.12_real32
  real(real32), parameter :: start_range_max = 5.12_real32
  real(real32), parameter :: omega = 0.5_real32
  real(real32), parameter :: c1 = 1.5_real32
  real(real32), parameter :: c2 = 1.5_real32
  real(real32), parameter :: phi = 3.1415_real32
  real(real32), parameter :: rand_max = 2147483647.0_real32
  real(real32), parameter :: tolerance = 1.0e-3_real32

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

  character(len=256) :: arg0, arg1, arg2
  integer(int32) :: p, r, i, size_total
  real(real32), allocatable :: positions(:), velocities(:), p_bests(:), p_bests_ref(:), g_best(:)
  real(real64) :: avg_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of particles> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) p
  read(arg2, *) r

  if (p <= 0_int32 .or. r <= 0_int32) then
    print '(A)', 'FAIL'
    stop 1
  end if

  print '(A,I0)', 'Number of particles is ', p
  print '(A,I0)', 'Number of dimensions is ', dim

  size_total = p * dim
  allocate(positions(1:size_total), velocities(1:size_total))
  allocate(p_bests(1:size_total), p_bests_ref(1:size_total), g_best(1:dim))

  call c_srand(123_c_int)
  do i = 1, size_total
    positions(i) = get_random(start_range_min, start_range_max)
    p_bests(i) = positions(i)
    p_bests_ref(i) = positions(i)
    velocities(i) = 0.0_real32
  end do
  g_best = p_bests(1:dim)

  print '(A)'
  print '(A)', 'Execute PSO on a device'
  call gpu_pso(p, r, positions, velocities, p_bests, g_best, avg_us)
  print '(A,F0.6)', 'Result=', host_fitness_function(g_best)

  print '(A)'
  print '(A)', 'Execute PSO on a host. This may take a while for large problem size..'
  g_best = p_bests_ref(1:dim)
  call pso_host(p, r, positions, velocities, p_bests_ref, g_best)
  print '(A,F0.6)', 'Result=', host_fitness_function(g_best)

  ok = compare_first_particle(p_bests_ref, p_bests)
  print '(A)', merge('PASS', 'FAIL', ok)

contains

  pure real(real32) function f_map(x) result(value)
    real(real32), intent(in) :: x
    value = 1.0_real32 + (x - 1.0_real32) / 4.0_real32
  end function f_map

  real(real32) function host_fitness_function(x) result(res)
    real(real32), intent(in) :: x(1:dim)
    integer(int32) :: i
    real(real32) :: y1, yn, y, yp

    y1 = f_map(x(1))
    yn = f_map(x(dim))
    res = sin(phi * y1) ** 2 + (yn - 1.0_real32) ** 2

    do i = 1, dim - 1
      y = f_map(x(i))
      yp = f_map(x(i + 1))
      res = res + (y - 1.0_real32) ** 2 * (1.0_real32 + 10.0_real32 * sin(phi * yp) ** 2)
    end do
  end function host_fitness_function

  real(real32) function get_random(low, high) result(value)
    real(real32), intent(in) :: low, high
    value = low + (((high - low) + 1.0_real32) * real(c_rand(), real32) / (rand_max + 1.0_real32))
  end function get_random

  real(real32) function get_random_clamped(seed) result(value)
    integer(int32), intent(in) :: seed
    call c_srand(int(seed, c_int))
    value = real(c_rand(), real32) / rand_max
  end function get_random_clamped

  subroutine pso_host(p, r, positions, velocities, p_bests, g_best)
    integer(int32), intent(in) :: p, r
    real(real32), intent(inout) :: positions(1:), velocities(1:), p_bests(1:), g_best(1:dim)
    real(real32) :: temp_particle1(1:dim), temp_particle2(1:dim)
    real(real32) :: rp, rg
    integer(int32) :: iter, idx, j

    do iter = 0, r - 1
      rp = get_random_clamped(iter)
      rg = get_random_clamped(r - iter)

      do idx = 1, p * dim
        velocities(idx) = omega * velocities(idx) + &
                          c1 * rp * (p_bests(idx) - positions(idx)) + &
                          c2 * rg * (g_best(mod(idx - 1, dim) + 1) - positions(idx))
        positions(idx) = positions(idx) + velocities(idx)
      end do

      do idx = 1, p * dim, dim
        do j = 1, dim
          temp_particle1(j) = positions(idx + j - 1)
          temp_particle2(j) = p_bests(idx + j - 1)
        end do

        if (host_fitness_function(temp_particle1) < host_fitness_function(temp_particle2)) then
          do j = 1, dim
            p_bests(idx + j - 1) = temp_particle1(j)
          end do

          if (host_fitness_function(temp_particle1) < 130.0_real32) then
            do j = 1, dim
              g_best(j) = g_best(j) + temp_particle1(j)
            end do
          end if
        end if
      end do
    end do
  end subroutine pso_host

  subroutine gpu_pso(p, r, positions, velocities, p_bests, g_best, avg_us)
    integer(int32), intent(in) :: p, r
    real(real32), intent(inout) :: positions(1:), velocities(1:), p_bests(1:), g_best(1:dim)
    real(real64), intent(out) :: avg_us
    integer(int32) :: iter, idx, particle, j, base, size_total
    real(real32) :: rp, rg
    real(real32) :: temp_particle1(1:dim), temp_particle2(1:dim)
    real(real32) :: fit1, fit2, y1, yn, y, yp
    real(real64) :: start_time, end_time

    size_total = p * dim

    !$omp target data map(to: positions(1:size_total), velocities(1:size_total)) &
    !$omp& map(tofrom: g_best(1:dim), p_bests(1:size_total))
    start_time = omp_get_wtime()
    do iter = 0, r - 1
      rp = get_random_clamped(iter)
      rg = get_random_clamped(r - iter)

      !$omp target teams distribute parallel do thread_limit(256) firstprivate(p, rp, rg)
      do idx = 1, size_total
        velocities(idx) = omega * velocities(idx) + &
                          c1 * rp * (p_bests(idx) - positions(idx)) + &
                          c2 * rg * (g_best(mod(idx - 1, dim) + 1) - positions(idx))
        positions(idx) = positions(idx) + velocities(idx)
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do thread_limit(256) &
      !$omp& private(base, j, temp_particle1, temp_particle2, fit1, fit2, y1, yn, y, yp)
      do particle = 0, p - 1
        base = particle * dim + 1

        do j = 1, dim
          temp_particle1(j) = positions(base + j - 1)
          temp_particle2(j) = p_bests(base + j - 1)
        end do

        y1 = 1.0_real32 + (temp_particle1(1) - 1.0_real32) / 4.0_real32
        yn = 1.0_real32 + (temp_particle1(dim) - 1.0_real32) / 4.0_real32
        fit1 = sin(phi * y1) ** 2 + (yn - 1.0_real32) ** 2
        do j = 1, dim - 1
          y = 1.0_real32 + (temp_particle1(j) - 1.0_real32) / 4.0_real32
          yp = 1.0_real32 + (temp_particle1(j + 1) - 1.0_real32) / 4.0_real32
          fit1 = fit1 + (y - 1.0_real32) ** 2 * (1.0_real32 + 10.0_real32 * sin(phi * yp) ** 2)
        end do

        y1 = 1.0_real32 + (temp_particle2(1) - 1.0_real32) / 4.0_real32
        yn = 1.0_real32 + (temp_particle2(dim) - 1.0_real32) / 4.0_real32
        fit2 = sin(phi * y1) ** 2 + (yn - 1.0_real32) ** 2
        do j = 1, dim - 1
          y = 1.0_real32 + (temp_particle2(j) - 1.0_real32) / 4.0_real32
          yp = 1.0_real32 + (temp_particle2(j + 1) - 1.0_real32) / 4.0_real32
          fit2 = fit2 + (y - 1.0_real32) ** 2 * (1.0_real32 + 10.0_real32 * sin(phi * yp) ** 2)
        end do

        if (fit1 < fit2) then
          do j = 1, dim
            p_bests(base + j - 1) = temp_particle1(j)
          end do

          if (fit1 < 130.0_real32) then
            do j = 1, dim
              !$omp atomic update
              g_best(j) = g_best(j) + temp_particle1(j)
            end do
          end if
        end if
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data

    avg_us = ((end_time - start_time) * 1.0d6) / real(r, real64)
    print '(A,F0.6,A)', 'Average kernel execution time ', avg_us, ' (us)'
  end subroutine gpu_pso

  logical function compare_first_particle(p_bests_ref, p_bests) result(ok)
    real(real32), intent(in) :: p_bests_ref(1:), p_bests(1:)
    integer(int32) :: i

    ok = .true.
    do i = 1, dim
      if (abs(p_bests_ref(i) - p_bests(i)) > tolerance) then
        print '(A,I0,1X,F0.6,1X,F0.6)', '@', i - 1, p_bests_ref(i), p_bests(i)
        ok = .false.
        exit
      end if
    end do
  end function compare_first_particle

end program main
