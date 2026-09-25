! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
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

  integer, parameter :: grid_size = 21
  integer, parameter :: n_particles = 147456
  integer, parameter :: grid_cells = grid_size * grid_size
  real(real32), parameter :: radius = 0.5_real32
  integer :: n_iterations, n_repeat, random_count, map_size
  integer :: repeat, mismatches
  integer(int64), allocatable :: map(:), map_ref(:)
  real(real32), allocatable :: random_x(:), random_y(:), particle_x(:), particle_y(:)
  real(real64) :: sim_start, sim_end, kernel_total, kernel_start, kernel_end

  if (command_argument_count() /= 2) then
    print '(A)', ' Incorrect number of parameters '
    print '(A)', ' Usage: ./main <Number of iterations within the kernel> <Kernel execution count>'
    stop 1
  end if

  n_iterations = read_arg(1)
  n_repeat = read_arg(2)

  random_count = n_particles * n_iterations
  map_size = n_particles * grid_cells
  allocate(random_x(random_count), random_y(random_count), particle_x(n_particles), particle_y(n_particles), &
           map(map_size), map_ref(map_size))

  call initialize_random(random_x, random_y)
  particle_x = 10.0_real32
  particle_y = 10.0_real32
  map = 0_int64
  map_ref = 0_int64

  sim_start = omp_get_wtime()
  call motion_device(particle_x, particle_y, random_x, random_y, map, n_iterations, n_repeat, kernel_total)
  sim_end = omp_get_wtime()

  call motion_host(random_x, random_y, map_ref, n_iterations)
  mismatches = count(map /= map_ref)

  print '(A,I0)', ' The number of kernel execution is ', n_repeat
  print '(A,I0)', ' The number of particles is ', n_particles
  print '(A)'
  print '(A,F0.6,A)', 'Average kernel execution time: ', kernel_total / real(n_repeat, real64), ' (s)'
  print '(A)'
  print '(A,F0.6,A)', 'Simulation time: ', sim_end - sim_start, ' (s) '
  if (mismatches <= 2) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(random_x, random_y, particle_x, particle_y, map, map_ref)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_random(random_x, random_y)
    real(real32), intent(out) :: random_x(:), random_y(:)
    integer :: i
    call c_srand(17_c_int)
    do i = 1, size(random_x)
      random_x(i) = real(mod(c_rand(), 100_c_int), real32)
      random_y(i) = real(mod(c_rand(), 100_c_int), real32)
    end do
  end subroutine initialize_random

  subroutine motion_device(particle_x, particle_y, random_x, random_y, map, n_iterations, n_repeat, kernel_total)
    real(real32), intent(inout) :: particle_x(:), particle_y(:)
    real(real32), intent(in) :: random_x(:), random_y(:)
    integer(int64), intent(inout) :: map(:)
    integer, intent(in) :: n_iterations, n_repeat
    real(real64), intent(out) :: kernel_total
    integer :: repeat, idx, particle, iter
    integer :: ix, iy, map_base, map_index
    real(real32) :: px, py, rand_x, rand_y, displacement_x, displacement_y, dx, dy

    kernel_total = 0.0_real64
    !$omp target data map(to: random_x(1:size(random_x)), random_y(1:size(random_y))) &
    !$omp& map(alloc: particle_x(1:size(particle_x)), particle_y(1:size(particle_y))) map(from: map(1:size(map)))
    do repeat = 1, n_repeat
      !$omp target update to(particle_x(1:size(particle_x)))
      !$omp target update to(particle_y(1:size(particle_y)))
      !$omp target update to(map(1:size(map)))

      kernel_start = omp_get_wtime()
      !$omp target teams distribute parallel do thread_limit(256) &
      !$omp& private(iter, px, py, rand_x, rand_y, displacement_x, displacement_y, dx, dy, ix, iy, map_base, map_index)
      do particle = 1, n_particles
        px = particle_x(particle)
        py = particle_y(particle)
        map_base = (particle - 1) * grid_cells
        do iter = 1, n_iterations
          rand_x = random_x((iter - 1) * n_particles + particle)
          rand_y = random_y((iter - 1) * n_particles + particle)
          displacement_x = rand_x / 1000.0_real32 - 0.0495_real32
          displacement_y = rand_y / 1000.0_real32 - 0.0495_real32
          px = px + displacement_x
          py = py + displacement_y
          dx = px - real(int(px), real32)
          dy = py - real(int(py), real32)
          ix = int(floor(px))
          iy = int(floor(py))
          if (px < real(grid_size, real32) .and. py < real(grid_size, real32) .and. px >= 0.0_real32 .and. py >= 0.0_real32) then
            if (dx * dx + dy * dy <= radius * radius) then
              map_index = map_base + iy * grid_size + ix + 1
              map(map_index) = map(map_index) + 1_int64
            end if
          end if
        end do
        particle_x(particle) = px
        particle_y(particle) = py
      end do
      !$omp end target teams distribute parallel do
      kernel_end = omp_get_wtime()
      kernel_total = kernel_total + kernel_end - kernel_start
    end do
    !$omp end target data
  end subroutine motion_device

  subroutine motion_host(random_x, random_y, map, n_iterations)
    real(real32), intent(in) :: random_x(:), random_y(:)
    integer(int64), intent(inout) :: map(:)
    integer, intent(in) :: n_iterations
    integer :: particle, iter
    integer :: ix, iy, map_base, map_index
    real(real32) :: px, py, rand_x, rand_y, displacement_x, displacement_y, dx, dy

    map = 0_int64
    do particle = 1, n_particles
      px = 10.0_real32
      py = 10.0_real32
      map_base = (particle - 1) * grid_cells
      do iter = 1, n_iterations
        rand_x = random_x((iter - 1) * n_particles + particle)
        rand_y = random_y((iter - 1) * n_particles + particle)
        displacement_x = rand_x / 1000.0_real32 - 0.0495_real32
        displacement_y = rand_y / 1000.0_real32 - 0.0495_real32
        px = px + displacement_x
        py = py + displacement_y
        dx = px - real(int(px), real32)
        dy = py - real(int(py), real32)
        ix = int(floor(px))
        iy = int(floor(py))
        if (px < real(grid_size, real32) .and. py < real(grid_size, real32) .and. px >= 0.0_real32 .and. py >= 0.0_real32) then
          if (dx * dx + dy * dy <= radius * radius) then
            map_index = map_base + iy * grid_size + ix + 1
            map(map_index) = map(map_index) + 1_int64
          end if
        end if
      end do
    end do
  end subroutine motion_host

end program main
