! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  interface
    subroutine libc_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine libc_srand

    function libc_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function libc_rand
  end interface

  integer, parameter :: num_elements = 10000000
  integer, parameter :: projectile_fields = 5
  integer, parameter :: angle_field = 1
  integer, parameter :: velocity_field = 2
  integer, parameter :: range_field = 3
  integer, parameter :: total_time_field = 4
  integer, parameter :: max_height_field = 5
  integer, parameter :: block_size = 256
  real(real32), parameter :: k_pi_value = 3.1415_real32
  real(real32), parameter :: k_g_value = 9.81_real32

  character(len=256) :: arg0, arg
  integer :: repeat, i, iter, errors, obj
  real(real32), allocatable :: input_vect(:), out_parallel_vect(:), out_scalar_vect(:)
  real(real32) :: proj_angle, proj_vel, sin_value, cos_value
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) stop 1

  allocate(input_vect(projectile_fields * num_elements))
  allocate(out_parallel_vect(projectile_fields * num_elements))
  allocate(out_scalar_vect(projectile_fields * num_elements))

  input_vect = 0.0_real32
  call libc_srand(2_c_int)
  do i = 1, num_elements
    obj = projectile_offset(i)
    input_vect(obj + angle_field) = real(mod(libc_rand(), 90_c_int) + 10_c_int, real32)
    input_vect(obj + velocity_field) = real(mod(libc_rand(), 400_c_int) + 10_c_int, real32)
    input_vect(obj + range_field) = 1.0_real32
    input_vect(obj + total_time_field) = 1.0_real32
    input_vect(obj + max_height_field) = 1.0_real32
  end do

  out_parallel_vect = 0.0_real32
  out_scalar_vect = 0.0_real32

  !$omp target data map(to: input_vect(1:projectile_fields*num_elements)) &
  !$omp& map(from: out_parallel_vect(1:projectile_fields*num_elements))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 1, num_elements
      obj = (i - 1) * projectile_fields
      proj_angle = input_vect(obj + angle_field)
      proj_vel = input_vect(obj + velocity_field)
      sin_value = sin(proj_angle * k_pi_value / 180.0_real32)
      cos_value = cos(proj_angle * k_pi_value / 180.0_real32)
      out_parallel_vect(obj + total_time_field) = abs(2.0_real32 * proj_vel * sin_value) / k_g_value
      out_parallel_vect(obj + range_field) = abs(proj_vel * out_parallel_vect(obj + total_time_field) * cos_value)
      out_parallel_vect(obj + angle_field) = proj_angle
      out_parallel_vect(obj + velocity_field) = proj_vel
      out_parallel_vect(obj + max_height_field) = (proj_vel * proj_vel * sin_value * sin_value) / 2.0_real32 * k_g_value
    end do
    !$omp end target teams distribute parallel do
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F8.6,A)') 'Average kernel execution time: ', elapsed / real(repeat, real64), ' (s)'

  do i = 1, num_elements
    obj = projectile_offset(i)
    call compute_projectile(input_vect(obj + angle_field), input_vect(obj + velocity_field), &
                            out_scalar_vect(obj + 1:obj + projectile_fields))
  end do

  errors = 0
  do i = 1, num_elements
    obj = projectile_offset(i)
    if (projectile_differs(out_parallel_vect(obj + 1:obj + projectile_fields), &
                           out_scalar_vect(obj + 1:obj + projectile_fields))) then
      errors = errors + 1
      exit
    end if
  end do

  if (errors == 0) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(input_vect, out_parallel_vect, out_scalar_vect)

contains

  integer function projectile_offset(index)
    integer, intent(in) :: index

    projectile_offset = (index - 1) * projectile_fields
  end function projectile_offset

  subroutine compute_projectile(proj_angle, proj_vel, pObj)
    real(real32), intent(in) :: proj_angle, proj_vel
    real(real32), intent(out) :: pObj(projectile_fields)
    real(real32) :: ref_sin_value, ref_cos_value

    ref_sin_value = sin(proj_angle * k_pi_value / 180.0_real32)
    ref_cos_value = cos(proj_angle * k_pi_value / 180.0_real32)
    pObj(total_time_field) = abs(2.0_real32 * proj_vel * ref_sin_value) / k_g_value
    pObj(range_field) = abs(proj_vel * pObj(total_time_field) * ref_cos_value)
    pObj(angle_field) = proj_angle
    pObj(velocity_field) = proj_vel
    pObj(max_height_field) = (proj_vel * proj_vel * ref_sin_value * ref_sin_value) / 2.0_real32 * k_g_value
  end subroutine compute_projectile

  logical function projectile_differs(a, b)
    real(real32), intent(in) :: a(projectile_fields), b(projectile_fields)

    projectile_differs = abs(a(angle_field) - b(angle_field)) > 1.0_real32 .or. &
                         abs(a(velocity_field) - b(velocity_field)) > 1.0_real32 .or. &
                         abs(a(range_field) - b(range_field)) > 1.0_real32 .or. &
                         abs(a(total_time_field) - b(total_time_field)) > 1.0_real32 .or. &
                         abs(a(max_height_field) - b(max_height_field)) > 1.0_real32
  end function projectile_differs

end program main
