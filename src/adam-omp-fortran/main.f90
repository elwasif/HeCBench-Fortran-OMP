! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: vector_size, time_step, repeat
  real(real32), allocatable :: m(:), v(:), g(:), p(:), r(:), m_ref(:), v_ref(:)
  real(real32), parameter :: step_size = 1.0e-3_real32
  real(real32), parameter :: decay = 0.5_real32
  real(real32), parameter :: beta1 = 0.9_real32
  real(real32), parameter :: beta2 = 0.999_real32
  real(real32), parameter :: eps = 1.0e-10_real32
  real(real32), parameter :: grad_scale = 256.0_real32
  real(real32), parameter :: rand_max = 2147483647.0_real32
  real(real64) :: start_time, end_time, cr, cp
  integer :: i
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <vector size> <number of time steps> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) vector_size
  read(arg2, *) time_step
  read(arg3, *) repeat

  allocate(m(vector_size), v(vector_size), g(vector_size), p(vector_size), r(vector_size))
  allocate(m_ref(vector_size), v_ref(vector_size))
  call c_srand(19937_c_int)
  call initialize_vectors(vector_size, m, v, g, p)
  r = p
  m_ref = m
  v_ref = v

  !$omp target data map(to: m(1:vector_size), v(1:vector_size), g(1:vector_size)) map(tofrom: p(1:vector_size))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call adam_kernel(p, m, v, g, beta1, beta2, eps, grad_scale, step_size, time_step, vector_size, decay)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time ', &
    (end_time - start_time) * 1.0e3_real64 / real(repeat, real64), ' (ms)'
  !$omp end target data

  call reference(repeat, r, m_ref, v_ref, g, beta1, beta2, eps, grad_scale, step_size, &
    time_step, vector_size, decay)

  ok = .true.
  cr = 0.0_real64
  cp = 0.0_real64
  do i = 1, vector_size
    if (abs(r(i) - p(i)) > 1.0e-3_real32) ok = .false.
    cr = cr + real(r(i), real64)
    cp = cp + real(p(i), real64)
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if
  write(*,'(A,F0.6,1X,F0.6)') 'Checksum: ', cr / real(vector_size, real64), cp / real(vector_size, real64)

  deallocate(m, v, g, p, r, m_ref, v_ref)

contains

  subroutine initialize_vectors(vector_size, m, v, g, p)
    integer, intent(in) :: vector_size
    real(real32), intent(out) :: m(:), v(:), g(:), p(:)
    integer :: i

    do i = 1, vector_size
      m(i) = random_uniform_unit()
      v(i) = random_uniform_unit()
      g(i) = random_uniform_unit()
      p(i) = random_uniform_unit()
    end do
  end subroutine initialize_vectors

  function random_uniform_unit() result(value)
    real(real32) :: value

    value = real(c_rand(), real32) / rand_max
  end function random_uniform_unit

  subroutine adam_kernel(p, m, v, g, b1, b2, eps, grad_scale, step_size, time_step, vector_size, decay)
    real(real32), intent(inout) :: p(:), m(:), v(:)
    real(real32), intent(in) :: g(:), b1, b2, eps, grad_scale, step_size, decay
    integer, intent(in) :: time_step, vector_size
    integer :: j, t
    real(real32) :: scaled_grad, m_corrected, v_corrected, denom, update

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(j, t, scaled_grad, m_corrected, v_corrected, denom, update)
    do j = 1, vector_size
      do t = 1, time_step
        scaled_grad = g(j) / grad_scale
        m(j) = b1 * m(j) + (1.0_real32 - b1) * scaled_grad
        v(j) = b2 * v(j) + (1.0_real32 - b2) * scaled_grad * scaled_grad
        m_corrected = m(j) / (1.0_real32 - b1 ** t)
        v_corrected = v(j) / (1.0_real32 - b2 ** t)
        denom = sqrt(v_corrected + eps)
        update = (m_corrected / denom) + (decay * p(j))
        p(j) = p(j) - (step_size * update)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine adam_kernel

  subroutine reference(repeat, p, m, v, g, b1, b2, eps, grad_scale, step_size, time_step, vector_size, decay)
    integer, intent(in) :: repeat, time_step, vector_size
    real(real32), intent(inout) :: p(:), m(:), v(:)
    real(real32), intent(in) :: g(:), b1, b2, eps, grad_scale, step_size, decay
    integer :: i, j, t
    real(real32) :: scaled_grad, m_corrected, v_corrected, denom, update

    do i = 1, repeat
      do j = 1, vector_size
        do t = 1, time_step
          scaled_grad = g(j) / grad_scale
          m(j) = b1 * m(j) + (1.0_real32 - b1) * scaled_grad
          v(j) = b2 * v(j) + (1.0_real32 - b2) * scaled_grad * scaled_grad
          m_corrected = m(j) / (1.0_real32 - b1 ** t)
          v_corrected = v(j) / (1.0_real32 - b2 ** t)
          denom = sqrt(v_corrected + eps)
          update = (m_corrected / denom) + (decay * p(j))
          p(j) = p(j) - (step_size * update)
        end do
      end do
    end do
  end subroutine reference

end program main
