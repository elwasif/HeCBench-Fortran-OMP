! SPDX-License-Identifier: CC0-1.0
module perplexity_mod
  use iso_fortran_env, only: real32, real64
  use iso_c_binding, only: c_int
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name="rand")
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

contains

  subroutine fill_random(values)
    real(real32), intent(out) :: values(0:)
    integer :: i

    call c_srand(123_c_int)
    do i = 0, size(values) - 1
      values(i) = real(c_rand(), real32) / 2147483647.0_real32
    end do
  end subroutine fill_random

  subroutine perplexity_reference(distances, p, perplexity, epochs, tol, n, k)
    real(real32), intent(in) :: distances(0:)
    real(real32), intent(out) :: p(0:)
    real(real32), intent(in) :: perplexity, tol
    integer, intent(in) :: epochs, n, k
    integer :: i, j, step, ik
    real(real32) :: desired_entropy, beta_min, beta_max, beta, sum_pi
    real(real32) :: sum_disti_pi, div_value, entropy, entropy_diff

    desired_entropy = log(perplexity)
    do i = 0, n - 1
      beta_min = -huge(1.0_real32)
      beta_max = huge(1.0_real32)
      beta = 1.0_real32
      ik = i * k

      do step = 0, epochs - 1
        sum_pi = epsilon(1.0_real32)
        do j = 0, k - 1
          p(ik + j) = exp(-distances(ik + j) * beta)
          sum_pi = sum_pi + p(ik + j)
        end do

        sum_disti_pi = 0.0_real32
        div_value = 1.0_real32 / sum_pi
        do j = 0, k - 1
          p(ik + j) = p(ik + j) * div_value
          sum_disti_pi = sum_disti_pi + distances(ik + j) * p(ik + j)
        end do

        entropy = log(sum_pi) + beta * sum_disti_pi
        entropy_diff = entropy - desired_entropy
        if (abs(entropy_diff) <= tol) exit

        if (entropy_diff > 0.0_real32) then
          beta_min = beta
          if (beta_max == huge(1.0_real32)) then
            beta = beta * 2.0_real32
          else
            beta = (beta + beta_max) * 0.5_real32
          end if
        else
          beta_max = beta
          if (beta_min == -huge(1.0_real32)) then
            beta = beta * 0.5_real32
          else
            beta = (beta + beta_min) * 0.5_real32
          end if
        end if
      end do
    end do
  end subroutine perplexity_reference

  subroutine perplexity_search(distances, p, perplexity, epochs, tol, n, k)
    real(real32), intent(in) :: distances(0:)
    real(real32), intent(out) :: p(0:)
    real(real32), intent(in) :: perplexity, tol
    integer, intent(in) :: epochs, n, k
    integer :: i, j, step, ik
    real(real32) :: desired_entropy, beta_min, beta_max, beta, sum_pi
    real(real32) :: sum_disti_pi, div_value, entropy, entropy_diff

    desired_entropy = log(perplexity)

    !$omp target teams distribute parallel do thread_limit(256) private(j, step, ik, beta_min, beta_max, beta, sum_pi, sum_disti_pi, div_value, entropy, entropy_diff)
    do i = 0, n - 1
      beta_min = -huge(1.0_real32)
      beta_max = huge(1.0_real32)
      beta = 1.0_real32
      ik = i * k

      do step = 0, epochs - 1
        sum_pi = epsilon(1.0_real32)
        do j = 0, k - 1
          p(ik + j) = exp(-distances(ik + j) * beta)
          sum_pi = sum_pi + p(ik + j)
        end do

        sum_disti_pi = 0.0_real32
        div_value = 1.0_real32 / sum_pi
        do j = 0, k - 1
          p(ik + j) = p(ik + j) * div_value
          sum_disti_pi = sum_disti_pi + distances(ik + j) * p(ik + j)
        end do

        entropy = log(sum_pi) + beta * sum_disti_pi
        entropy_diff = entropy - desired_entropy
        if (abs(entropy_diff) <= tol) exit

        if (entropy_diff > 0.0_real32) then
          beta_min = beta
          if (beta_max == huge(1.0_real32)) then
            beta = beta * 2.0_real32
          else
            beta = (beta + beta_max) * 0.5_real32
          end if
        else
          beta_max = beta
          if (beta_min == -huge(1.0_real32)) then
            beta = beta * 0.5_real32
          else
            beta = (beta + beta_min) * 0.5_real32
          end if
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine perplexity_search

end module perplexity_mod

program main
  use iso_fortran_env, only: real32, real64
  use omp_lib
  use perplexity_mod
  implicit none

  integer :: argc, n, perplexity_i, repeat, n_nbrs, max_iter, status, total_size, iter, i
  real(real32), parameter :: tol = 1.0e-8_real32
  real(real32), allocatable :: data(:), h_data(:), distance(:)
  real(real64) :: time_total, start_time
  logical :: ok
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 3) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage: ', trim(prog), ' <number of points> <perplexity> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) n
  if (status /= 0) stop 1
  call get_command_argument(2, arg)
  read(arg, *, iostat=status) perplexity_i
  if (status /= 0) stop 1
  call get_command_argument(3, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  n_nbrs = 4 * perplexity_i
  max_iter = 100
  total_size = n * n_nbrs

  allocate(data(0:total_size - 1), h_data(0:total_size - 1), distance(0:total_size - 1))
  call fill_random(distance)
  data = 0.0_real32
  h_data = 0.0_real32
  time_total = 0.0_real64

  !$omp target data map(from: data(0:total_size - 1)) map(to: distance(0:total_size - 1))
  do iter = 1, repeat
    start_time = omp_get_wtime()
    call perplexity_search(distance, data, real(perplexity_i, real32), max_iter, tol, n, n_nbrs)
    time_total = time_total + (omp_get_wtime() - start_time)
  end do
  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', time_total / real(repeat, real64), ' (s)'
  !$omp end target data

  call perplexity_reference(distance, h_data, real(perplexity_i, real32), max_iter, tol, n, n_nbrs)

  ok = .true.
  do i = 0, total_size - 1
    if (abs(data(i) - h_data(i)) > 1.0e-3_real32) then
      write(*,'(I0,1X,F0.6,1X,F0.6)') i, data(i), h_data(i)
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(data, h_data, distance)
end program main
