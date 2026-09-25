! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: m_size = 512_int32 * 8_int32
  integer(int32), parameter :: m_dim = m_size / 8_int32
  integer(int32), parameter :: n_dim = m_size / 4_int32
  integer(int32), parameter :: k_dim = m_size / 2_int32
  real(real32), parameter :: value_tolerance = 1.0e-5_real32

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

  character(len=256) :: arg0, arg1
  integer(int32) :: repeat, i, j, k, m, iter, c_idx
  real(real32), allocatable :: a_host(:), b_host(:), c_host(:), c_back(:)
  real(real32) :: p, one_over_p, sum_value
  real(real64) :: start_time, end_time, avg_seconds

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0_int32) then
    stop 1
  end if

  allocate(a_host(1:m_dim * n_dim))
  allocate(b_host(1:n_dim * k_dim))
  allocate(c_host(1:m_dim * k_dim))
  allocate(c_back(1:m_dim * k_dim))

  a_host = 1.0_real32 / real(n_dim, real32)

  call c_srand(123_c_int)
  do i = 1, n_dim
    do j = 1, k_dim
      b_host((i - 1_int32) * k_dim + j) = real(modulo(c_rand(), 256_c_int), real32)
    end do
  end do

  call normalize_columns(b_host, n_dim, k_dim)

  write(*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'Problem size: c(', m_dim, ',', k_dim, ') = a(', m_dim, ',', &
                                                 n_dim, ') * b(', n_dim, ',', k_dim, ')'

  !$omp target data map(to: a_host(1:m_dim * n_dim), b_host(1:n_dim * k_dim)) &
  !$omp& map(alloc: c_back(1:m_dim * k_dim))
  do m = 1, 4
    print '(A,I0)', 'Minkowski distance with p = ', m
    p = real(m, real32)
    one_over_p = 1.0_real32 / p

    start_time = omp_get_wtime()
    do iter = 1, repeat
      !$omp target teams distribute parallel do collapse(2) thread_limit(256) private(j, k, c_idx, sum_value)
      do i = 1, m_dim
        do j = 1, k_dim
          sum_value = 0.0_real32
          do k = 1, n_dim
            sum_value = sum_value + abs(a_host((i - 1_int32) * n_dim + k) - &
                                          b_host((k - 1_int32) * k_dim + j)) ** p
          end do
          c_idx = (i - 1_int32) * k_dim + j
          c_back(c_idx) = sum_value ** one_over_p
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()

    avg_seconds = (end_time - start_time) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average kernel execution time: ', avg_seconds, ' (s)'

    !$omp target update from(c_back(1:m_dim * k_dim))
    call verify_result(a_host, b_host, c_host, c_back, p, one_over_p)
  end do
  !$omp end target data

contains

  subroutine normalize_columns(matrix, row_count, col_count)
    real(real32), intent(inout) :: matrix(1:)
    integer(int32), intent(in) :: row_count, col_count
    integer(int32) :: row, col, idx
    real(real32) :: sum_value

    do col = 1, col_count
      sum_value = 0.0_real32
      do row = 1, row_count
        idx = (row - 1_int32) * col_count + col
        sum_value = sum_value + matrix(idx)
      end do
      do row = 1, row_count
        idx = (row - 1_int32) * col_count + col
        matrix(idx) = matrix(idx) / sum_value
      end do
    end do
  end subroutine normalize_columns

  logical function value_same(lhs, rhs) result(ok)
    real(real32), intent(in) :: lhs, rhs
    ok = abs(lhs - rhs) <= value_tolerance
  end function value_same

  subroutine verify_result(a_matrix, b_matrix, c_matrix, c_device, p, one_over_p)
    real(real32), intent(in) :: a_matrix(1:), b_matrix(1:)
    real(real32), intent(inout) :: c_matrix(1:)
    real(real32), intent(in) :: c_device(1:)
    real(real32), intent(in) :: p, one_over_p
    integer(int32) :: i, j, k, c_idx, print_count
    logical :: mismatch_found

    c_matrix = 0.0_real32
    do i = 1, m_dim
      do k = 1, n_dim
        do j = 1, k_dim
          c_idx = (i - 1_int32) * k_dim + j
          c_matrix(c_idx) = c_matrix(c_idx) + abs(a_matrix((i - 1_int32) * n_dim + k) - &
                                                 b_matrix((k - 1_int32) * k_dim + j)) ** p
        end do
      end do
    end do

    do i = 1, m_dim
      do j = 1, k_dim
        c_idx = (i - 1_int32) * k_dim + j
        c_matrix(c_idx) = c_matrix(c_idx) ** one_over_p
      end do
    end do

    mismatch_found = .false.
    print_count = 0

    do i = 1, m_dim
      do j = 1, k_dim
        c_idx = (i - 1_int32) * k_dim + j
        if (.not. value_same(c_device(c_idx), c_matrix(c_idx))) then
          write(*, '(A,I0,A,I0,A,F0.6,A,F0.6)') 'Fail - The result is incorrect for element: [', i - 1, ', ', &
                                                 j - 1, '], expected: ', c_matrix(c_idx), ', but found: ', c_device(c_idx)
          mismatch_found = .true.
          print_count = print_count + 1
          if (print_count == 5) exit
        end if
      end do
      if (print_count == 5) exit
    end do

    if (.not. mismatch_found) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine verify_result

end program main
