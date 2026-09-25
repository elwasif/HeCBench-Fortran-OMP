! SPDX-License-Identifier: CC0-1.0
program gaussian_elim
  use iso_fortran_env, only: int64, real32, real64, output_unit
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: sp = real32
  integer, parameter :: block_size_0 = 256

  real(sp), allocatable :: a(:), b(:), m(:), final_vec(:)
  real(sp), allocatable :: a_host(:), b_host(:), m_host(:), final_vec_host(:)
  character(len=200) :: filename
  integer :: size, quiet, timing
  logical :: parse_error, ok
  real(real64) :: offload_start, offload_end

  filename = ""
  quiet = 0
  timing = 0
  size = -1

  parse_error = parse_commandline(filename, quiet, timing, size)
  if (parse_error) then
    call print_usage()
    stop 0
  end if

  if (size < 1) then
    call read_input_file(trim(filename), a, b, size)
  else
    allocate(a(size * size), b(size))
    call init_matrix(a, size)
    b = 1.0_sp
  end if

  if (quiet == 0) then
    write(*,'(A)') "The input matrix a is:"
    call print_mat(a, size, size, size)
    write(*,'(A)') "The input array b is:"
    call print_ary(b, size)
  end if

  allocate(m(size * size), final_vec(size))
  m = 0.0_sp

  allocate(a_host(size * size), b_host(size), m_host(size * size), final_vec_host(size))
  a_host = a
  b_host = b
  m_host = m
  call gaussian_reference(a_host, b_host, m_host, final_vec_host, size)

  offload_start = omp_get_wtime()
  call forward_sub(a, b, m, size, timing)
  offload_end = omp_get_wtime()

  if (timing /= 0) then
    write(*,'(A,I0,A)') "Device offloading time ", int((offload_end - offload_start) * 1000000.0_real64, int64), " (us)"
    write(*,*)
  end if

  call back_sub(a, b, final_vec, size)

  if (quiet == 0) then
    write(*,'(A)') "The result of array a is after forwardsub: "
    call print_mat(a, size, size, size)
    write(*,'(A)') "The result of array b is after forwardsub: "
    call print_ary(b, size)
    write(*,'(A)') "The result of matrix m is after forwardsub: "
    call print_mat(m, size, size, size)
    write(*,'(A)') "The solution is: "
    call print_ary(final_vec, size)
  end if

  write(*,'(A)') "Checking the results.."
  ok = check_results(final_vec, final_vec_host, size)
  if (ok) then
    write(*,'(A)') "PASS"
  else
    write(*,'(A)') "FAIL"
  end if

contains

  integer pure function idx(row, col, n) result(pos)
    integer, intent(in) :: row, col, n
    pos = row * n + col + 1
  end function idx

  logical function parse_commandline(filename, quiet, timing, size) result(has_error)
    character(len=*), intent(out) :: filename
    integer, intent(out) :: quiet, timing, size
    integer :: i, argc
    character(len=200) :: arg

    argc = command_argument_count()
    has_error = argc < 1
    if (has_error) return

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      if (len_trim(arg) > 0 .and. arg(1:1) == "-") then
        select case (arg(2:2))
        case ("s")
          i = i + 1
          call get_command_argument(i, arg)
          read(arg, *) size
          write(*,'(A,I0,A,I0,A)') "Create a square matrix (", size, " x ", size, ") internally"
        case ("f")
          i = i + 1
          call get_command_argument(i, filename)
          write(*,'(A,A,A)') "Read file from ", trim(filename), " "
        case ("h")
          has_error = .true.
          return
        case ("q")
          quiet = 1
        case ("t")
          timing = 1
        end select
      end if
      i = i + 1
    end do
  end function parse_commandline

  subroutine init_matrix(mat, n)
    real(sp), intent(out) :: mat(n * n)
    integer, intent(in) :: n
    real(sp), allocatable :: coe(:)
    real(sp) :: lamda, coe_i
    integer :: i, j

    allocate(coe(2 * n - 1))
    lamda = -0.01_sp
    do i = 0, n - 1
      coe_i = 10.0_sp * exp(lamda * real(i, sp))
      coe(n + i) = coe_i
      coe(n - i) = coe_i
    end do

    do i = 0, n - 1
      do j = 0, n - 1
        mat(idx(i,j,n)) = coe(n - i + j)
      end do
    end do
    deallocate(coe)
  end subroutine init_matrix

  subroutine read_input_file(path, mat, rhs, n)
    character(len=*), intent(in) :: path
    real(sp), allocatable, intent(out) :: mat(:), rhs(:)
    integer, intent(out) :: n
    integer :: unit, i

    open(newunit=unit, file=path, status="old", action="read")
    read(unit, *) n
    allocate(mat(n * n), rhs(n))
    do i = 1, n * n
      read(unit, *) mat(i)
    end do
    do i = 1, n
      read(unit, *) rhs(i)
    end do
    close(unit)
  end subroutine read_input_file

  subroutine gaussian_reference(a, b, m, final_vec, n)
    real(sp), intent(inout) :: a(n * n), b(n), m(n * n)
    real(sp), intent(out) :: final_vec(n)
    integer, intent(in) :: n
    integer :: t, i, x, y

    do t = 0, n - 2
      do i = 0, n - 2 - t
        m(idx(i + t + 1, t, n)) = a(idx(i + t + 1, t, n)) / a(idx(t, t, n))
      end do
      do x = 0, n - 2 - t
        do y = 0, n - 1 - t
          a(idx(x + t + 1, y + t, n)) = a(idx(x + t + 1, y + t, n)) - &
              m(idx(x + t + 1, t, n)) * a(idx(t, y + t, n))
          if (y == 0) then
            b(x + t + 2) = b(x + t + 2) - m(idx(x + t + 1, y + t, n)) * b(t + 1)
          end if
        end do
      end do
    end do

    call back_sub(a, b, final_vec, n)
  end subroutine gaussian_reference

  subroutine forward_sub(a, b, m, n, timing)
    real(sp), intent(inout) :: a(n * n), b(n), m(n * n)
    integer, intent(in) :: n, timing
    integer :: t, i, x, y
    real(real64) :: start_time, end_time

    !$omp target data map(tofrom: a(1:n*n), b(1:n), m(1:n*n))
    start_time = omp_get_wtime()

    do t = 0, n - 2
      !$omp target teams distribute parallel do thread_limit(block_size_0)
      do i = 0, n - 2 - t
        m(idx(i + t + 1, t, n)) = a(idx(i + t + 1, t, n)) / a(idx(t, t, n))
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(2) thread_limit(block_size_0)
      do x = 0, n - 2 - t
        do y = 0, n - 1 - t
          a(idx(x + t + 1, y + t, n)) = a(idx(x + t + 1, y + t, n)) - &
              m(idx(x + t + 1, t, n)) * a(idx(t, y + t, n))
          if (y == 0) then
            b(x + t + 2) = b(x + t + 2) - m(idx(x + t + 1, y + t, n)) * b(t + 1)
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end do

    end_time = omp_get_wtime()
    if (timing /= 0) then
      write(*,'(A,I0,A)') "Total kernel execution time ", int((end_time - start_time) * 1000000.0_real64, int64), " (us)"
    end if
    !$omp end target data
  end subroutine forward_sub

  subroutine back_sub(a, b, final_vec, n)
    real(sp), intent(in) :: a(n * n), b(n)
    real(sp), intent(out) :: final_vec(n)
    integer, intent(in) :: n
    integer :: i, j, row

    do i = 0, n - 1
      row = n - i - 1
      final_vec(row + 1) = b(row + 1)
      do j = 0, i - 1
        final_vec(row + 1) = final_vec(row + 1) - &
            a(idx(row, n - j - 1, n)) * final_vec(n - j)
      end do
      final_vec(row + 1) = final_vec(row + 1) / a(idx(row, row, n))
    end do
  end subroutine back_sub

  logical function check_results(final_vec, final_vec_host, n) result(ok)
    real(sp), intent(in) :: final_vec(n), final_vec_host(n)
    integer, intent(in) :: n
    integer :: i

    ok = .true.
    do i = 1, n
      if (abs(final_vec(i) - final_vec_host(i)) > 1.0e-3_sp) then
        ok = .false.
        write(*,'(A,I0,A,F10.6,A,F10.6,A)') "Result mismatch at index ", i - 1, ": ", &
            final_vec(i), "(device)  ", final_vec_host(i), "(host)"
      end if
    end do
  end function check_results

  subroutine print_usage()
    write(*,'(A)') "Gaussian Elimination Usage"
    write(*,*)
    write(*,'(A)') "gaussianElimination -f [filename] [-hqt]"
    write(*,*)
    write(*,'(A)') "example:"
    write(*,'(A)') "$ ./gaussianElimination matrix4.txt"
    write(*,*)
    write(*,'(A)') "filename     the filename that holds the matrix data"
    write(*,*)
    write(*,'(A)') "-h           Display the help file"
    write(*,'(A)') "-q           Quiet mode. Suppress all text output."
    write(*,'(A)') "-t           Print timing information."
    write(*,'(A)') "-s           Specifiy the matrix size when the path to a matrix data file is not set."
    write(*,*)
    write(*,*)
    write(*,'(A)') "Notes: 1. The filename is required as the first parameter."
    write(*,'(A)') "       2. If you declare either the device or the platform,"
    write(*,'(A)') "          you must declare both."
    write(*,*)
  end subroutine print_usage

  subroutine print_mat(mat, size, nrow, ncol)
    real(sp), intent(in) :: mat(size * size)
    integer, intent(in) :: size, nrow, ncol
    integer :: i, j

    do i = 0, nrow - 1
      do j = 0, ncol - 1
        write(*,'(ES8.2,1X)', advance="no") mat(idx(i,j,size))
      end do
      write(*,*)
    end do
    write(*,*)
  end subroutine print_mat

  subroutine print_ary(vec, n)
    real(sp), intent(in) :: vec(n)
    integer, intent(in) :: n
    integer :: i

    do i = 1, n
      write(*,'(ES8.2,1X)', advance="no") vec(i)
    end do
    write(*,*)
    write(*,*)
  end subroutine print_ary

end program gaussian_elim
