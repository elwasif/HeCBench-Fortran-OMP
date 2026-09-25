! SPDX-License-Identifier: CC0-1.0
module simple_spmv_mod
  use iso_c_binding, only: c_double, c_int, c_long
  use iso_fortran_env, only: int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name="rand")
      import :: c_int
    end function c_rand

    subroutine c_srand48(seed) bind(C, name="srand48")
      import :: c_long
      integer(c_long), value :: seed
    end subroutine c_srand48

    real(c_double) function c_drand48() bind(C, name="drand48")
      import :: c_double
    end function c_drand48
  end interface

contains

  subroutine init_vector(vector, m)
    real(real32), intent(out) :: vector(0:)
    integer, intent(in) :: m
    integer :: i

    do i = 0, m - 1
      vector(i) = real(c_drand48(), real32)
    end do
  end subroutine init_vector

  subroutine init_matrix(matrix, num_rows, nnz)
    real(real32), intent(out) :: matrix(0:)
    integer, intent(in) :: num_rows, nnz
    real(real32), allocatable :: d(:)
    integer :: n, i, j, a, b
    real(real32) :: tmp

    n = num_rows * num_rows
    allocate(d(0:n - 1))

    call c_srand(123_c_int)
    do i = 0, n - 1
      d(i) = real(i, real32)
    end do

    do i = n, 1, -1
      a = i - 1
      b = mod(c_rand(), i)
      if (a /= b) then
        tmp = d(a)
        d(a) = d(b)
        d(b) = tmp
      end if
    end do

    do i = 0, num_rows - 1
      do j = 0, num_rows - 1
        if (d(i * num_rows + j) >= real(nnz, real32)) then
          matrix(i * num_rows + j) = 0.0_real32
        else
          matrix(i * num_rows + j) = real(c_drand48(), real32) + 1.0_real32
        end if
      end do
    end do

    deallocate(d)
  end subroutine init_matrix

  subroutine init_csr(row_indices, values, col_indices, matrix, num_rows, nnz)
    integer, intent(out) :: row_indices(0:), col_indices(0:)
    real(real32), intent(out) :: values(0:)
    real(real32), intent(in) :: matrix(0:)
    integer, intent(in) :: num_rows, nnz
    integer, allocatable :: non_zero_elements(:)
    integer :: i, j, tmp, nnz_per_row

    allocate(non_zero_elements(0:num_rows - 1))
    row_indices(num_rows) = nnz
    row_indices(0) = 0
    tmp = 0

    do i = 0, num_rows - 1
      nnz_per_row = 0
      do j = 0, num_rows - 1
        if (matrix(i * num_rows + j) /= 0.0_real32) then
          values(tmp) = matrix(i * num_rows + j)
          col_indices(tmp) = j
          tmp = tmp + 1
          nnz_per_row = nnz_per_row + 1
        end if
      end do
      non_zero_elements(i) = nnz_per_row
    end do

    if (tmp /= nnz) error stop "CSR nonzero count mismatch"

    do i = 1, num_rows - 1
      row_indices(i) = row_indices(i - 1) + non_zero_elements(i - 1)
    end do

    deallocate(non_zero_elements)
  end subroutine init_csr

  subroutine mv_csr_serial(num_rows, row_indices, col_indices, values, x, y)
    integer, intent(in) :: num_rows
    integer, intent(in) :: row_indices(0:), col_indices(0:)
    real(real32), intent(in) :: values(0:), x(0:)
    real(real32), intent(out) :: y(0:)
    integer :: row, row_start, row_end, i
    real(real32) :: dot

    do row = 0, num_rows - 1
      row_start = row_indices(row)
      row_end = row_indices(row + 1)
      dot = 0.0_real32
      do i = row_start, row_end - 1
        dot = dot + values(i) * x(col_indices(i))
      end do
      y(row) = dot
    end do
  end subroutine mv_csr_serial

  real(real32) function check_rate(a, b, n) result(rate)
    real(real32), intent(in) :: a(0:), b(0:)
    integer, intent(in) :: n
    integer :: i
    real(real64) :: diff_sum, sum_b

    diff_sum = 0.0_real64
    sum_b = 0.0_real64
    do i = 0, n - 1
      diff_sum = diff_sum + abs(real(a(i), real64) - real(b(i), real64))
      sum_b = sum_b + abs(real(b(i), real64))
    end do
    rate = real(diff_sum / sum_b, real32)
  end function check_rate

  real(real64) function mv_dense_parallel(repeat, bs, num_rows, x, matrix, y) result(time_ns)
    integer, intent(in) :: repeat, bs, num_rows
    real(real32), intent(in) :: x(0:)
    real(real32), intent(inout) :: matrix(0:)
    real(real32), intent(out) :: y(0:)
    integer :: n, i, j, num_elems
    real(real32) :: temp
    real(real64) :: start_time

    num_elems = num_rows * num_rows
    !$omp target data map(to: matrix(0:num_elems - 1), x(0:num_rows - 1)) map(from: y(0:num_rows - 1))
    start_time = omp_get_wtime()
    do n = 1, repeat
      !$omp target teams distribute parallel do num_threads(bs) private(i, j, temp)
      do i = 0, num_rows - 1
        temp = 0.0_real32
        do j = 0, num_rows - 1
          if (matrix(i * num_rows + j) /= 0.0_real32) then
            temp = temp + matrix(i * num_rows + j) * x(j)
          end if
        end do
        y(i) = temp
      end do
      !$omp end target teams distribute parallel do
    end do
    time_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    !$omp end target data
  end function mv_dense_parallel

  real(real64) function mv_csr_parallel(repeat, bs, num_rows, row_indices, col_indices, values, x, nnz, matrix, y) result(time_ns)
    integer, intent(in) :: repeat, bs, num_rows, nnz
    integer, intent(in) :: row_indices(0:), col_indices(0:)
    real(real32), intent(in) :: values(0:), x(0:)
    real(real32), intent(inout) :: matrix(0:)
    real(real32), intent(out) :: y(0:)
    integer :: n, i, j, row_start, row_end
    real(real32) :: temp
    real(real64) :: start_time

    !$omp target data map(to: row_indices(0:num_rows), col_indices(0:nnz - 1), values(0:nnz - 1), x(0:num_rows - 1)) &
    !$omp& map(from: y(0:num_rows - 1))
    start_time = omp_get_wtime()
    do n = 1, repeat
      !$omp target teams distribute parallel do num_threads(bs) private(i, j, row_start, row_end, temp)
      do i = 0, num_rows - 1
        row_start = row_indices(i)
        row_end = row_indices(i + 1)
        temp = 0.0_real32
        do j = row_start, row_end - 1
          temp = temp + values(j) * x(col_indices(j))
        end do
        y(i) = temp
      end do
      !$omp end target teams distribute parallel do
    end do
    time_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    !$omp end target data
  end function mv_csr_parallel

  integer function prev_power_of_2(v) result(value)
    integer, intent(in) :: v
    integer :: work

    work = v - 1
    work = ior(work, ishft(work, -1))
    work = ior(work, ishft(work, -2))
    work = ior(work, ishft(work, -4))
    work = ior(work, ishft(work, -8))
    work = ior(work, ishft(work, -16))
    work = work + 1
    value = ishft(work, -1)
  end function prev_power_of_2

  real(real64) function vector_mv_csr_parallel( &
      repeat, bs, num_rows, row_indices, col_indices, values, x, nnz, matrix, y) result(time_ns)
    integer, intent(in) :: repeat, bs, num_rows, nnz
    integer, intent(in) :: row_indices(0:), col_indices(0:)
    real(real32), intent(in) :: values(0:), x(0:)
    real(real32), intent(inout) :: matrix(0:)
    real(real32), intent(out) :: y(0:)
    integer :: nnz_per_row, threads_per_row, warp_size, rows_per_block, num_blocks
    integer :: n, i, j, row_start, row_end
    real(real32) :: temp
    real(real64) :: start_time

    nnz_per_row = nnz / num_rows
    threads_per_row = prev_power_of_2(nnz_per_row)
    warp_size = 32
    if (threads_per_row > warp_size) threads_per_row = warp_size
    rows_per_block = bs / threads_per_row
    if (rows_per_block == 0) rows_per_block = 1
    num_blocks = (num_rows + rows_per_block - 1) / rows_per_block

    !$omp target data map(to: row_indices(0:num_rows), col_indices(0:nnz - 1), values(0:nnz - 1), x(0:num_rows - 1)) &
    !$omp& map(from: y(0:num_rows - 1))
    start_time = omp_get_wtime()
    do n = 1, repeat
      !$omp target teams distribute num_teams(num_blocks * rows_per_block) private(i, j, row_start, row_end, temp)
      do i = 0, num_rows - 1
        row_start = row_indices(i)
        row_end = row_indices(i + 1)
        temp = 0.0_real32
        !$omp parallel do num_threads(threads_per_row) reduction(+:temp)
        do j = row_start, row_end - 1
          temp = temp + values(j) * x(col_indices(j))
        end do
        !$omp end parallel do
        y(i) = temp
      end do
      !$omp end target teams distribute
    end do
    time_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    !$omp end target data
  end function vector_mv_csr_parallel

end module simple_spmv_mod

program main
  use iso_c_binding, only: c_long
  use iso_fortran_env, only: int64, real32, real64
  use simple_spmv_mod
  implicit none

  integer :: argc, status, nnz, num_rows, repeat, num_elems
  integer :: vector_size, bs, i
  integer, allocatable :: row_indices(:), col_indices(:)
  real(real32), allocatable :: values(:), x(:), matrix(:), y0(:), y1(:), y2(:), y3(:)
  real(real64) :: elapsed(3), sparsity
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 3) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage ', trim(prog), ' <number of non-zero elements> <number of rows in a square matrix> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) nnz
  if (status /= 0) stop 1
  call get_command_argument(2, arg)
  read(arg, *, iostat=status) num_rows
  if (status /= 0) stop 1
  call get_command_argument(3, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  num_elems = num_rows * num_rows
  if (nnz <= 0 .or. num_rows <= 0 .or. nnz > num_elems) error stop "invalid benchmark dimensions"

  vector_size = num_rows
  allocate(row_indices(0:num_rows), col_indices(0:nnz - 1), values(0:nnz - 1))
  allocate(x(0:vector_size - 1), y0(0:vector_size - 1), y1(0:vector_size - 1))
  allocate(y2(0:vector_size - 1), y3(0:vector_size - 1), matrix(0:num_elems - 1))

  call c_srand48(int(ishft(1, 12), c_long))
  call init_matrix(matrix, num_rows, nnz)
  call init_vector(x, num_rows)
  call init_csr(row_indices, values, col_indices, matrix, num_rows, nnz)
  call mv_csr_serial(num_rows, row_indices, col_indices, values, x, y0)

  write(*,'(A,I0)') 'Number of non-zero elements: ', nnz
  write(*,'(A,I0)') 'Number of rows in a square matrix: ', num_rows
  sparsity = real(num_elems - nnz, real64) / real(num_elems, real64) * 100.0_real64
  write(*,'(A,F0.6,A)') 'Sparsity: ', sparsity, '%'

  do bs = 32, 1024, 32
    if (bs /= 32 .and. iand(bs, bs - 1) /= 0) cycle
    write(*,'(/,A,I0)') 'Thread block size: ', bs
    elapsed(1) = mv_dense_parallel(repeat, bs, num_rows, x, matrix, y1)
    elapsed(2) = mv_csr_parallel(repeat, bs, num_rows, row_indices, col_indices, values, x, nnz, matrix, y2)
    elapsed(3) = vector_mv_csr_parallel(repeat, bs, num_rows, row_indices, col_indices, values, x, nnz, matrix, y3)

    write(*,'(A)', advance='no') 'Average dense, sparse, and vector sparse kernel execution time (ms):'
    do i = 1, 3
      write(*,'(1X,F0.6)', advance='no') elapsed(i) * 1.0e-6_real64 / real(repeat, real64)
    end do
    write(*,*)

    write(*,'(A)', advance='no') 'Error rate:'
    write(*,'(1X,F8.6)', advance='no') check_rate(y0, y1, num_rows)
    write(*,'(1X,F8.6)', advance='no') check_rate(y0, y2, num_rows)
    write(*,'(1X,F8.6)', advance='no') check_rate(y0, y3, num_rows)
    write(*,*)
  end do

  deallocate(row_indices, col_indices, values, x, matrix, y0, y1, y2, y3)
end program main
