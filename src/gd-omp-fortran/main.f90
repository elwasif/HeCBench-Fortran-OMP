! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  type :: classification_data_crs
    integer(int64) :: nzmax = 0_int64
    integer :: m = 0
    integer :: n = 0
    real(real32), allocatable :: values(:)
    integer, allocatable :: col_index(:)
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: y_label(:)
  end type classification_data_crs

  type(classification_data_crs) :: a
  character(len=512) :: arg0, file_path, lambda_arg, alpha_arg, repeat_arg
  real(real32) :: lambda, alpha, obj_val, train_error
  real(real64) :: repeat_value
  real(real32), allocatable :: x(:), grad(:)
  integer :: iters, ios
  integer :: m, n, k, correct(1)
  real(real32) :: total_obj_val(1), l2_norm
  real(real64) :: train_start, train_end

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <path to file> <lambda> <alpha> <repeat>'
    stop 1
  end if

  call get_command_argument(1, file_path)
  call get_command_argument(2, lambda_arg)
  call get_command_argument(3, alpha_arg)
  call get_command_argument(4, repeat_arg)
  read(lambda_arg, *) lambda
  read(alpha_arg, *) alpha
  read(repeat_arg, *, iostat=ios) repeat_value
  if (ios /= 0) repeat_value = 0.0_real64
  iters = int(repeat_value)

  call get_crsm_from_svm(a, trim(file_path))
  m = a%m
  n = a%n

  allocate(x(n), grad(n))
  x = 0.0_real32
  grad = 0.0_real32

  !$omp target data map(to: a%row_ptr(1:size(a%row_ptr)), a%values(1:size(a%values)), &
  !$omp& a%col_index(1:size(a%col_index)), a%y_label(1:size(a%y_label))) &
  !$omp& map(tofrom: x(1:n)) map(alloc: grad(1:n), total_obj_val(1:1), correct(1:1))
  train_start = omp_get_wtime()

  do k = 1, iters
    total_obj_val(1) = 0.0_real32
    correct(1) = 0
    l2_norm = 0.0_real32

    !$omp target update to(total_obj_val(1:1))
    !$omp target update to(correct(1:1))

    grad = 0.0_real32
    !$omp target update to(grad(1:n))

    call compute_kernel(x, grad, a%row_ptr, a%col_index, a%values, a%y_label, &
      total_obj_val, correct, m)
    call update_kernel(x, grad, m, n, lambda, alpha, l2_norm)
  end do

  train_end = omp_get_wtime()
  write(*,'(A,F0.6,A,I0,A)') 'Training time takes ', train_end - train_start, &
    ' (s) for ', iters, ' iterations'
  write(*,*)

  !$omp target update from(total_obj_val(1:1))
  !$omp target update from(correct(1:1))

  obj_val = total_obj_val(1) / real(m, real32) + 0.5_real32 * lambda * l2_norm
  train_error = 1.0_real32 - (real(correct(1), real32) / real(m, real32))
  !$omp end target data

  write(*,'(A,F0.6,A,F0.6)') 'object value = ', obj_val, ' train_error = ', train_error

  call reference(a, n, iters, alpha, lambda, obj_val, train_error)

  deallocate(x, grad, a%values, a%col_index, a%row_ptr, a%y_label)

contains

  subroutine compute_kernel(x, grad, row_ptr, col_index, values, y_label, total_obj_val, correct, m)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: grad(:)
    integer, intent(in) :: row_ptr(:), col_index(:), y_label(:)
    real(real32), intent(in) :: values(:)
    real(real32), intent(inout) :: total_obj_val(:)
    integer, intent(inout) :: correct(:)
    integer, intent(in) :: m
    integer :: i, j, t
    real(real32) :: xp, v, prediction, accum, temp

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(i, j, t, xp, v, prediction, accum, temp)
    do i = 1, m
      xp = 0.0_real32
      do j = row_ptr(i), row_ptr(i + 1) - 1
        xp = xp + values(j) * x(col_index(j))
      end do

      v = log(1.0_real32 + exp(-xp * real(y_label(i), real32)))
      !$omp atomic update
      total_obj_val(1) = total_obj_val(1) + v

      prediction = 1.0_real32 / (1.0_real32 + exp(-xp))
      if (prediction >= 0.5_real32) then
        t = 1
      else
        t = -1
      end if
      if (y_label(i) == t) then
        !$omp atomic update
        correct(1) = correct(1) + 1
      end if

      accum = exp(-real(y_label(i), real32) * xp)
      accum = accum / (1.0_real32 + accum)
      do j = row_ptr(i), row_ptr(i + 1) - 1
        temp = -accum * values(j) * real(y_label(i), real32)
        !$omp atomic update
        grad(col_index(j)) = grad(col_index(j)) + temp
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_kernel

  subroutine update_kernel(x, grad, m, n, lambda, alpha, l2_norm)
    real(real32), intent(inout) :: x(:)
    real(real32), intent(in) :: grad(:)
    integer, intent(in) :: m, n
    real(real32), intent(in) :: lambda, alpha
    real(real32), intent(inout) :: l2_norm
    integer :: i
    real(real32) :: g

    !$omp target teams distribute parallel do reduction(+:l2_norm) thread_limit(256) private(i, g)
    do i = 1, n
      l2_norm = l2_norm + x(i) * x(i)
      g = grad(i) / real(m, real32) + lambda * x(i)
      x(i) = x(i) - alpha * g
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_kernel

  subroutine reference(a, n, iters, alpha, lambda, obj_val, train_error)
    type(classification_data_crs), intent(in) :: a
    integer, intent(in) :: n, iters
    real(real32), intent(in) :: alpha, lambda, obj_val, train_error
    real(real32), allocatable :: h_x(:), h_grad(:)
    real(real32) :: h_obj_val, h_train_error, total_obj_val, l2_norm
    integer :: correct, k
    logical :: ok

    allocate(h_x(n), h_grad(n))
    h_x = 0.0_real32

    do k = 1, iters
      total_obj_val = 0.0_real32
      l2_norm = 0.0_real32
      correct = 0
      h_grad = 0.0_real32

      call compute_ref(h_x, h_grad, a%row_ptr, a%col_index, a%values, a%y_label, &
        total_obj_val, correct, a%m)
      call l2_norm_ref(h_x, l2_norm, n)

      h_obj_val = total_obj_val / real(a%m, real32) + 0.5_real32 * lambda * l2_norm
      h_train_error = 1.0_real32 - (real(correct, real32) / real(a%m, real32))

      call update_ref(h_x, h_grad, a%m, n, lambda, alpha)
    end do

    ok = abs(obj_val - h_obj_val) < 1.0e-3_real32 .and. &
      abs(train_error - h_train_error) < 1.0e-3_real32
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if
    deallocate(h_x, h_grad)
  end subroutine reference

  subroutine compute_ref(x, grad, row_ptr, col_index, values, y_label, total_obj_val, correct, m)
    real(real32), intent(in) :: x(:)
    real(real32), intent(inout) :: grad(:)
    integer, intent(in) :: row_ptr(:), col_index(:), y_label(:)
    real(real32), intent(in) :: values(:)
    real(real32), intent(inout) :: total_obj_val
    integer, intent(inout) :: correct
    integer, intent(in) :: m
    integer :: i, j, t
    real(real32) :: xp, v, prediction, accum, temp

    do i = 1, m
      xp = 0.0_real32
      do j = row_ptr(i), row_ptr(i + 1) - 1
        xp = xp + values(j) * x(col_index(j))
      end do

      v = log(1.0_real32 + exp(-xp * real(y_label(i), real32)))
      total_obj_val = total_obj_val + v

      prediction = 1.0_real32 / (1.0_real32 + exp(-xp))
      if (prediction >= 0.5_real32) then
        t = 1
      else
        t = -1
      end if
      if (y_label(i) == t) correct = correct + 1

      accum = exp(-real(y_label(i), real32) * xp)
      accum = accum / (1.0_real32 + accum)
      do j = row_ptr(i), row_ptr(i + 1) - 1
        temp = -accum * values(j) * real(y_label(i), real32)
        grad(col_index(j)) = grad(col_index(j)) + temp
      end do
    end do
  end subroutine compute_ref

  subroutine l2_norm_ref(x, l2_norm, n)
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: l2_norm
    integer, intent(in) :: n
    integer :: i

    l2_norm = 0.0_real32
    do i = 1, n
      l2_norm = l2_norm + x(i) * x(i)
    end do
  end subroutine l2_norm_ref

  subroutine update_ref(x, grad, m, n, lambda, alpha)
    real(real32), intent(inout) :: x(:)
    real(real32), intent(in) :: grad(:)
    integer, intent(in) :: m, n
    real(real32), intent(in) :: lambda, alpha
    integer :: i
    real(real32) :: g

    do i = 1, n
      g = grad(i) / real(m, real32) + lambda * x(i)
      x(i) = x(i) - alpha * g
    end do
  end subroutine update_ref

  subroutine get_crsm_from_svm(a, file_path)
    type(classification_data_crs), intent(inout) :: a
    character(len=*), intent(in) :: file_path
    character(len=4096) :: line
    integer :: unit, ios, rows, nnz, max_col

    rows = 0
    nnz = 0
    max_col = 0
    open(newunit=unit, file=file_path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Could not find the SMV file, check again!'
      stop 1
    end if

    write(*,'(A)') 'Processing the SVM file'
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call scan_svm_line(trim(line), rows, nnz, max_col, a, .false.)
    end do
    close(unit)

    a%m = rows
    a%n = max_col
    a%nzmax = int(nnz, int64)
    allocate(a%row_ptr(rows + 1), a%col_index(nnz), a%values(nnz), a%y_label(rows))
    a%row_ptr(1) = 1

    rows = 0
    nnz = 0
    open(newunit=unit, file=file_path, status='old', action='read', iostat=ios)
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call scan_svm_line(trim(line), rows, nnz, max_col, a, .true.)
    end do
    close(unit)

    call normalize_rows(a)
    write(*,'(A,I0,A,I0,A,I0)') 'Finished processing the LIBSVM file. ', a%m, &
      ' observations and ', a%n, ' features were read. The total number of non-zero elements are: ', a%nzmax
  end subroutine get_crsm_from_svm

  subroutine scan_svm_line(line, rows, nnz, max_col, a, fill)
    character(len=*), intent(in) :: line
    integer, intent(inout) :: rows, nnz, max_col
    type(classification_data_crs), intent(inout) :: a
    logical, intent(in) :: fill
    character(len=256) :: token
    integer :: pos, next, colon, label, col
    real(real32) :: value

    if (len_trim(line) == 0) return
    rows = rows + 1
    pos = 1
    call next_token(line, pos, token)
    read(token, *) label
    if (fill) a%y_label(rows) = label

    do
      call next_token(line, pos, token)
      if (len_trim(token) == 0) exit
      colon = index(token, ':')
      if (colon <= 1) cycle
      read(token(:colon - 1), *) col
      read(token(colon + 1:), *) value
      nnz = nnz + 1
      max_col = max(max_col, col)
      if (fill) then
        next = nnz
        a%col_index(next) = col
        a%values(next) = value
      end if
    end do
    if (fill) a%row_ptr(rows + 1) = nnz + 1
  end subroutine scan_svm_line

  subroutine next_token(line, pos, token)
    character(len=*), intent(in) :: line
    integer, intent(inout) :: pos
    character(len=*), intent(out) :: token
    integer :: start, finish, n

    token = ''
    n = len_trim(line)
    do while (pos <= n .and. (line(pos:pos) == ' ' .or. line(pos:pos) == char(9)))
      pos = pos + 1
    end do
    if (pos > n) return

    start = pos
    do while (pos <= n .and. line(pos:pos) /= ' ' .and. line(pos:pos) /= char(9))
      pos = pos + 1
    end do
    finish = pos - 1
    token = line(start:finish)
  end subroutine next_token

  subroutine normalize_rows(a)
    type(classification_data_crs), intent(inout) :: a
    integer :: i, j
    real(real32) :: norm_sqrd, norm

    do i = 1, a%m
      norm_sqrd = 0.0_real32
      do j = a%row_ptr(i), a%row_ptr(i + 1) - 1
        norm_sqrd = norm_sqrd + a%values(j) ** 2
      end do
      norm = sqrt(norm_sqrd)
      if (norm > 0.0_real32) then
        do j = a%row_ptr(i), a%row_ptr(i + 1) - 1
          a%values(j) = a%values(j) / norm
        end do
      end if
    end do
  end subroutine normalize_rows

end program main
