! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_char, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only : int64, real64
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

    function c_atoi(str) bind(C, name="atoi") result(value)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: str(*)
      integer(c_int) :: value
    end function c_atoi
  end interface

  character(len=256) :: arg0, arg
  integer :: nelem, repeat, i
  integer, allocatable :: grad_in(:), grad_out(:)
  integer(int64) :: idx_dim, size, step, grad_in_dim_stride
  integer(int64) :: grad_in_last_dim_stride, grad_in_dim_size
  integer(int64) :: grad_out_dim_stride
  real(real64) :: start_time, elapsed
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  nelem = atoi_argument(arg)
  call get_command_argument(2, arg)
  repeat = atoi_argument(arg)

  size = 2_int64
  step = 1_int64
  grad_in_dim_stride = 1_int64
  grad_in_last_dim_stride = 1_int64
  grad_in_dim_size = int(nelem, int64)
  grad_out_dim_stride = 1_int64
  allocate(grad_in(nelem), grad_out(nelem))

  call c_srand(123_c_int)
  do i = 1, nelem
    grad_in(i) = modulo(c_rand(), 256_c_int)
  end do

  idx_dim = 0_int64

  !$omp target data map(to: grad_in(1:nelem), idx_dim) map(from: grad_out(1:nelem))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call unfold_backward_internal_kernel(grad_out, grad_in, idx_dim, size, step, &
      grad_in_dim_stride, grad_in_last_dim_stride, grad_in_dim_size, &
      grad_out_dim_stride)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average execution time of unfold backward kernel: ', &
    (elapsed * 1000000.0_real64) / real(repeat, real64), ' (us)'

  ok = .true.
  do i = 1, nelem
    if (repeat * grad_in(i) /= grad_out(i)) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(grad_in, grad_out)

contains

  integer function atoi_argument(value) result(parsed)
    character(len=*), intent(in) :: value
    character(kind=c_char), allocatable :: c_value(:)
    integer :: pos, trimmed_len

    trimmed_len = len_trim(value)
    allocate(c_value(trimmed_len + 1))
    do pos = 1, trimmed_len
      c_value(pos) = value(pos:pos)
    end do
    c_value(trimmed_len + 1) = c_null_char
    parsed = int(c_atoi(c_value), kind(parsed))
    deallocate(c_value)
  end function atoi_argument

  subroutine unfold_backward_internal_kernel(grad_out, grad_in, idx_dim, size, step, &
      grad_in_dim_stride, grad_in_last_dim_stride, grad_in_dim_size, &
      grad_out_dim_stride)
    integer, intent(inout) :: grad_out(:)
    integer, intent(in) :: grad_in(:)
    integer(int64), intent(in) :: idx_dim, size, step
    integer(int64), intent(in) :: grad_in_dim_stride, grad_in_last_dim_stride
    integer(int64), intent(in) :: grad_in_dim_size, grad_out_dim_stride
    integer, parameter :: n_threads = 64
    integer, parameter :: n_elems_per_thread = 4
    integer, parameter :: total_work_block = n_threads * n_elems_per_thread
    integer :: grid

    grid = int((grad_in_dim_size + total_work_block - 1_int64) / total_work_block)
    call unfold_backward_elementwise_kernel(grid, int(grad_in_dim_size), grad_out, &
      grad_in, idx_dim, size, step, grad_in_dim_stride, grad_in_last_dim_stride, &
      grad_in_dim_size, grad_out_dim_stride)
  end subroutine unfold_backward_internal_kernel

  subroutine unfold_backward_elementwise_kernel(grid, total_n_elems, grad_out, &
      grad_in, idx_dim, size, step, grad_in_dim_stride, grad_in_last_dim_stride, &
      grad_in_dim_size, grad_out_dim_stride)
    integer, intent(in) :: grid, total_n_elems
    integer, intent(inout) :: grad_out(:)
    integer, intent(in) :: grad_in(:)
    integer(int64), intent(in) :: idx_dim, size, step
    integer(int64), intent(in) :: grad_in_dim_stride, grad_in_last_dim_stride
    integer(int64), intent(in) :: grad_in_dim_size, grad_out_dim_stride
    integer, parameter :: n_threads = 64
    integer, parameter :: n_elems_per_thread = 4
    integer, parameter :: total_work_block = n_threads * n_elems_per_thread
    integer :: bid, tid, elem, j
    integer(int64) :: left_fold_idx, right_fold_idx, fold_idx, idx_last_dim
    integer(int64) :: grad_out_index, grad_in_index

    !$omp target teams distribute parallel do collapse(2) num_teams(grid) thread_limit(n_threads) &
    !$omp& private(elem, j, left_fold_idx, right_fold_idx, fold_idx, idx_last_dim, &
    !$omp& grad_out_index, grad_in_index)
    do bid = 0, grid - 1
      do tid = 0, n_threads - 1
        elem = bid * total_work_block + tid
        do j = 1, n_elems_per_thread
          if (elem < total_n_elems) then
            grad_out_index = int(elem, int64) * grad_out_dim_stride
            left_fold_idx = 0_int64
            if (idx_dim > size) left_fold_idx = (idx_dim - size) / step
            if (.not. (left_fold_idx * step <= idx_dim .and. &
                idx_dim < left_fold_idx * step + size)) then
              left_fold_idx = left_fold_idx + 1_int64
            end if

            right_fold_idx = idx_dim / step
            if (right_fold_idx >= grad_in_dim_size) then
              right_fold_idx = grad_in_dim_size - 1_int64
            end if

            do fold_idx = left_fold_idx, right_fold_idx
              idx_last_dim = idx_dim - fold_idx * step
              grad_in_index = int(elem, int64) + fold_idx * grad_in_dim_stride + &
                idx_last_dim * grad_in_last_dim_stride
              grad_out(grad_out_index + 1_int64) = grad_out(grad_out_index + 1_int64) + &
                grad_in(grad_in_index + 1_int64)
            end do
            elem = elem + n_threads
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine unfold_backward_elementwise_kernel

end program main
