! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
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
  end interface

  integer(int32), parameter :: num_threads = 256_int32

  character(len=256) :: arg0, arg1, arg2, arg3, arg4
  integer(int32) :: nrows, ndims, top_k, repeat, data_size
  integer(int32), allocatable :: label(:)
  real(real32), allocatable :: data(:)
  integer(int32) :: count_ref, ngrid, iter, count(0:0)
  integer(int32) :: row, col, label_data, ngt
  real(real32) :: label_pred, pred
  real(real64) :: start_time, end_time, avg_us

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of rows> <number of columns> <top K> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  read(arg1, *) nrows
  read(arg2, *) ndims
  read(arg3, *) top_k
  read(arg4, *) repeat
  if (nrows <= 0_int32 .or. ndims <= 0_int32 .or. top_k <= 0_int32 .or. repeat <= 0_int32) stop 1

  data_size = nrows * ndims
  allocate(label(0:nrows - 1), data(0:data_size - 1))
  call initialize_inputs(nrows, ndims, label, data)

  count_ref = reference_count(nrows, ndims, top_k, data, label)

  !$omp target data map(to: label(0:nrows - 1), data(0:data_size - 1)) map(alloc: count(0:0))
  do ngrid = nrows / 4, nrows, nrows / 4
    print '(A,I0)', 'Grid size is ', ngrid
    start_time = omp_get_wtime()
    do iter = 1, repeat
      count(0) = 0_int32
      !$omp target update to(count(0:0))

      !$omp target teams distribute num_teams(ngrid) private(label_data, label_pred, ngt, col, pred)
      do row = 0, nrows - 1
        label_data = label(row)
        label_pred = data(row * ndims + label_data)
        ngt = 0_int32
        !$omp parallel do reduction(+:ngt) num_threads(num_threads) private(pred)
        do col = 0, ndims - 1
          pred = data(row * ndims + col)
          if (pred > label_pred .or. (pred == label_pred .and. col <= label_data)) ngt = ngt + 1_int32
        end do
        !$omp end parallel do
        if (ngt <= top_k) then
          !$omp atomic update
          count(0) = count(0) + 1_int32
        end if
      end do
      !$omp end target teams distribute
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0d6) / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of accuracy kernel: ', avg_us, ' (us)'
    !$omp target update from(count(0:0))
    if (count(0) == count_ref) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
      stop 1
    end if
  end do
  !$omp end target data

  deallocate(label, data)

contains

  subroutine initialize_inputs(rows, dims, label, data)
    integer(int32), intent(in) :: rows, dims
    integer(int32), intent(out) :: label(0:)
    real(real32), intent(out) :: data(0:)
    integer(int32) :: row, idx
    integer(int64) :: g

    call c_srand(123_c_int)
    do row = 0, rows - 1
      label(row) = int(mod(c_rand(), int(dims, c_int)), int32)
    end do

    g = 123_int64
    do idx = 0, rows * dims - 1
      data(idx) = uniform_real_distribution(g)
    end do
  end subroutine initialize_inputs

  real(real32) function uniform_real_distribution(g) result(value)
    integer(int64), intent(inout) :: g
    integer(int64), parameter :: a = 16807_int64
    integer(int64), parameter :: m = 2147483647_int64
    integer(int64), parameter :: range = 2147483646_int64

    g = mod(a * g, m)
    value = real(real(g - 1_int64, real64) / real(range, real64), real32)
  end function uniform_real_distribution

  integer(int32) function reference_count(rows, dims, topk, data, label) result(count)
    integer(int32), intent(in) :: rows, dims, topk
    real(real32), intent(in) :: data(0:)
    integer(int32), intent(in) :: label(0:)
    integer(int32) :: row, col, label_data, ngt
    real(real32) :: label_pred, pred

    count = 0_int32
    do row = 0, rows - 1
      label_data = label(row)
      label_pred = data(row * dims + label_data)
      ngt = 0_int32
      do col = 0, dims - 1
        pred = data(row * dims + col)
        if (pred > label_pred .or. (pred == label_pred .and. col <= label_data)) ngt = ngt + 1_int32
      end do
      if (ngt <= topk) count = count + 1_int32
    end do
  end function reference_count

end program main
