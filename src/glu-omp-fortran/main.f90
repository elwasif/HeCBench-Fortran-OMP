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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: ndims, dim_size, repeat
  integer, allocatable :: xshape(:), yshape(:)
  integer(int64) :: nelems
  real(real32), allocatable :: x(:), y(:), y_ref(:)
  integer :: i, input_dim, split_index, split_dim_size, m, n
  real(real64) :: start_time, end_time, elapsed_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <number of dimensions> <size of each dimension> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) ndims
  read(arg2, *) dim_size
  read(arg3, *) repeat

  allocate(xshape(ndims), yshape(ndims))
  xshape = dim_size
  yshape = dim_size

  write(*, '(A)', advance='no') 'Shape of input tensor: ( '
  do i = 1, ndims
    write(*, '(I0,A)', advance='no') xshape(i), ' '
  end do
  print '(A)', ')'

  nelems = size_from_dim(1, xshape)
  allocate(x(nelems), y(nelems), y_ref(nelems))

  call initialize_input(x)
  y = 0.0_real32
  y_ref = 0.0_real32

  !$omp target data map(to: x(1:nelems)) map(from: y(1:nelems))
  do input_dim = -1, 3 * (ndims - 1) - 1
    if (input_dim == -1) then
      split_index = ndims
    else
      split_index = mod(input_dim, ndims) + 1
    end if

    if (mod(yshape(split_index), 2) /= 0) then
      print '(A,I0,A)', 'Split dimension ', yshape(split_index), ' should be divided by two. Skip'
      cycle
    end if

    split_dim_size = yshape(split_index) / 2
    m = int(size_to_dim(split_index, xshape))
    n = int(size_from_dim(split_index + 1, xshape))

    call compute_glu_ref(m, split_dim_size, n, x, y_ref)

    start_time = omp_get_wtime()
    do i = 1, repeat
      call glu_kernel(m, split_dim_size, n, x, y)
    end do
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A,I0,A,F0.6,A)', 'Average execution time of GLU kernel (split dimension = ', &
        split_index - 1, '): ', elapsed_us, ' (us)'

    !$omp target update from(y(1:nelems))

    ok = .true.
    do i = 1, int(nelems / 2_int64)
      if (abs(y(i) - y_ref(i)) > 1.0e-3_real32) then
        ok = .false.
        exit
      end if
    end do

    if (ok) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end do
  !$omp end target data

  deallocate(xshape, yshape, x, y, y_ref)

contains

  subroutine initialize_input(x)
    real(real32), intent(out) :: x(:)
    integer(c_int), parameter :: rand_max = huge(0_c_int)
    integer :: idx

    call c_srand(123_c_int)
    do idx = 1, size(x)
      x(idx) = real(12.0_real64 * real(c_rand(), real64) / real(rand_max, real64) - 6.0_real64, real32)
    end do
  end subroutine initialize_input

  integer(int64) function size_from_dim(k, dims)
    integer, intent(in) :: k
    integer, intent(in) :: dims(:)
    integer :: idx

    size_from_dim = 1_int64
    do idx = k, size(dims)
      size_from_dim = size_from_dim * int(dims(idx), int64)
    end do
  end function size_from_dim

  integer(int64) function size_to_dim(k, dims)
    integer, intent(in) :: k
    integer, intent(in) :: dims(:)
    integer :: idx

    size_to_dim = 1_int64
    do idx = 1, k - 1
      size_to_dim = size_to_dim * int(dims(idx), int64)
    end do
  end function size_to_dim

  real(real32) function sigmoid(x)
    real(real32), intent(in) :: x
    real(real32) :: exp_x

    if (x >= 0.0_real32) then
      sigmoid = 1.0_real32 / (1.0_real32 + exp(-x))
    else
      exp_x = exp(x)
      sigmoid = exp_x / (1.0_real32 + exp_x)
    end if
  end function sigmoid

  subroutine compute_glu_ref(m, split_dim_size, n, x, y)
    integer, intent(in) :: m, split_dim_size, n
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: y(:)
    integer :: i, j, k, idx, idy, jn, jdx1, jdx2, jdy, x_stride, y_stride

    y_stride = split_dim_size * n
    x_stride = 2 * y_stride
    do i = 0, m - 1
      idx = i * x_stride
      idy = i * y_stride
      do j = 0, split_dim_size - 1
        jn = j * n
        jdx1 = idx + jn
        jdx2 = idx + (j + split_dim_size) * n
        jdy = idy + jn
        do k = 0, n - 1
          y(jdy + k + 1) = x(jdx1 + k + 1) * sigmoid(x(jdx2 + k + 1))
        end do
      end do
    end do
  end subroutine compute_glu_ref

  subroutine glu_kernel(m, split_dim_size, n, x, y)
    integer, intent(in) :: m, split_dim_size, n
    real(real32), intent(in) :: x(:)
    real(real32), intent(out) :: y(:)
    integer :: index, i, j, k, x_offset, y_offset
    real(real32) :: x1, x2

    x_offset = 2 * split_dim_size * n
    y_offset = split_dim_size * n

    !$omp target teams distribute parallel do thread_limit(256) private(i, j, k, x1, x2)
    do index = 0, m * split_dim_size * n - 1
      i = index / (split_dim_size * n)
      j = mod(index / n, split_dim_size)
      k = mod(index, n)
      x1 = x(i * x_offset + j * n + k + 1)
      x2 = x(i * x_offset + (j + split_dim_size) * n + k + 1)
      y(i * y_offset + j * n + k + 1) = x1 * (1.0_real32 / (1.0_real32 + exp(-x2)))
    end do
    !$omp end target teams distribute parallel do
  end subroutine glu_kernel

end program main
