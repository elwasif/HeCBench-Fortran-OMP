! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: GPU_THREADS = 256
  integer, parameter :: fill_val = -1
  integer :: m, n, b, repeat

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name='rand')
      import :: c_int
    end function c_rand
  end interface

  if (command_argument_count() /= 4) then
    print '(A)', 'Usage: ./main <sequence length> <sequence length> <batch size> <repeat>'
    stop 1
  end if

  m = read_arg(1)
  n = read_arg(2)
  b = read_arg(3)
  repeat = read_arg(4)
  if (m <= 0 .or. n <= 0 .or. repeat <= 0) then
    print '(A)', 'invalid arguments'
    stop 1
  end if

  call eval_mask(m, n, b, repeat)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine eval_mask(m, n, b, repeat)
    integer, intent(in) :: m, n, b, repeat
    integer :: batch_dim, data_size, radius, idx
    integer(int32), allocatable :: h_in(:), h_out(:), out_ref(:), seq_len(:), window(:)

    batch_dim = merge(1, b, b <= 0)
    radius = m / 4
    data_size = n * m * batch_dim
    print *
    write(*, '(A,I0,A,I0,A,I0)') 'M = ', m, ', N = ', n, ', B = ', batch_dim

    allocate(h_in(data_size), h_out(data_size), out_ref(data_size), seq_len(n), window(n))
    call c_srand(123_c_int)
    do idx = 1, n
      seq_len(idx) = int(mod(c_rand(), max(1, m / 2)), int32)
    end do
    do idx = 1, n
      window(idx) = int(mod(c_rand(), m), int32)
    end do
    do idx = 1, data_size
      h_in(idx) = int(mod(c_rand(), m * n), int32)
    end do

    !$omp target data map(to: h_in(1:data_size), seq_len(1:n), window(1:n)) map(alloc: h_out(1:data_size))
      call run_one_mask('sequenceMask', 1, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
      call run_one_mask('windowMask', 2, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
      call run_one_mask('upperMask', 3, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
      call run_one_mask('lowerMask', 4, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
      call run_one_mask('upperDiagMask', 5, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
      call run_one_mask('lowerDiagMask', 6, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
    !$omp end target data

    deallocate(h_in, h_out, out_ref, seq_len, window)
  end subroutine eval_mask

  subroutine run_one_mask(name, mode, m, n, b, batch_dim, repeat, radius, h_in, h_out, out_ref, seq_len, window, data_size)
    character(len=*), intent(in) :: name
    integer, intent(in) :: mode, m, n, b, batch_dim, repeat, radius, data_size
    integer(int32), intent(in) :: h_in(:), seq_len(:), window(:)
    integer(int32), intent(inout) :: h_out(:)
    integer(int32), intent(out) :: out_ref(:)
    integer :: r
    real(real64) :: start_time, end_time

    call mask_reference(mode, m, n, b, batch_dim, radius, h_in, seq_len, window, out_ref)
    start_time = omp_get_wtime()
    do r = 1, repeat
      select case (mode)
      case (1)
        call sequenceMaskKernel(n, m, batch_dim, h_in, seq_len, fill_val, h_out)
      case (2)
        call windowMaskKernel(n, m, batch_dim, h_in, window, radius, fill_val, h_out)
      case (3)
        call upperMaskKernel(n, m, batch_dim, h_in, fill_val, h_out)
      case (4)
        call lowerMaskKernel(n, m, batch_dim, h_in, fill_val, h_out)
      case (5)
        call upperDiagMaskKernel(n, m, batch_dim, h_in, fill_val, h_out)
      case default
        call lowerDiagMaskKernel(n, m, batch_dim, h_in, fill_val, h_out)
      end select
    end do
    end_time = omp_get_wtime()
    write(*, '(A,A,A,F0.6,A)') 'Average execution time of ', trim(name), ' kernel: ', &
      ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64), ' (us)'
    !$omp target update from(h_out(1:data_size))
    call print_mask_ratio(h_out, out_ref, data_size)
  end subroutine run_one_mask

  subroutine sequenceMaskKernel(n, m, b, h_in, seq_lengths, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, fill_val_arg
    integer(int32), intent(in) :: h_in(:), seq_lengths(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), k >= seq_lengths(j + 1))
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), j >= seq_lengths(i + 1))
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine sequenceMaskKernel

  subroutine windowMaskKernel(n, m, b, h_in, window_centers, radius, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, radius, fill_val_arg
    integer(int32), intent(in) :: h_in(:), window_centers(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), &
          k < window_centers(j + 1) - radius .or. k > window_centers(j + 1) + radius)
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), &
          j < window_centers(i + 1) - radius .or. j > window_centers(i + 1) + radius)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine windowMaskKernel

  subroutine upperMaskKernel(n, m, b, h_in, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, fill_val_arg
    integer(int32), intent(in) :: h_in(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), k > j)
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), j > i)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine upperMaskKernel

  subroutine lowerMaskKernel(n, m, b, h_in, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, fill_val_arg
    integer(int32), intent(in) :: h_in(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), k < j)
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), j < i)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine lowerMaskKernel

  subroutine upperDiagMaskKernel(n, m, b, h_in, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, fill_val_arg
    integer(int32), intent(in) :: h_in(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), k >= j)
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), j >= i)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine upperDiagMaskKernel

  subroutine lowerDiagMaskKernel(n, m, b, h_in, fill_val_arg, h_out)
    integer, intent(in) :: n, m, b, fill_val_arg
    integer(int32), intent(in) :: h_in(:)
    integer(int32), intent(inout) :: h_out(:)
    integer :: index, i, j, k, ind

    if (b >= 0) then
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n) thread_limit(GPU_THREADS) private(index, i, j, k, ind)
      do index = 0, b * n * m - 1
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        h_out(ind) = merge(fill_val_arg, h_in(ind), k <= j)
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do &
      !$omp& num_teams(m * n / GPU_THREADS) thread_limit(GPU_THREADS) private(index, i, j)
      do index = 0, n * m - 1
        i = index / m
        j = mod(index, m)
        h_out(index + 1) = merge(fill_val_arg, h_in(index + 1), j <= i)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine lowerDiagMaskKernel

  subroutine mask_reference(mode, m, n, b, batch_dim, radius, h_in, seq_len, window, out_ref)
    integer, intent(in) :: mode, m, n, b, batch_dim, radius
    integer(int32), intent(in) :: h_in(:), seq_len(:), window(:)
    integer(int32), intent(out) :: out_ref(:)
    integer :: index, i, j, k, ind, total

    total = n * m * batch_dim
    do index = 0, total - 1
      if (b >= 0) then
        k = mod(index, m)
        j = mod((index - k) / m, n)
        i = (index - m * j - k) / (n * m)
        ind = n * m * i + m * j + k + 1
        out_ref(ind) = masked_value(mode, k, j, h_in(ind), seq_len(j + 1), window(j + 1), radius)
      else
        i = index / m
        j = mod(index, m)
        ind = index + 1
        out_ref(ind) = masked_value(mode, j, i, h_in(ind), seq_len(i + 1), window(i + 1), radius)
      end if
    end do
  end subroutine mask_reference

  integer(int32) function masked_value(mode, col, row, in_value, seq_value, window_value, radius)
    integer, intent(in) :: mode, col, row, seq_value, window_value, radius
    integer(int32), intent(in) :: in_value
    select case (mode)
    case (1)
      masked_value = merge(fill_val, in_value, col >= seq_value)
    case (2)
      masked_value = merge(fill_val, in_value, col < window_value - radius .or. col > window_value + radius)
    case (3)
      masked_value = merge(fill_val, in_value, col > row)
    case (4)
      masked_value = merge(fill_val, in_value, col < row)
    case (5)
      masked_value = merge(fill_val, in_value, col >= row)
    case default
      masked_value = merge(fill_val, in_value, col <= row)
    end select
  end function masked_value

  subroutine print_mask_ratio(h_out, out_ref, data_size)
    integer(int32), intent(in) :: h_out(:), out_ref(:)
    integer, intent(in) :: data_size
    integer :: idx, cnt_fill, errors
    cnt_fill = 0
    errors = 0
    do idx = 1, data_size
      if (h_out(idx) == fill_val) cnt_fill = cnt_fill + 1
      if (h_out(idx) /= out_ref(idx)) errors = errors + 1
    end do
    if (errors == 0) then
      write(*, '(A,F8.6)') 'PASS, Mask ratio: ', real(cnt_fill, real64) / real(data_size, real64)
    else
      write(*, '(A,F8.6)') 'FAIL, Mask ratio: ', real(cnt_fill, real64) / real(data_size, real64)
    end if
  end subroutine print_mask_ratio

end program main
