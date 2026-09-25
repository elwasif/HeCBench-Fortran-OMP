! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: n_contractions = 18
  character(len=256) :: arg
  integer :: max_n, max_channels, repeat

  if (command_argument_count() /= 2) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <dimension> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) max_n
  call get_command_argument(2, arg)
  read(arg, *) repeat
  if (max_n <= 0 .or. repeat <= 0) stop 1

  max_channels = n_contractions
  call contract_real32(max_n, max_channels, repeat)
  call contract_real64(max_n, max_channels, repeat)

contains

  subroutine contract_real32(max_n, max_channels, repeat)
    integer, intent(in) :: max_n, max_channels, repeat
    integer(int64) :: tensor_size, adj_size, output_size
    real(real32), allocatable :: tensor_value(:), adj_value(:), value(:)
    real(real64) :: start_time, end_time, checksum, min_value, max_value
    integer(int64) :: i
    integer :: iter

    tensor_size = int(max_n, int64) * max_n * max_n * max_channels
    adj_size = int(max_n, int64) * max_n
    output_size = int(max_n, int64) * max_n * max_channels * n_contractions

    allocate(tensor_value(0:tensor_size - 1), adj_value(0:adj_size - 1), value(0:output_size - 1))
    do i = 0_int64, tensor_size - 1_int64
      tensor_value(i) = 1.0_real32
    end do
    do i = 0_int64, adj_size - 1_int64
      adj_value(i) = 1.0_real32
    end do

    !$omp target data map(to: tensor_value(0:tensor_size - 1), adj_value(0:adj_size - 1)) &
    !$omp& map(from: value(0:output_size - 1))
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call contraction_real32(tensor_value, adj_value, value, output_size, max_n, max_channels)
    end do
    end_time = omp_get_wtime()
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', &
      (end_time - start_time) / real(repeat, real64), ' (s)'
    !$omp end target data

    call summarize_real32(value, output_size, checksum, min_value, max_value)
    write(*,'(A,F0.6,A,F0.6,A,F0.6)') 'Checksum: ', checksum, ' min:', min_value, ' max:', max_value
    deallocate(value, tensor_value, adj_value)
  end subroutine contract_real32

  subroutine contract_real64(max_n, max_channels, repeat)
    integer, intent(in) :: max_n, max_channels, repeat
    integer(int64) :: tensor_size, adj_size, output_size
    real(real64), allocatable :: tensor_value(:), adj_value(:), value(:)
    real(real64) :: start_time, end_time, checksum, min_value, max_value
    integer(int64) :: i
    integer :: iter

    tensor_size = int(max_n, int64) * max_n * max_n * max_channels
    adj_size = int(max_n, int64) * max_n
    output_size = int(max_n, int64) * max_n * max_channels * n_contractions

    allocate(tensor_value(0:tensor_size - 1), adj_value(0:adj_size - 1), value(0:output_size - 1))
    do i = 0_int64, tensor_size - 1_int64
      tensor_value(i) = 1.0_real64
    end do
    do i = 0_int64, adj_size - 1_int64
      adj_value(i) = 1.0_real64
    end do

    !$omp target data map(to: tensor_value(0:tensor_size - 1), adj_value(0:adj_size - 1)) &
    !$omp& map(from: value(0:output_size - 1))
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call contraction_real64(tensor_value, adj_value, value, output_size, max_n, max_channels)
    end do
    end_time = omp_get_wtime()
    write(*,'(A,F0.6,A)') 'Average kernel execution time ', &
      (end_time - start_time) / real(repeat, real64), ' (s)'
    !$omp end target data

    call summarize_real64(value, output_size, checksum, min_value, max_value)
    write(*,'(A,F0.6,A,F0.6,A,F0.6)') 'Checksum: ', checksum, ' min:', min_value, ' max:', max_value
    deallocate(value, tensor_value, adj_value)
  end subroutine contract_real64

  subroutine contraction_real32(tensor, adj, value, output_size, n, channels)
    real(real32), intent(in) :: tensor(0:), adj(0:)
    real(real32), intent(out) :: value(0:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: n, channels
    integer(int64) :: tid, a_stride, b_stride, c_stride, y_size
    integer :: f, case_id, x, y, a, b, c, d, e
    real(real32) :: sum_value, adj_value

    c_stride = int(channels, int64)
    b_stride = int(n, int64) * c_stride
    a_stride = int(n, int64) * b_stride
    y_size = int(channels * n_contractions, int64)

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(tid, f, case_id, x, y, a, b, c, d, e, sum_value, adj_value)
    do tid = 0_int64, output_size - 1_int64
      f = int(mod(mod(tid, y_size), int(channels, int64)))
      case_id = int(mod(tid, y_size) / int(channels, int64)) + 1
      y = int(mod(tid / y_size, int(n, int64)))
      x = int((tid / y_size) / int(n, int64))
      sum_value = 0.0_real32

      select case (case_id)
      case (1)
        a = x
        b = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real32) then
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end if
          end do
        end do
      case (2)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            do b = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end if
        end do
      case (3)
        b = x
        c = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real32) then
              do a = 0, n - 1
                sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end if
          end do
        end do
      case (4)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            do a = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end if
        end do
      case (5)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real32) then
          do a = 0, n - 1
            do b = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end do
        end if
      case (6)
        a = x
        b = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            c = d
            sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end do
        end do
      case (7)
        a = x
        b = y
        do d = 0, n - 1
          e = d
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (8)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            do b = 0, n - 1
              c = b
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (9)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            b = e
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (10)
        b = x
        c = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real32) then
              a = d
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end if
          end do
        end do
      case (11)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            do a = 0, n - 1
              c = a
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (12)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            a = e
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (13)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            c = e
            do a = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (14)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real32) then
          do a = 0, n - 1
            b = a
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end do
        end if
      case (15)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real32) then
          do b = 0, n - 1
            c = b
            do a = 0, n - 1
              sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end do
        end if
      case (16)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            b = e
            c = e
            sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end if
        end do
      case (17)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real32) then
            a = e
            c = e
            sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end if
        end do
      case (18)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real32) then
          do a = 0, n - 1
            b = a
            c = a
            sum_value = sum_value + tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end do
        end if
      end select

      value(tid) = sum_value
    end do
    !$omp end target teams distribute parallel do
  end subroutine contraction_real32

  subroutine contraction_real64(tensor, adj, value, output_size, n, channels)
    real(real64), intent(in) :: tensor(0:), adj(0:)
    real(real64), intent(out) :: value(0:)
    integer(int64), intent(in) :: output_size
    integer, intent(in) :: n, channels
    integer(int64) :: tid, a_stride, b_stride, c_stride, y_size
    integer :: f, case_id, x, y, a, b, c, d, e
    real(real64) :: sum_value, adj_value

    c_stride = int(channels, int64)
    b_stride = int(n, int64) * c_stride
    a_stride = int(n, int64) * b_stride
    y_size = int(channels * n_contractions, int64)

    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(tid, f, case_id, x, y, a, b, c, d, e, sum_value, adj_value)
    do tid = 0_int64, output_size - 1_int64
      f = int(mod(mod(tid, y_size), int(channels, int64)))
      case_id = int(mod(tid, y_size) / int(channels, int64)) + 1
      y = int(mod(tid / y_size, int(n, int64)))
      x = int((tid / y_size) / int(n, int64))
      sum_value = 0.0_real64

      select case (case_id)
      case (1)
        a = x
        b = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real64) then
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end if
          end do
        end do
      case (2)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            do b = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end if
        end do
      case (3)
        b = x
        c = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real64) then
              do a = 0, n - 1
                sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end if
          end do
        end do
      case (4)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            do a = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end if
        end do
      case (5)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real64) then
          do a = 0, n - 1
            do b = 0, n - 1
              do c = 0, n - 1
                sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
              end do
            end do
          end do
        end if
      case (6)
        a = x
        b = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            c = d
            sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end do
        end do
      case (7)
        a = x
        b = y
        do d = 0, n - 1
          e = d
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (8)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            do b = 0, n - 1
              c = b
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (9)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            b = e
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (10)
        b = x
        c = y
        do d = 0, n - 1
          do e = 0, n - 1
            adj_value = adj(d * n + e)
            if (adj_value > 0.0_real64) then
              a = d
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end if
          end do
        end do
      case (11)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            do a = 0, n - 1
              c = a
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (12)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            a = e
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (13)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            c = e
            do a = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end if
        end do
      case (14)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real64) then
          do a = 0, n - 1
            b = a
            do c = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end do
        end if
      case (15)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real64) then
          do b = 0, n - 1
            c = b
            do a = 0, n - 1
              sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
            end do
          end do
        end if
      case (16)
        a = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            b = e
            c = e
            sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end if
        end do
      case (17)
        b = x
        d = y
        do e = 0, n - 1
          adj_value = adj(d * n + e)
          if (adj_value > 0.0_real64) then
            a = e
            c = e
            sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end if
        end do
      case (18)
        d = x
        e = y
        adj_value = adj(d * n + e)
        if (adj_value > 0.0_real64) then
          do a = 0, n - 1
            b = a
            c = a
            sum_value = sum_value + tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) * adj_value
          end do
        end if
      end select

      value(tid) = sum_value
    end do
    !$omp end target teams distribute parallel do
  end subroutine contraction_real64

  real(real32) function tensor_index_real32(tensor, a, b, c, f, a_stride, b_stride, c_stride) result(value)
    real(real32), intent(in) :: tensor(0:)
    integer, intent(in) :: a, b, c, f
    integer(int64), intent(in) :: a_stride, b_stride, c_stride

    value = tensor(int(a, int64) * a_stride + int(b, int64) * b_stride + int(c, int64) * c_stride + f)
  end function tensor_index_real32

  real(real64) function tensor_index_real64(tensor, a, b, c, f, a_stride, b_stride, c_stride) result(value)
    real(real64), intent(in) :: tensor(0:)
    integer, intent(in) :: a, b, c, f
    integer(int64), intent(in) :: a_stride, b_stride, c_stride

    value = tensor(int(a, int64) * a_stride + int(b, int64) * b_stride + int(c, int64) * c_stride + f)
  end function tensor_index_real64

  subroutine summarize_real32(value, output_size, checksum, min_value, max_value)
    real(real32), intent(in) :: value(0:)
    integer(int64), intent(in) :: output_size
    real(real64), intent(out) :: checksum, min_value, max_value
    integer(int64) :: i

    checksum = 0.0_real64
    min_value = real(value(0), real64)
    max_value = real(value(0), real64)
    do i = 0_int64, output_size - 1_int64
      checksum = checksum + real(value(i), real64)
      min_value = min(min_value, real(value(i), real64))
      max_value = max(max_value, real(value(i), real64))
    end do
  end subroutine summarize_real32

  subroutine summarize_real64(value, output_size, checksum, min_value, max_value)
    real(real64), intent(in) :: value(0:)
    integer(int64), intent(in) :: output_size
    real(real64), intent(out) :: checksum, min_value, max_value
    integer(int64) :: i

    checksum = 0.0_real64
    min_value = value(0)
    max_value = value(0)
    do i = 0_int64, output_size - 1_int64
      checksum = checksum + value(i)
      min_value = min(min_value, value(i))
      max_value = max(max_value, value(i))
    end do
  end subroutine summarize_real64

end program main
