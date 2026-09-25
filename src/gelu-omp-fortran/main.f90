! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int16, int32, int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: batch_size, seq_len, hidden_dim, repeat, block_size, i
  integer(int64) :: src_size
  integer(int16), allocatable :: output(:), output_ref(:), bias(:)
  real(real64) :: start_time, end_time, elapsed_ms
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <batch> <sequence length> <hidden dimension> <repeat>'
    write(*,'(A)') 'The hidden dimension is a multiple of two'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) batch_size
  call get_command_argument(2, arg)
  read(arg, *) seq_len
  call get_command_argument(3, arg)
  read(arg, *) hidden_dim
  call get_command_argument(4, arg)
  read(arg, *) repeat
  if (batch_size <= 0 .or. seq_len <= 0 .or. hidden_dim <= 0 .or. repeat <= 0) stop 1
  if (mod(hidden_dim, 2) /= 0) stop 1

  src_size = int(batch_size, int64) * int(seq_len, int64) * int(hidden_dim, int64)
  allocate(output(0:src_size - 1), output_ref(0:src_size - 1), bias(0:hidden_dim - 1))
  call initialize_inputs(output, output_ref, bias, src_size, hidden_dim)

  if (hidden_dim >= 4096) then
    block_size = 512
  else if (hidden_dim >= 2048) then
    block_size = 256
  else
    block_size = 128
  end if

  call gelu_bias_loop_cpu(output_ref, bias, batch_size, hidden_dim, seq_len)

  !$omp target data map(to: bias(0:hidden_dim - 1), output(0:src_size - 1))
  call gelu_bias_loop(output, bias, batch_size, hidden_dim, seq_len, block_size)
  !$omp target update from(output(0:src_size - 1))

  ok = .true.
  do i = 0, int(src_size - 1_int64)
    if (abs(half_to_real(output_ref(i)) - half_to_real(output(i))) > 1.0e-3_real32) then
      ok = .false.
      exit
    end if
  end do
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  start_time = omp_get_wtime()
  do i = 1, repeat
    call gelu_bias_loop(output, bias, batch_size, hidden_dim, seq_len, block_size)
  end do
  end_time = omp_get_wtime()
  elapsed_ms = (end_time - start_time) * 1.0e3_real64 / real(repeat, real64)
  call print_time('Average execution time of vectorized kernel ', elapsed_ms)

  start_time = omp_get_wtime()
  do i = 1, repeat
    call gelu_bias_loop_base(output, bias, batch_size, hidden_dim, seq_len, block_size)
  end do
  end_time = omp_get_wtime()
  elapsed_ms = (end_time - start_time) * 1.0e3_real64 / real(repeat, real64)
  call print_time('Average execution time of baseline kernel ', elapsed_ms)
  !$omp end target data

  deallocate(output, output_ref, bias)

contains

  subroutine initialize_inputs(output, output_ref, bias, src_size, hidden_dim)
    integer(int16), intent(out) :: output(0:), output_ref(0:), bias(0:)
    integer(int64), intent(in) :: src_size
    integer, intent(in) :: hidden_dim
    integer(int64) :: i
    integer :: j

    call c_srand(123_c_int)
    do i = 0_int64, src_size - 1_int64
      output(i) = real_to_half(real(c_rand(), real32) / real(huge(0_c_int), real32))
      output_ref(i) = output(i)
    end do
    do j = 0, hidden_dim - 1
      bias(j) = real_to_half(real(-6_c_int + mod(c_rand(), 12_c_int), real32))
    end do
  end subroutine initialize_inputs

  subroutine gelu_bias_loop(src, bias, batch_size, width, height, block_size)
    integer(int16), intent(inout) :: src(0:)
    integer(int16), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height, block_size
    integer :: batch, x, y
    integer(int64) :: index, left, right
    real(real32) :: tx, ty

    !$omp target teams distribute collapse(2) num_teams(batch_size * height) private(batch, x, y, index, left, right, tx, ty)
    do batch = 0, batch_size - 1
      do x = 0, height - 1
        index = (int(batch, int64) * width * height + int(x, int64) * width) / 2_int64
        !$omp parallel do num_threads(block_size) private(y, left, right, tx, ty)
        do y = 0, width / 2 - 1
          left = 2_int64 * (index + int(y, int64))
          right = left + 1_int64
          tx = half_to_real(src(left)) + half_to_real(bias(2 * y))
          ty = half_to_real(src(right)) + half_to_real(bias(2 * y + 1))
          src(left) = real_to_half(gelu_value(tx))
          src(right) = real_to_half(gelu_value(ty))
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine gelu_bias_loop

  subroutine gelu_bias_loop_base(src, bias, batch_size, width, height, block_size)
    integer(int16), intent(inout) :: src(0:)
    integer(int16), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height, block_size
    integer :: batch, x, y
    integer(int64) :: base, idx
    real(real32) :: t

    !$omp target teams distribute collapse(2) num_teams(batch_size * height) private(batch, x, y, base, idx, t)
    do batch = 0, batch_size - 1
      do x = 0, height - 1
        base = int(batch, int64) * width * height + int(x, int64) * width
        !$omp parallel do num_threads(block_size) private(y, idx, t)
        do y = 0, width - 1
          idx = base + int(y, int64)
          t = half_to_real(src(idx)) + half_to_real(bias(y))
          src(idx) = real_to_half(gelu_value(t))
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine gelu_bias_loop_base

  subroutine gelu_bias_loop_cpu(src, bias, batch_size, width, height)
    integer(int16), intent(inout) :: src(0:)
    integer(int16), intent(in) :: bias(0:)
    integer, intent(in) :: batch_size, width, height
    integer :: batch, x, y
    integer(int64) :: idx
    real(real32) :: t

    do batch = 0, batch_size - 1
      do x = 0, height - 1
        do y = 0, width - 1
          idx = int(batch, int64) * width * height + int(x, int64) * width + int(y, int64)
          t = half_to_real(src(idx)) + half_to_real(bias(y))
          src(idx) = real_to_half(gelu_value(t))
        end do
      end do
    end do
  end subroutine gelu_bias_loop_cpu

  pure real(real32) function gelu_value(t)
    real(real32), intent(in) :: t

    gelu_value = 0.5_real32 * t * (1.0_real32 + tanh(0.79788456_real32 * (t + 0.044715_real32 * t * t * t)))
  end function gelu_value

  subroutine print_time(label, elapsed_ms)
    character(len=*), intent(in) :: label
    real(real64), intent(in) :: elapsed_ms

    if (elapsed_ms >= 0.0_real64 .and. elapsed_ms < 1.0_real64) then
      write(*,'(A,A,F0.6,A)') label, '0', elapsed_ms, ' (ms)'
    else
      write(*,'(A,F0.6,A)') label, elapsed_ms, ' (ms)'
    end if
  end subroutine print_time

  pure real(real32) function half_to_real(h)
    integer(int16), intent(in) :: h
    integer(int32) :: bits, sign, exponent, fraction
    real(real32) :: mantissa

    bits = iand(int(h, int32), int(z'0000ffff', int32))
    sign = iand(bits, int(z'00008000', int32))
    exponent = iand(shiftr(bits, 10), 31_int32)
    fraction = iand(bits, int(z'000003ff', int32))

    if (exponent == 0) then
      if (fraction == 0) then
        half_to_real = 0.0_real32
      else
        half_to_real = scale(real(fraction, real32) / 1024.0_real32, -14)
      end if
    else if (exponent == 31) then
      half_to_real = huge(half_to_real)
    else
      mantissa = 1.0_real32 + real(fraction, real32) / 1024.0_real32
      half_to_real = scale(mantissa, exponent - 15)
    end if
    if (sign /= 0) half_to_real = -half_to_real
  end function half_to_real

  pure integer(int16) function real_to_half(value)
    real(real32), intent(in) :: value
    integer(int32) :: sign, exponent, fraction, bits
    real(real32) :: ax, mantissa

    if (value < 0.0_real32) then
      sign = int(z'00008000', int32)
      ax = -value
    else
      sign = 0_int32
      ax = value
    end if

    if (ax == 0.0_real32) then
      bits = sign
    else if (ax >= 65504.0_real32) then
      bits = sign + int(z'00007bff', int32)
    else if (ax < scale(1.0_real32, -14)) then
      fraction = nint(ax * scale(1.0_real32, 24))
      if (fraction > 1023) fraction = 1023
      bits = sign + fraction
    else
      exponent = floor(log(ax) / log(2.0_real32))
      mantissa = ax / scale(1.0_real32, exponent) - 1.0_real32
      fraction = nint(mantissa * 1024.0_real32)
      if (fraction == 1024) then
        exponent = exponent + 1
        fraction = 0
      end if
      bits = sign + shiftl(exponent + 15, 10) + fraction
    end if

    real_to_half = int(iand(bits, int(z'0000ffff', int32)), int16)
  end function real_to_half

end program main
