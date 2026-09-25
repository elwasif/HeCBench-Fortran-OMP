! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use omp_lib
  implicit none

  integer, parameter :: block_sizes(6) = [32, 64, 128, 256, 512, 1024]
  real(real32), parameter :: eps = 1.0e-5_real32

  interface
    subroutine hecbench_srand(seed) bind(C, name='hecbench_srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine hecbench_srand

    real(c_float) function hecbench_rand_float() bind(C, name='hecbench_rand_float')
      import :: c_float
    end function hecbench_rand_float
  end interface

  character(len=256) :: arg0, arg
  integer :: bsz, tsz, csz, repeat, i, block_size
  integer(int64) :: n_btc, n_bt, memory_ops
  real(real32), allocatable :: out(:), d_out(:), mean(:), d_mean(:), rstd(:), d_rstd(:)
  real(real32), allocatable :: inp(:), weight(:), bias(:)
  real(real64) :: start_time, elapsed_ms, bandwidth
  character(len=32) :: elapsed_text, bandwidth_text

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <batch size> <sequence length> <channel length> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) bsz
  call get_command_argument(2, arg); read(arg, *) tsz
  call get_command_argument(3, arg); read(arg, *) csz
  call get_command_argument(4, arg); read(arg, *) repeat
  if (bsz <= 0 .or. tsz <= 0 .or. csz <= 0 .or. repeat <= 0) stop 1

  n_btc = int(bsz, int64) * int(tsz, int64) * int(csz, int64)
  n_bt = int(bsz, int64) * int(tsz, int64)
  allocate(out(n_btc), d_out(n_btc), mean(n_bt), d_mean(n_bt), rstd(n_bt), d_rstd(n_bt))
  allocate(inp(n_btc), weight(csz), bias(csz))

  call hecbench_srand(0_c_int)
  call fill_random(inp)
  call fill_random(weight)
  call fill_random(bias)

  !$omp target data map(to: inp(1:n_btc), weight(1:csz), bias(1:csz)) &
  !$omp& map(tofrom: d_out(1:n_btc), d_mean(1:n_bt), d_rstd(1:n_bt))
  call layernorm_forward_cpu(out, mean, rstd, inp, weight, bias, bsz, tsz, csz)

  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    write(*,'(A,I0,A)') 'Checking block size ', block_size, '.'

    call layernorm_forward_kernel(d_out, d_mean, d_rstd, inp, weight, bias, bsz, tsz, csz, block_size)
    !$omp target update from(d_out(1:n_btc), d_mean(1:n_bt), d_rstd(1:n_bt))

    call validate_result(d_out, out, 'out', n_btc, 1.0e-5_real32)
    call validate_result(d_mean, mean, 'mean', n_bt, 1.0e-5_real32)
    call validate_result(d_rstd, rstd, 'rstd', n_bt, 1.0e-5_real32)
  end do

  write(*,'(A)') 'All results match. Starting benchmarks.'
  write(*,*)

  do i = 1, size(block_sizes)
    block_size = block_sizes(i)
    start_time = omp_get_wtime()
    call run_repeated(d_out, d_mean, d_rstd, inp, weight, bias, bsz, tsz, csz, block_size, repeat)
    elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
    memory_ops = 2_int64 * n_btc * 4_int64
    bandwidth = real(memory_ops, real64) / elapsed_ms / 1.0e6_real64
    write(elapsed_text, '(F16.4)') elapsed_ms
    write(bandwidth_text, '(F16.2)') bandwidth
    write(*,'(A,I4,A,A,A,A,A)') 'block_size ', block_size, ' | time ', trim(adjustl(elapsed_text)), &
      ' ms | bandwidth ', trim(adjustl(bandwidth_text)), ' GB/s'
  end do
  !$omp end target data

  deallocate(out, d_out, mean, d_mean, rstd, d_rstd, inp, weight, bias)

contains

  subroutine fill_random(values)
    real(real32), intent(out) :: values(:)
    integer :: i

    do i = 1, size(values)
      values(i) = real(hecbench_rand_float(), real32)
    end do
  end subroutine fill_random

  subroutine layernorm_forward_kernel(out, mean, rstd, inp, weight, bias, bsz, tsz, csz, block_size)
    real(real32), intent(out) :: out(:), mean(:), rstd(:)
    real(real32), intent(in) :: inp(:), weight(:), bias(:)
    integer, intent(in) :: bsz, tsz, csz, block_size
    integer :: b, t, i, base, row
    real(real32) :: m, v, s, xshift, nval

    !$omp target teams distribute collapse(2) num_teams(bsz * tsz) private(i, base, row, m, v, s, xshift, nval)
    do b = 0, bsz - 1
      do t = 0, tsz - 1
        base = b * tsz * csz + t * csz
        row = b * tsz + t + 1
        m = 0.0_real32
        !$omp parallel do reduction(+:m) num_threads(block_size)
        do i = 1, csz
          m = m + inp(base + i)
        end do
        !$omp end parallel do
        m = m / real(csz, real32)

        v = 0.0_real32
        !$omp parallel do private(xshift) reduction(+:v) num_threads(block_size)
        do i = 1, csz
          xshift = inp(base + i) - m
          v = v + xshift * xshift
        end do
        !$omp end parallel do
        v = v / real(csz, real32)
        s = 1.0_real32 / sqrt(v + eps)

        !$omp parallel do private(nval) num_threads(block_size)
        do i = 1, csz
          nval = s * (inp(base + i) - m)
          out(base + i) = nval * weight(i) + bias(i)
        end do
        !$omp end parallel do
        mean(row) = m
        rstd(row) = s
      end do
    end do
    !$omp end target teams distribute
  end subroutine layernorm_forward_kernel

  subroutine run_repeated(out, mean, rstd, inp, weight, bias, bsz, tsz, csz, block_size, repeat)
    real(real32), intent(out) :: out(:), mean(:), rstd(:)
    real(real32), intent(in) :: inp(:), weight(:), bias(:)
    integer, intent(in) :: bsz, tsz, csz, block_size, repeat
    integer :: iter

    do iter = 1, repeat
      call layernorm_forward_kernel(out, mean, rstd, inp, weight, bias, bsz, tsz, csz, block_size)
    end do
  end subroutine run_repeated

  subroutine layernorm_forward_cpu(out, mean, rstd, inp, weight, bias, bsz, tsz, csz)
    real(real32), intent(out) :: out(:), mean(:), rstd(:)
    real(real32), intent(in) :: inp(:), weight(:), bias(:)
    integer, intent(in) :: bsz, tsz, csz
    integer :: b, t, i, base, row
    real(real32) :: m, v, s, xshift, nval

    do b = 0, bsz - 1
      do t = 0, tsz - 1
        base = b * tsz * csz + t * csz
        row = b * tsz + t + 1
        m = 0.0_real32
        do i = 1, csz
          m = m + inp(base + i)
        end do
        m = m / real(csz, real32)

        v = 0.0_real32
        do i = 1, csz
          xshift = inp(base + i) - m
          v = v + xshift * xshift
        end do
        v = v / real(csz, real32)
        s = 1.0_real32 / sqrt(v + eps)

        do i = 1, csz
          nval = s * (inp(base + i) - m)
          out(base + i) = nval * weight(i) + bias(i)
        end do
        mean(row) = m
        rstd(row) = s
      end do
    end do
  end subroutine layernorm_forward_cpu

  subroutine validate_result(device_result, cpu_reference, name, num_elements, tolerance)
    real(real32), intent(in) :: device_result(:), cpu_reference(:), tolerance
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: num_elements
    integer(int64) :: i
    integer :: nfaults
    real(real32) :: fp_epsilon, t_eff

    nfaults = 0
    fp_epsilon = epsilon(1.0_real32)
    do i = 1, num_elements
      if (.not. ieee_is_finite(cpu_reference(i))) cycle
      t_eff = tolerance + abs(cpu_reference(i)) * fp_epsilon
      if (abs(cpu_reference(i) - device_result(i)) > t_eff) then
        write(*,'(A,A,A,I0,A,F0.6,A,F0.6)') 'Mismatch of ', trim(name), ' at ', i - 1, &
          ': CPU_ref: ', cpu_reference(i), ' vs GPU: ', device_result(i)
        nfaults = nfaults + 1
        if (nfaults >= 10) stop 1
      end if
    end do
    if (nfaults > 0) stop 1
  end subroutine validate_result

  logical function ieee_is_finite(value)
    real(real32), intent(in) :: value

    ieee_is_finite = value == value .and. abs(value) <= huge(value)
  end function ieee_is_finite

end program main
