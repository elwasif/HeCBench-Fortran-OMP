! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: tile_width = 16
  character(len=256) :: arg0, arg
  integer :: n_batch, channels, maps, win, hin, kernel, repeat

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

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 7) then
    write(*,'(A,A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)', advance='no') ' <batch size:N> <input channels:C> <output feature maps:M>'
    write(*,'(A)') ' <input width:Win> <input height:Hin> <kernel size:K> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) n_batch
  call get_command_argument(2, arg); read(arg, *) channels
  call get_command_argument(3, arg); read(arg, *) maps
  call get_command_argument(4, arg); read(arg, *) win
  call get_command_argument(5, arg); read(arg, *) hin
  call get_command_argument(6, arg); read(arg, *) kernel
  call get_command_argument(7, arg); read(arg, *) repeat
  if (n_batch <= 0 .or. channels <= 0 .or. maps <= 0 .or. win <= 0 .or. hin <= 0 .or. &
      kernel <= 0 .or. repeat <= 0 .or. kernel > win .or. kernel > hin) stop 1

  write(*,'(A)') '3D convolution (FP32)'
  write(*,*)
  write(*,'(A)') '========== Warmup start =========='
  call conv3d(n_batch, channels, maps, win, hin, kernel, 1000)
  write(*,*)
  write(*,'(A)') '========== Warmup done =========='
  call conv3d(n_batch, channels, maps, win, hin, kernel, repeat)

contains

  subroutine conv3d(n_batch, channels, maps, win, hin, kernel, repeat)
    integer, intent(in) :: n_batch, channels, maps, win, hin, kernel, repeat
    integer :: hout, wout, x_size, w_size, y_size
    real(real32), allocatable :: x(:), weights(:), y(:), y_ref(:)
    integer :: idx
    real(real64) :: start_time, elapsed
    logical :: ok

    hout = hin - kernel + 1
    wout = win - kernel + 1
    x_size = n_batch * channels * hin * win
    w_size = maps * channels * kernel * kernel
    y_size = n_batch * maps * hout * wout

    allocate(x(x_size), weights(w_size), y(y_size), y_ref(y_size))

    call c_srand(123_c_int)

    do idx = 1, w_size
      weights(idx) = real(mod(c_rand(), 31_c_int), real32)
    end do
    do idx = 1, x_size
      x(idx) = real(mod(c_rand(), 13_c_int), real32)
    end do
    y = -1.0_real32
    y_ref = -1.0_real32

    write(*,'(A,I0,A,I0,A,I0)') 'input dimensions: C=', channels, ' Win=', win, ' Hin=', hin
    write(*,'(A,I0,A,I0,A,I0)') 'output dimensions: M=', maps, ' Wout=', wout, ' Hout=', hout

    !$omp target data map(to: x(1:x_size), weights(1:w_size)) map(tofrom: y(1:y_size))
    start_time = omp_get_wtime()
    do idx = 1, repeat
      call conv3d_kernel(x, weights, y, n_batch, channels, maps, win, hin, kernel, hout, wout)
    end do
    elapsed = omp_get_wtime() - start_time
    !$omp end target data

    write(*,'(A,F0.6,A)') 'Average kernel execution time of conv3d kernel: ', &
      elapsed * 1.0e6_real64 / real(repeat, real64), ' (us)'

    call reference_conv3d(x, weights, y_ref, n_batch, channels, maps, win, hin, kernel, hout, wout)
    ok = all(abs(y - y_ref) <= 1.0e-3_real32)
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    deallocate(x, weights, y, y_ref)
  end subroutine conv3d

  subroutine conv3d_kernel(x, weights, y, n_batch, channels, maps, win, hin, kernel, hout, wout)
    real(real32), intent(in) :: x(:), weights(:)
    real(real32), intent(out) :: y(:)
    integer, intent(in) :: n_batch, channels, maps, win, hin, kernel, hout, wout
    integer :: n, m, h, w, c, p, q
    real(real32) :: total

    !$omp target teams distribute parallel do collapse(4) thread_limit(tile_width * tile_width)
    do n = 0, n_batch - 1
      do m = 0, maps - 1
        do h = 0, hout - 1
          do w = 0, wout - 1
            total = 0.0_real32
            do c = 0, channels - 1
              do p = 0, kernel - 1
                do q = 0, kernel - 1
                  total = total + x(input_index(n, c, h + p, w + q, channels, hin, win)) * &
                    weights(weight_index(m, c, p, q, channels, kernel))
                end do
              end do
            end do
            y(output_index(n, m, h, w, maps, hout, wout)) = total
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine conv3d_kernel

  subroutine reference_conv3d(x, weights, y, n_batch, channels, maps, win, hin, kernel, hout, wout)
    real(real32), intent(in) :: x(:), weights(:)
    real(real32), intent(out) :: y(:)
    integer, intent(in) :: n_batch, channels, maps, win, hin, kernel, hout, wout
    integer :: n, m, h, w, c, p, q
    real(real32) :: total

    do n = 0, n_batch - 1
      do m = 0, maps - 1
        do h = 0, hout - 1
          do w = 0, wout - 1
            total = 0.0_real32
            do c = 0, channels - 1
              do p = 0, kernel - 1
                do q = 0, kernel - 1
                  total = total + x(input_index(n, c, h + p, w + q, channels, hin, win)) * &
                    weights(weight_index(m, c, p, q, channels, kernel))
                end do
              end do
            end do
            y(output_index(n, m, h, w, maps, hout, wout)) = total
          end do
        end do
      end do
    end do
  end subroutine reference_conv3d

  integer function input_index(n, c, h, w, channels, hin, win)
    integer, intent(in) :: n, c, h, w, channels, hin, win
    input_index = n * channels * hin * win + c * hin * win + h * win + w + 1
  end function input_index

  integer function weight_index(m, c, h, w, channels, kernel)
    integer, intent(in) :: m, c, h, w, channels, kernel
    weight_index = m * channels * kernel * kernel + c * kernel * kernel + h * kernel + w + 1
  end function weight_index

  integer function output_index(n, m, h, w, maps, hout, wout)
    integer, intent(in) :: n, m, h, w, maps, hout, wout
    output_index = n * maps * hout * wout + m * hout * wout + h * wout + w + 1
  end function output_index

end program main
