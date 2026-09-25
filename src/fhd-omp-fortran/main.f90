! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  real(real32), parameter :: pi = 3.14159265358979323846_real32
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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: samples, voxels, verify
  real(real32), allocatable :: h_rmu(:), h_imu(:), h_kx(:), h_ky(:), h_kz(:)
  real(real32), allocatable :: h_rfhd(:), h_ifhd(:), h_x(:), h_y(:), h_z(:)
  real(real32), allocatable :: rfhd(:), ifhd(:)
  real(real64) :: start_time, end_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <#samples> <#voxels> <verify>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) samples
  read(arg2, *) voxels
  read(arg3, *) verify

  allocate(h_rmu(voxels), h_imu(voxels), h_kx(voxels), h_ky(voxels), h_kz(voxels))
  allocate(h_rfhd(samples), h_ifhd(samples), h_x(samples), h_y(samples), h_z(samples))
  allocate(rfhd(samples), ifhd(samples))

  call initialize_inputs(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, &
    h_rfhd, h_ifhd, h_x, h_y, h_z, rfhd, ifhd)

  write(*,'(A)') 'Run FHd on a device'

  !$omp target data map(to: h_rmu(1:voxels), h_imu(1:voxels), h_kx(1:voxels), &
  !$omp& h_ky(1:voxels), h_kz(1:voxels), h_x(1:samples), h_y(1:samples), h_z(1:samples)) &
  !$omp& map(tofrom: rfhd(1:samples), ifhd(1:samples))
  start_time = omp_get_wtime()
  call fhd_kernel(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, h_x, h_y, h_z, rfhd, ifhd)
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Device execution time ', end_time - start_time, ' (s)'
  !$omp end target data

  if (verify /= 0) then
    call verify_host(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, h_x, h_y, h_z, &
      h_rfhd, h_ifhd, rfhd, ifhd)
  end if

  deallocate(h_rmu, h_imu, h_kx, h_ky, h_kz, h_rfhd, h_ifhd, h_x, h_y, h_z, rfhd, ifhd)

contains

  subroutine initialize_inputs(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, &
      h_rfhd, h_ifhd, h_x, h_y, h_z, rfhd, ifhd)
    integer, intent(in) :: samples, voxels
    real(real32), intent(out) :: h_rmu(:), h_imu(:), h_kx(:), h_ky(:), h_kz(:)
    real(real32), intent(out) :: h_rfhd(:), h_ifhd(:), h_x(:), h_y(:), h_z(:)
    real(real32), intent(out) :: rfhd(:), ifhd(:)
    integer :: idx

    ! Preserve the C++ original's srand(2) / rand()%2 initialization sequence.
    call c_srand(2_c_int)

    do idx = 1, samples
      rfhd(idx) = real(idx - 1, real32) / real(samples, real32)
      ifhd(idx) = rfhd(idx)
      h_rfhd(idx) = rfhd(idx)
      h_ifhd(idx) = ifhd(idx)
      h_x(idx) = 0.3_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
      h_y(idx) = 0.2_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
      h_z(idx) = 0.1_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
    end do

    do idx = 1, voxels
      h_rmu(idx) = real(idx - 1, real32) / real(voxels, real32)
      h_imu(idx) = h_rmu(idx)
      h_kx(idx) = 0.1_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
      h_ky(idx) = 0.2_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
      h_kz(idx) = 0.3_real32 + merge(0.1_real32, -0.1_real32, rand_bit())
    end do
  end subroutine initialize_inputs

  logical function rand_bit()
    rand_bit = mod(c_rand(), 2_c_int) /= 0_c_int
  end function rand_bit

  subroutine fhd_kernel(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, h_x, h_y, h_z, rfhd, ifhd)
    integer, intent(in) :: samples, voxels
    real(real32), intent(in) :: h_rmu(:), h_imu(:), h_kx(:), h_ky(:), h_kz(:), h_x(:), h_y(:), h_z(:)
    real(real32), intent(inout) :: rfhd(:), ifhd(:)
    integer :: n, m
    real(real32) :: r, im, xn, yn, zn, e, c, s

    !$omp target teams distribute parallel do private(n, m, r, im, xn, yn, zn, e, c, s)
    do n = 1, samples
      r = rfhd(n)
      im = ifhd(n)
      xn = h_x(n)
      yn = h_y(n)
      zn = h_z(n)
      !$omp parallel do simd reduction(+:r, im) private(e, c, s)
      do m = 1, voxels
        e = 2.0_real32 * pi * (h_kx(m) * xn + h_ky(m) * yn + h_kz(m) * zn)
        c = cos(e)
        s = sin(e)
        r = r + h_rmu(m) * c - h_imu(m) * s
        im = im + h_imu(m) * c + h_rmu(m) * s
      end do
      !$omp end parallel do simd
      rfhd(n) = r
      ifhd(n) = im
    end do
    !$omp end target teams distribute parallel do
  end subroutine fhd_kernel

  subroutine verify_host(samples, voxels, h_rmu, h_imu, h_kx, h_ky, h_kz, h_x, h_y, h_z, &
      h_rfhd, h_ifhd, rfhd, ifhd)
    integer, intent(in) :: samples, voxels
    real(real32), intent(in) :: h_rmu(:), h_imu(:), h_kx(:), h_ky(:), h_kz(:), h_x(:), h_y(:), h_z(:)
    real(real32), intent(inout) :: h_rfhd(:), h_ifhd(:)
    real(real32), intent(in) :: rfhd(:), ifhd(:)
    integer :: n, m
    real(real32) :: r, im, e, c, s, err

    write(*,'(A)') 'Computing root mean square error between host and device results.'
    write(*,'(A)') 'This will take a while..'

    !$omp parallel do private(n, m, r, im, e, c, s)
    do n = 1, samples
      r = h_rfhd(n)
      im = h_ifhd(n)
      !$omp parallel do simd reduction(+:r, im) private(e, c, s)
      do m = 1, voxels
        e = 2.0_real32 * pi * (h_kx(m) * h_x(n) + h_ky(m) * h_y(n) + h_kz(m) * h_z(n))
        c = cos(e)
        s = sin(e)
        r = r + h_rmu(m) * c - h_imu(m) * s
        im = im + h_imu(m) * c + h_rmu(m) * s
      end do
      !$omp end parallel do simd
      h_rfhd(n) = r
      h_ifhd(n) = im
    end do
    !$omp end parallel do

    err = 0.0_real32
    do n = 1, samples
      err = err + (h_rfhd(n) - rfhd(n)) ** 2 + (h_ifhd(n) - ifhd(n)) ** 2
    end do
    write(*,'(A,F0.6)') 'RMSE = ', sqrt(err / real(2 * samples, real32))
  end subroutine verify_host

end program main
