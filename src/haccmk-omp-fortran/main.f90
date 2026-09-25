! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer :: repeat, n1, n2, i, error
  real(real32) :: fsrrmax2, mp_rsm2, fcoeff, dx1, dy1, dz1
  real(real32) :: dx2, dy2, dz2, eps
  real(real32), allocatable :: xx(:), yy(:), zz(:), mass(:)
  real(real32), allocatable :: vx2(:), vy2(:), vz2(:), vx2_hw(:), vy2_hw(:), vz2_hw(:)

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <repeat>'
    stop 1
  end if
  repeat = read_arg(1)

  n1 = 784
  n2 = 15000
  print '(A,I0)', 'Outer loop count is set ', n1
  print '(A,I0)', 'Inner loop count is set ', n2

  allocate(xx(n2), yy(n2), zz(n2), mass(n2))
  allocate(vx2(n2), vy2(n2), vz2(n2), vx2_hw(n2), vy2_hw(n2), vz2_hw(n2))

  fcoeff = 0.23_real32
  fsrrmax2 = 0.5_real32
  mp_rsm2 = 0.03_real32
  dx1 = 1.0_real32 / real(n2, real32)
  dy1 = 2.0_real32 / real(n2, real32)
  dz1 = 3.0_real32 / real(n2, real32)
  xx(1) = 0.0_real32
  yy(1) = 0.0_real32
  zz(1) = 0.0_real32
  mass(1) = 2.0_real32

  do i = 2, n2
    xx(i) = xx(i - 1) + dx1
    yy(i) = yy(i - 1) + dy1
    zz(i) = zz(i - 1) + dz1
    mass(i) = real(i - 1, real32) * 0.01_real32 + xx(i)
  end do

  vx2 = 0.0_real32
  vy2 = 0.0_real32
  vz2 = 0.0_real32
  vx2_hw = 0.0_real32
  vy2_hw = 0.0_real32
  vz2_hw = 0.0_real32

  do i = 1, n1
    call haccmk_gold(n2, xx(i), yy(i), zz(i), fsrrmax2, mp_rsm2, xx, yy, zz, mass, dx2, dy2, dz2)
    vx2(i) = vx2(i) + dx2 * fcoeff
    vy2(i) = vy2(i) + dy2 * fcoeff
    vz2(i) = vz2(i) + dz2 * fcoeff
  end do

  call haccmk(repeat, n1, n2, fsrrmax2, mp_rsm2, fcoeff, xx, yy, zz, mass, vx2_hw, vy2_hw, vz2_hw)

  error = 0
  eps = 1.0_real32
  do i = 1, n2
    if (abs(vx2(i) - vx2_hw(i)) > eps) then
      print '(A,I0,A,F0.6,1X,F0.6)', 'error at vx2[', i - 1, '] ', vx2(i), vx2_hw(i)
      error = 1
      exit
    end if
    if (abs(vy2(i) - vy2_hw(i)) > eps) then
      print '(A,I0,A,F0.6,1X,F0.6)', 'error at vy2[', i - 1, ']: ', vy2(i), vy2_hw(i)
      error = 1
      exit
    end if
    if (abs(vz2(i) - vz2_hw(i)) > eps) then
      print '(A,I0,A,F0.6,1X,F0.6)', 'error at vz2[', i - 1, ']: ', vz2(i), vz2_hw(i)
      error = 1
      exit
    end if
  end do

  if (error /= 0) then
    print '(A)', 'FAIL'
  else
    print '(A)', 'PASS'
  end if

  deallocate(xx, yy, zz, mass, vx2, vy2, vz2, vx2_hw, vy2_hw, vz2_hw)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine haccmk(repeat, n, ilp, fsrrmax, mp_rsm, fcoeff, xx, yy, zz, mass, vx2, vy2, vz2)
    integer, intent(in) :: repeat, n, ilp
    real(real32), intent(in) :: fsrrmax, mp_rsm, fcoeff
    real(real32), intent(in) :: xx(:), yy(:), zz(:), mass(:)
    real(real32), intent(inout) :: vx2(:), vy2(:), vz2(:)
    integer :: rep, idx, j
    real(real32) :: ma0, ma1, ma2, ma3, ma4, ma5
    real(real32) :: dxc, dyc, dzc, m, r2, f, xi, yi, zi, xxi, yyi, zzi
    real(real64) :: start_time, end_time, total_time

    ma0 = 0.269327_real32
    ma1 = -0.0750978_real32
    ma2 = 0.0114808_real32
    ma3 = -0.00109313_real32
    ma4 = 0.0000605491_real32
    ma5 = -0.00000147177_real32
    total_time = 0.0_real64

    !$omp target data map(to: xx(1:ilp), yy(1:ilp), zz(1:ilp), mass(1:ilp)) &
    !$omp& map(from: vx2(1:n), vy2(1:n), vz2(1:n))
    do rep = 1, repeat
      !$omp target update to(vx2(1:n))
      !$omp target update to(vy2(1:n))
      !$omp target update to(vz2(1:n))
      start_time = omp_get_wtime()
      !$omp target teams distribute parallel do private(dxc, dyc, dzc, m, r2, f, xi, yi, zi, xxi, yyi, zzi, j)
      do idx = 1, n
        xi = 0.0_real32
        yi = 0.0_real32
        zi = 0.0_real32
        xxi = xx(idx)
        yyi = yy(idx)
        zzi = zz(idx)
        do j = 1, ilp
          dxc = xx(j) - xxi
          dyc = yy(j) - yyi
          dzc = zz(j) - zzi
          r2 = dxc * dxc + dyc * dyc + dzc * dzc
          if (r2 < fsrrmax) then
            m = mass(j)
          else
            m = 0.0_real32
          end if
          f = r2 + mp_rsm
          f = m * (1.0_real32 / (f * sqrt(f)) - &
              (ma0 + r2 * (ma1 + r2 * (ma2 + r2 * (ma3 + r2 * (ma4 + r2 * ma5))))))
          xi = xi + f * dxc
          yi = yi + f * dyc
          zi = zi + f * dzc
        end do
        vx2(idx) = vx2(idx) + xi * fcoeff
        vy2(idx) = vy2(idx) + yi * fcoeff
        vz2(idx) = vz2(idx) + zi * fcoeff
      end do
      !$omp end target teams distribute parallel do
      end_time = omp_get_wtime()
      total_time = total_time + end_time - start_time
    end do
    !$omp end target data

    print '(A,F8.6,A)', 'Average kernel execution time ', total_time / real(repeat, real64), ' (s)'
  end subroutine haccmk

  subroutine haccmk_gold(count1, xxi, yyi, zzi, fsrrmax2, mp_rsm2, xx1, yy1, zz1, mass1, dxi, dyi, dzi)
    integer, intent(in) :: count1
    real(real32), intent(in) :: xxi, yyi, zzi, fsrrmax2, mp_rsm2
    real(real32), intent(in) :: xx1(:), yy1(:), zz1(:), mass1(:)
    real(real32), intent(out) :: dxi, dyi, dzi
    integer :: j
    real(real32) :: ma0, ma1, ma2, ma3, ma4, ma5
    real(real32) :: dxc, dyc, dzc, m, r2, f, xi, yi, zi

    ma0 = 0.269327_real32
    ma1 = -0.0750978_real32
    ma2 = 0.0114808_real32
    ma3 = -0.00109313_real32
    ma4 = 0.0000605491_real32
    ma5 = -0.00000147177_real32
    xi = 0.0_real32
    yi = 0.0_real32
    zi = 0.0_real32

    do j = 1, count1
      dxc = xx1(j) - xxi
      dyc = yy1(j) - yyi
      dzc = zz1(j) - zzi
      r2 = dxc * dxc + dyc * dyc + dzc * dzc
      if (r2 < fsrrmax2) then
        m = mass1(j)
      else
        m = 0.0_real32
      end if
      f = r2 + mp_rsm2
      f = m * (1.0_real32 / (f * sqrt(f)) - &
          (ma0 + r2 * (ma1 + r2 * (ma2 + r2 * (ma3 + r2 * (ma4 + r2 * ma5))))))
      xi = xi + f * dxc
      yi = yi + f * dyc
      zi = zi + f * dzc
    end do
    dxi = xi
    dyi = yi
    dzi = zi
  end subroutine haccmk_gold

end program main
