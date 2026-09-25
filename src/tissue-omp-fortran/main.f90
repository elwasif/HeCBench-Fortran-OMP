! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  interface
    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand

    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand
  end interface

  character(len=256) :: arg0, arg
  integer, parameter :: nnt_dev = 32 * 32 * 32
  integer, parameter :: nsp = 2
  integer, parameter :: step = 4
  integer, parameter :: isp = 2
  integer :: dim, repeat, nnt, i, iter
  real(real32), parameter :: c_rand_max = 2147483647.0_real32
  integer, allocatable :: tisspoints(:)
  real(real32), allocatable :: gtt(:), gbartt(:), ct(:), ctprev(:), qt(:), ct_gold(:)
  real(real64) :: start_time, elapsed
  character(len=32) :: time_text
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <dimension of a 3D grid> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) dim
  if (dim > 32) then
    write(*,'(A)') 'Maximum dimension is 32'
    stop 1
  end if
  call get_command_argument(2, arg); read(arg, *) repeat
  if (dim <= 0 .or. repeat <= 0) stop 1

  nnt = dim * dim * dim
  allocate(tisspoints(3 * nnt_dev))
  allocate(gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev))
  allocate(ct(nnt_dev), ctprev(nnt_dev), qt(nnt_dev), ct_gold(nnt_dev))

  call c_srand(1_c_int)
  do i = 1, 3 * nnt_dev
    tisspoints(i) = modulo(c_rand(), nnt_dev / 3)
  end do
  do i = 1, nsp * nnt_dev
    gtt(i) = real(c_rand(), real32) / c_rand_max
    gbartt(i) = real(c_rand(), real32) / c_rand_max
  end do
  do i = 1, nnt_dev
    ct(i) = 0.0_real32
    ct_gold(i) = 0.0_real32
    ctprev(i) = real(c_rand(), real32) / c_rand_max
    qt(i) = real(c_rand(), real32) / c_rand_max
  end do

  !$omp target data map(to: tisspoints(1:3*nnt_dev), gtt(1:nsp*nnt_dev), gbartt(1:nsp*nnt_dev), &
  !$omp& ctprev(1:nnt_dev), qt(1:nnt_dev)) map(tofrom: ct(1:nnt_dev))
  do i = 1, 2
    call tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, step, isp)
  end do

  do i = 1, 2
    call reference(tisspoints, gtt, gbartt, ct_gold, ctprev, qt, nnt, nnt_dev, step, isp)
  end do

  !$omp target update from(ct(1:nnt_dev))
  ok = .true.
  do i = 1, nnt_dev
    if (abs(ct(i) - ct_gold(i)) > 1.0e-1_real32) then
      write(*,'(A,I0,A,F0.6,1X,F0.6)') '@', i - 1, ': ', ct(i), ct_gold(i)
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
  do iter = 1, repeat
    call tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, step, isp)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(time_text,'(F0.6)') elapsed / real(repeat, real64)
  if (time_text(1:1) == '.') time_text = '0' // trim(time_text)
  write(*,'(A,A,A)') 'Average kernel execution time: ', trim(time_text), ' (s)'

  deallocate(tisspoints, gtt, gbartt, ct, ctprev, qt, ct_gold)

contains

  subroutine tissue(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, step, isp)
    integer, intent(in) :: nnt, nnt_dev, step, isp
    integer, intent(in) :: tisspoints(3 * nnt_dev)
    real(real32), intent(in) :: gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev), ctprev(nnt_dev), qt(nnt_dev)
    real(real32), intent(inout) :: ct(nnt_dev)
    integer :: i, itp, itp1, jtp, ix, iy, iz, jx, jy, jz, ixyz, istep
    real(real32) :: p

    !$omp target teams distribute parallel do private(jtp, ix, iy, iz, jx, jy, jz, ixyz, p, itp, itp1, istep) thread_limit(256)
    do i = 0, step * nnt - 1
      itp = i / step
      itp1 = modulo(i, step)
      p = 0.0_real32
      if (itp < nnt) then
        ix = tisspoints(itp + 1)
        iy = tisspoints(itp + nnt + 1)
        iz = tisspoints(itp + 2 * nnt + 1)
        do jtp = itp1, nnt - 1, step
          jx = tisspoints(jtp + 1)
          jy = tisspoints(jtp + nnt + 1)
          jz = tisspoints(jtp + 2 * nnt + 1)
          ixyz = abs(jx - ix) + abs(jy - iy) + abs(jz - iz) + (isp - 1) * nnt_dev
          p = p + gtt(ixyz + 1) * ctprev(jtp + 1) + gbartt(ixyz + 1) * qt(jtp + 1)
        end do
        if (itp1 == 0) ct(itp + 1) = p
      end if
      do istep = 1, step - 1
        if (itp1 == istep .and. itp < nnt) ct(itp + 1) = ct(itp + 1) + p
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine tissue

  subroutine reference(tisspoints, gtt, gbartt, ct, ctprev, qt, nnt, nnt_dev, step, isp)
    integer, intent(in) :: nnt, nnt_dev, step, isp
    integer, intent(in) :: tisspoints(3 * nnt_dev)
    real(real32), intent(in) :: gtt(nsp * nnt_dev), gbartt(nsp * nnt_dev), ctprev(nnt_dev), qt(nnt_dev)
    real(real32), intent(inout) :: ct(nnt_dev)
    integer :: i, itp, itp1, jtp, ix, iy, iz, jx, jy, jz, ixyz, istep
    real(real32) :: p

    do i = 0, step * nnt - 1
      itp = i / step
      itp1 = modulo(i, step)
      p = 0.0_real32
      if (itp < nnt) then
        ix = tisspoints(itp + 1)
        iy = tisspoints(itp + nnt + 1)
        iz = tisspoints(itp + 2 * nnt + 1)
        do jtp = itp1, nnt - 1, step
          jx = tisspoints(jtp + 1)
          jy = tisspoints(jtp + nnt + 1)
          jz = tisspoints(jtp + 2 * nnt + 1)
          ixyz = abs(jx - ix) + abs(jy - iy) + abs(jz - iz) + (isp - 1) * nnt_dev
          p = p + gtt(ixyz + 1) * ctprev(jtp + 1) + gbartt(ixyz + 1) * qt(jtp + 1)
        end do
        if (itp1 == 0) ct(itp + 1) = p
      end if
      do istep = 1, step - 1
        if (itp1 == istep .and. itp < nnt) ct(itp + 1) = ct(itp + 1) + p
      end do
    end do
  end subroutine reference

end program main
