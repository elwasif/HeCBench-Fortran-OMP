! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: natom = 5877
  integer, parameter :: ngrid = 134918
  integer, parameter :: ngadj = ngrid + (512 - iand(ngrid, 511))
  real(real32), parameter :: c_rand_max = 2147483647.0_real32
  real(real32), parameter :: pre1 = 4.46184985145e19_real32
  real(real32), parameter :: xkappa = 0.0735516324639_real32

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name='rand')
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer :: itmax, wgsize
  real(real32), allocatable :: ax(:), ay(:), az(:), gx(:), gy(:), gz(:)
  real(real32), allocatable :: charge(:), atom_size(:), val_cpu(:), val_gpu(:)
  real(real64) :: timer_start, cpu_start, gpu_start, cpu_elapsed, gpu_elapsed

  itmax = 100
  wgsize = 256
  call getargs(itmax, wgsize)

  allocate(ax(natom), ay(natom), az(natom), charge(natom), atom_size(natom))
  allocate(gx(ngadj), gy(ngadj), gz(ngadj), val_cpu(ngadj), val_gpu(ngadj))

  call gendata(ax, ay, az, gx, gy, gz, charge, atom_size)
  val_cpu = 0.0_real32
  val_gpu = 0.0_real32

  cpu_start = omp_get_wtime()
  call run_cpu_kernel(itmax, ax, ay, az, gx, gy, gz, charge, atom_size, val_cpu)
  cpu_elapsed = omp_get_wtime() - cpu_start
  if (cpu_elapsed < 1.0_real64) then
    write(*,'(A,F0.12,A,I0,A)') 'CPU Time: 0', cpu_elapsed, ' (Number of tests = ', itmax, ')'
  else
    write(*,'(A,F0.12,A,I0,A)') 'CPU Time: ', cpu_elapsed, ' (Number of tests = ', itmax, ')'
  end if
  write(*,*)

  gpu_start = omp_get_wtime()
  call run_gpu_kernel(wgsize, itmax, ax, ay, az, gx, gy, gz, charge, atom_size, val_gpu)
  gpu_elapsed = omp_get_wtime() - gpu_start
  if (gpu_elapsed < 1.0_real64) then
    write(*,'(A,F0.12,A,I0,A)') 'GPU Time: 0', gpu_elapsed, ' (Number of tests = ', itmax, ')'
  else
    write(*,'(A,F0.12,A,I0,A)') 'GPU Time: ', gpu_elapsed, ' (Number of tests = ', itmax, ')'
  end if
  write(*,*)

  call compare_results(val_cpu, val_gpu)

  deallocate(ax, ay, az, charge, atom_size, gx, gy, gz, val_cpu, val_gpu)

contains

  subroutine getargs(itmax, wgsize)
    integer, intent(inout) :: itmax, wgsize
    integer :: i, argc
    character(len=128) :: arg, next_arg

    argc = command_argument_count()
    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      if (trim(arg) == '-itmax' .and. i + 1 <= argc) then
        call get_command_argument(i + 1, next_arg)
        read(next_arg, *) itmax
        i = i + 1
      else if (trim(arg) == '-wgsize' .and. i + 1 <= argc) then
        call get_command_argument(i + 1, next_arg)
        read(next_arg, *) wgsize
        i = i + 1
      end if
      i = i + 1
    end do

    write(*,'(A)') 'Run parameters:'
    write(*,'(A,I0)') '  kernel loop count: ', itmax
    write(*,'(A,I0)') '     workgroup size: ', wgsize
  end subroutine getargs

  subroutine gendata(ax, ay, az, gx, gy, gz, charge, atom_size)
    real(real32), intent(out) :: ax(:), ay(:), az(:), gx(:), gy(:), gz(:), charge(:), atom_size(:)
    integer :: i

    write(*,'(A)') 'Generating Data.. '
    call c_srand(1_c_int)
    do i = 1, natom
      ax(i) = next_rand()
      ay(i) = next_rand()
      az(i) = next_rand()
      charge(i) = next_rand()
      atom_size(i) = real(natom, real32)
    end do
    gx = 0.0_real32
    gy = 0.0_real32
    gz = 0.0_real32
    do i = 1, ngrid
      gx(i) = next_rand()
      gy(i) = next_rand()
      gz(i) = next_rand()
    end do
    write(*,'(A)') 'Done generating inputs.'
    write(*,*)
  end subroutine gendata

  real(real32) function next_rand() result(value)
    value = real(c_rand(), real32) / c_rand_max
  end function next_rand

  subroutine run_cpu_kernel(itmax, ax, ay, az, gx, gy, gz, charge, atom_size, val)
    integer, intent(in) :: itmax
    real(real32), intent(in) :: ax(:), ay(:), az(:), gx(:), gy(:), gz(:), charge(:), atom_size(:)
    real(real32), intent(out) :: val(:)
    integer :: n, igrid, iatom
    real(real32) :: sum_value, l_gx, l_gy, l_gz, dist
    real(real64) :: start_time, elapsed

    start_time = omp_get_wtime()
    do n = 1, itmax
      !$omp parallel do private(igrid, sum_value, l_gx, l_gy, l_gz)
      do igrid = 1, ngadj
        sum_value = 0.0_real32
        l_gx = gx(igrid)
        l_gy = gy(igrid)
        l_gz = gz(igrid)
        !$omp parallel do simd reduction(+:sum_value) private(iatom, dist)
        do iatom = 1, natom
          dist = sqrt((l_gx - ax(iatom)) * (l_gx - ax(iatom)) + &
                      (l_gy - ay(iatom)) * (l_gy - ay(iatom)) + &
                      (l_gz - az(iatom)) * (l_gz - az(iatom)))
          sum_value = sum_value + pre1 * (charge(iatom) / dist) * &
              exp(-xkappa * (dist - atom_size(iatom))) / (1.0_real32 + xkappa * atom_size(iatom))
        end do
        !$omp end parallel do simd
        val(igrid) = sum_value
      end do
      !$omp end parallel do
    end do
    elapsed = (omp_get_wtime() - start_time) / real(itmax, real64)
    write(*,'(A,G0.12)') 'Average kernel execution time: ', elapsed
  end subroutine run_cpu_kernel

  subroutine run_gpu_kernel(wgsize, itmax, ax, ay, az, gx, gy, gz, charge, atom_size, val)
    integer, intent(in) :: wgsize, itmax
    real(real32), intent(in) :: ax(:), ay(:), az(:), gx(:), gy(:), gz(:), charge(:), atom_size(:)
    real(real32), intent(out) :: val(:)
    integer :: n, igrid, iatom
    real(real32) :: sum_value, l_gx, l_gy, l_gz, dist
    real(real64) :: start_time, elapsed

    !$omp target data map(to: ax(1:natom), ay(1:natom), az(1:natom), charge(1:natom), atom_size(1:natom), &
    !$omp& gx(1:ngadj), gy(1:ngadj), gz(1:ngadj)) map(alloc: val(1:ngadj))
    start_time = omp_get_wtime()
    do n = 1, itmax
      !$omp target teams distribute thread_limit(wgsize) private(igrid, sum_value, l_gx, l_gy, l_gz)
      do igrid = 1, ngrid
        sum_value = 0.0_real32
        l_gx = gx(igrid)
        l_gy = gy(igrid)
        l_gz = gz(igrid)
        !$omp parallel do reduction(+:sum_value) private(iatom, dist)
        do iatom = 1, natom
          dist = sqrt((l_gx - ax(iatom)) * (l_gx - ax(iatom)) + &
                      (l_gy - ay(iatom)) * (l_gy - ay(iatom)) + &
                      (l_gz - az(iatom)) * (l_gz - az(iatom)))
          sum_value = sum_value + pre1 * (charge(iatom) / dist) * &
              exp(-xkappa * (dist - atom_size(iatom))) / (1.0_real32 + xkappa * atom_size(iatom))
        end do
        !$omp end parallel do
        val(igrid) = sum_value
      end do
      !$omp end target teams distribute
    end do
    elapsed = (omp_get_wtime() - start_time) / real(itmax, real64)
    write(*,'(A,G0.12)') 'Average kernel time on the device: ', elapsed
    !$omp target update from(val(1:ngrid))
    !$omp end target data
  end subroutine run_gpu_kernel

  subroutine compare_results(arr, arr2)
    real(real32), intent(in) :: arr(:), arr2(:)
    integer :: i
    logical :: ok
    ok = .true.
    do i = 1, ngrid
      if (abs(arr(i) - arr2(i)) > 1.0e-3_real32) then
        ok = .false.
        exit
      end if
    end do
    write(*,'(A)') merge('PASS', 'FAIL', ok)
  end subroutine compare_results

end program main
