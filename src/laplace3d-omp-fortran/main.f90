! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_x = 32
  integer, parameter :: block_y = 8

  character(len=256) :: arg
  integer :: nx, ny, nz, repeat, verify
  integer(int64) :: grid_size
  integer :: i, j, k, iter, ind
  real(real32), allocatable :: u1(:), u2(:), u3(:), final_gpu(:), final_ref(:)
  real(real64) :: start_time, elapsed, err
  logical :: ok

  if (command_argument_count() /= 5) then
    call print_help()
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) nx
  call get_command_argument(2, arg); read(arg, *) ny
  call get_command_argument(3, arg); read(arg, *) nz
  call get_command_argument(4, arg); read(arg, *) repeat
  call get_command_argument(5, arg); read(arg, *) verify

  if (nx <= 0 .or. mod(nx, 32) /= 0 .or. ny <= 0 .or. nz <= 0 .or. repeat <= 0) stop 1

  write(*,*)
  write(*,'(A,I0,A,I0,A,I0)') 'Grid dimensions: ', nx, ' x ', ny, ' x ', nz
  if (verify /= 0) then
    write(*,'(A)') 'Result verification enabled '
  else
    write(*,'(A)') 'Result verification disabled '
  end if

  grid_size = int(nx, int64) * int(ny, int64) * int(nz, int64)
  allocate(u1(grid_size), u2(grid_size), u3(grid_size), final_gpu(grid_size), final_ref(grid_size))

  do k = 0, nz - 1
    do j = 0, ny - 1
      do i = 0, nx - 1
        ind = i + j * nx + k * nx * ny + 1
        if (i == 0 .or. i == nx - 1 .or. j == 0 .or. j == ny - 1 .or. k == 0 .or. k == nz - 1) then
          u1(ind) = 1.0_real32
        else
          u1(ind) = 0.0_real32
        end if
      end do
    end do
  end do
  u2 = u1

  !$omp target data map(tofrom: u1(1:grid_size), u2(1:grid_size))
  call laplace3d(nx, ny, nz, u1, u2)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    if (mod(iter, 2) == 1) then
      call laplace3d(nx, ny, nz, u1, u2)
    else
      call laplace3d(nx, ny, nz, u2, u1)
    end if
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed / real(repeat, real64), ' (s)'

  if (verify /= 0) then
    if (mod(repeat, 2) == 1) then
      final_gpu = u2
    else
      final_gpu = u1
    end if

    u2 = final_ref_initial(nx, ny, nz, grid_size)
    do iter = 1, repeat
      if (mod(iter, 2) == 1) then
        call reference(nx, ny, nz, u2, u3)
      else
        call reference(nx, ny, nz, u3, u2)
      end if
    end do
    if (mod(repeat, 2) == 1) then
      final_ref = u3
    else
      final_ref = u2
    end if

    err = 0.0_real64
    ok = .true.
    do i = 1, grid_size
      err = err + real((final_gpu(i) - final_ref(i)) * (final_gpu(i) - final_ref(i)), real64)
      if (abs(final_gpu(i) - final_ref(i)) > 1.0e-3_real32) ok = .false.
    end do
    write(*,*)
    write(*,'(A,F0.6,A)') ' RMS error = ', sqrt(err / real(nx, real64) * real(ny, real64) * real(nz, real64)), ' '
    if (ok) then
      write(*,'(A)') ' PASS'
    else
      write(*,'(A)') ' FAIL'
    end if
  end if

  deallocate(u1, u2, u3, final_gpu, final_ref)

contains

  function final_ref_initial(nx, ny, nz, grid_size) result(values)
    integer, intent(in) :: nx, ny, nz
    integer(int64), intent(in) :: grid_size
    real(real32) :: values(grid_size)
    integer :: i, j, k, ind

    do k = 0, nz - 1
      do j = 0, ny - 1
        do i = 0, nx - 1
          ind = i + j * nx + k * nx * ny + 1
          if (i == 0 .or. i == nx - 1 .or. j == 0 .or. j == ny - 1 .or. k == 0 .or. k == nz - 1) then
            values(ind) = 1.0_real32
          else
            values(ind) = 0.0_real32
          end if
        end do
      end do
    end do
  end function final_ref_initial

  subroutine laplace3d(nx, ny, nz, u1, u2)
    integer, intent(in) :: nx, ny, nz
    real(real32), intent(in) :: u1(:)
    real(real32), intent(out) :: u2(:)
    integer :: i, j, k, ind

    !$omp target teams distribute parallel do collapse(3) thread_limit(block_x * block_y) private(ind)
    do k = 0, nz - 1
      do j = 0, ny - 1
        do i = 0, nx - 1
          ind = i + j * nx + k * nx * ny + 1
          if (i == 0 .or. i == nx - 1 .or. j == 0 .or. j == ny - 1 .or. k == 0 .or. k == nz - 1) then
            u2(ind) = u1(ind)
          else
            u2(ind) = (u1(ind - 1) + u1(ind + 1) + u1(ind - nx) + u1(ind + nx) + &
              u1(ind - nx * ny) + u1(ind + nx * ny)) * (1.0_real32 / 6.0_real32)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine laplace3d

  subroutine reference(nx, ny, nz, u1, u2)
    integer, intent(in) :: nx, ny, nz
    real(real32), intent(in) :: u1(:)
    real(real32), intent(out) :: u2(:)
    integer :: i, j, k, ind

    do k = 0, nz - 1
      do j = 0, ny - 1
        do i = 0, nx - 1
          ind = i + j * nx + k * nx * ny + 1
          if (i == 0 .or. i == nx - 1 .or. j == 0 .or. j == ny - 1 .or. k == 0 .or. k == nz - 1) then
            u2(ind) = u1(ind)
          else
            u2(ind) = (u1(ind - 1) + u1(ind + 1) + u1(ind - nx) + u1(ind + nx) + &
              u1(ind - nx * ny) + u1(ind + nx * ny)) * (1.0_real32 / 6.0_real32)
          end if
        end do
      end do
    end do
  end subroutine reference

  subroutine print_help()
    write(*,'(A)') 'Usage:  laplace3d [OPTION]...'
    write(*,'(A)') '6-point stencil 3D Laplace test '
    write(*,*)
    write(*,'(A)') 'Example: run 100 iterations on a 256x128x128 grid'
    write(*,'(A)') './main 256 128 128 100 1'
    write(*,*)
    write(*,'(A)') 'Options:'
    write(*,'(A)') 'Grid width'
    write(*,'(A)') 'Grid height'
    write(*,'(A)') 'Grid depth'
    write(*,'(A)') 'Number of repetitions'
    write(*,'(A)') 'verify the result'
  end subroutine print_help

end program main
