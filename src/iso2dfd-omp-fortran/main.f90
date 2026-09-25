! SPDX-License-Identifier: CC0-1.0
program iso2dfd
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  real(real32), parameter :: dt = 0.002_real32
  real(real32), parameter :: dxy = 20.0_real32
  integer, parameter :: half_length = 1

  real(real32), allocatable :: prev_base(:), next_base(:), next_cpu(:), vel_base(:)
  integer :: nrows, ncols, niterations, nsize
  character(len=64) :: arg
  real(real32) :: dtdivdxy
  real(real64) :: kstart, ktime, cpu_start, cpu_time_ms
  logical :: error
  integer :: k

  if (command_argument_count() < 3) then
    call usage()
  end if

  call get_command_argument(1, arg)
  read(arg, *) nrows
  call get_command_argument(2, arg)
  read(arg, *) ncols
  call get_command_argument(3, arg)
  read(arg, *) niterations

  nsize = nrows * ncols
  allocate(prev_base(0:nsize-1), next_base(0:nsize-1), next_cpu(0:nsize-1), vel_base(0:nsize-1))

  dtdivdxy = (dt * dt) / (dxy * dxy)
  call initialize(prev_base, next_base, vel_base, nrows, ncols)

  write(*,'("Grid Sizes: ",I0," ",I0)') nrows, ncols
  write(*,'("Iterations: ",I0)') niterations
  write(*,*)
  write(*,'("Computing wavefield in device ..")')

  !$omp target data map(tofrom: next_base, prev_base) map(to: vel_base)
  kstart = omp_get_wtime()
  do k = 0, niterations - 1
    if (mod(k, 2) == 1) then
      call iso_2dfd_kernel(prev_base, next_base, vel_base, dtdivdxy, nrows, ncols)
    else
      call iso_2dfd_kernel(next_base, prev_base, vel_base, dtdivdxy, nrows, ncols)
    end if
  end do
  ktime = omp_get_wtime() - kstart
  !$omp end target data

  write(*,'("Total kernel execution time ",F0.6," (ms)")') ktime * 1.0e3_real64
  write(*,'("Average kernel execution time ",F0.6," (us)")') ktime * 1.0e6_real64 / real(niterations, real64)

  call write_binary("wavefield_snapshot.bin", next_base)

  write(*,'("Computing wavefield in CPU ..")')
  call initialize(prev_base, next_cpu, vel_base, nrows, ncols)

  cpu_start = omp_get_wtime()
  call iso_2dfd_iteration_cpu(next_cpu, prev_base, vel_base, dtdivdxy, nrows, ncols, niterations)
  cpu_time_ms = (omp_get_wtime() - cpu_start) * 1.0e3_real64

  write(*,'("CPU time: ",I0," ms")') int(cpu_time_ms)
  write(*,*)
  write(*,'("Check difference between final wavefields computed in device and host")')
  error = within_epsilon(next_base, next_cpu, nrows, ncols, half_length, 0.1_real32)
  if (error) then
    write(*,'("FAIL")')
  else
    write(*,'("PASS")')
  end if

  call write_binary("wavefield_snapshot_cpu.bin", next_cpu)
  write(*,'("Final wavefields (from device and CPU) written to disk")')
  write(*,'("Finished.  ")')

  if (error) stop 1

contains

  subroutine usage()
    write(*,'(" Incorrect parameters ")')
    write(*,'(" Usage: ./main n1 n2 Iterations ")')
    write(*,*)
    write(*,'(" n1 n2      : Grid sizes for the stencil ")')
    write(*,'(" Iterations : No. of timesteps. ")')
    stop 1
  end subroutine usage

  subroutine initialize(ptr_prev, ptr_next, ptr_vel, nr, nc)
    real(real32), intent(out) :: ptr_prev(0:), ptr_next(0:), ptr_vel(0:)
    integer, intent(in) :: nr, nc
    real(real32), parameter :: wavelet(0:11) = [ &
      0.016387336_real32, -0.041464937_real32, -0.067372555_real32, 0.386110067_real32, &
      0.812723635_real32, 0.416998396_real32, 0.076488599_real32, -0.059434419_real32, &
      0.023680172_real32, 0.005611435_real32, 0.001823209_real32, -0.000720549_real32]
    integer :: i, j, s, offset

    write(*,'("Initializing ... ")')
    ptr_prev = 0.0_real32
    ptr_next = 0.0_real32
    ptr_vel = 2250000.0_real32

    do s = 11, 0, -1
      do i = nr / 2 - s, nr / 2 + s - 1
        offset = i * nc
        do j = nc / 2 - s, nc / 2 + s - 1
          ptr_prev(offset + j) = wavelet(s)
        end do
      end do
    end do
  end subroutine initialize

  subroutine iso_2dfd_iteration_cpu(next_arr, prev_arr, vel, dtdiv, nr, nc, niters)
    real(real32), intent(inout) :: next_arr(0:), prev_arr(0:)
    real(real32), intent(in) :: vel(0:), dtdiv
    integer, intent(in) :: nr, nc, niters
    integer :: iter, i, j, gid
    real(real32) :: value
    logical :: write_next

    do iter = 0, niters - 1
      write_next = mod(iter, 2) == 0
      do i = 1, nr - half_length - 1
        do j = 1, nc - half_length - 1
          gid = j + i * nc
          if (write_next) then
            value = 0.0_real32
            value = value + prev_arr(gid + 1) - 2.0_real32 * prev_arr(gid) + prev_arr(gid - 1)
            value = value + prev_arr(gid + nc) - 2.0_real32 * prev_arr(gid) + prev_arr(gid - nc)
            value = value * dtdiv * vel(gid)
            next_arr(gid) = 2.0_real32 * prev_arr(gid) - next_arr(gid) + value
          else
            value = 0.0_real32
            value = value + next_arr(gid + 1) - 2.0_real32 * next_arr(gid) + next_arr(gid - 1)
            value = value + next_arr(gid + nc) - 2.0_real32 * next_arr(gid) + next_arr(gid - nc)
            value = value * dtdiv * vel(gid)
            prev_arr(gid) = 2.0_real32 * next_arr(gid) - prev_arr(gid) + value
          end if
        end do
      end do
    end do
  end subroutine iso_2dfd_iteration_cpu

  subroutine iso_2dfd_kernel(next_arr, prev_arr, vel, dtdiv, nr, nc)
    real(real32), intent(inout) :: next_arr(0:)
    real(real32), intent(in) :: prev_arr(0:), vel(0:), dtdiv
    integer, intent(in) :: nr, nc
    integer :: gid_row, gid_col, gid
    real(real32) :: value

    !$omp target teams distribute parallel do simd collapse(2) thread_limit(256)
    do gid_row = 0, nr - 1
      do gid_col = 0, nc - 1
        gid = gid_row * nc + gid_col
        if ((gid_col >= half_length .and. gid_col < nc - half_length) .and. &
            (gid_row >= half_length .and. gid_row < nr - half_length)) then
          value = 0.0_real32
          value = value + prev_arr(gid + 1) - 2.0_real32 * prev_arr(gid) + prev_arr(gid - 1)
          value = value + prev_arr(gid + nc) - 2.0_real32 * prev_arr(gid) + prev_arr(gid - nc)
          value = value * dtdiv * vel(gid)
          next_arr(gid) = 2.0_real32 * prev_arr(gid) - next_arr(gid) + value
        end if
      end do
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine iso_2dfd_kernel

  logical function within_epsilon(output, reference, dimx, dimy, radius, delta)
    real(real32), intent(in) :: output(0:), reference(0:), delta
    integer, intent(in) :: dimx, dimy, radius
    integer :: ix, iy, gid, unit
    real(real32) :: difference
    real(real64) :: norm2

    within_epsilon = .false.
    norm2 = 0.0_real64
    open(newunit=unit, file="error_diff.txt", status="replace", action="write")
    do iy = 0, dimy - 1
      do ix = 0, dimx - 1
        gid = ix + iy * dimx
        if (ix >= radius .and. ix < dimx - radius .and. iy >= radius .and. iy < dimy - radius) then
          difference = abs(reference(gid) - output(gid))
          norm2 = norm2 + real(difference * difference, real64)
          if (difference > delta) then
            within_epsilon = .true.
            write(unit,'(" ERROR: (",I0,",",I0,")",A,ES12.6E2," instead of ",ES12.6E2," (|e|=",ES12.6E2,")")') &
              ix, iy, char(9), output(gid), reference(gid), difference
          end if
        end if
      end do
    end do
    close(unit)
    if (within_epsilon) then
      write(*,'("error (Euclidean norm): ",ES15.9E2)') sqrt(norm2)
    end if
  end function within_epsilon

  subroutine write_binary(path, data)
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: data(0:)
    integer :: unit

    open(newunit=unit, file=path, access="stream", form="unformatted", status="replace", action="write")
    write(unit) data
    close(unit)
  end subroutine write_binary

end program iso2dfd
