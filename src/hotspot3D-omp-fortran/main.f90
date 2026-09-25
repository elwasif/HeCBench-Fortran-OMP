! SPDX-License-Identifier: CC0-1.0
program hotspot3d
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_double, c_float, c_int
  use omp_lib
  implicit none

  interface
    subroutine compute_temp_cpu_c(pin, tin, tout, nx, ny, nz, cap, rx, ry, rz, dt, amb_temp, numiter) bind(C)
      import :: c_float, c_int
      real(c_float), intent(in) :: pin(*), tin(*)
      real(c_float), intent(out) :: tout(*)
      integer(c_int), value :: nx, ny, nz, numiter
      real(c_float), value :: cap, rx, ry, rz, dt, amb_temp
    end subroutine compute_temp_cpu_c

    subroutine hotspot3d_target_c(tin, pin, tout, num_cols, num_rows, layers, iterations, ce, cw, cn, cs, &
                                  ct, cb, cc, step_div_cap, kernel_time) bind(C)
      import :: c_double, c_float, c_int
      real(c_float), intent(inout) :: tin(*), tout(*)
      real(c_float), intent(in) :: pin(*)
      integer(c_int), value :: num_cols, num_rows, layers, iterations
      real(c_float), value :: ce, cw, cn, cs, ct, cb, cc, step_div_cap
      real(c_double), intent(out) :: kernel_time
    end subroutine hotspot3d_target_c

    function accuracy_c(arr1, arr2, len) result(rms) bind(C)
      import :: c_float, c_int
      real(c_float), intent(in) :: arr1(*), arr2(*)
      integer(c_int), value :: len
      real(c_float) :: rms
    end function accuracy_c

    subroutine print_rms_c(rms) bind(C)
      import :: c_float
      real(c_float), value :: rms
    end subroutine print_rms_c
  end interface

  real(real32), parameter :: max_pd = 3.0e6_real32
  real(real32), parameter :: tol = 0.001_real32
  real(real32), parameter :: precision = 0.001_real32
  real(real32), parameter :: spec_heat_si = 1.75e6_real32
  real(real32), parameter :: k_si = 100.0_real32
  real(real32), parameter :: factor_chip = 0.5_real32
  real(real32), parameter :: t_chip = 0.0005_real32
  real(real32), parameter :: chip_height = 0.016_real32
  real(real32), parameter :: chip_width = 0.016_real32
  real(real32), parameter :: amb_temp = 80.0_real32

  integer :: num_cols, num_rows, layers, iterations, size
  character(len=256) :: arg, pfile, tfile, ofile
  real(real32) :: dx, dy, dz, cap, rx, ry, rz, max_slope, dt, rms
  real(real32) :: ce, cw, cn, cs, ct, cb, cc, step_div_cap
  real(real32), allocatable :: tin(:), pin(:), tcopy(:), tout(:), answer(:)
  real(real64) :: start_time, stop_time, kernel_time
  logical :: sel_is_tin

  if (command_argument_count() /= 6) then
    call usage()
  end if

  call get_command_argument(1, arg)
  read(arg, *) num_cols
  num_rows = num_cols
  call get_command_argument(2, arg)
  read(arg, *) layers
  call get_command_argument(3, arg)
  read(arg, *) iterations
  call get_command_argument(4, pfile)
  call get_command_argument(5, tfile)
  call get_command_argument(6, ofile)

  dx = chip_height / real(num_rows, real32)
  dy = chip_width / real(num_cols, real32)
  dz = t_chip / real(layers, real32)
  cap = factor_chip * spec_heat_si * t_chip * dx * dy
  rx = dy / (2.0_real32 * k_si * t_chip * dx)
  ry = dx / (2.0_real32 * k_si * t_chip * dy)
  rz = dz / (k_si * dx * dy)
  max_slope = max_pd / (factor_chip * t_chip * spec_heat_si)
  dt = precision / max_slope
  step_div_cap = dt / cap
  ce = step_div_cap / rx
  cw = ce
  cn = step_div_cap / ry
  cs = cn
  ct = step_div_cap / rz
  cb = ct
  cc = 1.0_real32 - (2.0_real32 * ce + 2.0_real32 * cn + 3.0_real32 * ct)

  size = num_cols * num_rows * layers
  allocate(tin(0:size-1), pin(0:size-1), tcopy(0:size-1), tout(0:size-1), answer(0:size-1))
  tin = 0.0_real32
  pin = 0.0_real32
  tout = 0.0_real32
  answer = 0.0_real32

  call read_input(tin, num_rows, num_cols, layers, trim(tfile))
  call read_input(pin, num_rows, num_cols, layers, trim(pfile))
  tcopy = tin

  start_time = omp_get_wtime()
  call hotspot3d_target_c(tin, pin, tout, num_cols, num_rows, layers, iterations, ce, cw, cn, cs, ct, cb, cc, &
                          step_div_cap, kernel_time)
  sel_is_tin = .false.

  stop_time = omp_get_wtime()

  call compute_temp_cpu_c(pin, tcopy, answer, num_cols, num_rows, layers, cap, rx, ry, rz, dt, amb_temp, iterations)

  write(*,'("Average kernel execution time ",F0.6," (us)")') kernel_time * 1.0e6_real64 / real(iterations, real64)
  write(*,'("Device offloading time: ",F5.3," (s)")') stop_time - start_time
  if (sel_is_tin) then
    rms = accuracy_c(tin, answer, size)
  else
    rms = accuracy_c(tout, answer, size)
  end if
  call print_rms_c(rms)
  call write_output(tout, num_rows, num_cols, layers, trim(ofile))

contains

  subroutine usage()
    write(0,'("Usage: ./main <rows/cols> <layers> <iterations> <powerFile> <tempFile> <outputFile>")')
    stop 1
  end subroutine usage

  subroutine read_input(vect, grid_rows, grid_cols, nz, path)
    real(real32), intent(out) :: vect(0:)
    integer, intent(in) :: grid_rows, grid_cols, nz
    character(len=*), intent(in) :: path
    integer :: i, j, k, unit, ios
    real(real32) :: val

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'("Error: The file was not opened")')
      stop 1
    end if

    do i = 0, grid_rows - 1
      do j = 0, grid_cols - 1
        do k = 0, nz - 1
          read(unit, *, iostat=ios) val
          if (ios /= 0) then
            write(*,'("Error: invalid file format")')
            stop 1
          end if
          vect(i * grid_cols + j + k * grid_rows * grid_cols) = val
        end do
      end do
    end do
    close(unit)
  end subroutine read_input

  subroutine write_output(vect, grid_rows, grid_cols, nz, path)
    real(real32), intent(in) :: vect(0:)
    integer, intent(in) :: grid_rows, grid_cols, nz
    character(len=*), intent(in) :: path
    integer :: i, j, k, idx, unit

    open(newunit=unit, file=path, status='replace', action='write')
    idx = 0
    do i = 0, grid_rows - 1
      do j = 0, grid_cols - 1
        do k = 0, nz - 1
          write(unit,'(I0,A,G0)') idx, char(9), vect(i * grid_cols + j + k * grid_rows * grid_cols)
          idx = idx + 1
        end do
      end do
    end do
    close(unit)
  end subroutine write_output

  subroutine swap_arrays(a, b)
    real(real32), allocatable, intent(inout) :: a(:), b(:)
    real(real32), allocatable :: tmp(:)

    call move_alloc(a, tmp)
    call move_alloc(b, a)
    call move_alloc(tmp, b)
  end subroutine swap_arrays

  subroutine hotspot_step(tin, pin, tout, nx, ny, nz, ce, cw, cn, cs, ct, cb, cc, step_div_cap)
    real(real32), intent(in) :: tin(0:), pin(0:)
    real(real32), intent(out) :: tout(0:)
    integer, intent(in) :: nx, ny, nz
    real(real32), intent(in) :: ce, cw, cn, cs, ct, cb, cc, step_div_cap
    integer :: x, y, k, c, w, e, n, s, xy
    real(real32) :: temp1, temp2, temp3

    !$omp target teams distribute parallel do collapse(2) thread_limit(256)
    do y = 0, ny - 1
      do x = 0, nx - 1
        c = x + y * nx
        xy = nx * ny

        w = merge(c, c - 1, x == 0)
        e = merge(c, c + 1, x == nx - 1)
        n = merge(c, c - nx, y == 0)
        s = merge(c, c + nx, y == ny - 1)

        temp1 = tin(c)
        temp2 = tin(c)
        temp3 = tin(c + xy)
        tout(c) = ((((((((cc * temp2 + cw * tin(w)) + ce * tin(e)) + cs * tin(s)) + &
                  cn * tin(n)) + cb * temp1) + ct * temp3) + step_div_cap * pin(c)) + &
                  ct * amb_temp)
        c = c + xy
        w = w + xy
        e = e + xy
        n = n + xy
        s = s + xy

        do k = 1, nz - 2
          temp1 = temp2
          temp2 = temp3
          temp3 = tin(c + xy)
          tout(c) = ((((((((cc * temp2 + cw * tin(w)) + ce * tin(e)) + cs * tin(s)) + &
                    cn * tin(n)) + cb * temp1) + ct * temp3) + step_div_cap * pin(c)) + &
                    ct * amb_temp)
          c = c + xy
          w = w + xy
          e = e + xy
          n = n + xy
          s = s + xy
        end do

        temp1 = temp2
        temp2 = temp3
        tout(c) = ((((((((cc * temp2 + cw * tin(w)) + ce * tin(e)) + cs * tin(s)) + &
                  cn * tin(n)) + cb * temp1) + ct * temp3) + step_div_cap * pin(c)) + &
                  ct * amb_temp)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine hotspot_step

  subroutine compute_temp_cpu(pin, tin0, answer, nx, ny, nz, cap, rx, ry, rz, dt, ambient, numiter)
    real(real32), intent(in) :: pin(0:), tin0(0:)
    real(real32), intent(out) :: answer(0:)
    integer, intent(in) :: nx, ny, nz, numiter
    real(real32), intent(in) :: cap, rx, ry, rz, dt, ambient
    real(real32), allocatable :: current(:), next(:)
    real(real32) :: lce, lcw, lcn, lcs, lct, lcb, lcc, lstep
    integer :: iter, x, y, z, c, w, e, n, s, b, t, len

    len = nx * ny * nz
    allocate(current(0:len-1), next(0:len-1))
    current = tin0
    next = 0.0_real32
    lstep = dt / cap
    lce = lstep / rx
    lcw = lce
    lcn = lstep / ry
    lcs = lcn
    lct = lstep / rz
    lcb = lct
    lcc = 1.0_real32 - (2.0_real32 * lce + 2.0_real32 * lcn + 3.0_real32 * lct)

    do iter = 1, numiter
      do z = 0, nz - 1
        do y = 0, ny - 1
          do x = 0, nx - 1
            c = x + y * nx + z * nx * ny
            w = merge(c, c - 1, x == 0)
            e = merge(c, c + 1, x == nx - 1)
            n = merge(c, c - nx, y == 0)
            s = merge(c, c + nx, y == ny - 1)
            b = merge(c, c - nx * ny, z == 0)
            t = merge(c, c + nx * ny, z == nz - 1)
            next(c) = ((((((((current(c) * lcc + current(n) * lcn) + current(s) * lcs) + &
                      current(e) * lce) + current(w) * lcw) + current(t) * lct) + &
                      current(b) * lcb) + lstep * pin(c)) + lct * ambient)
          end do
        end do
      end do
      call swap_arrays(current, next)
    end do

    if (mod(numiter, 2) == 1) then
      answer = current
    else
      answer = next
    end if
    deallocate(current, next)
  end subroutine compute_temp_cpu

  real(real32) function accuracy(arr1, arr2, len)
    real(real32), intent(in) :: arr1(0:), arr2(0:)
    integer, intent(in) :: len
    integer :: i
    real(real32) :: err

    err = 0.0_real32
    do i = 0, len - 1
      err = err + (arr1(i) - arr2(i)) * (arr1(i) - arr2(i))
    end do
    accuracy = sqrt(err / real(len, real32))
  end function accuracy

end program hotspot3d
