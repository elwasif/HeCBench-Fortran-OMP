! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use omp_lib
  implicit none

  type :: float3
    real(real32) :: x, y, z, pad
  end type float3

  type :: float4
    real(real32) :: x, y, z, w
  end type float4

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

  character(len=256) :: arg0, arg
  integer :: width, height, repeat, num_pts, i
  type(float3), allocatable :: points(:)
  type(float4), allocatable :: normal_points(:), ref_normal_points(:)
  real(real32) :: sx, sy, sz, sw
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <width> <height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) repeat

  num_pts = width * height
  allocate(points(num_pts))
  allocate(normal_points(num_pts))
  allocate(ref_normal_points(num_pts))

  call c_srand(123_c_int)
  do i = 1, num_pts
    points(i)%x = real(modulo(c_rand(), width), real32)
    points(i)%y = real(modulo(c_rand(), height), real32)
    points(i)%z = real(modulo(c_rand(), 256_c_int), real32)
    points(i)%pad = 0.0_real32
  end do
  normal_points%x = 0.0_real32; normal_points%y = 0.0_real32
  normal_points%z = 0.0_real32; normal_points%w = 0.0_real32
  ref_normal_points%x = 0.0_real32; ref_normal_points%y = 0.0_real32
  ref_normal_points%z = 0.0_real32; ref_normal_points%w = 0.0_real32

  !$omp target data map(to: points(1:num_pts)) map(from: normal_points(1:num_pts))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call estimate_normals(points, normal_points, width, height, num_pts)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', &
    elapsed / real(repeat, real64), ' (s)'

  call estimate_normals_cpu(points, ref_normal_points, width, height, num_pts)
  if (maxval(abs(normal_points%x - ref_normal_points%x)) > 1.0e-5_real32 .or. &
      maxval(abs(normal_points%y - ref_normal_points%y)) > 1.0e-5_real32 .or. &
      maxval(abs(normal_points%z - ref_normal_points%z)) > 1.0e-5_real32 .or. &
      maxval(abs(normal_points%w - ref_normal_points%w)) > 1.0e-5_real32) then
    error stop 'Fortran reference validation failed'
  end if

  sx = sum(normal_points%x)
  sy = sum(normal_points%y)
  sz = sum(normal_points%z)
  sw = sum(normal_points%w)
  write(*,'(A,F0.6,A,F0.6,A,F0.6,A,F0.6)') 'Checksum: x=', sx, ' y=', sy, &
    ' z=', sz, ' w=', sw

  deallocate(points, normal_points, ref_normal_points)

contains

  subroutine estimate_normals(points, normal_points, width, height, num_pts)
    integer, intent(in) :: width, height, num_pts
    type(float3), intent(in) :: points(:)
    type(float4), intent(out) :: normal_points(:)
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256)
    do idx = 0, num_pts - 1
      call normal_estimate(points, idx, width, height, normal_points(idx + 1))
    end do
    !$omp end target teams distribute parallel do
  end subroutine estimate_normals

  subroutine estimate_normals_cpu(points, normal_points, width, height, num_pts)
    integer, intent(in) :: width, height, num_pts
    type(float3), intent(in) :: points(:)
    type(float4), intent(out) :: normal_points(:)
    integer :: idx

    do idx = 0, num_pts - 1
      call normal_estimate(points, idx, width, height, normal_points(idx + 1))
    end do
  end subroutine estimate_normals_cpu

  subroutine normal_estimate(points, idx, width, height, normal_point)
    type(float3), intent(in) :: points(:)
    integer, intent(in) :: idx, width, height
    type(float4), intent(out) :: normal_point
    integer :: p, x_idx, y_idx
    logical :: west_valid, east_valid, north_valid, south_valid
    real(real32) :: qx, qy, qz, hx, hy, hz, vx, vy, vz
    real(real32) :: nx0, ny0, nz0, curvature, len, inv_len, dotv

    p = idx + 1
    qx = points(p)%x; qy = points(p)%y; qz = points(p)%z
    if (ieee_is_nan(qz)) then
      normal_point%x = 0.0_real32; normal_point%y = 0.0_real32
      normal_point%z = 0.0_real32; normal_point%w = 0.0_real32
      return
    end if

    x_idx = modulo(idx, width)
    y_idx = idx / width

    west_valid = .false.
    east_valid = .false.
    north_valid = .false.
    south_valid = .false.
    if (x_idx > 1) west_valid = (.not. ieee_is_nan(points(p - 1)%z)) .and. &
      abs(points(p - 1)%z - qz) < 200.0_real32
    if (x_idx < width - 1) east_valid = (.not. ieee_is_nan(points(p + 1)%z)) .and. &
      abs(points(p + 1)%z - qz) < 200.0_real32
    if (y_idx > 1) north_valid = (.not. ieee_is_nan(points(p - width)%z)) .and. &
      abs(points(p - width)%z - qz) < 200.0_real32
    if (y_idx < height - 1) south_valid = (.not. ieee_is_nan(points(p + width)%z)) .and. &
      abs(points(p + width)%z - qz) < 200.0_real32

    if (west_valid .and. east_valid) then
      hx = points(p + 1)%x - points(p - 1)%x
      hy = points(p + 1)%y - points(p - 1)%y
      hz = points(p + 1)%z - points(p - 1)%z
    else if (west_valid .and. .not. east_valid) then
      hx = qx - points(p - 1)%x
      hy = qy - points(p - 1)%y
      hz = qz - points(p - 1)%z
    else if (.not. west_valid .and. east_valid) then
      hx = points(p + 1)%x - qx
      hy = points(p + 1)%y - qy
      hz = points(p + 1)%z - qz
    else
      normal_point%x = 0.0_real32; normal_point%y = 0.0_real32
      normal_point%z = 0.0_real32; normal_point%w = 1.0_real32
      return
    end if

    if (south_valid .and. north_valid) then
      vx = points(p - width)%x - points(p + width)%x
      vy = points(p - width)%y - points(p + width)%y
      vz = points(p - width)%z - points(p + width)%z
    else if (south_valid .and. .not. north_valid) then
      vx = qx - points(p + width)%x
      vy = qy - points(p + width)%y
      vz = qz - points(p + width)%z
    else if (.not. south_valid .and. north_valid) then
      vx = points(p - width)%x - qx
      vy = points(p - width)%y - qy
      vz = points(p - width)%z - qz
    else
      normal_point%x = 0.0_real32; normal_point%y = 0.0_real32
      normal_point%z = 0.0_real32; normal_point%w = 1.0_real32
      return
    end if

    nx0 = hy * vz - hz * vy
    ny0 = hz * vx - hx * vz
    nz0 = hx * vy - hy * vx
    len = sqrt(nx0 * nx0 + ny0 * ny0 + nz0 * nz0)
    curvature = merge(1.0_real32, 0.0_real32, abs(hz) > 0.04_real32 .or. &
      abs(vz) > 0.04_real32 .or. .not. west_valid .or. .not. east_valid .or. &
      .not. north_valid .or. .not. south_valid)

    inv_len = 1.0_real32 / len
    normal_point%x = nx0 * inv_len
    normal_point%y = ny0 * inv_len
    normal_point%z = nz0 * inv_len
    dotv = qx * normal_point%x + qy * normal_point%y + qz * normal_point%z
    if (dotv > 0.0_real32) then
      normal_point%x = -normal_point%x
      normal_point%y = -normal_point%y
      normal_point%z = -normal_point%z
    end if
    normal_point%w = curvature
  end subroutine normal_estimate

end program main
