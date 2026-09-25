! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: width = 256
  integer, parameter :: height = 256
  integer, parameter :: nsubsamples = 2
  integer, parameter :: nao_samples = 8

  type :: vec
    real(real32) :: x
    real(real32) :: y
    real(real32) :: z
  end type vec

  type :: isect
    real(real32) :: t
    type(vec) :: p
    type(vec) :: n
    integer(int32) :: hit
  end type isect

  type :: sphere
    type(vec) :: center
    real(real32) :: radius
  end type sphere

  type :: plane
    type(vec) :: p
    type(vec) :: n
  end type plane

  type :: ray
    type(vec) :: org
    type(vec) :: dir
  end type ray

  type :: rng_state
    integer(int32) :: x
  end type rng_state

  integer :: argc, ios, loopmax, iter, mismatches
  integer(int8), allocatable :: img(:), ref_img(:)
  type(sphere) :: spheres(0:2)
  type(plane) :: scene_plane
  real(real64) :: total_time
  character(len=256) :: arg, program_name

  argc = command_argument_count()
  if (argc /= 1) then
    call get_command_argument(0, program_name)
    write(*, '(A,A,A)') 'Usage: ', trim(program_name), ' <iterations>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=ios) loopmax
  if (ios /= 0 .or. loopmax <= 0) stop 1

  allocate(img(0:width * height * 3 - 1), ref_img(0:width * height * 3 - 1))
  img = 0_int8
  ref_img = 0_int8

  call init_scene(spheres, scene_plane)

  total_time = 0.0_real64
  do iter = 1, loopmax
    total_time = total_time + render_device(img, width, height, nsubsamples, spheres, scene_plane)
  end do

  call render_host(ref_img, width, height, nsubsamples, spheres, scene_plane)
  mismatches = count(img /= ref_img)
  if (mismatches /= 0) then
    write(*, '(A,I0)') 'FAIL: image mismatches=', mismatches
    stop 1
  end if

  write(*, '(A,F0.6,A)') 'Average kernel time: ', total_time * 1.0e6_real64 / real(loopmax, real64), ' usec.'

  call saveppm('ao.ppm', width, height, img)

  deallocate(img, ref_img)

contains

  real(real32) function vdot(v0, v1)
    type(vec), intent(in) :: v0, v1
    vdot = v0%x * v1%x + v0%y * v1%y + v0%z * v1%z
  end function vdot

  subroutine vcross(c, v0, v1)
    type(vec), intent(out) :: c
    type(vec), intent(in) :: v0, v1
    c%x = v0%y * v1%z - v0%z * v1%y
    c%y = v0%z * v1%x - v0%x * v1%z
    c%z = v0%x * v1%y - v0%y * v1%x
  end subroutine vcross

  subroutine vnormalize(c)
    type(vec), intent(inout) :: c
    real(real32) :: length

    length = sqrt(vdot(c, c))
    if (abs(length) > 1.0e-17_real32) then
      c%x = c%x / length
      c%y = c%y / length
      c%z = c%z / length
    end if
  end subroutine vnormalize

  subroutine ray_sphere_intersect(hit_rec, in_ray, in_sphere)
    type(isect), intent(inout) :: hit_rec
    type(ray), intent(in) :: in_ray
    type(sphere), intent(in) :: in_sphere
    type(vec) :: rs
    real(real32) :: b, c, d, t

    rs%x = in_ray%org%x - in_sphere%center%x
    rs%y = in_ray%org%y - in_sphere%center%y
    rs%z = in_ray%org%z - in_sphere%center%z

    b = vdot(rs, in_ray%dir)
    c = vdot(rs, rs) - in_sphere%radius * in_sphere%radius
    d = b * b - c

    if (d > 0.0_real32) then
      t = -b - sqrt(d)
      if (t > 0.0_real32 .and. t < hit_rec%t) then
        hit_rec%t = t
        hit_rec%hit = 1_int32
        hit_rec%p%x = in_ray%org%x + in_ray%dir%x * t
        hit_rec%p%y = in_ray%org%y + in_ray%dir%y * t
        hit_rec%p%z = in_ray%org%z + in_ray%dir%z * t
        hit_rec%n%x = hit_rec%p%x - in_sphere%center%x
        hit_rec%n%y = hit_rec%p%y - in_sphere%center%y
        hit_rec%n%z = hit_rec%p%z - in_sphere%center%z
        call vnormalize(hit_rec%n)
      end if
    end if
  end subroutine ray_sphere_intersect

  subroutine ray_plane_intersect(hit_rec, in_ray, in_plane)
    type(isect), intent(inout) :: hit_rec
    type(ray), intent(in) :: in_ray
    type(plane), intent(in) :: in_plane
    real(real32) :: d, v, t

    d = -vdot(in_plane%p, in_plane%n)
    v = vdot(in_ray%dir, in_plane%n)
    if (abs(v) < 1.0e-17_real32) return

    t = -(vdot(in_ray%org, in_plane%n) + d) / v
    if (t > 0.0_real32 .and. t < hit_rec%t) then
      hit_rec%t = t
      hit_rec%hit = 1_int32
      hit_rec%p%x = in_ray%org%x + in_ray%dir%x * t
      hit_rec%p%y = in_ray%org%y + in_ray%dir%y * t
      hit_rec%p%z = in_ray%org%z + in_ray%dir%z * t
      hit_rec%n = in_plane%n
    end if
  end subroutine ray_plane_intersect

  subroutine ortho_basis(basis, n)
    type(vec), intent(out) :: basis(0:2)
    type(vec), intent(in) :: n

    basis(2) = n
    basis(1)%x = 0.0_real32
    basis(1)%y = 0.0_real32
    basis(1)%z = 0.0_real32

    if (n%x < 0.6_real32 .and. n%x > -0.6_real32) then
      basis(1)%x = 1.0_real32
    else if (n%y < 0.6_real32 .and. n%y > -0.6_real32) then
      basis(1)%y = 1.0_real32
    else if (n%z < 0.6_real32 .and. n%z > -0.6_real32) then
      basis(1)%z = 1.0_real32
    else
      basis(1)%x = 1.0_real32
    end if

    call vcross(basis(0), basis(1), basis(2))
    call vnormalize(basis(0))
    call vcross(basis(1), basis(2), basis(0))
    call vnormalize(basis(1))
  end subroutine ortho_basis

  integer(int32) function rng_next(rng)
    type(rng_state), intent(inout) :: rng

    rng%x = ieor(rng%x, shiftr(rng%x, 6))
    rng%x = ieor(rng%x, shiftl(rng%x, 17))
    rng%x = ieor(rng%x, shiftr(rng%x, 9))
    rng_next = rng%x
  end function rng_next

  real(real32) function rng_uniform(rng)
    type(rng_state), intent(inout) :: rng
    integer(int32), parameter :: fmask = int(z'007fffff', int32)
    integer(int32), parameter :: one_bits = int(z'3f800000', int32)
    integer(int32) :: bits

    bits = ior(iand(rng_next(rng), fmask), one_bits)
    rng_uniform = transfer(bits, rng_uniform) - 1.0_real32
  end function rng_uniform

  subroutine ambient_occlusion(col, hit_rec, spheres, scene_plane, rng)
    type(vec), intent(out) :: col
    type(isect), intent(in) :: hit_rec
    type(sphere), intent(in) :: spheres(0:2)
    type(plane), intent(in) :: scene_plane
    type(rng_state), intent(inout) :: rng
    integer :: i, j
    real(real32), parameter :: eps = 0.0001_real32
    real(real32), parameter :: pi = 3.14159265358979323846_real32
    real(real32) :: occlusion, theta, phi, lx, ly, lz, rx, ry, rz
    type(vec) :: p, basis(0:2)
    type(ray) :: occ_ray
    type(isect) :: occ_hit

    p%x = hit_rec%p%x + eps * hit_rec%n%x
    p%y = hit_rec%p%y + eps * hit_rec%n%y
    p%z = hit_rec%p%z + eps * hit_rec%n%z

    call ortho_basis(basis, hit_rec%n)

    occlusion = 0.0_real32
    do j = 0, nao_samples - 1
      do i = 0, nao_samples - 1
        theta = sqrt(rng_uniform(rng))
        phi = 2.0_real32 * pi * rng_uniform(rng)
        lx = cos(phi) * theta
        ly = sin(phi) * theta
        lz = sqrt(1.0_real32 - theta * theta)

        rx = lx * basis(0)%x + ly * basis(1)%x + lz * basis(2)%x
        ry = lx * basis(0)%y + ly * basis(1)%y + lz * basis(2)%y
        rz = lx * basis(0)%z + ly * basis(1)%z + lz * basis(2)%z

        occ_ray%org = p
        occ_ray%dir%x = rx
        occ_ray%dir%y = ry
        occ_ray%dir%z = rz

        occ_hit%t = 1.0e17_real32
        occ_hit%hit = 0_int32
        call ray_sphere_intersect(occ_hit, occ_ray, spheres(0))
        call ray_sphere_intersect(occ_hit, occ_ray, spheres(1))
        call ray_sphere_intersect(occ_hit, occ_ray, spheres(2))
        call ray_plane_intersect(occ_hit, occ_ray, scene_plane)

        if (occ_hit%hit /= 0_int32) occlusion = occlusion + 1.0_real32
      end do
    end do

    occlusion = real(nao_samples * nao_samples, real32) - occlusion
    occlusion = occlusion / real(nao_samples * nao_samples, real32)
    col%x = occlusion
    col%y = occlusion
    col%z = occlusion
  end subroutine ambient_occlusion

  integer(int8) function my_clamp(f)
    real(real32), intent(in) :: f
    integer(int32) :: i

    i = int(f * 255.5_real32, int32)
    if (i < 0_int32) i = 0_int32
    if (i > 255_int32) i = 255_int32
    my_clamp = pack_u8(i)
  end function my_clamp

  integer(int8) function pack_u8(value)
    integer(int32), intent(in) :: value
    integer(int32) :: masked

    masked = iand(value, 255_int32)
    if (masked > 127_int32) then
      pack_u8 = int(masked - 256_int32, int8)
    else
      pack_u8 = int(masked, int8)
    end if
  end function pack_u8

  subroutine init_scene(spheres, scene_plane)
    type(sphere), intent(out) :: spheres(0:2)
    type(plane), intent(out) :: scene_plane

    spheres(0)%center%x = -2.0_real32
    spheres(0)%center%y = 0.0_real32
    spheres(0)%center%z = -3.5_real32
    spheres(0)%radius = 0.5_real32

    spheres(1)%center%x = -0.5_real32
    spheres(1)%center%y = 0.0_real32
    spheres(1)%center%z = -3.0_real32
    spheres(1)%radius = 0.5_real32

    spheres(2)%center%x = 1.0_real32
    spheres(2)%center%y = 0.0_real32
    spheres(2)%center%z = -2.2_real32
    spheres(2)%radius = 0.5_real32

    scene_plane%p%x = 0.0_real32
    scene_plane%p%y = -0.5_real32
    scene_plane%p%z = 0.0_real32

    scene_plane%n%x = 0.0_real32
    scene_plane%n%y = 1.0_real32
    scene_plane%n%z = 0.0_real32
  end subroutine init_scene

  subroutine saveppm(fname, w, h, img)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: w, h
    integer(int8), intent(in) :: img(0:)
    integer :: unit, ios
    character(len=64) :: header

    open(newunit=unit, file=fname, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A)') 'Failed to open the file ', trim(fname)
      stop 1
    end if

    write(header, '("P6",A,I0,1X,I0,A,"255",A)') achar(10), w, h, achar(10), achar(10)
    write(unit) header(1:len_trim(header))
    write(unit) img(0:w * h * 3 - 1)
    close(unit)
  end subroutine saveppm

  real(real64) function render_device(img, w, h, samples, spheres, scene_plane)
    integer(int8), intent(inout) :: img(0:)
    integer, intent(in) :: w, h, samples
    type(sphere), intent(in) :: spheres(0:2)
    type(plane), intent(in) :: scene_plane
    integer :: x, y, u, v, idx
    real(real32) :: px, py, s0, s1, s2
    real(real64) :: start_time, end_time
    type(ray) :: primary_ray
    type(isect) :: hit_rec
    type(vec) :: col
    type(rng_state) :: rng

    !$omp target data map(from: img(0:w*h*3-1)) map(to: spheres(0:2), scene_plane)
    start_time = omp_get_wtime()
    !$omp target teams distribute parallel do simd collapse(2) thread_limit(256) &
    !$omp& private(u, v, idx, px, py, s0, s1, s2, primary_ray, hit_rec, col, rng)
    do x = 0, w - 1
      do y = 0, h - 1
        rng%x = int(y * w + x, int32)
        s0 = 0.0_real32
        s1 = 0.0_real32
        s2 = 0.0_real32

        do v = 0, samples - 1
          do u = 0, samples - 1
            px = (real(x, real32) + real(u, real32) / real(samples, real32) - real(w, real32) / 2.0_real32) / &
                (real(w, real32) / 2.0_real32)
            py = -(real(y, real32) + real(v, real32) / real(samples, real32) - real(h, real32) / 2.0_real32) / &
                (real(h, real32) / 2.0_real32)

            primary_ray%org%x = 0.0_real32
            primary_ray%org%y = 0.0_real32
            primary_ray%org%z = 0.0_real32
            primary_ray%dir%x = px
            primary_ray%dir%y = py
            primary_ray%dir%z = -1.0_real32
            call vnormalize(primary_ray%dir)

            hit_rec%t = 1.0e17_real32
            hit_rec%hit = 0_int32
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(0))
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(1))
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(2))
            call ray_plane_intersect(hit_rec, primary_ray, scene_plane)

            if (hit_rec%hit /= 0_int32) then
              call ambient_occlusion(col, hit_rec, spheres, scene_plane, rng)
              s0 = s0 + col%x
              s1 = s1 + col%y
              s2 = s2 + col%z
            end if
          end do
        end do

        idx = 3 * (y * w + x)
        img(idx) = my_clamp(s0 / real(samples * samples, real32))
        img(idx + 1) = my_clamp(s1 / real(samples * samples, real32))
        img(idx + 2) = my_clamp(s2 / real(samples * samples, real32))
      end do
    end do
    !$omp end target teams distribute parallel do simd
    end_time = omp_get_wtime()
    !$omp end target data
    render_device = end_time - start_time
  end function render_device

  subroutine render_host(img, w, h, samples, spheres, scene_plane)
    integer(int8), intent(out) :: img(0:)
    integer, intent(in) :: w, h, samples
    type(sphere), intent(in) :: spheres(0:2)
    type(plane), intent(in) :: scene_plane
    integer :: x, y, u, v, idx
    real(real32) :: px, py, s0, s1, s2
    type(ray) :: primary_ray
    type(isect) :: hit_rec
    type(vec) :: col
    type(rng_state) :: rng

    do x = 0, w - 1
      do y = 0, h - 1
        rng%x = int(y * w + x, int32)
        s0 = 0.0_real32
        s1 = 0.0_real32
        s2 = 0.0_real32

        do v = 0, samples - 1
          do u = 0, samples - 1
            px = (real(x, real32) + real(u, real32) / real(samples, real32) - real(w, real32) / 2.0_real32) / &
                (real(w, real32) / 2.0_real32)
            py = -(real(y, real32) + real(v, real32) / real(samples, real32) - real(h, real32) / 2.0_real32) / &
                (real(h, real32) / 2.0_real32)

            primary_ray%org%x = 0.0_real32
            primary_ray%org%y = 0.0_real32
            primary_ray%org%z = 0.0_real32
            primary_ray%dir%x = px
            primary_ray%dir%y = py
            primary_ray%dir%z = -1.0_real32
            call vnormalize(primary_ray%dir)

            hit_rec%t = 1.0e17_real32
            hit_rec%hit = 0_int32
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(0))
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(1))
            call ray_sphere_intersect(hit_rec, primary_ray, spheres(2))
            call ray_plane_intersect(hit_rec, primary_ray, scene_plane)

            if (hit_rec%hit /= 0_int32) then
              call ambient_occlusion(col, hit_rec, spheres, scene_plane, rng)
              s0 = s0 + col%x
              s1 = s1 + col%y
              s2 = s2 + col%z
            end if
          end do
        end do

        idx = 3 * (y * w + x)
        img(idx) = my_clamp(s0 / real(samples * samples, real32))
        img(idx + 1) = my_clamp(s1 / real(samples * samples, real32))
        img(idx + 2) = my_clamp(s2 / real(samples * samples, real32))
      end do
    end do
  end subroutine render_host

end program main
