! SPDX-License-Identifier: CC0-1.0
program bilateral_main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  real(real32), parameter :: pi = 3.14159265358979323846_real32
  integer :: width, height, img_size, repeat, i
  real(real32) :: variance_i, variance_spatial, a_square
  real(real32), allocatable :: src(:), dst(:), ref(:)
  real(real64) :: start_time, elapsed_ms
  character(len=64) :: arg
  logical :: ok

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

  if (command_argument_count() /= 5) then
    write(*,'("Usage: ./main <image width> <image height> <intensity> <spatial> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) width
  call get_command_argument(2, arg)
  read(arg, *) height
  call get_command_argument(3, arg)
  read(arg, *) variance_i
  call get_command_argument(4, arg)
  read(arg, *) variance_spatial
  call get_command_argument(5, arg)
  read(arg, *) repeat

  img_size = width * height
  a_square = 0.5_real32 / (variance_i * pi)

  allocate(src(0:img_size-1), dst(0:img_size-1), ref(0:img_size-1))
  call c_srand(123_c_int)
  do i = 0, img_size - 1
    src(i) = real(mod(c_rand(), 256_c_int), real32)
  end do

  ok = .true.
  !$omp target data map(to: src) map(alloc: dst)
  start_time = omp_get_wtime()
  do i = 1, repeat
    call bilateral_filter(3, src, dst, width, height, a_square, variance_i, variance_spatial)
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (3x3) ",F0.6," (ms)")') elapsed_ms
  !$omp target update from(dst)
  call bilateral_reference(3, src, ref, width, height, a_square, variance_i, variance_spatial)
  ok = ok .and. check_values(dst, ref, img_size)

  start_time = omp_get_wtime()
  do i = 1, repeat
    call bilateral_filter(6, src, dst, width, height, a_square, variance_i, variance_spatial)
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (6x6) ",F0.6," (ms)")') elapsed_ms
  !$omp target update from(dst)
  call bilateral_reference(6, src, ref, width, height, a_square, variance_i, variance_spatial)
  ok = ok .and. check_values(dst, ref, img_size)

  start_time = omp_get_wtime()
  do i = 1, repeat
    call bilateral_filter(9, src, dst, width, height, a_square, variance_i, variance_spatial)
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64 / real(repeat, real64)
  write(*,'("Average kernel execution time (9x9) ",F0.6," (ms)")') elapsed_ms
  !$omp target update from(dst)
  !$omp end target data
  call bilateral_reference(9, src, ref, width, height, a_square, variance_i, variance_spatial)
  ok = ok .and. check_values(dst, ref, img_size)

  if (ok) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

contains

  subroutine bilateral_filter(radius, input, output, width, height, a_square, variance_i, variance_spatial)
    integer, intent(in) :: radius, width, height
    real(real32), intent(in) :: input(0:), a_square, variance_i, variance_spatial
    real(real32), intent(out) :: output(0:)
    integer :: idx, idy, ioff, joff, id, idk, idl, id_w
    real(real32) :: center, neighbor, range_term, spatial_term, weight
    real(real32) :: result, normalization

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) &
    !$omp& private(ioff, joff, id, idk, idl, id_w, center, neighbor, range_term, spatial_term, weight, result, normalization)
    do idy = 0, height - 1
      do idx = 0, width - 1
        id = idy * width + idx
        center = input(id)
        result = 0.0_real32
        normalization = 0.0_real32
        do ioff = -radius, radius
          do joff = -radius, radius
            idk = idx + ioff
            idl = idy + joff
            if (idk < 0) idk = -idk
            if (idl < 0) idl = -idl
            if (idk > width - 1) idk = width - 1 - ioff
            if (idl > height - 1) idl = height - 1 - joff
            id_w = idl * width + idk
            neighbor = input(id_w)
            range_term = -((center - neighbor) * (center - neighbor)) / (2.0_real32 * variance_i)
            spatial_term = -real((idk - idx) * (idk - idx) + (idl - idy) * (idl - idy), real32) / &
                (2.0_real32 * variance_spatial)
            weight = a_square * exp(spatial_term + range_term)
            normalization = normalization + weight
            result = result + neighbor * weight
          end do
        end do
        output(id) = result / normalization
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine bilateral_filter

  subroutine bilateral_reference(radius, input, output, width, height, a_square, variance_i, variance_spatial)
    integer, intent(in) :: radius, width, height
    real(real32), intent(in) :: input(0:), a_square, variance_i, variance_spatial
    real(real32), intent(out) :: output(0:)
    integer :: idx, idy, ioff, joff, id, idk, idl, id_w
    real(real32) :: center, neighbor, range_term, spatial_term, weight
    real(real32) :: result, normalization

    !$omp parallel do collapse(2) private(ioff, joff, id, idk, idl, id_w, center, neighbor, range_term, spatial_term, weight, result, normalization)
    do idx = 0, width - 1
      do idy = 0, height - 1
        id = idy * width + idx
        center = input(id)
        result = 0.0_real32
        normalization = 0.0_real32
        do ioff = -radius, radius
          do joff = -radius, radius
            idk = idx + ioff
            idl = idy + joff
            if (idk < 0) idk = -idk
            if (idl < 0) idl = -idl
            if (idk > width - 1) idk = width - 1 - ioff
            if (idl > height - 1) idl = height - 1 - joff
            id_w = idl * width + idk
            neighbor = input(id_w)
            range_term = -((center - neighbor) * (center - neighbor)) / (2.0_real32 * variance_i)
            spatial_term = -real((idk - idx) * (idk - idx) + (idl - idy) * (idl - idy), real32) / &
                (2.0_real32 * variance_spatial)
            weight = a_square * exp(spatial_term + range_term)
            normalization = normalization + weight
            result = result + neighbor * weight
          end do
        end do
        output(id) = result / normalization
      end do
    end do
    !$omp end parallel do
  end subroutine bilateral_reference

  logical function check_values(device_values, reference_values, size)
    real(real32), intent(in) :: device_values(0:), reference_values(0:)
    integer, intent(in) :: size
    integer :: idx

    check_values = .true.
    do idx = 0, size - 1
      if (abs(device_values(idx) - reference_values(idx)) > 1.0e-3_real32) then
        check_values = .false.
        exit
      end if
    end do
  end function check_values

end program bilateral_main
