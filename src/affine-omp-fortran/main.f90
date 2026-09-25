! SPDX-License-Identifier: CC0-1.0
program affine_main
  use, intrinsic :: iso_fortran_env, only : int16, int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: x_size = 512
  integer, parameter :: y_size = 512
  integer, parameter :: image_size = x_size * y_size
  real(real32), parameter :: pi = 3.14159265359_real32
  integer(int16), parameter :: white = 1_int16
  integer(int16), allocatable :: input_image(:), output_image(:), output_image_ref(:)
  integer :: iterations
  integer :: x, y, max_error, bytes_count
  real(real64) :: start_time, elapsed_s
  character(len=512) :: input_filename, output_filename, arg

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <input image> <output image> <iterations>")')
    stop 1
  end if

  call get_command_argument(1, input_filename)
  call get_command_argument(2, output_filename)
  call get_command_argument(3, arg)
  read(arg, *) iterations

  allocate(input_image(0:image_size-1), output_image(0:image_size-1), output_image_ref(0:image_size-1))

  write(*,'("Reading input image...")')
  write(*,*)
  write(*,'("   Reading RAW Image")')
  call read_raw_image(trim(input_filename), input_image)
  bytes_count = image_size * 2
  write(*,'("   Bytes read = ",I0)') bytes_count
  write(*,*)

  !$omp target data map(to: input_image) map(from: output_image)
  start_time = omp_get_wtime()
  do x = 1, iterations
    call affine_kernel(input_image, output_image)
  end do
  elapsed_s = (omp_get_wtime() - start_time) / real(iterations, real64)
  !$omp end target data
  write(*,'("   Average kernel execution time ",F11.9," (s)")') elapsed_s

  call affine_reference(input_image, output_image_ref)
  max_error = 0
  do y = 0, y_size - 1
    do x = 0, x_size - 1
      max_error = max(max_error, abs(unpack_u16(output_image(y * x_size + x)) - &
          unpack_u16(output_image_ref(y * x_size + x))))
    end do
  end do
  write(*,'("   Max output error is ",I0)') max_error
  write(*,*)

  write(*,'("   Writing RAW Image")')
  call write_raw_image(trim(output_filename), output_image)
  write(*,'("   Bytes written = ",I0)') bytes_count
  write(*,*)

contains

  subroutine read_raw_image(path, image)
    character(len=*), intent(in) :: path
    integer(int16), intent(out) :: image(0:)
    integer :: unit
    logical :: exists

    inquire(file=path, exist=exists)
    if (.not. exists) then
      write(*,'("Error: Unable to open input image file ",A,"!")') trim(path)
      stop 1
    end if

    open(newunit=unit, file=path, access="stream", form="unformatted", status="old", action="read")
    read(unit) image
    close(unit)
  end subroutine read_raw_image

  subroutine write_raw_image(path, image)
    character(len=*), intent(in) :: path
    integer(int16), intent(in) :: image(0:)
    integer :: unit

    open(newunit=unit, file=path, access="stream", form="unformatted", status="replace", action="write")
    write(unit) image
    close(unit)
  end subroutine write_raw_image

  integer(int16) function pack_u16(value)
    integer(int32), intent(in) :: value
    integer(int32) :: masked

    masked = iand(value, 65535_int32)
    pack_u16 = transfer(masked, pack_u16)
  end function pack_u16

  integer(int32) function unpack_u16(value)
    integer(int16), intent(in) :: value

    unpack_u16 = iand(int(value, int32), 65535_int32)
  end function unpack_u16

  subroutine affine_kernel(src, dst)
    integer(int16), intent(in) :: src(0:)
    integer(int16), intent(out) :: dst(0:)
    integer :: x, y

    !$omp target teams distribute parallel do collapse(2) thread_limit(256)
    do y = 0, y_size - 1
      do x = 0, x_size - 1
        dst(y * x_size + x) = affine_pixel(src, x, y)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine affine_kernel

  subroutine affine_reference(src, dst)
    integer(int16), intent(in) :: src(0:)
    integer(int16), intent(out) :: dst(0:)
    integer :: x, y

    do y = 0, y_size - 1
      do x = 0, x_size - 1
        dst(y * x_size + x) = affine_pixel(src, x, y)
      end do
    end do
  end subroutine affine_reference

  integer(int16) function affine_pixel(src, x, y)
    integer(int16), intent(in) :: src(0:)
    integer, intent(in) :: x, y
    real(real32), parameter :: lx_rot = 30.0_real32
    real(real32), parameter :: ly_rot = 0.0_real32
    real(real32), parameter :: lx_expan = 0.5_real32
    real(real32), parameter :: ly_expan = 0.5_real32
    real(real32) :: affine00, affine01, affine10, affine11
    real(real32) :: i_affine00, i_affine01, i_affine10, i_affine11
    real(real32) :: beta0, beta1, i_beta0, i_beta1, det
    real(real32) :: x_new, y_new, x_frac, y_frac, gray_new
    integer :: m, n

    affine00 = lx_expan * cos(lx_rot * pi / 180.0_real32)
    affine01 = ly_expan * sin(ly_rot * pi / 180.0_real32)
    affine10 = lx_expan * sin(lx_rot * pi / 180.0_real32)
    affine11 = ly_expan * cos(ly_rot * pi / 180.0_real32)
    beta0 = 0.0_real32
    beta1 = 0.0_real32

    det = affine00 * affine11 - affine01 * affine10
    if (det == 0.0_real32) then
      i_affine00 = 1.0_real32
      i_affine01 = 0.0_real32
      i_affine10 = 0.0_real32
      i_affine11 = 1.0_real32
      i_beta0 = -beta0
      i_beta1 = -beta1
    else
      i_affine00 = affine11 / det
      i_affine01 = -affine01 / det
      i_affine10 = -affine10 / det
      i_affine11 = affine00 / det
      i_beta0 = -i_affine00 * beta0 - i_affine01 * beta1
      i_beta1 = -i_affine10 * beta0 - i_affine11 * beta1
    end if

    x_new = i_beta0 + i_affine00 * (real(x, real32) - real(x_size, real32) / 2.0_real32) + &
        i_affine01 * (real(y, real32) - real(y_size, real32) / 2.0_real32) + real(x_size, real32) / 2.0_real32
    y_new = i_beta1 + i_affine10 * (real(x, real32) - real(x_size, real32) / 2.0_real32) + &
        i_affine11 * (real(y, real32) - real(y_size, real32) / 2.0_real32) + real(y_size, real32) / 2.0_real32

    m = floor(x_new)
    n = floor(y_new)
    x_frac = x_new - real(m, real32)
    y_frac = y_new - real(n, real32)

    if (m >= 0 .and. m + 1 < x_size .and. n >= 0 .and. n + 1 < y_size) then
      gray_new = (1.0_real32 - y_frac) * ((1.0_real32 - x_frac) * &
          real(unpack_u16(src(n * x_size + m)), real32) + &
          x_frac * real(unpack_u16(src(n * x_size + m + 1)), real32)) + &
          y_frac * ((1.0_real32 - x_frac) * &
          real(unpack_u16(src((n + 1) * x_size + m)), real32) + &
          x_frac * real(unpack_u16(src((n + 1) * x_size + m + 1)), real32))
      affine_pixel = pack_u16(int(gray_new, int32))
    else if (((m + 1 == x_size) .and. n >= 0 .and. n < y_size) .or. &
        ((n + 1 == y_size) .and. m >= 0 .and. m < x_size)) then
      affine_pixel = src(n * x_size + m)
    else
      affine_pixel = white
    end if
  end function affine_pixel

end program affine_main
