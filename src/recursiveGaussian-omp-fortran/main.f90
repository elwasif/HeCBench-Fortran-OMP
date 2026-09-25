! SPDX-License-Identifier: CC0-1.0
program recursive_gaussian_omp_fortran
  use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
  use omp_lib, only: omp_get_team_num, omp_get_thread_num, omp_get_wtime
  implicit none

  type :: gauss_parms
    real(real32) :: nsigma, alpha, ema, ema2
    real(real32) :: b1, b2, a0, a1, a2, a3, coefp, coefn
  end type gauss_parms

  integer, parameter :: max_image_width = 1920
  integer, parameter :: max_image_height = 1080
  real(real32), parameter :: sigma = 10.0_real32
  integer, parameter :: order = 0

  character(len=512) :: image_path, arg, program_name
  integer :: argc, ios, cycles, width, height
  integer(int32), allocatable :: input(:), tmp(:), output(:), golden(:)
  type(gauss_parms) :: gp
  logical :: status, match
  integer :: n_pixels, iter
  real(real64) :: elapsed, warmup_time

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, program_name)
    write(*, '(A,A,A)') 'Usage: ', trim(program_name), ' <path to image> <repeat>'
    stop 1
  end if

  call get_command_argument(1, image_path)
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) cycles
  if (ios /= 0) stop 1

  call load_ppm4ub(trim(image_path), input, width, height, status)

  write(*, '(A,I0,A,I0,A,I0)') 'Image Width = ', width, ', Height = ', height, ', bpp = ', 32
  write(*, *)

  if (width > max_image_width .or. height > max_image_height) then
    write(*, '(A)', advance='no') 'Error: Image Dimensions exceed the maximum values'
    status = .false.
  end if
  if (.not. status) then
    if (allocated(input)) deallocate(input)
    stop 1
  end if

  n_pixels = width * height
  allocate(tmp(n_pixels), output(n_pixels), golden(n_pixels))
  tmp = 0_int32
  output = 0_int32
  golden = 0_int32
  write(*, '(A)') 'Allocate Host Image Buffers...'

  call preprocess_gauss_parms(sigma, order, gp)

  !$omp target data map(alloc: input(1:n_pixels), tmp(1:n_pixels), output(1:n_pixels))
    warmup_time = gpu_gaussian_filter_rgba(input, tmp, output, width, height, gp)

    write(*, *)
    write(*, '(A,I0,A)') 'Running GPUGaussianFilterRGBA for ', cycles, ' cycles...'
    write(*, *)

    elapsed = 0.0_real64
    do iter = 1, cycles
      elapsed = elapsed + gpu_gaussian_filter_rgba(input, tmp, output, width, height, gp)
    end do
    write(*, '(A,F8.6,A)') 'Average execution time of kernels: ', elapsed / real(cycles, real64), ' (s)'
  !$omp end target data

  call host_recursive_gaussian_rgba(input, tmp, golden, width, height, gp)

  write(*, '(A)') 'Comparing GPU Result to CPU Result...'
  match = compare_uint_threshold(golden, output, n_pixels, 1.0_real32, 0.01_real32)
  write(*, *)
  if (match) then
    write(*, '(A)') 'GPU Result matches CPU Result within tolerance...'
  else
    write(*, '(A)') "GPU Result DOESN'T match CPU Result within tolerance..."
  end if

  deallocate(input, tmp, output, golden)

contains

  subroutine load_ppm4ub(path, rgba, width, height, ok)
    character(len=*), intent(in) :: path
    integer(int32), allocatable, intent(out) :: rgba(:)
    integer, intent(out) :: width, height
    logical, intent(out) :: ok
    integer(int8), allocatable :: bytes(:)
    integer :: unit, ios, file_size, pos, maxval, pixel, base
    character(len=64) :: magic, token

    ok = .false.
    width = 0
    height = 0
    inquire(file=path, size=file_size)
    if (file_size <= 0) return

    allocate(bytes(file_size))
    open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      deallocate(bytes)
      return
    end if
    read(unit, iostat=ios) bytes
    close(unit)
    if (ios /= 0) then
      deallocate(bytes)
      return
    end if

    pos = 1
    call next_token(bytes, pos, magic)
    call next_token(bytes, pos, token)
    read(token, *, iostat=ios) width
    if (ios /= 0) then
      deallocate(bytes)
      return
    end if
    call next_token(bytes, pos, token)
    read(token, *, iostat=ios) height
    if (ios /= 0) then
      deallocate(bytes)
      return
    end if
    call next_token(bytes, pos, token)
    read(token, *, iostat=ios) maxval
    if (ios /= 0) then
      deallocate(bytes)
      return
    end if

    if (trim(magic) /= 'P6' .or. width <= 0 .or. height <= 0 .or. maxval /= 255) then
      deallocate(bytes)
      return
    end if
    if (pos <= file_size .and. is_ws(byte_value(bytes(pos)))) pos = pos + 1
    if (pos + width * height * 3 - 1 > file_size) then
      deallocate(bytes)
      return
    end if

    allocate(rgba(width * height))
    do pixel = 1, width * height
      base = pos + (pixel - 1) * 3
      rgba(pixel) = ior(ior(byte_value(bytes(base)), ishft(byte_value(bytes(base + 1)), 8)), &
                        ishft(byte_value(bytes(base + 2)), 16))
    end do

    deallocate(bytes)
    ok = .true.
  end subroutine load_ppm4ub

  subroutine next_token(bytes, pos, token)
    integer(int8), intent(in) :: bytes(:)
    integer, intent(inout) :: pos
    character(len=*), intent(out) :: token
    integer :: len_token, c

    token = ''
    call skip_ws_and_comments(bytes, pos)
    len_token = 0
    do while (pos <= size(bytes))
      c = byte_value(bytes(pos))
      if (is_ws(c) .or. c == iachar('#')) exit
      len_token = len_token + 1
      if (len_token <= len(token)) token(len_token:len_token) = achar(c)
      pos = pos + 1
    end do
  end subroutine next_token

  subroutine skip_ws_and_comments(bytes, pos)
    integer(int8), intent(in) :: bytes(:)
    integer, intent(inout) :: pos
    integer :: c

    do while (pos <= size(bytes))
      c = byte_value(bytes(pos))
      if (is_ws(c)) then
        pos = pos + 1
      else if (c == iachar('#')) then
        do while (pos <= size(bytes) .and. byte_value(bytes(pos)) /= 10)
          pos = pos + 1
        end do
      else
        exit
      end if
    end do
  end subroutine skip_ws_and_comments

  logical function is_ws(c)
    integer, intent(in) :: c
    is_ws = (c == 9 .or. c == 10 .or. c == 13 .or. c == 32)
  end function is_ws

  integer(int32) function byte_value(b)
    integer(int8), intent(in) :: b
    integer :: tmp
    tmp = int(b)
    if (tmp < 0) tmp = tmp + 256
    byte_value = int(tmp, int32)
  end function byte_value

  subroutine preprocess_gauss_parms(f_sigma, i_order, gp)
    real(real32), intent(in) :: f_sigma
    integer, intent(in) :: i_order
    type(gauss_parms), intent(out) :: gp
    real(real32) :: k, ea, kn

    gp%nsigma = f_sigma
    gp%alpha = 1.695_real32 / gp%nsigma
    gp%ema = exp(-gp%alpha)
    gp%ema2 = exp(-2.0_real32 * gp%alpha)
    gp%b1 = -2.0_real32 * gp%ema
    gp%b2 = gp%ema2
    gp%a0 = 0.0_real32
    gp%a1 = 0.0_real32
    gp%a2 = 0.0_real32
    gp%a3 = 0.0_real32
    gp%coefp = 0.0_real32
    gp%coefn = 0.0_real32

    select case (i_order)
    case (0)
      k = (1.0_real32 - gp%ema) * (1.0_real32 - gp%ema) / &
          (1.0_real32 + 2.0_real32 * gp%alpha * gp%ema - gp%ema2)
      gp%a0 = k
      gp%a1 = k * (gp%alpha - 1.0_real32) * gp%ema
      gp%a2 = k * (gp%alpha + 1.0_real32) * gp%ema
      gp%a3 = -k * gp%ema2
    case (1)
      gp%a0 = (1.0_real32 - gp%ema) * (1.0_real32 - gp%ema)
      gp%a2 = -gp%a0
    case (2)
      ea = exp(-gp%alpha)
      k = -(gp%ema2 - 1.0_real32) / (2.0_real32 * gp%alpha * gp%ema)
      kn = -2.0_real32 * (-1.0_real32 + 3.0_real32 * ea - 3.0_real32 * ea * ea + ea * ea * ea)
      kn = kn / (3.0_real32 * ea + 1.0_real32 + 3.0_real32 * ea * ea + ea * ea * ea)
      gp%a0 = kn
      gp%a1 = -kn * (1.0_real32 + k * gp%alpha) * gp%ema
      gp%a2 = kn * (1.0_real32 - k * gp%alpha) * gp%ema
      gp%a3 = -kn * gp%ema2
    end select
    gp%coefp = (gp%a0 + gp%a1) / (1.0_real32 + gp%b1 + gp%b2)
    gp%coefn = (gp%a2 + gp%a3) / (1.0_real32 + gp%b1 + gp%b2)
  end subroutine preprocess_gauss_parms

  real(real64) function gpu_gaussian_filter_rgba(input, tmp, output, width, height, gp)
    integer(int32), intent(inout) :: input(:), tmp(:), output(:)
    integer, intent(in) :: width, height
    type(gauss_parms), intent(in) :: gp
    real(real64) :: start_time, end_time

    !$omp target update to(input(1:width * height))
    start_time = omp_get_wtime()
    call recursive_rgba_device(input, tmp, width, height, gp%a0, gp%a1, gp%a2, gp%a3, gp%b1, gp%b2, gp%coefp, gp%coefn)
    call transpose_device(tmp, output, width, height)
    call recursive_rgba_device(output, tmp, height, width, gp%a0, gp%a1, gp%a2, gp%a3, gp%b1, gp%b2, gp%coefp, gp%coefn)
    call transpose_device(tmp, output, height, width)
    end_time = omp_get_wtime()
    !$omp target update from(output(1:width * height))
    gpu_gaussian_filter_rgba = end_time - start_time
  end function gpu_gaussian_filter_rgba

  subroutine recursive_rgba_device(data_in, data_out, width, height, a0, a1, a2, a3, b1, b2, coefp, coefn)
    integer(int32), intent(in) :: data_in(:)
    integer(int32), intent(inout) :: data_out(:)
    integer, intent(in) :: width, height
    real(real32), intent(in) :: a0, a1, a2, a3, b1, b2, coefp, coefn
    integer :: x
    integer :: szGaussLocalWork, szGaussGlobalWork, szTeams

    szGaussLocalWork = 256
    szGaussGlobalWork = ((width + szGaussLocalWork - 1) / szGaussLocalWork) * szGaussLocalWork
    szTeams = szGaussGlobalWork / szGaussLocalWork

    !$omp target teams distribute parallel do num_teams(szTeams) thread_limit(szGaussLocalWork)
    do x = 1, width
      block
        integer :: y, idx, c
        real(real32) :: xp0, xp1, xp2, xp3, yp0, yp1, yp2, yp3
        real(real32) :: yb0, yb1, yb2, yb3, xc0, xc1, xc2, xc3
        real(real32) :: yc0, yc1, yc2, yc3
        real(real32) :: xn0, xn1, xn2, xn3, xa0, xa1, xa2, xa3
        real(real32) :: yn0, yn1, yn2, yn3, ya0, ya1, ya2, ya3
        real(real32) :: out0, out1, out2, out3

        c = data_in(x)
        xp0 = real(iand(c, int(z'000000ff', int32)), real32)
        xp1 = real(iand(ishft(c, -8), int(z'000000ff', int32)), real32)
        xp2 = real(iand(ishft(c, -16), int(z'000000ff', int32)), real32)
        xp3 = real(iand(ishft(c, -24), int(z'000000ff', int32)), real32)
        yb0 = xp0 * coefp
        yb1 = xp1 * coefp
        yb2 = xp2 * coefp
        yb3 = xp3 * coefp
        yp0 = yb0
        yp1 = yb1
        yp2 = yb2
        yp3 = yb3

        do y = 1, height
          idx = (y - 1) * width + x
          c = data_in(idx)
          xc0 = real(iand(c, int(z'000000ff', int32)), real32)
          xc1 = real(iand(ishft(c, -8), int(z'000000ff', int32)), real32)
          xc2 = real(iand(ishft(c, -16), int(z'000000ff', int32)), real32)
          xc3 = real(iand(ishft(c, -24), int(z'000000ff', int32)), real32)
          yc0 = xc0 * a0 + xp0 * a1 - yp0 * b1 - yb0 * b2
          yc1 = xc1 * a0 + xp1 * a1 - yp1 * b1 - yb1 * b2
          yc2 = xc2 * a0 + xp2 * a1 - yp2 * b1 - yb2 * b2
          yc3 = xc3 * a0 + xp3 * a1 - yp3 * b1 - yb3 * b2
          data_out(idx) = pack_rgba_device(yc0, yc1, yc2, yc3)
          xp0 = xc0
          xp1 = xc1
          xp2 = xc2
          xp3 = xc3
          yb0 = yp0
          yb1 = yp1
          yb2 = yp2
          yb3 = yp3
          yp0 = yc0
          yp1 = yc1
          yp2 = yc2
          yp3 = yc3
        end do

        idx = (height - 1) * width + x
        c = data_in(idx)
        xn0 = real(iand(c, int(z'000000ff', int32)), real32)
        xn1 = real(iand(ishft(c, -8), int(z'000000ff', int32)), real32)
        xn2 = real(iand(ishft(c, -16), int(z'000000ff', int32)), real32)
        xn3 = real(iand(ishft(c, -24), int(z'000000ff', int32)), real32)
        xa0 = xn0
        xa1 = xn1
        xa2 = xn2
        xa3 = xn3
        yn0 = xn0 * coefn
        yn1 = xn1 * coefn
        yn2 = xn2 * coefn
        yn3 = xn3 * coefn
        ya0 = yn0
        ya1 = yn1
        ya2 = yn2
        ya3 = yn3

        do y = height, 1, -1
          idx = (y - 1) * width + x
          c = data_in(idx)
          xc0 = real(iand(c, int(z'000000ff', int32)), real32)
          xc1 = real(iand(ishft(c, -8), int(z'000000ff', int32)), real32)
          xc2 = real(iand(ishft(c, -16), int(z'000000ff', int32)), real32)
          xc3 = real(iand(ishft(c, -24), int(z'000000ff', int32)), real32)
          yc0 = xn0 * a2 + xa0 * a3 - yn0 * b1 - ya0 * b2
          yc1 = xn1 * a2 + xa1 * a3 - yn1 * b1 - ya1 * b2
          yc2 = xn2 * a2 + xa2 * a3 - yn2 * b1 - ya2 * b2
          yc3 = xn3 * a2 + xa3 * a3 - yn3 * b1 - ya3 * b2
          xa0 = xn0
          xa1 = xn1
          xa2 = xn2
          xa3 = xn3
          xn0 = xc0
          xn1 = xc1
          xn2 = xc2
          xn3 = xc3
          ya0 = yn0
          ya1 = yn1
          ya2 = yn2
          ya3 = yn3
          yn0 = yc0
          yn1 = yc1
          yn2 = yc2
          yn3 = yc3
          c = data_out(idx)
          out0 = real(iand(c, int(z'000000ff', int32)), real32) + yc0
          out1 = real(iand(ishft(c, -8), int(z'000000ff', int32)), real32) + yc1
          out2 = real(iand(ishft(c, -16), int(z'000000ff', int32)), real32) + yc2
          out3 = real(iand(ishft(c, -24), int(z'000000ff', int32)), real32) + yc3
          data_out(idx) = pack_rgba_device(out0, out1, out2, out3)
        end do
      end block
    end do
    !$omp end target teams distribute parallel do
  end subroutine recursive_rgba_device

  subroutine transpose_device(data_in, data_out, width, height)
    integer(int32), intent(in) :: data_in(:)
    integer(int32), intent(inout) :: data_out(:)
    integer, intent(in) :: width, height
    integer :: szTransposeGlobalWork(2), szTransposeLocalWork(2)
    integer :: numTeamsX, numTeamsY, numTeams

    szTransposeLocalWork = [16, 16]
    szTransposeGlobalWork(1) = ((width + szTransposeLocalWork(1) - 1) / szTransposeLocalWork(1)) * szTransposeLocalWork(1)
    szTransposeGlobalWork(2) = ((height + szTransposeLocalWork(2) - 1) / szTransposeLocalWork(2)) * szTransposeLocalWork(2)
    numTeamsX = szTransposeGlobalWork(1) / szTransposeLocalWork(1)
    numTeamsY = szTransposeGlobalWork(2) / szTransposeLocalWork(2)
    numTeams = numTeamsX * numTeamsY

    !$omp target teams num_teams(numTeams) thread_limit(256)
    block
      integer(int32) :: uiLocalBuff(16 * 17)

      !$omp parallel
      block
        integer :: lidX, lidY, tidX, tidY, xIndex, yIndex

        lidX = mod(omp_get_thread_num(), 16)
        lidY = omp_get_thread_num() / 16
        tidX = mod(omp_get_team_num(), numTeamsX)
        tidY = omp_get_team_num() / numTeamsX

        xIndex = tidX * 16 + lidX
        yIndex = tidY * 16 + lidY

        if (xIndex < width .and. yIndex < height) then
          uiLocalBuff(lidY * 17 + lidX + 1) = data_in(yIndex * width + xIndex + 1)
        end if

        !$omp barrier

        xIndex = tidY * 16 + lidX
        yIndex = tidX * 16 + lidY

        if (xIndex < height .and. yIndex < width) then
          data_out(yIndex * height + xIndex + 1) = uiLocalBuff(lidX * 17 + lidY + 1)
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine transpose_device

  integer(int32) function pack_rgba_device(r, g, b, a)
    real(real32), intent(in) :: r, g, b, a
    pack_rgba_device = ior(ior(iand(int(r, int32), int(z'000000ff', int32)), &
                              ishft(iand(int(g, int32), int(z'000000ff', int32)), 8)), &
                          ior(ishft(iand(int(b, int32), int(z'000000ff', int32)), 16), &
                              ishft(iand(int(a, int32), int(z'000000ff', int32)), 24)))
  end function pack_rgba_device

  subroutine host_recursive_gaussian_rgba(input, tmp, output, width, height, gp)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: tmp(:), output(:)
    integer, intent(in) :: width, height
    type(gauss_parms), intent(in) :: gp

    call recursive_rgba_host(input, tmp, width, height, gp%a0, gp%a1, gp%a2, gp%a3, gp%b1, gp%b2, gp%coefp, gp%coefn)
    call transpose_host(tmp, output, width, height)
    call recursive_rgba_host(output, tmp, height, width, gp%a0, gp%a1, gp%a2, gp%a3, gp%b1, gp%b2, gp%coefp, gp%coefn)
    call transpose_host(tmp, output, height, width)
  end subroutine host_recursive_gaussian_rgba

  subroutine recursive_rgba_host(data_in, data_out, width, height, a0, a1, a2, a3, b1, b2, coefp, coefn)
    integer(int32), intent(in) :: data_in(:)
    integer(int32), intent(inout) :: data_out(:)
    integer, intent(in) :: width, height
    real(real32), intent(in) :: a0, a1, a2, a3, b1, b2, coefp, coefn
    integer :: x, y, idx, c
    real(real32) :: xp(4), yp(4), yb(4), xc(4), yc(4)
    real(real32) :: xn(4), xa(4), yn(4), ya(4), outv(4)

    do x = 1, width
      call unpack_rgba(data_in(x), xp)
      yb = xp * coefp
      yp = yb
      do y = 1, height
        idx = (y - 1) * width + x
        call unpack_rgba(data_in(idx), xc)
        yc = xc * a0 + xp * a1 - yp * b1 - yb * b2
        data_out(idx) = pack_rgba_host(yc)
        xp = xc
        yb = yp
        yp = yc
      end do

      idx = (height - 1) * width + x
      call unpack_rgba(data_in(idx), xn)
      xa = xn
      yn = xn * coefn
      ya = yn
      do y = height, 1, -1
        idx = (y - 1) * width + x
        call unpack_rgba(data_in(idx), xc)
        yc = xn * a2 + xa * a3 - yn * b1 - ya * b2
        xa = xn
        xn = xc
        ya = yn
        yn = yc
        c = data_out(idx)
        call unpack_rgba(c, outv)
        outv = outv + yc
        data_out(idx) = pack_rgba_host(outv)
      end do
    end do
  end subroutine recursive_rgba_host

  subroutine transpose_host(data_in, data_out, width, height)
    integer(int32), intent(in) :: data_in(:)
    integer(int32), intent(inout) :: data_out(:)
    integer, intent(in) :: width, height
    integer :: x, y

    do y = 1, height
      do x = 1, width
        data_out((x - 1) * height + y) = data_in((y - 1) * width + x)
      end do
    end do
  end subroutine transpose_host

  subroutine unpack_rgba(pixel, rgba)
    integer(int32), intent(in) :: pixel
    real(real32), intent(out) :: rgba(4)
    rgba(1) = real(iand(pixel, int(z'000000ff', int32)), real32)
    rgba(2) = real(iand(ishft(pixel, -8), int(z'000000ff', int32)), real32)
    rgba(3) = real(iand(ishft(pixel, -16), int(z'000000ff', int32)), real32)
    rgba(4) = real(iand(ishft(pixel, -24), int(z'000000ff', int32)), real32)
  end subroutine unpack_rgba

  integer(int32) function pack_rgba_host(rgba)
    real(real32), intent(inout) :: rgba(4)
    integer :: i

    do i = 1, 4
      if (rgba(i) < 0.0_real32) rgba(i) = 0.0_real32
    end do
    pack_rgba_host = ior(ior(iand(int(rgba(1), int32), int(z'000000ff', int32)), &
                            ishft(iand(int(rgba(2), int32), int(z'000000ff', int32)), 8)), &
                        ior(ishft(iand(int(rgba(3), int32), int(z'000000ff', int32)), 16), &
                            ishft(iand(int(rgba(4), int32), int(z'000000ff', int32)), 24)))
  end function pack_rgba_host

  logical function compare_uint_threshold(reference, data, n, epsilon, threshold)
    integer(int32), intent(in) :: reference(:), data(:)
    integer, intent(in) :: n
    real(real32), intent(in) :: epsilon, threshold
    integer :: i, error_count
    real(real32) :: max_error, diff

    max_error = max(epsilon, 1.0e-3_real32)
    error_count = 0
    do i = 1, n
      diff = abs(real(reference(i), real32) - real(data(i), real32))
      if (diff >= max_error) error_count = error_count + 1
    end do
    if (threshold == 0.0_real32) then
      compare_uint_threshold = (error_count == 0)
    else
      compare_uint_threshold = (real(n, real32) * threshold > real(error_count, real32))
    end if
  end function compare_uint_threshold

end program recursive_gaussian_omp_fortran
