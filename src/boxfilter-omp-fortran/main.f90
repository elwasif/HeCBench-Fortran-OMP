! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: radius = 10
  real(real32), parameter :: scale = 1.0_real32 / real(2 * radius + 1, real32)

  character(len=512) :: image_path
  integer :: repeat, width, height, n_pixels, i, error
  integer(int32), allocatable :: input(:), tmp(:), dev_output(:), host_output(:)

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage ./main <PPM image> <repeat>'
    stop 1
  end if

  call get_command_argument(1, image_path)
  repeat = read_arg(2)

  call load_ppm(trim(image_path), input, width, height)
  n_pixels = width * height
  allocate(tmp(n_pixels), dev_output(n_pixels), host_output(n_pixels))
  tmp = 0_int32
  dev_output = 0_int32
  host_output = 0_int32

  print '(A,I0,A,I0,A,I0,A,I0)', 'Image Width = ', width, ', Height = ', height, &
      ', bpp = ', 32, ', Mask Radius = ', radius
  print '(A)', 'Using Local Memory for Row Processing'
  print '(A)'

  !$omp target data map(to: input(1:n_pixels)) map(tofrom: tmp(1:n_pixels), dev_output(1:n_pixels))
  print '(A)', 'Warmup..'
  call box_filter_device(input, tmp, dev_output, width, height, repeat)
  print '(A)'
  print '(A,I0,A)', 'Running BoxFilterGPU for ', repeat, ' cycles...'
  print '(A)'
  call box_filter_device(input, tmp, dev_output, width, height, repeat)
  !$omp end target data

  call box_filter_host(input, tmp, host_output, width, height)

  error = 0
  do i = radius * width + 1, n_pixels - radius * width
    if (dev_output(i) /= host_output(i)) then
      write(*, '(I0,1X,Z8.8,1X,Z8.8)') i - 1, dev_output(i), host_output(i)
      error = 1
      exit
    end if
  end do
  if (error == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(input, tmp, dev_output, host_output)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine load_ppm(path, image, width, height)
    character(len=*), intent(in) :: path
    integer(int32), allocatable, intent(out) :: image(:)
    integer, intent(out) :: width, height
    integer :: unit, file_size, pos, maxval, pixel, base
    character(len=1), allocatable :: bytes(:)
    character(len=64) :: token

    inquire(file=path, size=file_size)
    if (file_size <= 0) error stop 'invalid PPM input'
    allocate(bytes(file_size))
    open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read')
    read(unit) bytes
    close(unit)

    pos = 1
    call next_token(bytes, file_size, pos, token)
    if (trim(token) /= 'P6') error stop 'expected binary P6 PPM'
    call next_token(bytes, file_size, pos, token)
    read(token, *) width
    call next_token(bytes, file_size, pos, token)
    read(token, *) height
    call next_token(bytes, file_size, pos, token)
    read(token, *) maxval
    if (maxval /= 255) error stop 'expected 8-bit PPM'

    allocate(image(width * height))
    do pixel = 1, width * height
      base = pos + (pixel - 1) * 3
      image(pixel) = pack_rgba(byte_value(bytes(base)), byte_value(bytes(base + 1)), byte_value(bytes(base + 2)), 0)
    end do
    deallocate(bytes)
  end subroutine load_ppm

  subroutine next_token(bytes, file_size, pos, token)
    character(len=1), intent(in) :: bytes(:)
    integer, intent(in) :: file_size
    integer, intent(inout) :: pos
    character(len=*), intent(out) :: token
    integer :: out_pos, ch
    token = ''
    do while (pos <= file_size)
      ch = iachar(bytes(pos))
      if (ch == iachar('#')) then
        do while (pos <= file_size .and. iachar(bytes(pos)) /= 10)
          pos = pos + 1
        end do
      else if (.not. is_space(ch)) then
        exit
      end if
      pos = pos + 1
    end do
    out_pos = 1
    do while (pos <= file_size)
      ch = iachar(bytes(pos))
      if (is_space(ch) .or. ch == iachar('#')) exit
      if (out_pos <= len(token)) token(out_pos:out_pos) = bytes(pos)
      out_pos = out_pos + 1
      pos = pos + 1
    end do
    do while (pos <= file_size .and. is_space(iachar(bytes(pos))))
      pos = pos + 1
    end do
  end subroutine next_token

  logical function is_space(ch)
    integer, intent(in) :: ch
    is_space = ch == 9 .or. ch == 10 .or. ch == 13 .or. ch == 32
  end function is_space

  integer function byte_value(ch)
    character(len=1), intent(in) :: ch
    byte_value = iachar(ch)
  end function byte_value

  integer(int32) function pack_rgba(r, g, b, a)
    integer, intent(in) :: r, g, b, a
    pack_rgba = int(ior(ior(ior(iand(r, 255), ishft(iand(g, 255), 8)), &
        ishft(iand(b, 255), 16)), ishft(iand(a, 255), 24)), int32)
  end function pack_rgba

  subroutine unpack_rgba(pixel, r, g, b, a)
    integer(int32), intent(in) :: pixel
    real(real32), intent(out) :: r, g, b, a
    r = real(iand(pixel, int(z'000000ff', int32)), real32)
    g = real(iand(ishft(pixel, -8), int(z'000000ff', int32)), real32)
    b = real(iand(ishft(pixel, -16), int(z'000000ff', int32)), real32)
    a = real(iand(ishft(pixel, -24), int(z'000000ff', int32)), real32)
  end subroutine unpack_rgba

  integer(int32) function pack_scaled(r, g, b, a)
    real(real32), intent(in) :: r, g, b, a
    pack_scaled = pack_rgba(int(r * scale), int(g * scale), int(b * scale), int(a * scale))
  end function pack_scaled

  subroutine box_filter_device(input, tmp, output, width, height, repeat)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: tmp(:), output(:)
    integer, intent(in) :: width, height, repeat
    integer, parameter :: szMaxWorkgroupSize = 256
    integer :: cycle, globalPosX
    integer :: iRadiusAligned, uiNumOutputPix, uiBlockWidth, numTeams, blockSize
    integer :: lid, gidx, gidy, globalPosY, iGlobalOffset, iOffsetX, iLimit
    integer :: scratchOffset, y, uiInputOffset, uiOutputOffset
    real(real32), allocatable :: uc4LocalDataR(:), uc4LocalDataG(:), uc4LocalDataB(:), uc4LocalDataA(:)
    real(real32) :: f4SumR, f4SumG, f4SumB, f4SumA
    real(real32) :: r, g, b, a
    real(real32) :: topR, topG, topB, topA, botR, botG, botB, botA
    real(real64) :: start_time, end_time, avg_us

    iRadiusAligned = ((radius + 15) / 16) * 16
    uiNumOutputPix = 64
    if (szMaxWorkgroupSize < (iRadiusAligned + uiNumOutputPix + radius)) then
      uiNumOutputPix = szMaxWorkgroupSize - iRadiusAligned - radius
    end if
    uiBlockWidth = (width + uiNumOutputPix - 1) / uiNumOutputPix
    numTeams = height * uiBlockWidth
    blockSize = iRadiusAligned + uiNumOutputPix + radius
    allocate(uc4LocalDataR(numTeams * 90), uc4LocalDataG(numTeams * 90), &
        uc4LocalDataB(numTeams * 90), uc4LocalDataA(numTeams * 90))

    start_time = omp_get_wtime()
    do cycle = 1, repeat
      !$omp target teams num_teams(numTeams) thread_limit(blockSize) &
      !$omp& map(alloc: uc4LocalDataR(1:numTeams * 90), uc4LocalDataG(1:numTeams * 90), &
      !$omp& uc4LocalDataB(1:numTeams * 90), uc4LocalDataA(1:numTeams * 90))
        !$omp parallel private(lid, gidx, gidy, globalPosX, globalPosY, iGlobalOffset, scratchOffset, &
        !$omp& iOffsetX, iLimit, f4SumR, f4SumG, f4SumB, f4SumA, r, g, b, a)
        lid = omp_get_thread_num()
        gidx = mod(omp_get_team_num(), uiBlockWidth)
        gidy = omp_get_team_num() / uiBlockWidth

        globalPosX = gidx * uiNumOutputPix + lid - iRadiusAligned
        globalPosY = gidy
        iGlobalOffset = globalPosY * width + globalPosX + 1
        scratchOffset = omp_get_team_num() * 90 + lid + 1

        if (globalPosX >= 0 .and. globalPosX < width) then
          call unpack_rgba(input(iGlobalOffset), r, g, b, a)
          uc4LocalDataR(scratchOffset) = r
          uc4LocalDataG(scratchOffset) = g
          uc4LocalDataB(scratchOffset) = b
          uc4LocalDataA(scratchOffset) = a
        else
          uc4LocalDataR(scratchOffset) = 0.0_real32
          uc4LocalDataG(scratchOffset) = 0.0_real32
          uc4LocalDataB(scratchOffset) = 0.0_real32
          uc4LocalDataA(scratchOffset) = 0.0_real32
        end if

        !$omp barrier

        if (globalPosX >= 0 .and. globalPosX < width .and. lid >= iRadiusAligned .and. &
            lid < iRadiusAligned + uiNumOutputPix) then
          f4SumR = 0.0_real32
          f4SumG = 0.0_real32
          f4SumB = 0.0_real32
          f4SumA = 0.0_real32
          iOffsetX = lid - radius
          iLimit = iOffsetX + (2 * radius) + 1
          do while (iOffsetX < iLimit)
            scratchOffset = omp_get_team_num() * 90 + iOffsetX + 1
            f4SumR = f4SumR + uc4LocalDataR(scratchOffset)
            f4SumG = f4SumG + uc4LocalDataG(scratchOffset)
            f4SumB = f4SumB + uc4LocalDataB(scratchOffset)
            f4SumA = f4SumA + uc4LocalDataA(scratchOffset)
            iOffsetX = iOffsetX + 1
          end do
          tmp(iGlobalOffset) = pack_scaled(f4SumR, f4SumG, f4SumB, f4SumA)
        end if
        !$omp end parallel
      !$omp end target teams

      !$omp target teams distribute parallel do thread_limit(64) private(y, uiInputOffset, uiOutputOffset, &
      !$omp& f4SumR, f4SumG, f4SumB, f4SumA, topR, topG, topB, topA, botR, botG, botB, botA, r, g, b, a)
      do globalPosX = 0, width - 1
        call unpack_rgba(tmp(globalPosX + 1), topR, topG, topB, topA)
        call unpack_rgba(tmp((height - 1) * width + globalPosX + 1), botR, botG, botB, botA)

        f4SumR = topR * real(radius, real32)
        f4SumG = topG * real(radius, real32)
        f4SumB = topB * real(radius, real32)
        f4SumA = topA * real(radius, real32)
        do y = 0, radius
          uiInputOffset = y * width + globalPosX + 1
          call unpack_rgba(tmp(uiInputOffset), r, g, b, a)
          f4SumR = f4SumR + r
          f4SumG = f4SumG + g
          f4SumB = f4SumB + b
          f4SumA = f4SumA + a
        end do
        output(globalPosX + 1) = pack_scaled(f4SumR, f4SumG, f4SumB, f4SumA)

        do y = 1, radius
          uiInputOffset = (y + radius) * width + globalPosX + 1
          call unpack_rgba(tmp(uiInputOffset), r, g, b, a)
          f4SumR = f4SumR + r - topR
          f4SumG = f4SumG + g - topG
          f4SumB = f4SumB + b - topB
          f4SumA = f4SumA + a - topA
          uiOutputOffset = y * width + globalPosX + 1
          output(uiOutputOffset) = pack_scaled(f4SumR, f4SumG, f4SumB, f4SumA)
        end do

        do y = radius + 1, height - radius - 1
          uiInputOffset = (y + radius) * width + globalPosX + 1
          call unpack_rgba(tmp(uiInputOffset), r, g, b, a)
          f4SumR = f4SumR + r
          f4SumG = f4SumG + g
          f4SumB = f4SumB + b
          f4SumA = f4SumA + a
          uiInputOffset = ((y - radius) * width) + globalPosX + 1 - width
          call unpack_rgba(tmp(uiInputOffset), r, g, b, a)
          f4SumR = f4SumR - r
          f4SumG = f4SumG - g
          f4SumB = f4SumB - b
          f4SumA = f4SumA - a
          uiOutputOffset = y * width + globalPosX + 1
          output(uiOutputOffset) = pack_scaled(f4SumR, f4SumG, f4SumB, f4SumA)
        end do

        do y = height - radius, height - 1
          f4SumR = f4SumR + botR
          f4SumG = f4SumG + botG
          f4SumB = f4SumB + botB
          f4SumA = f4SumA + botA
          uiInputOffset = ((y - radius) * width) + globalPosX + 1 - width
          call unpack_rgba(tmp(uiInputOffset), r, g, b, a)
          f4SumR = f4SumR - r
          f4SumG = f4SumG - g
          f4SumB = f4SumB - b
          f4SumA = f4SumA - a
          uiOutputOffset = y * width + globalPosX + 1
          output(uiOutputOffset) = pack_scaled(f4SumR, f4SumG, f4SumB, f4SumA)
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    avg_us = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
    write(*, '(A,F0.6,A)') 'Average kernel execution time ', avg_us, ' (us)'
    deallocate(uc4LocalDataR, uc4LocalDataG, uc4LocalDataB, uc4LocalDataA)
  end subroutine box_filter_device

  subroutine box_filter_host(input, tmp, output, width, height)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: tmp(:), output(:)
    integer, intent(in) :: width, height
    integer :: idx, x, y, offset, yy
    real(real32) :: rs, gs, bs, as, r, g, b, a
    do y = 1, height
      do x = 1, width
        rs = 0.0_real32
        gs = 0.0_real32
        bs = 0.0_real32
        as = 0.0_real32
        do offset = -radius, radius
          if (x + offset >= 1 .and. x + offset <= width) then
            call unpack_rgba(input((y - 1) * width + x + offset), r, g, b, a)
            rs = rs + r
            gs = gs + g
            bs = bs + b
            as = as + a
          end if
        end do
        idx = (y - 1) * width + x
        tmp(idx) = pack_scaled(rs, gs, bs, as)
      end do
    end do
    do y = 1, height
      do x = 1, width
        rs = 0.0_real32
        gs = 0.0_real32
        bs = 0.0_real32
        as = 0.0_real32
        do offset = -radius, radius
          yy = min(height, max(1, y + offset))
          call unpack_rgba(tmp((yy - 1) * width + x), r, g, b, a)
          rs = rs + r
          gs = gs + g
          bs = bs + b
          as = as + a
        end do
        idx = (y - 1) * width + x
        output(idx) = pack_scaled(rs, gs, bs, as)
      end do
    end do
  end subroutine box_filter_host

end program main
