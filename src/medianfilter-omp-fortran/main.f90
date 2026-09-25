! SPDX-License-Identifier: CC0-1.0
program medianfilter_omp_fortran
  use, intrinsic :: iso_fortran_env, only: int8, int32, real32, real64
  use omp_lib, only: omp_get_wtime, omp_get_team_num, omp_get_thread_num
  implicit none

  integer, parameter :: max_image_width = 1920
  integer, parameter :: max_image_height = 1080

  character(len=512) :: image_path, arg, program_name
  integer :: argc, cycles, ios, width, height, n_pixels, iter
  integer(int32), allocatable :: input(:), output(:), golden(:)
  logical :: status, match
  real(real64) :: total_time, warmup_time

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, program_name)
    write(*, '(A,A,A)') 'Usage: ', trim(program_name), ' <image file> <repeat>'
    stop 1
  end if

  call get_command_argument(1, image_path)
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) cycles
  if (ios /= 0) cycles = 0

  call load_ppm4ub(trim(image_path), input, width, height, status)

  write(*, '(A,A)') 'Image File', char(9)//' = '//trim(image_path)
  write(*, '(A,I0,A,I0,A,I0,A)') 'Image Dimensions = ', width, ' w x ', height, ' h x ', 32, ' bpp'
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
  allocate(output(n_pixels), golden(n_pixels))
  output = 0_int32
  golden = 0_int32

  !$omp target data map(to: input(1:n_pixels)) map(from: output(1:n_pixels))
    warmup_time = median_filter_gpu(input, output, width, height)

    write(*, *)
    write(*, '(A,I0,A)') 'Running MedianFilterGPU for ', cycles, ' cycles...'
    write(*, *)
    total_time = 0.0_real64
    do iter = 1, cycles
      total_time = total_time + median_filter_gpu(input, output, width, height)
    end do
    write(*, '(A,F8.6,A)') 'Average kernel execution time: ', total_time / real(cycles, real64), ' (s)'
    write(*, *)
  !$omp end target data

  call median_filter_host(input, golden, width, height)

  write(*, '(A)') 'Comparing GPU Result to CPU Result...'
  match = compare_uint_threshold(golden, output, n_pixels, 1.0_real32, 0.0001_real32)
  write(*, *)
  if (match) then
    write(*, '(A)') 'GPU Result matches CPU Result within tolerance...'
    write(*, '(A)') 'PASS'
  else
    write(*, '(A)') "GPU Result DOESN'T match CPU Result within tolerance..."
    write(*, '(A)') 'FAIL'
  end if

  deallocate(input, output, golden)

contains

  real(real64) function median_filter_gpu(input, output, width, height)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(inout) :: output(:)
    integer, intent(in) :: width, height
    integer, parameter :: iBlockDimX = 16
    integer, parameter :: iBlockDimY = 4
    integer, parameter :: iLocalPixPitch = iBlockDimX + 2
    integer :: szLocalWorkSize(2), szGlobalWorkSize(2)
    integer :: iTeamX, iTeamY, iNumTeams, iNumThreads
    real(real64) :: start_time, end_time

    szLocalWorkSize(1) = iBlockDimX
    szLocalWorkSize(2) = iBlockDimY
    szGlobalWorkSize(1) = round_up(szLocalWorkSize(1), width)
    szGlobalWorkSize(2) = round_up(szLocalWorkSize(2), height)

    iTeamX = szGlobalWorkSize(1) / szLocalWorkSize(1)
    iTeamY = szGlobalWorkSize(2) / szLocalWorkSize(2)
    iNumTeams = iTeamX * iTeamY
    iNumThreads = iBlockDimX * iBlockDimY

    start_time = omp_get_wtime()
    !$omp target teams num_teams(iNumTeams) thread_limit(iNumThreads)
    block
      integer(int32) :: uc4LocalData(iLocalPixPitch * (iBlockDimY + 2))

      !$omp parallel
      block
        integer :: iLocalIdX, iLocalIdY, iGroupIdX, iGroupIdY
        integer :: iBlockX, iBlockY, iImagePosX, iDevYPrime, iImageX
        integer :: iDevGMEMOffset, iLocalPixOffset, iSearch, channel
        real(real32) :: fMedianEstimate(3), fMinBound(3), fMaxBound(3)
        integer :: uiHighCount(3)
        integer(int32) :: uiPackedPix

        iLocalIdX = mod(omp_get_thread_num(), iBlockDimX)
        iLocalIdY = omp_get_thread_num() / iBlockDimX
        iGroupIdX = mod(omp_get_team_num(), iTeamX)
        iGroupIdY = omp_get_team_num() / iTeamX
        iBlockX = iBlockDimX
        iBlockY = iBlockDimY
        iImagePosX = iGroupIdX * iBlockX + iLocalIdX
        iDevYPrime = iGroupIdY * iBlockY + iLocalIdY - 1
        iImageX = width

        iDevGMEMOffset = iDevYPrime * iImageX + iImagePosX
        iLocalPixOffset = iLocalIdY * iLocalPixPitch + iLocalIdX + 1

        if ((iDevYPrime > -1) .and. (iDevYPrime < height) .and. (iImagePosX < width)) then
          uc4LocalData(iLocalPixOffset + 1) = input(iDevGMEMOffset + 1)
        else
          uc4LocalData(iLocalPixOffset + 1) = 0_int32
        end if

        if (iLocalIdY < 2) then
          iLocalPixOffset = iLocalPixOffset + iBlockY * iLocalPixPitch
          if (((iDevYPrime + iBlockY) < height) .and. (iImagePosX < width)) then
            uc4LocalData(iLocalPixOffset + 1) = input(iDevGMEMOffset + iBlockY * iImageX + 1)
          else
            uc4LocalData(iLocalPixOffset + 1) = 0_int32
          end if
        end if

        if (iLocalIdX == (iBlockX - 1)) then
          iLocalPixOffset = iLocalIdY * iLocalPixPitch

          if ((iDevYPrime > -1) .and. (iDevYPrime < height) .and. (iGroupIdX > 0)) then
            uc4LocalData(iLocalPixOffset + 1) = input(iDevYPrime * iImageX + iGroupIdX * iBlockX)
          else
            uc4LocalData(iLocalPixOffset + 1) = 0_int32
          end if

          if (iLocalIdY < 2) then
            iLocalPixOffset = iLocalPixOffset + iBlockY * iLocalPixPitch

            if (((iDevYPrime + iBlockY) < height) .and. (iGroupIdX > 0)) then
              uc4LocalData(iLocalPixOffset + 1) = input((iDevYPrime + iBlockY) * iImageX + &
                                                        iGroupIdX * iBlockX)
            else
              uc4LocalData(iLocalPixOffset + 1) = 0_int32
            end if
          end if
        else if (iLocalIdX == 0) then
          iLocalPixOffset = (iLocalIdY + 1) * iLocalPixPitch - 1

          if ((iDevYPrime > -1) .and. (iDevYPrime < height) .and. &
              ((iGroupIdX + 1) * iBlockX < width)) then
            uc4LocalData(iLocalPixOffset + 1) = input(iDevYPrime * iImageX + &
                                                      (iGroupIdX + 1) * iBlockX + 1)
          else
            uc4LocalData(iLocalPixOffset + 1) = 0_int32
          end if

          if (iLocalIdY < 2) then
            iLocalPixOffset = iLocalPixOffset + iBlockY * iLocalPixPitch

            if (((iDevYPrime + iBlockY) < height) .and. &
                ((iGroupIdX + 1) * iBlockX < width)) then
              uc4LocalData(iLocalPixOffset + 1) = input((iDevYPrime + iBlockY) * iImageX + &
                                                        (iGroupIdX + 1) * iBlockX + 1)
            else
              uc4LocalData(iLocalPixOffset + 1) = 0_int32
            end if
          end if
        end if

        !$omp barrier

        fMedianEstimate = 128.0_real32
        fMinBound = 0.0_real32
        fMaxBound = 255.0_real32

        do iSearch = 1, 8
          uiHighCount = 0

          iLocalPixOffset = iLocalIdY * iLocalPixPitch + iLocalIdX
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)

          iLocalPixOffset = iLocalPixOffset + (iLocalPixPitch - 2)
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)

          iLocalPixOffset = iLocalPixOffset + (iLocalPixPitch - 2)
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)
          iLocalPixOffset = iLocalPixOffset + 1
          call count_high_pixel(uc4LocalData(iLocalPixOffset + 1), fMedianEstimate, uiHighCount)

          do channel = 1, 3
            if (uiHighCount(channel) > 4) then
              fMinBound(channel) = fMedianEstimate(channel)
            else
              fMaxBound(channel) = fMedianEstimate(channel)
            end if
            fMedianEstimate(channel) = 0.5_real32 * (fMaxBound(channel) + fMinBound(channel))
          end do
        end do

        uiPackedPix = pack_pixel(fMedianEstimate)

        if ((iDevYPrime < height) .and. (iImagePosX < width)) then
          output(iDevGMEMOffset + iImageX + 1) = uiPackedPix
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
    end_time = omp_get_wtime()
    median_filter_gpu = end_time - start_time
  end function median_filter_gpu

  integer function round_up(group_size, global_size)
    integer, intent(in) :: group_size, global_size
    integer :: remainder

    remainder = mod(global_size, group_size)
    if (remainder == 0) then
      round_up = global_size
    else
      round_up = global_size + group_size - remainder
    end if
  end function round_up

  subroutine count_high_pixel(pixel, fMedianEstimate, uiHighCount)
    integer(int32), intent(in) :: pixel
    real(real32), intent(in) :: fMedianEstimate(3)
    integer, intent(inout) :: uiHighCount(3)

    if (fMedianEstimate(1) < real(iand(pixel, int(z'000000ff', int32)), real32)) uiHighCount(1) = uiHighCount(1) + 1
    if (fMedianEstimate(2) < real(iand(ishft(pixel, -8), int(z'000000ff', int32)), real32)) uiHighCount(2) = uiHighCount(2) + 1
    if (fMedianEstimate(3) < real(iand(ishft(pixel, -16), int(z'000000ff', int32)), real32)) uiHighCount(3) = uiHighCount(3) + 1
  end subroutine count_high_pixel

  integer(int32) function pack_pixel(fMedianEstimate)
    real(real32), intent(in) :: fMedianEstimate(3)

    pack_pixel = ior(ior(iand(int(fMedianEstimate(1) + 0.5_real32, int32), int(z'000000ff', int32)), &
                         ishft(iand(int(fMedianEstimate(2) + 0.5_real32, int32), int(z'000000ff', int32)), 8)), &
                     ishft(iand(int(fMedianEstimate(3) + 0.5_real32, int32), int(z'000000ff', int32)), 16))
  end function pack_pixel

  subroutine median_filter_host(input, output, width, height)
    integer(int32), intent(in) :: input(:)
    integer(int32), intent(out) :: output(:)
    integer, intent(in) :: width, height
    integer :: x, y

    do y = 1, height
      do x = 1, width
        output((y - 1) * width + x) = median_pixel(input, width, height, x, y)
      end do
    end do
  end subroutine median_filter_host

  integer(int32) function median_pixel(input, width, height, x, y)
    integer(int32), intent(in) :: input(:)
    integer, intent(in) :: width, height, x, y
    real(real32) :: estimate(3), min_bound(3), max_bound(3)
    integer :: search, row, col, channel, count_high(3), pixel, idx

    estimate = 128.0_real32
    min_bound = 0.0_real32
    max_bound = 255.0_real32

    do search = 1, 8
      count_high = 0
      do row = y - 1, y + 1
        do col = x - 1, x + 1
          pixel = 0_int32
          if (col >= 1 .and. col <= width .and. row >= 1 .and. row <= height) then
            idx = (row - 1) * width + col
            pixel = input(idx)
          end if
          if (estimate(1) < real(iand(pixel, int(z'000000ff', int32)), real32)) count_high(1) = count_high(1) + 1
          if (estimate(2) < real(iand(ishft(pixel, -8), int(z'000000ff', int32)), real32)) count_high(2) = count_high(2) + 1
          if (estimate(3) < real(iand(ishft(pixel, -16), int(z'000000ff', int32)), real32)) count_high(3) = count_high(3) + 1
        end do
      end do

      do channel = 1, 3
        if (count_high(channel) > 4) then
          min_bound(channel) = estimate(channel)
        else
          max_bound(channel) = estimate(channel)
        end if
        estimate(channel) = 0.5_real32 * (max_bound(channel) + min_bound(channel))
      end do
    end do

    median_pixel = ior(ior(iand(int(estimate(1) + 0.5_real32, int32), int(z'000000ff', int32)), &
                           ishft(iand(int(estimate(2) + 0.5_real32, int32), int(z'000000ff', int32)), 8)), &
                       ishft(iand(int(estimate(3) + 0.5_real32, int32), int(z'000000ff', int32)), 16))
  end function median_pixel

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
      if (error_count > 0) then
        write(*, '(A,F4.2,A,I0,A)') '    ', real(error_count, real32) * 100.0_real32 / real(n, real32), &
          '(%) of bytes mismatched (count=', error_count, ')'
      end if
      compare_uint_threshold = (real(n, real32) * threshold > real(error_count, real32))
    end if
  end function compare_uint_threshold

end program medianfilter_omp_fortran
