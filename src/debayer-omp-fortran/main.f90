! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_signed_char
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: tile_rows = 5
  integer, parameter :: tile_cols = 32
  integer, parameter :: kernel_size = 5
  integer, parameter :: half_ksize = kernel_size / 2
  integer, parameter :: apron_rows = tile_rows + kernel_size - 1
  integer, parameter :: apron_cols = tile_cols + kernel_size - 1
  integer, parameter :: n_tile_pixels = tile_rows * tile_cols
  integer, parameter :: n_apron_fill_tasks = apron_rows * apron_cols
  integer, parameter :: rggb = 0

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
  integer :: width, height, repeat, num_pix, input_image_pitch, output_image_pitch
  integer :: teamX, teamY
  integer(c_signed_char), allocatable :: input(:), output(:), reference(:)
  integer :: i, rand_value
  integer(int64) :: checksum
  real(real64) :: start_time, elapsed

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <width> <height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) width
  call get_command_argument(2, arg); read(arg, *) height
  call get_command_argument(3, arg); read(arg, *) repeat
  if (width <= kernel_size .or. height <= kernel_size .or. repeat <= 0) stop 1

  input_image_pitch = width
  output_image_pitch = width * 4
  num_pix = width * height
  allocate(input(num_pix), output(4 * num_pix), reference(4 * num_pix))

  call c_srand(123_c_int)
  do i = 1, num_pix
    rand_value = modulo(c_rand(), 256_c_int)
    input(i) = byte_from_uchar(rand_value)
  end do
  teamX = (width + tile_cols - 1) / tile_cols
  teamY = (height + tile_rows - 1) / tile_rows
  output = byte_from_uchar(0)
  reference = byte_from_uchar(0)

  !$omp target data map(to: input(1:num_pix)) map(from: output(1:4*num_pix))
  start_time = omp_get_wtime()
  do i = 1, repeat
    output = byte_from_uchar(0)
    !$omp target update to(output(1:4*num_pix))
    call malvar_he_cutler_demosaic(teamX, teamY, height, width, input, input_image_pitch, output, &
      output_image_pitch, rggb)
  end do
  elapsed = omp_get_wtime() - start_time
  !$omp end target data

  call reference_demosaic(height, width, input, input_image_pitch, reference, output_image_pitch, rggb)
  if (any(output /= reference)) error stop 'Fortran reference validation failed'

  write(*,'(A,F0.6,A)') 'Average kernel execution time ', elapsed / real(repeat, real64), ' (s)'

  checksum = 0_int64
  do i = 1, num_pix
    checksum = checksum + int(unsigned_byte(output(i)), int64)
  end do
  write(*,'(A,I0)') 'Checksum: ', checksum

  deallocate(input, output, reference)

contains

  subroutine malvar_he_cutler_demosaic(teamX, teamY, height, width, input_image, input_pitch, output_image, &
      output_pitch, bayer_pattern)
    integer, intent(in) :: teamX, teamY, height, width, input_pitch, output_pitch, bayer_pattern
    integer(c_signed_char), intent(in) :: input_image(:)
    integer(c_signed_char), intent(inout) :: output_image(:)

    !$omp target teams num_teams(teamX * teamY) thread_limit(tile_cols * tile_rows)
    block
      integer :: apron(apron_rows * apron_cols)

      !$omp parallel
      block
        integer :: tile_col_blocksize, tile_row_blocksize
        integer :: tile_col_block, tile_row_block, tile_col, tile_row
        integer :: g_c, g_r, tile_flat_id, apron_fill_task_id
        integer :: apron_read_row, apron_read_col, ag_c, ag_r, a_c, a_r
        logical :: valid_pixel_task

        tile_col_blocksize = tile_cols
        tile_row_blocksize = tile_rows
        tile_col_block = mod(omp_get_team_num(), teamX)
        tile_row_block = omp_get_team_num() / teamX
        tile_col = mod(omp_get_thread_num(), tile_cols)
        tile_row = omp_get_thread_num() / tile_cols
        g_c = tile_col_blocksize * tile_col_block + tile_col
        g_r = tile_row_blocksize * tile_row_block + tile_row
        valid_pixel_task = (g_r < height) .and. (g_c < width)

        tile_flat_id = tile_row * tile_cols + tile_col
        do apron_fill_task_id = tile_flat_id, n_apron_fill_tasks - 1, n_tile_pixels
          apron_read_row = apron_fill_task_id / apron_cols
          apron_read_col = mod(apron_fill_task_id, apron_cols)
          ag_c = apron_read_col + tile_col_block * tile_col_blocksize - half_ksize
          ag_r = apron_read_row + tile_row_block * tile_row_blocksize - half_ksize
          apron(apron_read_row * apron_cols + apron_read_col + 1) = &
            sample_pixel(input_image, height, width, input_pitch, ag_r, ag_c)
        end do

        !$omp barrier

        a_c = tile_col + half_ksize
        a_r = tile_row + half_ksize
        if (valid_pixel_task) then
          call demosaic_pixel(g_r, g_c, a_r, a_c, apron, output_image, output_pitch, bayer_pattern)
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine malvar_he_cutler_demosaic

  subroutine reference_demosaic(height, width, input_image, input_pitch, output_image, &
      output_pitch, bayer_pattern)
    integer, intent(in) :: height, width, input_pitch, output_pitch, bayer_pattern
    integer(c_signed_char), intent(in) :: input_image(:)
    integer(c_signed_char), intent(inout) :: output_image(:)
    integer :: row, col, apron_row, apron_col
    integer :: apron(apron_rows * apron_cols)

    do row = 0, height - 1
      do col = 0, width - 1
        apron = 0
        do apron_row = 0, kernel_size - 1
          do apron_col = 0, kernel_size - 1
            apron(apron_row * apron_cols + apron_col + 1) = &
              sample_pixel(input_image, height, width, input_pitch, &
                row + apron_row - half_ksize, col + apron_col - half_ksize)
          end do
        end do
        call demosaic_pixel(row, col, half_ksize, half_ksize, apron, output_image, &
          output_pitch, bayer_pattern)
      end do
    end do
  end subroutine reference_demosaic

  subroutine demosaic_pixel(g_r, g_c, a_r, a_c, apron, output_image, output_pitch, bayer_pattern)
    integer, intent(in) :: g_r, g_c, a_r, a_c, output_pitch, bayer_pattern
    integer, intent(in) :: apron(:)
    integer(c_signed_char), intent(inout) :: output_image(:)
    integer :: f_ij, r1, r2, r3, r4
    integer :: green_at_red_or_blue, red_at_green_in_red, red_at_green_in_blue
    integer :: blue_at_green_in_red, blue_at_green_in_blue, red_at_blue, blue_at_red
    integer :: r_mod_2, c_mod_2, red_col, red_row, blue_col, blue_row
    logical :: in_red_row, in_blue_row, is_red_pixel, is_blue_pixel, is_green_pixel
    integer :: red_value, green_value, blue_value, out_idx

    f_ij = apron_pixel(apron, a_r, a_c)

    r1 = (4 * ap(apron, a_r, a_c) + &
      2 * (ap(apron, a_r, a_c - 1) + ap(apron, a_r - 1, a_c) + &
           ap(apron, a_r, a_c + 1) + ap(apron, a_r + 1, a_c)) - &
      ap(apron, a_r, a_c - 2) - ap(apron, a_r, a_c + 2) - &
      ap(apron, a_r - 2, a_c) - ap(apron, a_r + 2, a_c)) / 8

    r2 = (8 * (ap(apron, a_r, a_c - 1) + ap(apron, a_r, a_c + 1)) + &
      10 * ap(apron, a_r, a_c) + ap(apron, a_r - 2, a_c) + ap(apron, a_r + 2, a_c) - &
      2 * (ap(apron, a_r - 1, a_c - 1) + ap(apron, a_r - 1, a_c + 1) + &
           ap(apron, a_r + 1, a_c - 1) + ap(apron, a_r + 1, a_c + 1) + &
           ap(apron, a_r, a_c - 2) + ap(apron, a_r, a_c + 2))) / 16

    r3 = (8 * (ap(apron, a_r - 1, a_c) + ap(apron, a_r + 1, a_c)) + &
      10 * ap(apron, a_r, a_c) + ap(apron, a_r, a_c - 2) + ap(apron, a_r, a_c + 2) - &
      2 * (ap(apron, a_r - 1, a_c - 1) + ap(apron, a_r - 1, a_c + 1) + &
           ap(apron, a_r + 1, a_c - 1) + ap(apron, a_r + 1, a_c + 1) + &
           ap(apron, a_r - 2, a_c) + ap(apron, a_r + 2, a_c))) / 16

    r4 = (12 * ap(apron, a_r, a_c) - &
      3 * (ap(apron, a_r, a_c - 2) + ap(apron, a_r, a_c + 2) + &
           ap(apron, a_r - 2, a_c) + ap(apron, a_r + 2, a_c)) + &
      4 * (ap(apron, a_r - 1, a_c - 1) + ap(apron, a_r - 1, a_c + 1) + &
           ap(apron, a_r + 1, a_c - 1) + ap(apron, a_r + 1, a_c + 1))) / 16

    green_at_red_or_blue = r1
    red_at_green_in_red = r2
    blue_at_green_in_blue = r2
    red_at_green_in_blue = r3
    blue_at_green_in_red = r3
    red_at_blue = r4
    blue_at_red = r4

    r_mod_2 = iand(g_r, 1)
    c_mod_2 = iand(g_c, 1)
    red_col = merge(1, 0, bayer_pattern == 1 .or. bayer_pattern == 3)
    red_row = merge(1, 0, bayer_pattern == 2 .or. bayer_pattern == 3)
    blue_col = 1 - red_col
    blue_row = 1 - red_row

    in_red_row = r_mod_2 == red_row
    in_blue_row = r_mod_2 == blue_row
    is_red_pixel = (r_mod_2 == red_row) .and. (c_mod_2 == red_col)
    is_blue_pixel = (r_mod_2 == blue_row) .and. (c_mod_2 == blue_col)
    is_green_pixel = .not. (is_red_pixel .or. is_blue_pixel)

    red_value = bool_int(is_red_pixel) * f_ij + bool_int(is_blue_pixel) * red_at_blue + &
      bool_int(is_green_pixel .and. in_red_row) * red_at_green_in_red + &
      bool_int(is_green_pixel .and. in_blue_row) * red_at_green_in_blue
    blue_value = bool_int(is_blue_pixel) * f_ij + bool_int(is_red_pixel) * blue_at_red + &
      bool_int(is_green_pixel .and. in_red_row) * blue_at_green_in_red + &
      bool_int(is_green_pixel .and. in_blue_row) * blue_at_green_in_blue
    green_value = bool_int(is_green_pixel) * f_ij + bool_int(.not. is_green_pixel) * green_at_red_or_blue

    out_idx = g_r * output_pitch + g_c * 4 + 1
    output_image(out_idx) = byte_from_uchar(saturate_uchar(red_value))
    output_image(out_idx + 1) = byte_from_uchar(saturate_uchar(green_value))
    output_image(out_idx + 2) = byte_from_uchar(saturate_uchar(blue_value))
    output_image(out_idx + 3) = byte_from_uchar(0)
  end subroutine demosaic_pixel

  integer function ap(apron, row, col) result(value)
    integer, intent(in) :: apron(:), row, col
    value = apron_pixel(apron, row, col)
  end function ap

  integer function apron_pixel(apron, row, col) result(value)
    integer, intent(in) :: apron(:), row, col
    value = apron(row * apron_cols + col + 1)
  end function apron_pixel

  integer function sample_pixel(input_image, height, width, pitch, row, col) result(value)
    integer(c_signed_char), intent(in) :: input_image(:)
    integer, intent(in) :: height, width, pitch, row, col
    integer :: rr, cc

    rr = reflect_exclusive(row, height)
    cc = reflect_exclusive(col, width)
    value = unsigned_byte(input_image(rr * pitch + cc + 1))
  end function sample_pixel

  integer function reflect_exclusive(coord, limit) result(value)
    integer, intent(in) :: coord, limit

    value = coord
    if (value < 0) then
      value = -value
    end if
    if (value >= limit) then
      value = limit - (value - limit) - 2
    end if
  end function reflect_exclusive

  integer function saturate_uchar(value) result(out)
    integer, intent(in) :: value

    if (value > 255) then
      out = 255
    else if (value < 0) then
      out = 0
    else
      out = value
    end if
  end function saturate_uchar

  integer function bool_int(value) result(out)
    logical, intent(in) :: value
    out = merge(1, 0, value)
  end function bool_int

  integer(c_signed_char) function byte_from_uchar(value) result(out)
    integer, intent(in) :: value
    integer :: wrapped

    wrapped = modulo(value, 256)
    if (wrapped >= 128) wrapped = wrapped - 256
    out = int(wrapped, c_signed_char)
  end function byte_from_uchar

  integer function unsigned_byte(value) result(out)
    integer(c_signed_char), intent(in) :: value

    out = int(value)
    if (out < 0) out = out + 256
  end function unsigned_byte

end program main
