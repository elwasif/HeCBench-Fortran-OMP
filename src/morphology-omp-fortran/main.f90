! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer(int32), parameter :: black = 0_int32
  integer(int32), parameter :: white = 255_int32
  character(len=256) :: arg0, arg1, arg2, arg3, arg4, arg5
  integer(int32) :: hsize, vsize, width, height, repeat
  integer(int32), allocatable :: src_img(:), tmp_img(:)
  integer(int32) :: i, j, n, total
  real(real64) :: dilate_time, erode_time

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 5) then
    print '(3A)', 'Usage: ', trim(arg0), ' <kernel width> <kernel height> ', &
      '<image width> <image height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  read(arg1, *) hsize
  read(arg2, *) vsize
  read(arg3, *) width
  read(arg4, *) height
  read(arg5, *) repeat

  if (hsize <= 0_int32 .or. vsize <= 0_int32 .or. width <= 0_int32 .or. height <= 0_int32) stop 1
  if (repeat <= 0_int32) stop 1

  allocate(src_img(width * height), tmp_img(width * height))
  src_img = black
  tmp_img = black

  do i = 1, height
    do j = 1, width
      if (i == height / 2 .and. j == width / 2) then
        src_img(index_1d(i, j, width)) = white
      end if
    end do
  end do

  dilate_time = 0.0_real64
  erode_time = 0.0_real64

  !$omp target data map(tofrom: src_img(1:width * height)) map(alloc: tmp_img(1:width * height))
  do n = 1, repeat
    dilate_time = dilate_time + morphology_pass(src_img, tmp_img, width, height, hsize, vsize, .true.)
    erode_time = erode_time + morphology_pass(src_img, tmp_img, width, height, hsize, vsize, .false.)
  end do
  !$omp end target data

  print '(A,F8.6,A)', 'Average kernel execution time (dilate): ', dilate_time / real(repeat, real64), ' (s)'
  print '(A,F8.6,A)', 'Average kernel execution time (erode): ', erode_time / real(repeat, real64), ' (s)'

  total = sum(src_img)
  if (total == white) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    stop 1
  end if

  deallocate(src_img, tmp_img)

contains

  integer(int32) function index_1d(row, col, width_value) result(idx)
    integer(int32), intent(in) :: row, col, width_value

    idx = (row - 1_int32) * width_value + col
  end function index_1d

  real(real64) function morphology_pass(img, tmp, width, height, hsize, vsize, is_dilate) result(elapsed)
    integer(int32), intent(inout) :: img(:), tmp(:)
    integer(int32), intent(in) :: width, height, hsize, vsize
    logical, intent(in) :: is_dilate
    real(real64) :: start_time, end_time

    call clear_image(tmp, width, height)
    start_time = omp_get_wtime()
    call horizontal_pass(img, tmp, width, height, hsize, is_dilate)
    call vertical_pass(tmp, img, width, height, vsize, is_dilate)
    end_time = omp_get_wtime()
    elapsed = end_time - start_time
  end function morphology_pass

  subroutine clear_image(img, width, height)
    integer(int32), intent(inout) :: img(:)
    integer(int32), intent(in) :: width, height
    integer(int32) :: p

    !$omp target teams distribute parallel do thread_limit(256)
    do p = 1, width * height
      img(p) = black
    end do
    !$omp end target teams distribute parallel do
  end subroutine clear_image

  pure integer(int32) function round_up(x, y) result(value)
    integer(int32), intent(in) :: x, y

    value = (x + y - 1_int32) / y
  end function round_up

  pure integer(int32) function element_op(a, b, is_dilate) result(value)
    integer(int32), intent(in) :: a, b
    logical, intent(in) :: is_dilate

    if (is_dilate) then
      value = max(a, b)
    else
      value = min(a, b)
    end if
  end function element_op

  pure integer(int32) function border_value(is_dilate) result(value)
    logical, intent(in) :: is_dilate

    value = merge(white, black, is_dilate)
  end function border_value

  subroutine two_way_scan(sMem, selSize, tid, is_dilate)
    integer(int32), intent(inout) :: sMem(128)
    integer(int32), intent(in) :: selSize, tid
    logical, intent(in) :: is_dilate
    integer(int32) :: offset

    sMem(tid + 2_int32 * selSize + 1_int32) = sMem(tid + 1_int32)
    sMem(tid + 3_int32 * selSize + 1_int32) = sMem(tid + selSize + 1_int32)
    !$omp barrier

    offset = 1_int32
    do while (offset < selSize)
      if (tid >= offset) then
        sMem(tid + 3_int32 * selSize) = element_op(sMem(tid + 3_int32 * selSize), &
          sMem(tid + 3_int32 * selSize - offset), is_dilate)
      end if
      if (tid <= selSize - 1_int32 - offset) then
        sMem(tid + 2_int32 * selSize + 1_int32) = element_op( &
          sMem(tid + 2_int32 * selSize + 1_int32), &
          sMem(tid + 2_int32 * selSize + 1_int32 + offset), is_dilate)
      end if
      !$omp barrier
      offset = offset * 2_int32
    end do
  end subroutine two_way_scan

  subroutine horizontal_pass(src, dst, width, height, hsize, is_dilate)
    integer(int32), intent(in) :: src(:)
    integer(int32), intent(inout) :: dst(:)
    integer(int32), intent(in) :: width, height, hsize
    logical, intent(in) :: is_dilate
    integer(int32) :: blockSize_x_h, blockSize_y_h, gridSize_x_h, gridSize_y_h

    blockSize_x_h = hsize
    blockSize_y_h = 1_int32
    gridSize_x_h = round_up(width, blockSize_x_h)
    gridSize_y_h = round_up(height, blockSize_y_h)

    !$omp target teams num_teams(gridSize_x_h * gridSize_y_h) thread_limit(blockSize_x_h * blockSize_y_h)
    block
      integer(int32) :: sMem(128)

      !$omp parallel
      block
        integer(int32) :: bx, by, tx, tidx, tidy

        bx = mod(omp_get_team_num(), gridSize_x_h)
        by = omp_get_team_num() / gridSize_x_h
        tx = omp_get_thread_num()

        tidx = tx + bx * blockSize_x_h
        tidy = by * blockSize_y_h

        if (tidx < width .and. tidy < height) then
          sMem(tx + 1_int32) = src(tidy * width + tidx + 1_int32)
          if (tidx + hsize < width) then
            sMem(tx + hsize + 1_int32) = src(tidy * width + tidx + hsize + 1_int32)
          end if
        end if
        !$omp barrier

        if (tidx < width .and. tidy < height) then
          call two_way_scan(sMem, hsize, tx, is_dilate)
        end if

        if (tidx < width .and. tidy < height) then
          if (tidx + hsize / 2_int32 < width - hsize / 2_int32) then
            dst(tidy * width + tidx + hsize / 2_int32 + 1_int32) = element_op( &
              sMem(tx + 2_int32 * hsize + 1_int32), sMem(tx + 3_int32 * hsize), is_dilate)
          end if
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine horizontal_pass

  subroutine vertical_pass(src, dst, width, height, vsize, is_dilate)
    integer(int32), intent(in) :: src(:)
    integer(int32), intent(inout) :: dst(:)
    integer(int32), intent(in) :: width, height, vsize
    logical, intent(in) :: is_dilate
    integer(int32) :: blockSize_x_v, blockSize_y_v, gridSize_x_v, gridSize_y_v

    blockSize_x_v = 1_int32
    blockSize_y_v = vsize
    gridSize_x_v = round_up(width, blockSize_x_v)
    gridSize_y_v = round_up(height, blockSize_y_v)

    !$omp target teams num_teams(gridSize_x_v * gridSize_y_v) thread_limit(blockSize_x_v * blockSize_y_v)
    block
      integer(int32) :: sMem(128)

      !$omp parallel
      block
        integer(int32) :: bx, by, ty, tidx, tidy

        bx = mod(omp_get_team_num(), gridSize_x_v)
        by = omp_get_team_num() / gridSize_x_v
        ty = omp_get_thread_num()

        tidx = bx * blockSize_x_v
        tidy = ty + by * blockSize_y_v

        if (tidx < width .and. tidy < height) then
          sMem(ty + 1_int32) = src(tidy * width + tidx + 1_int32)
          if (tidy + vsize < height) then
            sMem(ty + vsize + 1_int32) = src((tidy + vsize) * width + tidx + 1_int32)
          end if
        end if
        !$omp barrier

        if (tidx < width .and. tidy < height) then
          call two_way_scan(sMem, vsize, ty, is_dilate)
        end if

        if (tidx < width .and. tidy < height) then
          if (tidy + vsize / 2_int32 < height - vsize / 2_int32) then
            dst((tidy + vsize / 2_int32) * width + tidx + 1_int32) = element_op( &
              sMem(ty + 2_int32 * vsize + 1_int32), sMem(ty + 3_int32 * vsize), is_dilate)
          end if
          if (tidy < vsize / 2_int32 .or. tidy >= height - vsize / 2_int32) then
            dst(tidy * width + tidx + 1_int32) = border_value(is_dilate)
          end if
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine vertical_pass

end program main
