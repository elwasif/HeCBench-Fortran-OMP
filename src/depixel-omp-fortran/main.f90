! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use omp_lib
  implicit none

  interface
    subroutine depixel_seed_rng(seed) bind(C, name='depixel_seed_rng')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine depixel_seed_rng

    function depixel_next_random() bind(C, name='depixel_next_random') result(value)
      import :: c_float
      real(c_float) :: value
    end function depixel_next_random
  end interface

  integer, parameter :: nthreads = 256
  type, bind(C) :: float3
    real(c_float) :: x
    real(c_float) :: y
    real(c_float) :: z
    real(c_float) :: pad
  end type float3

  integer :: width, height, repeat, size, n, i, errors
  integer(int32), allocatable :: out(:), tmp(:), ref_tmp(:), ref_out(:)
  type(float3), allocatable :: img(:)
  real(real32) :: sum_value, lsum
  real(real64) :: start_time, total_time
  character(len=256) :: arg0

  if (command_argument_count() /= 3) then
    call get_command_argument(0, arg0)
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <image width> <image height> <repeat>'
    stop 1
  end if

  width = read_int_arg(1)
  height = read_int_arg(2)
  repeat = read_int_arg(3)
  if (width <= 0 .or. height <= 0 .or. repeat <= 0) stop 1

  size = width * height
  allocate(img(size), tmp(size), out(size), ref_tmp(size), ref_out(size))
  tmp = 0_int32
  out = 0_int32
  ref_tmp = 0_int32
  ref_out = 0_int32
  call depixel_seed_rng(19937_c_int)
  sum_value = 0.0_real32
  total_time = 0.0_real64
  errors = 0

  !$omp target data map(alloc: img(1:size), tmp(1:size)) map(from: out(1:size))
  do n = 1, repeat
    do i = 1, size
      img(i)%x = depixel_next_random()
      img(i)%y = depixel_next_random()
      img(i)%z = depixel_next_random()
      img(i)%pad = 0.0_c_float
    end do

    !$omp target update to(img(1:size))
    start_time = omp_get_wtime()
    call check_connect_device(img, tmp, width, height, size)
    call eliminate_crosses_device(tmp, out, width, height, size)
    total_time = total_time + (omp_get_wtime() - start_time)
    !$omp target update from(out(1:size))

    call check_connect_host(img, ref_tmp, width, height, size)
    call eliminate_crosses_host(ref_tmp, ref_out, width, height, size)
    do i = 1, size
      if (.not. pixel_matches(out(i), ref_out(i))) then
        if (errors == 0) then
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'First mismatch: iteration ', n, &
              ', index ', i, ', row ', (i - 1) / width, ', column ', mod(i - 1, width)
          write(*,'(A,Z8.8,A,Z8.8)') 'Device out=0x', out(i), ' reference=0x', ref_out(i)
        end if
        errors = errors + 1
        exit
      end if
    end do

    lsum = 0.0_real32
    do i = 1, size
      lsum = lsum + real(iand(out(i), int(z'000000ff', int32)), real32) / 256.0_real32 &
          + real(iand(shiftr(out(i), 8), int(z'000000ff', int32)), real32) / 256.0_real32 &
          + real(iand(shiftr(out(i), 16), int(z'000000ff', int32)), real32) / 256.0_real32 &
          + real(iand(shiftr(out(i), 24), int(z'000000ff', int32)), real32) / 256.0_real32
    end do
    sum_value = sum_value + lsum / real(size, real32)
  end do
  !$omp end target data

  write(*,'(A,I0,A,I0,A)') 'Image size: ', width, ' (width) x ', height, ' (height)'
  write(*,'(A,F0.6)') 'checkSum: ', sum_value
  write(*,'(A,I0,A,F0.6,A)') 'Average kernel time over ', repeat, ' iterations: ', &
      total_time / real(repeat, real64), ' (s)'
  if (errors /= 0) write(*,'(A)') 'FAIL'

  deallocate(img, tmp, out, ref_tmp, ref_out)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  logical function pixel_matches(actual, expected) result(ok)
    integer(int32), intent(in) :: actual, expected
    integer, parameter :: byte_tol = 1
    integer :: shift, actual_byte, expected_byte
    ! Connectivity is exact; packed Y/U/V payload bytes can differ by one count from host/device rounding.
    ok = iand(actual, int(z'000000ff', int32)) == iand(expected, int(z'000000ff', int32))
    if (.not. ok) return
    do shift = 8, 24, 8
      actual_byte = int(iand(shiftr(actual, shift), int(z'000000ff', int32)))
      expected_byte = int(iand(shiftr(expected, shift), int(z'000000ff', int32)))
      if (abs(actual_byte - expected_byte) > byte_tol) then
        ok = .false.
        return
      end if
    end do
  end function pixel_matches

  real(real32) function saturatef(v) result(out)
    real(real32), intent(in) :: v
    if (v < 0.0_real32) then
      out = 0.0_real32
    else if (v > 1.0_real32) then
      out = 1.0_real32
    else
      out = v
    end if
  end function saturatef

  integer(int32) function rgb_to_yuv(rgba) result(packed)
    type(float3), intent(in) :: rgba
    real(real32) :: y, u, v
    integer(int32) :: yi, ui, vi
    y = 0.299_real32 * rgba%x + 0.587_real32 * rgba%y + 0.114_real32 * rgba%z
    u = 0.713_real32 * (rgba%x - y) + 0.5_real32
    v = 0.564_real32 * (rgba%z - y) + 0.5_real32
    yi = int(saturatef(y) * 255.0_real32, int32)
    ui = int(saturatef(u) * 255.0_real32, int32)
    vi = int(saturatef(v) * 255.0_real32, int32)
    packed = ior(ior(ior(shiftl(int(255, int32), 24), shiftl(vi, 16)), shiftl(ui, 8)), yi)
  end function rgb_to_yuv

  logical function is_connected(lnode, rnode) result(ok)
    integer(int32), intent(in) :: lnode, rnode
    integer :: ly, lu, lv, ry, ru, rv
    ly = int(iand(lnode, int(z'000000ff', int32)))
    lu = int(iand(shiftr(lnode, 8), int(z'000000ff', int32)))
    lv = int(iand(shiftr(lnode, 16), int(z'000000ff', int32)))
    ry = int(iand(rnode, int(z'000000ff', int32)))
    ru = int(iand(shiftr(rnode, 8), int(z'000000ff', int32)))
    rv = int(iand(shiftr(rnode, 16), int(z'000000ff', int32)))
    ok = .not. ((abs(ly - ry) > 48) .or. (abs(lu - ru) > 7) .or. (abs(lv - rv) > 6))
  end function is_connected

  integer function bit_count32(v) result(count)
    integer(int32), intent(in) :: v
    integer(int32) :: x
    count = 0
    x = v
    do while (x /= 0_int32)
      x = iand(x, x - 1_int32)
      count = count + 1
    end do
  end function bit_count32

  integer function idx0(row, column, width) result(idx)
    integer, intent(in) :: row, column, width
    idx = row * width + column + 1
  end function idx0

  subroutine check_connect_device(rgba, connect, w, h, size)
    type(float3), intent(in) :: rgba(:)
    integer(int32), intent(inout) :: connect(:)
    integer, intent(in) :: w, h, size
    integer :: center, row, column, nr, nc
    integer(int32) :: yuv_c, yuv_n, con

    !$omp target teams distribute parallel do thread_limit(nthreads) private(row, column, nr, nc, yuv_c, yuv_n, con)
    do center = 1, size
      row = (center - 1) / w
      column = mod(center - 1, w)
      con = 0_int32
      yuv_c = rgb_to_yuv(rgba(center))

      nr = merge(row - 1, row, row > 0 .and. column > 0)
      nc = merge(column - 1, column, column > 0 .and. row > 0)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + 1_int32

      nr = merge(row - 1, row, row > 0)
      nc = column
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 1)

      nr = merge(row - 1, row, row > 0 .and. column < w - 1)
      nc = merge(column + 1, column, column < w - 1 .and. row > 0)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 2)

      nr = row
      nc = merge(column + 1, column, column < w - 1)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 3)

      nr = merge(row + 1, row, row < h - 1 .and. column < w - 1)
      nc = merge(column + 1, column, column < w - 1 .and. row < h - 1)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 4)

      nr = merge(row + 1, row, row < h - 1)
      nc = column
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 5)

      nr = merge(row + 1, row, row < h - 1 .and. column > 0)
      nc = merge(column - 1, column, column > 0 .and. row < h - 1)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 6)

      nr = row
      nc = merge(column - 1, column, column > 0)
      yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
      if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, 7)

      connect(center) = ior(ior(ior(shiftl(iand(shiftr(yuv_c, 16), int(z'000000ff', int32)), 24), &
          shiftl(iand(shiftr(yuv_c, 8), int(z'000000ff', int32)), 16)), &
          shiftl(iand(yuv_c, int(z'000000ff', int32)), 8)), con)
    end do
    !$omp end target teams distribute parallel do
  end subroutine check_connect_device

  subroutine check_connect_host(rgba, connect, w, h, size)
    type(float3), intent(in) :: rgba(:)
    integer(int32), intent(inout) :: connect(:)
    integer, intent(in) :: w, h, size
    integer :: center, row, column, nr, nc
    integer(int32) :: yuv_c, yuv_n, con
    do center = 1, size
      row = (center - 1) / w
      column = mod(center - 1, w)
      con = 0_int32
      yuv_c = rgb_to_yuv(rgba(center))
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row - 1, row, row > 0 .and. column > 0), merge(column - 1, column, column > 0 .and. row > 0), 0)
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row - 1, row, row > 0), column, 1)
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row - 1, row, row > 0 .and. column < w - 1), merge(column + 1, column, column < w - 1 .and. row > 0), 2)
      call add_neighbor(con, yuv_c, rgba, w, row, column, row, merge(column + 1, column, column < w - 1), 3)
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row + 1, row, row < h - 1 .and. column < w - 1), merge(column + 1, column, column < w - 1 .and. row < h - 1), 4)
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row + 1, row, row < h - 1), column, 5)
      call add_neighbor(con, yuv_c, rgba, w, row, column, merge(row + 1, row, row < h - 1 .and. column > 0), merge(column - 1, column, column > 0 .and. row < h - 1), 6)
      call add_neighbor(con, yuv_c, rgba, w, row, column, row, merge(column - 1, column, column > 0), 7)
      connect(center) = ior(ior(ior(shiftl(iand(shiftr(yuv_c, 16), int(z'000000ff', int32)), 24), &
          shiftl(iand(shiftr(yuv_c, 8), int(z'000000ff', int32)), 16)), &
          shiftl(iand(yuv_c, int(z'000000ff', int32)), 8)), con)
    end do
  end subroutine check_connect_host

  subroutine add_neighbor(con, yuv_c, rgba, w, row, column, nr, nc, bit)
    integer(int32), intent(inout) :: con
    integer(int32), intent(in) :: yuv_c
    type(float3), intent(in) :: rgba(:)
    integer, intent(in) :: w, row, column, nr, nc, bit
    integer(int32) :: yuv_n
    yuv_n = rgb_to_yuv(rgba(idx0(nr, nc, w)))
    if (.not. (row == nr .and. column == nc) .and. is_connected(yuv_c, yuv_n)) con = con + shiftl(1_int32, bit)
  end subroutine add_neighbor

  subroutine eliminate_crosses_device(id, od, w, h, size)
    integer(int32), intent(in) :: id(:)
    integer(int32), intent(inout) :: od(:)
    integer, intent(in) :: w, h, size
    integer :: center
    !$omp target teams distribute parallel do thread_limit(nthreads)
    do center = 1, size
      od(center) = eliminate_one(id, center, w, h)
    end do
    !$omp end target teams distribute parallel do
  end subroutine eliminate_crosses_device

  subroutine eliminate_crosses_host(id, od, w, h, size)
    integer(int32), intent(in) :: id(:)
    integer(int32), intent(inout) :: od(:)
    integer, intent(in) :: w, h, size
    integer :: center
    do center = 1, size
      od(center) = eliminate_one(id, center, w, h)
    end do
  end subroutine eliminate_crosses_host

  integer(int32) function eliminate_one(id, center, w, h) result(pixel)
    integer(int32), intent(in) :: id(:)
    integer, intent(in) :: center, w, h
    integer :: row, column, start_row, start_column, end_row, end_column
    integer :: weight_l, weight_r, sum_l, sum_r, i, j, c_row, c_column
    integer(int32) :: curve_l, curve_r, edge_l, edge_r
    logical :: remove_cross

    row = (center - 1) / w
    column = mod(center - 1, w)
    start_row = merge(row - 3, 0, row > 2)
    start_column = merge(column - 3, 0, column > 2)
    end_row = merge(row + 4, w - 1, row < w - 4)
    end_column = merge(column + 4, h - 1, column < h - 4)
    weight_l = 0
    weight_r = 0
    pixel = 0_int32
    remove_cross = .false.

    if ((row < h - 1) .and. (column < w - 1)) then
      pixel = ior(ior(ior(shiftr(iand(id(center), int(z'00000008', int32)), 3), &
          shiftl(shiftr(iand(id(center + w + 1), int(z'00000002', int32)), 1), 1)), &
          shiftl(shiftr(iand(id(center + w + 1), int(z'00000080', int32)), 7), 2)), &
          shiftl(shiftr(iand(id(center), int(z'00000020', int32)), 5), 3))

      if (iand(id(center), int(z'00000010', int32)) /= 0_int32 .and. &
          iand(id(center + 1), int(z'00000040', int32)) /= 0_int32) then
        if (iand(id(center), int(z'00000028', int32)) /= 0_int32 .and. &
            iand(id(center + 1), int(z'000000a0', int32)) /= 0_int32) then
          pixel = ior(shiftl(iand(shiftr(id(center), 8), int(z'00ffffff', int32)), 8), pixel)
          remove_cross = .true.
        else
          if (id(center) == int(z'00000010', int32)) weight_l = weight_l + 5
          if (id(center + 1) == int(z'00000040', int32)) weight_r = weight_r + 5
          sum_l = 0
          sum_r = 0
          do i = start_row, end_row
            do j = start_column, end_column
              if (i * w + j + 1 /= center .and. i * w + j + 1 /= center + 1) then
                if (is_connected(shiftr(id(center), 8), shiftr(id(idx0(i, j, w)), 8))) sum_l = sum_l + 1
                if (is_connected(shiftr(id(center + 1), 8), shiftr(id(idx0(i, j, w)), 8))) sum_r = sum_r + 1
              end if
            end do
          end do
          weight_r = weight_r + merge(sum_l - sum_r, 0, sum_l > sum_r)
          weight_l = weight_l + merge(sum_r - sum_l, 0, sum_l < sum_r)

          call trace_left(id, w, h, row, column, sum_l)
          call trace_right(id, w, h, row, column, sum_r)
          weight_l = weight_l + merge(sum_l - sum_r, 0, sum_l > sum_r)
          weight_r = weight_r + merge(sum_r - sum_l, 0, sum_l < sum_r)

          if (weight_l > weight_r) then
            pixel = ior(pixel, int(z'00000010', int32))
            pixel = ior(shiftl(iand(shiftr(id(center), 8), int(z'00ffffff', int32)), 8), pixel)
            remove_cross = .true.
          else if (weight_r > weight_l) then
            pixel = ior(pixel, int(z'00000020', int32))
            pixel = ior(shiftl(iand(shiftr(id(center), 8), int(z'00ffffff', int32)), 8), pixel)
            remove_cross = .true.
          end if
        end if
      end if
      if (.not. remove_cross) then
        pixel = ior(pixel, ior(shiftl(shiftr(iand(id(center), int(z'00000010', int32)), 4), 4), &
            shiftl(shiftr(iand(id(center + 1), int(z'00000040', int32)), 6), 5)))
      end if
    end if
    if (.not. remove_cross) pixel = ior(shiftl(iand(shiftr(id(center), 8), int(z'00ffffff', int32)), 8), pixel)
  end function eliminate_one

  subroutine trace_left(id, w, h, row, column, total)
    integer(int32), intent(in) :: id(:)
    integer, intent(in) :: w, h, row, column
    integer, intent(out) :: total
    integer :: c_row, c_column
    integer(int32) :: curve, edge
    total = 1
    c_row = row
    c_column = column
    curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
    edge = int(z'00000010', int32)
    do while (bit_count32(curve) == 2 .and. total < w * h)
      edge = curve - edge
      call step_edge(edge, c_row, c_column)
      if (c_row < 0 .or. c_row >= h .or. c_column < 0 .or. c_column >= w) exit
      edge = merge(shiftr(edge, 4), shiftl(edge, 4), edge > 8_int32)
      curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
      total = total + 1
    end do
    c_row = row + 1
    c_column = column + 1
    curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
    edge = int(z'00000001', int32)
    do while (bit_count32(curve) == 2 .and. total < w * h)
      edge = curve - edge
      call step_edge(edge, c_row, c_column)
      if (c_row < 0 .or. c_row >= h .or. c_column < 0 .or. c_column >= w) exit
      edge = merge(shiftr(edge, 4), shiftl(edge, 4), edge > 8_int32)
      curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
      total = total + 1
    end do
  end subroutine trace_left

  subroutine trace_right(id, w, h, row, column, total)
    integer(int32), intent(in) :: id(:)
    integer, intent(in) :: w, h, row, column
    integer, intent(out) :: total
    integer :: c_row, c_column
    integer(int32) :: curve, edge
    total = 1
    c_row = row
    c_column = column + 1
    curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
    edge = int(z'00000040', int32)
    do while (bit_count32(curve) == 2 .and. total < w * h)
      edge = curve - edge
      if (edge == 64_int32) then
        c_row = c_row + 1
        c_column = c_column - 1
        c_row = c_row - 1
        c_column = c_column - 1
      else
        call step_edge(edge, c_row, c_column)
      end if
      if (c_row < 0 .or. c_row >= h .or. c_column < 0 .or. c_column >= w) exit
      edge = merge(shiftr(edge, 4), shiftl(edge, 4), edge > 8_int32)
      curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
      total = total + 1
    end do
    c_row = row + 1
    c_column = column
    curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
    edge = int(z'00000004', int32)
    do while (bit_count32(curve) == 2 .and. total < w * h)
      edge = curve - edge
      call step_edge(edge, c_row, c_column)
      if (c_row < 0 .or. c_row >= h .or. c_column < 0 .or. c_column >= w) exit
      edge = merge(shiftr(edge, 4), shiftl(edge, 4), edge > 8_int32)
      curve = iand(id(idx0(c_row, c_column, w)), int(z'000000ff', int32))
      total = total + 1
    end do
  end subroutine trace_right

  subroutine step_edge(edge, c_row, c_column)
    integer(int32), intent(in) :: edge
    integer, intent(inout) :: c_row, c_column
    select case (edge)
    case (1_int32)
      c_row = c_row - 1
      c_column = c_column - 1
    case (2_int32)
      c_row = c_row - 1
    case (4_int32)
      c_row = c_row - 1
      c_column = c_column + 1
    case (8_int32)
      c_column = c_column + 1
    case (16_int32)
      c_row = c_row + 1
      c_column = c_column + 1
    case (32_int32)
      c_row = c_row + 1
    case (64_int32)
      c_row = c_row + 1
      c_column = c_column - 1
    case (128_int32)
      c_column = c_column - 1
    end select
  end subroutine step_edge

end program main
