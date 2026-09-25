! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: max_detections = 4096
  integer, parameter :: n_partitions = 32

  type :: float4
    real(real32) :: x
    real(real32) :: y
    real(real32) :: z
    real(real32) :: w
  end type float4

  character(len=512) :: input_file, output_file
  integer :: repeat, ndetections, limit, threads, totaldets
  type(float4), allocatable :: points(:)
  integer(int8), allocatable :: pointsbitmap(:), nmsbitmap(:)
  real(real64) :: start_time, end_time

  if (command_argument_count() /= 3) then
    call print_help()
    stop
  end if

  call get_command_argument(1, input_file)
  call get_command_argument(2, output_file)
  repeat = read_arg(3)

  allocate(points(max_detections))
  allocate(pointsbitmap(max_detections), nmsbitmap(max_detections * max_detections))
  points%x = 0.0_real32
  points%y = 0.0_real32
  points%z = 0.0_real32
  points%w = 0.0_real32
  pointsbitmap = 0_int8
  nmsbitmap = 1_int8

  call read_points(trim(input_file), points, ndetections)
  print '(A,A,A,I0)', 'Number of detections read from input file (', trim(input_file), '): ', ndetections

  limit = get_upper_limit(ndetections, 16)
  threads = get_optimal_dim(limit) * get_optimal_dim(limit)

  !$omp target data map(to: points(1:max_detections), &
  !$omp& nmsbitmap(1:max_detections * max_detections)) map(tofrom: pointsbitmap(1:max_detections))
  start_time = omp_get_wtime()
  call generate_nms_bitmap(points, nmsbitmap, limit, repeat, threads)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'Average kernel execution time (generate_nms_bitmap): ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  start_time = omp_get_wtime()
  call reduce_nms_bitmap(nmsbitmap, pointsbitmap, ndetections, repeat)
  end_time = omp_get_wtime()
  print '(A,F0.6,A)', 'Average kernel execution time (reduce_nms_bitmap): ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'
  !$omp end target data

  call write_points(trim(output_file), points, pointsbitmap, ndetections, totaldets)
  print '(A,I0)', 'Detections after NMS: ', totaldets

  deallocate(points, pointsbitmap, nmsbitmap)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine print_help()
    print '(A)', ''
    print '(A)', 'Usage: nmstest  <detections.txt>  <output.txt>'
    print '(A)', ''
    print '(A)', '               detections.txt -> Input file containing the coordinates, width, and scores of detected objects'
    print '(A)', '               output.txt     -> Output file after performing NMS'
    print '(A)', '               repeat         -> Kernel execution count'
    print '(A)', ''
  end subroutine print_help

  integer function get_optimal_dim(val)
    integer, intent(in) :: val
    integer :: div, neg, cntneg, cntpos, i
    neg = 1
    div = 16
    cntneg = div
    cntpos = div
    do i = 1, 5
      if (mod(val, div) == 0) then
        get_optimal_dim = div
        return
      end if
      if (neg /= 0) then
        cntneg = cntneg - 1
        div = cntneg
        neg = 0
      else
        cntpos = cntpos + 1
        div = cntpos
        neg = 1
      end if
    end do
    get_optimal_dim = 16
  end function get_optimal_dim

  integer function get_upper_limit(val, mul)
    integer, intent(in) :: val, mul
    integer :: cnt
    cnt = mul
    do while (cnt < val)
      cnt = cnt + mul
    end do
    if (cnt > max_detections) cnt = max_detections
    get_upper_limit = cnt
  end function get_upper_limit

  subroutine read_points(path, points, ndetections)
    character(len=*), intent(in) :: path
    type(float4), intent(inout) :: points(:)
    integer, intent(out) :: ndetections
    character(len=256) :: line
    integer :: unit, ios, k, x, y, w
    real(real32) :: score

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(A,A,A)', 'Error: Unable to open file ', path, ' for input detection coordinates.'
      stop 1
    end if

    ndetections = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      do k = 1, len_trim(line)
        if (line(k:k) == ',') line(k:k) = ' '
      end do
      read(line, *, iostat=ios) x, y, w, score
      if (ios /= 0) then
        print '(A,I0,A,A)', 'Error: Invalid file format in line ', ndetections, ' when reading ', path
        stop 1
      end if
      ndetections = ndetections + 1
      if (ndetections > max_detections) exit
      points(ndetections)%x = real(x, real32)
      points(ndetections)%y = real(y, real32)
      points(ndetections)%z = real(w, real32)
      points(ndetections)%w = score
    end do
    close(unit)
  end subroutine read_points

  subroutine generate_nms_bitmap(points, nmsbitmap, limit, repeat, threads)
    type(float4), intent(in) :: points(:)
    integer(int8), intent(inout) :: nmsbitmap(:)
    integer, intent(in) :: limit, repeat, threads
    integer :: rep, i, j, idx
    real(real32) :: area, overlap_w, overlap_h

    do rep = 1, repeat
      !$omp target teams distribute parallel do collapse(2) thread_limit(threads) private(idx, area, overlap_w, overlap_h)
      do i = 1, limit
        do j = 1, limit
          if (points(i)%w < points(j)%w) then
            area = (points(j)%z + 1.0_real32) * (points(j)%z + 1.0_real32)
            overlap_w = max(0.0_real32, min(points(i)%x + points(i)%z, points(j)%x + points(j)%z) - &
              max(points(i)%x, points(j)%x) + 1.0_real32)
            overlap_h = max(0.0_real32, min(points(i)%y + points(i)%z, points(j)%y + points(j)%z) - &
              max(points(i)%y, points(j)%y) + 1.0_real32)
            idx = (i - 1) * max_detections + j
            if (((overlap_w * overlap_h) / area) < 0.3_real32 .and. points(j)%z /= 0.0_real32) then
              nmsbitmap(idx) = 1_int8
            else
              nmsbitmap(idx) = 0_int8
            end if
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
  end subroutine generate_nms_bitmap

  subroutine reduce_nms_bitmap(nmsbitmap, pointsbitmap, ndetections, repeat)
    integer(int8), intent(in) :: nmsbitmap(:)
    integer(int8), intent(inout) :: pointsbitmap(:)
    integer, intent(in) :: ndetections, repeat
    integer :: rep, bid, lid, idx, part, s

    do rep = 1, repeat
      !$omp target teams num_teams(ndetections) thread_limit(max_detections / n_partitions) private(bid, lid, idx, part, s)
      !$omp parallel private(bid, lid, idx, part)
      bid = omp_get_team_num()
      lid = omp_get_thread_num()
      idx = bid * max_detections + lid + 1

      if (lid == 0) s = 1
      !$omp barrier

      !$omp atomic update
      s = iand(s, int(nmsbitmap(idx)))
      !$omp barrier

      do part = 1, n_partitions - 1
        idx = idx + max_detections / n_partitions
        !$omp atomic update
        s = iand(s, int(nmsbitmap(idx)))
        !$omp barrier
      end do
      pointsbitmap(bid + 1) = int(s, int8)
      !$omp end parallel
      !$omp end target teams
    end do
  end subroutine reduce_nms_bitmap

  subroutine write_points(path, points, pointsbitmap, ndetections, totaldets)
    character(len=*), intent(in) :: path
    type(float4), intent(in) :: points(:)
    integer(int8), intent(in) :: pointsbitmap(:)
    integer, intent(in) :: ndetections
    integer, intent(out) :: totaldets
    integer :: unit, ios, i

    open(newunit=unit, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print '(A,A,A)', 'Error: Unable to open file ', path, ' for detection outcome.'
      stop 1
    end if
    totaldets = 0
    do i = 1, ndetections
      if (pointsbitmap(i) /= 0_int8) then
        write(unit,'(I0,A,I0,A,I0,A,F0.6)') int(points(i)%x), ',', int(points(i)%y), ',', int(points(i)%z), ',', points(i)%w
        totaldets = totaldets + 1
      end if
    end do
    close(unit)
  end subroutine write_points

end program main
