! SPDX-License-Identifier: CC0-1.0
program kmeans_omp_fortran
  use iso_fortran_env, only: real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: block_size2 = 256

  character(len=:), allocatable :: filename
  logical :: is_binary, is_rmse, is_output
  real(real32) :: threshold
  integer :: min_nclusters, max_nclusters, nloops
  real(real32), allocatable :: features(:)
  real(real32), allocatable :: cluster_centres(:)
  integer :: npoints, nfeatures, best_nclusters, index
  real(real32) :: rmse
  real(real64) :: start_time, elapsed_ms
  integer :: i, j
  character(len=64) :: coord_text

  print '(A,I0,A,I0,A)', "WG size of kernel_swap = ", block_size, &
      ", WG size of kernel_kmeans = ", block_size2, " "

  call parse_args(filename, is_binary, threshold, min_nclusters, max_nclusters, &
      nloops, is_rmse, is_output)
  call read_features(filename, is_binary, features, npoints, nfeatures)

  print *
  print '(A)', "I/O completed"
  print *
  print '(A,I0)', "Number of objects: ", npoints
  print '(A,I0)', "Number of features: ", nfeatures

  if (npoints < min_nclusters) then
    print '(A,I0,A,I0,A)', "Error: min_nclusters(", min_nclusters, &
        ") > npoints(", npoints, ") -- cannot proceed"
    stop 0
  end if

  start_time = omp_get_wtime()
  call cluster(npoints, nfeatures, features, min_nclusters, max_nclusters, &
      threshold, best_nclusters, cluster_centres, rmse, is_rmse, nloops, index)
  elapsed_ms = (omp_get_wtime() - start_time) * 1000.0_real64
  print '(A,I0,A)', "Kmeans core timing: ", int(elapsed_ms), " ms"

  if ((min_nclusters == max_nclusters) .and. is_output) then
    print *
    print '(A)', "================= Centroid Coordinates ================="
    do i = 1, max_nclusters
      write (*,'(I0,A)', advance='no') i - 1, ":"
      do j = 1, nfeatures
        call format_fixed2(cluster_centres((i - 1) * nfeatures + j), coord_text)
        write (*,'(A,A)', advance='no') " ", trim(coord_text)
      end do
      print *
      print *
    end do
  end if

  print '(A,I0)', "Number of Iteration: ", nloops

  if (min_nclusters /= max_nclusters) then
    print '(A,I0)', "Best number of clusters is ", best_nclusters
  else
    if (nloops /= 1) then
      if (is_rmse) then
        print '(A,F0.3,A,I0)', "Number of trials to approach the best RMSE of ", &
            rmse, " is ", index + 1
      end if
    else
      if (is_rmse) then
        print '(A,F0.3)', "Root Mean Squared Error: ", rmse
      end if
    end if
  end if

contains

  subroutine parse_args(filename, is_binary, threshold, min_nclusters, max_nclusters, &
      nloops, is_rmse, is_output)
    character(len=:), allocatable, intent(out) :: filename
    logical, intent(out) :: is_binary, is_rmse, is_output
    real(real32), intent(out) :: threshold
    integer, intent(out) :: min_nclusters, max_nclusters, nloops

    character(len=4096) :: arg, value
    integer :: argc, pos

    threshold = 0.001_real32
    max_nclusters = 5
    min_nclusters = 5
    nloops = 1
    is_binary = .false.
    is_rmse = .false.
    is_output = .false.
    filename = ""

    argc = command_argument_count()
    pos = 1
    do while (pos <= argc)
      call get_command_argument(pos, arg)
      select case (trim(arg))
      case ("-i")
        call require_value(pos, argc, value)
        filename = trim(value)
        pos = pos + 2
      case ("-b")
        is_binary = .true.
        pos = pos + 1
      case ("-t")
        call require_value(pos, argc, value)
        read (value, *) threshold
        pos = pos + 2
      case ("-m")
        call require_value(pos, argc, value)
        read (value, *) max_nclusters
        pos = pos + 2
      case ("-n")
        call require_value(pos, argc, value)
        read (value, *) min_nclusters
        pos = pos + 2
      case ("-l")
        call require_value(pos, argc, value)
        read (value, *) nloops
        pos = pos + 2
      case ("-r")
        is_rmse = .true.
        pos = pos + 1
      case ("-o")
        is_output = .true.
        pos = pos + 1
      case default
        call usage()
      end select
    end do

    if (len_trim(filename) == 0) call usage()
  end subroutine parse_args

  subroutine require_value(pos, argc, value)
    integer, intent(in) :: pos, argc
    character(len=*), intent(out) :: value

    if (pos + 1 > argc) call usage()
    call get_command_argument(pos + 1, value)
  end subroutine require_value

  subroutine usage()
    character(len=4096) :: argv0

    call get_command_argument(0, argv0)
    write (*,'(A)') ""
    write (*,'(A,A)') "Usage: ", trim(argv0) // " [switches] -i filename"
    write (*,'(A)') ""
    write (*,'(A)') "    -i filename      :file containing data to be clustered"
    write (*,'(A)') "    -m max_nclusters :maximum number of clusters allowed    [default=5]"
    write (*,'(A)') "    -n min_nclusters :minimum number of clusters allowed    [default=5]"
    write (*,'(A)') "    -t threshold     :threshold value                       [default=0.001]"
    write (*,'(A)') "    -l nloops        :iteration for each number of clusters [default=1]"
    write (*,'(A)') "    -b               :input file is in binary format"
    write (*,'(A)') "    -r               :calculate RMSE                        [default=off]"
    write (*,'(A)') "    -o               :output cluster center coordinates     [default=off]"
    stop 1
  end subroutine usage

  subroutine format_fixed2(value, text)
    real(real32), intent(in) :: value
    character(len=*), intent(out) :: text
    character(len=64) :: raw

    write (raw, '(F0.2)') value
    raw = adjustl(raw)
    if (raw(1:1) == ".") then
      text = "0" // trim(raw)
    else if (len_trim(raw) >= 2 .and. raw(1:2) == "-.") then
      text = "-0" // trim(raw(2:))
    else
      text = trim(raw)
    end if
  end subroutine format_fixed2

  subroutine read_features(filename, is_binary, features, npoints, nfeatures)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: is_binary
    real(real32), allocatable, intent(out) :: features(:)
    integer, intent(out) :: npoints, nfeatures

    if (is_binary) then
      call read_binary_features(filename, features, npoints, nfeatures)
    else
      call read_ascii_features(filename, features, npoints, nfeatures)
    end if
  end subroutine read_features

  subroutine read_binary_features(filename, features, npoints, nfeatures)
    character(len=*), intent(in) :: filename
    real(real32), allocatable, intent(out) :: features(:)
    integer, intent(out) :: npoints, nfeatures

    integer :: unit, ios

    open (newunit=unit, file=filename, access="stream", form="unformatted", &
        status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write (*,'(A,A,A)') "Error: no such file (", trim(filename), ")"
      stop 1
    end if
    read (unit) npoints
    read (unit) nfeatures
    allocate(features(npoints * nfeatures))
    read (unit) features
    close (unit)
  end subroutine read_binary_features

  subroutine read_ascii_features(filename, features, npoints, nfeatures)
    character(len=*), intent(in) :: filename
    real(real32), allocatable, intent(out) :: features(:)
    integer, intent(out) :: npoints, nfeatures

    character(len=16384) :: line
    integer :: unit, ios, count, row, col
    real(real32), allocatable :: values(:)

    open (newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write (*,'(A,A,A)') "Error: no such file (", trim(filename), ")"
      stop 1
    end if

    npoints = 0
    nfeatures = 0
    do
      read (unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      count = count_fields(line)
      if (count > 0) then
        npoints = npoints + 1
        if (nfeatures == 0) nfeatures = count - 1
      end if
    end do

    if (npoints <= 0 .or. nfeatures <= 0) then
      write (*,'(A,A,A)') "Error: empty or invalid file (", trim(filename), ")"
      stop 1
    end if

    allocate(features(npoints * nfeatures))
    allocate(values(nfeatures + 1))
    rewind(unit)
    row = 0
    do
      read (unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      count = count_fields(line)
      if (count == 0) cycle
      read (line, *) values(1:count)
      row = row + 1
      do col = 1, nfeatures
        features((row - 1) * nfeatures + col) = values(col + 1)
      end do
    end do

    close(unit)
  end subroutine read_ascii_features

  integer function count_fields(line)
    character(len=*), intent(in) :: line
    integer :: i, n
    logical :: in_field

    n = 0
    in_field = .false.
    do i = 1, len_trim(line)
      select case (line(i:i))
      case (" ", achar(9), ",")
        in_field = .false.
      case default
        if (.not. in_field) then
          n = n + 1
          in_field = .true.
        end if
      end select
    end do
    count_fields = n
  end function count_fields

  subroutine cluster(npoints, nfeatures, features, min_nclusters, max_nclusters, &
      threshold, best_nclusters, cluster_centres, min_rmse, is_rmse, nloops, return_index)
    integer, intent(in) :: npoints, nfeatures, min_nclusters, max_nclusters, nloops
    real(real32), intent(in) :: features(:), threshold
    integer, intent(out) :: best_nclusters, return_index
    real(real32), allocatable, intent(out) :: cluster_centres(:)
    real(real32), intent(out) :: min_rmse
    logical, intent(in) :: is_rmse

    integer, allocatable :: membership(:), membership_device(:)
    real(real32), allocatable :: feature_swap(:), clusters(:)
    real(real32), allocatable :: best_clusters(:)
    integer :: nclusters, lp, i, j, point_id, loop_count
    integer :: selected, initial_points, n, temp
    integer, allocatable :: initial(:), new_centers_len(:)
    real(real32), allocatable :: new_centers(:)
    real(real32) :: delta, trial_rmse
    integer :: trial_rmse_int

    allocate(membership(npoints), membership_device(npoints))
    allocate(feature_swap(npoints * nfeatures))
    min_rmse = huge(1.0_real32)
    best_nclusters = 0
    return_index = 0

    !$omp target data map(to: features(1:npoints*nfeatures)) &
    !$omp& map(alloc: feature_swap(1:npoints*nfeatures), membership_device(1:npoints))
    do nclusters = min_nclusters, max_nclusters
      if (nclusters > npoints) exit

      !$omp target teams distribute parallel do thread_limit(block_size) private(j)
      do point_id = 1, npoints
        do j = 1, nfeatures
          feature_swap((j - 1) * npoints + point_id) = features((point_id - 1) * nfeatures + j)
        end do
      end do
      !$omp end target teams distribute parallel do

      allocate(clusters(nclusters * nfeatures))
      allocate(initial(npoints))
      do i = 1, npoints
        initial(i) = i
      end do
      initial_points = npoints

      do lp = 1, nloops
        n = 1
        do i = 1, nclusters
          if (initial_points < 1) exit
          selected = initial(n)
          do j = 1, nfeatures
            clusters((i - 1) * nfeatures + j) = features((selected - 1) * nfeatures + j)
          end do
          temp = initial(n)
          initial(n) = initial(initial_points)
          initial(initial_points) = temp
          initial_points = initial_points - 1
          n = n + 1
        end do

        membership = -1
        allocate(new_centers_len(nclusters), new_centers(nclusters * nfeatures))
        new_centers_len = 0
        new_centers = 0.0_real32
        loop_count = 0

        do
          delta = 0.0_real32
          !$omp target data map(to: clusters(1:nclusters*nfeatures))
          !$omp target teams distribute parallel do thread_limit(block_size2)
          do point_id = 1, npoints
            membership_device(point_id) = nearest_cluster_transposed( &
                point_id, npoints, nfeatures, nclusters, feature_swap, clusters)
          end do
          !$omp end target teams distribute parallel do
          !$omp end target data
          !$omp target update from(membership_device(1:npoints))

          do i = 1, npoints
            selected = membership_device(i)
            new_centers_len(selected) = new_centers_len(selected) + 1
            if (membership_device(i) /= membership(i)) then
              delta = delta + 1.0_real32
              membership(i) = membership_device(i)
            end if
            do j = 1, nfeatures
              new_centers((selected - 1) * nfeatures + j) = &
                  new_centers((selected - 1) * nfeatures + j) + &
                  features((i - 1) * nfeatures + j)
            end do
          end do

          do i = 1, nclusters
            do j = 1, nfeatures
              if (new_centers_len(i) > 0) then
                clusters((i - 1) * nfeatures + j) = &
                    new_centers((i - 1) * nfeatures + j) / real(new_centers_len(i), real32)
              end if
              new_centers((i - 1) * nfeatures + j) = 0.0_real32
            end do
            new_centers_len(i) = 0
          end do

          if (.not. (delta > threshold .and. loop_count < 500)) exit
          loop_count = loop_count + 1
        end do

        deallocate(new_centers, new_centers_len)

        if (is_rmse) then
          trial_rmse = rms_err(features, nfeatures, npoints, clusters, nclusters)
          trial_rmse_int = int(trial_rmse)
          if (real(trial_rmse_int, real32) < min_rmse) then
            min_rmse = real(trial_rmse_int, real32)
            best_nclusters = nclusters
            return_index = lp - 1
          end if
        end if
      end do

      if (allocated(best_clusters)) deallocate(best_clusters)
      allocate(best_clusters(nclusters * nfeatures))
      best_clusters = clusters
      deallocate(initial, clusters)
    end do
    !$omp end target data

    if (allocated(best_clusters)) then
      allocate(cluster_centres(size(best_clusters)))
      cluster_centres = best_clusters
      deallocate(best_clusters)
    else
      allocate(cluster_centres(0))
    end if
    deallocate(feature_swap, membership_device, membership)
  end subroutine cluster

  integer function nearest_cluster_transposed(point_id, npoints, nfeatures, nclusters, feature_swap, clusters)
    !$omp declare target
    integer, intent(in) :: point_id, npoints, nfeatures, nclusters
    real(real32), intent(in) :: feature_swap(:), clusters(:)
    integer :: i, l
    real(real32) :: min_dist, dist, diff

    min_dist = huge(1.0_real32)
    nearest_cluster_transposed = 1
    do i = 1, nclusters
      dist = 0.0_real32
      do l = 1, nfeatures
        diff = feature_swap((l - 1) * npoints + point_id) - &
            clusters((i - 1) * nfeatures + l)
        dist = dist + diff * diff
      end do
      if (dist < min_dist) then
        min_dist = dist
        nearest_cluster_transposed = i
      end if
    end do
  end function nearest_cluster_transposed

  real(real32) function rms_err(features, nfeatures, npoints, cluster_centres, nclusters)
    real(real32), intent(in) :: features(:), cluster_centres(:)
    integer, intent(in) :: nfeatures, npoints, nclusters
    integer :: i, nearest
    real(real32) :: sum_euclid

    sum_euclid = 0.0_real32
    !$omp parallel do default(shared) private(i, nearest) reduction(+:sum_euclid) schedule(static)
    do i = 1, npoints
      nearest = nearest_cluster_row(features, i, nfeatures, cluster_centres, nclusters)
      sum_euclid = sum_euclid + euclid_dist_2(features, i, cluster_centres, nearest, nfeatures)
    end do
    !$omp end parallel do
    rms_err = sqrt(sum_euclid / real(npoints, real32))
  end function rms_err

  integer function nearest_cluster_row(features, point_id, nfeatures, cluster_centres, nclusters)
    real(real32), intent(in) :: features(:), cluster_centres(:)
    integer, intent(in) :: point_id, nfeatures, nclusters
    integer :: i
    real(real32) :: min_dist, dist

    min_dist = huge(1.0_real32)
    nearest_cluster_row = 1
    do i = 1, nclusters
      dist = euclid_dist_2(features, point_id, cluster_centres, i, nfeatures)
      if (dist < min_dist) then
        min_dist = dist
        nearest_cluster_row = i
      end if
    end do
  end function nearest_cluster_row

  real(real32) function euclid_dist_2(features, point_id, cluster_centres, cluster_id, nfeatures)
    real(real32), intent(in) :: features(:), cluster_centres(:)
    integer, intent(in) :: point_id, cluster_id, nfeatures
    integer :: j
    real(real32) :: diff

    euclid_dist_2 = 0.0_real32
    do j = 1, nfeatures
      diff = features((point_id - 1) * nfeatures + j) - &
          cluster_centres((cluster_id - 1) * nfeatures + j)
      euclid_dist_2 = euclid_dist_2 + diff * diff
    end do
  end function euclid_dist_2

end program kmeans_omp_fortran
