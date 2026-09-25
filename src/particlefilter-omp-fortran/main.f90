! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real32, real64
  use omp_lib
  implicit none

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 256
#endif

  integer(int32), parameter :: block_size = BLOCK_SIZE
  real(real32), parameter :: pi = 3.1415926535897932_real32
  integer(int64), parameter :: lcg_a = 1103515245_int64
  integer(int64), parameter :: lcg_c = 12345_int64
  integer(int64), parameter :: lcg_m = 2147483647_int64
  real(real32), parameter :: scale_factor = 300.0_real32

  character(len=256) :: arg
  integer(int32) :: iszx, iszy, nfr, nparticles
  integer(int32), allocatable :: seed(:)
  integer(int8), allocatable :: image(:)
  real(real64) :: start_time, end_video_sequence, end_particle_filter

  if (command_argument_count() /= 8) then
    call print_usage()
    stop 0
  end if

  call get_command_argument(1, arg)
  if (trim(arg) /= '-x') then
    call print_usage()
    stop 0
  end if
  call get_command_argument(3, arg)
  if (trim(arg) /= '-y') then
    call print_usage()
    stop 0
  end if
  call get_command_argument(5, arg)
  if (trim(arg) /= '-z') then
    call print_usage()
    stop 0
  end if
  call get_command_argument(7, arg)
  if (trim(arg) /= '-np') then
    call print_usage()
    stop 0
  end if

  iszx = read_int_arg(2)
  if (iszx <= 0_int32) then
    write(*,'(A)') 'dimX must be > 0'
    stop 0
  end if
  iszy = read_int_arg(4)
  if (iszy <= 0_int32) then
    write(*,'(A)') 'dimY must be > 0'
    stop 0
  end if
  nfr = read_int_arg(6)
  if (nfr <= 0_int32) then
    write(*,'(A)') 'number of frames must be > 0'
    stop 0
  end if
  nparticles = read_int_arg(8)
  if (nparticles <= 0_int32) then
    write(*,'(A)') 'Number of particles must be > 0'
    stop 0
  end if

#ifdef DEBUG
  write(*,'(A,I0,A,I0,A,I0,A,I0)') 'dimX=', iszx, ' dimY=', iszy, ' Nfr=', nfr, ' Nparticles=', nparticles
#endif

  allocate(seed(0:nparticles - 1))
  call initialize_seed(seed)

  allocate(image(0:iszx * iszy * nfr - 1))
  image = to_byte(0_int32)

  start_time = wall_time_seconds()
  call video_sequence(image, iszx, iszy, nfr, seed)
  end_video_sequence = wall_time_seconds()
  write(*,'(A,F8.6,A)') 'VIDEO SEQUENCE TOOK ', real(end_video_sequence - start_time, real64), ' (s)'

  call particle_filter(image, iszx, iszy, nfr, seed, nparticles)
  end_particle_filter = wall_time_seconds()
  write(*,'(A,F8.6,A)') 'PARTICLE FILTER TOOK ', real(end_particle_filter - end_video_sequence, real64), ' (s)'
  write(*,'(A,F8.6,A)') 'ENTIRE PROGRAM TOOK ', real(end_particle_filter - start_time, real64), ' (s)'

  deallocate(seed, image)

contains

  subroutine print_usage()
    write(*,'(A)') './main -x <dimX> -y <dimY> -z <Nfr> -np <Nparticles>'
  end subroutine print_usage

  integer(int32) function read_int_arg(position) result(value)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) value
  end function read_int_arg

  real(real64) function wall_time_seconds() result(value)
    integer(int64) :: count, rate

    call system_clock(count, rate)
    value = real(count, real64) / real(rate, real64)
  end function wall_time_seconds

  subroutine initialize_seed(seed)
    integer(int32), intent(out) :: seed(0:)
    integer :: i

    do i = 0, size(seed) - 1
      seed(i) = int(i + 1, int32)
    end do
  end subroutine initialize_seed

  real(real32) function randu(seed, index) result(value)
    integer(int32), intent(inout) :: seed(0:)
    integer, intent(in) :: index

    seed(index) = lcg_next(seed(index))
    value = abs(real(seed(index), real32) / real(lcg_m, real32))
  end function randu

  real(real32) function randn(seed, index) result(value)
    integer(int32), intent(inout) :: seed(0:)
    integer, intent(in) :: index
    real(real32) :: u, v

    u = randu(seed, index)
    v = randu(seed, index)
    if (u <= 0.0_real32) u = tiny(u)
    value = sqrt(-2.0_real32 * log(u)) * cos(2.0_real32 * pi * v)
  end function randn

  integer(int32) function lcg_next(old_seed) result(new_seed)
    integer(int32), intent(in) :: old_seed
    integer(int64) :: raw, quotient, remainder

    raw = modulo(lcg_a * int(old_seed, int64) + lcg_c, 4294967296_int64)
    if (raw >= 2147483648_int64) raw = raw - 4294967296_int64
    quotient = raw / lcg_m
    remainder = raw - quotient * lcg_m
    new_seed = int(remainder, int32)
  end function lcg_next

  integer(int32) function round_float(value) result(new_value)
    real(real32), intent(in) :: value
    integer(int32) :: truncated

    truncated = int(value, int32)
    if (value - real(truncated, real32) < 0.5_real32) then
      new_value = truncated
    else
      new_value = truncated + 1_int32
    end if
  end function round_float

  integer(int32) function to_u8_from_real(value) result(byte)
    real(real32), intent(in) :: value
    integer(int32) :: truncated

    truncated = int(value, int32)
    byte = modulo(truncated, 256_int32)
  end function to_u8_from_real

  integer(int8) function to_byte(value) result(byte)
    integer(int32), intent(in) :: value
    integer(int32) :: wrapped

    wrapped = modulo(value, 256_int32)
    if (wrapped > 127_int32) wrapped = wrapped - 256_int32
    byte = int(wrapped, int8)
  end function to_byte

  integer(int32) function byte_as_int(byte) result(value)
    integer(int8), intent(in) :: byte

    value = int(byte, int32)
    if (value < 0_int32) value = value + 256_int32
  end function byte_as_int

  subroutine set_if(test_value, new_value, array3d, dimx, dimy, dimz)
    integer(int32), intent(in) :: test_value, new_value, dimx, dimy, dimz
    integer(int8), intent(inout) :: array3d(0:)
    integer(int32) :: x, y, z, idx

    do x = 0, dimx - 1
      do y = 0, dimy - 1
        do z = 0, dimz - 1
          idx = x * dimy * dimz + y * dimz + z
          if (byte_as_int(array3d(idx)) == test_value) array3d(idx) = to_byte(new_value)
        end do
      end do
    end do
  end subroutine set_if

  subroutine add_noise(array3d, dimx, dimy, dimz, seed)
    integer(int32), intent(in) :: dimx, dimy, dimz
    integer(int8), intent(inout) :: array3d(0:)
    integer(int32), intent(inout) :: seed(0:)
    integer(int32) :: x, y, z, idx, noise

    do x = 0, dimx - 1
      do y = 0, dimy - 1
        do z = 0, dimz - 1
          idx = x * dimy * dimz + y * dimz + z
          noise = to_u8_from_real(5.0_real32 * randn(seed, 0))
          array3d(idx) = to_byte(byte_as_int(array3d(idx)) + noise)
        end do
      end do
    end do
  end subroutine add_noise

  subroutine strel_disk(disk, radius)
    integer(int32), intent(out) :: disk(0:)
    integer(int32), intent(in) :: radius
    integer(int32) :: diameter, x, y
    real(real32) :: distance

    diameter = radius * 2_int32 - 1_int32
    do x = 0, diameter - 1
      do y = 0, diameter - 1
        distance = sqrt(real((x - radius + 1_int32) * (x - radius + 1_int32) + &
          (y - radius + 1_int32) * (y - radius + 1_int32), real32))
        if (distance < real(radius, real32)) then
          disk(x * diameter + y) = 1_int32
        else
          disk(x * diameter + y) = 0_int32
        end if
      end do
    end do
  end subroutine strel_disk

  subroutine dilate_matrix(matrix, posx, posy, posz, dimx, dimy, dimz, error_radius)
    integer(int8), intent(inout) :: matrix(0:)
    integer(int32), intent(in) :: posx, posy, posz, dimx, dimy, dimz, error_radius
    integer(int32) :: startx, starty, endx, endy, x, y
    real(real32) :: distance

    startx = posx - error_radius
    do while (startx < 0_int32)
      startx = startx + 1_int32
    end do
    starty = posy - error_radius
    do while (starty < 0_int32)
      starty = starty + 1_int32
    end do
    endx = posx + error_radius
    do while (endx > dimx)
      endx = endx - 1_int32
    end do
    endy = posy + error_radius
    do while (endy > dimy)
      endy = endy - 1_int32
    end do

    do x = startx, endx - 1
      do y = starty, endy - 1
        distance = sqrt(real((x - posx) * (x - posx) + (y - posy) * (y - posy), real32))
        if (distance < real(error_radius, real32)) matrix(x * dimy * dimz + y * dimz + posz) = to_byte(1_int32)
      end do
    end do
  end subroutine dilate_matrix

  subroutine imdilate_disk(matrix, dimx, dimy, dimz, error_radius, new_matrix)
    integer(int8), intent(in) :: matrix(0:)
    integer(int8), intent(inout) :: new_matrix(0:)
    integer(int32), intent(in) :: dimx, dimy, dimz, error_radius
    integer(int32) :: x, y, z

    do z = 0, dimz - 1
      do x = 0, dimx - 1
        do y = 0, dimy - 1
          if (byte_as_int(matrix(x * dimy * dimz + y * dimz + z)) == 1_int32) then
            call dilate_matrix(new_matrix, x, y, z, dimx, dimy, dimz, error_radius)
          end if
        end do
      end do
    end do
  end subroutine imdilate_disk

  subroutine getneighbors(se, num_ones, objxy, radius)
    integer(int32), intent(in) :: se(0:), num_ones, radius
    integer(int32), intent(out) :: objxy(0:)
    integer(int32) :: x, y, neighy, center, diameter

    neighy = 0_int32
    center = radius - 1_int32
    diameter = radius * 2_int32 - 1_int32
    do x = 0, diameter - 1
      do y = 0, diameter - 1
        if (se(x * diameter + y) /= 0_int32) then
          objxy(neighy * 2_int32) = y - center
          objxy(neighy * 2_int32 + 1_int32) = x - center
          neighy = neighy + 1_int32
        end if
      end do
    end do
  end subroutine getneighbors

  subroutine video_sequence(image, iszx, iszy, nfr, seed)
    integer(int32), intent(in) :: iszx, iszy, nfr
    integer(int8), intent(inout) :: image(0:)
    integer(int32), intent(inout) :: seed(0:)
    integer(int8), allocatable :: new_matrix(:)
    integer(int32) :: k, max_size, x0, y0, xk, yk, pos, x, y

    max_size = iszx * iszy * nfr
    x0 = round_float(real(iszy, real32) / 2.0_real32)
    y0 = round_float(real(iszx, real32) / 2.0_real32)
    image(x0 * iszy * nfr + y0 * nfr) = to_byte(1_int32)

    do k = 1, nfr - 1
      xk = abs(x0 + (k - 1_int32))
      yk = abs(y0 - 2_int32 * (k - 1_int32))
      pos = yk * iszy * nfr + xk * nfr + k
      if (pos >= max_size) pos = 0_int32
      image(pos) = to_byte(1_int32)
    end do

    allocate(new_matrix(0:max_size - 1))
    new_matrix = to_byte(0_int32)
    call imdilate_disk(image, iszx, iszy, nfr, 5_int32, new_matrix)
    do x = 0, iszx - 1
      do y = 0, iszy - 1
        do k = 0, nfr - 1
          image(x * iszy * nfr + y * nfr + k) = new_matrix(x * iszy * nfr + y * nfr + k)
        end do
      end do
    end do
    deallocate(new_matrix)

    call set_if(0_int32, 100_int32, image, iszx, iszy, nfr)
    call set_if(1_int32, 228_int32, image, iszx, iszy, nfr)
    call add_noise(image, iszx, iszy, nfr, seed)
  end subroutine video_sequence

  subroutine particle_filter(image, iszx, iszy, nfr, seed, nparticles)
    integer(int8), intent(in) :: image(0:)
    integer(int32), intent(in) :: iszx, iszy, nfr, nparticles
    integer(int32), intent(inout) :: seed(0:)
    integer(int32), allocatable :: disk(:), objxy(:), ind(:)
    real(real32), allocatable :: weights(:), likelihood(:), partial_sums(:), arrayx(:), arrayy(:)
    real(real32), allocatable :: xj(:), yj(:), cdf(:), u(:)
    integer(int32) :: radius, diameter, count_ones, x, y, k, num_blocks, max_size
    real(real32) :: xe, ye, distance
    real(real64) :: offload_start, start_time, end_time, offload_end

    max_size = iszx * iszy * nfr
    xe = real(round_float(real(iszy, real32) / 2.0_real32), real32)
    ye = real(round_float(real(iszx, real32) / 2.0_real32), real32)

    radius = 5_int32
    diameter = radius * 2_int32 - 1_int32
    allocate(disk(0:diameter * diameter - 1))
    call strel_disk(disk, radius)
    count_ones = 0_int32
    do x = 0, diameter - 1
      do y = 0, diameter - 1
        if (disk(x * diameter + y) == 1_int32) count_ones = count_ones + 1_int32
      end do
    end do

    allocate(objxy(0:2 * count_ones - 1))
    call getneighbors(disk, count_ones, objxy, radius)

    allocate(weights(0:nparticles - 1), likelihood(0:nparticles), partial_sums(0:nparticles))
    allocate(arrayx(0:nparticles - 1), arrayy(0:nparticles - 1))
    allocate(xj(0:nparticles - 1), yj(0:nparticles - 1), cdf(0:nparticles - 1))
    allocate(ind(0:count_ones * nparticles - 1), u(0:nparticles - 1))

    weights = 1.0_real32 / real(nparticles, real32)
    do x = 0, nparticles - 1
      xj(x) = xe
      yj(x) = ye
    end do

    num_blocks = (nparticles + block_size - 1_int32) / block_size
    offload_start = wall_time_seconds()

    !$omp target data map(alloc: likelihood(0:nparticles), ind(0:count_ones*nparticles-1), &
    !$omp& u(0:nparticles-1), partial_sums(0:nparticles), cdf(0:nparticles-1)) &
    !$omp& map(from: arrayx(0:nparticles-1), arrayy(0:nparticles-1)) &
    !$omp& map(tofrom: weights(0:nparticles-1)) &
    !$omp& map(to: xj(0:nparticles-1), yj(0:nparticles-1), seed(0:nparticles-1), &
    !$omp& image(0:iszx*iszy*nfr-1), objxy(0:2*count_ones-1))

    start_time = wall_time_seconds()
    do k = 1, nfr - 1
      call likelihood_kernel(arrayx, arrayy, xj, yj, ind, objxy, likelihood, image, &
        weights, seed, partial_sums, nparticles, count_ones, iszy, nfr, k, max_size, num_blocks)
      call sum_kernel(partial_sums, nparticles, num_blocks)
      call normalize_kernel(weights, partial_sums, cdf, u, seed, nparticles, num_blocks)
      call find_index_kernel(arrayx, arrayy, cdf, u, xj, yj, nparticles)
    end do
    end_time = wall_time_seconds()
    if (nfr > 1_int32) then
      write(*,'(A,F8.6,A)') 'Average execution time of kernels: ', &
        real((end_time - start_time) / real(nfr - 1_int32, real64), real64), ' (s)'
    else
      write(*,'(A,F8.6,A)') 'Average execution time of kernels: ', 0.0_real64, ' (s)'
    end if

    !$omp end target data

    offload_end = wall_time_seconds()
    write(*,'(A,F8.6,A)') 'Device offloading time: ', real(offload_end - offload_start, real64), ' (s)'

    xe = 0.0_real32
    ye = 0.0_real32
    do x = 0, nparticles - 1
      xe = xe + arrayx(x) * weights(x)
      ye = ye + arrayy(x) * weights(x)
    end do
    distance = sqrt((xe - real(round_float(real(iszy, real32) / 2.0_real32), real32)) ** 2 + &
      (ye - real(round_float(real(iszx, real32) / 2.0_real32), real32)) ** 2)

    open(unit=10, file='output.txt', status='replace', action='write')
    write(10,'(A,F0.6)') 'XE: ', real(xe, real64)
    write(10,'(A,F0.6)') 'YE: ', real(ye, real64)
    write(10,'(A,F0.6)') 'distance: ', real(distance, real64)
    close(10)

    deallocate(disk, objxy, weights, likelihood, partial_sums, arrayx, arrayy, xj, yj, cdf, ind, u)
  end subroutine particle_filter

  subroutine likelihood_kernel(arrayx, arrayy, xj, yj, ind, objxy, likelihood, image, &
      weights, seed, partial_sums, nparticles, count_ones, iszy, nfr, frame, max_size, num_blocks)
    integer(int32), intent(in) :: nparticles, count_ones, iszy, nfr, frame, max_size, num_blocks
    real(real32), intent(inout) :: arrayx(0:), arrayy(0:), weights(0:), likelihood(0:), partial_sums(0:)
    real(real32), intent(in) :: xj(0:), yj(0:)
    integer(int32), intent(inout) :: ind(0:), seed(0:)
    integer(int32), intent(in) :: objxy(0:)
    integer(int8), intent(in) :: image(0:)
    integer(int32) :: i, y, ix, iy, rnd_ix, rnd_iy, indx, indy
    integer(int32) :: team_id, thread_id, block_dim, s
    real(real32) :: ur, vr, likelihood_sum

    !$omp target teams num_teams(num_blocks) thread_limit(block_size)
    block
      real(real32) :: weights_local(block_size)

      !$omp parallel private(team_id, thread_id, block_dim, i, y, ix, iy, rnd_ix, rnd_iy, indx, indy, ur, vr, likelihood_sum, s)
      team_id = omp_get_team_num()
      thread_id = omp_get_thread_num()
      block_dim = omp_get_num_threads()
      i = team_id * block_dim + thread_id

      if (i < nparticles) then
        arrayx(i) = xj(i)
        arrayy(i) = yj(i)
        weights(i) = 1.0_real32 / real(nparticles, real32)

        seed(i) = lcg_next(seed(i))
        ur = abs(real(seed(i), real32) / real(lcg_m, real32))
        seed(i) = lcg_next(seed(i))
        vr = abs(real(seed(i), real32) / real(lcg_m, real32))
        arrayx(i) = arrayx(i) + 1.0_real32 + 5.0_real32 * (sqrt(-2.0_real32 * log(ur)) * cos(2.0_real32 * pi * vr))

        seed(i) = lcg_next(seed(i))
        ur = abs(real(seed(i), real32) / real(lcg_m, real32))
        seed(i) = lcg_next(seed(i))
        vr = abs(real(seed(i), real32) / real(lcg_m, real32))
        arrayy(i) = arrayy(i) - 2.0_real32 + 2.0_real32 * (sqrt(-2.0_real32 * log(ur)) * cos(2.0_real32 * pi * vr))
      end if

      !$omp barrier

      if (i < nparticles) then
        ix = int(arrayx(i), int32)
        iy = int(arrayy(i), int32)
        if (arrayx(i) - real(ix, real32) < 0.5_real32) then
          rnd_ix = ix
        else
          rnd_ix = ix
        end if
        if (arrayy(i) - real(iy, real32) < 0.5_real32) then
          rnd_iy = iy
        else
          rnd_iy = iy
        end if

        do y = 0, count_ones - 1
          indx = rnd_ix + objxy(y * 2 + 1)
          indy = rnd_iy + objxy(y * 2)
          ind(i * count_ones + y) = abs(indx * iszy * nfr + indy * nfr + frame)
          if (ind(i * count_ones + y) >= max_size) ind(i * count_ones + y) = 0_int32
        end do

        likelihood_sum = 0.0_real32
        do y = 0, count_ones - 1
          likelihood_sum = likelihood_sum + &
            real((byte_as_int(image(ind(i * count_ones + y))) - 100_int32) * &
                 (byte_as_int(image(ind(i * count_ones + y))) - 100_int32) - &
                 (byte_as_int(image(ind(i * count_ones + y))) - 228_int32) * &
                 (byte_as_int(image(ind(i * count_ones + y))) - 228_int32), real32) / 50.0_real32
        end do
        likelihood(i) = likelihood_sum / real(count_ones, real32) - scale_factor
        weights(i) = weights(i) * exp(likelihood(i))
      end if

      weights_local(thread_id + 1) = 0.0_real32
      if (i < nparticles) weights_local(thread_id + 1) = weights(i)

      !$omp barrier

      s = block_dim / 2_int32
      do while (s > 0_int32)
        if (thread_id < s) weights_local(thread_id + 1) = weights_local(thread_id + 1) + weights_local(thread_id + s + 1)
        !$omp barrier
        s = s / 2_int32
      end do

      if (thread_id == 0_int32) partial_sums(team_id) = weights_local(1)
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine likelihood_kernel

  subroutine sum_kernel(partial_sums, nparticles, num_blocks)
    real(real32), intent(inout) :: partial_sums(0:)
    integer(int32), intent(in) :: nparticles, num_blocks
    integer(int32) :: x
    real(real32) :: sum_value

    !$omp target
    sum_value = 0.0_real32
    do x = 0, num_blocks - 1
      sum_value = sum_value + partial_sums(x)
    end do
    partial_sums(0) = sum_value
    !$omp end target
  end subroutine sum_kernel

  subroutine normalize_kernel(weights, partial_sums, cdf, u, seed, nparticles, num_blocks)
    real(real32), intent(inout) :: weights(0:), cdf(0:), u(0:)
    real(real32), intent(in) :: partial_sums(0:)
    integer(int32), intent(inout) :: seed(0:)
    integer(int32), intent(in) :: nparticles, num_blocks
    integer(int32) :: i, local_id, x
    real(real32) :: sum_weights, p, q, u1

    !$omp target teams num_teams(num_blocks) thread_limit(block_size)
      !$omp parallel private(local_id, i, x, p, q)
      local_id = omp_get_thread_num()
      i = omp_get_team_num() * omp_get_num_threads() + local_id

      if (local_id == 0_int32) sum_weights = partial_sums(0)

      !$omp barrier

      if (i < nparticles) weights(i) = weights(i) / sum_weights

      !$omp barrier

      if (i == 0_int32) then
        cdf(0) = weights(0)
        do x = 1, nparticles - 1
          cdf(x) = weights(x) + cdf(x - 1)
        end do

        seed(0) = lcg_next(seed(0))
        p = abs(real(seed(0), real32) / real(lcg_m, real32))
        seed(0) = lcg_next(seed(0))
        q = abs(real(seed(0), real32) / real(lcg_m, real32))
        u(0) = (1.0_real32 / real(nparticles, real32)) * (sqrt(-2.0_real32 * log(p)) * cos(2.0_real32 * pi * q))
      end if

      !$omp barrier

      if (local_id == 0_int32) u1 = u(0)

      !$omp barrier

      if (i < nparticles) u(i) = u1 + real(i, real32) / real(nparticles, real32)
      !$omp end parallel
    !$omp end target teams
  end subroutine normalize_kernel

  subroutine find_index_kernel(arrayx, arrayy, cdf, u, xj, yj, nparticles)
    real(real32), intent(in) :: arrayx(0:), arrayy(0:), cdf(0:), u(0:)
    real(real32), intent(inout) :: xj(0:), yj(0:)
    integer(int32), intent(in) :: nparticles
    integer(int32) :: i, x, index

    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, nparticles - 1
      index = -1_int32
      do x = 0, nparticles - 1
        if (cdf(x) >= u(i)) then
          index = x
          exit
        end if
      end do
      if (index == -1_int32) index = nparticles - 1
      xj(i) = arrayx(index)
      yj(i) = arrayy(index)
    end do
    !$omp end target teams distribute parallel do
  end subroutine find_index_kernel

end program main
