! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : error_unit, int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: tilesize = 128
  integer, parameter :: max_threads = 256
  integer(int64), parameter :: lcg_m = 2147483648_int64
  integer(int64), parameter :: lcg_a = 26757677_int64
  integer(int64), parameter :: lcg_c = 1_int64

  character(len=512) :: input_file
  integer :: cities, restarts, repeat, threads
  integer(int32), allocatable :: buf(:, :)
  integer(int32) :: climbs(1), best(1)
  integer(int64) :: moves
  character(len=64) :: metric
  real(real32), allocatable :: posx(:), posy(:), pxglob(:, :), pyglob(:, :)
  real(real64) :: ktime, kstart, kend
  integer :: rep

  write(*, '(A)') '2-opt TSP OpenMP target offloading GPU code v2.3'
  write(*, '(A)') 'Copyright (c) 2014-2020, Texas State University. All rights reserved.'

  if (command_argument_count() /= 3) then
    write(error_unit, '(A)') ''
    write(error_unit, '(A)') 'arguments: <input_file> <restart_count> <repeat>'
    stop 1
  end if

  call get_command_argument(1, input_file)
  restarts = read_int_arg(2)
  if (restarts < 1) then
    write(error_unit, '(A,I0)') 'restart_count is too small: ', restarts
    stop 1
  end if
  repeat = read_int_arg(3)

  call read_tsp_file(trim(input_file), cities, posx, posy)

  write(*, '(A,I0,A,I0,A,A,A)') 'configuration: ', cities, ' cities, ', restarts, &
    ' restarts, ', trim(input_file), ' input'

  climbs(1) = 0_int32
  best(1) = huge(best(1))
  allocate(buf(0:cities - 1, 0:restarts - 1))
  allocate(pxglob(0:cities, 0:restarts - 1), pyglob(0:cities, 0:restarts - 1))

  threads = best_thread_count(cities)
  write(*, '(A,I0)') 'number of threads per team: ', threads
  ktime = 0.0_real64

  !$omp target data map(to: posx(0:cities - 1), posy(0:cities - 1)) &
  !$omp& map(alloc: buf(0:cities - 1, 0:restarts - 1), pxglob(0:cities, 0:restarts - 1), &
  !$omp& pyglob(0:cities, 0:restarts - 1), climbs(1:1), best(1:1))
  do rep = 0, repeat - 1
    !$omp target update to(climbs(1:1), best(1:1))

    kstart = omp_get_wtime()
    call two_opt_kernel(cities, restarts, threads, posx, posy, buf, pxglob, pyglob, climbs, best)
    kend = omp_get_wtime()

    if (rep > 0) ktime = ktime + kend - kstart
  end do

  !$omp target update from(climbs(1:1), best(1:1))
  !$omp end target data

  moves = int(climbs(1), int64) * int(cities - 2, int64) * int(cities - 1, int64) / 2_int64

  write(metric, '(F20.4)') ktime / real(repeat, real64)
  write(*, '(A,A,A)') 'Average kernel time: ', trim(adjustl(metric)), ' s'
  if (ktime > 0.0_real64) then
    write(metric, '(F20.3)') real(moves * int(repeat, int64), real64) / ktime / 1.0e9_real64
  else
    write(metric, '(F20.3)') 0.0_real64
  end if
  write(*, '(A,A)') trim(adjustl(metric)), ' Gmoves/s'
  write(*, '(A,I0,A,I0,A)') 'Best found tour length is ', best(1), ' with ', climbs(1), ' climbers'

  if (best(1) < 38000 .and. best(1) >= 35002) then
    write(*, '(A)') 'PASS'
  else
    write(*, '(A)') 'FAIL'
  end if

  deallocate(posx, posy, buf, pxglob, pyglob)

contains

  integer function read_int_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_int_arg
  end function read_int_arg

  subroutine read_tsp_file(path, cities, posx, posy)
    character(len=*), intent(in) :: path
    integer, intent(out) :: cities
    real(real32), allocatable, intent(out) :: posx(:), posy(:)

    character(len=512) :: line, token
    integer :: unit, ios, colon, idx, count
    real(real32) :: x, y
    logical :: saw_coords

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit, '(A,A)') 'could not open file ', trim(path)
      stop 1
    end if

    cities = -1
    saw_coords = .false.
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      colon = index(line, ':')
      if (colon > 0) then
        token = adjustl(line(1:colon - 1))
        if (trim(token) == 'DIMENSION') then
          read(line(colon + 1:), *) cities
          if (cities < 100) then
            write(error_unit, '(A)') 'the problem size must be at least 100 for this version of the code'
            close(unit)
            stop 1
          end if
        end if
      end if
      if (trim(adjustl(line)) == 'NODE_COORD_SECTION') then
        saw_coords = .true.
        exit
      end if
    end do

    if (.not. saw_coords .or. cities < 0) then
      write(error_unit, '(A)') 'wrong file format'
      close(unit)
      stop 1
    end if

    allocate(posx(0:cities - 1), posy(0:cities - 1))
    count = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(adjustl(line)) == 'EOF') exit
      if (len_trim(line) == 0) cycle

      read(line, *, iostat=ios) idx, x, y
      if (ios /= 0) exit

      if (count >= cities) then
        write(error_unit, '(A)') 'input too long'
      else
        posx(count) = x
        posy(count) = y
      end if
      count = count + 1
      if (count /= idx) write(error_unit, '(A,I0,A,I0)') 'input line mismatch: expected ', count, ' instead of ', idx
    end do

    if (count /= cities) write(error_unit, '(A,I0,A,I0,A)') 'read ', count, ' instead of ', cities, ' cities'

    close(unit)
  end subroutine read_tsp_file

  integer function best_thread_count(cities)
    integer, intent(in) :: cities
    integer :: max_count, best_perf, threads, smem, blocks, rounded_threads, perf

    max_count = cities - 2
    if (max_count > 256) max_count = 256
    best_perf = 0
    best_thread_count = 4
    do threads = 1, max_count
      smem = 4 * threads + 2 * 4 * tilesize + 4 * tilesize
      blocks = (16384 * 2) / smem
      if (blocks > 16) blocks = 16
      rounded_threads = ((threads + 31) / 32) * 32
      do while (blocks * rounded_threads > 2048)
        blocks = blocks - 1
      end do
      perf = threads * blocks
      if (perf > best_perf) then
        best_perf = perf
        best_thread_count = threads
      end if
    end do
  end function best_thread_count

  subroutine two_opt_kernel(cities, restarts, threads, posx, posy, buf, pxglob, pyglob, climbs, best)
    integer, intent(in) :: cities, restarts, threads
    real(real32), intent(in), contiguous :: posx(0:), posy(0:)
    integer(int32), intent(inout), contiguous :: buf(0:, 0:)
    real(real32), intent(inout), contiguous :: pxglob(0:, 0:), pyglob(0:, 0:)
    integer(int32), intent(inout) :: climbs(1), best(1)

    real(real32) :: px_s(0:tilesize - 1), py_s(0:tilesize - 1)
    integer(int32) :: bf_s(0:tilesize - 1), buf_s(0:max_threads - 1)
    integer :: lid, bid, dim, i, ii, jj, k, j, jm, bound, lower
    integer :: minchange, mini, minj, change, tmp, sumidx, term
    integer(int64) :: seed
    real(real32) :: pxi0, pyi0, pxi1, pyi1, pxj0, pyj0, pxj1, pyj1, swap_value

    !$omp target teams num_teams(restarts) thread_limit(threads) private(px_s, py_s, bf_s, buf_s)
      px_s = 0.0_real32
      py_s = 0.0_real32
      bf_s = 0_int32
      buf_s = 0_int32
      !$omp parallel private(lid, bid, dim, i, ii, jj, k, j, jm, bound, lower, minchange, mini, minj, &
      !$omp& change, tmp, sumidx, term, seed, pxi0, pyi0, pxi1, pyi1, pxj0, pyj0, pxj1, pyj1, swap_value) &
      !$omp& shared(px_s, py_s, bf_s, buf_s, posx, posy, buf, pxglob, pyglob, climbs, best)
        lid = omp_get_thread_num()
        bid = omp_get_team_num()
        dim = omp_get_num_threads()

        do i = lid, cities - 1, dim
          pxglob(i, bid) = posx(i)
          pyglob(i, bid) = posy(i)
        end do
        !$omp barrier

        if (lid == 0) then
          seed = int(bid, int64)
          do i = 1, cities - 1
            seed = modulo(lcg_a * seed + lcg_c, lcg_m)
            j = int(real(seed, real32) / real(lcg_m, real32) * real(cities - 1, real32)) + 1
            swap_value = pxglob(i, bid)
            pxglob(i, bid) = pxglob(j, bid)
            pxglob(j, bid) = swap_value
            swap_value = pyglob(i, bid)
            pyglob(i, bid) = pyglob(j, bid)
            pyglob(j, bid) = swap_value
          end do
          pxglob(cities, bid) = pxglob(0, bid)
          pyglob(cities, bid) = pyglob(0, bid)
        end if
        !$omp barrier

        do
          do i = lid, cities - 1, dim
            buf(i, bid) = -int(sqrt((pxglob(i, bid) - pxglob(i + 1, bid)) * &
              (pxglob(i, bid) - pxglob(i + 1, bid)) + &
              (pyglob(i, bid) - pyglob(i + 1, bid)) * (pyglob(i, bid) - pyglob(i + 1, bid))), int32)
          end do
          !$omp barrier

          minchange = 0
          mini = 1
          minj = 0
          do ii = 0, cities - 3, dim
            i = ii + lid
            if (i < cities - 2) then
              minchange = minchange - buf(i, bid)
              pxi0 = pxglob(i, bid)
              pyi0 = pyglob(i, bid)
              pxi1 = pxglob(i + 1, bid)
              pyi1 = pyglob(i + 1, bid)
              pxj1 = pxglob(cities, bid)
              pyj1 = pyglob(cities, bid)
            end if

            jj = cities - 1
            do while (jj >= ii + 2)
              bound = jj - tilesize + 1
              do k = lid, tilesize - 1, dim
                if (k + bound >= ii + 2) then
                  px_s(k) = pxglob(k + bound, bid)
                  py_s(k) = pyglob(k + bound, bid)
                  bf_s(k) = buf(k + bound, bid)
                end if
              end do
              !$omp barrier

              lower = bound
              if (lower < i + 2) lower = i + 2
              do j = jj, lower, -1
                jm = j - bound
                pxj0 = px_s(jm)
                pyj0 = py_s(jm)
                change = bf_s(jm) + int(sqrt((pxi0 - pxj0) * (pxi0 - pxj0) + &
                  (pyi0 - pyj0) * (pyi0 - pyj0)), int32) + &
                  int(sqrt((pxi1 - pxj1) * (pxi1 - pxj1) + &
                  (pyi1 - pyj1) * (pyi1 - pyj1)), int32)
                pxj1 = pxj0
                pyj1 = pyj0
                if (minchange > change) then
                  minchange = change
                  mini = i
                  minj = j
                end if
              end do
              !$omp barrier
              jj = jj - tilesize
            end do

            if (i < cities - 2) minchange = minchange + buf(i, bid)
          end do
          !$omp barrier

          change = minchange
          buf_s(lid) = minchange
          if (lid == 0) then
            !$omp atomic update
            climbs(1) = climbs(1) + 1_int32
          end if
          !$omp barrier

          j = dim
          do
            k = (j + 1) / 2
            if (lid + k < j) then
              tmp = buf_s(lid + k)
              if (change > tmp) change = tmp
              buf_s(lid) = change
            end if
            j = k
            !$omp barrier
            if (j <= 1) exit
          end do

          if (minchange == buf_s(0)) buf_s(1) = lid
          !$omp barrier

          if (lid == buf_s(1)) then
            buf_s(2) = mini + 1
            buf_s(3) = minj
          end if
          !$omp barrier

          minchange = buf_s(0)
          mini = buf_s(2)
          sumidx = buf_s(3) + mini
          i = lid
          do while (i + i < sumidx)
            if (mini <= i) then
              j = sumidx - i
              swap_value = pxglob(i, bid)
              pxglob(i, bid) = pxglob(j, bid)
              pxglob(j, bid) = swap_value
              swap_value = pyglob(i, bid)
              pyglob(i, bid) = pyglob(j, bid)
              pyglob(j, bid) = swap_value
            end if
            i = i + dim
          end do
          !$omp barrier

          if (minchange >= 0) exit
        end do

        term = 0
        do i = lid, cities - 1, dim
          term = term + int(sqrt((pxglob(i, bid) - pxglob(i + 1, bid)) * &
            (pxglob(i, bid) - pxglob(i + 1, bid)) + &
            (pyglob(i, bid) - pyglob(i + 1, bid)) * (pyglob(i, bid) - pyglob(i + 1, bid))), int32)
        end do
        buf_s(lid) = term
        !$omp barrier

        j = dim
        do
          k = (j + 1) / 2
          if (lid + k < j) term = term + buf_s(lid + k)
          !$omp barrier
          if (lid + k < j) buf_s(lid) = term
          j = k
          !$omp barrier
          if (j <= 1) exit
        end do

        if (lid == 0) then
          !$omp atomic update
          best(1) = min(best(1), term)
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine two_opt_kernel

end program main
