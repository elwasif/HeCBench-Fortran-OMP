! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use omp_lib
  implicit none

  integer, parameter :: m_seed = 9
  integer, parameter :: halo = 1
  integer, parameter :: lws = 250
  integer :: rows, cols, pyramid_height, total, t, iteration
  integer :: active_parity, mismatches, borderCols, gws
  integer(int32), allocatable :: wall(:), gpuSrc(:), gpuResult(:), reference(:), outputBuffer(:)
  logical :: src_is_gpuSrc
  real(real64) :: offload_start, offload_end, kernel_start, kernel_end

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name="rand")
      import :: c_int
    end function c_rand
  end interface

  if (command_argument_count() == 3) then
    cols = read_arg(1)
    rows = read_arg(2)
    pyramid_height = read_arg(3)
  else
    print '(A)', 'Usage: ./main <column length> <row length> <pyramid_height>'
    stop 0
  end if

  if (cols <= 0 .or. rows <= 0 .or. pyramid_height <= 0) then
    error stop 'pathfinder dimensions and pyramid_height must be positive'
  end if

  total = rows * cols
  allocate(wall(total), gpuSrc(cols), gpuResult(cols), reference(cols), outputBuffer(16384))
  call initialize_wall(wall, rows, cols)
  gpuSrc = wall(1:cols)
  gpuResult = 0_int32
  outputBuffer = 0_int32
  reference = gpuSrc
  borderCols = pyramid_height * halo
  gws = total / lws
  src_is_gpuSrc = .true.

  offload_start = omp_get_wtime()
  kernel_start = omp_get_wtime()

  !$omp target data map(to: gpuSrc(1:cols), wall(cols + 1:total)) &
  !$omp& map(alloc: gpuResult(1:cols)) map(from: outputBuffer(1:16384))
  do t = 0, rows - 2, pyramid_height
    if (t == pyramid_height) then
      kernel_start = omp_get_wtime()
    end if

    iteration = min(pyramid_height, rows - t - 1)
    if (src_is_gpuSrc) then
      call launch_kernel(gpuSrc, gpuResult, wall, outputBuffer, cols, iteration, t, borderCols, gws)
    else
      call launch_kernel(gpuResult, gpuSrc, wall, outputBuffer, cols, iteration, t, borderCols, gws)
    end if
    src_is_gpuSrc = .not. src_is_gpuSrc
  end do
  if (src_is_gpuSrc) then
    !$omp target update from(gpuSrc(1:cols))
  else
    !$omp target update from(gpuResult(1:cols))
  end if
  !$omp end target data

  kernel_end = omp_get_wtime()
  offload_end = omp_get_wtime()

  call cpu_reference(wall, rows, cols, reference)
  active_parity = merge(1, 0, src_is_gpuSrc)
  if (active_parity == 1) then
    mismatches = count(gpuSrc /= reference)
  else
    mismatches = count(gpuResult /= reference)
  end if
  if (mismatches /= 0) then
    write(*, '(A,I0)') 'FAIL: mismatches=', mismatches
    error stop 'pathfinder validation failed'
  end if

  write(*, '(A,F8.6,A)') 'Total kernel execution time: ', kernel_end - kernel_start, ' (s)'
  write(*, '(A,F8.6,A)') 'Device offloading time = ', offload_end - offload_start, ' (s)'

  deallocate(wall, gpuSrc, gpuResult, reference, outputBuffer)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_wall(wall, rows, cols)
    integer(int32), intent(out) :: wall(:)
    integer, intent(in) :: rows, cols
    integer :: i, j, idx
    call c_srand(int(m_seed, c_int))
    do j = 1, rows
      do i = 1, cols
        idx = (j - 1) * cols + i
        wall(idx) = int(mod(c_rand(), 10_c_int), int32)
      end do
    end do
  end subroutine initialize_wall

  subroutine cpu_reference(wall, rows, cols, result)
    integer(int32), intent(in) :: wall(:)
    integer, intent(in) :: rows, cols
    integer(int32), intent(out) :: result(:)
    integer(int32), allocatable :: scratch(:)
    integer :: row, col, left_col, right_col

    allocate(scratch(cols))
    result = wall(1:cols)
    scratch = 0_int32
    do row = 2, rows
      do col = 1, cols
        left_col = max(1, col - 1)
        right_col = min(cols, col + 1)
        scratch(col) = wall((row - 1) * cols + col) + min(result(left_col), min(result(col), result(right_col)))
      end do
      result = scratch
    end do
    deallocate(scratch)
  end subroutine cpu_reference

  subroutine launch_kernel(gpuSrc, gpuResult, wall, outputBuffer, cols, iteration, t, borderCols, gws)
    integer(int32), intent(in) :: gpuSrc(:), wall(:)
    integer(int32), intent(inout) :: outputBuffer(:)
    integer(int32), intent(out) :: gpuResult(:)
    integer, intent(in) :: cols, iteration, t, borderCols, gws
    integer :: theHalo

    theHalo = halo
    !$omp target teams num_teams(gws) thread_limit(lws)
    block
      integer(int32) :: prev(0:lws - 1), result(0:lws - 1)
      !$omp parallel
      block
        integer :: BLOCK_SIZE, bx, tx, small_block_cols, blkX, blkXmax
        integer :: xidx, validXmin, validXmax, W, E, i, index, bufIndex
        integer(int32) :: left, up, right, shortest
        logical :: isValid, computed

        BLOCK_SIZE = omp_get_num_threads()
        bx = omp_get_team_num()
        tx = omp_get_thread_num()

        small_block_cols = BLOCK_SIZE - (iteration * theHalo * 2)
        blkX = (small_block_cols * bx) - borderCols
        blkXmax = blkX + BLOCK_SIZE - 1
        xidx = blkX + tx

        validXmin = merge(-blkX, 0, blkX < 0)
        validXmax = merge(BLOCK_SIZE - 1 - (blkXmax - cols + 1), BLOCK_SIZE - 1, blkXmax > cols - 1)

        W = max(validXmin, tx - 1)
        E = min(validXmax, tx + 1)
        isValid = tx >= validXmin .and. tx <= validXmax

        if (xidx >= 0 .and. xidx <= cols - 1) then
          prev(tx) = gpuSrc(xidx + 1)
        end if

        !$omp barrier

        do i = 0, iteration - 1
          computed = .false.

          if (tx >= i + 1 .and. tx <= BLOCK_SIZE - i - 2 .and. isValid) then
            computed = .true.
            left = prev(W)
            up = prev(tx)
            right = prev(E)
            shortest = min(left, up)
            shortest = min(shortest, right)

            index = cols * (t + i + 1) + xidx + 1
            result(tx) = shortest + wall(index)

            if (tx == 11 .and. i == 0) then
              bufIndex = gpuSrc(xidx + 1)
              if (bufIndex >= 0 .and. bufIndex <= 16383) then
                outputBuffer(bufIndex + 1) = 1_int32
              end if
            end if
          end if

          !$omp barrier

          if (i == iteration - 1) then
            exit
          end if

          if (computed) then
            prev(tx) = result(tx)
          end if
          !$omp barrier
        end do

        if (computed) then
          gpuResult(xidx + 1) = result(tx)
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine launch_kernel

end program main
