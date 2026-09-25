! SPDX-License-Identifier: CC0-1.0
program nw_omp_fortran
  use, intrinsic :: iso_c_binding, only: c_int
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use omp_lib, only: omp_get_wtime, omp_get_team_num, omp_get_thread_num
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer, parameter :: block_size = 16
  integer, parameter :: warmup = 100
  integer, parameter :: blosum62(0:23,0:23) = reshape([ &
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4, &
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4, &
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4, &
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4, &
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, &
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4, &
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, &
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4, &
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4, &
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4, &
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4, &
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4, &
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4, &
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4, &
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4, &
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4, &
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4, &
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4, &
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4, &
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4, &
    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4, &
    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, &
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4, &
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1], &
    [24, 24], order=[2, 1])

  integer :: argc, dim, penalty, repeat, ios
  integer :: max_rows, max_cols, n_items, i, j
  integer(int32), allocatable :: reference(:), host_items(:), device_items(:)
  character(len=128) :: arg, program_name
  real(real64) :: start_time, end_time, total_time
  logical :: ok

  write(*, '(A,I0,A)') 'WG size of kernel = ', block_size, ' '

  argc = command_argument_count()
  if (argc /= 3) call usage()

  call get_command_argument(1, arg)
  read(arg, *, iostat=ios) dim
  if (ios /= 0) call usage()
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) penalty
  if (ios /= 0) call usage()
  call get_command_argument(3, arg)
  read(arg, *, iostat=ios) repeat
  if (ios /= 0) call usage()

  if (mod(dim, block_size) /= 0) then
    write(0, '(A)') 'The dimension values must be a multiple of 16'
    stop 1
  end if

  max_rows = dim + 1
  max_cols = dim + 1
  n_items = max_rows * max_cols
  allocate(reference(0:n_items - 1), host_items(0:n_items - 1), device_items(0:n_items - 1))
  reference = 0_int32
  host_items = 0_int32
  device_items = 0_int32

  call init_inputs(reference, host_items, device_items, max_rows, max_cols, penalty)

  !$omp target data map(tofrom: device_items(0:n_items - 1)) map(to: reference(0:n_items - 1))
    do i = 1, warmup + repeat
      if (i == warmup + 1) start_time = omp_get_wtime()
      call nw_device(device_items, reference, max_cols, penalty)
    end do
    end_time = omp_get_wtime()
  !$omp end target data

  total_time = (end_time - start_time) / real(repeat, real64)
  write(*, '(A,F8.6,A)') 'Total kernel execution time: ', total_time, ' (s)'

  call nw_host(host_items, reference, max_cols, penalty)
  ok = all(device_items == host_items)
  if (ok) then
    write(*, '(A)') 'PASS'
  else
    write(*, '(A)') 'FAIL'
    stop 1
  end if

  deallocate(reference, host_items, device_items)

contains

  subroutine usage()
    call get_command_argument(0, program_name)
    write(0, '(A,A,A)') 'Usage: ', trim(program_name), ' <max_rows/max_cols> <penalty> <repeat>'
    write(0, '(A)') char(9)//'<dimension>  - x and y dimensions'
    write(0, '(A)') char(9)//'<penalty> - penalty(positive integer)'
    write(0, '(A)') char(9)//'<repeat> - the number of kernel executions'
    write(0, '(A)') char(9)//'<file> - filename'
    stop 1
  end subroutine usage

  subroutine init_inputs(reference, host_items, device_items, max_rows, max_cols, penalty)
    integer(int32), intent(out) :: reference(0:), host_items(0:), device_items(0:)
    integer, intent(in) :: max_rows, max_cols, penalty
    integer :: i, j

    call c_srand(7_c_int)
    reference = 0_int32
    host_items = 0_int32
    device_items = 0_int32

    do i = 1, max_rows - 1
      host_items(i * max_cols) = int(mod(c_rand(), 10_c_int) + 1_c_int, int32)
      device_items(i * max_cols) = host_items(i * max_cols)
    end do

    do j = 1, max_cols - 1
      host_items(j) = int(mod(c_rand(), 10_c_int) + 1_c_int, int32)
      device_items(j) = host_items(j)
    end do

    do i = 1, max_cols - 1
      do j = 1, max_rows - 1
        reference(i * max_cols + j) = blosum62(device_items(i * max_cols), device_items(j))
      end do
    end do

    do i = 1, max_rows - 1
      host_items(i * max_cols) = -i * penalty
      device_items(i * max_cols) = host_items(i * max_cols)
    end do
    do j = 1, max_cols - 1
      host_items(j) = -j * penalty
      device_items(j) = host_items(j)
    end do
  end subroutine init_inputs

  subroutine nw_device(items, reference, max_cols, penalty)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty
    integer :: blk, block_width

    block_width = (max_cols - 1) / block_size

    do blk = 1, block_width
      !$omp target teams num_teams(blk) thread_limit(block_size)
      block
        integer(int32) :: input_itemsets_l(0:(block_size + 1) * (block_size + 1) - 1)
        integer(int32) :: reference_l(0:block_size * block_size - 1)

        !$omp parallel
        block
          integer :: bx, tx, ty, m, base, b_index_x, b_index_y
          integer :: index, index_n, index_w, index_nw, t_index_x, t_index_y

          bx = omp_get_team_num()
          tx = omp_get_thread_num()

          base = 0
          b_index_x = bx
          b_index_y = blk - 1 - bx

          index = base + max_cols * block_size * b_index_y + block_size * b_index_x + tx + (max_cols + 1)
          index_n = base + max_cols * block_size * b_index_y + block_size * b_index_x + tx + 1
          index_w = base + max_cols * block_size * b_index_y + block_size * b_index_x + max_cols
          index_nw = base + max_cols * block_size * b_index_y + block_size * b_index_x

          if (tx == 0) input_itemsets_l(tx * (block_size + 1)) = items(index_nw + tx)

          do ty = 0, block_size - 1
            reference_l(ty * block_size + tx) = reference(index + max_cols * ty)
          end do

          input_itemsets_l((tx + 1) * (block_size + 1)) = items(index_w + max_cols * tx)
          input_itemsets_l(tx + 1) = items(index_n)

          !$omp barrier

          do m = 0, block_size - 1
            if (tx <= m) then
              t_index_x = tx + 1
              t_index_y = m - tx + 1

              input_itemsets_l(t_index_y * (block_size + 1) + t_index_x) = max( &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x - 1) + &
                  reference_l((t_index_y - 1) * block_size + t_index_x - 1), &
                input_itemsets_l(t_index_y * (block_size + 1) + t_index_x - 1) - penalty, &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x) - penalty)
            end if
            !$omp barrier
          end do

          do m = block_size - 2, 0, -1
            if (tx <= m) then
              t_index_x = tx + block_size - m
              t_index_y = block_size - tx

              input_itemsets_l(t_index_y * (block_size + 1) + t_index_x) = max( &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x - 1) + &
                  reference_l((t_index_y - 1) * block_size + t_index_x - 1), &
                input_itemsets_l(t_index_y * (block_size + 1) + t_index_x - 1) - penalty, &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x) - penalty)
            end if
            !$omp barrier
          end do

          do ty = 0, block_size - 1
            items(index + max_cols * ty) = input_itemsets_l((ty + 1) * (block_size + 1) + tx + 1)
          end do
        end block
        !$omp end parallel
      end block
      !$omp end target teams
    end do

    do blk = block_width - 1, 1, -1
      !$omp target teams num_teams(blk) thread_limit(block_size)
      block
        integer(int32) :: input_itemsets_l(0:(block_size + 1) * (block_size + 1) - 1)
        integer(int32) :: reference_l(0:block_size * block_size - 1)

        !$omp parallel
        block
          integer :: bx, tx, ty, m, base, b_index_x, b_index_y
          integer :: index, index_n, index_w, index_nw, t_index_x, t_index_y

          bx = omp_get_team_num()
          tx = omp_get_thread_num()

          base = 0
          b_index_x = bx + block_width - blk
          b_index_y = block_width - bx - 1

          index = base + max_cols * block_size * b_index_y + block_size * b_index_x + tx + (max_cols + 1)
          index_n = base + max_cols * block_size * b_index_y + block_size * b_index_x + tx + 1
          index_w = base + max_cols * block_size * b_index_y + block_size * b_index_x + max_cols
          index_nw = base + max_cols * block_size * b_index_y + block_size * b_index_x

          if (tx == 0) input_itemsets_l(tx * (block_size + 1)) = items(index_nw)

          do ty = 0, block_size - 1
            reference_l(ty * block_size + tx) = reference(index + max_cols * ty)
          end do

          input_itemsets_l((tx + 1) * (block_size + 1)) = items(index_w + max_cols * tx)
          input_itemsets_l(tx + 1) = items(index_n)

          !$omp barrier

          do m = 0, block_size - 1
            if (tx <= m) then
              t_index_x = tx + 1
              t_index_y = m - tx + 1

              input_itemsets_l(t_index_y * (block_size + 1) + t_index_x) = max( &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x - 1) + &
                  reference_l((t_index_y - 1) * block_size + t_index_x - 1), &
                input_itemsets_l(t_index_y * (block_size + 1) + t_index_x - 1) - penalty, &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x) - penalty)
            end if
            !$omp barrier
          end do

          do m = block_size - 2, 0, -1
            if (tx <= m) then
              t_index_x = tx + block_size - m
              t_index_y = block_size - tx

              input_itemsets_l(t_index_y * (block_size + 1) + t_index_x) = max( &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x - 1) + &
                  reference_l((t_index_y - 1) * block_size + t_index_x - 1), &
                input_itemsets_l(t_index_y * (block_size + 1) + t_index_x - 1) - penalty, &
                input_itemsets_l((t_index_y - 1) * (block_size + 1) + t_index_x) - penalty)
            end if
            !$omp barrier
          end do

          do ty = 0, block_size - 1
            items(index + ty * max_cols) = input_itemsets_l((ty + 1) * (block_size + 1) + tx + 1)
          end do
        end block
        !$omp end parallel
      end block
      !$omp end target teams
    end do
  end subroutine nw_device

  subroutine nw_host(items, reference, max_cols, penalty)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty
    integer :: blk, bx, block_width

    block_width = (max_cols - 1) / block_size
    do blk = 1, block_width
      do bx = 0, blk - 1
        call compute_block(items, reference, max_cols, penalty, bx, blk - 1 - bx)
      end do
    end do
    do blk = 2, block_width
      do bx = blk - 1, block_width - 1
        call compute_block(items, reference, max_cols, penalty, bx, block_width + blk - 2 - bx)
      end do
    end do
  end subroutine nw_host

  subroutine compute_block(items, reference, max_cols, penalty, b_index_x, b_index_y)
    integer(int32), intent(inout) :: items(0:)
    integer(int32), intent(in) :: reference(0:)
    integer, intent(in) :: max_cols, penalty, b_index_x, b_index_y
    integer :: ii, jj, global_i, global_j, base
    integer(int32) :: local_items(0:block_size,0:block_size)
    integer(int32) :: local_ref(1:block_size,1:block_size)

    base = max_cols * block_size * b_index_y + block_size * b_index_x
    do ii = 0, block_size
      do jj = 0, block_size
        local_items(ii, jj) = items(base + ii * max_cols + jj)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        global_i = b_index_y * block_size + ii
        global_j = b_index_x * block_size + jj
        local_ref(ii, jj) = reference(global_i * max_cols + global_j)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        local_items(ii, jj) = max(local_items(ii - 1, jj - 1) + local_ref(ii, jj), &
                                  local_items(ii, jj - 1) - penalty, &
                                  local_items(ii - 1, jj) - penalty)
      end do
    end do

    do ii = 1, block_size
      do jj = 1, block_size
        items(base + ii * max_cols + jj) = local_items(ii, jj)
      end do
    end do
  end subroutine compute_block

end program nw_omp_fortran
