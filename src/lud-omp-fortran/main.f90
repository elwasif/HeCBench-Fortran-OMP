! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 16
  integer :: matrix_dim
  logical :: do_verify
  character(len=512) :: input_file
  logical :: have_input
  real(real32), allocatable :: m(:), mm(:)
  real(real64) :: start_total, end_total, start_kernel, end_kernel

  matrix_dim = 32
  do_verify = .false.
  input_file = ''
  have_input = .false.
  call parse_args(matrix_dim, do_verify, input_file, have_input)

  if (have_input) then
    write(*, '(A,A)') 'Reading matrix from file ', trim(input_file)
    call create_matrix_from_file(m, input_file, matrix_dim)
  else
    write(*, '(A,I0)') 'Creating matrix internally size=', matrix_dim
    call create_matrix(m, matrix_dim)
  end if

  if (do_verify) then
    print '(A)', 'Before LUD'
    allocate(mm(size(m)))
    mm = m
  end if

  print '(A,I0,A,I0)', 'WG size of kernel = ', block_size, ' X ', block_size
  start_total = omp_get_wtime()
  call lud_offload(m, matrix_dim)
  end_total = omp_get_wtime()
  write(*, '(A,F8.6)') 'Device offloading time (s): ', end_total - start_total

  if (do_verify) then
    print '(A)', 'After LUD'
    print '(A)', '>>>Verify<<<<'
    call lud_verify(mm, m, matrix_dim)
    deallocate(mm)
  end if

  deallocate(m)

contains

  subroutine parse_args(matrix_dim, do_verify, input_file, have_input)
    integer, intent(inout) :: matrix_dim
    logical, intent(out) :: do_verify, have_input
    character(len=*), intent(out) :: input_file
    integer :: argc, idx
    character(len=512) :: arg, value

    do_verify = .false.
    have_input = .false.
    argc = command_argument_count()
    if (argc == 0) then
      call usage_stop()
    end if

    idx = 1
    do while (idx <= argc)
      call get_command_argument(idx, arg)
      select case (trim(arg))
      case ('-v', '--verify')
        do_verify = .true.
        idx = idx + 1
      case ('-s', '--size')
        if (idx == argc) then
          print '(A)', 'missing argument'
          stop 1
        end if
        call get_command_argument(idx + 1, value)
        read(value, *) matrix_dim
        if (matrix_dim <= 0) then
          print '(A)', 'Matrix dimension must be positive!'
          stop 1
        end if
        if (mod(matrix_dim, block_size) /= 0) then
          write(*, '(A,I0,A)') 'Matrix dimension of ', matrix_dim, ' not supported by the benchmark'
          stop 1
        end if
        write(*, '(A,I0)') 'Generate input matrix internally, size =', matrix_dim
        idx = idx + 2
      case ('-i', '--input')
        if (idx == argc) then
          print '(A)', 'missing argument'
          stop 1
        end if
        call get_command_argument(idx + 1, input_file)
        have_input = .true.
        idx = idx + 2
      case default
        call usage_stop()
      end select
    end do
  end subroutine parse_args

  subroutine usage_stop()
    write(*, '(A)') 'Usage: ./main [-v] [-s matrix_size|-i input_file]'
    stop 1
  end subroutine usage_stop

  subroutine create_matrix(m, n)
    real(real32), allocatable, intent(out) :: m(:)
    integer, intent(in) :: n
    real(real32), allocatable :: coe(:)
    real(real32) :: lamda, coe_i
    integer :: i, j

    allocate(m(n * n), coe(2 * n - 1))
    lamda = -0.001_real32
    do i = 0, n - 1
      coe_i = 10.0_real32 * exp(lamda * real(i, real32))
      coe(n + i) = coe_i
      coe(n - i) = coe_i
    end do

    do i = 0, n - 1
      do j = 0, n - 1
        m(i * n + j + 1) = coe(n - i + j)
      end do
    end do
    deallocate(coe)
  end subroutine create_matrix

  subroutine create_matrix_from_file(m, filename, n)
    real(real32), allocatable, intent(out) :: m(:)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n
    integer :: unit, ios, i, j

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*, '(A,A)') 'error create matrix from file ', trim(filename)
      stop 1
    end if
    read(unit, *, iostat=ios) n
    if (ios /= 0 .or. n <= 0) then
      write(*, '(A,A)') 'error create matrix from file ', trim(filename)
      close(unit)
      stop 1
    end if
    allocate(m(n * n))
    do i = 0, n - 1
      read(unit, *, iostat=ios) (m(i * n + j + 1), j = 0, n - 1)
      if (ios /= 0) then
        write(*, '(A,A)') 'error create matrix from file ', trim(filename)
        close(unit)
        stop 1
      end if
    end do
    close(unit)
  end subroutine create_matrix_from_file

  subroutine lud_offload(m, n)
    real(real32), intent(inout), contiguous :: m(:)
    integer, intent(in) :: n
    integer :: offset, tile_count

    !$omp target data map(tofrom: m(1:n*n))
      start_kernel = omp_get_wtime()
      offset = 0
      do while (offset < n - block_size)
        call diagonal_kernel(m, n, offset)
        tile_count = (n - offset) / block_size - 1
        call perimeter_kernel(m, n, offset, tile_count)
        call internal_kernel(m, n, offset, tile_count)
        offset = offset + block_size
      end do
      call diagonal_kernel(m, n, offset)
      end_kernel = omp_get_wtime()
      write(*, '(A,F8.6,A)') 'Total kernel execution time : ', end_kernel - start_kernel, ' (s)'
    !$omp end target data
  end subroutine lud_offload

  subroutine diagonal_kernel(m, n, offset)
    real(real32), intent(inout), contiguous :: m(:)
    integer, intent(in) :: n, offset
    integer :: tx, i, j, array_offset
    real(real32) :: shadow(block_size * block_size)

    !$omp target teams num_teams(1) thread_limit(block_size) private(shadow)
      shadow = 0.0_real32
      !$omp parallel private(tx, i, j, array_offset) shared(shadow, m)
        tx = omp_get_thread_num()
        array_offset = offset * n + offset
        do i = 0, block_size - 1
          shadow(i * block_size + tx + 1) = m(array_offset + tx + 1)
          array_offset = array_offset + n
        end do
        !$omp barrier
        do i = 0, block_size - 2
          if (tx > i) then
            do j = 0, i - 1
              shadow(tx * block_size + i + 1) = shadow(tx * block_size + i + 1) - &
                shadow(tx * block_size + j + 1) * shadow(j * block_size + i + 1)
            end do
            shadow(tx * block_size + i + 1) = shadow(tx * block_size + i + 1) / shadow(i * block_size + i + 1)
          end if
          !$omp barrier
          if (tx > i) then
            do j = 0, i
              shadow((i + 1) * block_size + tx + 1) = shadow((i + 1) * block_size + tx + 1) - &
                shadow((i + 1) * block_size + j + 1) * shadow(j * block_size + tx + 1)
            end do
          end if
          !$omp barrier
        end do
        array_offset = (offset + 1) * n + offset
        do i = 1, block_size - 1
          m(array_offset + tx + 1) = shadow(i * block_size + tx + 1)
          array_offset = array_offset + n
        end do
      !$omp end parallel
    !$omp end target teams
  end subroutine diagonal_kernel

  subroutine perimeter_kernel(m, n, offset, tile_count)
    real(real32), intent(inout), contiguous :: m(:)
    integer, intent(in) :: n, offset, tile_count
    integer :: bx, tx, idx, i, j, array_offset
    real(real32) :: dia(block_size * block_size), peri_row(block_size * block_size), peri_col(block_size * block_size)

    !$omp target teams num_teams(tile_count) thread_limit(2 * block_size) private(dia, peri_row, peri_col)
      dia = 0.0_real32
      peri_row = 0.0_real32
      peri_col = 0.0_real32
      !$omp parallel private(bx, tx, idx, i, j, array_offset) shared(dia, peri_row, peri_col, m)
        bx = omp_get_team_num()
        tx = omp_get_thread_num()
        if (tx < block_size) then
          idx = tx
          array_offset = offset * n + offset
          do i = 0, block_size / 2 - 1
            dia(i * block_size + idx + 1) = m(array_offset + idx + 1)
            array_offset = array_offset + n
          end do
          array_offset = offset * n + offset
          do i = 0, block_size - 1
            peri_row(i * block_size + idx + 1) = m(array_offset + (bx + 1) * block_size + idx + 1)
            array_offset = array_offset + n
          end do
        else
          idx = tx - block_size
          array_offset = (offset + block_size / 2) * n + offset
          do i = block_size / 2, block_size - 1
            dia(i * block_size + idx + 1) = m(array_offset + idx + 1)
            array_offset = array_offset + n
          end do
          array_offset = (offset + (bx + 1) * block_size) * n + offset
          do i = 0, block_size - 1
            peri_col(i * block_size + idx + 1) = m(array_offset + idx + 1)
            array_offset = array_offset + n
          end do
        end if
        !$omp barrier
        if (tx < block_size) then
          idx = tx
          do i = 1, block_size - 1
            do j = 0, i - 1
              peri_row(i * block_size + idx + 1) = peri_row(i * block_size + idx + 1) - &
                dia(i * block_size + j + 1) * peri_row(j * block_size + idx + 1)
            end do
          end do
        else
          idx = tx - block_size
          do i = 0, block_size - 1
            do j = 0, i - 1
              peri_col(idx * block_size + i + 1) = peri_col(idx * block_size + i + 1) - &
                peri_col(idx * block_size + j + 1) * dia(j * block_size + i + 1)
            end do
            peri_col(idx * block_size + i + 1) = peri_col(idx * block_size + i + 1) / dia(i * block_size + i + 1)
          end do
        end if
        !$omp barrier
        if (tx < block_size) then
          idx = tx
          array_offset = (offset + 1) * n + offset
          do i = 1, block_size - 1
            m(array_offset + (bx + 1) * block_size + idx + 1) = peri_row(i * block_size + idx + 1)
            array_offset = array_offset + n
          end do
        else
          idx = tx - block_size
          array_offset = (offset + (bx + 1) * block_size) * n + offset
          do i = 0, block_size - 1
            m(array_offset + idx + 1) = peri_col(i * block_size + idx + 1)
            array_offset = array_offset + n
          end do
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine perimeter_kernel

  subroutine internal_kernel(m, n, offset, tile_count)
    real(real32), intent(inout), contiguous :: m(:)
    integer, intent(in) :: n, offset, tile_count
    integer :: team, bx, by, tx, ty, i, global_row_id, global_col_id
    real(real32) :: sum
    real(real32) :: peri_row(block_size * block_size), peri_col(block_size * block_size)

    !$omp target teams num_teams(tile_count * tile_count) thread_limit(block_size * block_size) private(peri_row, peri_col)
      peri_row = 0.0_real32
      peri_col = 0.0_real32
      !$omp parallel private(team, bx, by, tx, ty, i, global_row_id, global_col_id, sum) shared(peri_row, peri_col, m)
        team = omp_get_team_num()
        bx = mod(team, tile_count)
        by = team / tile_count
        tx = mod(omp_get_thread_num(), block_size)
        ty = omp_get_thread_num() / block_size
        global_row_id = offset + (by + 1) * block_size
        global_col_id = offset + (bx + 1) * block_size
        peri_row(ty * block_size + tx + 1) = m((offset + ty) * n + global_col_id + tx + 1)
        peri_col(ty * block_size + tx + 1) = m((global_row_id + ty) * n + offset + tx + 1)
        !$omp barrier
        sum = 0.0_real32
        do i = 0, block_size - 1
          sum = sum + peri_col(ty * block_size + i + 1) * peri_row(i * block_size + tx + 1)
        end do
        m((global_row_id + ty) * n + global_col_id + tx + 1) = &
          m((global_row_id + ty) * n + global_col_id + tx + 1) - sum
      !$omp end parallel
    !$omp end target teams
  end subroutine internal_kernel

  subroutine lud_verify(original, lu, n)
    real(real32), intent(in) :: original(:), lu(:)
    integer, intent(in) :: n
    real(real32), allocatable :: tmp(:)
    real(real32) :: sum, l_value, u_value
    integer :: i, j, k

    allocate(tmp(n * n))
    do i = 0, n - 1
      do j = 0, n - 1
        sum = 0.0_real32
        do k = 0, min(i, j)
          if (i == k) then
            l_value = 1.0_real32
          else
            l_value = lu(i * n + k + 1)
          end if
          u_value = lu(k * n + j + 1)
          sum = sum + l_value * u_value
        end do
        tmp(i * n + j + 1) = sum
      end do
    end do

    do i = 0, n - 1
      do j = 0, n - 1
        if (abs(original(i * n + j + 1) - tmp(i * n + j + 1)) > 1.0e-4_real32) then
          write(*, '(A,I0,A,I0,A,F0.6,A,F0.6)') 'mismatch at (', i, ', ', j, '): (o)', &
            original(i * n + j + 1), ' (n)', tmp(i * n + j + 1)
        end if
      end do
    end do
    deallocate(tmp)
  end subroutine lud_verify

end program main
