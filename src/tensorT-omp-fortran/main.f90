! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: TILE_SIZE_CONST = 5900
  integer, parameter :: NTHREADS = 256
  integer, parameter :: d1 = 41, d2 = 13, d3 = 11, d4 = 9, d5 = 76, d6 = 50
  integer(int64), parameter :: data_size = int(d1, int64) * int(d2, int64) * &
    int(d3, int64) * int(d4, int64) * int(d5, int64) * int(d6, int64)
  character(len=256) :: arg0, arg1
  integer :: repeat, iter
  integer(int64) :: idx
  real(real64), allocatable :: input(:), output(:)
  real(real64) :: start_time, end_time
  integer :: nblocks, tile_size, dim_output, dim_input
  integer :: shape_output(0:2), shape_input(0:2)
  real(real64) :: shape_output_r(0:2), shape_input_r(0:2)
  integer :: stride_output_local(0:2), stride_output_global(0:2), stride_input(0:2)

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat
  if (repeat <= 0) stop 1

  allocate(input(0:data_size - 1), output(0:data_size - 1))
  !$omp parallel do
  do idx = 0, data_size - 1
    input(idx) = real(idx, real64)
    output(idx) = 0.0_real64
  end do
  !$omp end parallel do

  shape_output = [d2, d3, d1]
  shape_input = [d4, d5, d6]
  shape_output_r = [1.0_real64 / real(d2, real64), 1.0_real64 / real(d3, real64), &
    1.0_real64 / real(d1, real64)]
  shape_input_r = [1.0_real64 / real(d4, real64), 1.0_real64 / real(d5, real64), &
    1.0_real64 / real(d6, real64)]
  stride_output_local = [d1, d1 * d2, 1]
  stride_output_global = [1, d2, d2 * d3 * d4 * d6]
  stride_input = [d2 * d3, d2 * d3 * d4 * d6 * d1, d2 * d3 * d4]
  nblocks = d4 * d5 * d6
  tile_size = d1 * d2 * d3
  dim_output = 3
  dim_input = 3

  !$omp target data map(to: input(0:data_size - 1), shape_input(0:dim_input - 1), &
  !$omp& shape_input_r(0:dim_input - 1), shape_output(0:dim_output - 1), &
  !$omp& shape_output_r(0:dim_output - 1), stride_input(0:dim_input - 1), &
  !$omp& stride_output_local(0:dim_output - 1), stride_output_global(0:dim_output - 1)) &
  !$omp& map(from: output(0:data_size - 1))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call tensor_transpose(input, output, nblocks, tile_size, dim_input, dim_output, &
      shape_input, shape_input_r, shape_output, shape_output_r, stride_input, &
      stride_output_local, stride_output_global)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', &
    (end_time - start_time) * 1.0e3_real64 / real(repeat, real64), ' (ms)'
  !$omp end target data

  call verify(input, output)
  deallocate(input, output)

contains

  subroutine tensor_transpose(input, output, nblocks, tile_size, dim_input, dim_output, &
      shape_input, shape_input_r, shape_output, shape_output_r, stride_input, &
      stride_output_local, stride_output_global)
    real(real64), intent(in) :: input(0:)
    real(real64), intent(out) :: output(0:)
    integer, intent(in) :: nblocks, tile_size, dim_input, dim_output
    integer, intent(in) :: shape_input(0:), shape_output(0:)
    real(real64), intent(in) :: shape_input_r(0:), shape_output_r(0:)
    integer, intent(in) :: stride_input(0:), stride_output_local(0:), stride_output_global(0:)
    real(real64) :: tile(0:TILE_SIZE_CONST - 1)
    integer :: block_idx, i, j, it, im, offset1, offset2, local_offset, tmp

    !$omp target teams num_teams(nblocks) thread_limit(NTHREADS) private(tile, block_idx, i, j, it, im, offset1, offset2, local_offset, tmp)
    !$omp parallel private(i, j, it, im, offset1, offset2, local_offset, tmp)
    do block_idx = omp_get_team_num(), nblocks - 1, omp_get_num_teams()
      it = block_idx
      im = 0
      offset1 = 0
      do i = 0, dim_input - 1
        im = int(real(it, real64) * shape_input_r(i))
        offset1 = offset1 + stride_input(i) * (it - im * shape_input(i))
        it = im
      end do

      do i = omp_get_thread_num(), tile_size - 1, omp_get_num_threads()
        tile(i) = input(i + block_idx * tile_size)
      end do

      !$omp barrier

      do i = omp_get_thread_num(), tile_size - 1, omp_get_num_threads()
        it = i
        offset2 = 0
        local_offset = 0
        do j = 0, dim_output - 1
          im = int(real(it, real64) * shape_output_r(j))
          tmp = it - im * shape_output(j)
          offset2 = offset2 + stride_output_global(j) * tmp
          local_offset = local_offset + stride_output_local(j) * tmp
          it = im
        end do
        output(offset1 + offset2) = tile(local_offset)
      end do
      !$omp barrier
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine tensor_transpose

  subroutine verify(input, output)
    real(real64), intent(in) :: input(0:), output(0:)
    integer(int64) :: input_offset, output_offset, step_input, step_output
    integer :: i
    logical :: error

    input_offset = 2_int64 + int(d1, int64) * (2_int64 + int(d2, int64) * &
      (2_int64 + int(d3, int64) * (2_int64 + int(d4, int64) * &
      (0_int64 + 2_int64 * int(d5, int64)))))
    output_offset = 2_int64 + int(d2, int64) * (2_int64 + int(d3, int64) * &
      (2_int64 + int(d4, int64) * (2_int64 + int(d6, int64) * &
      (2_int64 + 0_int64 * int(d1, int64)))))
    step_input = int(d1, int64) * int(d2, int64) * int(d3, int64) * int(d4, int64)
    step_output = int(d2, int64) * int(d3, int64) * int(d4, int64) * int(d6, int64) * int(d1, int64)

    error = .false.
    do i = 0, d5 - 1
      if (input(input_offset + int(i, int64) * step_input) /= &
          output(output_offset + int(i, int64) * step_output)) then
        write(*,'(A)') 'FAIL'
        error = .true.
        exit
      end if
    end do
    if (.not. error) write(*,'(A)') 'PASS'
  end subroutine verify

end program main
