! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: nelems, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <number of elements> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) nelems
  read(arg2, *) repeat
  if (nelems <= 0 .or. repeat <= 0) stop 1

  call test_block(64, nelems, repeat)
  call test_block(128, nelems, repeat)
  call test_block(256, nelems, repeat)
  call test_block(512, nelems, repeat)
  call test_block(1024, nelems, repeat)

contains

  subroutine test_block(block_threads, requested_items, repeat)
    integer, intent(in) :: block_threads, requested_items, repeat
    integer, parameter :: items_per_thread = 4
    integer :: items_per_block, num_items, grid_size, i, rep
    integer, allocatable :: h_in(:), h_out(:), r_out(:)
    real(real64) :: start_time, elapsed
    logical :: ok

    items_per_block = block_threads * items_per_thread
    num_items = ((requested_items + items_per_block - 1) / items_per_block) * items_per_block
    grid_size = num_items / items_per_block

    allocate(h_in(num_items), h_out(num_items), r_out(num_items))
    do i = 1, num_items
      h_in(i) = mod(i - 1, 17)
    end do

    !$omp target data map(to: h_in(1:num_items)) map(alloc: h_out(1:num_items))
    do rep = 1, repeat
      call adjacent_kernel(h_in, h_out, .true., num_items, items_per_block, grid_size)
    end do
    !$omp target update from(h_out(1:num_items))
    call reference_adjacent(h_in, r_out, .true., num_items, items_per_block, grid_size)
    ok = all(h_out == r_out)
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    do rep = 1, repeat
      call adjacent_kernel(h_in, h_out, .false., num_items, items_per_block, grid_size)
    end do
    !$omp target update from(h_out(1:num_items))
    call reference_adjacent(h_in, r_out, .false., num_items, items_per_block, grid_size)
    ok = all(h_out == r_out)
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    start_time = omp_get_wtime()
    do rep = 1, repeat
      call adjacent_kernel(h_in, h_out, .true., num_items, items_per_block, grid_size)
      call adjacent_kernel(h_out, h_out, .false., num_items, items_per_block, grid_size)
    end do
    elapsed = omp_get_wtime() - start_time
    write(*,'(A,I4,A,F0.6,A)') 'Average execution time of the kernels (thread block size = ', &
      block_threads, '): ', elapsed * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp end target data

    deallocate(h_in, h_out, r_out)
  end subroutine test_block

  subroutine adjacent_kernel(input, output, subtract_left, num_items, items_per_block, grid_size)
    integer, intent(in) :: input(:), num_items, items_per_block, grid_size
    integer, intent(out) :: output(:)
    logical, intent(in) :: subtract_left
    integer :: b, i, idx

    !$omp target teams distribute private(i, idx)
    do b = 0, grid_size - 1
      !$omp parallel do private(idx)
      do i = 0, items_per_block - 1
        idx = b * items_per_block + i + 1
        if (subtract_left) then
          if (i == 0) then
            output(idx) = input(idx)
          else
            output(idx) = input(idx) - input(idx - 1)
          end if
        else
          if (i == items_per_block - 1) then
            output(idx) = input(idx)
          else
            output(idx) = input(idx) - input(idx + 1)
          end if
        end if
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine adjacent_kernel

  subroutine reference_adjacent(input, output, subtract_left, num_items, items_per_block, grid_size)
    integer, intent(in) :: input(:), num_items, items_per_block, grid_size
    integer, intent(out) :: output(:)
    logical, intent(in) :: subtract_left
    integer :: b, i, idx

    do b = 0, grid_size - 1
      do i = 0, items_per_block - 1
        idx = b * items_per_block + i + 1
        if (subtract_left) then
          if (i == 0) then
            output(idx) = input(idx)
          else
            output(idx) = input(idx) - input(idx - 1)
          end if
        else
          if (i == items_per_block - 1) then
            output(idx) = input(idx)
          else
            output(idx) = input(idx) - input(idx + 1)
          end if
        end if
      end do
    end do
  end subroutine reference_adjacent

end program main
