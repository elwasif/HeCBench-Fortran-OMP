! SPDX-License-Identifier: CC0-1.0
program dense_embedding_main
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_embed_dims = 3
  integer, parameter :: embed_dims(num_embed_dims) = [768, 2048, 12288]
  integer :: nrows, batch_size, repeat
  character(len=64) :: arg

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand

    subroutine dense_embedding_random_fill(dense, input, dense_size, input_size) bind(C, name="dense_embedding_random_fill")
      import :: c_float, c_int
      real(c_float) :: dense(*), input(*)
      integer(c_int), value :: dense_size, input_size
    end subroutine dense_embedding_random_fill
  end interface

  if (command_argument_count() /= 3) then
    write(*,'("Usage: ./main <number of rows> <batch size> <repeat>")')
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nrows
  call get_command_argument(2, arg)
  read(arg, *) batch_size
  call get_command_argument(3, arg)
  read(arg, *) repeat

  if (nrows <= batch_size * batch_size) then
    stop 1
  end if

  write(*,'("Number of rows in the embedding table: ",I0)') nrows
  write(*,'("Batch size: ",I0)') batch_size
  call run_sweep(nrows, batch_size, repeat)

contains

  subroutine run_sweep(nrows, batch_size, repeat)
    integer, intent(in) :: nrows, batch_size, repeat
    integer :: dim_index, ncols, input_size, dense_size, block_size, i
    integer, allocatable :: offset(:)
    real(real32), allocatable :: input(:), dense(:), output_k1(:), output_k2(:), output_k3(:), output_ref(:)
    real(real64) :: start_time, elapsed_us
    logical :: ok

    do dim_index = 1, num_embed_dims
      ncols = embed_dims(dim_index)
      write(*,*)
      write(*,'("Embedding dimension: ",I0)') ncols

      input_size = nrows * ncols
      dense_size = batch_size * ncols
      allocate(input(0:input_size-1), dense(0:dense_size-1))
      allocate(output_k1(0:input_size-1), output_k2(0:input_size-1), output_k3(0:input_size-1), output_ref(0:input_size-1))
      allocate(offset(0:batch_size))

      call c_srand(123_c_int)
      offset(0) = 0
      do i = 1, batch_size
        offset(i) = offset(i-1) + (mod(c_rand(), batch_size) + 1) * ncols
      end do

      call dense_embedding_random_fill(dense, input, int(dense_size, c_int), int(input_size, c_int))
      do i = 0, input_size - 1
        output_k1(i) = 0.0_real32
        output_k2(i) = 0.0_real32
        output_k3(i) = 0.0_real32
        output_ref(i) = 0.0_real32
      end do

      call dense_reference(input, dense, output_ref, ncols, batch_size, offset)

      !$omp target data map(to: input, dense, offset, output_k1, output_k2, output_k3)
      do block_size = 128, 1024, 128
        if (block_size /= 128 .and. block_size /= 256 .and. block_size /= 512 .and. block_size /= 1024) cycle
        write(*,'("block size: ",I0)') block_size

        start_time = omp_get_wtime()
        do i = 1, repeat
          call dense_kernel_k1(input, dense, output_k1, ncols, batch_size, offset, block_size)
        end do
        elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
        write(*,'("Average execution time of dense embedding kernel (k1): ",F0.6," (us)")') elapsed_us
        !$omp target update from(output_k1)

        start_time = omp_get_wtime()
        do i = 1, repeat
          call dense_kernel_k2(input, dense, output_k2, ncols, batch_size, offset, block_size)
        end do
        elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
        write(*,'("Average execution time of dense embedding kernel (k2): ",F0.6," (us)")') elapsed_us
        !$omp target update from(output_k2)

        start_time = omp_get_wtime()
        do i = 1, repeat
          call dense_kernel_flat(input, dense, output_k3, ncols, batch_size, offset, block_size)
        end do
        elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat, real64)
        write(*,'("Average execution time of dense embedding kernel (k3): ",F0.6," (us)")') elapsed_us
        !$omp target update from(output_k3)

        ok = .true.
        do i = 0, input_size - 1
          if (abs(output_k1(i) - output_ref(i)) > 1.0e-3_real32 .or. &
              abs(output_k2(i) - output_ref(i)) > 1.0e-3_real32 .or. &
              abs(output_k3(i) - output_ref(i)) > 1.0e-3_real32) then
            ok = .false.
            exit
          end if
        end do
        write(*,'(A)') merge("PASS", "FAIL", ok)
      end do
      !$omp end target data

      deallocate(input, dense, output_k1, output_k2, output_k3, output_ref, offset)
    end do
  end subroutine run_sweep

  subroutine dense_reference(input, dense, output, embedding_dim, batch_size, offset)
    real(real32), intent(in) :: input(0:), dense(0:)
    real(real32), intent(inout) :: output(0:)
    integer, intent(in) :: embedding_dim, batch_size, offset(0:)
    integer :: batch_idx, idx, nested_idx, range
    real(real32) :: dense_elem

    do batch_idx = 0, batch_size - 1
      range = offset(batch_idx + 1) - offset(batch_idx)
      do idx = 0, embedding_dim - 1
        dense_elem = dense(batch_idx * embedding_dim + idx)
        do nested_idx = idx, range - 1, embedding_dim
          output(offset(batch_idx) + nested_idx) = input(offset(batch_idx) + nested_idx) + dense_elem
        end do
      end do
    end do
  end subroutine dense_reference

  subroutine dense_kernel_k1(input, dense, output, embedding_dim, batch_size, offset, block_size)
    real(real32), intent(in) :: input(0:), dense(0:)
    real(real32), intent(inout) :: output(0:)
    integer, intent(in) :: embedding_dim, batch_size, offset(0:), block_size
    integer :: batch_idx, grain_size, idx, nested_idx, range, tid
    real(real32) :: dense_elem

    !$omp target teams num_teams(batch_size) private(batch_idx, grain_size, idx, nested_idx, range, tid, dense_elem)
    !$omp parallel num_threads(block_size)
    batch_idx = omp_get_team_num()
    grain_size = omp_get_num_threads()
    tid = omp_get_thread_num()
    range = offset(batch_idx + 1) - offset(batch_idx)
    do idx = tid, embedding_dim - 1, grain_size
      dense_elem = dense(batch_idx * embedding_dim + idx)
      do nested_idx = idx, range - 1, embedding_dim
        output(offset(batch_idx) + nested_idx) = input(offset(batch_idx) + nested_idx) + dense_elem
      end do
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine dense_kernel_k1

  subroutine dense_kernel_k2(input, dense, output, embedding_dim, batch_size, offset, block_size)
    real(real32), intent(in) :: input(0:), dense(0:)
    real(real32), intent(inout) :: output(0:)
    integer, intent(in) :: embedding_dim, batch_size, offset(0:), block_size
    integer :: batch_idx, idx, nested_idx, range, start_idx
    real(real32) :: dense_elem

    !$omp target teams num_teams(batch_size) private(batch_idx, idx, nested_idx, range, start_idx, dense_elem)
    !$omp parallel num_threads(block_size)
    batch_idx = omp_get_team_num()
    start_idx = offset(batch_idx)
    range = offset(batch_idx + 1) - start_idx
    do idx = omp_get_thread_num(), embedding_dim - 1, omp_get_num_threads()
      dense_elem = dense(batch_idx * embedding_dim + idx)
      do nested_idx = idx, range - 1, embedding_dim
        output(start_idx + nested_idx) = input(start_idx + nested_idx) + dense_elem
      end do
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine dense_kernel_k2

  subroutine dense_kernel_flat(input, dense, output, embedding_dim, batch_size, offset, block_size)
    real(real32), intent(in) :: input(0:), dense(0:)
    real(real32), intent(inout) :: output(0:)
    integer, intent(in) :: embedding_dim, batch_size, offset(0:), block_size
    integer :: batch_idx, idx, start_idx, range

    !$omp target teams distribute num_teams(batch_size) private(start_idx, range, idx)
    do batch_idx = 0, batch_size - 1
      start_idx = offset(batch_idx)
      range = offset(batch_idx + 1) - start_idx
      !$omp parallel do num_threads(block_size)
      do idx = 0, range - 1
        output(start_idx + idx) = input(start_idx + idx) + dense(batch_idx * embedding_dim + mod(idx, embedding_dim))
      end do
      !$omp end parallel do
    end do
    !$omp end target teams distribute
  end subroutine dense_kernel_flat

end program dense_embedding_main
