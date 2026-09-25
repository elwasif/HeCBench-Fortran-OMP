! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: vector_size = 1024 * 1024
  integer, parameter :: total_iterations = 1024
  integer, parameter :: block_size = 256
  character(len=256) :: arg0, arg
  integer :: repeat
  real(real32), allocatable :: g_data(:, :)
  real(real64), allocatable :: c(:)
  real(real64) :: time_ms, time_ns, checksum
  integer(int64) :: datasize, operations_bytes, operations_128bit

  call get_command_argument(0, arg0)
  write(*,'(A)') 'Shared memory bandwidth microbenchmark'

  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) stop 1

  datasize = int(vector_size, 8) * 8_8
  write(*,'(A,I0,A)') 'Buffer sizes: ', datasize / (1024_8 * 1024_8), 'MB'

  allocate(g_data(4, vector_size / 2))
  g_data = 0.0_real32

  call shmembenchGPU(g_data, int(vector_size, int64), repeat, time_ns)

  time_ms = time_ns * 1.0e-6_real64
  write(*,'(A,F0.6,A)') 'Average kernel execution time : ', time_ms, ' (ms)'

  allocate(c(vector_size))
  c = transfer(g_data, c, vector_size)
  checksum = sum(c)
  if (checksum /= 21256458760384741137729978368.00_real64) then
    write(*,'(A)') 'checksum failed'
  end if

  write(*,'(A)') 'Memory throughput'

  operations_bytes = (6_8 + 4_8 * 5_8 * int(total_iterations, 8) + 6_8) * &
    int(vector_size, 8) * 4_8
  operations_128bit = (6_8 + 4_8 * 5_8 * int(total_iterations, 8) + 6_8) * &
    int(vector_size, 8) / 4_8

  write(*,'(A,F8.2,A,F6.2,A)') achar(9)//'using 128bit operations : ', &
    real(operations_bytes, real64) / time_ns, ' GB/sec (', &
    real(operations_128bit, real64) / time_ns, ' billion accesses/sec)'

  deallocate(g_data)
  deallocate(c)

contains

  subroutine shmembenchGPU(g_data, size, n, time_shmem_128b)
    real(real32), intent(out) :: g_data(:, :)
    integer(int64), intent(in) :: size
    integer, intent(in) :: n
    real(real64), intent(out) :: time_shmem_128b
    integer :: total_blocks, iter, tid, blk, gid, globaltid, j
    real(real64) :: start_time
    real(real32) :: shm_buffer(4, block_size * 6)
    real(real32) :: tmp(4)

    total_blocks = int(size / int(block_size, int64))

    !$omp target data map(from: g_data)
    start_time = omp_get_wtime()
    do iter = 1, n
      !$omp target teams num_teams(total_blocks / 4) thread_limit(block_size) private(shm_buffer, tmp)
      !$omp parallel private(tid, blk, gid, globaltid, j, tmp)
      tid = omp_get_thread_num()
      blk = omp_get_num_threads()
      gid = omp_get_team_num()
      globaltid = gid * blk + tid + 1

      call init_val(tid, shm_buffer(:, tid + 0 * blk + 1))
      call init_val(tid + 1, shm_buffer(:, tid + 1 * blk + 1))
      call init_val(tid + 3, shm_buffer(:, tid + 2 * blk + 1))
      call init_val(tid + 7, shm_buffer(:, tid + 3 * blk + 1))
      call init_val(tid + 13, shm_buffer(:, tid + 4 * blk + 1))
      call init_val(tid + 17, shm_buffer(:, tid + 5 * blk + 1))

      !$omp barrier

      do j = 1, total_iterations
        call shmem_swap(shm_buffer(:, tid + 0 * blk + 1), shm_buffer(:, tid + 1 * blk + 1), tmp)
        call shmem_swap(shm_buffer(:, tid + 2 * blk + 1), shm_buffer(:, tid + 3 * blk + 1), tmp)
        call shmem_swap(shm_buffer(:, tid + 4 * blk + 1), shm_buffer(:, tid + 5 * blk + 1), tmp)

        !$omp barrier

        call shmem_swap(shm_buffer(:, tid + 1 * blk + 1), shm_buffer(:, tid + 2 * blk + 1), tmp)
        call shmem_swap(shm_buffer(:, tid + 3 * blk + 1), shm_buffer(:, tid + 4 * blk + 1), tmp)

        !$omp barrier
      end do

      call reduce_vector(shm_buffer(:, tid + 0 * blk + 1), &
        shm_buffer(:, tid + 1 * blk + 1), shm_buffer(:, tid + 2 * blk + 1), &
        shm_buffer(:, tid + 3 * blk + 1), shm_buffer(:, tid + 4 * blk + 1), &
        shm_buffer(:, tid + 5 * blk + 1), g_data(:, globaltid))
      !$omp end parallel
      !$omp end target teams
    end do
    time_shmem_128b = ((omp_get_wtime() - start_time) * 1.0e9_real64) / real(n, real64)
    !$omp end target data
  end subroutine shmembenchGPU

  subroutine shmem_swap(v1, v2, tmp)
    real(real32), intent(inout) :: v1(4), v2(4)
    real(real32), intent(out) :: tmp(4)
    tmp = v2
    v2 = v1
    v1 = tmp
  end subroutine shmem_swap

  subroutine init_val(i, v)
    integer, intent(in) :: i
    real(real32), intent(out) :: v(4)

    v(1) = real(i, real32)
    v(2) = real(i + 11, real32)
    v(3) = real(i + 19, real32)
    v(4) = real(i + 23, real32)
  end subroutine init_val

  subroutine reduce_vector(v1, v2, v3, v4, v5, v6, v)
    real(real32), intent(in) :: v1(4), v2(4), v3(4), v4(4), v5(4), v6(4)
    real(real32), intent(out) :: v(4)

    v(1) = v1(1) + v2(1) + v3(1) + v4(1) + v5(1) + v6(1)
    v(2) = v1(2) + v2(2) + v3(2) + v4(2) + v5(2) + v6(2)
    v(3) = v1(3) + v2(3) + v3(3) + v4(3) + v5(3) + v6(3)
    v(4) = v1(4) + v2(4) + v3(4) + v4(4) + v5(4) + v6(4)
  end subroutine reduce_vector
end program main
