! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: data_len = 1024
  integer, parameter :: n_threads = 3 * 4 * 7 * 8 * 9 * block_size
  integer :: repeat

  if (command_argument_count() /= 1) then
    print '(A)', 'Usage: ./main <repeat>'
    stop 1
  end if

  repeat = read_arg(1)

  print '(A)', ''
  print '(A)', 'FP64 atomic add'
  call atomic_perf_real64(n_threads, data_len, repeat)

  print '(A)', ''
  print '(A)', 'INT32 atomic add'
  call atomic_perf_int32(n_threads, data_len, repeat)

  print '(A)', ''
  print '(A)', 'FP32 atomic add'
  call atomic_perf_real32(n_threads, data_len, repeat)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine atomic_perf_real64(n, t, repeat)
    integer, intent(in) :: n, t, repeat
    real(real64), allocatable :: data(:), h_data(:), r_data(:)
    integer :: i, rep, offset
    real(real64) :: start_time, end_time

    allocate(data(t), h_data(t), r_data(t))
    do i = 1, t
      h_data(i) = real(mod(i - 1, 1024) + 1, real64)
      data(i) = h_data(i)
    end do

    !$omp target data map(alloc: data(1:t))
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call block_global_real64(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of BlockRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call block_global_ref_real64(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call warp_global_real64(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of WarpRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call warp_global_ref_real64(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_real64(data, offset, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of SingleRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_ref_real64(r_data, offset, n)
    end do
    call print_status(all(data == r_data))

    call time_shared_real64(data, h_data, t, n, repeat, 'BlockRangeAtomicOnSharedMem', 1)
    call time_shared_real64(data, h_data, t, n, repeat, 'WarpRangeAtomicOnSharedMem', 2)
    call time_shared_real64(data, h_data, t, n, repeat, 'SingleRangeAtomicOnSharedMem', 3)
    !$omp end target data
    deallocate(data, h_data, r_data)
  end subroutine atomic_perf_real64

  subroutine atomic_perf_real32(n, t, repeat)
    integer, intent(in) :: n, t, repeat
    real(real32), allocatable :: data(:), h_data(:), r_data(:)
    integer :: i, rep, offset
    real(real64) :: start_time, end_time

    allocate(data(t), h_data(t), r_data(t))
    do i = 1, t
      h_data(i) = real(mod(i - 1, 1024) + 1, real32)
      data(i) = h_data(i)
    end do

    !$omp target data map(alloc: data(1:t))
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call block_global_real32(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of BlockRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call block_global_ref_real32(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call warp_global_real32(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of WarpRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call warp_global_ref_real32(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_real32(data, offset, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of SingleRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_ref_real32(r_data, offset, n)
    end do
    call print_status(all(data == r_data))

    call time_shared_real32(data, h_data, t, n, repeat, 'BlockRangeAtomicOnSharedMem', 1)
    call time_shared_real32(data, h_data, t, n, repeat, 'WarpRangeAtomicOnSharedMem', 2)
    call time_shared_real32(data, h_data, t, n, repeat, 'SingleRangeAtomicOnSharedMem', 3)
    !$omp end target data
    deallocate(data, h_data, r_data)
  end subroutine atomic_perf_real32

  subroutine atomic_perf_int32(n, t, repeat)
    integer, intent(in) :: n, t, repeat
    integer(int32), allocatable :: data(:), h_data(:), r_data(:)
    integer :: i, rep, offset
    real(real64) :: start_time, end_time

    allocate(data(t), h_data(t), r_data(t))
    do i = 1, t
      h_data(i) = int(mod(i - 1, 1024) + 1, int32)
      data(i) = h_data(i)
    end do

    !$omp target data map(alloc: data(1:t))
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call block_global_int32(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of BlockRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call block_global_ref_int32(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call warp_global_int32(data, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of WarpRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      call warp_global_ref_int32(r_data, n)
    end do
    call print_status(all(data == r_data))

    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_int32(data, offset, n)
    end do
    end_time = omp_get_wtime()
    print '(A,F0.6,A)', 'Average execution time of SingleRangeAtomicOnGlobalMem: ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    r_data = h_data
    do rep = 1, repeat
      offset = mod(rep - 1, block_size)
      call single_global_ref_int32(r_data, offset, n)
    end do
    call print_status(all(data == r_data))

    call time_shared_int32(data, h_data, t, n, repeat, 'BlockRangeAtomicOnSharedMem', 1)
    call time_shared_int32(data, h_data, t, n, repeat, 'WarpRangeAtomicOnSharedMem', 2)
    call time_shared_int32(data, h_data, t, n, repeat, 'SingleRangeAtomicOnSharedMem', 3)
    !$omp end target data
    deallocate(data, h_data, r_data)
  end subroutine atomic_perf_int32

  subroutine block_global_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_global_real64

  subroutine warp_global_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_global_real64

  subroutine single_global_real64(data, offset, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_global_real64

  subroutine block_global_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_global_real32

  subroutine warp_global_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_global_real32

  subroutine single_global_real32(data, offset, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      data(idx) = data(idx) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_global_real32

  subroutine block_global_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_global_int32

  subroutine warp_global_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      data(idx) = data(idx) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_global_int32

  subroutine single_global_int32(data, offset, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    !$omp target teams distribute parallel do thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      data(idx) = data(idx) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_global_int32

  subroutine time_shared_real64(data, h_data, t, n, repeat, label, mode)
    real(real64), intent(inout) :: data(:)
    real(real64), intent(in) :: h_data(:)
    integer, intent(in) :: t, n, repeat, mode
    character(len=*), intent(in) :: label
    integer :: rep, offset
    real(real64) :: start_time, end_time
    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      select case (mode)
      case (1)
        call block_shared_real64(data, n)
      case (2)
        call warp_shared_real64(data, n)
      case (3)
        offset = mod(rep - 1, block_size)
        call single_shared_real64(data, offset, n)
      end select
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    call print_status(all(data == h_data))
  end subroutine time_shared_real64

  subroutine time_shared_real32(data, h_data, t, n, repeat, label, mode)
    real(real32), intent(inout) :: data(:)
    real(real32), intent(in) :: h_data(:)
    integer, intent(in) :: t, n, repeat, mode
    character(len=*), intent(in) :: label
    integer :: rep, offset
    real(real64) :: start_time, end_time
    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      select case (mode)
      case (1)
        call block_shared_real32(data, n)
      case (2)
        call warp_shared_real32(data, n)
      case (3)
        offset = mod(rep - 1, block_size)
        call single_shared_real32(data, offset, n)
      end select
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    call print_status(all(data == h_data))
  end subroutine time_shared_real32

  subroutine time_shared_int32(data, h_data, t, n, repeat, label, mode)
    integer(int32), intent(inout) :: data(:)
    integer(int32), intent(in) :: h_data(:)
    integer, intent(in) :: t, n, repeat, mode
    character(len=*), intent(in) :: label
    integer :: rep, offset
    real(real64) :: start_time, end_time
    data = h_data
    !$omp target update to(data(1:t))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      select case (mode)
      case (1)
        call block_shared_int32(data, n)
      case (2)
        call warp_shared_int32(data, n)
      case (3)
        offset = mod(rep - 1, block_size)
        call single_shared_int32(data, offset, n)
      end select
    end do
    end_time = omp_get_wtime()
    print '(A,A,A,F0.6,A)', 'Average execution time of ', label, ': ', &
      (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    !$omp target update from(data(1:t))
    call print_status(all(data == h_data))
  end subroutine time_shared_int32

  subroutine block_shared_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    real(real64) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_shared_real64

  subroutine warp_shared_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    real(real64) :: smem_data(32)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_shared_real64

  subroutine single_shared_real64(data, offset, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i
    real(real64) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      smem_data(offset + 1) = smem_data(offset + 1) + 1.0_real64
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_shared_real64

  subroutine block_shared_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    real(real32) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_shared_real32

  subroutine warp_shared_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    real(real32) :: smem_data(32)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_shared_real32

  subroutine single_shared_real32(data, offset, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i
    real(real32) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      smem_data(offset + 1) = smem_data(offset + 1) + 1.0_real32
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_shared_real32

  subroutine block_shared_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    integer(int32) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = mod(i, block_size) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine block_shared_int32

  subroutine warp_shared_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i, idx
    integer(int32) :: smem_data(32)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size) private(idx)
    do i = 0, n - 1
      idx = iand(i, 31) + 1
      !$omp atomic update
      smem_data(idx) = smem_data(idx) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine warp_shared_int32

  subroutine single_shared_int32(data, offset, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i
    integer(int32) :: smem_data(block_size)
    !$omp target teams distribute parallel do num_teams(n / block_size) thread_limit(block_size)
    do i = 0, n - 1
      !$omp atomic update
      smem_data(offset + 1) = smem_data(offset + 1) + 1_int32
    end do
    !$omp end target teams distribute parallel do
  end subroutine single_shared_int32

  subroutine block_global_ref_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(mod(i, block_size) + 1) = data(mod(i, block_size) + 1) + 1.0_real64
    end do
  end subroutine block_global_ref_real64

  subroutine warp_global_ref_real64(data, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(iand(i, 31) + 1) = data(iand(i, 31) + 1) + 1.0_real64
    end do
  end subroutine warp_global_ref_real64

  subroutine single_global_ref_real64(data, offset, n)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    do i = 0, n - 1
      data(idx) = data(idx) + 1.0_real64
    end do
  end subroutine single_global_ref_real64

  subroutine block_global_ref_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(mod(i, block_size) + 1) = data(mod(i, block_size) + 1) + 1.0_real32
    end do
  end subroutine block_global_ref_real32

  subroutine warp_global_ref_real32(data, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(iand(i, 31) + 1) = data(iand(i, 31) + 1) + 1.0_real32
    end do
  end subroutine warp_global_ref_real32

  subroutine single_global_ref_real32(data, offset, n)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    do i = 0, n - 1
      data(idx) = data(idx) + 1.0_real32
    end do
  end subroutine single_global_ref_real32

  subroutine block_global_ref_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(mod(i, block_size) + 1) = data(mod(i, block_size) + 1) + 1_int32
    end do
  end subroutine block_global_ref_int32

  subroutine warp_global_ref_int32(data, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n - 1
      data(iand(i, 31) + 1) = data(iand(i, 31) + 1) + 1_int32
    end do
  end subroutine warp_global_ref_int32

  subroutine single_global_ref_int32(data, offset, n)
    integer(int32), intent(inout) :: data(:)
    integer, intent(in) :: offset, n
    integer :: i, idx
    idx = offset + 1
    do i = 0, n - 1
      data(idx) = data(idx) + 1_int32
    end do
  end subroutine single_global_ref_int32

  subroutine print_status(pass)
    logical, intent(in) :: pass
    if (pass) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine print_status

end program main
