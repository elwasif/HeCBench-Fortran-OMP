! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() result(value) bind(C, name='rand')
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  integer, parameter :: max_distance = 200
  integer :: num_nodes, num_iterations, block_size, block_threads
  integer :: matrix_size, iter, k, x, y, idx, yk_idx, kx_idx
  integer :: mismatches
  integer(int32), allocatable :: pathDistanceMatrix(:), pathMatrix(:)
  integer(int32), allocatable :: verificationPathDistanceMatrix(:), verificationPathMatrix(:)
  integer(int32) :: distance_yx, distance_yk, distance_kx, indirect
  real(real64) :: total_time, start_time, end_time
  character(len=64) :: time_text

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <number of nodes> <iterations> <block size>'
    stop 1
  end if

  num_nodes = read_arg(1)
  num_iterations = read_arg(2)
  block_size = read_arg(3)
  if (num_nodes <= 0 .or. num_iterations <= 0 .or. block_size <= 0) then
    error stop 'floydwarshall arguments must be positive'
  end if

  if (mod(num_nodes, block_size) /= 0) then
    num_nodes = (num_nodes / block_size + 1) * block_size
  end if
  block_threads = block_size * block_size
  matrix_size = num_nodes * num_nodes

  allocate(pathDistanceMatrix(matrix_size), pathMatrix(matrix_size))
  allocate(verificationPathDistanceMatrix(matrix_size), verificationPathMatrix(matrix_size))
  call initialize_matrix(pathDistanceMatrix, num_nodes)
  call initialize_path_matrix(pathMatrix, num_nodes)
  verificationPathDistanceMatrix = pathDistanceMatrix
  verificationPathMatrix = pathMatrix

  total_time = 0.0_real64
  !$omp target data map(alloc: pathDistanceMatrix(1:matrix_size), pathMatrix(1:matrix_size))
  do iter = 1, num_iterations
    !$omp target update to(pathDistanceMatrix(1:matrix_size))
    start_time = omp_get_wtime()
    do k = 1, num_nodes
      !$omp target teams distribute parallel do collapse(2) thread_limit(block_threads) &
      !$omp& private(idx, yk_idx, kx_idx, distance_yx, distance_yk, distance_kx, indirect)
      do y = 1, num_nodes
        do x = 1, num_nodes
          idx = (y - 1) * num_nodes + x
          yk_idx = (y - 1) * num_nodes + k
          kx_idx = (k - 1) * num_nodes + x
          distance_yx = pathDistanceMatrix(idx)
          distance_yk = pathDistanceMatrix(yk_idx)
          distance_kx = pathDistanceMatrix(kx_idx)
          indirect = distance_yk + distance_kx
          if (indirect < distance_yx) then
            pathDistanceMatrix(idx) = indirect
            pathMatrix(idx) = k - 1
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    total_time = total_time + end_time - start_time
  end do
  !$omp target update from(pathDistanceMatrix(1:matrix_size))
  !$omp end target data

  call floyd_warshall_cpu(verificationPathDistanceMatrix, verificationPathMatrix, num_nodes)
  mismatches = count(pathDistanceMatrix /= verificationPathDistanceMatrix)

  write(time_text, '(F0.6)') total_time / real(num_iterations, real64)
  if (time_text(1:1) == '.') time_text = '0' // trim(time_text)
  write(*, '(A,A,A)') 'Average kernel execution time ', trim(time_text), ' (s)'
  if (mismatches == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
    if (num_nodes <= 8) then
      call print_debug(verificationPathDistanceMatrix, pathDistanceMatrix, num_nodes)
    end if
  end if

  deallocate(pathDistanceMatrix, pathMatrix)
  deallocate(verificationPathDistanceMatrix, verificationPathMatrix)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_matrix(matrix, num_nodes)
    integer(int32), intent(out) :: matrix(:)
    integer, intent(in) :: num_nodes
    integer :: i, j, idx
    call c_srand(2_c_int)
    do i = 1, num_nodes
      do j = 1, num_nodes
        idx = (i - 1) * num_nodes + j
        matrix(idx) = int(mod(int(c_rand(), int64), int(max_distance + 1, int64)), int32)
      end do
    end do
    do i = 1, num_nodes
      matrix((i - 1) * num_nodes + i) = 0_int32
    end do
  end subroutine initialize_matrix

  subroutine initialize_path_matrix(pathMatrix, num_nodes)
    integer(int32), intent(out) :: pathMatrix(:)
    integer, intent(in) :: num_nodes
    integer :: i, j
    do i = 1, num_nodes
      do j = 1, i - 1
        pathMatrix((i - 1) * num_nodes + j) = int(i - 1, int32)
        pathMatrix((j - 1) * num_nodes + i) = int(j - 1, int32)
      end do
      pathMatrix((i - 1) * num_nodes + i) = int(i - 1, int32)
    end do
  end subroutine initialize_path_matrix

  subroutine floyd_warshall_cpu(pathDistanceMatrix, pathMatrix, num_nodes)
    integer(int32), intent(inout) :: pathDistanceMatrix(:), pathMatrix(:)
    integer, intent(in) :: num_nodes
    integer :: k, x, y, idx, yk_idx, kx_idx
    integer(int32) :: indirect
    do k = 1, num_nodes
      do y = 1, num_nodes
        yk_idx = (y - 1) * num_nodes + k
        do x = 1, num_nodes
          idx = (y - 1) * num_nodes + x
          kx_idx = (k - 1) * num_nodes + x
          indirect = pathDistanceMatrix(yk_idx) + pathDistanceMatrix(kx_idx)
          if (indirect < pathDistanceMatrix(idx)) then
            pathDistanceMatrix(idx) = indirect
            pathMatrix(idx) = k - 1
          end if
        end do
      end do
    end do
  end subroutine floyd_warshall_cpu

  subroutine print_debug(reference, device, num_nodes)
    integer(int32), intent(in) :: reference(:), device(:)
    integer, intent(in) :: num_nodes
    integer :: i, j
    do i = 1, num_nodes
      do j = 1, num_nodes
        write(*, '(A,I0,1X)', advance='no') 'host: ', reference((i - 1) * num_nodes + j)
      end do
      print '(A)'
    end do
    do i = 1, num_nodes
      do j = 1, num_nodes
        write(*, '(A,I0,1X)', advance='no') 'device: ', device((i - 1) * num_nodes + j)
      end do
      print '(A)'
    end do
  end subroutine print_debug

end program main
