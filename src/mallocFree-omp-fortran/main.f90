! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_null_ptr, c_ptr, c_size_t, c_sizeof
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none
  !$omp requires unified_shared_memory

  integer, parameter :: num_size = 19
  integer, parameter :: num_iter = 500
  character(len=256) :: arg0, arg1
  integer(c_size_t) :: total_global_mem
  integer(c_size_t) :: sizes(num_size)
  type(c_ptr) :: ad(num_iter)
  integer(c_int), allocatable :: a_host(:)
  integer :: num, device_num, i, j
  real(real64) :: start_time, end_time, elapsed_us

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(3A)') 'Usage: ', trim(arg0), ' <total global memory size in bytes>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) total_global_mem

  sizes = 0_c_size_t
  ad = c_null_ptr
  num = num_size
  call setup(sizes, num, a_host, total_global_mem)

  device_num = 0
  call test_init(sizes(1), device_num)

  do i = 1, num
    start_time = omp_get_wtime()
    do j = 1, num_iter
      ad(j) = omp_target_alloc(sizes(i), device_num)
    end do
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(num_iter, real64)
    write(*,'(A,I0,A,F0.6,A)') 'omp_target_alloc(', int(sizes(i), int64), ') takes ', elapsed_us, ' us'

    start_time = omp_get_wtime()
    do j = 1, num_iter
      call omp_target_free(ad(j), device_num)
      ad(j) = c_null_ptr
    end do
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(num_iter, real64)
    write(*,'(A,I0,A,F0.6,A)') 'omp_target_free(', int(sizes(i), int64), ') takes ', elapsed_us, ' us'
  end do

  deallocate(a_host)

contains

  subroutine setup(sizes, num, a_host, total_global_mem)
    integer(c_size_t), intent(out) :: sizes(:)
    integer, intent(inout) :: num
    integer(c_int), allocatable, intent(out) :: a_host(:)
    integer(c_size_t), intent(in) :: total_global_mem
    integer :: i, n

    do i = 1, num
      sizes(i) = int(ishft(1_int64, i + 5), c_size_t)
      if (int(num_iter + 1, c_size_t) * sizes(i) > total_global_mem) then
        num = i - 1
        exit
      end if
    end do
    if (num <= 0) stop 1

    n = int(sizes(num) / c_sizeof(0_c_int))
    allocate(a_host(n))
    a_host = 1_c_int
  end subroutine setup

  subroutine test_init(size_bytes, device_num)
    integer(c_size_t), intent(in) :: size_bytes
    integer, intent(in) :: device_num
    type(c_ptr) :: ad
    real(real64) :: start_time, end_time, elapsed_us

    write(*,'(A)') 'Initial allocation and deallocation'

    start_time = omp_get_wtime()
    ad = omp_target_alloc(size_bytes, device_num)
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64
    write(*,'(A,I0,A,F0.6,A)') 'omp_target_alloc(', int(size_bytes, int64), ') takes ', elapsed_us, ' us'

    start_time = omp_get_wtime()
    call omp_target_free(ad, device_num)
    end_time = omp_get_wtime()
    elapsed_us = (end_time - start_time) * 1.0e6_real64
    write(*,'(A,I0,A,F0.6,A)') 'omp_target_free(', int(size_bytes, int64), ') takes ', elapsed_us, ' us'
    write(*,'(A)') ''
  end subroutine test_init

end program main
