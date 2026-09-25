! SPDX-License-Identifier: CC0-1.0
module f16atomic_mod
  use iso_fortran_env, only: int16, int32, int64, real16, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256

contains

  subroutine f16_atomic_on_global_mem(result, n, num_teams)
    real(real16), intent(inout) :: result(0:)
    integer, intent(in) :: n, num_teams
    integer :: tid, i
    real(real16), parameter :: zero_fp16 = 0.0_real16
    real(real16), parameter :: one_fp16 = 1.0_real16

    !$omp target teams distribute parallel do num_teams(num_teams) thread_limit(block_size) private(i)
    do tid = 0, n - 1
      i = mod(tid, block_size)
      !$omp atomic update
      result(2 * i) = result(2 * i) + zero_fp16
      !$omp atomic update
      result(2 * i + 1) = result(2 * i + 1) + one_fp16
    end do
    !$omp end target teams distribute parallel do
  end subroutine f16_atomic_on_global_mem

  integer(int32) function fp16_bits(x) result(bits)
    real(real16), intent(in) :: x
    integer(int16) :: raw

    raw = transfer(x, raw)
    bits = iand(int(raw, int32), int(z'ffff', int32))
  end function fp16_bits

  subroutine atomic_cost(nelems, repeat)
    integer, intent(in) :: nelems, repeat
    integer, parameter :: result_size = block_size * 2
    real(real16), allocatable :: result(:)
    integer :: i, num_teams
    real(real64) :: start_time, elapsed_ns

    allocate(result(0:result_size - 1))

    !$omp target data map(alloc: result(0:result_size - 1))
    !$omp target teams distribute parallel do
    do i = 0, result_size - 1
      result(i) = 0.0_real16
    end do
    !$omp end target teams distribute parallel do

    num_teams = (nelems / 2 + block_size - 1) / block_size

    call f16_atomic_on_global_mem(result, nelems / 2, num_teams)
    !$omp target update from(result(0:result_size - 1))

    write(*,'(A,Z4.4,1X,A,Z4.4)') 'Print the first two elements in HEX: 0x', &
      fp16_bits(result(0)), '0x', fp16_bits(result(1))
    write(*,'(A,F8.6,1X,F0.6,/)') 'Print the first two elements in FLOAT32: ', &
      real(result(0), real32), real(result(1), real32)

    start_time = omp_get_wtime()
    do i = 1, repeat
      call f16_atomic_on_global_mem(result, nelems / 2, num_teams)
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,F12.6,A)') &
      'Average execution time of 16-bit floating-point atomic add on global memory: ', &
      elapsed_ns * 1.0e-3_real64 / real(repeat, real64), ' (us)'
    !$omp end target data

    deallocate(result)
  end subroutine atomic_cost

end module f16atomic_mod

program main
  use f16atomic_mod
  implicit none

  integer :: argc, nelems, repeat, status
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage: ', trim(prog), ' <N> <repeat>'
    write(*,'(A)') 'N: total number of elements (a multiple of 2)'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) nelems
  if (status /= 0) stop 1
  call get_command_argument(2, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  if (nelems <= 0 .or. mod(nelems, 2) /= 0) stop 1

  write(*,'(/,A)') 'FP16 atomic add'
  call atomic_cost(nelems, repeat)
end program main
