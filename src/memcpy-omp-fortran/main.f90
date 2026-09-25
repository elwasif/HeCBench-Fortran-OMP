! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: num_size = 16
  integer(int64) :: sizes(num_size)
  character(len=256) :: arg0, arg1
  integer :: argc
  integer(int32) :: repeat

  argc = command_argument_count()
  call get_command_argument(0, arg0)
  if (argc /= 1) then
    write(*,'(3A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *) repeat

  call setup(sizes)

  call run_benchmark(sizes, repeat)

contains

  subroutine setup(sizes)
    integer(int64), intent(out) :: sizes(num_size)
    integer :: i

    do i = 1, num_size
      sizes(i) = ishft(1_int64, i + 5)
    end do
  end subroutine setup

  subroutine val_set(a, val)
    integer(int32), intent(out) :: a(:)
    integer(int32), intent(in) :: val

    a = val
  end subroutine val_set

  subroutine run_benchmark(sizes, repeat)
    integer(int64), intent(in) :: sizes(num_size)
    integer(int32), intent(in) :: repeat
    integer :: i, j, len
    character(len=32) :: gap_label
    integer(int32), allocatable :: a(:)
    real(real64) :: start_time, end_time, time_h2d_ns, time_d2h_ns, gap_ns_per_byte

    do i = 1, num_size
      len = int(sizes(i) / 4_int64)
      allocate(a(len))
      call val_set(a, 1_int32)

      !$omp target data map(alloc: a(1:len))
      do j = 1, repeat
        !$omp target update to(a(1:len))
      end do

      start_time = omp_get_wtime()
      do j = 1, repeat
        !$omp target update to(a(1:len))
      end do
      end_time = omp_get_wtime()
      time_h2d_ns = (end_time - start_time) * 1.0d9
      write(*,'("Copy ",I0," bytes from host to device takes ",F0.6," us")') &
        sizes(i), (time_h2d_ns * 1.0d-3) / real(repeat, kind=real64)

      do j = 1, repeat
        !$omp target update from(a(1:len))
      end do

      start_time = omp_get_wtime()
      do j = 1, repeat
        !$omp target update from(a(1:len))
      end do
      end_time = omp_get_wtime()
      time_d2h_ns = (end_time - start_time) * 1.0d9
      write(*,'("Copy ",I0," bytes from device to host takes ",F0.6," us")') &
        sizes(i), (time_d2h_ns * 1.0d-3) / real(repeat, kind=real64)

      gap_ns_per_byte = abs(time_h2d_ns - time_d2h_ns) / &
        (real(repeat, kind=real64) * real(sizes(i), kind=real64))
      write(gap_label, '(F12.6)') gap_ns_per_byte
      gap_label = adjustl(gap_label)
      write(*,'("Timing gap in nanoseconds per byte: ",A)') trim(gap_label)

      !$omp end target data

      deallocate(a)
      write(*,*)
    end do
  end subroutine run_benchmark

end program main
