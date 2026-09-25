! SPDX-License-Identifier: CC0-1.0
program main
  use iso_fortran_env, only: real32, real64, int64
  use iso_c_binding, only: c_double, c_float, c_int64_t
  use omp_lib, only: omp_get_wtime
  implicit none

  integer :: argc, repeat
  integer(int64) :: num_elements
  character(len=64) :: arg

  interface
    subroutine dp_fill_real32(srcA, srcB, iNumElements, src_size, dst_ref) bind(C, name="dp_fill_real32")
      import :: c_float, c_int64_t
      real(c_float) :: srcA(*), srcB(*)
      integer(c_int64_t), value :: iNumElements, src_size
      real(c_float) :: dst_ref
    end subroutine dp_fill_real32

    subroutine dp_fill_real64(srcA, srcB, iNumElements, src_size, dst_ref) bind(C, name="dp_fill_real64")
      import :: c_double, c_int64_t
      real(c_double) :: srcA(*), srcB(*)
      integer(c_int64_t), value :: iNumElements, src_size
      real(c_double) :: dst_ref
    end subroutine dp_fill_real64
  end interface

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, arg)
    write(*,'("Usage: ",A," <number of elements> <repeat>")') trim(arg)
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) num_elements
  call get_command_argument(2, arg)
  read(arg, *) repeat

  write(*,'("------------- Data type is Float32 ---------------")')
  call dot_real32(num_elements, repeat)
  write(*,'("------------- Data type is Float64 ---------------")')
  call dot_real64(num_elements, repeat)

contains

  integer(int64) function round_up(local_size, elements) result(global_size)
    integer(int64), intent(in) :: local_size, elements

    global_size = ((elements + local_size - 1_int64) / local_size) * local_size
  end function round_up

  subroutine print_average_kernel_time(milliseconds)
    real(real64), intent(in) :: milliseconds
    character(len=64) :: value

    write(value, '(F0.6)') milliseconds
    if (value(1:1) == '.') then
      write(*,'("Average kernel execution time 0",A," (ms)")') trim(value)
    else
      write(*,'("Average kernel execution time ",A," (ms)")') trim(value)
    end if
  end subroutine print_average_kernel_time

  subroutine dot_real32(num_elements, repeat)
    integer(int64), intent(in) :: num_elements
    integer, intent(in) :: repeat
    integer, parameter :: local_work_size = 256
    integer(int64) :: global_work_size, src_size, i, gid, offset
    real(real32), allocatable :: src_a(:), src_b(:)
    real(real32) :: dst, dst_ref
    integer :: iter
    real(real64) :: start_time, end_time

    global_work_size = round_up(int(local_work_size, int64), num_elements)
    src_size = global_work_size
    write(*,'("Global Work Size ",A,A,"= ",I0)') achar(9), achar(9), global_work_size
    write(*,'("Local Work Size ",A,A,"= ",I0)') achar(9), achar(9), local_work_size

    allocate(src_a(src_size), src_b(src_size))
    call dp_fill_real32(src_a, src_b, num_elements, src_size, dst_ref)

    !$omp target data map(to: src_a(1:src_size), src_b(1:src_size))
      do iter = 1, 100
        dst = 0.0_real32
        !$omp target teams distribute parallel do reduction(+:dst) thread_limit(local_work_size) private(offset)
        do gid = 0, src_size / 4 - 1
          offset = gid * 4 + 1
          dst = dst + src_a(offset) * src_b(offset) + src_a(offset + 1) * src_b(offset + 1) + &
                      src_a(offset + 2) * src_b(offset + 2) + src_a(offset + 3) * src_b(offset + 3)
        end do
        !$omp end target teams distribute parallel do
      end do

      start_time = omp_get_wtime()
      do iter = 1, repeat
        dst = 0.0_real32
        !$omp target teams distribute parallel do reduction(+:dst) thread_limit(local_work_size) private(offset)
        do gid = 0, src_size / 4 - 1
          offset = gid * 4 + 1
          dst = dst + src_a(offset) * src_b(offset) + src_a(offset + 1) * src_b(offset + 1) + &
                      src_a(offset + 2) * src_b(offset + 2) + src_a(offset + 3) * src_b(offset + 3)
        end do
        !$omp end target teams distribute parallel do
      end do
      end_time = omp_get_wtime()
    !$omp end target data

    call print_average_kernel_time((end_time - start_time) * 1000.0_real64 / repeat)
    if (abs(real(dst, real64) - real(dst_ref, real64)) <= 0.0_real64) then
      write(*,'("PASS",/)')
    else
      write(*,'("FAIL",/)')
    end if
    deallocate(src_a, src_b)
  end subroutine dot_real32

  subroutine dot_real64(num_elements, repeat)
    integer(int64), intent(in) :: num_elements
    integer, intent(in) :: repeat
    integer, parameter :: local_work_size = 256
    integer(int64) :: global_work_size, src_size, i, gid, offset
    real(real64), allocatable :: src_a(:), src_b(:)
    real(real64) :: dst, dst_ref, start_time, end_time
    integer :: iter

    global_work_size = round_up(int(local_work_size, int64), num_elements)
    src_size = global_work_size
    write(*,'("Global Work Size ",A,A,"= ",I0)') achar(9), achar(9), global_work_size
    write(*,'("Local Work Size ",A,A,"= ",I0)') achar(9), achar(9), local_work_size

    allocate(src_a(src_size), src_b(src_size))
    call dp_fill_real64(src_a, src_b, num_elements, src_size, dst_ref)

    !$omp target data map(to: src_a(1:src_size), src_b(1:src_size))
      do iter = 1, 100
        dst = 0.0_real64
        !$omp target teams distribute parallel do reduction(+:dst) thread_limit(local_work_size) private(offset)
        do gid = 0, src_size / 4 - 1
          offset = gid * 4 + 1
          dst = dst + src_a(offset) * src_b(offset) + src_a(offset + 1) * src_b(offset + 1) + &
                      src_a(offset + 2) * src_b(offset + 2) + src_a(offset + 3) * src_b(offset + 3)
        end do
        !$omp end target teams distribute parallel do
      end do

      start_time = omp_get_wtime()
      do iter = 1, repeat
        dst = 0.0_real64
        !$omp target teams distribute parallel do reduction(+:dst) thread_limit(local_work_size) private(offset)
        do gid = 0, src_size / 4 - 1
          offset = gid * 4 + 1
          dst = dst + src_a(offset) * src_b(offset) + src_a(offset + 1) * src_b(offset + 1) + &
                      src_a(offset + 2) * src_b(offset + 2) + src_a(offset + 3) * src_b(offset + 3)
        end do
        !$omp end target teams distribute parallel do
      end do
      end_time = omp_get_wtime()
    !$omp end target data

    call print_average_kernel_time((end_time - start_time) * 1000.0_real64 / repeat)
    if (dst == dst_ref) then
      write(*,'("PASS",/)')
    else
      write(*,'("FAIL",/)')
    end if
    deallocate(src_a, src_b)
  end subroutine dot_real64

end program main
