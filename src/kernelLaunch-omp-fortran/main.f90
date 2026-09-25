! SPDX-License-Identifier: CC0-1.0
module kernel_launch_kernels
  use, intrinsic :: iso_c_binding, only : c_signed_char
  use omp_lib
  implicit none

  type :: SmallKernelArgs
    integer(c_signed_char) :: args(16)
  end type SmallKernelArgs

  type :: MediumKernelArgs
    integer(c_signed_char) :: args(256)
  end type MediumKernelArgs

  type :: LargeKernelArgs
    integer(c_signed_char) :: args(4096)
  end type LargeKernelArgs

contains

  subroutine KernelWithSmallArgs(args)
    type(SmallKernelArgs), intent(in) :: args
    integer :: i

    i = omp_get_num_teams() * omp_get_num_threads() + omp_get_thread_num()
    if (i < 0) call consume_small(args)
  end subroutine KernelWithSmallArgs

  subroutine KernelWithMediumArgs(args)
    type(MediumKernelArgs), intent(in) :: args
    integer :: i

    i = omp_get_num_teams() * omp_get_num_threads() + omp_get_thread_num()
    if (i < 0) call consume_medium(args)
  end subroutine KernelWithMediumArgs

  subroutine KernelWithLargeArgs(args)
    type(LargeKernelArgs), intent(in) :: args
    integer :: i

    i = omp_get_num_teams() * omp_get_num_threads() + omp_get_thread_num()
    if (i < 0) call consume_large(args)
  end subroutine KernelWithLargeArgs

  subroutine consume_small(args)
    type(SmallKernelArgs), intent(in) :: args
    integer(c_signed_char) :: tmp

    tmp = args%args(1)
  end subroutine consume_small

  subroutine consume_medium(args)
    type(MediumKernelArgs), intent(in) :: args
    integer(c_signed_char) :: tmp

    tmp = args%args(1)
  end subroutine consume_medium

  subroutine consume_large(args)
    type(LargeKernelArgs), intent(in) :: args
    integer(c_signed_char) :: tmp

    tmp = args%args(1)
  end subroutine consume_large

end module kernel_launch_kernels

program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  use kernel_launch_kernels
  implicit none

  character(len=256) :: arg0, arg1
  integer :: repeat, parse_status
  type(SmallKernelArgs) :: small_kernel_args
  type(MediumKernelArgs) :: medium_kernel_args
  type(LargeKernelArgs) :: large_kernel_args

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    print '(3A)', 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  read(arg1, *, iostat=parse_status) repeat
  if (parse_status /= 0) repeat = 0

  call run_small_kernel_args(small_kernel_args, repeat)
  call run_medium_kernel_args(medium_kernel_args, repeat)
  call run_large_kernel_args(large_kernel_args, repeat)

contains

  subroutine run_small_kernel_args(small_kernel_args, repeat)
    type(SmallKernelArgs), intent(in) :: small_kernel_args
    integer, intent(in) :: repeat
    integer :: i
    real(real64) :: start_time, end_time, elapsed_us

    do i = 1, repeat
      !$omp target map(to: small_kernel_args)
      call KernelWithSmallArgs(small_kernel_args)
      !$omp end target
    end do

    start_time = omp_get_wtime()
    do i = 1, repeat
      !$omp target map(to: small_kernel_args)
      call KernelWithSmallArgs(small_kernel_args)
      !$omp end target
    end do
    end_time = omp_get_wtime()

    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of kernelWithSmallArgs: ', elapsed_us, ' (us)'
  end subroutine run_small_kernel_args

  subroutine run_medium_kernel_args(medium_kernel_args, repeat)
    type(MediumKernelArgs), intent(in) :: medium_kernel_args
    integer, intent(in) :: repeat
    integer :: i
    real(real64) :: start_time, end_time, elapsed_us

    do i = 1, repeat
      !$omp target map(to: medium_kernel_args)
      call KernelWithMediumArgs(medium_kernel_args)
      !$omp end target
    end do

    start_time = omp_get_wtime()
    do i = 1, repeat
      !$omp target map(to: medium_kernel_args)
      call KernelWithMediumArgs(medium_kernel_args)
      !$omp end target
    end do
    end_time = omp_get_wtime()

    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of kernelWithMediumArgs: ', elapsed_us, ' (us)'
  end subroutine run_medium_kernel_args

  subroutine run_large_kernel_args(large_kernel_args, repeat)
    type(LargeKernelArgs), intent(in) :: large_kernel_args
    integer, intent(in) :: repeat
    integer :: i
    real(real64) :: start_time, end_time, elapsed_us

    do i = 1, repeat
      !$omp target map(to: large_kernel_args)
      call KernelWithLargeArgs(large_kernel_args)
      !$omp end target
    end do

    start_time = omp_get_wtime()
    do i = 1, repeat
      !$omp target map(to: large_kernel_args)
      call KernelWithLargeArgs(large_kernel_args)
      !$omp end target
    end do
    end_time = omp_get_wtime()

    elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    print '(A,F0.6,A)', 'Average execution time of kernelWithLargeArgs: ', elapsed_us, ' (us)'
  end subroutine run_large_kernel_args

end program main
