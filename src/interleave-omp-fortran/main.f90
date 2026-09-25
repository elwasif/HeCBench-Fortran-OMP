! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: num_elements = 4096
  integer, parameter :: count = 4096
  integer, parameter :: members = 16

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name='rand') result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: repeat_count
  integer(int32), allocatable :: interleaved_src(:,:), interleaved_dst(:,:)
  integer(int32), allocatable :: non_interleaved_src(:,:), non_interleaved_dst(:,:)

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat_count
  if (repeat_count <= 0) stop 1

  allocate(interleaved_src(members, num_elements), interleaved_dst(members, num_elements))
  allocate(non_interleaved_src(num_elements, members), non_interleaved_dst(num_elements, members))

  call c_srand(1_c_int)
  call initialize(interleaved_src, interleaved_dst, non_interleaved_src, non_interleaved_dst)
  call add_test_non_interleaved(non_interleaved_dst, non_interleaved_src, repeat_count)
  call add_test_interleaved(interleaved_dst, interleaved_src, repeat_count)
  call verify(interleaved_dst, non_interleaved_dst)

  deallocate(interleaved_src, interleaved_dst, non_interleaved_src, non_interleaved_dst)

contains

  subroutine initialize(inter_src, inter_dst, non_src, non_dst)
    integer(int32), intent(out) :: inter_src(:,:), inter_dst(:,:), non_src(:,:), non_dst(:,:)
    integer :: i, field, value

    do i = 1, num_elements
      do field = 1, members
        ! Preserve the C++ original's implicit srand(1)/rand() input stream.
        value = modulo(c_rand(), 16)
        inter_src(field, i) = value
        non_src(i, field) = value
        inter_dst(field, i) = 0
        non_dst(i, field) = 0
      end do
    end do
  end subroutine initialize

  subroutine add_test_interleaved(dst, src, repeat_count)
    integer(int32), intent(inout) :: dst(:,:)
    integer(int32), intent(in) :: src(:,:)
    integer, intent(in) :: repeat_count
    integer :: rep, tid, iter, field
    real(real64) :: start_time, elapsed_s

    !$omp target data map(to: src(1:members,1:num_elements)) map(tofrom: dst(1:members,1:num_elements))
    start_time = omp_get_wtime()
    do rep = 1, repeat_count
      !$omp target teams distribute parallel do thread_limit(256) private(iter, field)
      do tid = 1, num_elements
        do iter = 1, count
          do field = 1, members
            dst(field, tid) = dst(field, tid) + src(field, tid)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    elapsed_s = (omp_get_wtime() - start_time) / real(repeat_count, real64)
    write(*,'(A,F0.6,A)') 'Average kernel (interleaved) execution time ', elapsed_s, ' (s)'
    !$omp end target data
  end subroutine add_test_interleaved

  subroutine add_test_non_interleaved(dst, src, repeat_count)
    integer(int32), intent(inout) :: dst(:,:)
    integer(int32), intent(in) :: src(:,:)
    integer, intent(in) :: repeat_count
    integer :: rep, tid, iter, field
    real(real64) :: start_time, elapsed_s

    !$omp target data map(to: src(1:num_elements,1:members)) map(tofrom: dst(1:num_elements,1:members))
    start_time = omp_get_wtime()
    do rep = 1, repeat_count
      !$omp target teams distribute parallel do thread_limit(256) private(iter, field)
      do tid = 1, num_elements
        do iter = 1, count
          do field = 1, members
            dst(tid, field) = dst(tid, field) + src(tid, field)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    elapsed_s = (omp_get_wtime() - start_time) / real(repeat_count, real64)
    write(*,'(A,F0.6,A)') 'Average kernel (non-interleaved) execution time ', elapsed_s, ' (s)'
    !$omp end target data
  end subroutine add_test_non_interleaved

  subroutine verify(inter_dst, non_dst)
    integer(int32), intent(in) :: inter_dst(:,:), non_dst(:,:)
    integer :: i, field

    do i = 1, num_elements
      do field = 1, members
        if (inter_dst(field, i) /= non_dst(i, field)) then
          write(*,'(A,I0,A,I0)') 'Mismatch at element ', i - 1, ' field ', field - 1
          stop 1
        end if
      end do
    end do
  end subroutine verify

end program main
