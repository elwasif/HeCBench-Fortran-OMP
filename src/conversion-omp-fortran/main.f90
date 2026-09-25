! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_signed_char
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer :: nelems, niters

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(3A)') 'Usage: ', trim(arg0), ' <number of elements> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  nelems = atoi_compat(arg1)
  niters = atoi_compat(arg2)

  write(*,'(A)') 'float -> float'
  call convert_r4_r4(nelems, niters, 4, 4)
  write(*,'(A)') 'float -> int'
  call convert_r4_i4(nelems, niters, 4, 4)
  write(*,'(A)') 'float -> char'
  call convert_r4_i1(nelems, niters, 4, 1)
  write(*,'(A)') 'float -> uchar'
  call convert_r4_u1(nelems, niters, 4, 1)

  write(*,'(A)') 'int -> int'
  call convert_i4_i4(nelems, niters, 4, 4)
  write(*,'(A)') 'int -> float'
  call convert_i4_r4(nelems, niters, 4, 4)
  write(*,'(A)') 'int -> char'
  call convert_i4_i1(nelems, niters, 4, 1)
  write(*,'(A)') 'int -> uchar'
  call convert_i4_u1(nelems, niters, 4, 1)

  write(*,'(A)') 'char -> int'
  call convert_i1_i4(nelems, niters, 1, 4)
  write(*,'(A)') 'char -> float'
  call convert_i1_r4(nelems, niters, 1, 4)
  write(*,'(A)') 'char -> char'
  call convert_i1_i1(nelems, niters, 1, 1)
  write(*,'(A)') 'char -> uchar'
  call convert_i1_u1(nelems, niters, 1, 1)

  write(*,'(A)') 'uchar -> int'
  call convert_u1_i4(nelems, niters, 1, 4)
  write(*,'(A)') 'uchar -> float'
  call convert_u1_r4(nelems, niters, 1, 4)
  write(*,'(A)') 'uchar -> char'
  call convert_u1_i1(nelems, niters, 1, 1)
  write(*,'(A)') 'uchar -> uchar'
  call convert_u1_u1(nelems, niters, 1, 1)

contains

  integer function atoi_compat(arg)
    character(len=*), intent(in) :: arg
    integer :: ios

    read(arg, *, iostat=ios) atoi_compat
    if (ios /= 0) atoi_compat = 0
  end function atoi_compat

  pure integer(c_signed_char) function to_uchar_i4(value)
    integer(int32), intent(in) :: value

    to_uchar_i4 = int(iand(value, int(z'000000ff', int32)), c_signed_char)
  end function to_uchar_i4

  pure integer(c_signed_char) function to_uchar_r4(value)
    real(real32), intent(in) :: value

    to_uchar_r4 = to_uchar_i4(int(value, int32))
  end function to_uchar_r4

  pure integer(int32) function from_uchar(value)
    integer(c_signed_char), intent(in) :: value

    from_uchar = iand(int(value, int32), int(z'000000ff', int32))
  end function from_uchar

  subroutine print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    integer, intent(in) :: src_bytes, dst_bytes, nelems, niters
    real(real64), intent(in) :: start_time, end_time
    real(real64) :: elapsed_sec, size_gb, bandwidth

    elapsed_sec = (end_time - start_time) / real(niters, real64)
    size_gb = real(src_bytes + dst_bytes, real64) * real(nelems, real64) / 1.0e9_real64
    if (elapsed_sec > 0.0_real64) then
      bandwidth = size_gb / elapsed_sec
    else
      bandwidth = 0.0_real64
    end if
    write(*,'(A,F0.2,A,F0.6,A,F0.6)') 'size(GB):', size_gb, &
      ', average time(sec):', elapsed_sec, ', BW:', bandwidth
  end subroutine print_timing

  subroutine convert_r4_r4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    real(real32), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = src(i)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_r4_r4

  subroutine convert_r4_i4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    real(real32), allocatable :: src(:)
    integer(int32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = int(src(i), int32)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = int(src(i), int32)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_r4_i4

  subroutine convert_r4_i1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    real(real32), allocatable :: src(:)
    integer(c_signed_char), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = int(src(i), c_signed_char)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = int(src(i), c_signed_char)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_r4_i1

  subroutine convert_r4_u1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    real(real32), allocatable :: src(:)
    integer(c_signed_char), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = to_uchar_r4(src(i))
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = to_uchar_r4(src(i))
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_r4_u1

  subroutine convert_i4_i4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(int32), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = src(i)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i4_i4

  subroutine convert_i4_r4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(int32), allocatable :: src(:)
    real(real32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = real(src(i), real32)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = real(src(i), real32)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i4_r4

  subroutine convert_i4_i1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(int32), allocatable :: src(:)
    integer(c_signed_char), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = int(src(i), c_signed_char)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = int(src(i), c_signed_char)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i4_i1

  subroutine convert_i4_u1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(int32), allocatable :: src(:)
    integer(c_signed_char), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = to_uchar_i4(src(i))
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = to_uchar_i4(src(i))
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i4_u1

  subroutine convert_i1_i4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:)
    integer(int32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = int(src(i), int32)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = int(src(i), int32)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i1_i4

  subroutine convert_i1_r4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:)
    real(real32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = real(src(i), real32)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = real(src(i), real32)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i1_r4

  subroutine convert_u1_i4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:)
    integer(int32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = from_uchar(src(i))
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = from_uchar(src(i))
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_u1_i4

  subroutine convert_u1_r4(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:)
    real(real32), allocatable :: dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = real(from_uchar(src(i)), real32)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = real(from_uchar(src(i)), real32)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_u1_r4

  subroutine convert_i1_u1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = to_uchar_i4(int(src(i), int32))
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = to_uchar_i4(int(src(i), int32))
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i1_u1

  subroutine convert_u1_i1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = int(from_uchar(src(i)), c_signed_char)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = int(from_uchar(src(i)), c_signed_char)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_u1_i1

  subroutine convert_u1_u1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = src(i)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_u1_u1

  subroutine convert_i1_i1(nelems, niters, src_bytes, dst_bytes)
    integer, intent(in) :: nelems, niters, src_bytes, dst_bytes
    integer(c_signed_char), allocatable :: src(:), dst(:)
    integer :: iter, i, ls, gs
    real(real64) :: start_time, end_time

    allocate(src(nelems), dst(nelems))
    ls = min(nelems, 256)
    gs = (nelems + ls - 1) / ls
    !$omp target data map(alloc: src(1:nelems), dst(1:nelems))
    !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
    do i = 1, nelems
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do
    start_time = omp_get_wtime()
    do iter = 1, niters
      !$omp target teams distribute parallel do num_teams(gs) num_threads(ls)
      do i = 1, nelems
        dst(i) = src(i)
      end do
      !$omp end target teams distribute parallel do
    end do
    end_time = omp_get_wtime()
    !$omp end target data
    call print_timing(src_bytes, dst_bytes, nelems, niters, start_time, end_time)
    deallocate(src, dst)
  end subroutine convert_i1_i1

end program main
