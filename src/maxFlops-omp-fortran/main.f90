! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_double, c_long
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer, parameter :: num_floats = 2 * 1024 * 1024
  integer :: repeat

  interface
    subroutine c_srand48(seed) bind(C, name="srand48")
      import :: c_long
      integer(c_long), value :: seed
    end subroutine c_srand48

    real(c_double) function c_drand48() bind(C, name="drand48")
      import :: c_double
    end function c_drand48
  end interface

  if (command_argument_count() /= 1) then
    call print_usage()
    stop 1
  end if

  repeat = read_int_arg(1)
  if (repeat < 0) stop 1

  write(*,'(A)') '=== Single-precision floating-point kernels ==='
  call test_sp(repeat, num_floats)
  write(*,'(A)') '=== Double-precision floating-point kernels ==='
  call test_dp(repeat, num_floats)

contains

  subroutine print_usage()
    character(len=256) :: arg0
    call get_command_argument(0, arg0)
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
  end subroutine print_usage

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  subroutine test_sp(repeat, n)
    integer, intent(in) :: repeat, n
    real(real32), allocatable :: data(:)
    integer :: j
    allocate(data(n))
    call c_srand48(123_c_long)
    do j = 1, n / 2
      data(j) = real(c_drand48() * 10.0_c_double, real32)
      data(n - j + 1) = data(j)
    end do
    !$omp target data map(alloc: data(1:n))
    do j = 1, 4
      call add1_sp(data, n, repeat)
      call add2_sp(data, n, repeat)
      call add4_sp(data, n, repeat)
      call add8_sp(data, n, repeat)
    end do
    call timed_add1_sp(data, n, repeat, 'Add1')
    call timed_add2_sp(data, n, repeat, 'Add2')
    call timed_add4_sp(data, n, repeat, 'Add4')
    call timed_add8_sp(data, n, repeat, 'Add8')
    do j = 1, 4
      call mul1_sp(data, n, repeat)
      call mul2_sp(data, n, repeat)
      call mul4_sp(data, n, repeat)
      call mul8_sp(data, n, repeat)
    end do
    call timed_mul1_sp(data, n, repeat, 'Mul1')
    call timed_mul2_sp(data, n, repeat, 'Mul2')
    call timed_mul4_sp(data, n, repeat, 'Mul4')
    call timed_mul8_sp(data, n, repeat, 'Mul8')
    do j = 1, 4
      call madd1_sp(data, n, repeat)
      call madd2_sp(data, n, repeat)
      call madd4_sp(data, n, repeat)
      call madd8_sp(data, n, repeat)
    end do
    call timed_madd1_sp(data, n, repeat, 'MAdd1')
    call timed_madd2_sp(data, n, repeat, 'MAdd2')
    call timed_madd4_sp(data, n, repeat, 'MAdd4')
    call timed_madd8_sp(data, n, repeat, 'MAdd8')
    do j = 1, 4
      call mulmadd1_sp(data, n, repeat)
      call mulmadd2_sp(data, n, repeat)
      call mulmadd4_sp(data, n, repeat)
      call mulmadd8_sp(data, n, repeat)
    end do
    call timed_mulmadd1_sp(data, n, repeat, 'MulMAdd1')
    call timed_mulmadd2_sp(data, n, repeat, 'MulMAdd2')
    call timed_mulmadd4_sp(data, n, repeat, 'MulMAdd4')
    call timed_mulmadd8_sp(data, n, repeat, 'MulMAdd8')
    !$omp end target data
    deallocate(data)
  end subroutine test_sp

  subroutine timed_add1_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add1_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add1_sp

  subroutine timed_add2_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add2_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add2_sp

  subroutine timed_add4_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add4_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add4_sp

  subroutine timed_add8_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add8_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add8_sp

  subroutine timed_mul1_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul1_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul1_sp

  subroutine timed_mul2_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul2_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul2_sp

  subroutine timed_mul4_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul4_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul4_sp

  subroutine timed_mul8_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul8_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul8_sp

  subroutine timed_madd1_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd1_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd1_sp

  subroutine timed_madd2_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd2_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd2_sp

  subroutine timed_madd4_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd4_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd4_sp

  subroutine timed_madd8_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd8_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd8_sp

  subroutine timed_mulmadd1_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd1_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd1_sp

  subroutine timed_mulmadd2_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd2_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd2_sp

  subroutine timed_mulmadd4_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd4_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd4_sp

  subroutine timed_mulmadd8_sp(data, n, repeat, label)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd8_sp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd8_sp

  subroutine add1_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 240
          s = 10.0_real32 - s
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine add1_sp

  subroutine add2_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      do j = 1, repeat
        do k = 1, 120
          s = 10.0_real32 - s
          s2 = 10.0_real32 - s2
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine add2_sp

  subroutine add4_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      do j = 1, repeat
        do k = 1, 60
          s = 10.0_real32 - s
          s2 = 10.0_real32 - s2
          s3 = 10.0_real32 - s3
          s4 = 10.0_real32 - s4
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine add4_sp

  subroutine add8_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      s5 = 8.0_real32 - s
      s6 = 8.0_real32 - s2
      s7 = 7.0_real32 - s
      s8 = 7.0_real32 - s2
      do j = 1, repeat
        do k = 1, 30
          s = 10.0_real32 - s
          s2 = 10.0_real32 - s2
          s3 = 10.0_real32 - s3
          s4 = 10.0_real32 - s4
          s5 = 10.0_real32 - s5
          s6 = 10.0_real32 - s6
          s7 = 10.0_real32 - s7
          s8 = 10.0_real32 - s8
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine add8_sp

  subroutine mul1_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real32
      do j = 1, repeat
        do k = 1, 200
          s = s * s * 1.01_real32
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul1_sp

  subroutine mul2_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real32
      s2 = s - 0.0001_real32
      do j = 1, repeat
        do k = 1, 100
          s = s * s * 1.01_real32
          s2 = s2 * s2 * 1.01_real32
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul2_sp

  subroutine mul4_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real32
      s2 = s - 0.0001_real32
      s3 = s - 0.0002_real32
      s4 = s - 0.0003_real32
      do j = 1, repeat
        do k = 1, 50
          s = s * s * 1.01_real32
          s2 = s2 * s2 * 1.01_real32
          s3 = s3 * s3 * 1.01_real32
          s4 = s4 * s4 * 1.01_real32
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul4_sp

  subroutine mul8_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real32
      s2 = s - 0.0001_real32
      s3 = s - 0.0002_real32
      s4 = s - 0.0003_real32
      s5 = s - 0.0004_real32
      s6 = s - 0.0005_real32
      s7 = s - 0.0006_real32
      s8 = s - 0.0007_real32
      do j = 1, repeat
        do k = 1, 25
          s = s * s * 1.01_real32
          s2 = s2 * s2 * 1.01_real32
          s3 = s3 * s3 * 1.01_real32
          s4 = s4 * s4 * 1.01_real32
          s5 = s5 * s5 * 1.01_real32
          s6 = s6 * s6 * 1.01_real32
          s7 = s7 * s7 * 1.01_real32
          s8 = s8 * s8 * 1.01_real32
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul8_sp

  subroutine madd1_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 240
          s = 10.0_real32 - s * 0.9899_real32
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd1_sp

  subroutine madd2_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      do j = 1, repeat
        do k = 1, 120
          s = 10.0_real32 - s * 0.9899_real32
          s2 = 10.0_real32 - s2 * 0.9899_real32
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd2_sp

  subroutine madd4_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      do j = 1, repeat
        do k = 1, 60
          s = 10.0_real32 - s * 0.9899_real32
          s2 = 10.0_real32 - s2 * 0.9899_real32
          s3 = 10.0_real32 - s3 * 0.9899_real32
          s4 = 10.0_real32 - s4 * 0.9899_real32
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd4_sp

  subroutine madd8_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      s5 = 8.0_real32 - s
      s6 = 8.0_real32 - s2
      s7 = 7.0_real32 - s
      s8 = 7.0_real32 - s2
      do j = 1, repeat
        do k = 1, 30
          s = 10.0_real32 - s * 0.9899_real32
          s2 = 10.0_real32 - s2 * 0.9899_real32
          s3 = 10.0_real32 - s3 * 0.9899_real32
          s4 = 10.0_real32 - s4 * 0.9899_real32
          s5 = 10.0_real32 - s5 * 0.9899_real32
          s6 = 10.0_real32 - s6 * 0.9899_real32
          s7 = 10.0_real32 - s7 * 0.9899_real32
          s8 = 10.0_real32 - s8 * 0.9899_real32
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd8_sp

  subroutine mulmadd1_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 160
          s = (3.75_real32 - 0.355_real32 * s) * s
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd1_sp

  subroutine mulmadd2_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      do j = 1, repeat
        do k = 1, 80
          s = (3.75_real32 - 0.355_real32 * s) * s
          s2 = (3.75_real32 - 0.355_real32 * s2) * s2
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd2_sp

  subroutine mulmadd4_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      do j = 1, repeat
        do k = 1, 40
          s = (3.75_real32 - 0.355_real32 * s) * s
          s2 = (3.75_real32 - 0.355_real32 * s2) * s2
          s3 = (3.75_real32 - 0.355_real32 * s3) * s3
          s4 = (3.75_real32 - 0.355_real32 * s4) * s4
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd4_sp

  subroutine mulmadd8_sp(data, n, repeat)
    real(real32), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real32) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real32 - s
      s3 = 9.0_real32 - s
      s4 = 9.0_real32 - s2
      s5 = 8.0_real32 - s
      s6 = 8.0_real32 - s2
      s7 = 7.0_real32 - s
      s8 = 7.0_real32 - s2
      do j = 1, repeat
        do k = 1, 20
          s = (3.75_real32 - 0.355_real32 * s) * s
          s2 = (3.75_real32 - 0.355_real32 * s2) * s2
          s3 = (3.75_real32 - 0.355_real32 * s3) * s3
          s4 = (3.75_real32 - 0.355_real32 * s4) * s4
          s5 = (3.75_real32 - 0.355_real32 * s5) * s5
          s6 = (3.75_real32 - 0.355_real32 * s6) * s6
          s7 = (3.75_real32 - 0.355_real32 * s7) * s7
          s8 = (3.75_real32 - 0.355_real32 * s8) * s8
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd8_sp

  subroutine test_dp(repeat, n)
    integer, intent(in) :: repeat, n
    real(real64), allocatable :: data(:)
    integer :: j
    allocate(data(n))
    call c_srand48(123_c_long)
    do j = 1, n / 2
      data(j) = real(c_drand48() * 10.0_c_double, real64)
      data(n - j + 1) = data(j)
    end do
    !$omp target data map(alloc: data(1:n))
    do j = 1, 4
      call add1_dp(data, n, repeat)
      call add2_dp(data, n, repeat)
      call add4_dp(data, n, repeat)
      call add8_dp(data, n, repeat)
    end do
    call timed_add1_dp(data, n, repeat, 'Add1')
    call timed_add2_dp(data, n, repeat, 'Add2')
    call timed_add4_dp(data, n, repeat, 'Add4')
    call timed_add8_dp(data, n, repeat, 'Add8')
    do j = 1, 4
      call mul1_dp(data, n, repeat)
      call mul2_dp(data, n, repeat)
      call mul4_dp(data, n, repeat)
      call mul8_dp(data, n, repeat)
    end do
    call timed_mul1_dp(data, n, repeat, 'Mul1')
    call timed_mul2_dp(data, n, repeat, 'Mul2')
    call timed_mul4_dp(data, n, repeat, 'Mul4')
    call timed_mul8_dp(data, n, repeat, 'Mul8')
    do j = 1, 4
      call madd1_dp(data, n, repeat)
      call madd2_dp(data, n, repeat)
      call madd4_dp(data, n, repeat)
      call madd8_dp(data, n, repeat)
    end do
    call timed_madd1_dp(data, n, repeat, 'MAdd1')
    call timed_madd2_dp(data, n, repeat, 'MAdd2')
    call timed_madd4_dp(data, n, repeat, 'MAdd4')
    call timed_madd8_dp(data, n, repeat, 'MAdd8')
    do j = 1, 4
      call mulmadd1_dp(data, n, repeat)
      call mulmadd2_dp(data, n, repeat)
      call mulmadd4_dp(data, n, repeat)
      call mulmadd8_dp(data, n, repeat)
    end do
    call timed_mulmadd1_dp(data, n, repeat, 'MulMAdd1')
    call timed_mulmadd2_dp(data, n, repeat, 'MulMAdd2')
    call timed_mulmadd4_dp(data, n, repeat, 'MulMAdd4')
    call timed_mulmadd8_dp(data, n, repeat, 'MulMAdd8')
    !$omp end target data
    deallocate(data)
  end subroutine test_dp

  subroutine timed_add1_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add1_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add1_dp

  subroutine timed_add2_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add2_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add2_dp

  subroutine timed_add4_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add4_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add4_dp

  subroutine timed_add8_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call add8_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_add8_dp

  subroutine timed_mul1_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul1_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul1_dp

  subroutine timed_mul2_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul2_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul2_dp

  subroutine timed_mul4_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul4_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul4_dp

  subroutine timed_mul8_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mul8_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mul8_dp

  subroutine timed_madd1_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd1_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd1_dp

  subroutine timed_madd2_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd2_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd2_dp

  subroutine timed_madd4_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd4_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd4_dp

  subroutine timed_madd8_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call madd8_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_madd8_dp

  subroutine timed_mulmadd1_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd1_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd1_dp

  subroutine timed_mulmadd2_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd2_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd2_dp

  subroutine timed_mulmadd4_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd4_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd4_dp

  subroutine timed_mulmadd8_dp(data, n, repeat, label)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    character(len=*), intent(in) :: label
    real(real64) :: t0, elapsed
    !$omp target update to(data(1:n))
    t0 = omp_get_wtime()
    call mulmadd8_dp(data, n, repeat)
    elapsed = omp_get_wtime() - t0
    write(*,'(A,A,A,F0.6,A)') 'kernel execution time (', trim(label), '): ', elapsed, ' (s)'
  end subroutine timed_mulmadd8_dp

  subroutine add1_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 240
          s = 10.0_real64 - s
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine add1_dp

  subroutine add2_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      do j = 1, repeat
        do k = 1, 120
          s = 10.0_real64 - s
          s2 = 10.0_real64 - s2
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine add2_dp

  subroutine add4_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      do j = 1, repeat
        do k = 1, 60
          s = 10.0_real64 - s
          s2 = 10.0_real64 - s2
          s3 = 10.0_real64 - s3
          s4 = 10.0_real64 - s4
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine add4_dp

  subroutine add8_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      s5 = 8.0_real64 - s
      s6 = 8.0_real64 - s2
      s7 = 7.0_real64 - s
      s8 = 7.0_real64 - s2
      do j = 1, repeat
        do k = 1, 30
          s = 10.0_real64 - s
          s2 = 10.0_real64 - s2
          s3 = 10.0_real64 - s3
          s4 = 10.0_real64 - s4
          s5 = 10.0_real64 - s5
          s6 = 10.0_real64 - s6
          s7 = 10.0_real64 - s7
          s8 = 10.0_real64 - s8
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine add8_dp

  subroutine mul1_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real64
      do j = 1, repeat
        do k = 1, 200
          s = s * s * 1.01_real64
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul1_dp

  subroutine mul2_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real64
      s2 = s - 0.0001_real64
      do j = 1, repeat
        do k = 1, 100
          s = s * s * 1.01_real64
          s2 = s2 * s2 * 1.01_real64
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul2_dp

  subroutine mul4_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real64
      s2 = s - 0.0001_real64
      s3 = s - 0.0002_real64
      s4 = s - 0.0003_real64
      do j = 1, repeat
        do k = 1, 50
          s = s * s * 1.01_real64
          s2 = s2 * s2 * 1.01_real64
          s3 = s3 * s3 * 1.01_real64
          s4 = s4 * s4 * 1.01_real64
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul4_dp

  subroutine mul8_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid) - data(gid) + 0.999_real64
      s2 = s - 0.0001_real64
      s3 = s - 0.0002_real64
      s4 = s - 0.0003_real64
      s5 = s - 0.0004_real64
      s6 = s - 0.0005_real64
      s7 = s - 0.0006_real64
      s8 = s - 0.0007_real64
      do j = 1, repeat
        do k = 1, 25
          s = s * s * 1.01_real64
          s2 = s2 * s2 * 1.01_real64
          s3 = s3 * s3 * 1.01_real64
          s4 = s4 * s4 * 1.01_real64
          s5 = s5 * s5 * 1.01_real64
          s6 = s6 * s6 * 1.01_real64
          s7 = s7 * s7 * 1.01_real64
          s8 = s8 * s8 * 1.01_real64
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine mul8_dp

  subroutine madd1_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 240
          s = 10.0_real64 - s * 0.9899_real64
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd1_dp

  subroutine madd2_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      do j = 1, repeat
        do k = 1, 120
          s = 10.0_real64 - s * 0.9899_real64
          s2 = 10.0_real64 - s2 * 0.9899_real64
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd2_dp

  subroutine madd4_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      do j = 1, repeat
        do k = 1, 60
          s = 10.0_real64 - s * 0.9899_real64
          s2 = 10.0_real64 - s2 * 0.9899_real64
          s3 = 10.0_real64 - s3 * 0.9899_real64
          s4 = 10.0_real64 - s4 * 0.9899_real64
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd4_dp

  subroutine madd8_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      s5 = 8.0_real64 - s
      s6 = 8.0_real64 - s2
      s7 = 7.0_real64 - s
      s8 = 7.0_real64 - s2
      do j = 1, repeat
        do k = 1, 30
          s = 10.0_real64 - s * 0.9899_real64
          s2 = 10.0_real64 - s2 * 0.9899_real64
          s3 = 10.0_real64 - s3 * 0.9899_real64
          s4 = 10.0_real64 - s4 * 0.9899_real64
          s5 = 10.0_real64 - s5 * 0.9899_real64
          s6 = 10.0_real64 - s6 * 0.9899_real64
          s7 = 10.0_real64 - s7 * 0.9899_real64
          s8 = 10.0_real64 - s8 * 0.9899_real64
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine madd8_dp

  subroutine mulmadd1_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s)
    do gid = 1, n
      s = data(gid)
      do j = 1, repeat
        do k = 1, 160
          s = (3.75_real64 - 0.355_real64 * s) * s
        end do
      end do
      data(gid) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd1_dp

  subroutine mulmadd2_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      do j = 1, repeat
        do k = 1, 80
          s = (3.75_real64 - 0.355_real64 * s) * s
          s2 = (3.75_real64 - 0.355_real64 * s2) * s2
        end do
      end do
      data(gid) = s + s2
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd2_dp

  subroutine mulmadd4_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      do j = 1, repeat
        do k = 1, 40
          s = (3.75_real64 - 0.355_real64 * s) * s
          s2 = (3.75_real64 - 0.355_real64 * s2) * s2
          s3 = (3.75_real64 - 0.355_real64 * s3) * s3
          s4 = (3.75_real64 - 0.355_real64 * s4) * s4
        end do
      end do
      data(gid) = (s + s2) + (s3 + s4)
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd4_dp

  subroutine mulmadd8_dp(data, n, repeat)
    real(real64), intent(inout) :: data(:)
    integer, intent(in) :: n, repeat
    integer :: gid, j, k
    real(real64) :: s, s2, s3, s4, s5, s6, s7, s8
    !$omp target teams distribute parallel do thread_limit(block_size) private(j,k,s,s2,s3,s4,s5,s6,s7,s8)
    do gid = 1, n
      s = data(gid)
      s2 = 10.0_real64 - s
      s3 = 9.0_real64 - s
      s4 = 9.0_real64 - s2
      s5 = 8.0_real64 - s
      s6 = 8.0_real64 - s2
      s7 = 7.0_real64 - s
      s8 = 7.0_real64 - s2
      do j = 1, repeat
        do k = 1, 20
          s = (3.75_real64 - 0.355_real64 * s) * s
          s2 = (3.75_real64 - 0.355_real64 * s2) * s2
          s3 = (3.75_real64 - 0.355_real64 * s3) * s3
          s4 = (3.75_real64 - 0.355_real64 * s4) * s4
          s5 = (3.75_real64 - 0.355_real64 * s5) * s5
          s6 = (3.75_real64 - 0.355_real64 * s6) * s6
          s7 = (3.75_real64 - 0.355_real64 * s7) * s7
          s8 = (3.75_real64 - 0.355_real64 * s8) * s8
        end do
      end do
      data(gid) = ((s + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8))
    end do
    !$omp end target teams distribute parallel do
  end subroutine mulmadd8_dp

end program main
