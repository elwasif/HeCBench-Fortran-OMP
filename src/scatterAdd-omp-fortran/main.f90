! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: max_threads_per_block = 512

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand") result(value)
      import :: c_int
      integer(c_int) :: value
    end function c_rand
  end interface

  character(len=256) :: arg0, arg
  integer :: batch_size, output_size, vector_dim, repeat

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 4) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <batch size> <output size> <vector dimension> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg); read(arg, *) batch_size
  call get_command_argument(2, arg); read(arg, *) output_size
  call get_command_argument(3, arg); read(arg, *) vector_dim
  call get_command_argument(4, arg); read(arg, *) repeat
  if (batch_size <= 0 .or. output_size <= 0 .or. vector_dim <= 0 .or. repeat <= 0) stop 1

  write(*,'(A,I0)') 'batch_size: ', batch_size
  write(*,'(A,I0)') 'output_size (range of index values): ', output_size
  write(*,'(A,I0)') 'vector_dimension: ', vector_dim

  call index_accumulate(batch_size, output_size, vector_dim, repeat)

contains

  subroutine index_accumulate(batch_size, output_size, vector_dim, repeat)
    integer, intent(in) :: batch_size, output_size, vector_dim, repeat
    integer :: i, iter
    integer, allocatable :: index(:)
    real(real32), allocatable :: source(:), output(:), output_ref(:)
    real(real64) :: start_time, elapsed
    logical :: ok

    allocate(index(batch_size))
    allocate(source(batch_size * vector_dim))
    allocate(output(output_size * vector_dim))
    allocate(output_ref(output_size * vector_dim))

    call c_srand(2_c_int)
    do i = 1, batch_size
      index(i) = modulo(int(c_rand()), output_size)
    end do
    source = -1.0_real32
    output = 0.0_real32
    output_ref = 0.0_real32
    call scatter_add_reference(batch_size, vector_dim, output_ref, index, source)

    !$omp target data map(to: source(1:batch_size*vector_dim), index(1:batch_size)) &
    !$omp& map(tofrom: output(1:output_size*vector_dim))
    do iter = 1, 10
      output = 0.0_real32
      !$omp target update to(output(1:output_size*vector_dim))
      call scatter_add2_kernel(index, source, output, batch_size, vector_dim)
    end do

    !$omp target update from(output(1:output_size*vector_dim))
    ok = .true.
    do i = 1, output_size * vector_dim
      if (abs(output(i) - output_ref(i)) > 1.0e-3_real32) then
        write(*,'(A,I0,A,F0.6,1X,F0.6)') 'output ', i - 1, ': ', output(i), output_ref(i)
        ok = .false.
        exit
      end if
    end do
    if (ok) then
      write(*,'(A)') 'PASS'
    else
      write(*,'(A)') 'FAIL'
    end if

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call scatter_add_kernel(index, source, output, batch_size, vector_dim)
    end do
    elapsed = omp_get_wtime() - start_time
    write(*,'(A,F0.6,A)') 'Average execution time of kernel1: ', (elapsed * 1.0e6_real64) / real(repeat, real64), ' (us)'

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call scatter_add2_kernel(index, source, output, batch_size, vector_dim)
    end do
    elapsed = omp_get_wtime() - start_time
    !$omp end target data
    write(*,'(A,F0.6,A)') 'Average execution time of kernel2: ', (elapsed * 1.0e6_real64) / real(repeat, real64), ' (us)'

    deallocate(index, source, output, output_ref)
  end subroutine index_accumulate

  subroutine scatter_add_reference(batch_size, vector_dim, output, index, source)
    integer, intent(in) :: batch_size, vector_dim
    integer, intent(in) :: index(batch_size)
    real(real32), intent(in) :: source(batch_size * vector_dim)
    real(real32), intent(inout) :: output(:)
    integer :: d, i

    do d = 0, vector_dim - 1
      do i = 0, batch_size - 1
        output(index(i + 1) * vector_dim + d + 1) = &
          output(index(i + 1) * vector_dim + d + 1) + source(i * vector_dim + d + 1)
      end do
    end do
  end subroutine scatter_add_reference

  subroutine scatter_add_kernel(index, source, output, batch_size, vector_dim)
    integer, intent(in) :: batch_size, vector_dim
    integer, intent(in) :: index(batch_size)
    real(real32), intent(in) :: source(batch_size * vector_dim)
    real(real32), intent(inout) :: output(:)
    integer :: d, i, out_pos, src_pos

    !$omp target teams distribute parallel do collapse(2) private(out_pos, src_pos) &
    !$omp& thread_limit(max_threads_per_block)
    do d = 0, vector_dim - 1
      do i = 0, batch_size - 1
        out_pos = index(i + 1) * vector_dim + d + 1
        src_pos = i * vector_dim + d + 1
        !$omp atomic update
        output(out_pos) = output(out_pos) + source(src_pos)
        !$omp end atomic
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine scatter_add_kernel

  subroutine scatter_add2_kernel(index, source, output, batch_size, vector_dim)
    integer, intent(in) :: batch_size, vector_dim
    integer, intent(in) :: index(batch_size)
    real(real32), intent(in) :: source(batch_size * vector_dim)
    real(real32), intent(inout) :: output(:)
    integer :: i, d, out_pos, src_pos

    !$omp target teams distribute parallel do collapse(2) private(out_pos, src_pos) &
    !$omp& thread_limit(max_threads_per_block)
    do i = 0, batch_size - 1
      do d = 0, vector_dim - 1
        out_pos = index(i + 1) * vector_dim + d + 1
        src_pos = i * vector_dim + d + 1
        !$omp atomic update
        output(out_pos) = output(out_pos) + source(src_pos)
        !$omp end atomic
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine scatter_add2_kernel

end program main
