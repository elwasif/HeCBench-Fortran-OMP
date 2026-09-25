! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_long_long
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand_mod(modulus) bind(C, name="concat_rand_mod") result(value)
      import :: c_int
      integer(c_int), value :: modulus
      integer(c_int) :: value
    end function c_rand_mod

    subroutine c_fill_rand_float(values, n) bind(C, name="concat_fill_rand_float")
      import :: c_float, c_long_long
      real(c_float), intent(out) :: values(*)
      integer(c_long_long), value :: n
    end subroutine c_fill_rand_float
  end interface

  integer, parameter :: seq_len = 1024
  integer, parameter :: batch_size = 8
  integer, parameter :: beam_size = 8
  integer, parameter :: head_dim = 128
  integer :: repeat
  character(len=256) :: arg0, arg

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) repeat
  if (repeat <= 0) stop 1

  call run_case(6, repeat)
  call run_case(12, repeat)
  call run_case(24, repeat)
  call run_case(48, repeat)

contains

  subroutine run_case(nhead, repeat)
    integer, intent(in) :: nhead, repeat
    integer :: hidden_dim, sz0, sz2, sl1, sl2, i
    integer(int64) :: inp1_size, inp2_size, outp_size
    real(real32), allocatable :: inp1(:), inp2(:), outp(:), outp_ref(:)
    real(real64) :: size_bytes_gb, start_time, end_time, avg_time_us
    logical :: ok

    hidden_dim = nhead * 128
    sz0 = batch_size * beam_size * nhead
    sz2 = hidden_dim / nhead
    call c_srand(int(nhead, c_int))
    sl1 = int(c_rand_mod(int(seq_len - 1, c_int))) + 1
    sl2 = seq_len - sl1

    write(*,*)
    write(*,'(A,I0,A)', advance='no') 'num_head = ', nhead, char(9)
    write(*,'(A,I0,A)', advance='no') 'seq_len = ', seq_len, char(9)
    write(*,'(A,I0,A)', advance='no') 'batch_size = ', batch_size, char(9)
    write(*,'(A,I0,A)', advance='no') 'hidden_dimension = ', hidden_dim, char(9)
    write(*,'(A,I0)') 'beam_size = ', beam_size

    inp1_size = int(batch_size, int64) * beam_size * hidden_dim * sl1
    inp2_size = int(batch_size, int64) * beam_size * hidden_dim * sl2
    outp_size = int(batch_size, int64) * beam_size * hidden_dim * seq_len
    size_bytes_gb = real(2_int64 * outp_size * 4_int64, real64) * 1.0e-9_real64
    write(*,'(A,F4.2)') 'Total device memory usage (GB) = ', size_bytes_gb

    allocate(inp1(0:inp1_size - 1), inp2(0:inp2_size - 1), outp(0:outp_size - 1), outp_ref(0:outp_size - 1))
    call c_fill_rand_float(inp1, int(inp1_size, c_long_long))
    call c_fill_rand_float(inp2, int(inp2_size, c_long_long))
    outp = -1.0_real32
    outp_ref = -2.0_real32

    !$omp target data map(to: inp1(0:inp1_size - 1), inp2(0:inp2_size - 1)) map(alloc: outp(0:outp_size - 1))
    call concat_kernel(inp1, inp2, outp, sz0, sz2, sl1, sl2)
    !$omp target update from(outp(0:outp_size - 1))

    call concat_cpu(inp1, inp2, outp_ref, sz0, sz2, sl1, sl2)
    ok = .true.
    do i = 0, int(outp_size - 1_int64)
      if (outp(i) /= outp_ref(i)) then
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
    do i = 1, repeat
      call concat_kernel(inp1, inp2, outp, sz0, sz2, sl1, sl2)
    end do
    end_time = omp_get_wtime()
    avg_time_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
    write(*,'(A,F0.6,A)') 'Average kernel execution time: ', avg_time_us, ' (us)'
    write(*,'(A,F0.6,A)') 'Average kernel throughput : ', size_bytes_gb / (avg_time_us * 1.0e-6_real64), ' (GB/s)'
    !$omp end target data

    deallocate(inp1, inp2, outp, outp_ref)
  end subroutine run_case

  subroutine concat_kernel(inp1, inp2, output, sz0, sz2, sz1_1, sz1_2)
    real(real32), intent(in) :: inp1(0:), inp2(0:)
    real(real32), intent(inout) :: output(0:)
    integer, intent(in) :: sz0, sz2, sz1_1, sz1_2
    integer(int64) :: idx, nele, idx_y, idx0, src_index, dst_index
    integer :: idx_x, idx1, sz1

    nele = int(sz0, int64) * sz2 * (sz1_1 + sz1_2)
    !$omp target teams distribute parallel do thread_limit(256) private(idx, idx_x, idx_y, idx1, idx0, sz1, src_index, dst_index)
    do idx = 0_int64, nele - 1_int64
      dst_index = idx
      idx_x = int(mod(idx, int(sz2, int64)))
      idx_y = idx / sz2
      idx1 = int(mod(idx_y, int(sz1_1 + sz1_2, int64)))
      idx0 = idx_y / (sz1_1 + sz1_2)
      if (idx1 < sz1_1) then
        sz1 = sz1_1
        src_index = flat_3dim(idx0, int(idx1, int64), int(idx_x, int64), sz1, sz2)
        output(dst_index) = inp1(src_index)
      else
        idx1 = idx1 - sz1_1
        sz1 = sz1_2
        src_index = flat_3dim(idx0, int(idx1, int64), int(idx_x, int64), sz1, sz2)
        output(dst_index) = inp2(src_index)
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine concat_kernel

  subroutine concat_cpu(inp1, inp2, output, sz0, sz2, sz1_1, sz1_2)
    real(real32), intent(in) :: inp1(0:), inp2(0:)
    real(real32), intent(out) :: output(0:)
    integer, intent(in) :: sz0, sz2, sz1_1, sz1_2
    integer :: idx1, idx2, sz1
    integer(int64) :: idx0, src_index, dst_index

    do idx0 = 0_int64, int(sz0 - 1, int64)
      do idx1 = 0, sz1_1 + sz1_2 - 1
        do idx2 = 0, sz2 - 1
          dst_index = flat_3dim(idx0, int(idx1, int64), int(idx2, int64), sz1_1 + sz1_2, sz2)
          if (idx1 < sz1_1) then
            sz1 = sz1_1
            src_index = flat_3dim(idx0, int(idx1, int64), int(idx2, int64), sz1, sz2)
            output(dst_index) = inp1(src_index)
          else
            sz1 = sz1_2
            src_index = flat_3dim(idx0, int(idx1 - sz1_1, int64), int(idx2, int64), sz1, sz2)
            output(dst_index) = inp2(src_index)
          end if
        end do
      end do
    end do
  end subroutine concat_cpu

  pure integer(int64) function flat_3dim(id1, id2, id3, dim2, dim3)
    integer(int64), intent(in) :: id1, id2, id3
    integer, intent(in) :: dim2, dim3

    flat_3dim = id1 * dim2 * dim3 + id2 * dim3 + id3
  end function flat_3dim

end program main
