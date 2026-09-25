! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  integer, parameter :: chunk_size = 512
  integer, parameter :: value_count = 16
  integer, parameter :: warp_size = 32

  type :: uint4
    integer(int32) :: x, y, z, w
  end type uint4

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
  integer :: nkeys, repeat_count
  integer(int32), allocatable :: keys(:), out(:)
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of keys> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) nkeys
  call get_command_argument(2, arg)
  read(arg, *) repeat_count
  if (nkeys <= 0 .or. repeat_count <= 0 .or. modulo(nkeys, chunk_size) /= 0) stop 1

  allocate(keys(nkeys), out(nkeys))
  call c_srand(512_c_int)
  call initialize(keys)
  out = keys

  call split_sort(out, repeat_count)
  ok = verify(out, keys)
  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
    stop 1
  end if

  deallocate(keys, out)

contains

  subroutine initialize(keys)
    integer(int32), intent(out) :: keys(:)
    integer :: i

    do i = 1, size(keys)
      ! Preserve the C++ original's srand(512)/rand() input stream.
      keys(i) = int(modulo(c_rand(), value_count), int32)
    end do
  end subroutine initialize

  subroutine split_sort(out, repeat_count)
    integer(int32), intent(inout) :: out(:)
    integer, intent(in) :: repeat_count
    integer :: iter, nkeys
    real(real64) :: start_time, elapsed_us

    nkeys = size(out)
    !$omp target data map(tofrom: out(1:nkeys))
    start_time = omp_get_wtime()
    do iter = 1, repeat_count
      call sort_chunks(out, nkeys)
    end do
    elapsed_us = (omp_get_wtime() - start_time) * 1.0e6_real64 / real(repeat_count, real64)
    write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed_us, ' (us)'
    !$omp end target data
  end subroutine split_sort

  function scanwarp(val, sData, maxlevel) result(scan)
    integer(int32), value :: val
    integer(int32), volatile, intent(inout) :: sData(:)
    integer, value :: maxlevel
    integer(int32) :: scan
    integer :: localId, idx

    localId = omp_get_thread_num()
    idx = 2 * localId - iand(localId, warp_size - 1)
    sData(idx + 1) = 0_int32
    idx = idx + warp_size
    sData(idx + 1) = val

    if (0 <= maxlevel) sData(idx + 1) = sData(idx + 1) + sData(idx)
    if (1 <= maxlevel) sData(idx + 1) = sData(idx + 1) + sData(idx - 1)
    if (2 <= maxlevel) sData(idx + 1) = sData(idx + 1) + sData(idx - 3)
    if (3 <= maxlevel) sData(idx + 1) = sData(idx + 1) + sData(idx - 7)
    if (4 <= maxlevel) sData(idx + 1) = sData(idx + 1) + sData(idx - 15)

    scan = sData(idx + 1) - val
  end function scanwarp

  function scan4(idata, ptr) result(val4)
    type(uint4), value :: idata
    integer(int32), volatile, intent(inout) :: ptr(:)
    type(uint4) :: val4
    integer(int32) :: sum0, sum1, sum2, val
    integer :: idx

    idx = omp_get_thread_num()
    val4 = idata
    sum0 = val4%x
    sum1 = val4%y + sum0
    sum2 = val4%z + sum1

    val = val4%w + sum2
    val = scanwarp(val, ptr, 4)
    !$omp barrier

    if (iand(idx, warp_size - 1) == warp_size - 1) then
      ptr(ishft(idx, -5) + 1) = val + val4%w + sum2
    end if
    !$omp barrier

    if (idx < warp_size) ptr(idx + 1) = scanwarp(ptr(idx + 1), ptr, 2)
    !$omp barrier

    val = val + ptr(ishft(idx, -5) + 1)

    val4%x = val
    val4%y = val + sum0
    val4%z = val + sum1
    val4%w = val + sum2
  end function scan4

  function rank4(preds, sMem, numtrue) result(rank)
    type(uint4), value :: preds
    integer(int32), volatile, intent(inout) :: sMem(:), numtrue(:)
    type(uint4) :: address, rank
    integer :: localId, localSize, idx

    localId = omp_get_thread_num()
    localSize = omp_get_num_threads()
    address = scan4(preds, sMem)

    if (localId == localSize - 1) numtrue(1) = address%w + preds%w
    !$omp barrier

    idx = localId * 4
    if (preds%x /= 0) then
      rank%x = address%x
    else
      rank%x = numtrue(1) + idx - address%x
    end if
    if (preds%y /= 0) then
      rank%y = address%y
    else
      rank%y = numtrue(1) + idx + 1 - address%y
    end if
    if (preds%z /= 0) then
      rank%z = address%z
    else
      rank%z = numtrue(1) + idx + 2 - address%z
    end if
    if (preds%w /= 0) then
      rank%w = address%w
    else
      rank%w = numtrue(1) + idx + 3 - address%w
    end if
  end function rank4

  subroutine sort_chunks(out, nkeys)
    integer(int32), intent(inout) :: out(:)
    integer, intent(in) :: nkeys
    integer, parameter :: startbit = 0, nbits = 4, threads = 128
    integer :: teams

    teams = nkeys / 4 / threads
    !$omp target teams num_teams(teams) thread_limit(threads)
    block
      integer(int32), volatile :: numtrue(1), sMem(chunk_size)

      !$omp parallel
      block
        integer :: localId, localSize, globalId, shift
        type(uint4) :: key, lsb, r

        localId = omp_get_thread_num()
        localSize = omp_get_num_threads()
        globalId = omp_get_team_num() * localSize + localId

        key%x = out(globalId * 4 + 1)
        key%y = out(globalId * 4 + 2)
        key%z = out(globalId * 4 + 3)
        key%w = out(globalId * 4 + 4)

        do shift = startbit, startbit + nbits - 1
          lsb%x = merge(1_int32, 0_int32, iand(ishft(key%x, -shift), 1_int32) == 0_int32)
          lsb%y = merge(1_int32, 0_int32, iand(ishft(key%y, -shift), 1_int32) == 0_int32)
          lsb%z = merge(1_int32, 0_int32, iand(ishft(key%z, -shift), 1_int32) == 0_int32)
          lsb%w = merge(1_int32, 0_int32, iand(ishft(key%w, -shift), 1_int32) == 0_int32)

          r = rank4(lsb, sMem, numtrue)

          sMem(iand(r%x, 3_int32) * localSize + ishft(r%x, -2) + 1) = key%x
          sMem(iand(r%y, 3_int32) * localSize + ishft(r%y, -2) + 1) = key%y
          sMem(iand(r%z, 3_int32) * localSize + ishft(r%z, -2) + 1) = key%z
          sMem(iand(r%w, 3_int32) * localSize + ishft(r%w, -2) + 1) = key%w
          !$omp barrier

          key%x = sMem(localId + 1)
          key%y = sMem(localId + localSize + 1)
          key%z = sMem(localId + 2 * localSize + 1)
          key%w = sMem(localId + 3 * localSize + 1)
          !$omp barrier
        end do

        out(globalId * 4 + 1) = key%x
        out(globalId * 4 + 2) = key%y
        out(globalId * 4 + 3) = key%z
        out(globalId * 4 + 4) = key%w
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine sort_chunks

  logical function verify(sorted_keys, keys)
    integer(int32), intent(in) :: sorted_keys(:), keys(:)
    integer :: nkeys, base, i, value
    integer :: hist_original(0:value_count - 1), hist_sorted(0:value_count - 1)

    nkeys = size(sorted_keys)
    verify = .true.
    do base = 1, nkeys, chunk_size
      do i = 0, chunk_size - 2
        if (sorted_keys(base + i) > sorted_keys(base + i + 1)) then
          verify = .false.
          return
        end if
      end do
    end do
    do i = 1, nkeys
      if (sorted_keys(i) < 0 .or. sorted_keys(i) >= value_count) then
        verify = .false.
        return
      end if
    end do
    do base = 1, nkeys, chunk_size
      hist_original = 0
      hist_sorted = 0
      do i = 0, chunk_size - 1
        value = keys(base + i)
        hist_original(value) = hist_original(value) + 1
        value = sorted_keys(base + i)
        hist_sorted(value) = hist_sorted(value) + 1
      end do
      do value = 0, value_count - 1
        if (hist_original(value) /= hist_sorted(value)) then
          verify = .false.
          return
        end if
      end do
    end do
  end function verify

end program main
