! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  logical :: test_result

  call run_perf()
  test_result = run_test()

  if (test_result) then
    write(*,'(A)') 'Test assert completed, returned OK'
  else
    write(*,'(A)') 'Test assert completed, returned ERROR!'
    stop 1
  end if

contains

  subroutine test_kernel(num_teams, num_threads, n)
    integer, intent(in) :: num_teams, num_threads, n
    integer :: gid

    !$omp target teams distribute parallel do num_teams(num_teams) num_threads(num_threads)
    do gid = 0, n - 1
      if (.not. (gid < n)) error stop 'assert(gid < n)'
    end do
    !$omp end target teams distribute parallel do
  end subroutine test_kernel

  subroutine perf_kernel(num_teams, num_threads)
    integer, intent(in) :: num_teams, num_threads
    integer :: gid, n, s

    !$omp target teams num_teams(num_teams)
    !$omp parallel private(gid, n, s) num_threads(num_threads)
    gid = omp_get_team_num() * omp_get_num_threads() + omp_get_thread_num()
    if (.not. (gid <= omp_get_num_threads() * omp_get_num_teams())) error stop 'assert(gid <= omp_get_num_threads() * omp_get_num_teams())'

    s = 0
    do n = 1, gid
      s = s + 1
      if (.not. (s <= gid)) error stop 'assert(s <= gid)'
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine perf_kernel

  subroutine perf_kernel2(num_teams, num_threads)
    integer, intent(in) :: num_teams, num_threads
    integer :: gid, n, s

    !$omp target teams num_teams(num_teams)
    !$omp parallel private(gid, n, s) num_threads(num_threads)
    gid = omp_get_team_num() * omp_get_num_threads() + omp_get_thread_num()
    s = 0
    do n = 1, gid
      s = s + 1
      if (.not. (s <= gid)) error stop 'assert(s <= gid)'
    end do
    !$omp end parallel
    !$omp end target teams
  end subroutine perf_kernel2

  logical function run_test() result(ok)
    integer :: nblocks, nthreads

    nblocks = 2
    nthreads = 32

    write(*,'(A)') ''
    write(*,'(A)') 'Launch kernel to generate assertion failures'
    write(*,'(A)') ''
    write(*,'(A)') '-- Begin assert output'
    write(*,'(A)') ''

    call test_kernel(nblocks, nthreads, 60)

    write(*,'(A)') ''
    write(*,'(A)') '-- End assert output'
    write(*,'(A)') ''

    ok = .true.
  end function run_test

  subroutine run_perf()
    integer :: nblocks, nthreads
    real(real64) :: start_time, end_time, elapsed

    nblocks = 1000
    nthreads = 256

    write(*,'(A)') ''
    write(*,'(A)') 'Launch kernel to evaluate the impact of assertion on performance '
    write(*,'(A)') 'Each thread in the kernel executes threadID + 1 assertions'
    start_time = omp_get_wtime()
    call perf_kernel(nblocks, nthreads)
    end_time = omp_get_wtime()
    elapsed = end_time - start_time
    write(*,'(A,F8.6)') 'Kernel time : ', elapsed

    write(*,'(A)') 'Each thread in the kernel executes threadID assertions'
    call perf_kernel2(nblocks, nthreads)
    end_time = omp_get_wtime()
    elapsed = end_time - start_time
    write(*,'(A,F8.6)') 'Kernel time : ', elapsed
  end subroutine run_perf

end program main
