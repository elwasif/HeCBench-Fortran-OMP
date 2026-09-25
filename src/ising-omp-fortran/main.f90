! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int64, real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

  real(real32), parameter :: tcrit = 2.26918531421_real32
  integer, parameter :: threads = 128

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

  integer(int64) :: nx, ny, seed
  integer :: nwarmup, niters, total_cells, i
  real(real32) :: alpha, inv_temp
  real(real64) :: t0, t1, duration_us, updates_per_ns
  real(real32), allocatable :: randvals(:)
  integer(int8), allocatable :: lattice_b(:), lattice_w(:), lattice_b_ref(:), lattice_w_ref(:)
  logical :: ok

  nx = 5120_int64
  ny = 5120_int64
  alpha = 0.1_real32
  nwarmup = 100
  niters = 1000
  seed = 1234_int64

  call parse_args(nx, ny, alpha, seed, nwarmup, niters)

  if (mod(nx, 2_int64) /= 0 .or. mod(ny, 2_int64) /= 0) then
    write(0,'(A)') 'ERROR: Lattice dimensions must be even values.'
    stop 1
  end if
  if (nx * ny / 2_int64 > int(huge(total_cells), int64)) then
    write(0,'(A)') 'ERROR: Lattice dimensions are too large for this Fortran port.'
    stop 1
  end if

  inv_temp = 1.0_real32 / (alpha * tcrit)
  total_cells = int(nx * ny / 2_int64)
  allocate(randvals(total_cells), lattice_b(total_cells), lattice_w(total_cells))
  allocate(lattice_b_ref(total_cells), lattice_w_ref(total_cells))

  call fill_randvals(randvals, seed)

  !$omp target data map(to: randvals(1:total_cells)) map(alloc: lattice_b(1:total_cells), lattice_w(1:total_cells))
  call init_spins(lattice_b, randvals, int(nx), int(ny / 2_int64))
  call init_spins(lattice_w, randvals, int(nx), int(ny / 2_int64))

  print '(A)', 'Starting warmup...'
  do i = 1, nwarmup
    call update(lattice_b, lattice_w, randvals, inv_temp, int(nx), int(ny))
  end do

  print '(A)', 'Starting trial iterations...'
  t0 = omp_get_wtime()
  do i = 1, niters
    call update(lattice_b, lattice_w, randvals, inv_temp, int(nx), int(ny))
  end do
  t1 = omp_get_wtime()
  duration_us = (t1 - t0) * 1.0e6_real64
  !$omp target update from(lattice_b(1:total_cells), lattice_w(1:total_cells))
  !$omp end target data

  print '(A)', 'REPORT:'
  print '(A,I0)', char(9)//'nGPUs: ', 1
  print '(A,F8.6,A,F8.6)', char(9)//'temperature: ', alpha, ' * ', tcrit
  print '(A,I0)', char(9)//'seed: ', seed
  print '(A,I0)', char(9)//'warmup iterations: ', nwarmup
  print '(A,I0)', char(9)//'trial iterations: ', niters
  print '(A,I0,A,I0)', char(9)//'lattice dimensions: ', nx, ' x ', ny
  print '(A,F0.6,A)', char(9)//'elapsed time: ', duration_us * 1.0e-6_real64, ' sec'
  updates_per_ns = real(nx * ny, real64) * real(niters, real64) / duration_us * 1.0e-3_real64
  print '(A,F0.6)', char(9)//'updates per ns: ', updates_per_ns

  print '(A)', 'Starting verification iterations ...'
  call init_spins_ref(lattice_b_ref, randvals, int(nx), int(ny / 2_int64))
  call init_spins_ref(lattice_w_ref, randvals, int(nx), int(ny / 2_int64))
  do i = 1, nwarmup + niters
    call update_ref(lattice_b_ref, lattice_w_ref, randvals, inv_temp, int(nx), int(ny))
  end do

  ok = all(lattice_b == lattice_b_ref) .and. all(lattice_w == lattice_w_ref)
  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(randvals, lattice_b, lattice_w, lattice_b_ref, lattice_w_ref)

contains

  subroutine parse_args(nx, ny, alpha, seed, nwarmup, niters)
    integer(int64), intent(inout) :: nx, ny, seed
    real(real32), intent(inout) :: alpha
    integer, intent(inout) :: nwarmup, niters
    integer :: arg, argc
    character(len=256) :: key, value

    argc = command_argument_count()
    arg = 1
    do while (arg <= argc)
      call get_command_argument(arg, key)
      if (trim(key) == '-h' .or. trim(key) == '--help') call usage()
      if (arg == argc) then
        write(0,'(A,A)') 'option requires an argument: ', trim(key)
        stop 1
      end if
      call get_command_argument(arg + 1, value)
      select case (trim(key))
      case ('-x', '--lattice-n')
        read(value, *) nx
      case ('-y', '--lattice-m')
        read(value, *) ny
      case ('--alpha')
        read(value, *) ny
      case ('-a')
        read(value, *) alpha
      case ('-s', '--seed')
        read(value, *) seed
      case ('-w', '--nwarmup')
        read(value, *) nwarmup
      case ('-n', '--niters')
        read(value, *) niters
      case default
        write(0,'(A,A)') 'unknown option: ', trim(key)
        stop 1
      end select
      arg = arg + 2
    end do
  end subroutine parse_args

  subroutine usage()
    print '(A)', 'Usage: ./main [options]'
    print '(A)', 'options:'
    print '(A)', char(9)//'-x|--lattice-n <LATTICE_N>'
    print '(A)', char(9)//'-y|--lattice_m <LATTICE_M>'
    print '(A)', char(9)//'-w|--nwarmup <NWARMUP>'
    print '(A)', char(9)//'-n|--niters <NITERS>'
    print '(A)', char(9)//'-a|--alpha <ALPHA>'
    print '(A)', char(9)//'-s|--seed <SEED>'
    stop
  end subroutine usage

  subroutine fill_randvals(randvals, seed)
    real(real32), intent(out) :: randvals(:)
    integer(int64), intent(in) :: seed
    integer :: i
    call c_srand(int(seed, c_int))
    do i = 1, size(randvals)
      randvals(i) = real(c_rand(), real32) / 2147483647.0_real32
    end do
  end subroutine fill_randvals

  subroutine init_spins(lattice, randvals, nx, half_ny)
    integer(int8), intent(inout) :: lattice(:)
    real(real32), intent(in) :: randvals(:)
    integer, intent(in) :: nx, half_ny
    integer :: tid
    !$omp target teams distribute parallel do simd thread_limit(threads)
    do tid = 1, nx * half_ny
      if (randvals(tid) < 0.5_real32) then
        lattice(tid) = -1_int8
      else
        lattice(tid) = 1_int8
      end if
    end do
    !$omp end target teams distribute parallel do simd
  end subroutine init_spins

  subroutine update(lattice_b, lattice_w, randvals, inv_temp, nx, ny)
    integer(int8), intent(inout) :: lattice_b(:), lattice_w(:)
    real(real32), intent(in) :: randvals(:), inv_temp
    integer, intent(in) :: nx, ny
    call update_lattice(lattice_b, lattice_w, randvals, inv_temp, nx, ny / 2, .true.)
    call update_lattice(lattice_w, lattice_b, randvals, inv_temp, nx, ny / 2, .false.)
  end subroutine update

  subroutine update_lattice(lattice, op_lattice, randvals, inv_temp, nx, half_ny, is_black)
    integer(int8), intent(inout) :: lattice(:)
    integer(int8), intent(in) :: op_lattice(:)
    real(real32), intent(in) :: randvals(:), inv_temp
    integer, intent(in) :: nx, half_ny
    logical, intent(in) :: is_black
    integer :: i, j, ipp, inn, jpp, jnn, joff, idx
    integer :: nn_sum
    integer(int8) :: lij
    real(real32) :: acceptance_ratio

    !$omp target teams distribute parallel do collapse(2) thread_limit(threads) &
    !$omp& private(ipp, inn, jpp, jnn, joff, idx, nn_sum, lij, acceptance_ratio)
    do i = 0, nx - 1
      do j = 0, half_ny - 1
        ipp = merge(i + 1, 0, i + 1 < nx)
        inn = merge(i - 1, nx - 1, i - 1 >= 0)
        jpp = merge(j + 1, 0, j + 1 < half_ny)
        jnn = merge(j - 1, half_ny - 1, j - 1 >= 0)
        if (is_black) then
          joff = merge(jpp, jnn, mod(i, 2) /= 0)
        else
          joff = merge(jnn, jpp, mod(i, 2) /= 0)
        end if
        idx = i * half_ny + j + 1
        nn_sum = int(op_lattice(inn * half_ny + j + 1)) + int(op_lattice(idx)) + &
                 int(op_lattice(ipp * half_ny + j + 1)) + int(op_lattice(i * half_ny + joff + 1))
        lij = lattice(idx)
        acceptance_ratio = exp(-2.0_real32 * inv_temp * real(nn_sum, real32) * real(lij, real32))
        if (randvals(idx) < acceptance_ratio) lattice(idx) = -lij
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine update_lattice

  subroutine init_spins_ref(lattice, randvals, nx, half_ny)
    integer(int8), intent(out) :: lattice(:)
    real(real32), intent(in) :: randvals(:)
    integer, intent(in) :: nx, half_ny
    integer :: tid
    do tid = 1, nx * half_ny
      lattice(tid) = merge(-1_int8, 1_int8, randvals(tid) < 0.5_real32)
    end do
  end subroutine init_spins_ref

  subroutine update_ref(lattice_b, lattice_w, randvals, inv_temp, nx, ny)
    integer(int8), intent(inout) :: lattice_b(:), lattice_w(:)
    real(real32), intent(in) :: randvals(:), inv_temp
    integer, intent(in) :: nx, ny
    call update_lattice_ref(lattice_b, lattice_w, randvals, inv_temp, nx, ny / 2, .true.)
    call update_lattice_ref(lattice_w, lattice_b, randvals, inv_temp, nx, ny / 2, .false.)
  end subroutine update_ref

  subroutine update_lattice_ref(lattice, op_lattice, randvals, inv_temp, nx, half_ny, is_black)
    integer(int8), intent(inout) :: lattice(:)
    integer(int8), intent(in) :: op_lattice(:)
    real(real32), intent(in) :: randvals(:), inv_temp
    integer, intent(in) :: nx, half_ny
    logical, intent(in) :: is_black
    integer :: i, j, ipp, inn, jpp, jnn, joff, idx, nn_sum
    integer(int8) :: lij
    real(real32) :: acceptance_ratio

    do i = 0, nx - 1
      do j = 0, half_ny - 1
        ipp = merge(i + 1, 0, i + 1 < nx)
        inn = merge(i - 1, nx - 1, i - 1 >= 0)
        jpp = merge(j + 1, 0, j + 1 < half_ny)
        jnn = merge(j - 1, half_ny - 1, j - 1 >= 0)
        if (is_black) then
          joff = merge(jpp, jnn, mod(i, 2) /= 0)
        else
          joff = merge(jnn, jpp, mod(i, 2) /= 0)
        end if
        idx = i * half_ny + j + 1
        nn_sum = int(op_lattice(inn * half_ny + j + 1)) + int(op_lattice(idx)) + &
                 int(op_lattice(ipp * half_ny + j + 1)) + int(op_lattice(i * half_ny + joff + 1))
        lij = lattice(idx)
        acceptance_ratio = exp(-2.0_real32 * inv_temp * real(nn_sum, real32) * real(lij, real32))
        if (randvals(idx) < acceptance_ratio) lattice(idx) = -lij
      end do
    end do
  end subroutine update_lattice_ref

end program main
