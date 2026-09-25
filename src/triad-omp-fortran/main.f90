! SPDX-License-Identifier: CC0-1.0
program main
  use iso_fortran_env, only: real32, real64
  use iso_c_binding, only: c_double, c_long
  use omp_lib
  implicit none

  integer, parameter :: n_sizes = 9
  integer, parameter :: block_sizes(n_sizes) = [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
  integer, parameter :: mem_size = 16384
  integer, parameter :: num_max_floats = 1024 * mem_size / 4
  integer, parameter :: half_num_floats = num_max_floats / 2
  integer, parameter :: max_block_size = 16384 * 1024
  integer, parameter :: block_size = 128
  real(real32), parameter :: scalar = 1.75_real32
  integer(c_long), parameter :: triad_seed = 8650341_c_long

  interface
    subroutine c_srand48(seedval) bind(C, name='srand48')
      import :: c_long
      integer(c_long), value :: seedval
    end subroutine c_srand48

    function c_drand48() bind(C, name='drand48') result(value)
      import :: c_double
      real(c_double) :: value
    end function c_drand48
  end interface

  logical :: verbose
  integer :: n_passes

  call parse_args(verbose, n_passes)
  call run_benchmark(verbose, n_passes)

contains

  subroutine parse_args(verbose, n_passes)
    logical, intent(out) :: verbose
    integer, intent(out) :: n_passes
    integer :: argc, i, n
    character(len=256) :: arg
    logical :: ok, help_requested

    verbose = .false.
    n_passes = 10
    ok = .true.
    help_requested = .false.
    argc = command_argument_count()
    i = 1
    do while (i <= argc .and. ok)
      call get_command_argument(i, arg, length=n)
      select case (arg(:n))
      case ('-v', '--verbose')
        verbose = .true.
      case ('-n', '--passes')
        i = i + 1
        if (i <= argc) then
          call get_command_argument(i, arg, length=n)
          call parse_passes_value(arg(:n), n_passes, ok)
        else
          print '(a)', 'failure, option: --passes with no value'
          print '(a)', 'Ignoring remaining options'
          ok = .false.
        end if
      case ('-c', '--configFile')
        i = i + 1
        if (i <= argc) then
          call get_command_argument(i, arg, length=n)
          ok = parse_config_file(arg(:n), verbose, n_passes, help_requested)
        else
          print '(a)', 'failure, option: --configFile with no value'
          print '(a)', 'Ignoring remaining options'
          ok = .false.
        end if
      case ('-h', '--help')
        help_requested = .true.
        ok = .false.
      case default
        ok = parse_unknown_or_short_option(arg(:n), i, argc, verbose, n_passes, help_requested)
      end select
      i = i + 1
    end do
    if (.not. ok) then
      call usage()
      if (help_requested) stop
      stop 1
    end if
  end subroutine parse_args

  subroutine parse_passes_value(text, n_passes, ok)
    character(len=*), intent(in) :: text
    integer, intent(out) :: n_passes
    logical, intent(out) :: ok
    integer :: ios

    read(text, *, iostat=ios) n_passes
    ok = ios == 0
    if (.not. ok) then
      print '(a)', 'failure, option: --passes with malformed value'
      print '(a)', 'Ignoring remaining options'
    end if
  end subroutine parse_passes_value

  recursive logical function parse_config_file(file_name, verbose, n_passes, help_requested) result(ok)
    character(len=*), intent(in) :: file_name
    logical, intent(inout) :: verbose
    integer, intent(inout) :: n_passes
    logical, intent(inout) :: help_requested
    character(len=512) :: line, key, value
    integer :: unit, ios, split_at

    ok = .true.
    open(newunit=unit, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print '(a)', 'Bad config file'
      ok = .false.
      return
    end if

    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      split_at = scan(line, ' ')
      if (split_at == 0) then
        key = trim(line)
        value = ''
      else
        key = trim(line(:split_at - 1))
        value = adjustl(line(split_at + 1:))
      end if

      select case (trim(key))
      case ('verbose')
        verbose = .true.
      case ('passes')
        call parse_passes_value(trim(value), n_passes, ok)
      case ('help')
        help_requested = .true.
        ok = .false.
      case ('configFile')
        ok = parse_config_file(trim(value), verbose, n_passes, help_requested)
      case default
        print '(a,a)', 'Option not recognized: --', trim(key)
        print '(a)', 'Ignoring remaining options'
        ok = .false.
      end select
      if (.not. ok) exit
    end do
    close(unit)
  end function parse_config_file

  logical function parse_unknown_or_short_option(arg, i, argc, verbose, n_passes, help_requested) result(ok)
    character(len=*), intent(in) :: arg
    integer, intent(inout) :: i
    integer, intent(in) :: argc
    logical, intent(inout) :: verbose
    integer, intent(inout) :: n_passes
    logical, intent(inout) :: help_requested
    character(len=256) :: value
    integer :: n, p, nopts

    ok = .true.
    n = len_trim(arg)
    if (n == 0 .or. arg(1:1) /= '-') then
      print '(a,a)', 'failure, no leading - in option: ', trim(arg)
      print '(a)', 'Ignoring remaining options'
      ok = .false.
      return
    end if

    if (n >= 2 .and. arg(1:2) == '--') then
      print '(a,a)', 'Option not recognized: ', trim(arg)
      print '(a)', 'Ignoring remaining options'
      ok = .false.
      return
    end if

    nopts = n - 1
    do p = 1, nopts
      select case (arg(p + 1:p + 1))
      case ('v')
        verbose = .true.
      case ('n')
        if (i + 1 > argc .or. p < nopts) then
          print '(a)', 'failure, option: -n with no value'
          print '(a)', 'Ignoring remaining options'
          ok = .false.
          return
        end if
        i = i + 1
        call get_command_argument(i, value, length=n)
        call parse_passes_value(value(:n), n_passes, ok)
        if (.not. ok) return
      case ('c')
        if (i + 1 > argc .or. p < nopts) then
          print '(a)', 'failure, option: -c with no value'
          print '(a)', 'Ignoring remaining options'
          ok = .false.
          return
        end if
        i = i + 1
        call get_command_argument(i, value, length=n)
        ok = parse_config_file(value(:n), verbose, n_passes, help_requested)
        if (.not. ok) return
      case ('h')
        help_requested = .true.
        ok = .false.
        return
      case default
        print '(a,a,a)', 'Option: ', trim(arg), ' not recognized.'
        print '(a)', 'Ignoring remaining options'
        ok = .false.
        return
      end select
    end do
  end function parse_unknown_or_short_option

  subroutine usage()
    print '(a)', 'Usage: ./main [options]'
    print '(a)', '  -c, --configFile <file> specify configuration file'
    print '(a)', '  -v, --verbose        enable verbose output'
    print '(a)', '  -n, --passes <N>     specify number of passes'
    print '(a)', '  -h, --help           print this usage'
  end subroutine usage

  subroutine run_benchmark(verbose, n_passes)
    logical, intent(in) :: verbose
    integer, intent(in) :: n_passes
    real(real32), allocatable :: h_mem(:)
    real(real32), allocatable :: a0(:), b0(:), c0(:), a1(:), b1(:), c1(:)
    integer :: i, j, pass, elems_in_block, crt_idx, block_idx
    integer :: gid
    logical :: curr_stream, ok
    real(real64) :: start_time, elapsed, triad, bdwth

    allocate(h_mem(0:num_max_floats - 1))
    allocate(a0(0:max_block_size - 1), b0(0:max_block_size - 1), c0(0:max_block_size - 1))
    allocate(a1(0:max_block_size - 1), b1(0:max_block_size - 1), c1(0:max_block_size - 1))

    call c_srand48(triad_seed)

    !$omp target data map(alloc: a0(0:max_block_size - 1), b0(0:max_block_size - 1), &
    !$omp& c0(0:max_block_size - 1), a1(0:max_block_size - 1), b1(0:max_block_size - 1), &
    !$omp& c1(0:max_block_size - 1))
    do i = 1, n_sizes
      do j = 0, num_max_floats - 1
        c0(j) = 0.0_real32
        c1(j) = 0.0_real32
      end do

      do j = 0, half_num_floats - 1
        a0(j) = real(c_drand48() * 10.0_c_double, real32)
        a0(half_num_floats + j) = a0(j)
        b0(j) = a0(j)
        b0(half_num_floats + j) = a0(j)
        a1(j) = a0(j)
        a1(half_num_floats + j) = a0(j)
        b1(j) = a0(j)
        b1(half_num_floats + j) = a0(j)
      end do

      elems_in_block = block_sizes(i) * 1024 / 4
      if (verbose) then
        write(*, '(a,i0,a,i0,a)') '>> Executing Triad with vectors of length ', &
            num_max_floats, ' and block size of ', elems_in_block, ' elements.'
        write(*, '(a,i0,a)') 'Block: ', block_sizes(i), 'KB'
      end if

      crt_idx = 0
      start_time = omp_get_wtime()
      do pass = 0, n_passes - 1
        !$omp target update to(a0(0:elems_in_block - 1)) nowait
        !$omp target update to(b0(0:elems_in_block - 1)) nowait

        !$omp target teams distribute parallel do thread_limit(block_size) nowait
        do gid = 0, elems_in_block - 1
          c0(gid) = a0(gid) + scalar * b0(gid)
        end do
        !$omp end target teams distribute parallel do

        if (elems_in_block < num_max_floats) then
          !$omp target update to(a1(elems_in_block:2 * elems_in_block - 1)) nowait
          !$omp target update to(b1(elems_in_block:2 * elems_in_block - 1)) nowait
        end if

        block_idx = 1
        curr_stream = .true.
        do while (crt_idx < num_max_floats)
          curr_stream = iand(block_idx, 1) == 1
          if (curr_stream) then
            !$omp target update from(c0(crt_idx:crt_idx + elems_in_block - 1)) nowait
          else
            !$omp target update from(c1(crt_idx:crt_idx + elems_in_block - 1)) nowait
          end if

          crt_idx = crt_idx + elems_in_block
          if (crt_idx < num_max_floats) then
            if (curr_stream) then
              !$omp target teams distribute parallel do thread_limit(block_size) nowait
              do gid = 0, elems_in_block - 1
                c1(crt_idx + gid) = a1(crt_idx + gid) + scalar * b1(crt_idx + gid)
              end do
              !$omp end target teams distribute parallel do
            else
              !$omp target teams distribute parallel do thread_limit(block_size) nowait
              do gid = 0, elems_in_block - 1
                c0(crt_idx + gid) = a0(crt_idx + gid) + scalar * b0(crt_idx + gid)
              end do
              !$omp end target teams distribute parallel do
            end if
          end if

          if (crt_idx + elems_in_block < num_max_floats) then
            if (curr_stream) then
              !$omp target update to(a0(crt_idx + elems_in_block:crt_idx + 2 * elems_in_block - 1)) nowait
              !$omp target update to(b0(crt_idx + elems_in_block:crt_idx + 2 * elems_in_block - 1)) nowait
            else
              !$omp target update to(a1(crt_idx + elems_in_block:crt_idx + 2 * elems_in_block - 1)) nowait
              !$omp target update to(b1(crt_idx + elems_in_block:crt_idx + 2 * elems_in_block - 1)) nowait
            end if
          end if
          block_idx = block_idx + 1
        end do
      end do

      elapsed = omp_get_wtime() - start_time
      if (verbose) then
        triad = (real(num_max_floats, real64) * 2.0_real64 * real(n_passes, real64)) / (elapsed * 1.0e9_real64)
        bdwth = (real(num_max_floats, real64) * 4.0_real64 * 3.0_real64 * real(n_passes, real64)) / &
            (elapsed * 1000.0_real64 * 1000.0_real64 * 1000.0_real64)
        write(*, '(a,f0.6,a)') 'Average TriadFlops ', triad, ' GFLOPS/s'
        write(*, '(a,f0.6,a)') 'Average TriadBdwth ', bdwth, ' GB/s'
      end if

      ok = .true.
      do j = 0, num_max_floats - 1, elems_in_block
        if (iand(j / elems_in_block, 1) == 0) then
          h_mem(j:j + elems_in_block - 1) = c0(j:j + elems_in_block - 1)
        else
          h_mem(j:j + elems_in_block - 1) = c1(j:j + elems_in_block - 1)
        end if
      end do

      do j = 0, half_num_floats - 1
        if (h_mem(j) /= h_mem(j + half_num_floats)) then
          write(*, '(a,i0,a,es13.6,a,i0,a,es13.6,a)') 'hostMem[', j, ']=', h_mem(j), &
              ' is different from its twin element hostMem[', j + half_num_floats, ']: ', &
              h_mem(j + half_num_floats), 'stopping check'
          ok = .false.
          exit
        end if
      end do

      if (ok) then
        print '(a)', 'PASS'
      else
        print '(a)', 'FAIL'
      end if
    end do
    !$omp end target data

    deallocate(h_mem, a0, b0, c0, a1, b1, c1)
  end subroutine run_benchmark

end program main
