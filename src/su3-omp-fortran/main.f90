! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_signed_char, c_sizeof
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  use omp_lib
  implicit none

  integer(c_signed_char), parameter :: EVEN = 2_c_signed_char, ODD = 1_c_signed_char

  type, bind(C) :: Complx
    real(c_double) :: real
    real(c_double) :: imag
  end type Complx

  type, bind(C) :: su3_matrix
    real(c_double) :: e(2,3,3)
  end type su3_matrix

  type, bind(C) :: site
    real(c_double) :: link(2,3,3,4)
    integer(c_int) :: x, y, z, t
    integer(c_int) :: index
    integer(c_signed_char) :: parity
    integer(c_int) :: pad(10)
  end type site

  integer(int64) :: iterations, ldim, threads_per_group, warmups
  integer(int32) :: verbose, device
  integer :: argc

  iterations = 100_int64
  ldim = 32_int64
  threads_per_group = 128_int64
  warmups = 1_int64
  verbose = 1_int32
  device = -1_int32

  argc = command_argument_count()
  call parse_args(argc, iterations, ldim, threads_per_group, verbose, device, warmups)
  call run_benchmark(iterations, ldim, threads_per_group, verbose, warmups)

contains

  subroutine parse_args(argc, iterations, ldim, threads_per_group, verbose, device, warmups)
    integer, intent(in) :: argc
    integer(int64), intent(inout) :: iterations, ldim, threads_per_group, warmups
    integer(int32), intent(inout) :: verbose, device
    character(len=256) :: arg, value
    integer :: i

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('-h')
        write(*,'(A)') 'Usage: ./main [-i iterations] [-l lattice dimension] [-t threads per workgroup] [-d device] [-v verbosity level [0,1,2,3]] [-w warmups]'
        stop 1
      case ('-i')
        call get_command_argument(i + 1, value)
        read(value, *) iterations
        i = i + 1
      case ('-l')
        call get_command_argument(i + 1, value)
        read(value, *) ldim
        i = i + 1
      case ('-t')
        call get_command_argument(i + 1, value)
        read(value, *) threads_per_group
        i = i + 1
      case ('-v')
        call get_command_argument(i + 1, value)
        read(value, *) verbose
        i = i + 1
      case ('-d')
        call get_command_argument(i + 1, value)
        read(value, *) device
        i = i + 1
      case ('-w')
        call get_command_argument(i + 1, value)
        read(value, *) warmups
        i = i + 1
      case ('-n')
        i = i + 1
      end select
      i = i + 1
    end do
  end subroutine parse_args

  subroutine run_benchmark(iterations, ldim, threads_per_group, verbose, warmups)
    integer(int64), intent(in) :: iterations, ldim, warmups
    integer(int64), intent(inout) :: threads_per_group
    integer(int32), intent(in) :: verbose
    integer(int64) :: total_sites
    type(site), allocatable :: a(:), c(:)
    type(su3_matrix), allocatable :: b(:)
    real(real64) :: ttotal, tflop, memory_usage, checksum

    if (threads_per_group == 0_int64) threads_per_group = 36_int64
    total_sites = ldim * ldim * ldim * ldim

    allocate(a(total_sites), c(total_sites), b(4))

    call make_lattice(a, ldim, Complx(1.0_c_double, 0.0_c_double))
    call init_link(b, Complx(1.0_c_double / 3.0_c_double, 0.0_c_double))
    c = a

    if (verbose >= 1_int32) then
      write(*,'(A,I0,A)') 'Number of sites = ', ldim, '^4'
      write(*,'(A,I0,A,I0,A)') 'Executing ', iterations, ' iterations with ', warmups, ' warmups'
      if (threads_per_group /= 0_int64) write(*,'(A,I0)') 'Threads per group = ', threads_per_group
    end if

    ttotal = su3_mat_nn(a, b, c, total_sites, iterations, &
      threads_per_group, warmups, verbose, checksum)

    if (verbose >= 1_int32) write(*,'(A,F8.6,A)') 'Total kernel execution time = ', ttotal, ' (s)'
    tflop = real(iterations, real64) * real(total_sites, real64) * 864.0_real64
    write(*,'(A,F0.3)') 'Total GFLOP/s = ', tflop / ttotal / 1.0e9_real64
    memory_usage = real(c_sizeof(a(1)), real64) * real(size(a) + size(c), real64) + &
      real(c_sizeof(b(1)), real64) * real(size(b), real64)
    write(*,'(A,F0.3)') 'Total GByte/s (GPU memory)  = ', &
      real(iterations, real64) * memory_usage / ttotal / 1.0e9_real64

    if (abs(checksum - 36.0_real64) > 1.0e-6_real64) error stop 'verification failed'

    deallocate(b, c, a)
  end subroutine run_benchmark

  subroutine init_link(s, val)
    type(su3_matrix), intent(inout) :: s(:)
    type(Complx), intent(in) :: val
    integer :: j, k, l

    do j = 1, size(s)
      do k = 1, 3
        do l = 1, 3
          s(j)%e(1,l,k) = val%real
          s(j)%e(2,l,k) = val%imag
        end do
      end do
    end do
  end subroutine init_link

  subroutine init_site_link(link, val)
    real(c_double), intent(inout) :: link(2,3,3,4)
    type(Complx), intent(in) :: val
    integer :: j, k, l

    do j = 1, 4
      do k = 1, 3
        do l = 1, 3
          link(1,l,k,j) = val%real
          link(2,l,k,j) = val%imag
        end do
      end do
    end do
  end subroutine init_site_link

  subroutine make_lattice(s, n, val)
    type(site), intent(inout) :: s(:)
    integer(int64), intent(in) :: n
    type(Complx), intent(in) :: val
    integer(int64) :: nx, ny, nz, nt, x, y, z, t, i

    nx = n
    ny = n
    nz = n
    nt = n
    do t = 0_int64, nt - 1_int64
      i = t * nz * ny * nx + 1_int64
      do z = 0_int64, nz - 1_int64
        do y = 0_int64, ny - 1_int64
          do x = 0_int64, nx - 1_int64
            s(i)%x = int(x, c_int)
            s(i)%y = int(y, c_int)
            s(i)%z = int(z, c_int)
            s(i)%t = int(t, c_int)
            s(i)%index = int(x + nx * (y + ny * (z + nz * t)), c_int)
            if (mod(x + y + z + t, 2_int64) == 0_int64) then
              s(i)%parity = EVEN
            else
              s(i)%parity = ODD
            end if
            call init_site_link(s(i)%link, val)
            i = i + 1_int64
          end do
        end do
      end do
    end do
  end subroutine make_lattice

  real(real64) function su3_mat_nn(a, b, c, total_sites, iterations, &
      threads_per_group, warmups, verbose, checksum) result(ttotal)
    integer(int64), intent(in) :: total_sites, iterations, threads_per_group, warmups
    type(site), intent(in) :: a(total_sites)
    type(su3_matrix), intent(in) :: b(4)
    type(site), intent(inout) :: c(total_sites)
    integer(int32), intent(in) :: verbose
    real(real64), intent(out) :: checksum
    integer :: num_work_items, n_sites, id, i, j, k, l, m, iters
    real(real64) :: cc_real, cc_imag, tstart

    n_sites = int(total_sites)
    num_work_items = n_sites * int(threads_per_group)
    if (verbose >= 1_int32) then
      write(*,'(A,I0)') 'Number of teams = ', total_sites
      write(*,'(A,I0)') 'Threads per team = ', threads_per_group
      write(*,'(A,I0)') 'Number of work items = ', num_work_items
    end if

    !$omp target data map(to: a(1:n_sites), b(1:4)) map(from: c(1:n_sites))
    do iters = 0, int(iterations + warmups) - 1
      if (iters == int(warmups)) tstart = omp_get_wtime()
      !$omp target teams distribute parallel do num_teams(total_sites) thread_limit(threads_per_group) &
      !$omp& private(i,j,k,l,m,cc_real,cc_imag)
      do id = 0, num_work_items - 1
        i = id / 36 + 1
        if (i <= n_sites) then
          j = mod(id, 36) / 9 + 1
          k = mod(id, 9) / 3 + 1
          l = mod(id, 3) + 1
          cc_real = 0.0_real64
          cc_imag = 0.0_real64
          do m = 1, 3
            cc_real = cc_real + a(i)%link(1,m,k,j) * b(j)%e(1,l,m) - &
              a(i)%link(2,m,k,j) * b(j)%e(2,l,m)
            cc_imag = cc_imag + a(i)%link(1,m,k,j) * b(j)%e(2,l,m) + &
              a(i)%link(2,m,k,j) * b(j)%e(1,l,m)
          end do
          c(i)%link(1,l,k,j) = cc_real
          c(i)%link(2,l,k,j) = cc_imag
        end if
      end do
      !$omp end target teams distribute parallel do
    end do
    !$omp end target data

    ttotal = omp_get_wtime() - tstart

    checksum = 0.0_real64
    do i = 1, n_sites
      do j = 1, 4
        do k = 1, 3
          do l = 1, 3
            cc_real = 0.0_real64
            cc_imag = 0.0_real64
            do m = 1, 3
              cc_real = cc_real + a(i)%link(1,m,k,j) * b(j)%e(1,l,m) - &
                a(i)%link(2,m,k,j) * b(j)%e(2,l,m)
              cc_imag = cc_imag + a(i)%link(1,m,k,j) * b(j)%e(2,l,m) + &
                a(i)%link(2,m,k,j) * b(j)%e(1,l,m)
            end do
            if (abs(c(i)%link(1,l,k,j) - cc_real) > 1.0e-6_real64 .or. &
                abs(c(i)%link(2,l,k,j) - cc_imag) > 1.0e-6_real64) error stop 'verification failed'
            checksum = checksum + c(i)%link(1,l,k,j)
          end do
        end do
      end do
    end do
    checksum = checksum / real(total_sites, real64)

    if (abs(checksum - 36.0_real64) < 1.0e-6_real64) then
      write(*,'(A,I0)') 'Checksum SUCCESS... though please be diligent and check the following value is not NaN: checksum=', nint(checksum)
    else
      write(*,'(A)') 'Checksum FAILURE'
    end if
  end function su3_mat_nn

end program main
