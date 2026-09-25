! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
  implicit none

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

  integer :: m, nsystems, block_size, repeat
  integer(int64) :: matrix_size
  real(real64), allocatable :: u_seq(:), u_thomas(:), u_input(:)
  real(real64), allocatable :: d_seq(:), d_thomas(:), d_input(:)
  real(real64), allocatable :: l_seq(:), l_thomas(:), l_input(:)
  real(real64), allocatable :: rhs_seq(:), rhs_thomas(:), rhs_input(:)
  real(real64), allocatable :: rhs_seq_output(:), rhs_seq_interleave(:)
  real(real64), allocatable :: params_u(:), params_l(:), params_d(:), params_rhs(:)
  real(real64) :: start_time, end_time, elapsed, error
  integer :: i, j, iter

  if (command_argument_count() /= 4) then
    write(*,'(A)') 'Usage: %s [system size] [#systems] [thread block size] [repeat]'
    stop -1
  end if

  m = read_int_arg(1)
  nsystems = read_int_arg(2)
  block_size = read_int_arg(3)
  repeat = read_int_arg(4)
  matrix_size = int(m, int64) * int(nsystems, int64)

  allocate(params_u(m), params_l(m), params_d(m), params_rhs(m))
  call load_thomas_matrix_syn(m, params_u, params_l, params_d, params_rhs)

  allocate(u_seq(matrix_size), u_thomas(matrix_size), u_input(matrix_size))
  allocate(d_seq(matrix_size), d_thomas(matrix_size), d_input(matrix_size))
  allocate(l_seq(matrix_size), l_thomas(matrix_size), l_input(matrix_size))
  allocate(rhs_seq(matrix_size), rhs_thomas(matrix_size), rhs_input(matrix_size))
  allocate(rhs_seq_output(matrix_size), rhs_seq_interleave(matrix_size))

  call initialize_systems(m, nsystems, params_u, params_l, params_d, params_rhs, &
      u_seq, u_input, d_seq, d_input, l_seq, l_input, rhs_seq, rhs_input)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call solve_seq(l_seq, d_seq, u_seq, rhs_seq, m, nsystems)
  end do
  end_time = omp_get_wtime()
  elapsed = (end_time - start_time) * 1000.0_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average serial execution time: ', elapsed, ' (ms)'

  rhs_seq_output = rhs_seq

  call initialize_systems(m, nsystems, params_u, params_l, params_d, params_rhs, &
      u_seq, u_input, d_seq, d_input, l_seq, l_input, rhs_seq, rhs_input)

  do i = 1, m
    do j = 1, nsystems
      u_thomas((i - 1) * nsystems + j) = u_input((j - 1) * m + i)
      l_thomas((i - 1) * nsystems + j) = l_input((j - 1) * m + i)
      d_thomas((i - 1) * nsystems + j) = d_input((j - 1) * m + i)
      rhs_thomas((i - 1) * nsystems + j) = rhs_input((j - 1) * m + i)
      rhs_seq_interleave((i - 1) * nsystems + j) = rhs_seq_output((j - 1) * m + i)
    end do
  end do

  !$omp target data map(to: l_thomas(1:matrix_size), d_thomas(1:matrix_size), u_thomas(1:matrix_size)) &
  !$omp& map(tofrom: rhs_thomas(1:matrix_size))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call solve_device(l_thomas, d_thomas, u_thomas, rhs_thomas, m, nsystems, block_size)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  elapsed = (end_time - start_time) * 1000.0_real64 / real(repeat, real64)
  write(*,'(A,F0.6,A)') 'Average kernel execution time: ', elapsed, ' (ms)'

  error = calc_error(rhs_seq_interleave, rhs_thomas)
  write(*,'(A,ES12.6)') 'Maximum error: ', error

  deallocate(params_u, params_l, params_d, params_rhs)
  deallocate(u_seq, u_thomas, u_input, d_seq, d_thomas, d_input)
  deallocate(l_seq, l_thomas, l_input, rhs_seq, rhs_thomas, rhs_input)
  deallocate(rhs_seq_output, rhs_seq_interleave)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  subroutine load_thomas_matrix_syn(size, upper, lower, diagonal, rhs)
    integer, intent(in) :: size
    real(real64), intent(out) :: upper(:), lower(:), diagonal(:), rhs(:)
  integer :: idx
    call c_srand(1_c_int)
    do idx = 1, size
      upper(idx) = frand(-2.0_real64, 2.0_real64)
      lower(idx) = frand(-2.0_real64, 2.0_real64)
      diagonal(idx) = frand(5.0_real64, 10.0_real64)
      rhs(idx) = frand(-2.0_real64, 2.0_real64)
    end do
  end subroutine load_thomas_matrix_syn

  real(real64) function frand(fmin, fmax) result(value)
    real(real64), intent(in) :: fmin, fmax
    value = fmin + real(c_rand(), real64) / 2147483647.0_real64 * (fmax - fmin)
  end function frand

  subroutine initialize_systems(m, nsystems, params_u, params_l, params_d, params_rhs, &
      u_seq, u_input, d_seq, d_input, l_seq, l_input, rhs_seq, rhs_input)
    integer, intent(in) :: m, nsystems
    real(real64), intent(in) :: params_u(:), params_l(:), params_d(:), params_rhs(:)
    real(real64), intent(out) :: u_seq(:), u_input(:), d_seq(:), d_input(:)
    real(real64), intent(out) :: l_seq(:), l_input(:), rhs_seq(:), rhs_input(:)
    integer :: i, j, idx
    do i = 1, nsystems
      do j = 1, m
        idx = (i - 1) * m + j
        u_seq(idx) = params_u(j)
        u_input(idx) = params_u(j)
        d_seq(idx) = params_d(j)
        d_input(idx) = params_d(j)
        l_seq(idx) = params_l(j)
        l_input(idx) = params_l(j)
        rhs_seq(idx) = params_rhs(j)
        rhs_input(idx) = params_rhs(j)
      end do
    end do
  end subroutine initialize_systems

  subroutine solve_seq(lower, diagonal, upper, rhs, m, nsystems)
    real(real64), intent(in) :: lower(:), diagonal(:)
    real(real64), intent(inout) :: upper(:), rhs(:)
    integer, intent(in) :: m, nsystems
    integer :: sys, first, last, idx
    do sys = 1, nsystems
      first = (sys - 1) * m + 1
      last = first + m - 1
      upper(first) = upper(first) / diagonal(first)
      rhs(first) = rhs(first) / diagonal(first)
      do idx = first + 1, last - 1
        upper(idx) = upper(idx) / (diagonal(idx) - lower(idx) * upper(idx - 1))
        rhs(idx) = (rhs(idx) - lower(idx) * rhs(idx - 1)) / (diagonal(idx) - lower(idx) * upper(idx - 1))
      end do
      rhs(last) = (rhs(last) - lower(last) * rhs(last - 1)) / (diagonal(last) - lower(last) * upper(last - 1))
      do idx = last - 1, first, -1
        rhs(idx) = rhs(idx) - upper(idx) * rhs(idx + 1)
      end do
    end do
  end subroutine solve_seq

  subroutine solve_device(lower, diagonal, upper, rhs, m, nsystems, block_size)
    real(real64), intent(in) :: lower(:), diagonal(:)
    real(real64), intent(inout) :: upper(:), rhs(:)
    integer, intent(in) :: m, nsystems, block_size
    integer :: tid, first, last, idx
    !$omp target teams distribute parallel do thread_limit(block_size) private(first, last, idx)
    do tid = 1, nsystems
      first = tid
      last = nsystems * (m - 1) + tid
      upper(first) = upper(first) / diagonal(first)
      rhs(first) = rhs(first) / diagonal(first)
      do idx = first + nsystems, last - nsystems, nsystems
        upper(idx) = upper(idx) / (diagonal(idx) - lower(idx) * upper(idx - nsystems))
        rhs(idx) = (rhs(idx) - lower(idx) * rhs(idx - nsystems)) / &
            (diagonal(idx) - lower(idx) * upper(idx - nsystems))
      end do
      rhs(last) = (rhs(last) - lower(last) * rhs(last - nsystems)) / &
          (diagonal(last) - lower(last) * upper(last - nsystems))
      do idx = last - nsystems, first, -nsystems
        rhs(idx) = rhs(idx) - upper(idx) * rhs(idx + nsystems)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine solve_device

  real(real64) function calc_error(src, dst) result(error)
    real(real64), intent(in) :: src(:), dst(:)
    integer :: idx
    real(real64) :: diff
    error = 0.0_real64
    do idx = 1, size(src)
      diff = abs(abs(src(idx)) - abs(dst(idx)))
      if (error < diff) error = diff
    end do
  end function calc_error

end program main
