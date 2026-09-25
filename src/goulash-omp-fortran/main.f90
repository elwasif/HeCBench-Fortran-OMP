! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use, intrinsic :: iso_c_binding, only : c_char, c_double, c_long, c_null_char
  use omp_lib
  implicit none

  character(len=256) :: arg0, arg1, arg2
  integer(int64) :: iterations, n_cells, itime, i
  real(real64) :: kernel_mem_used, kernel_starttime, kernel_endtime, kernel_runtime
  real(real64), allocatable :: m_gate(:), m_gate_h(:), vm(:)
  logical :: ok

  interface
    function c_atol(str) bind(C, name='atol') result(value)
      import :: c_char, c_long
      character(kind=c_char), intent(in) :: str(*)
      integer(c_long) :: value
    end function c_atol

    function c_atof(str) bind(C, name='atof') result(value)
      import :: c_char, c_double
      character(kind=c_char), intent(in) :: str(*)
      real(c_double) :: value
    end function c_atof
  end interface

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    print '(3A)', 'Usage: ', trim(arg0), ' <Iterations> <Kernel_GBs_used>'
    print '(A)'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  iterations = int(c_atol(trim(arg1) // c_null_char), int64)
  kernel_mem_used = real(c_atof(trim(arg2) // c_null_char), real64)

  n_cells = int((kernel_mem_used * 1024.0_real64 * 1024.0_real64 * 1024.0_real64) / &
      (2.0_real64 * storage_size(0.0_real64) / 8.0_real64), int64)

  print '(A,I0)', 'Number of cells: ', n_cells

  allocate(m_gate(n_cells), m_gate_h(n_cells), vm(n_cells))
  m_gate = 0.0_real64
  m_gate_h = 0.0_real64
  vm = 0.0_real64
  kernel_starttime = 0.0_real64

  !$omp target data map(to: m_gate(1:n_cells), vm(1:n_cells))
  do itime = 0_int64, iterations
    if (itime == 1_int64) then
      !$omp target update from(m_gate(1:n_cells))
      kernel_starttime = omp_get_wtime()
    end if
    call gate(m_gate, n_cells, vm)
  end do
  kernel_endtime = omp_get_wtime()
  kernel_runtime = kernel_endtime - kernel_starttime
  print '(A,F0.6,A,I0,A)', 'total kernel time ', kernel_runtime, '(s) for ', iterations - 1_int64, ' iterations'
  !$omp end target data

  call reference(m_gate_h, n_cells, vm)

  ok = .true.
  do i = 1_int64, n_cells
    if (abs(m_gate(i) - m_gate_h(i)) > 1.0e-6_real64) then
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(m_gate, m_gate_h, vm)

contains

  subroutine gate(m_gate, n_cells, vm)
    real(real64), intent(inout) :: m_gate(:)
    integer(int64), intent(in) :: n_cells
    real(real64), intent(in) :: vm(:)
    integer(int64) :: idx
    real(real64) :: sum1, sum2, x, mhu, taur

    !$omp target teams distribute parallel do thread_limit(256) private(sum1, sum2, x, mhu, taur)
    do idx = 1_int64, n_cells
      x = vm(idx)
      call eval_gate_terms(x, mhu, taur)
      m_gate(idx) = m_gate(idx) + (mhu - m_gate(idx)) * (1.0_real64 - exp(-taur))
    end do
    !$omp end target teams distribute parallel do
  end subroutine gate

  subroutine reference(m_gate, n_cells, vm)
    real(real64), intent(inout) :: m_gate(:)
    integer(int64), intent(in) :: n_cells
    real(real64), intent(in) :: vm(:)
    integer(int64) :: idx
    real(real64) :: mhu, taur

    do idx = 1_int64, n_cells
      call eval_gate_terms(vm(idx), mhu, taur)
      m_gate(idx) = m_gate(idx) + (mhu - m_gate(idx)) * (1.0_real64 - exp(-taur))
    end do
  end subroutine reference

  subroutine eval_gate_terms(x, mhu, taur)
    real(real64), intent(in) :: x
    real(real64), intent(out) :: mhu, taur
    integer :: j
    real(real64) :: sum1, sum2
    real(real64), parameter :: mhu_a(15) = [ &
        9.9632117206253790e-01_real64, 4.0825738726469545e-02_real64, &
        6.3401613233199589e-04_real64, 4.4158436861700431e-06_real64, &
        1.1622058324043520e-08_real64, 1.0000000000000000e+00_real64, &
        4.0568375699663400e-02_real64, 6.4216825832642788e-04_real64, &
        4.2661664422410096e-06_real64, 1.3559930396321903e-08_real64, &
       -1.3573468728873069e-11_real64,-4.2594802366702580e-13_real64, &
        7.6779952208246166e-15_real64, 1.4260675804433780e-16_real64, &
       -2.6656212072499249e-18_real64]
    real(real64), parameter :: tau_a(18) = [ &
        1.7765862602413648e+01_real64*0.02_real64, 5.0010202770602419e-02_real64*0.02_real64, &
       -7.8002064070783474e-04_real64*0.02_real64,-6.9399661775931530e-05_real64*0.02_real64, &
        1.6936588308244311e-06_real64*0.02_real64, 5.4629017090963798e-07_real64*0.02_real64, &
       -1.3805420990037933e-08_real64*0.02_real64,-8.0678945216155694e-10_real64*0.02_real64, &
        1.6209833004622630e-11_real64*0.02_real64, 6.5130101230170358e-13_real64*0.02_real64, &
       -6.9931705949674988e-15_real64*0.02_real64,-3.1161210504114690e-16_real64*0.02_real64, &
        5.0166191902609083e-19_real64*0.02_real64, 7.8608831661430381e-20_real64*0.02_real64, &
        4.3936315597226053e-22_real64*0.02_real64,-7.0535966258003289e-24_real64*0.02_real64, &
       -9.0473475495087118e-26_real64*0.02_real64,-2.9878427692323621e-28_real64*0.02_real64]

    sum1 = 0.0_real64
    do j = 5, 1, -1
      sum1 = mhu_a(j) + x * sum1
    end do
    sum2 = 0.0_real64
    do j = 15, 6, -1
      sum2 = mhu_a(j) + x * sum2
    end do
    mhu = sum1 / sum2

    sum1 = 0.0_real64
    do j = 18, 1, -1
      sum1 = tau_a(j) + x * sum1
    end do
    taur = sum1
  end subroutine eval_gate_terms

end program main
