! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int8, int64, real32, real64
  use omp_lib
  implicit none

  type :: FloatComplex
    real(real32) :: x, y
  end type FloatComplex

  type :: DoubleComplex
    real(real64) :: x, y
  end type DoubleComplex

  character(len=256) :: arg0, arg1, arg2
  integer :: n, repeat, iter
  integer(int8), allocatable :: checksum(:)
  real(real64) :: start_time, end_time
  logical :: complex_float_check, complex_double_check
  integer(int64), parameter :: lcg_a = 2806196910506780709_int64
  integer(int64), parameter :: lcg_c = 1_int64
  integer(int64), parameter :: lcg_mask = huge(0_int64)
  real(real64), parameter :: lcg_scale = 1.0_real64 / 9223372036854775808.0_real64

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  read(arg1, *) n
  read(arg2, *) repeat
  if (n <= 0 .or. repeat <= 0) stop 1

  allocate(checksum(n))

  !$omp target data map(alloc: checksum(1:n))
  call complex_float(checksum, n)
  call complex_double(checksum, n)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call complex_float(checksum, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time (float) ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  !$omp target update from(checksum(1:n))
  complex_float_check = all(checksum == 5_int8)

  start_time = omp_get_wtime()
  do iter = 1, repeat
    call complex_double(checksum, n)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average kernel execution time (double) ', &
    (end_time - start_time) / real(repeat, real64), ' (s)'

  !$omp target update from(checksum(1:n))
  complex_double_check = all(checksum == 5_int8)
  !$omp end target data

  if (complex_float_check .and. complex_double_check) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(checksum)

contains

  subroutine complex_float(checksum, n)
    integer(int8), intent(out) :: checksum(n)
    integer, intent(in) :: n
    integer :: i
    integer(int8) :: s
    real(real32) :: r1, r2, r3, r4
    type(FloatComplex) :: z1, z2

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, n
      call random_float4(i - 1, r1, r2, r3, r4)
      z1 = make_FloatComplex(r1, r2)
      z2 = make_FloatComplex(r3, r4)

      s = 0_int8
      if (abs(Cabsf(Cmulf(z1, z2)) - Cabsf(z1) * Cabsf(z2)) < 1.0e-3_real32) s = s + 1_int8
      if (abs(Cabsf(Caddf(z1, z2)) * Cabsf(Caddf(z1, z2)) - &
          Crealf(Cmulf(Caddf(z1, z2), Caddf(Conjf(z1), Conjf(z2))))) < 1.0e-3_real32) s = s + 1_int8
      if (abs(Cabsf(Csubf(z1, z2)) * Cabsf(Csubf(z1, z2)) - &
          Crealf(Cmulf(Csubf(z1, z2), Csubf(Conjf(z1), Conjf(z2))))) < 1.0e-3_real32) s = s + 1_int8
      if (abs(Crealf(Caddf(Cmulf(z1, Conjf(z2)), Cmulf(z2, Conjf(z1)))) - &
          2.0_real32 * (Crealf(z1) * Crealf(z2) + Cimagf(z1) * Cimagf(z2))) < 1.0e-3_real32) &
          s = s + 1_int8
      if (abs(Cabsf(Cdivf(Conjf(z1), z2)) - Cabsf(Cdivf(Conjf(z1), Conjf(z2)))) < 1.0e-3_real32) &
          s = s + 1_int8

      checksum(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine complex_float

  subroutine complex_double(checksum, n)
    integer(int8), intent(out) :: checksum(n)
    integer, intent(in) :: n
    integer :: i
    integer(int8) :: s
    real(real64) :: r1, r2, r3, r4
    type(DoubleComplex) :: z1, z2

    !$omp target teams distribute parallel do thread_limit(256)
    do i = 1, n
      call random_double4(i - 1, r1, r2, r3, r4)
      z1 = make_DoubleComplex(r1, r2)
      z2 = make_DoubleComplex(r3, r4)

      s = 0_int8
      if (abs(Cabs(Cmul(z1, z2)) - Cabs(z1) * Cabs(z2)) < 1.0e-3_real64) s = s + 1_int8
      if (abs(Cabs(Cadd(z1, z2)) * Cabs(Cadd(z1, z2)) - &
          Creal(Cmul(Cadd(z1, z2), Cadd(Conj(z1), Conj(z2))))) < 1.0e-3_real64) s = s + 1_int8
      if (abs(Cabs(Csub(z1, z2)) * Cabs(Csub(z1, z2)) - &
          Creal(Cmul(Csub(z1, z2), Csub(Conj(z1), Conj(z2))))) < 1.0e-3_real64) s = s + 1_int8
      if (abs(Creal(Cadd(Cmul(z1, Conj(z2)), Cmul(z2, Conj(z1)))) - &
          2.0_real64 * (Creal(z1) * Creal(z2) + Cimag(z1) * Cimag(z2))) < 1.0e-3_real64) &
          s = s + 1_int8
      if (abs(Cabs(Cdiv(Conj(z1), z2)) - Cabs(Cdiv(Conj(z1), Conj(z2)))) < 1.0e-3_real64) &
          s = s + 1_int8

      checksum(i) = s
    end do
    !$omp end target teams distribute parallel do
  end subroutine complex_double

  real(real32) function Crealf(x)
    type(FloatComplex), intent(in) :: x
    Crealf = x%x
  end function Crealf

  real(real32) function Cimagf(x)
    type(FloatComplex), intent(in) :: x
    Cimagf = x%y
  end function Cimagf

  type(FloatComplex) function make_FloatComplex(r, i)
    real(real32), intent(in) :: r, i
    make_FloatComplex%x = r
    make_FloatComplex%y = i
  end function make_FloatComplex

  type(FloatComplex) function Conjf(x)
    type(FloatComplex), intent(in) :: x
    Conjf = make_FloatComplex(Crealf(x), -Cimagf(x))
  end function Conjf

  type(FloatComplex) function Caddf(x, y)
    type(FloatComplex), intent(in) :: x, y
    Caddf = make_FloatComplex(Crealf(x) + Crealf(y), Cimagf(x) + Cimagf(y))
  end function Caddf

  type(FloatComplex) function Csubf(x, y)
    type(FloatComplex), intent(in) :: x, y
    Csubf = make_FloatComplex(Crealf(x) - Crealf(y), Cimagf(x) - Cimagf(y))
  end function Csubf

  type(FloatComplex) function Cmulf(x, y)
    type(FloatComplex), intent(in) :: x, y
    Cmulf = make_FloatComplex((Crealf(x) * Crealf(y)) - (Cimagf(x) * Cimagf(y)), &
                              (Crealf(x) * Cimagf(y)) + (Cimagf(x) * Crealf(y)))
  end function Cmulf

  type(FloatComplex) function Cdivf(x, y)
    type(FloatComplex), intent(in) :: x, y
    real(real32) :: s, oos, ars, ais, brs, bis

    s = abs(Crealf(y)) + abs(Cimagf(y))
    oos = 1.0_real32 / s
    ars = Crealf(x) * oos
    ais = Cimagf(x) * oos
    brs = Crealf(y) * oos
    bis = Cimagf(y) * oos
    s = (brs * brs) + (bis * bis)
    oos = 1.0_real32 / s
    Cdivf = make_FloatComplex(((ars * brs) + (ais * bis)) * oos, &
                              ((ais * brs) - (ars * bis)) * oos)
  end function Cdivf

  real(real32) function Cabsf(x)
    type(FloatComplex), intent(in) :: x
    real(real32) :: a, b, v, w, t

    a = abs(Crealf(x))
    b = abs(Cimagf(x))
    if (a > b) then
      v = a
      w = b
    else
      v = b
      w = a
    end if
    t = w / v
    t = 1.0_real32 + t * t
    t = v * sqrt(t)
    if (v == 0.0_real32 .or. v > huge(0.0_real32) .or. w > huge(0.0_real32)) t = v + w
    Cabsf = t
  end function Cabsf

  real(real64) function Creal(x)
    type(DoubleComplex), intent(in) :: x
    Creal = x%x
  end function Creal

  real(real64) function Cimag(x)
    type(DoubleComplex), intent(in) :: x
    Cimag = x%y
  end function Cimag

  type(DoubleComplex) function make_DoubleComplex(r, i)
    real(real64), intent(in) :: r, i
    make_DoubleComplex%x = r
    make_DoubleComplex%y = i
  end function make_DoubleComplex

  type(DoubleComplex) function Conj(x)
    type(DoubleComplex), intent(in) :: x
    Conj = make_DoubleComplex(Creal(x), -Cimag(x))
  end function Conj

  type(DoubleComplex) function Cadd(x, y)
    type(DoubleComplex), intent(in) :: x, y
    Cadd = make_DoubleComplex(Creal(x) + Creal(y), Cimag(x) + Cimag(y))
  end function Cadd

  type(DoubleComplex) function Csub(x, y)
    type(DoubleComplex), intent(in) :: x, y
    Csub = make_DoubleComplex(Creal(x) - Creal(y), Cimag(x) - Cimag(y))
  end function Csub

  type(DoubleComplex) function Cmul(x, y)
    type(DoubleComplex), intent(in) :: x, y
    Cmul = make_DoubleComplex((Creal(x) * Creal(y)) - (Cimag(x) * Cimag(y)), &
                              (Creal(x) * Cimag(y)) + (Cimag(x) * Creal(y)))
  end function Cmul

  type(DoubleComplex) function Cdiv(x, y)
    type(DoubleComplex), intent(in) :: x, y
    real(real64) :: s, oos, ars, ais, brs, bis

    s = abs(Creal(y)) + abs(Cimag(y))
    oos = 1.0_real64 / s
    ars = Creal(x) * oos
    ais = Cimag(x) * oos
    brs = Creal(y) * oos
    bis = Cimag(y) * oos
    s = (brs * brs) + (bis * bis)
    oos = 1.0_real64 / s
    Cdiv = make_DoubleComplex(((ars * brs) + (ais * bis)) * oos, &
                              ((ais * brs) - (ars * bis)) * oos)
  end function Cdiv

  real(real64) function Cabs(x)
    type(DoubleComplex), intent(in) :: x
    real(real64) :: a, b, v, w, t

    a = abs(Creal(x))
    b = abs(Cimag(x))
    if (a > b) then
      v = a
      w = b
    else
      v = b
      w = a
    end if
    t = w / v
    t = 1.0_real64 + t * t
    t = v * sqrt(t)
    if (v == 0.0_real64 .or. v > huge(0.0_real64) .or. w > huge(0.0_real64)) t = v + w
    Cabs = t
  end function Cabs

  subroutine random_float4(i0, r1, r2, r3, r4)
    integer, intent(in) :: i0
    real(real32), intent(out) :: r1, r2, r3, r4
    integer(int64) :: seed

    seed = fast_forward_lcg(1_int64, int(i0, int64))
    r1 = real(lcg_random_double(seed), real32)
    r2 = real(lcg_random_double(seed), real32)
    r3 = real(lcg_random_double(seed), real32)
    r4 = real(lcg_random_double(seed), real32)
  end subroutine random_float4

  subroutine random_double4(i0, r1, r2, r3, r4)
    integer, intent(in) :: i0
    real(real64), intent(out) :: r1, r2, r3, r4
    integer(int64) :: seed

    seed = fast_forward_lcg(1_int64, int(i0, int64))
    r1 = lcg_random_double(seed)
    r2 = lcg_random_double(seed)
    r3 = lcg_random_double(seed)
    r4 = lcg_random_double(seed)
  end subroutine random_double4

  real(real64) function lcg_random_double(seed)
    integer(int64), intent(inout) :: seed

    seed = iand(lcg_a * seed + lcg_c, lcg_mask)
    lcg_random_double = real(seed, real64) * lcg_scale
  end function lcg_random_double

  integer(int64) function fast_forward_lcg(seed, n)
    integer(int64), intent(in) :: seed, n
    integer(int64) :: a, c, a_new, c_new, n_work

    a = lcg_a
    c = lcg_c
    a_new = 1_int64
    c_new = 0_int64
    n_work = iand(n, lcg_mask)

    do while (n_work > 0_int64)
      if (iand(n_work, 1_int64) /= 0_int64) then
        a_new = iand(a_new * a, lcg_mask)
        c_new = iand(c_new * a + c, lcg_mask)
      end if

      c = iand(c * (a + 1_int64), lcg_mask)
      a = iand(a * a, lcg_mask)
      n_work = shiftr(n_work, 1)
    end do

    fast_forward_lcg = iand(a_new * seed + c_new, lcg_mask)
  end function fast_forward_lcg

end program main
