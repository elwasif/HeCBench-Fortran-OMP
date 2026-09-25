! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int8, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256
  integer(c_int), parameter :: c_rand_max = 2147483647_c_int
  real(real64), parameter :: two_pi = 6.283185307179586_real64
  real(real32), parameter :: code(256) = [ &
    -0.992968738079071_real32, -0.9789062738418579_real32, -0.96484375_real32, -0.9507812261581421_real32,  &
    -0.936718761920929_real32, -0.922656238079071_real32, -0.9085937738418579_real32, -0.89453125_real32,  &
    -0.8804687261581421_real32, -0.866406261920929_real32, -0.852343738079071_real32, -0.8382812738418579_real32,  &
    -0.82421875_real32, -0.8101562261581421_real32, -0.796093761920929_real32, -0.782031238079071_real32,  &
    -0.7679687738418579_real32, -0.75390625_real32, -0.7398437261581421_real32, -0.725781261920929_real32,  &
    -0.7117187976837158_real32, -0.6976562738418579_real32, -0.68359375_real32, -0.6695312261581421_real32,  &
    -0.655468761920929_real32, -0.6414062976837158_real32, -0.6273437738418579_real32, -0.61328125_real32,  &
    -0.5992187261581421_real32, -0.585156261920929_real32, -0.5710937976837158_real32, -0.5570312738418579_real32,  &
    -0.54296875_real32, -0.5289062261581421_real32, -0.5148437023162842_real32, -0.500781238079071_real32,  &
    -0.48671871423721313_real32, -0.47265625_real32, -0.4585937261581421_real32, -0.44453126192092896_real32,  &
    -0.43046873807907104_real32, -0.4164062440395355_real32, -0.40234375_real32, -0.3882812261581421_real32,  &
    -0.37421876192092896_real32, -0.36015623807907104_real32, -0.3460937440395355_real32, -0.33203125_real32,  &
    -0.3179687261581421_real32, -0.30390626192092896_real32, -0.28984373807907104_real32, -0.2757812738418579_real32,  &
    -0.26171875_real32, -0.24765624105930328_real32, -0.23359374701976776_real32, -0.21953125298023224_real32,  &
    -0.20546874403953552_real32, -0.19140625_real32, -0.17734375596046448_real32, -0.16328124701976776_real32,  &
    -0.14921875298023224_real32, -0.13515624403953552_real32, -0.12109375_real32, -0.10703125596046448_real32,  &
    -0.09859374910593033_real32, -0.09578125923871994_real32, -0.09296875447034836_real32, -0.09015624970197678_real32,  &
    -0.08734375238418579_real32, -0.08453124761581421_real32, -0.08171875774860382_real32, -0.07890625298023224_real32,  &
    -0.07609374821186066_real32, -0.07328125089406967_real32, -0.07046874612569809_real32, -0.0676562562584877_real32,  &
    -0.06484375149011612_real32, -0.062031250447034836_real32, -0.05921875312924385_real32, -0.05640624836087227_real32,  &
    -0.053593751043081284_real32, -0.05078125_real32, -0.047968748956918716_real32, -0.04515625163912773_real32,  &
    -0.04234374687075615_real32, -0.039531249552965164_real32, -0.03671875223517418_real32, -0.033906251192092896_real32,  &
    -0.031093750149011612_real32, -0.028281250968575478_real32, -0.025468749925494194_real32, -0.02265625074505806_real32,  &
    -0.019843751564621925_real32, -0.017031250521540642_real32, -0.014218750409781933_real32, -0.011406250298023224_real32,  &
    -0.009718749672174454_real32, -0.009156249463558197_real32, -0.008593750186264515_real32, -0.008031249977648258_real32,  &
    -0.0074687497690320015_real32, -0.006906250026077032_real32, -0.006343749817460775_real32, -0.005781250074505806_real32,  &
    -0.0052187503315508366_real32, -0.0046562496572732925_real32, -0.004093749914318323_real32, -0.0035312497057020664_real32,  &
    -0.002968749962747097_real32, -0.002406249986961484_real32, -0.001843750011175871_real32, -0.001281249918974936_real32,  &
    -0.0009437500848434865_real32, -0.0008312499849125743_real32, -0.0007187500596046448_real32, -0.0006062500760890543_real32,  &
    -0.000493750034365803_real32, -0.0003812500217463821_real32, -0.0002687500382307917_real32, -0.00015625001105945557_real32,  &
    -8.874999912222847e-05_real32, -6.625000241911039e-05_real32, -4.374999844003469e-05_real32, -2.1249998098937795e-05_real32,  &
    -7.749999895168003e-06_real32, -3.250000190746505e-06_real32, -5.500000384017767e-07_real32, 0.0_real32,  &
    5.500000384017767e-07_real32, 3.250000190746505e-06_real32, 7.749999895168003e-06_real32, 2.1249998098937795e-05_real32,  &
    4.374999844003469e-05_real32, 6.625000241911039e-05_real32, 8.874999912222847e-05_real32, 0.00015625001105945557_real32,  &
    0.0002687500382307917_real32, 0.0003812500217463821_real32, 0.000493750034365803_real32, 0.0006062500760890543_real32,  &
    0.0007187500596046448_real32, 0.0008312499849125743_real32, 0.0009437500848434865_real32, 0.001281249918974936_real32,  &
    0.001843750011175871_real32, 0.002406249986961484_real32, 0.002968749962747097_real32, 0.0035312497057020664_real32,  &
    0.004093749914318323_real32, 0.0046562496572732925_real32, 0.0052187503315508366_real32, 0.005781250074505806_real32,  &
    0.006343749817460775_real32, 0.006906250026077032_real32, 0.0074687497690320015_real32, 0.008031249977648258_real32,  &
    0.008593750186264515_real32, 0.009156249463558197_real32, 0.009718749672174454_real32, 0.011406250298023224_real32,  &
    0.014218750409781933_real32, 0.017031250521540642_real32, 0.019843751564621925_real32, 0.02265625074505806_real32,  &
    0.025468749925494194_real32, 0.028281250968575478_real32, 0.031093750149011612_real32, 0.033906251192092896_real32,  &
    0.03671875223517418_real32, 0.039531249552965164_real32, 0.04234374687075615_real32, 0.04515625163912773_real32,  &
    0.047968748956918716_real32, 0.05078125_real32, 0.053593751043081284_real32, 0.05640624836087227_real32,  &
    0.05921875312924385_real32, 0.062031250447034836_real32, 0.06484375149011612_real32, 0.0676562562584877_real32,  &
    0.07046874612569809_real32, 0.07328125089406967_real32, 0.07609374821186066_real32, 0.07890625298023224_real32,  &
    0.08171875774860382_real32, 0.08453124761581421_real32, 0.08734375238418579_real32, 0.09015624970197678_real32,  &
    0.09296875447034836_real32, 0.09578125923871994_real32, 0.09859374910593033_real32, 0.10703125596046448_real32,  &
    0.12109375_real32, 0.13515624403953552_real32, 0.14921875298023224_real32, 0.16328124701976776_real32,  &
    0.17734375596046448_real32, 0.19140625_real32, 0.20546874403953552_real32, 0.21953125298023224_real32,  &
    0.23359374701976776_real32, 0.24765624105930328_real32, 0.26171875_real32, 0.2757812738418579_real32,  &
    0.28984373807907104_real32, 0.30390626192092896_real32, 0.3179687261581421_real32, 0.33203125_real32,  &
    0.3460937440395355_real32, 0.36015623807907104_real32, 0.37421876192092896_real32, 0.3882812261581421_real32,  &
    0.40234375_real32, 0.4164062440395355_real32, 0.43046873807907104_real32, 0.44453126192092896_real32,  &
    0.4585937261581421_real32, 0.47265625_real32, 0.48671871423721313_real32, 0.500781238079071_real32,  &
    0.5148437023162842_real32, 0.5289062261581421_real32, 0.54296875_real32, 0.5570312738418579_real32,  &
    0.5710937976837158_real32, 0.585156261920929_real32, 0.5992187261581421_real32, 0.61328125_real32,  &
    0.6273437738418579_real32, 0.6414062976837158_real32, 0.655468761920929_real32, 0.6695312261581421_real32,  &
    0.68359375_real32, 0.6976562738418579_real32, 0.7117187976837158_real32, 0.725781261920929_real32,  &
    0.7398437261581421_real32, 0.75390625_real32, 0.7679687738418579_real32, 0.782031238079071_real32,  &
    0.796093761920929_real32, 0.8101562261581421_real32, 0.82421875_real32, 0.8382812738418579_real32,  &
    0.852343738079071_real32, 0.866406261920929_real32, 0.8804687261581421_real32, 0.89453125_real32,  &
    0.9085937738418579_real32, 0.922656238079071_real32, 0.936718761920929_real32, 0.9507812261581421_real32,  &
    0.96484375_real32, 0.9789062738418579_real32, 0.992968738079071_real32, 1.0_real32 &
  ]

  character(len=256) :: arg0, arg1, arg2
  integer(int64) :: n, i
  integer :: repeat, iter, grid
  real(real32), allocatable :: a(:)
  integer(int8), allocatable :: out(:), ref(:)
  real(real64) :: start_time, end_time

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

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 2) then
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg0)
    write(*,'(A)') ' <number of elements> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  n = parse_atol(arg1)
  repeat = parse_atoi(arg2)

  allocate(a(n), out(n), ref(n))
  call c_srand(19937_c_int)
  do i = 1, n
    a(i) = c_normal_input()
    ref(i) = uint8_byte(dquantize(code, 0.0_real32, a(i)))
  end do

  grid = int((n + int(block_size - 1, int64)) / int(block_size, int64))
  !$omp target data map(to: a(1:n), code(1:256)) map(from: out(1:n))
  start_time = omp_get_wtime()
  do iter = 1, repeat
    !$omp target teams distribute parallel do num_teams(grid) num_threads(block_size)
    do i = 1, n
      out(i) = uint8_byte(dquantize(code, 0.0_real32, a(i)))
    end do
    !$omp end target teams distribute parallel do
  end do
  end_time = omp_get_wtime()
  write(*,'(A,I0,A,F0.6,A)') 'Average execution time of kQuantize kernel with block size ', &
    block_size, ': ', (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
  !$omp end target data

  if (all(out == ref)) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(a, out, ref)

contains

  real(real32) function c_normal_input()
    real(real64) :: u1, u2

    u1 = max(c_rand_unit(), tiny(1.0_real64))
    u2 = c_rand_unit()
    c_normal_input = real(sqrt(-2.0_real64 * log(u1)) * cos(two_pi * u2), real32)
  end function c_normal_input

  real(real64) function c_rand_unit()
    c_rand_unit = (real(c_rand(), real64) + 0.5_real64) / (real(c_rand_max, real64) + 1.0_real64)
  end function c_rand_unit

  integer(int64) function parse_atol(arg) result(value)
    character(len=*), intent(in) :: arg
    integer :: pos, sign, digit

    value = 0_int64
    pos = first_nonblank(arg)
    sign = 1
    if (pos <= len_trim(arg)) then
      if (arg(pos:pos) == '-') then
        sign = -1
        pos = pos + 1
      else if (arg(pos:pos) == '+') then
        pos = pos + 1
      end if
    end if

    do while (pos <= len_trim(arg))
      digit = iachar(arg(pos:pos)) - iachar('0')
      if (digit < 0 .or. digit > 9) exit
      value = value * 10_int64 + int(digit, int64)
      pos = pos + 1
    end do
    value = value * int(sign, int64)
  end function parse_atol

  integer function parse_atoi(arg) result(value)
    character(len=*), intent(in) :: arg
    value = int(parse_atol(arg))
  end function parse_atoi

  integer function first_nonblank(arg) result(pos)
    character(len=*), intent(in) :: arg
    do pos = 1, len_trim(arg)
      if (arg(pos:pos) /= ' ' .and. arg(pos:pos) /= achar(9)) return
    end do
    pos = len_trim(arg) + 1
  end function first_nonblank

  integer(int8) function uint8_byte(value) result(byte_value)
    integer, intent(in) :: value
    if (value > 127) then
      byte_value = int(value - 256, int8)
    else
      byte_value = int(value, int8)
    end if
  end function uint8_byte

  integer function dquantize(smem_code, rand, x)
    real(real32), intent(in) :: smem_code(:)
    real(real32), intent(in) :: rand, x
    integer :: pivot, upper_pivot, lower_pivot, step
    real(real32) :: lower, upper, val, midpoint

    pivot = 127
    upper_pivot = 255
    lower_pivot = 0
    lower = -1.0_real32
    upper = 1.0_real32
    val = smem_code(pivot + 1)
    step = 64
    do while (step > 0)
      if (x > val) then
        lower_pivot = pivot
        lower = val
        pivot = pivot + step
      else
        upper_pivot = pivot
        upper = val
        pivot = pivot - step
      end if
      val = smem_code(pivot + 1)
      step = step / 2
    end do

    if (upper_pivot == 255) upper = smem_code(upper_pivot + 1)
    if (lower_pivot == 0) lower = smem_code(lower_pivot + 1)

    if (x > val) then
      midpoint = (upper + val) * 0.5_real32
      if (x > midpoint) then
        dquantize = upper_pivot
      else
        dquantize = pivot
      end if
    else
      midpoint = (lower + val) * 0.5_real32
      if (x < midpoint) then
        dquantize = lower_pivot
      else
        dquantize = pivot
      end if
    end if
  end function dquantize

end program main
