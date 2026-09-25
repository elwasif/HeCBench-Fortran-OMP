! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int16, int32, real32, real64
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

  character(len=256) :: arg
  integer :: n, repeat
  real(real32), allocatable :: x(:), y(:), hf(:), rf(:)
  integer(int32), allocatable :: hi(:), ri(:)
  integer(int16), allocatable :: hs(:), rs(:)
  real(real64) :: start_time, end_time, error
  integer :: iter

  if (command_argument_count() /= 2) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <number of coordinates> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) n
  call get_command_argument(2, arg)
  read(arg, *) repeat
  if (n <= 0 .or. repeat <= 0) stop 1

  allocate(x(0:n - 1), y(0:n - 1), hf(0:n - 1), rf(0:n - 1))
  allocate(hi(0:n - 1), ri(0:n - 1), hs(0:n - 1), rs(0:n - 1))

  call initialize_inputs(n, x, y)

  !$omp target data map(to: x(0:n - 1), y(0:n - 1)) map(alloc: hf(0:n - 1), hi(0:n - 1), hs(0:n - 1))
  write(*,'(A)') ''
  write(*,'(A)') '======== output type is f32 ========'
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call compute_f(n, y, x, hf)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time: ', (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
  !$omp target update from(hf(0:n - 1))

  call reference_f(n, y, x, rf)
  error = rmse_f(n, rf, hf)
  write(*,'(A,F0.6)') 'RMSE: ', error

  write(*,'(A)') ''
  write(*,'(A)') '======== output type is i32 ========'
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call compute_i(n, y, x, hi)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time: ', (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
  !$omp target update from(hi(0:n - 1))

  call reference_i(n, y, x, ri)
  error = rmse_i(n, ri, hi)
  write(*,'(A,F0.6)') 'RMSE: ', error

  write(*,'(A)') ''
  write(*,'(A)') '======== output type is i16 ========'
  start_time = omp_get_wtime()
  do iter = 1, repeat
    call compute_s(n, y, x, hs)
  end do
  end_time = omp_get_wtime()
  write(*,'(A,F0.6,A)') 'Average execution time: ', (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
  !$omp target update from(hs(0:n - 1))

  call reference_s(n, y, x, rs)
  error = rmse_s(n, rs, hs)
  write(*,'(A,F0.6)') 'RMSE: ', error
  !$omp end target data

  deallocate(rs, hs, ri, hi, rf, hf, y, x)

contains

  subroutine initialize_inputs(n, x, y)
    integer, intent(in) :: n
    real(real32), intent(out) :: x(0:), y(0:)
    integer :: i

    call c_srand(123_c_int)
    do i = 0, n - 1
      x(i) = real(c_rand(), real32) / 2147483647.0_real32 + 1.57_real32
      y(i) = real(c_rand(), real32) / 2147483647.0_real32 + 1.57_real32
    end do
  end subroutine initialize_inputs

  subroutine compute_f(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    real(real32), intent(out) :: r(0:)
    integer :: i
    real(real32) :: vy, vx

    !$omp target teams distribute parallel do thread_limit(256) private(i, vy, vx)
    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      r(i) = safe_atan2f(3, vy, vx) + safe_atan2f(5, vy, vx) + safe_atan2f(7, vy, vx) + &
        safe_atan2f(9, vy, vx) + safe_atan2f(11, vy, vx) + safe_atan2f(13, vy, vx) + safe_atan2f(15, vy, vx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_f

  subroutine compute_i(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    integer(int32), intent(out) :: r(0:)
    integer :: i
    real(real32) :: vy, vx

    !$omp target teams distribute parallel do thread_limit(256) private(i, vy, vx)
    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      r(i) = unsafe_atan2i(3, vy, vx) + unsafe_atan2i(5, vy, vx) + unsafe_atan2i(7, vy, vx) + &
        unsafe_atan2i(9, vy, vx) + unsafe_atan2i(11, vy, vx) + unsafe_atan2i(13, vy, vx) + unsafe_atan2i(15, vy, vx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_i

  subroutine compute_s(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    integer(int16), intent(out) :: r(0:)
    integer :: i
    integer(int32) :: value
    real(real32) :: vy, vx

    !$omp target teams distribute parallel do thread_limit(256) private(i, value, vy, vx)
    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      value = int(unsafe_atan2s(3, vy, vx), int32) + int(unsafe_atan2s(5, vy, vx), int32) + &
        int(unsafe_atan2s(7, vy, vx), int32) + int(unsafe_atan2s(9, vy, vx), int32)
      r(i) = int(value, int16)
    end do
    !$omp end target teams distribute parallel do
  end subroutine compute_s

  subroutine reference_f(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    real(real32), intent(out) :: r(0:)
    integer :: i
    real(real32) :: vy, vx

    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      r(i) = safe_atan2f(3, vy, vx) + safe_atan2f(5, vy, vx) + safe_atan2f(7, vy, vx) + &
        safe_atan2f(9, vy, vx) + safe_atan2f(11, vy, vx) + safe_atan2f(13, vy, vx) + safe_atan2f(15, vy, vx)
    end do
  end subroutine reference_f

  subroutine reference_i(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    integer(int32), intent(out) :: r(0:)
    integer :: i
    real(real32) :: vy, vx

    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      r(i) = unsafe_atan2i(3, vy, vx) + unsafe_atan2i(5, vy, vx) + unsafe_atan2i(7, vy, vx) + &
        unsafe_atan2i(9, vy, vx) + unsafe_atan2i(11, vy, vx) + unsafe_atan2i(13, vy, vx) + unsafe_atan2i(15, vy, vx)
    end do
  end subroutine reference_i

  subroutine reference_s(n, x, y, r)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    integer(int16), intent(out) :: r(0:)
    integer :: i
    integer(int32) :: value
    real(real32) :: vy, vx

    do i = 0, n - 1
      vy = y(i)
      vx = x(i)
      value = int(unsafe_atan2s(3, vy, vx), int32) + int(unsafe_atan2s(5, vy, vx), int32) + &
        int(unsafe_atan2s(7, vy, vx), int32) + int(unsafe_atan2s(9, vy, vx), int32)
      r(i) = int(value, int16)
    end do
  end subroutine reference_s

  real(real64) function rmse_f(n, ref, got) result(value)
    integer, intent(in) :: n
    real(real32), intent(in) :: ref(0:), got(0:)
    integer :: i
    real(real64) :: error

    error = 0.0_real64
    do i = 0, n - 1
      if (abs(ref(i) - got(i)) > 1.0e-3_real32) then
        error = error + real(ref(i) - got(i), real64) * real(ref(i) - got(i), real64)
      end if
    end do
    value = sqrt(error / real(n, real64))
  end function rmse_f

  real(real64) function rmse_i(n, ref, got) result(value)
    integer, intent(in) :: n
    integer(int32), intent(in) :: ref(0:), got(0:)
    integer :: i
    real(real64) :: error, diff

    error = 0.0_real64
    do i = 0, n - 1
      if (abs(ref(i) - got(i)) > 0) then
        diff = real(ref(i) - got(i), real64)
        error = error + diff * diff
      end if
    end do
    value = sqrt(error / real(n, real64))
  end function rmse_i

  real(real64) function rmse_s(n, ref, got) result(value)
    integer, intent(in) :: n
    integer(int16), intent(in) :: ref(0:), got(0:)
    integer :: i
    real(real64) :: error, diff

    error = 0.0_real64
    do i = 0, n - 1
      if (abs(int(ref(i), int32) - int(got(i), int32)) > 0) then
        diff = real(int(ref(i), int32) - int(got(i), int32), real64)
        error = error + diff * diff
      end if
    end do
    value = sqrt(error / real(n, real64))
  end function rmse_s

  real(real32) function safe_atan2f(degree, y, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: y, x
    real(real32) :: use_x

    use_x = x
    if (y == 0.0_real32 .and. x == 0.0_real32) use_x = 0.2_real32
    value = unsafe_atan2f_impl(degree, y, use_x)
  end function safe_atan2f

  real(real32) function unsafe_atan2f_impl(degree, y, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: y, x
    real(real32), parameter :: pi4f = 0.7853981633974483_real32
    real(real32), parameter :: pi34f = 2.356194490192345_real32
    real(real32) :: r, angle

    r = (abs(x) - abs(y)) / (abs(x) + abs(y))
    if (x < 0.0_real32) r = -r
    if (x >= 0.0_real32) then
      angle = pi4f
    else
      angle = pi34f
    end if
    angle = angle + approx_atan2f_p(degree, r)
    if (y < 0.0_real32) then
      value = -angle
    else
      value = angle
    end if
  end function unsafe_atan2f_impl

  integer(int32) function unsafe_atan2i(degree, y, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: y, x
    integer(int32), parameter :: pi4 = 536870912_int32
    integer(int32), parameter :: pi34 = 1610612736_int32
    real(real32) :: r
    integer(int32) :: angle

    r = (abs(x) - abs(y)) / (abs(x) + abs(y))
    if (x < 0.0_real32) r = -r
    if (x >= 0.0_real32) then
      angle = pi4
    else
      angle = pi34
    end if
    angle = angle + int(approx_atan2i_p(degree, r), int32)
    if (y < 0.0_real32) then
      value = -angle
    else
      value = angle
    end if
  end function unsafe_atan2i

  integer(int16) function unsafe_atan2s(degree, y, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: y, x
    integer(int32), parameter :: pi4 = 8192_int32
    integer(int32), parameter :: pi34 = 24576_int32
    real(real32) :: r
    integer(int32) :: angle

    r = (abs(x) - abs(y)) / (abs(x) + abs(y))
    if (x < 0.0_real32) r = -r
    if (x >= 0.0_real32) then
      angle = pi4
    else
      angle = pi34
    end if
    angle = angle + int(approx_atan2s_p(degree, r), int32)
    if (y < 0.0_real32) angle = -angle
    value = int(angle, int16)
  end function unsafe_atan2s

  real(real32) function approx_atan2f_p(degree, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: x
    real(real32) :: z

    z = x * x
    select case (degree)
    case (3)
      value = x * (-0.9723931551_real32 + z * 0.1915416718_real32)
    case (5)
      value = x * (-0.9951133728_real32 + z * (0.2886903882_real32 + z * (-0.07933916_real32)))
    case (7)
      value = x * (-0.9992147684_real32 + z * (0.3211759329_real32 + z * (-0.1462625563_real32 + z * 0.0390017629_real32)))
    case (9)
      value = x * (-0.9998816252_real32 + z * (0.3302908242_real32 + z * (-0.1801744998_real32 + &
        z * (0.0851577073_real32 + z * (-0.0208360124_real32)))))
    case (11)
      value = x * (-0.9999772310_real32 + z * (0.3326220512_real32 + z * (-0.1930824816_real32 + &
        z * (0.1164280772_real32 + z * (-0.0526488191_real32 + z * 0.0117195472_real32)))))
    case (13)
      value = x * (-0.9999960661_real32 + z * (0.3331727982_real32 + z * (-0.1985587776_real32 + &
        z * (0.1326826215_real32 + z * (-0.0795908570_real32 + z * (0.0337525904_real32 + &
        z * (-0.0068054050_real32)))))))
    case default
      value = x * (-0.9999885559_real32 + z * (0.3330533803_real32 + z * (-0.1990792751_real32 + &
        z * (0.1382764578_real32 + z * (-0.0963287801_real32 + z * (0.0558212623_real32 + &
        z * (-0.0210717320_real32 + z * 0.0040353443_real32)))))))
    end select
  end function approx_atan2f_p

  real(real32) function approx_atan2i_p(degree, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: x
    real(real32) :: z

    z = x * x
    select case (degree)
    case (3)
      value = x * (-664694912.0_real32 + z * 131209024.0_real32)
    case (5)
      value = x * (-680392064.0_real32 + z * (197338400.0_real32 + z * (-54233256.0_real32)))
    case (7)
      value = x * (-683027840.0_real32 + z * (219543904.0_real32 + z * (-99981040.0_real32 + z * 26649684.0_real32)))
    case (9)
      value = x * (-683473920.0_real32 + z * (225785056.0_real32 + z * (-123151184.0_real32 + &
        z * (58210592.0_real32 + z * (-14249276.0_real32)))))
    case (11)
      value = x * (-683549696.0_real32 + z * (227369312.0_real32 + z * (-132297008.0_real32 + &
        z * (79584144.0_real32 + z * (-35987016.0_real32 + z * 8010488.0_real32)))))
    case default
      value = x * (-683562624.0_real32 + z * (227746080.0_real32 + z * (-135400128.0_real32 + &
        z * (90460848.0_real32 + z * (-54431464.0_real32 + z * (22973256.0_real32 + z * (-4657049.0_real32)))))))
    end select
  end function approx_atan2i_p

  real(real32) function approx_atan2s_p(degree, x) result(value)
    integer, intent(in) :: degree
    real(real32), intent(in) :: x
    real(real32) :: z

    z = x * x
    select case (degree)
    case (3)
      value = x * (-10142.439453125_real32 + z * 2002.0908203125_real32)
    case (5)
      value = x * (-10381.9609375_real32 + z * (3011.1513671875_real32 + z * (-827.538330078125_real32)))
    case (7)
      value = x * (-10422.177734375_real32 + z * (3349.97412109375_real32 + &
        z * (-1525.589599609375_real32 + z * 406.64190673828125_real32)))
    case default
      value = x * (-10428.984375_real32 + z * (3445.20654296875_real32 + z * (-1879.137939453125_real32 + &
        z * (888.22314453125_real32 + z * (-217.42669677734375_real32)))))
    end select
  end function approx_atan2s_p

end program main
