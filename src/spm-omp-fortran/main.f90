! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_bool, c_int, c_signed_char
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: num_threads = 128
  integer, parameter :: num_blocks = 256
  integer, parameter :: hist_size = 65536
  integer, parameter :: ran_count = 97
  real(real32), parameter :: ran(ran_count) = [ &
    0.656619_real32,0.891183_real32,0.488144_real32,0.992646_real32,0.373326_real32, &
    0.531378_real32,0.181316_real32,0.501944_real32,0.422195_real32,0.660427_real32, &
    0.673653_real32,0.95733_real32,0.191866_real32,0.111216_real32,0.565054_real32, &
    0.969166_real32,0.0237439_real32,0.870216_real32,0.0268766_real32,0.519529_real32, &
    0.192291_real32,0.715689_real32,0.250673_real32,0.933865_real32,0.137189_real32, &
    0.521622_real32,0.895202_real32,0.942387_real32,0.335083_real32,0.437364_real32, &
    0.471156_real32,0.14931_real32,0.135864_real32,0.532498_real32,0.725789_real32, &
    0.398703_real32,0.358419_real32,0.285279_real32,0.868635_real32,0.626413_real32, &
    0.241172_real32,0.978082_real32,0.640501_real32,0.229849_real32,0.681335_real32, &
    0.665823_real32,0.134718_real32,0.0224933_real32,0.262199_real32,0.116515_real32, &
    0.0693182_real32,0.85293_real32,0.180331_real32,0.0324186_real32,0.733926_real32, &
    0.536517_real32,0.27603_real32,0.368458_real32,0.0128863_real32,0.889206_real32, &
    0.866021_real32,0.254247_real32,0.569481_real32,0.159265_real32,0.594364_real32, &
    0.3311_real32,0.658613_real32,0.863634_real32,0.567623_real32,0.980481_real32, &
    0.791832_real32,0.152594_real32,0.833027_real32,0.191863_real32,0.638987_real32, &
    0.669_real32,0.772088_real32,0.379818_real32,0.441585_real32,0.48306_real32, &
    0.608106_real32,0.175996_real32,0.00202556_real32,0.790224_real32,0.513609_real32, &
    0.213229_real32,0.10345_real32,0.157337_real32,0.407515_real32,0.407757_real32, &
    0.0526927_real32,0.941815_real32,0.149972_real32,0.384374_real32,0.311059_real32, &
    0.168534_real32,0.896648_real32]

  integer :: v, repeat, data_size, vol_size, i, count_device, count_host, max_diff
  integer(c_signed_char), allocatable :: f(:), g(:), ivf(:), ivg(:), ivf_ref(:), ivg_ref(:)
  logical(c_bool), allocatable :: threshold(:), threshold_ref(:)
  integer(int32), allocatable :: hist_device(:), hist_host(:)
  real(real32) :: matrix(16)
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

  if (command_argument_count() /= 2) then
    print '(A)', 'Usage: ./main <dimension> <repeat>'
    stop 1
  end if

  v = read_arg(1)
  repeat = read_arg(2)
  if (v <= 2 .or. repeat <= 0) error stop 'spm arguments must be dimension > 2 and repeat > 0'

  data_size = (v + 1) * (v + 1) * (v + 5)
  vol_size = v * v * v
  allocate(f(data_size), g(data_size), ivf(vol_size), ivg(vol_size), ivf_ref(vol_size), ivg_ref(vol_size))
  allocate(threshold(vol_size), threshold_ref(vol_size), hist_device(hist_size), hist_host(hist_size))

  call initialize_inputs(matrix, f, g)
  ivf = 0
  ivg = 0
  threshold = .false.

  !$omp target data map(to: matrix(1:16), g(1:data_size), f(1:data_size)) &
  !$omp& map(from: ivf(1:vol_size), ivg(1:vol_size), threshold(1:vol_size))
  start_time = omp_get_wtime()
  do i = 1, repeat
    call spm_device(matrix, vol_size, g, f, v, ivf, ivg, threshold)
  end do
  end_time = omp_get_wtime()
  !$omp end target data
  write(*, '(A,F0.6,A)') 'Average kernel execution time: ', ((end_time - start_time) * 1.0e3_real64) / real(repeat, real64), ' (ms)'

  hist_device = 0
  count_device = build_histogram(ivf, ivg, threshold, hist_device)
  print '(A,I0)', 'Device count: ', count_device

  call spm_host(matrix, vol_size, g, f, v, ivf_ref, ivg_ref, threshold_ref)
  hist_host = 0
  count_host = build_histogram(ivf_ref, ivg_ref, threshold_ref, hist_host)
  print '(A,I0)', 'Host count: ', count_host

  max_diff = maxval(abs(hist_host - hist_device))
  print '(A,I0)', 'Maximum difference ', max_diff

  deallocate(f, g, ivf, ivg, ivf_ref, ivg_ref, threshold, threshold_ref, hist_device, hist_host)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_inputs(matrix, f, g)
    real(real32), intent(out) :: matrix(:)
    integer(c_signed_char), intent(out) :: f(:), g(:)
    integer :: i
    call c_srand(123_c_int)
    do i = 1, size(matrix)
      matrix(i) = real(c_rand(), real32) / 2147483647.0_real32
    end do
    do i = 1, size(f)
      f(i) = int(mod(c_rand(), 256_c_int), c_signed_char)
      g(i) = int(mod(c_rand(), 256_c_int), c_signed_char)
    end do
  end subroutine initialize_inputs

  subroutine spm_device(matrix, data_size, g, f, v, ivf, ivg, threshold)
    real(real32), intent(in) :: matrix(:)
    integer, intent(in) :: data_size, v
    integer(c_signed_char), intent(in) :: g(:), f(:)
    integer(c_signed_char), intent(out) :: ivf(:), ivg(:)
    logical(c_bool), intent(out) :: threshold(:)
    integer :: idx, x_datasize, y_datasize, ix, iy, iz, base, upper_base
    integer :: k111, k112, k121, k122, k211, k212, k221, k222
    real(real32) :: xx_temp, yy_temp, zz_temp, rx, ry, rz, xp, yp, zp
    real(real32) :: dx1, dy1, dz1, dx2, dy2, dz2, vf, vg, r

    x_datasize = v - 2
    y_datasize = v - 2
    !$omp target teams distribute parallel do num_teams(num_blocks) thread_limit(num_threads) &
    !$omp& private(xx_temp, yy_temp, zz_temp, rx, ry, rz, xp, yp, zp, r, ix, iy, iz, base, upper_base) &
    !$omp& private(dx1, dy1, dz1, dx2, dy2, dz2, k111, k112, k121, k122, k211, k212, k221, k222, vf, vg)
    do idx = 1, data_size
      xx_temp = real(mod(idx - 1, x_datasize), real32) + 1.0_real32
      yy_temp = real(mod(int(floor(real(idx - 1, real32) / real(x_datasize, real32))), y_datasize), real32) + 1.0_real32
      zz_temp = floor(real(idx - 1, real32) / real(x_datasize, real32)) / real(y_datasize, real32) + 1.0_real32
      r = ran(mod(idx - 1, ran_count) + 1)
      rx = xx_temp + r
      ry = yy_temp + r
      rz = zz_temp + r
      xp = matrix(1) * rx + matrix(5) * ry + matrix(9) * rz + matrix(13)
      yp = matrix(2) * rx + matrix(6) * ry + matrix(10) * rz + matrix(14)
      zp = matrix(3) * rx + matrix(7) * ry + matrix(11) * rz + matrix(15)
      if (zp >= 1.0_real32 .and. zp < real(v, real32) .and. yp >= 1.0_real32 .and. yp < real(v, real32) &
          .and. xp >= 1.0_real32 .and. xp < real(v, real32)) then
        call interp_inline(f, v, xp, yp, zp, vf)
        call interp_inline(g, v, rx, ry, rz, vg)
        ivf(idx) = int(floor(vf + 0.5_real32), c_signed_char)
        ivg(idx) = int(floor(vg + 0.5_real32), c_signed_char)
        threshold(idx) = .true.
      else
        ivf(idx) = 0
        ivg(idx) = 0
        threshold(idx) = .false.
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine spm_device

  subroutine spm_host(matrix, data_size, g, f, v, ivf, ivg, threshold)
    real(real32), intent(in) :: matrix(:)
    integer, intent(in) :: data_size, v
    integer(c_signed_char), intent(in) :: g(:), f(:)
    integer(c_signed_char), intent(out) :: ivf(:), ivg(:)
    logical(c_bool), intent(out) :: threshold(:)
    integer :: idx, x_datasize, y_datasize
    real(real32) :: xx_temp, yy_temp, zz_temp, rx, ry, rz, xp, yp, zp, vf, vg, r
    x_datasize = v - 2
    y_datasize = v - 2
    do idx = 1, data_size
      xx_temp = real(mod(idx - 1, x_datasize), real32) + 1.0_real32
      yy_temp = real(mod(int(floor(real(idx - 1, real32) / real(x_datasize, real32))), y_datasize), real32) + 1.0_real32
      zz_temp = floor(real(idx - 1, real32) / real(x_datasize, real32)) / real(y_datasize, real32) + 1.0_real32
      r = ran(mod(idx - 1, ran_count) + 1)
      rx = xx_temp + r
      ry = yy_temp + r
      rz = zz_temp + r
      xp = matrix(1) * rx + matrix(5) * ry + matrix(9) * rz + matrix(13)
      yp = matrix(2) * rx + matrix(6) * ry + matrix(10) * rz + matrix(14)
      zp = matrix(3) * rx + matrix(7) * ry + matrix(11) * rz + matrix(15)
      if (zp >= 1.0_real32 .and. zp < real(v, real32) .and. yp >= 1.0_real32 .and. yp < real(v, real32) &
          .and. xp >= 1.0_real32 .and. xp < real(v, real32)) then
        call interp_inline(f, v, xp, yp, zp, vf)
        call interp_inline(g, v, rx, ry, rz, vg)
        ivf(idx) = int(floor(vf + 0.5_real32), c_signed_char)
        ivg(idx) = int(floor(vg + 0.5_real32), c_signed_char)
        threshold(idx) = .true.
      else
        ivf(idx) = 0
        ivg(idx) = 0
        threshold(idx) = .false.
      end if
    end do
  end subroutine spm_host

  subroutine interp_inline(field, v, x, y, z, value)
    integer(c_signed_char), intent(in) :: field(:)
    integer, intent(in) :: v
    real(real32), intent(in) :: x, y, z
    real(real32), intent(out) :: value
    integer :: ix, iy, iz, base, upper_base
    integer :: k111, k112, k121, k122, k211, k212, k221, k222
    real(real32) :: dx1, dy1, dz1, dx2, dy2, dz2
    ix = int(floor(x))
    iy = int(floor(y))
    iz = int(floor(z))
    dx1 = x - real(ix, real32)
    dy1 = y - real(iy, real32)
    dz1 = z - real(iz, real32)
    dx2 = 1.0_real32 - dx1
    dy2 = 1.0_real32 - dy1
    dz2 = 1.0_real32 - dz1
    base = ix + v * ((iy - 1) + v * (iz - 1))
    k222 = byte_value(field(base))
    k122 = byte_value(field(base + 1))
    k212 = byte_value(field(base + v))
    k112 = byte_value(field(base + v + 1))
    upper_base = base + v * v
    k221 = byte_value(field(upper_base))
    k121 = byte_value(field(upper_base + 1))
    k211 = byte_value(field(upper_base + v))
    k111 = byte_value(field(upper_base + v + 1))
    value = (((real(k222, real32) * dx2 + real(k122, real32) * dx1) * dy2 + &
             (real(k212, real32) * dx2 + real(k112, real32) * dx1) * dy1)) * dz2 + &
            (((real(k221, real32) * dx2 + real(k121, real32) * dx1) * dy2 + &
             (real(k211, real32) * dx2 + real(k111, real32) * dx1) * dy1)) * dz1
  end subroutine interp_inline

  integer function build_histogram(ivf, ivg, threshold, histogram)
    integer(c_signed_char), intent(in) :: ivf(:), ivg(:)
    logical(c_bool), intent(in) :: threshold(:)
    integer(int32), intent(inout) :: histogram(:)
    integer :: idx, bin
    build_histogram = 0
    do idx = 1, size(threshold)
      if (threshold(idx)) then
        bin = byte_value(ivf(idx)) + byte_value(ivg(idx)) * 256 + 1
        histogram(bin) = histogram(bin) + 1
        build_histogram = build_histogram + 1
      end if
    end do
  end function build_histogram

  integer function byte_value(value)
    integer(c_signed_char), intent(in) :: value
    byte_value = iand(int(value, int32), 255)
  end function byte_value

end program main
