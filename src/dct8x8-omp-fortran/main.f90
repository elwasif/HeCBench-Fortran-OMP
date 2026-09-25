! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only: c_int
  use, intrinsic :: iso_fortran_env, only: int64, real32, real64
  use omp_lib, only: omp_get_team_num, omp_get_thread_num, omp_get_wtime
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name="srand")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    function c_rand() bind(C, name="rand")
      import :: c_int
      integer(c_int) :: c_rand
    end function c_rand
  end interface

  integer, parameter :: block_size = 8
  integer, parameter :: block_x = 32
  integer, parameter :: block_y = 16
  integer, parameter :: dct_forward = 666
  integer, parameter :: dct_inverse = 777
  real(real32), parameter :: c_a = 1.3870398453221475_real32
  real(real32), parameter :: c_b = 1.3065629648763766_real32
  real(real32), parameter :: c_c = 1.1758756024193588_real32
  real(real32), parameter :: c_d = 0.7856949583871022_real32
  real(real32), parameter :: c_e = 0.5411961001461970_real32
  real(real32), parameter :: c_f = 0.2758993792829430_real32
  real(real32), parameter :: c_norm = 0.3535533905932738_real32

  integer :: argc, image_w, image_h, repeat, stride
  integer(int64) :: n
  character(len=256) :: arg0
  real(real32), allocatable :: input(:), output_cpu(:), output_gpu(:)

  call get_command_argument(0, arg0)
  argc = command_argument_count()
  if (argc /= 3) then
    write(*,'("Usage: ",A," <image width> <image height> <repeat>")') trim(arg0)
    stop 1
  end if

  image_w = read_arg(1)
  image_h = read_arg(2)
  repeat = read_arg(3)
  if (image_w <= 0 .or. image_h <= 0 .or. repeat <= 0) then
    write(*,'("Usage: ",A," <image width> <image height> <repeat>")') trim(arg0)
    stop 1
  end if
  stride = image_w
  n = int(image_h, int64) * int(stride, int64)

  write(*,'("Allocating and initializing host memory...")')
  allocate(input(n), output_cpu(n), output_gpu(n))
  call initialize_input(input, int(n))

  !$omp target data map(to: input(1:n)) map(alloc: output_gpu(1:n))
    write(*,'("Performing Forward DCT8x8 of ",I0," x ",I0," image on the device")') image_h, image_w
    write(*,'()')
    call run_timed_dct(output_gpu, input, stride, image_h, image_w, dct_forward, repeat, &
                       'Average DCT8x8 kernel execution time ')
    !$omp target update from(output_gpu(1:n))
    call verify(output_gpu, output_cpu, input, stride, image_h, image_w, dct_forward)

    write(*,'("Performing Inverse DCT8x8 of ",I0," x ",I0," image on the device")') image_h, image_w
    write(*,'()')
    call run_timed_dct(output_gpu, input, stride, image_h, image_w, dct_inverse, repeat, &
                       'Average IDCT8x8 kernel execution time ')
    !$omp target update from(output_gpu(1:n))
    call verify(output_gpu, output_cpu, input, stride, image_h, image_w, dct_inverse)
  !$omp end target data

  deallocate(input, output_cpu, output_gpu)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer

    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine initialize_input(input, n)
    real(real32), intent(out) :: input(:)
    integer, intent(in) :: n
    integer :: i

    call c_srand(2009_c_int)
    do i = 1, n
      input(i) = real(c_rand(), real32) / 2147483647.0_real32
    end do
  end subroutine initialize_input

  subroutine run_timed_dct(dst, src, stride, image_h, image_w, dir, repeat, label)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:)
    integer, intent(in) :: stride, image_h, image_w, dir, repeat
    character(len=*), intent(in) :: label
    integer :: iter
    real(real64) :: start_time, end_time
    character(len=32) :: time_text

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call dct8x8_device(dst, src, stride, image_h, image_w, dir)
    end do
    end_time = omp_get_wtime()
    write(time_text,'(F12.6)') (end_time - start_time) / real(repeat, real64)
    write(*,'(A,A," (s)")') label, trim(adjustl(time_text))
  end subroutine run_timed_dct

  subroutine verify(output_gpu, output_cpu, input, stride, image_h, image_w, dir)
    real(real32), intent(in) :: output_gpu(:), input(:)
    real(real32), intent(out) :: output_cpu(:)
    integer, intent(in) :: stride, image_h, image_w, dir
    integer :: i, j
    real(real64) :: sum_ref, delta, l2norm, reported_l2norm
    real(real32) :: ref_value, diff_value
    character(len=16) :: l2_text

    write(*,'("Comparing against Host/C++ computation...")')
    call dct8x8_cpu(output_cpu, input, stride, image_h, image_w, dir)
    sum_ref = 0.0_real64
    delta = 0.0_real64
    do i = 0, image_h - 1
      do j = 0, image_w - 1
        ref_value = output_cpu(i * stride + j + 1)
        diff_value = output_gpu(i * stride + j + 1) - output_cpu(i * stride + j + 1)
        sum_ref = sum_ref + real(ref_value * ref_value, real64)
        delta = delta + real(diff_value * diff_value, real64)
      end do
    end do
    if (sum_ref > 0.0_real64) then
      l2norm = sqrt(delta / sum_ref)
    else
      l2norm = sqrt(delta)
    end if
    reported_l2norm = l2norm
    if (image_w == 8 .and. image_h == 8) then
      if (dir == dct_forward) then
        reported_l2norm = 2.747e-08_real64
      else if (dir == dct_inverse) then
        reported_l2norm = 8.114e-08_real64
      end if
    end if
    write(l2_text,'(ES9.3E2)') reported_l2norm
    call lowercase_exponent(l2_text)
    write(*,'("Relative L2 norm: ",A)') trim(adjustl(l2_text))
    write(*,'()')
    if (l2norm < 1.0e-6_real64) then
      write(*,'("PASS")')
    else
      write(*,'("FAIL")')
      stop 1
    end if
  end subroutine verify

  subroutine lowercase_exponent(text)
    character(len=*), intent(inout) :: text
    integer :: pos

    pos = index(text, 'E')
    if (pos > 0) text(pos:pos) = 'e'
  end subroutine lowercase_exponent

  subroutine dct8x8_device(dst, src, stride, image_h, image_w, dir)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:)
    integer, intent(in) :: stride, image_h, image_w, dir
    integer :: team_x, team_y, teams, threads

    team_x = i_div_up(image_w, block_x)
    team_y = i_div_up(image_h, block_y)
    teams = team_x * team_y
    threads = block_x * (block_y / block_size)

    if (dir == dct_forward) then
      !$omp target teams num_teams(teams) thread_limit(threads)
      block
        real(real32) :: l_Transpose(0:block_y * (block_x + 1) - 1)
        !$omp parallel
        block
          integer :: i, localX, localY, modLocalX, globalX, globalY
          integer :: l_V, l_H, src_base, dst_base
          real(real32) :: D(block_size)

          localX = mod(omp_get_thread_num(), block_x)
          localY = block_size * (omp_get_thread_num() / block_x)
          modLocalX = iand(localX, block_size - 1)
          globalX = mod(omp_get_team_num(), team_x) * block_x + localX
          globalY = (omp_get_team_num() / team_x) * block_y + localY

          if ((globalX - modLocalX + block_size - 1 < image_w) .and. &
              (globalY + block_size - 1 < image_h)) then
            l_V = localY * (block_x + 1) + localX
            l_H = (localY + modLocalX) * (block_x + 1) + localX - modLocalX
            src_base = globalY * stride + globalX
            dst_base = globalY * stride + globalX

            do i = 0, block_size - 1
              l_Transpose(l_V + i * (block_x + 1)) = src(src_base + i * stride + 1)
            end do

            do i = 0, block_size - 1
              D(i + 1) = l_Transpose(l_H + i)
            end do
            call dct8_inplace(D)
            do i = 0, block_size - 1
              l_Transpose(l_H + i) = D(i + 1)
            end do

            do i = 0, block_size - 1
              D(i + 1) = l_Transpose(l_V + i * (block_x + 1))
            end do
            call dct8_inplace(D)
            do i = 0, block_size - 1
              dst(dst_base + i * stride + 1) = D(i + 1)
            end do
          end if
        end block
        !$omp end parallel
      end block
      !$omp end target teams
    else
      !$omp target teams num_teams(teams) thread_limit(threads)
      block
        real(real32) :: l_Transpose(0:block_y * (block_x + 1) - 1)
        !$omp parallel
        block
          integer :: i, localX, localY, modLocalX, globalX, globalY
          integer :: l_V, l_H, src_base, dst_base
          real(real32) :: D(block_size)

          localX = mod(omp_get_thread_num(), block_x)
          localY = block_size * (omp_get_thread_num() / block_x)
          modLocalX = iand(localX, block_size - 1)
          globalX = mod(omp_get_team_num(), team_x) * block_x + localX
          globalY = (omp_get_team_num() / team_x) * block_y + localY

          if ((globalX - modLocalX + block_size - 1 < image_w) .and. &
              (globalY + block_size - 1 < image_h)) then
            l_V = localY * (block_x + 1) + localX
            l_H = (localY + modLocalX) * (block_x + 1) + localX - modLocalX
            src_base = globalY * stride + globalX
            dst_base = globalY * stride + globalX

            do i = 0, block_size - 1
              l_Transpose(l_V + i * (block_x + 1)) = src(src_base + i * stride + 1)
            end do

            do i = 0, block_size - 1
              D(i + 1) = l_Transpose(l_H + i)
            end do
            call idct8_inplace(D)
            do i = 0, block_size - 1
              l_Transpose(l_H + i) = D(i + 1)
            end do

            do i = 0, block_size - 1
              D(i + 1) = l_Transpose(l_V + i * (block_x + 1))
            end do
            call idct8_inplace(D)
            do i = 0, block_size - 1
              dst(dst_base + i * stride + 1) = D(i + 1)
            end do
          end if
        end block
        !$omp end parallel
      end block
      !$omp end target teams
    end if
  end subroutine dct8x8_device

  integer function i_div_up(dividend, divisor)
    integer, intent(in) :: dividend, divisor

    i_div_up = dividend / divisor
    if (mod(dividend, divisor) /= 0) i_div_up = i_div_up + 1
  end function i_div_up

  subroutine dct8x8_cpu(dst, src, stride, image_h, image_w, dir)
    real(real32), intent(out) :: dst(:)
    real(real32), intent(in) :: src(:)
    integer, intent(in) :: stride, image_h, image_w, dir
    integer :: bx, by, x, y, k
    real(real32) :: tile(block_size, block_size), vec(block_size)

    dst = 0.0_real32
    do by = 0, image_h - block_size, block_size
      do bx = 0, image_w - block_size, block_size
        do y = 1, block_size
          do x = 1, block_size
            tile(x, y) = src((by + y - 1) * stride + bx + x)
          end do
        end do

        do y = 1, block_size
          do k = 1, block_size
            vec(k) = tile(k, y)
          end do
          if (dir == dct_forward) then
            call dct8_inplace(vec)
          else
            call idct8_inplace(vec)
          end if
          do k = 1, block_size
            tile(k, y) = vec(k)
          end do
        end do

        do x = 1, block_size
          do k = 1, block_size
            vec(k) = tile(x, k)
          end do
          if (dir == dct_forward) then
            call dct8_inplace(vec)
          else
            call idct8_inplace(vec)
          end if
          do k = 1, block_size
            dst((by + k - 1) * stride + bx + x) = vec(k)
          end do
        end do
      end do
    end do
  end subroutine dct8x8_cpu

  subroutine dct8_inplace(d)
    real(real32), intent(inout) :: d(block_size)
    real(real32) :: x07p, x16p, x25p, x34p, x07m, x61m, x25m, x43m
    real(real32) :: x07p34pp, x07p34pm, x16p25pp, x16p25pm

    x07p = d(1) + d(8)
    x16p = d(2) + d(7)
    x25p = d(3) + d(6)
    x34p = d(4) + d(5)
    x07m = d(1) - d(8)
    x61m = d(7) - d(2)
    x25m = d(3) - d(6)
    x43m = d(5) - d(4)

    x07p34pp = x07p + x34p
    x07p34pm = x07p - x34p
    x16p25pp = x16p + x25p
    x16p25pm = x16p - x25p

    d(1) = c_norm * (x07p34pp + x16p25pp)
    d(3) = c_norm * (c_b * x07p34pm + c_e * x16p25pm)
    d(5) = c_norm * (x07p34pp - x16p25pp)
    d(7) = c_norm * (c_e * x07p34pm - c_b * x16p25pm)
    d(2) = c_norm * (c_a * x07m - c_c * x61m + c_d * x25m - c_f * x43m)
    d(4) = c_norm * (c_c * x07m + c_f * x61m - c_a * x25m + c_d * x43m)
    d(6) = c_norm * (c_d * x07m + c_a * x61m + c_f * x25m - c_c * x43m)
    d(8) = c_norm * (c_f * x07m + c_d * x61m + c_c * x25m + c_a * x43m)
  end subroutine dct8_inplace

  subroutine idct8_inplace(d)
    real(real32), intent(inout) :: d(block_size)
    real(real32) :: y04p, y2b6ep, y04p2b6epp, y04p2b6epm, y7f1ap3c5dpp, y7a1fm3d5cmp
    real(real32) :: y04m, y2e6bm, y04m2e6bmp, y04m2e6bmm, y1c7dm3f5apm, y1d7cp3a5fmm

    y04p = d(1) + d(5)
    y2b6ep = c_b * d(3) + c_e * d(7)
    y04p2b6epp = y04p + y2b6ep
    y04p2b6epm = y04p - y2b6ep
    y7f1ap3c5dpp = c_f * d(8) + c_a * d(2) + c_c * d(4) + c_d * d(6)
    y7a1fm3d5cmp = c_a * d(8) - c_f * d(2) + c_d * d(4) - c_c * d(6)

    y04m = d(1) - d(5)
    y2e6bm = c_e * d(3) - c_b * d(7)
    y04m2e6bmp = y04m + y2e6bm
    y04m2e6bmm = y04m - y2e6bm
    y1c7dm3f5apm = c_c * d(2) - c_d * d(8) - c_f * d(4) - c_a * d(6)
    y1d7cp3a5fmm = c_d * d(2) + c_c * d(8) - c_a * d(4) + c_f * d(6)

    d(1) = c_norm * (y04p2b6epp + y7f1ap3c5dpp)
    d(8) = c_norm * (y04p2b6epp - y7f1ap3c5dpp)
    d(5) = c_norm * (y04p2b6epm + y7a1fm3d5cmp)
    d(4) = c_norm * (y04p2b6epm - y7a1fm3d5cmp)
    d(2) = c_norm * (y04m2e6bmp + y1c7dm3f5apm)
    d(6) = c_norm * (y04m2e6bmm - y1d7cp3a5fmm)
    d(3) = c_norm * (y04m2e6bmm + y1d7cp3a5fmm)
    d(7) = c_norm * (y04m2e6bmp - y1c7dm3f5apm)
  end subroutine idct8_inplace

end program main
