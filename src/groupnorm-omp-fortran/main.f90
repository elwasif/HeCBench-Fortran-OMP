! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: tpb = 1024
  real(real32), parameter :: eps = 1.0e-5_real32

  integer :: bsz, channels, height, width, n_groups, repeat
  integer :: img_size, n, ng, total, i
  character(len=256) :: progname
  real(real32), allocatable :: x(:), weight(:), bias(:), dout(:)
  real(real32), allocatable :: out_ref(:), mean_ref(:), rstd_ref(:)
  real(real32), allocatable :: dx_ref(:), dweight_ref(:), dbias_ref(:)
  real(real32), allocatable :: out_dev(:), mean_dev(:), rstd_dev(:)
  real(real32), allocatable :: dx_dev(:), dweight_dev(:), dbias_dev(:)
  real(real64) :: elapsed

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name='rand')
      import :: c_int
    end function c_rand
  end interface

  if (command_argument_count() /= 6) then
    call get_command_argument(0, progname)
    print '(A,A,A)', 'Usage: ', trim(progname), ' <batch size> <number of channels> <height> <width> <number of groups> <repeat>'
    stop 1
  end if

  bsz = read_arg(1)
  channels = read_arg(2)
  height = read_arg(3)
  width = read_arg(4)
  n_groups = read_arg(5)
  repeat = read_arg(6)
  if (bsz <= 0 .or. channels <= 0 .or. height <= 0 .or. width <= 0 .or. n_groups <= 0 .or. repeat <= 0) error stop 'invalid arguments'
  if (mod(channels, n_groups) /= 0) error stop 'channels must be divisible by number of groups'

  img_size = height * width
  total = bsz * channels * img_size
  ng = bsz * n_groups

  allocate(x(total), weight(channels), bias(channels), dout(total))
  allocate(out_ref(total), mean_ref(ng), rstd_ref(ng))
  allocate(dx_ref(total), dweight_ref(channels), dbias_ref(channels))
  allocate(out_dev(total), mean_dev(ng), rstd_dev(ng))
  allocate(dx_dev(total), dweight_dev(channels), dbias_dev(channels))

  call c_srand(0_c_int)
  call fill_random(x)
  call fill_random(weight)
  call fill_random(bias)
  call fill_random(dout)
  out_ref = 0.0_real32
  mean_ref = 0.0_real32
  rstd_ref = 0.0_real32
  dx_ref = 0.0_real32
  dweight_ref = 0.0_real32
  dbias_ref = 0.0_real32
  out_dev = 0.0_real32
  mean_dev = 0.0_real32
  rstd_dev = 0.0_real32
  dx_dev = 0.0_real32
  dweight_dev = 0.0_real32
  dbias_dev = 0.0_real32

  !$omp target data map(to: x(1:total), weight(1:channels), bias(1:channels), dout(1:total)) &
  !$omp& map(tofrom: out_dev(1:total), mean_dev(1:ng), rstd_dev(1:ng), dx_dev(1:total), dweight_dev(1:channels), dbias_dev(1:channels))
    print '(A)', 'Checking forward pass'

    call groupnorm_forward_ref(x, weight, bias, out_ref, mean_ref, rstd_ref, bsz, channels, img_size, n_groups)
    call groupnorm_forward_dev(x, weight, bias, out_dev, mean_dev, rstd_dev, bsz, channels, img_size, n_groups)
    !$omp target update from(out_dev(1:total))
    call validate_result(out_dev, out_ref, 'out', 1.0e-2_real32)

    print '(A)', 'Checking backward pass'

    call groupnorm_backward_ref(dout, x, mean_ref, rstd_ref, weight, dx_ref, dweight_ref, dbias_ref, bsz, channels, img_size, n_groups)
    call groupnorm_backward_dev(dout, x, mean_dev, rstd_dev, weight, dx_dev, dweight_dev, dbias_dev, bsz, channels, img_size, n_groups)

    print '(A)', 'Checking dbias'
    !$omp target update from(dbias_dev(1:channels))
    call validate_result(dbias_dev, dbias_ref, 'dbias', 1.0e-2_real32)
    print '(A)', 'Checking dweight'
    !$omp target update from(dweight_dev(1:channels))
    call validate_result(dweight_dev, dweight_ref, 'dweight', 1.0e-2_real32)
    print '(A)', 'Checking dx'
    !$omp target update from(dx_dev(1:total))
    call validate_result(dx_dev, dx_ref, 'dx', 1.0_real32)
    print '(A)', ''
    print '(A)', '─────────────────────────────────────────────────────'

    print '(A)', 'Forward pass benchmarks'
    elapsed = benchmark_forward(repeat, x, weight, bias, out_dev, mean_dev, rstd_dev, bsz, channels, img_size, n_groups)
    write(*, '(A,F0.4,A)') 'time ', elapsed, ' us'

    print '(A)', 'Backward pass benchmarks'
    elapsed = benchmark_backward(repeat, dout, x, mean_dev, rstd_dev, weight, dx_dev, dweight_dev, dbias_dev, bsz, channels, img_size, n_groups)
    write(*, '(A,F0.4,A)') 'time ', elapsed, ' us'
  !$omp end target data

  deallocate(x, weight, bias, dout, out_ref, mean_ref, rstd_ref, dx_ref, dweight_ref, dbias_ref)
  deallocate(out_dev, mean_dev, rstd_dev, dx_dev, dweight_dev, dbias_dev)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  subroutine fill_random(a)
    real(real32), intent(out) :: a(:)
    integer :: i
    integer(c_int) :: raw
    real(real32), parameter :: rand_max = 2147483647.0_real32
    do i = 1, size(a)
      raw = c_rand()
      a(i) = real(raw, real32) / rand_max * 2.0_real32 - 1.0_real32
    end do
  end subroutine fill_random

  subroutine validate_result(actual, expected, name, tolerance)
    real(real32), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: name
    integer :: i, nfaults
    nfaults = 0
    do i = 1, size(actual)
      if (abs(expected(i) - actual(i)) > tolerance .and. expected(i) == expected(i)) then
        write(*, '(A,A,A,I0,A,F0.6,A,F0.6)') 'Mismatch of ', trim(name), ' at ', i - 1, ': CPU_ref: ', expected(i), ' vs GPU: ', actual(i)
        nfaults = nfaults + 1
        if (nfaults >= 10) exit
      end if
    end do
    if (nfaults == 0) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end subroutine validate_result

  subroutine groupnorm_forward_ref(x, weight, bias, out, mean, rstd, bsz, channels, img_size, n_groups)
    real(real32), intent(in) :: x(:), weight(:), bias(:)
    real(real32), intent(out) :: out(:), mean(:), rstd(:)
    integer, intent(in) :: bsz, channels, img_size, n_groups
    integer :: group_size, group_pixels, b, g, i, c, block_idx, base, chan
    real(real32) :: sumv, sum2, m, var, s, val, norm

    group_size = channels / n_groups
    group_pixels = img_size * group_size
    do b = 1, bsz
      do g = 1, n_groups
        block_idx = (b - 1) * n_groups + g
        base = (block_idx - 1) * group_pixels
        sumv = 0.0_real32
        sum2 = 0.0_real32
        do i = 1, group_pixels
          val = x(base + i)
          sumv = sumv + val
          sum2 = sum2 + val * val
        end do
        m = sumv / real(group_pixels, real32)
        var = sum2 / real(group_pixels, real32) - m * m
        s = 1.0_real32 / sqrt(var + eps)
        mean(block_idx) = m
        rstd(block_idx) = s
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          norm = s * (x(base + i) - m)
          out(base + i) = norm * weight(chan) + bias(chan)
        end do
      end do
    end do
  end subroutine groupnorm_forward_ref

  subroutine groupnorm_backward_ref(dout, x, mean, rstd, weight, dx, dweight, dbias, bsz, channels, img_size, n_groups)
    real(real32), intent(in) :: dout(:), x(:), mean(:), rstd(:), weight(:)
    real(real32), intent(out) :: dx(:), dweight(:), dbias(:)
    integer, intent(in) :: bsz, channels, img_size, n_groups
    integer :: group_size, group_pixels, b, g, i, c, block_idx, base, chan, pix
    real(real32) :: m, s, wdout_sum, wdout_norm_sum, wdout_block, wdout_norm_block
    real(real32) :: norm, wdout, dw, db

    group_size = channels / n_groups
    group_pixels = img_size * group_size
    dx = 0.0_real32
    dweight = 0.0_real32
    dbias = 0.0_real32
    do b = 1, bsz
      do g = 1, n_groups
        block_idx = (b - 1) * n_groups + g
        base = (block_idx - 1) * group_pixels
        m = mean(block_idx)
        s = rstd(block_idx)
        wdout_sum = 0.0_real32
        wdout_norm_sum = 0.0_real32
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          wdout = weight(chan) * dout(base + i)
          norm = (x(base + i) - m) * s
          wdout_sum = wdout_sum + wdout
          wdout_norm_sum = wdout_norm_sum + wdout * norm
        end do
        wdout_block = wdout_sum / real(group_pixels, real32)
        wdout_norm_block = wdout_norm_sum / real(group_pixels, real32)
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          norm = (x(base + i) - m) * s
          wdout = weight(chan) * dout(base + i)
          dx(base + i) = (wdout - wdout_block - norm * wdout_norm_block) * s
        end do
        do c = 1, group_size
          chan = (g - 1) * group_size + c
          dw = 0.0_real32
          db = 0.0_real32
          do pix = 1, img_size
            i = (c - 1) * img_size + pix
            norm = (x(base + i) - m) * s
            db = db + dout(base + i)
            dw = dw + dout(base + i) * norm
          end do
          dweight(chan) = dweight(chan) + dw
          dbias(chan) = dbias(chan) + db
        end do
      end do
    end do
  end subroutine groupnorm_backward_ref

  subroutine groupnorm_forward_dev(x, weight, bias, out, mean, rstd, bsz, channels, img_size, n_groups)
    real(real32), intent(in) :: x(:), weight(:), bias(:)
    real(real32), intent(inout) :: out(:), mean(:), rstd(:)
    integer, intent(in) :: bsz, channels, img_size, n_groups
    integer :: group_size, group_pixels, team_size, block_size
    integer :: b, g, i, c, block_idx, base, chan
    real(real32) :: sumv, sum2, m, var, s, val, norm

    group_size = channels / n_groups
    group_pixels = img_size * group_size
    team_size = bsz * n_groups
    block_size = max(min(tpb, group_pixels), 32)
    !$omp target teams distribute collapse(2) num_teams(team_size) &
    !$omp& private(block_idx, base, sumv, sum2, m, var, s, i, c, chan, val, norm)
    do b = 1, bsz
      do g = 1, n_groups
        block_idx = (b - 1) * n_groups + g
        base = (block_idx - 1) * group_pixels
        sumv = 0.0_real32
        sum2 = 0.0_real32
        !$omp parallel do reduction(+:sumv, sum2) num_threads(block_size) private(val)
        do i = 1, group_pixels
          val = x(base + i)
          sumv = sumv + val
          sum2 = sum2 + val * val
        end do
        !$omp end parallel do
        m = sumv / real(group_pixels, real32)
        var = sum2 / real(group_pixels, real32) - m * m
        s = 1.0_real32 / sqrt(var + eps)
        mean(block_idx) = m
        rstd(block_idx) = s
        !$omp parallel do num_threads(block_size) private(c, chan, norm)
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          norm = s * (x(base + i) - m)
          out(base + i) = norm * weight(chan) + bias(chan)
        end do
        !$omp end parallel do
      end do
    end do
    !$omp end target teams distribute
  end subroutine groupnorm_forward_dev

  subroutine groupnorm_backward_dev(dout, x, mean, rstd, weight, dx, dweight, dbias, bsz, channels, img_size, n_groups)
    real(real32), intent(in) :: dout(:), x(:), mean(:), rstd(:), weight(:)
    real(real32), intent(inout) :: dx(:), dweight(:), dbias(:)
    integer, intent(in) :: bsz, channels, img_size, n_groups
    integer :: group_size, group_pixels, team_size, block_size
    integer :: b, g, i, c, block_idx, base, chan, pix
    real(real32) :: m, s, wdout_sum, wdout_norm_sum, wdout_block, wdout_norm_block
    real(real32) :: norm, wdout, dw, db

    group_size = channels / n_groups
    group_pixels = img_size * group_size
    team_size = bsz * n_groups
    block_size = max(min(tpb, group_pixels), 32 * group_size)
    !$omp target teams distribute collapse(2) num_teams(team_size) &
    !$omp& private(block_idx, base, m, s, wdout_sum, wdout_norm_sum, wdout_block, wdout_norm_block, i, c, chan, pix, norm, wdout, dw, db)
    do b = 1, bsz
      do g = 1, n_groups
        block_idx = (b - 1) * n_groups + g
        base = (block_idx - 1) * group_pixels
        m = mean(block_idx)
        s = rstd(block_idx)
        wdout_sum = 0.0_real32
        wdout_norm_sum = 0.0_real32
        !$omp parallel do reduction(+:wdout_sum, wdout_norm_sum) num_threads(block_size) private(c, chan, norm, wdout)
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          wdout = weight(chan) * dout(base + i)
          norm = (x(base + i) - m) * s
          wdout_sum = wdout_sum + wdout
          wdout_norm_sum = wdout_norm_sum + wdout * norm
        end do
        !$omp end parallel do
        wdout_block = wdout_sum / real(group_pixels, real32)
        wdout_norm_block = wdout_norm_sum / real(group_pixels, real32)
        !$omp parallel do num_threads(block_size) private(c, chan, norm, wdout)
        do i = 1, group_pixels
          c = (i - 1) / img_size + 1
          chan = (g - 1) * group_size + c
          norm = (x(base + i) - m) * s
          wdout = weight(chan) * dout(base + i)
          dx(base + i) = (wdout - wdout_block - norm * wdout_norm_block) * s
        end do
        !$omp end parallel do
        do c = 1, group_size
          chan = (g - 1) * group_size + c
          dw = 0.0_real32
          db = 0.0_real32
          !$omp parallel do reduction(+:dw, db) num_threads(block_size) private(i, norm)
          do pix = 1, img_size
            i = (c - 1) * img_size + pix
            norm = (x(base + i) - m) * s
            db = db + dout(base + i)
            dw = dw + dout(base + i) * norm
          end do
          !$omp end parallel do
          !$omp atomic update
          dweight(chan) = dweight(chan) + dw
          !$omp atomic update
          dbias(chan) = dbias(chan) + db
        end do
      end do
    end do
    !$omp end target teams distribute
  end subroutine groupnorm_backward_dev

  real(real64) function benchmark_forward(repeat, x, weight, bias, out, mean, rstd, bsz, channels, img_size, n_groups)
    integer, intent(in) :: repeat, bsz, channels, img_size, n_groups
    real(real32), intent(in) :: x(:), weight(:), bias(:)
    real(real32), intent(inout) :: out(:), mean(:), rstd(:)
    integer :: i
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call groupnorm_forward_dev(x, weight, bias, out, mean, rstd, bsz, channels, img_size, n_groups)
    end do
    end_time = omp_get_wtime()
    benchmark_forward = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
  end function benchmark_forward

  real(real64) function benchmark_backward(repeat, dout, x, mean, rstd, weight, dx, dweight, dbias, bsz, channels, img_size, n_groups)
    integer, intent(in) :: repeat, bsz, channels, img_size, n_groups
    real(real32), intent(in) :: dout(:), x(:), mean(:), rstd(:), weight(:)
    real(real32), intent(inout) :: dx(:), dweight(:), dbias(:)
    integer :: i
    real(real64) :: start_time, end_time
    start_time = omp_get_wtime()
    do i = 1, repeat
      call groupnorm_backward_dev(dout, x, mean, rstd, weight, dx, dweight, dbias, bsz, channels, img_size, n_groups)
    end do
    end_time = omp_get_wtime()
    benchmark_backward = ((end_time - start_time) * 1.0e6_real64) / real(repeat, real64)
  end function benchmark_backward

end program main
