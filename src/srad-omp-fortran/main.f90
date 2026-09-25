! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  integer, parameter :: number_threads = 256
  real(real32), parameter :: tolerance = 5.0e-2_real32

  integer :: niter, nr, nc, ne, image_ori_rows, image_ori_cols, image_ori_elem
  integer :: iter, i, r1, r2, c1, c2, ne_roi, blocks_x, blocks_work_size
  integer, allocatable :: iN(:), iS(:), jW(:), jE(:)
  real(real32) :: lambda
  real(real32), allocatable :: image_ori(:), image(:)
  real(real32), allocatable :: dN(:), dS(:), dW(:), dE(:), c(:), sums(:), sums2(:)
  real(real64) :: t(0:12), stage_start, stage_end

  t = 0.0_real64
  t(0) = omp_get_wtime()
  t(1) = omp_get_wtime()

  if (command_argument_count() /= 4) then
    call print_usage()
    stop 1
  end if
  niter = read_int_arg(1)
  lambda = read_real_arg(2)
  nr = read_int_arg(3)
  nc = read_int_arg(4)
  if (niter < 0 .or. nr <= 0 .or. nc <= 0) stop 1
  t(2) = omp_get_wtime()

  image_ori_rows = 502
  image_ori_cols = 458
  image_ori_elem = image_ori_rows * image_ori_cols
  allocate(image_ori(image_ori_elem))
  call read_pgm('../data/srad/image.pgm', image_ori, image_ori_rows, image_ori_cols)
  t(3) = omp_get_wtime()

  ne = nr * nc
  blocks_x = ne / number_threads
  if (mod(ne, number_threads) /= 0) blocks_x = blocks_x + 1
  blocks_work_size = blocks_x
  allocate(image(ne))
  call resize_colmajor(image_ori, image_ori_rows, image_ori_cols, image, nr, nc)
  t(4) = omp_get_wtime()

  r1 = 0
  r2 = nr - 1
  c1 = 0
  c2 = nc - 1
  ne_roi = (r2 - r1 + 1) * (c2 - c1 + 1)
  allocate(iN(nr), iS(nr), jW(nc), jE(nc))
  do i = 1, nr
    iN(i) = i - 2
    iS(i) = i
  end do
  do i = 1, nc
    jW(i) = i - 2
    jE(i) = i
  end do
  iN(1) = 0
  iS(nr) = nr - 1
  jW(1) = 0
  jE(nc) = nc - 1
  allocate(dN(ne), dS(ne), dW(ne), dE(ne), c(ne), sums(ne), sums2(ne))
  t(5) = omp_get_wtime()

  !$omp target data map(to: iN(1:nr), iS(1:nr), jE(1:nc), jW(1:nc)) &
  !$omp& map(tofrom: image(1:ne)) &
  !$omp& map(alloc: dN(1:ne), dS(1:ne), dW(1:ne), dE(1:ne), c(1:ne), sums(1:ne), sums2(1:ne))
  t(6) = omp_get_wtime()

  !$omp target teams distribute parallel do num_teams(blocks_work_size) thread_limit(number_threads)
  do i = 1, ne
    image(i) = exp(image(i) / 255.0_real32)
  end do
  !$omp end target teams distribute parallel do
  t(7) = omp_get_wtime()

  stage_start = omp_get_wtime()
  do iter = 1, niter
    call srad_iteration_device(image, dN, dS, dW, dE, c, sums, sums2, iN, iS, jW, jE, &
                               nr, nc, ne, ne_roi, blocks_work_size, lambda)
  end do
  stage_end = omp_get_wtime()
  t(8) = t(7) + (stage_end - stage_start)

  !$omp target teams distribute parallel do num_teams(blocks_work_size) thread_limit(number_threads)
  do i = 1, ne
    image(i) = log(image(i)) * 255.0_real32
  end do
  !$omp end target teams distribute parallel do
  t(9) = t(8) + (omp_get_wtime() - stage_end)

  !$omp target update from(image(1:ne))
  t(10) = omp_get_wtime()
  !$omp end target data

  call write_pgm('image_out.pgm', image, nr, nc)
  t(11) = omp_get_wtime()

  deallocate(image_ori, image, iN, iS, jW, jE, dN, dS, dW, dE, c, sums, sums2)
  t(12) = omp_get_wtime()

  call print_timing(t, niter)

contains

  subroutine print_usage()
    character(len=256) :: arg0
    call get_command_argument(0, arg0)
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat> <lambda> <number of rows> <number of columns>'
  end subroutine print_usage

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  real(real32) function read_real_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_real_arg

  subroutine read_pgm(path, image, rows, cols)
    character(len=*), intent(in) :: path
    real(real32), intent(out) :: image(:)
    integer, intent(in) :: rows, cols
    character(len=16) :: magic
    integer :: unit, file_cols, file_rows, max_value, row, col, temp
    integer, allocatable :: row_values(:)
    open(newunit=unit, file=path, status='old', action='read', form='formatted')
    read(unit, *) magic
    read(unit, *) file_cols, file_rows
    read(unit, *) max_value
    if (trim(magic) /= 'P2' .or. file_rows /= rows .or. file_cols /= cols .or. max_value /= 255) then
      close(unit)
      error stop 'unexpected PGM header'
    end if
    allocate(row_values(cols))
    do row = 1, rows
      read(unit, *) row_values
      do col = 1, cols
        temp = row_values(col)
        image(row + (col - 1) * rows) = real(temp, real32)
      end do
    end do
    deallocate(row_values)
    close(unit)
  end subroutine read_pgm

  subroutine resize_colmajor(input, input_rows, input_cols, output, output_rows, output_cols)
    real(real32), intent(in) :: input(:)
    integer, intent(in) :: input_rows, input_cols, output_rows, output_cols
    real(real32), intent(out) :: output(:)
    integer :: row, col, row2, col2
    col2 = 1
    do col = 1, output_cols
      if (col2 > input_cols) col2 = col2 - input_cols
      row2 = 1
      do row = 1, output_rows
        if (row2 > input_rows) row2 = row2 - input_rows
        output(row + (col - 1) * output_rows) = input(row2 + (col2 - 1) * input_rows)
        row2 = row2 + 1
      end do
      col2 = col2 + 1
    end do
  end subroutine resize_colmajor

  subroutine write_pgm(path, image, rows, cols)
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: image(:)
    integer, intent(in) :: rows, cols
    integer :: unit, row, col
    open(newunit=unit, file=path, status='replace', action='write', form='formatted')
    write(unit,'(A)') 'P2'
    write(unit,'(I0,1X,I0)') cols, rows
    write(unit,'(I0)') 255
    do row = 1, rows
      do col = 1, cols
        write(unit,'(I0,1X)', advance='no') int(image(row + (col - 1) * rows))
      end do
      write(unit,*)
    end do
    close(unit)
  end subroutine write_pgm

  subroutine srad_iteration_device(image, dN, dS, dW, dE, c, sums, sums2, iN, iS, jW, jE, &
                                   nr, nc, ne, ne_roi, blocks_work_size, lambda)
    real(real32), intent(inout) :: image(:)
    real(real32), intent(inout) :: dN(:), dS(:), dW(:), dE(:), c(:), sums(:), sums2(:)
    integer, intent(in) :: iN(:), iS(:), jW(:), jE(:), nr, nc, ne, ne_roi, blocks_work_size
    real(real32), intent(in) :: lambda
    integer :: ei, row, col, blocks_work_size2, blocks_x, no, mul, bx, nf, j
    real(real32) :: psum, psum2
    real(real32) :: mean_roi, mean_roi2, var_roi, q0sqr
    real(real32) :: jc, n_loc, s_loc, w_loc, e_loc, g2, lap, num, den, qsqr, c_loc
    real(real32) :: cN, cS, cW, cE, div

    !$omp target teams distribute parallel do num_teams(blocks_work_size) thread_limit(number_threads)
    do ei = 1, ne
      sums(ei) = image(ei)
      sums2(ei) = image(ei) * image(ei)
    end do
    !$omp end target teams distribute parallel do

    blocks_work_size2 = blocks_work_size
    no = ne
    mul = 1
    do while (blocks_work_size2 /= 0)
      !$omp target teams distribute parallel do num_teams(blocks_work_size2) thread_limit(number_threads) &
      !$omp& private(nf, j, psum, psum2)
      do bx = 0, blocks_work_size2 - 1
        nf = number_threads
        if (bx == blocks_work_size2 - 1) nf = number_threads - (blocks_work_size2 * number_threads - no)
        psum = 0.0_real32
        psum2 = 0.0_real32
        do j = 1, nf
          ei = bx * number_threads + j
          psum = psum + sums((ei - 1) * mul + 1)
          psum2 = psum2 + sums2((ei - 1) * mul + 1)
        end do
        sums(bx * mul * number_threads + 1) = psum
        sums2(bx * mul * number_threads + 1) = psum2
      end do
      !$omp end target teams distribute parallel do

      no = blocks_work_size2
      if (blocks_work_size2 == 1) then
        blocks_work_size2 = 0
      else
        mul = mul * number_threads
        blocks_x = blocks_work_size2 / number_threads
        if (mod(blocks_work_size2, number_threads) /= 0) blocks_x = blocks_x + 1
        blocks_work_size2 = blocks_x
      end if
    end do

    !$omp target update from(sums(1:1))
    !$omp target update from(sums2(1:1))

    mean_roi = sums(1) / real(ne_roi, real32)
    mean_roi2 = mean_roi * mean_roi
    var_roi = (sums2(1) / real(ne_roi, real32)) - mean_roi2
    q0sqr = var_roi / mean_roi2

    !$omp target teams distribute parallel do num_teams(blocks_work_size) thread_limit(number_threads) &
    !$omp& private(row, col, jc, n_loc, s_loc, w_loc, e_loc, g2, lap, num, den, qsqr, c_loc)
    do ei = 1, ne
      row = mod(ei - 1, nr)
      col = (ei - 1) / nr
      jc = image(ei)
      n_loc = image(iN(row + 1) + nr * col + 1) - jc
      s_loc = image(iS(row + 1) + nr * col + 1) - jc
      w_loc = image(row + nr * jW(col + 1) + 1) - jc
      e_loc = image(row + nr * jE(col + 1) + 1) - jc
      g2 = (n_loc*n_loc + s_loc*s_loc + w_loc*w_loc + e_loc*e_loc) / (jc*jc)
      lap = (n_loc + s_loc + w_loc + e_loc) / jc
      num = 0.5_real32 * g2 - (1.0_real32 / 16.0_real32) * (lap * lap)
      den = 1.0_real32 + 0.25_real32 * lap
      qsqr = num / (den * den)
      den = (qsqr - q0sqr) / (q0sqr * (1.0_real32 + q0sqr))
      c_loc = 1.0_real32 / (1.0_real32 + den)
      if (c_loc < 0.0_real32) then
        c_loc = 0.0_real32
      else if (c_loc > 1.0_real32) then
        c_loc = 1.0_real32
      end if
      dN(ei) = n_loc
      dS(ei) = s_loc
      dW(ei) = w_loc
      dE(ei) = e_loc
      c(ei) = c_loc
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do thread_limit(number_threads) private(row, col, cN, cS, cW, cE, div)
    do ei = 1, ne
      row = mod(ei - 1, nr)
      col = (ei - 1) / nr
      cN = c(ei)
      cS = c(iS(row + 1) + nr * col + 1)
      cW = c(ei)
      cE = c(row + nr * jE(col + 1) + 1)
      div = cN * dN(ei) + cS * dS(ei) + cW * dW(ei) + cE * dE(ei)
      image(ei) = image(ei) + 0.25_real32 * lambda * div
    end do
    !$omp end target teams distribute parallel do
  end subroutine srad_iteration_device

  subroutine run_reference(image, iN, iS, jW, jE, nr, nc, ne, ne_roi, niter, lambda)
    real(real32), intent(inout) :: image(:)
    integer, intent(in) :: iN(:), iS(:), jW(:), jE(:), nr, nc, ne, ne_roi, niter
    real(real32), intent(in) :: lambda
    real(real32), allocatable :: dN(:), dS(:), dW(:), dE(:), c(:)
    integer :: iter, ei, row, col
    real(real32) :: sum1, sum2, mean_roi, mean_roi2, var_roi, q0sqr
    real(real32) :: jc, n_loc, s_loc, w_loc, e_loc, g2, lap, num, den, qsqr, c_loc
    real(real32) :: cN, cS, cW, cE, div
    allocate(dN(ne), dS(ne), dW(ne), dE(ne), c(ne))
    do ei = 1, ne
      image(ei) = exp(image(ei) / 255.0_real32)
    end do
    do iter = 1, niter
      sum1 = 0.0_real32
      sum2 = 0.0_real32
      do ei = 1, ne
        sum1 = sum1 + image(ei)
        sum2 = sum2 + image(ei) * image(ei)
      end do
      mean_roi = sum1 / real(ne_roi, real32)
      mean_roi2 = mean_roi * mean_roi
      var_roi = (sum2 / real(ne_roi, real32)) - mean_roi2
      q0sqr = var_roi / mean_roi2
      do ei = 1, ne
        row = mod(ei - 1, nr)
        col = (ei - 1) / nr
        jc = image(ei)
        n_loc = image(iN(row + 1) + nr * col + 1) - jc
        s_loc = image(iS(row + 1) + nr * col + 1) - jc
        w_loc = image(row + nr * jW(col + 1) + 1) - jc
        e_loc = image(row + nr * jE(col + 1) + 1) - jc
        g2 = (n_loc*n_loc + s_loc*s_loc + w_loc*w_loc + e_loc*e_loc) / (jc*jc)
        lap = (n_loc + s_loc + w_loc + e_loc) / jc
        num = 0.5_real32 * g2 - (1.0_real32 / 16.0_real32) * (lap * lap)
        den = 1.0_real32 + 0.25_real32 * lap
        qsqr = num / (den * den)
        den = (qsqr - q0sqr) / (q0sqr * (1.0_real32 + q0sqr))
        c_loc = 1.0_real32 / (1.0_real32 + den)
        if (c_loc < 0.0_real32) c_loc = 0.0_real32
        if (c_loc > 1.0_real32) c_loc = 1.0_real32
        dN(ei) = n_loc
        dS(ei) = s_loc
        dW(ei) = w_loc
        dE(ei) = e_loc
        c(ei) = c_loc
      end do
      do ei = 1, ne
        row = mod(ei - 1, nr)
        col = (ei - 1) / nr
        cN = c(ei)
        cS = c(iS(row + 1) + nr * col + 1)
        cW = c(ei)
        cE = c(row + nr * jE(col + 1) + 1)
        div = cN * dN(ei) + cS * dS(ei) + cW * dW(ei) + cE * dE(ei)
        image(ei) = image(ei) + 0.25_real32 * lambda * div
      end do
    end do
    do ei = 1, ne
      image(ei) = log(image(ei)) * 255.0_real32
    end do
    deallocate(dN, dS, dW, dE, c)
  end subroutine run_reference

  logical function compare_images(image, ref_image, ne) result(ok)
    real(real32), intent(in) :: image(:), ref_image(:)
    integer, intent(in) :: ne
    integer :: ei
    ok = .true.
    do ei = 1, ne
      if (abs(image(ei) - ref_image(ei)) > tolerance) then
        ok = .false.
        return
      end if
    end do
  end function compare_images

  subroutine print_timing(t, niter)
    real(real64), intent(in) :: t(0:12)
    integer, intent(in) :: niter
    real(real64) :: total
    total = t(12) - t(0)
    write(*,'(A)') 'Time spent in different stages of the application:'
    call print_stage(t(1)-t(0), total, 'SETUP VARIABLES')
    call print_stage(t(2)-t(1), total, 'READ COMMAND LINE PARAMETERS')
    call print_stage(t(3)-t(2), total, 'READ IMAGE FROM FILE')
    call print_stage(t(4)-t(3), total, 'RESIZE IMAGE')
    call print_stage(t(5)-t(4), total, 'GPU DRIVER INIT, CPU/GPU SETUP, MEMORY ALLOCATION')
    call print_stage(t(6)-t(5), total, 'COPY DATA TO CPU->GPU')
    call print_stage(t(7)-t(6), total, 'EXTRACT IMAGE')
    call print_stage_iter(t(8)-t(7), total, 'COMPUTE', niter)
    call print_stage(t(9)-t(8), total, 'COMPRESS IMAGE')
    call print_stage(t(10)-t(9), total, 'COPY DATA TO GPU->CPU')
    call print_stage(t(11)-t(10), total, 'SAVE IMAGE INTO FILE')
    call print_stage(t(12)-t(11), total, 'FREE MEMORY')
    write(*,'(A)') 'Total time:'
    write(*,'(F0.12,A)') total, ' s'
  end subroutine print_timing

  subroutine print_stage(seconds, total, label)
    real(real64), intent(in) :: seconds, total
    character(len=*), intent(in) :: label
    write(*,'(F15.12,A,F15.12,A,A)') seconds, ' s, ', seconds / total * 100.0_real64, ' % : ', trim(label)
  end subroutine print_stage

  subroutine print_stage_iter(seconds, total, label, niter)
    real(real64), intent(in) :: seconds, total
    character(len=*), intent(in) :: label
    integer, intent(in) :: niter
    write(*,'(F15.12,A,F15.12,A,A,A,I0,A)') seconds, ' s, ', seconds / total * 100.0_real64, &
        ' % : ', trim(label), ' (', niter, ' iterations)'
  end subroutine print_stage_iter

end program main
