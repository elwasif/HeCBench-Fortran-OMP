! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only: c_int
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use omp_lib, only: omp_get_wtime, omp_get_team_num, omp_get_thread_num
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

  integer, parameter :: kernel_radius = 8
  integer, parameter :: kernel_length = 2 * kernel_radius + 1
  integer, parameter :: rows_blockdim_x = 16
  integer, parameter :: rows_blockdim_y = 4
  integer, parameter :: rows_result_steps = 8
  integer, parameter :: rows_halo_steps = 1
  integer, parameter :: columns_blockdim_x = 16
  integer, parameter :: columns_blockdim_y = 8
  integer, parameter :: columns_result_steps = 8
  integer, parameter :: columns_halo_steps = 1

  integer :: argc, image_w, image_h, repeat
  integer(int64) :: n
  character(len=256) :: arg0
  real(real32), allocatable :: kernel(:), input(:), buffer(:), output_cpu(:), output_gpu(:)
  real(real64) :: sum_ref, delta, l2norm
  character(len=16) :: l2_text
  integer :: i

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
  if (mod(image_w, rows_result_steps * rows_blockdim_x) /= 0 .or. &
      mod(image_h, rows_blockdim_y) /= 0 .or. &
      mod(image_w, columns_blockdim_x) /= 0 .or. &
      mod(image_h, columns_result_steps * columns_blockdim_y) /= 0) then
    write(*,'("Image dimensions do not satisfy tiled convolution block geometry")')
    stop 1
  end if

  n = int(image_w, int64) * int(image_h, int64)
  if (n > huge(i)) then
    write(*,'("Image is too large")')
    stop 1
  end if

  allocate(kernel(kernel_length), input(n), buffer(n), output_cpu(n), output_gpu(n))
  call initialize_inputs(kernel, input, int(n))

  !$omp target data map(to: kernel(1:kernel_length), input(1:n)) &
  !$omp& map(alloc: buffer(1:n)) map(from: output_gpu(1:n))
    call convolution_rows_device(buffer, input, kernel, image_w, image_h)
    call convolution_columns_device(output_gpu, buffer, kernel, image_w, image_h)

    call run_timed_convolution(output_gpu, buffer, input, kernel, image_w, image_h, repeat)
  !$omp end target data

  write(*,'("Comparing against Host/C++ computation...")')
  call convolution_row_host(buffer, input, kernel, image_w, image_h)
  call convolution_column_host(output_cpu, buffer, kernel, image_w, image_h)

  delta = 0.0_real64
  sum_ref = 0.0_real64
  do i = 1, int(n)
    delta = delta + real(output_cpu(i) - output_gpu(i), real64) * real(output_cpu(i) - output_gpu(i), real64)
    sum_ref = sum_ref + real(output_cpu(i), real64) * real(output_cpu(i), real64)
  end do
  if (sum_ref > 0.0_real64) then
    l2norm = sqrt(delta / sum_ref)
  else
    l2norm = sqrt(delta)
  end if
  write(l2_text,'(ES9.3)') l2norm
  do i = 1, len(l2_text)
    if (l2_text(i:i) == 'E') l2_text(i:i) = 'e'
  end do
  write(*,'("Relative L2 norm: ",A)') trim(adjustl(l2_text))
  write(*,'()')

  if (l2norm < 1.0e-6_real64) then
    write(*,'("PASS")')
  else
    write(*,'("FAIL")')
    stop 1
  end if

  deallocate(kernel, input, buffer, output_cpu, output_gpu)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer_arg

    call get_command_argument(position, buffer_arg)
    read(buffer_arg, *) read_arg
  end function read_arg

  subroutine initialize_inputs(kernel, input, n)
    real(real32), intent(out) :: kernel(:), input(:)
    integer, intent(in) :: n
    integer :: idx

    call c_srand(2009_c_int)
    do idx = 1, kernel_length
      kernel(idx) = real(mod(c_rand(), 16_c_int), real32)
    end do

    do idx = 1, n
      input(idx) = real(mod(c_rand(), 16_c_int), real32)
    end do
  end subroutine initialize_inputs

  subroutine run_timed_convolution(output, buffer, input, kernel, image_w, image_h, repeat)
    real(real32), intent(inout) :: output(:), buffer(:)
    real(real32), intent(in) :: input(:), kernel(:)
    integer, intent(in) :: image_w, image_h, repeat
    integer :: iter
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    do iter = 1, repeat
      call convolution_rows_device(buffer, input, kernel, image_w, image_h)
      call convolution_columns_device(output, buffer, kernel, image_w, image_h)
    end do
    end_time = omp_get_wtime()

    write(*,'("Average kernel execution time ",F8.6," (s)")') (end_time - start_time) / real(repeat, real64)
  end subroutine run_timed_convolution

  subroutine convolution_rows_device(dst, src, kernel, image_w, image_h)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: team_x, team_y, num_teams
    real(real32), allocatable :: l_data(:, :, :)

    team_x = (image_w / rows_result_steps) / rows_blockdim_x
    team_y = image_h / rows_blockdim_y
    num_teams = team_x * team_y
    allocate(l_data(num_teams, rows_blockdim_y, (rows_result_steps + 2 * rows_halo_steps) * rows_blockdim_x))

    !$omp target teams num_teams(num_teams) thread_limit(rows_blockdim_y * rows_blockdim_x) map(alloc: l_data)
      !$omp parallel
      block
      integer :: gid_x, gid_y, lid_x, lid_y, base_x, base_y, i, j
      real(real32) :: accum

        gid_x = mod(omp_get_team_num(), team_x)
        gid_y = omp_get_team_num() / team_x
        lid_x = mod(omp_get_thread_num(), rows_blockdim_x)
        lid_y = omp_get_thread_num() / rows_blockdim_x
        base_x = (gid_x * rows_result_steps - rows_halo_steps) * rows_blockdim_x + lid_x
        base_y = gid_y * rows_blockdim_y + lid_y

        do i = rows_halo_steps, rows_halo_steps + rows_result_steps - 1
          l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + 1) = &
              src(base_y * image_w + base_x + i * rows_blockdim_x + 1)
        end do

        do i = 0, rows_halo_steps - 1
          if (base_x + i * rows_blockdim_x >= 0) then
            l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + 1) = &
                src(base_y * image_w + base_x + i * rows_blockdim_x + 1)
          else
            l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + 1) = 0.0_real32
          end if
        end do

        do i = rows_halo_steps + rows_result_steps, rows_halo_steps + rows_result_steps + rows_halo_steps - 1
          if (base_x + i * rows_blockdim_x < image_w) then
            l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + 1) = &
                src(base_y * image_w + base_x + i * rows_blockdim_x + 1)
          else
            l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + 1) = 0.0_real32
          end if
        end do

        !$omp barrier

        do i = rows_halo_steps, rows_halo_steps + rows_result_steps - 1
          accum = 0.0_real32
          do j = -kernel_radius, kernel_radius
            accum = accum + kernel(kernel_radius - j + 1) * &
                l_data(omp_get_team_num() + 1, lid_y + 1, lid_x + i * rows_blockdim_x + j + 1)
          end do
          dst(base_y * image_w + base_x + i * rows_blockdim_x + 1) = accum
        end do
      end block
      !$omp end parallel
    !$omp end target teams
    deallocate(l_data)
  end subroutine convolution_rows_device

  subroutine convolution_columns_device(dst, src, kernel, image_w, image_h)
    real(real32), intent(inout) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: team_x, team_y, num_teams
    real(real32), allocatable :: l_data(:, :, :)

    team_x = image_w / columns_blockdim_x
    team_y = image_h / columns_result_steps / columns_blockdim_y
    num_teams = team_x * team_y
    allocate(l_data(num_teams, columns_blockdim_x, (columns_result_steps + 2 * columns_halo_steps) * columns_blockdim_y + 1))

    !$omp target teams num_teams(num_teams) thread_limit(columns_blockdim_y * columns_blockdim_x) map(alloc: l_data)
      !$omp parallel
      block
      integer :: gid_x, gid_y, lid_x, lid_y, base_x, base_y, i, j
      real(real32) :: accum

        gid_x = mod(omp_get_team_num(), team_x)
        gid_y = omp_get_team_num() / team_x
        lid_x = mod(omp_get_thread_num(), columns_blockdim_x)
        lid_y = omp_get_thread_num() / columns_blockdim_x
        base_x = gid_x * columns_blockdim_x + lid_x
        base_y = (gid_y * columns_result_steps - columns_halo_steps) * columns_blockdim_y + lid_y

        do i = columns_halo_steps, columns_halo_steps + columns_result_steps - 1
          l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + 1) = &
              src((base_y + i * columns_blockdim_y) * image_w + base_x + 1)
        end do

        do i = 0, columns_halo_steps - 1
          if (base_y + i * columns_blockdim_y >= 0) then
            l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + 1) = &
                src((base_y + i * columns_blockdim_y) * image_w + base_x + 1)
          else
            l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + 1) = 0.0_real32
          end if
        end do

        do i = columns_halo_steps + columns_result_steps, columns_halo_steps + columns_result_steps + columns_halo_steps - 1
          if (base_y + i * columns_blockdim_y < image_h) then
            l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + 1) = &
                src((base_y + i * columns_blockdim_y) * image_w + base_x + 1)
          else
            l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + 1) = 0.0_real32
          end if
        end do

        !$omp barrier

        do i = columns_halo_steps, columns_halo_steps + columns_result_steps - 1
          accum = 0.0_real32
          do j = -kernel_radius, kernel_radius
            accum = accum + kernel(kernel_radius - j + 1) * &
                l_data(omp_get_team_num() + 1, lid_x + 1, lid_y + i * columns_blockdim_y + j + 1)
          end do
          dst((base_y + i * columns_blockdim_y) * image_w + base_x + 1) = accum
        end do
      end block
      !$omp end parallel
    !$omp end target teams
    deallocate(l_data)
  end subroutine convolution_columns_device

  subroutine convolution_row_host(dst, src, kernel, image_w, image_h)
    real(real32), intent(out) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d
    real(real64) :: accum

    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real64
        do k = -kernel_radius, kernel_radius
          d = x + k
          if (d >= 0 .and. d < image_w) then
            accum = accum + real(src(y * image_w + d + 1), real64) * real(kernel(kernel_radius - k + 1), real64)
          end if
        end do
        dst(y * image_w + x + 1) = real(accum, real32)
      end do
    end do
  end subroutine convolution_row_host

  subroutine convolution_column_host(dst, src, kernel, image_w, image_h)
    real(real32), intent(out) :: dst(:)
    real(real32), intent(in) :: src(:), kernel(:)
    integer, intent(in) :: image_w, image_h
    integer :: x, y, k, d
    real(real64) :: accum

    do y = 0, image_h - 1
      do x = 0, image_w - 1
        accum = 0.0_real64
        do k = -kernel_radius, kernel_radius
          d = y + k
          if (d >= 0 .and. d < image_h) then
            accum = accum + real(src(d * image_w + x + 1), real64) * real(kernel(kernel_radius - k + 1), real64)
          end if
        end do
        dst(y * image_w + x + 1) = real(accum, real32)
      end do
    end do
  end subroutine convolution_column_host

end program main
