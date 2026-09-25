! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: n_directions = 32
  integer, parameter :: max_m = 17
  integer, parameter :: max_dimensions = 10200
  real(real32), parameter :: l1error_tolerance = 1.0e-6_real32
  real(real32), parameter :: k_2powneg32 = 2.3283064e-10_real32

  character(len=256) :: arg0
  integer :: n_vectors, n_dimensions, repeat_count
  integer(int32), allocatable :: directions(:)
  real(real32), allocatable :: output_cpu(:), output_gpu(:)
  real(real32) :: l1norm_diff, l1norm_ref, l1error, ref
  real(real64) :: ktime
  integer :: d, v

  if (command_argument_count() /= 3) then
    call get_command_argument(0, arg0)
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <number of vectors> <number of dimensions> <repeat>'
    stop 1
  end if

  n_vectors = read_int_arg(1)
  n_dimensions = read_int_arg(2)
  repeat_count = read_int_arg(3)

  if (n_vectors < 1 .or. n_dimensions < 1 .or. repeat_count < 1) then
    write(*,'(A)') 'Error: arguments must be positive'
    stop 1
  end if

  if (n_dimensions > max_dimensions) then
    write(*,'(A,I0,A,I0)') 'Error: requested dimensions ', n_dimensions, ' exceed table size ', max_dimensions
    stop 1
  end if

  write(*,'(A)') 'Allocating CPU memory...'
  allocate(directions(0:n_dimensions * n_directions - 1))
  allocate(output_cpu(0:n_vectors * n_dimensions - 1))
  allocate(output_gpu(0:n_vectors * n_dimensions - 1))

  write(*,'(A)') 'Initializing direction numbers...'
  call init_sobol_direction_vectors(n_dimensions, directions)

  write(*,'(A)') 'Executing QRNG on GPU...'
  !$omp target data map(to: directions(0:n_dimensions * n_directions - 1)) &
  !$omp& map(from: output_gpu(0:n_dimensions * n_vectors - 1))
    ktime = sobol_gpu(repeat_count, n_vectors, n_dimensions, directions, output_gpu)
    write(*,'(A,ES12.5,A)') 'Average kernel execution time: ', (ktime * 1.0e-9_real64) / real(repeat_count, real64), ' (s)'
  !$omp end target data

  write(*,*)
  write(*,'(A)') 'Executing QRNG on CPU...'
  call sobol_cpu(n_vectors, n_dimensions, directions, output_cpu)

  write(*,'(A)') 'Checking results...'
  l1norm_diff = 0.0_real32
  l1norm_ref = 0.0_real32

  if (n_vectors == 1) then
    do d = 0, n_dimensions - 1
      v = 0
      ref = output_cpu(d * n_vectors + v)
      l1norm_diff = l1norm_diff + abs(output_gpu(d * n_vectors + v) - ref)
      l1norm_ref = l1norm_ref + abs(ref)
    end do
    l1error = l1norm_diff
    if (l1norm_ref /= 0.0_real32) then
      write(*,'(A)') 'Error: L1-Norm of the reference is not zero (for single vector), golden generator appears broken'
    else
      write(*,'(A,ES12.5)') 'L1-Error: ', l1error
    end if
  else
    do d = 0, n_dimensions - 1
      do v = 0, n_vectors - 1
        ref = output_cpu(d * n_vectors + v)
        l1norm_diff = l1norm_diff + abs(output_gpu(d * n_vectors + v) - ref)
        l1norm_ref = l1norm_ref + abs(ref)
      end do
    end do
    if (l1norm_ref == 0.0_real32) then
      l1error = huge(l1error)
      write(*,'(A)') 'Error: L1-Norm of the reference is zero, golden generator appears broken'
    else
      l1error = l1norm_diff / l1norm_ref
      write(*,'(A,ES12.5)') 'L1-Error: ', l1error
    end if
  end if

  write(*,'(A)') 'Shutting down...'
  if (l1error < l1error_tolerance) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(output_gpu, output_cpu, directions)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  integer function ffs32(x) result(pos)
    !$omp declare target
    integer(int32), intent(in) :: x
    integer :: bit

    pos = 0
    do bit = 0, 31
      if (iand(shiftr(x, bit), 1_int32) /= 0_int32) then
        pos = bit + 1
        return
      end if
    end do
  end function ffs32

  real(real32) function uint32_to_real32(x) result(value)
    !$omp declare target
    integer(int32), intent(in) :: x
    integer(int64), parameter :: mask32 = int(z'00000000ffffffff', int64)

    value = real(iand(int(x, int64), mask32), real32)
  end function uint32_to_real32

  subroutine init_sobol_direction_vectors(n_dimensions, directions)
    integer, intent(in) :: n_dimensions
    integer(int32), intent(out) :: directions(0:)
    integer(int32), allocatable :: degree(:), poly_a(:), primitive_m(:, :)
    integer :: dim, i, j, d, base

    allocate(degree(0:max_dimensions - 1))
    allocate(poly_a(0:max_dimensions - 1))
    allocate(primitive_m(0:max_m - 1, 0:max_dimensions - 1))
    include 'sobol_primitives.inc'

    do dim = 0, n_dimensions - 1
      base = dim * n_directions
      if (dim == 0) then
        do i = 0, n_directions - 1
          directions(base + i) = shiftl(1_int32, 31 - i)
        end do
      else
        d = degree(dim)
        do i = 0, d - 1
          directions(base + i) = shiftl(primitive_m(i, dim), 31 - i)
        end do
        do i = d, n_directions - 1
          directions(base + i) = ieor(directions(base + i - d), shiftr(directions(base + i - d), d))
          do j = 1, d - 1
            if (iand(shiftr(poly_a(dim), d - 1 - j), 1_int32) /= 0_int32) then
              directions(base + i) = ieor(directions(base + i), directions(base + i - j))
            end if
          end do
        end do
      end if
    end do

    deallocate(primitive_m, poly_a, degree)
  end subroutine init_sobol_direction_vectors

  subroutine sobol_cpu(n_vectors, n_dimensions, directions, output)
    integer, intent(in) :: n_vectors, n_dimensions
    integer(int32), intent(in) :: directions(0:)
    real(real32), intent(out) :: output(0:)
    integer :: dim, i, base_dir, base_out, c
    integer(int32) :: x

    do dim = 0, n_dimensions - 1
      base_dir = dim * n_directions
      base_out = dim * n_vectors
      x = 0_int32
      output(base_out) = 0.0_real32
      do i = 1, n_vectors - 1
        c = ffs32(not(int(i - 1, int32))) - 1
        x = ieor(x, directions(base_dir + c))
        output(base_out + i) = uint32_to_real32(x) * k_2powneg32
      end do
    end do
  end subroutine sobol_cpu

  real(real64) function sobol_gpu(repeat_count, n_vectors, n_dimensions, directions, output) result(elapsed_ns)
    integer, intent(in) :: repeat_count, n_vectors, n_dimensions
    integer(int32), intent(in) :: directions(0:)
    real(real32), intent(out) :: output(0:)
    integer, parameter :: threadsperblock = 128
    integer :: dim_grid_x, dim_grid_y, target_dim_grid_x, num_team
    integer :: rep, team_x, team_y, tid_x, thread_size_x
    integer :: i0, stride, k, i, dir_base, out_base, ffs_stride, team_id
    integer(int32) :: scratch_v(0:n_directions - 1)
    integer(int32) :: g, x, mask, v_log2stridem1, v_stridemask
    real(real64) :: t0, t1

    dim_grid_y = n_dimensions
    if (n_dimensions < 4 * 24) then
      dim_grid_x = 4 * 24
    else
      dim_grid_x = 1
    end if

    if (dim_grid_x > (n_vectors / threadsperblock)) then
      dim_grid_x = (n_vectors + threadsperblock - 1) / threadsperblock
    end if

    target_dim_grid_x = dim_grid_x
    dim_grid_x = 1
    do while (dim_grid_x < target_dim_grid_x)
      dim_grid_x = dim_grid_x * 2
    end do
    if (dim_grid_x < 1) dim_grid_x = 1

    num_team = dim_grid_x * dim_grid_y
    t0 = omp_get_wtime()
    do rep = 1, repeat_count
      !$omp target teams num_teams(num_team) thread_limit(threadsperblock) private(scratch_v)
        !$omp parallel private(team_id, team_x, team_y, tid_x, thread_size_x, dir_base, out_base, i0, stride, g, x, mask, k, ffs_stride, v_log2stridem1, v_stridemask, i)
          team_id = omp_get_team_num()
          team_x = mod(team_id, dim_grid_x)
          team_y = team_id / dim_grid_x
          tid_x = omp_get_thread_num()
          thread_size_x = omp_get_num_threads()
          dir_base = n_directions * team_y
          out_base = n_vectors * team_y

          if (tid_x < n_directions) then
            scratch_v(tid_x) = directions(dir_base + tid_x)
          end if

          !$omp barrier

          i0 = team_x * thread_size_x + tid_x
          stride = dim_grid_x * thread_size_x
          g = ieor(int(i0, int32), shiftr(int(i0, int32), 1))
          x = 0_int32
          ffs_stride = ffs32(int(stride, int32))

          do k = 0, ffs_stride - 2
            if (iand(g, 1_int32) /= 0_int32) then
              mask = not(0_int32)
            else
              mask = 0_int32
            end if
            x = ieor(x, iand(mask, scratch_v(k)))
            g = shiftr(g, 1)
          end do

          if (i0 < n_vectors) then
            output(out_base + i0) = uint32_to_real32(x) * k_2powneg32
          end if

          v_log2stridem1 = scratch_v(ffs_stride - 2)
          v_stridemask = int(stride - 1, int32)
          i = i0 + stride
          do while (i < n_vectors)
            x = ieor(x, ieor(v_log2stridem1, scratch_v(ffs32(not(ior(int(i - stride, int32), v_stridemask))) - 1)))
            output(out_base + i) = uint32_to_real32(x) * k_2powneg32
            i = i + stride
          end do
        !$omp end parallel
      !$omp end target teams
    end do
    t1 = omp_get_wtime()
    elapsed_ns = (t1 - t0) * 1.0e9_real64
  end function sobol_gpu

end program main
