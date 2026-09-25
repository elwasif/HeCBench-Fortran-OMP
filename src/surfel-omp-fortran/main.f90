! SPDX-License-Identifier: CC0-1.0
module surfel_mod
  use iso_c_binding, only: c_float, c_int
  use iso_fortran_env, only: real32, real64
  use omp_lib
  implicit none

  integer, parameter :: col_p_x = 0, col_p_y = 1, col_p_z = 2
  integer, parameter :: col_n_x = 3, col_n_y = 4, col_n_z = 5
  integer, parameter :: col_rsq = 6, col_dim = 7

  interface
    subroutine c_fill_src(src, n) bind(C, name="surfel_fill_src")
      import :: c_float, c_int
      real(c_float), intent(out) :: src(*)
      integer(c_int), value :: n
    end subroutine c_fill_src
  end interface

contains

  subroutine fill_src(src, n)
    real(real32), intent(out) :: src(0:)
    integer, intent(in) :: n

    call c_fill_src(src, int(n, c_int))
  end subroutine fill_src

  subroutine surfel_render(src, n, f, w, h, dst)
    real(real32), intent(in) :: src(0:)
    integer, intent(in) :: n, w, h
    real(real32), intent(in) :: f
    real(real32), intent(out) :: dst(0:)
    integer :: idx, idy, i
    real(real32) :: ray0, ray1, ray2, pt0, pt1, pt2
    real(real32) :: p0, p1, p2, n0, n1, n2, rsq_max
    real(real32) :: p_dot_n, ds_dot_ray, alpha, t, rsq, d_min

    !$omp target teams distribute parallel do collapse(2) thread_limit(256) &
    !$omp& private(ray0, ray1, ray2, pt0, pt1, pt2, p0, p1, p2, n0, n1, n2, rsq_max, p_dot_n, ds_dot_ray, alpha, t, rsq, d_min, i)
    do idy = 0, h - 1
      do idx = 0, w - 1
        ray0 = real(idx, real32) - real(w - 1, real32) * 0.5_real32
        ray1 = real(idy, real32) - real(h - 1, real32) * 0.5_real32
        ray2 = f
        d_min = 1.0e20_real32

        do i = 0, n - 1
          p0 = src(i * col_dim + col_p_x)
          p1 = src(i * col_dim + col_p_y)
          p2 = src(i * col_dim + col_p_z)
          n0 = src(i * col_dim + col_n_x)
          n1 = src(i * col_dim + col_n_y)
          n2 = src(i * col_dim + col_n_z)
          rsq_max = src(i * col_dim + col_rsq)
          p_dot_n = p0 * n0 + p1 * n1 + p2 * n2
          ds_dot_ray = ray0 * n0 + ray1 * n1 + ray2 * n2
          alpha = p_dot_n / ds_dot_ray
          pt0 = ray0 * alpha - p0
          pt1 = ray1 * alpha - p1
          pt2 = ray2 * alpha - p2
          t = ray2 * alpha
          rsq = pt0 * pt0 + pt1 * pt1 + pt2 * pt2
          if (rsq < rsq_max .and. d_min > t) d_min = t
        end do
        if (d_min > 100.0_real32) then
          dst(idy * w + idx) = 0.0_real32
        else
          dst(idy * w + idx) = d_min
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine surfel_render

  subroutine reference(src, n, f, w, h, dst)
    real(real32), intent(in) :: src(0:)
    integer, intent(in) :: n, w, h
    real(real32), intent(in) :: f
    real(real32), intent(out) :: dst(0:)
    integer :: idx, idy, i
    real(real32) :: ray0, ray1, ray2, pt0, pt1, pt2
    real(real32) :: p0, p1, p2, n0, n1, n2, rsq_max
    real(real32) :: p_dot_n, ds_dot_ray, alpha, t, rsq, d_min

    do idy = 0, h - 1
      do idx = 0, w - 1
        ray0 = real(idx, real32) - real(w - 1, real32) * 0.5_real32
        ray1 = real(idy, real32) - real(h - 1, real32) * 0.5_real32
        ray2 = f
        d_min = 1.0e20_real32
        do i = 0, n - 1
          p0 = src(i * col_dim + col_p_x)
          p1 = src(i * col_dim + col_p_y)
          p2 = src(i * col_dim + col_p_z)
          n0 = src(i * col_dim + col_n_x)
          n1 = src(i * col_dim + col_n_y)
          n2 = src(i * col_dim + col_n_z)
          rsq_max = src(i * col_dim + col_rsq)
          p_dot_n = p0 * n0 + p1 * n1 + p2 * n2
          ds_dot_ray = ray0 * n0 + ray1 * n1 + ray2 * n2
          alpha = p_dot_n / ds_dot_ray
          pt0 = ray0 * alpha - p0
          pt1 = ray1 * alpha - p1
          pt2 = ray2 * alpha - p2
          t = ray2 * alpha
          rsq = pt0 * pt0 + pt1 * pt1 + pt2 * pt2
          if (rsq < rsq_max .and. d_min > t) d_min = t
        end do
        dst(idy * w + idx) = merge(0.0_real32, d_min, d_min > 100.0_real32)
      end do
    end do
  end subroutine reference

end module surfel_mod

program main
  use iso_fortran_env, only: real32, real64
  use omp_lib
  use surfel_mod
  implicit none

  integer :: argc, n, w, h, repeat, status, f_idx, iter, i
  integer :: src_size, dst_size
  real(real32), allocatable :: src(:), h_dst(:), r_dst(:)
  real(real32), parameter :: inverse_focal_length(0:2) = [0.005_real32, 0.02_real32, 0.036_real32]
  real(real64) :: start_time, elapsed_ns
  logical :: ok
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 4) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') 'Usage: ', trim(prog), ' <number of surfels> <output width> <output height> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) n
  if (status /= 0) stop 1
  call get_command_argument(2, arg)
  read(arg, *, iostat=status) w
  if (status /= 0) stop 1
  call get_command_argument(3, arg)
  read(arg, *, iostat=status) h
  if (status /= 0) stop 1
  call get_command_argument(4, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  src_size = n * col_dim
  dst_size = w * h
  allocate(src(0:src_size - 1), h_dst(0:dst_size - 1), r_dst(0:dst_size - 1))
  call fill_src(src, n)
  ok = .true.

  write(*,'(A)') '-------------------------------------'
  write(*,'(A)') ' surfelRenderTest with type float32  '
  write(*,'(A)') '-------------------------------------'

  !$omp target data map(to: src(0:src_size - 1)) map(alloc: h_dst(0:dst_size - 1))
  do f_idx = 0, 2
    write(*,'(/,A,I0)') 'f = ', f_idx
    call reference(src, n, inverse_focal_length(f_idx), w, h, r_dst)
    start_time = omp_get_wtime()
    do iter = 1, repeat
      call surfel_render(src, n, inverse_focal_length(f_idx), w, h, h_dst)
    end do
    elapsed_ns = (omp_get_wtime() - start_time) * 1.0e9_real64
    write(*,'(A,F8.6,A)') 'Average kernel execution time: ', elapsed_ns * 1.0e-6_real64 / real(repeat, real64), ' (ms)'
    !$omp target update from(h_dst(0:dst_size - 1))
    do i = 0, dst_size - 1
      if (abs(h_dst(i) - r_dst(i)) > 1.0e-3_real32) then
        write(*,'(F0.6,1X,F0.6)') h_dst(i), r_dst(i)
        ok = .false.
        exit
      end if
    end do
    if (.not. ok) exit
  end do
  !$omp end target data

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(src, h_dst, r_dst)
end program main
