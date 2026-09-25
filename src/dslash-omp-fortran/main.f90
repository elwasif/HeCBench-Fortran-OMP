! SPDX-License-Identifier: CC0-1.0
program dslash_main
  use iso_c_binding, only: c_double, c_int
  use iso_fortran_env, only: int64, real64, output_unit
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: ldim = 32
  integer, parameter :: iterations = 100
  integer, parameter :: warmups = 1
  integer, parameter :: total_sites = ldim * ldim * ldim * ldim
  integer, parameter :: even_sites = total_sites / 2
  real(dp), parameter :: eps = 2.0e-6_dp

  type, bind(C) :: dcomplex
    real(c_double) :: real
    real(c_double) :: imag
  end type dcomplex

  type, bind(C) :: su3_vector
    type(dcomplex) :: c(3)
  end type su3_vector

  type, bind(C) :: su3_matrix
    type(dcomplex) :: e(3,3)
  end type su3_matrix

  interface
    subroutine dslash_rng_seed_random_device() bind(C)
    end subroutine dslash_rng_seed_random_device

    function dslash_rng_uniform_real() bind(C) result(value)
      import :: c_double
      real(c_double) :: value
    end function dslash_rng_uniform_real

    function dslash_maxrss_mb() bind(C) result(value)
      import :: c_double
      real(c_double) :: value
    end function dslash_maxrss_mb
  end interface

  integer :: workgroup_size, arg_status
  character(len=64) :: arg
  type(su3_vector), allocatable :: src(:), dst(:), chkdst(:)
  type(su3_matrix), allocatable :: fat(:,:), lng(:,:)
  type(su3_matrix), allocatable :: fatbck(:,:), lngbck(:,:)
  integer(int64), allocatable :: fwd(:,:), bck(:,:), fwd3(:,:), bck3(:,:)
  real(dp) :: total_time, tflop, memory_usage, memory_allocated
  real(dp) :: max_error
  character(len=32) :: number_text

  if (command_argument_count() < 1) then
    write(*,'(A)') "Usage <workgroup size>"
    stop 1
  end if

  call get_command_argument(1, arg, status=arg_status)
  if (arg_status /= 0) stop 1
  read(arg, *) workgroup_size

  allocate(src(total_sites), dst(total_sites), chkdst(total_sites))
  allocate(fat(4,total_sites), lng(4,total_sites))
  allocate(fatbck(4,total_sites), lngbck(4,total_sites))
  allocate(fwd(4,total_sites), bck(4,total_sites), fwd3(4,total_sites), bck3(4,total_sites))

  call zero_vectors(dst)
  call zero_vectors(chkdst)
  call zero_matrices(fatbck)
  call zero_matrices(lngbck)

  call set_neighbors(fwd, bck, fwd3, bck3)
  call make_data(src, fat, lng)

  write(*,'(A,I0,A)') "Number of sites = ", ldim, "^4"
  write(*,'(A,I0,A,I0,A)') "Executing ", iterations, " iterations with ", warmups, " warmups"
  if (workgroup_size /= 0) write(*,'(A,I0)') "Threads per group = ", workgroup_size
  flush(output_unit)

  total_time = dslash_fn(src, dst, fat, lng, fatbck, lngbck, &
                         fwd, bck, fwd3, bck3, workgroup_size)
  write(number_text,'(F12.6)') total_time
  write(*,'(A,A,A)') "Total execution time = ", trim(adjustl(number_text)), " secs"

  write(*,'(A)') "Validating the result"
  flush(output_unit)
  call dslash_host(src, chkdst, fat, lng, fatbck, lngbck, fwd, bck, fwd3, bck3)
  max_error = max_vector_error(dst, chkdst)
  if (max_error >= eps) then
    write(*,'(A,ES15.8)') "Validation failed, max error = ", max_error
    stop 2
  end if

  tflop = real(iterations, dp) * real(even_sites, dp) * 1182.0_dp
  write(*,'(A,F0.3)') "Total GFLOP/s = ", tflop / total_time / 1.0e9_dp

  memory_usage = real(even_sites, dp) * &
      (real(144 * 4 * 4, dp) + real(48 * 16, dp) + real(8 * 16, dp) + real(48, dp))
  write(*,'(A,F0.3)') "Total GByte/s (GPU memory) = ", &
      real(iterations, dp) * memory_usage / total_time / 1.0e9_dp

  memory_allocated = real(total_sites, dp) * &
      (real(144 * 4 * 4, dp) + real(48 * 2, dp) + real(8 * 4 * 4, dp))
  write(*,'(A,I0)') "Total allocation for matrices = ", nint(memory_allocated / 1048576.0_dp)
  write(*,'(A,F0.2)') "Approximate memory usage = ", dslash_maxrss_mb()

contains

  integer function node_index(x, y, z, t) result(idx)
    integer, intent(in) :: x, y, z, t
    integer :: xr, yr, zr, tr, linear

    xr = modulo(x + ldim, ldim)
    yr = modulo(y + ldim, ldim)
    zr = modulo(z + ldim, ldim)
    tr = modulo(t + ldim, ldim)
    linear = xr + ldim * (yr + ldim * (zr + ldim * tr))
    if (modulo(x + y + z + t, 2) == 0) then
      idx = linear / 2 + 1
    else
      idx = (linear + total_sites) / 2 + 1
    end if
  end function node_index

  subroutine set_neighbors(fwd, bck, fwd3, bck3)
    integer(int64), intent(out) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(out) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer :: x, y, z, t, idx

    do t = 0, ldim - 1
      do z = 0, ldim - 1
        do y = 0, ldim - 1
          do x = 0, ldim - 1
            idx = node_index(x, y, z, t)
            fwd(1,idx) = node_index(x + 1, y, z, t)
            bck(1,idx) = node_index(x - 1, y, z, t)
            fwd(2,idx) = node_index(x, y + 1, z, t)
            bck(2,idx) = node_index(x, y - 1, z, t)
            fwd(3,idx) = node_index(x, y, z + 1, t)
            bck(3,idx) = node_index(x, y, z - 1, t)
            fwd(4,idx) = node_index(x, y, z, t + 1)
            bck(4,idx) = node_index(x, y, z, t - 1)
            fwd3(1,idx) = node_index(x + 3, y, z, t)
            bck3(1,idx) = node_index(x - 3, y, z, t)
            fwd3(2,idx) = node_index(x, y + 3, z, t)
            bck3(2,idx) = node_index(x, y - 3, z, t)
            fwd3(3,idx) = node_index(x, y, z + 3, t)
            bck3(3,idx) = node_index(x, y, z - 3, t)
            fwd3(4,idx) = node_index(x, y, z, t + 3)
            bck3(4,idx) = node_index(x, y, z, t - 3)
          end do
        end do
      end do
    end do
  end subroutine set_neighbors

  subroutine zero_vectors(vec)
    type(su3_vector), intent(inout) :: vec(:)
    integer :: site, row

    do site = 1, size(vec)
      do row = 1, 3
        vec(site)%c(row)%real = 0.0_dp
        vec(site)%c(row)%imag = 0.0_dp
      end do
    end do
  end subroutine zero_vectors

  subroutine zero_matrices(mat)
    type(su3_matrix), intent(inout) :: mat(:,:)
    integer :: site, dir, row, col

    do site = 1, size(mat, 2)
      do dir = 1, size(mat, 1)
        do row = 1, 3
          do col = 1, 3
            mat(dir,site)%e(row,col)%real = 0.0_dp
            mat(dir,site)%e(row,col)%imag = 0.0_dp
          end do
        end do
      end do
    end do
  end subroutine zero_matrices

  subroutine make_data(src, fat, lng)
    type(su3_vector), intent(out) :: src(total_sites)
    type(su3_matrix), intent(out) :: fat(4,total_sites)
    type(su3_matrix), intent(out) :: lng(4,total_sites)
    real(dp) :: r, i
    integer :: site, dir

    call dslash_rng_seed_random_device()
    do site = 1, total_sites
      r = dslash_rng_uniform_real()
      i = dslash_rng_uniform_real()
      src(site)%c(:)%real = r
      src(site)%c(:)%imag = i
      do dir = 1, 4
        r = dslash_rng_uniform_real()
        i = dslash_rng_uniform_real()
        fat(dir,site)%e(:,:)%real = r
        fat(dir,site)%e(:,:)%imag = i
        r = dslash_rng_uniform_real()
        i = dslash_rng_uniform_real()
        lng(dir,site)%e(:,:)%real = r
        lng(dir,site)%e(:,:)%imag = i
      end do
    end do
  end subroutine make_data

  real(dp) function dslash_fn(src, dst, fat, lng, fatbck, lngbck, &
                              fwd, bck, fwd3, bck3, wgsize) result(ttotal)
    type(su3_vector), intent(in) :: src(total_sites)
    type(su3_vector), intent(inout) :: dst(total_sites)
    type(su3_matrix), intent(in) :: fat(4,total_sites), lng(4,total_sites)
    type(su3_matrix), intent(inout) :: fatbck(4,total_sites), lngbck(4,total_sites)
    integer(int64), intent(in) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(in) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer, intent(in) :: wgsize

    integer :: my_site, dir, row, col, k, iter
    integer(int64) :: src_site
    real(dp) :: start_time, stop_time
    real(dp) :: acc_r, acc_i
    type(dcomplex) :: v(3)

    !$omp target data map(to: src, fat, lng, fwd, bck, fwd3, bck3) &
    !$omp& map(tofrom: dst, fatbck, lngbck)

    !$omp target teams distribute parallel do thread_limit(1) &
    !$omp& private(my_site, dir, row, col, src_site)
    do my_site = 1, even_sites
      do dir = 1, 4
        src_site = bck(dir,my_site)
        do row = 1, 3
          do col = 1, 3
            fatbck(dir,my_site)%e(row,col)%real = fat(dir,src_site)%e(col,row)%real
            fatbck(dir,my_site)%e(row,col)%imag = -fat(dir,src_site)%e(col,row)%imag
          end do
        end do
        src_site = bck3(dir,my_site)
        do row = 1, 3
          do col = 1, 3
            lngbck(dir,my_site)%e(row,col)%real = lng(dir,src_site)%e(col,row)%real
            lngbck(dir,my_site)%e(row,col)%imag = -lng(dir,src_site)%e(col,row)%imag
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

    write(*,'(A)') "Running dslash loop"
    write(*,'(A,I0)') "Setting number of work items to ", even_sites
    write(*,'(A,I0)') "Setting workgroup size to ", wgsize
    flush(output_unit)

    start_time = omp_get_wtime()
    do iter = 0, iterations + warmups - 1
      if (iter == warmups) start_time = omp_get_wtime()
      !$omp target teams distribute parallel do thread_limit(wgsize) &
      !$omp& private(my_site, dir, row, col, k, src_site, acc_r, acc_i, v)
      do my_site = 1, even_sites
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          src_site = fwd(1,my_site)
          do col = 1, 3
            acc_r = acc_r + fat(1,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                            fat(1,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
            acc_i = acc_i + fat(1,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                            fat(1,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
          end do
          dst(my_site)%c(row)%real = acc_r
          dst(my_site)%c(row)%imag = acc_i
        end do
        do dir = 2, 4
          src_site = fwd(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + fat(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                              fat(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
              acc_i = acc_i + fat(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                              fat(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
            end do
            dst(my_site)%c(row)%real = dst(my_site)%c(row)%real + acc_r
            dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag + acc_i
          end do
        end do

        do row = 1, 3
          v(row)%real = 0.0_dp
          v(row)%imag = 0.0_dp
        end do
        do dir = 1, 4
          src_site = fwd3(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + lng(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                              lng(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
              acc_i = acc_i + lng(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                              lng(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
            end do
            v(row)%real = v(row)%real + acc_r
            v(row)%imag = v(row)%imag + acc_i
          end do
        end do
        do row = 1, 3
          dst(my_site)%c(row)%real = dst(my_site)%c(row)%real + v(row)%real
          dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag + v(row)%imag
        end do

        do row = 1, 3
          v(row)%real = 0.0_dp
          v(row)%imag = 0.0_dp
        end do
        do dir = 1, 4
          src_site = bck(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + fatbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                              fatbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
              acc_i = acc_i + fatbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                              fatbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
            end do
            v(row)%real = v(row)%real + acc_r
            v(row)%imag = v(row)%imag + acc_i
          end do
        end do
        do row = 1, 3
          dst(my_site)%c(row)%real = dst(my_site)%c(row)%real - v(row)%real
          dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag - v(row)%imag
        end do

        do row = 1, 3
          v(row)%real = 0.0_dp
          v(row)%imag = 0.0_dp
        end do
        do dir = 1, 4
          src_site = bck3(dir,my_site)
          do row = 1, 3
            acc_r = 0.0_dp
            acc_i = 0.0_dp
            do col = 1, 3
              acc_r = acc_r + lngbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                              lngbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
              acc_i = acc_i + lngbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                              lngbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
            end do
            v(row)%real = v(row)%real + acc_r
            v(row)%imag = v(row)%imag + acc_i
          end do
        end do
        do row = 1, 3
          dst(my_site)%c(row)%real = dst(my_site)%c(row)%real - v(row)%real
          dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag - v(row)%imag
        end do
      end do
      !$omp end target teams distribute parallel do
    end do
    stop_time = omp_get_wtime()

    !$omp end target data
    ttotal = stop_time - start_time
  end function dslash_fn

  subroutine dslash_host(src, dst, fat, lng, fatbck, lngbck, fwd, bck, fwd3, bck3)
    type(su3_vector), intent(in) :: src(total_sites)
    type(su3_vector), intent(out) :: dst(total_sites)
    type(su3_matrix), intent(in) :: fat(4,total_sites), lng(4,total_sites)
    type(su3_matrix), intent(in) :: fatbck(4,total_sites), lngbck(4,total_sites)
    integer(int64), intent(in) :: fwd(4,total_sites), bck(4,total_sites)
    integer(int64), intent(in) :: fwd3(4,total_sites), bck3(4,total_sites)
    integer :: my_site, dir, row, col
    integer(int64) :: src_site
    real(dp) :: acc_r, acc_i
    type(dcomplex) :: v(3)

    call zero_vectors(dst)
    !$omp parallel do private(my_site, dir, row, col, src_site, acc_r, acc_i, v)
    do my_site = 1, even_sites
      do row = 1, 3
        acc_r = 0.0_dp
        acc_i = 0.0_dp
        src_site = fwd(1,my_site)
        do col = 1, 3
          acc_r = acc_r + fat(1,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                          fat(1,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
          acc_i = acc_i + fat(1,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                          fat(1,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
        end do
        dst(my_site)%c(row)%real = acc_r
        dst(my_site)%c(row)%imag = acc_i
      end do
      do dir = 2, 4
        src_site = fwd(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + fat(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                            fat(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
            acc_i = acc_i + fat(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                            fat(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
          end do
          dst(my_site)%c(row)%real = dst(my_site)%c(row)%real + acc_r
          dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag + acc_i
        end do
      end do

      v(:)%real = 0.0_dp
      v(:)%imag = 0.0_dp
      do dir = 1, 4
        src_site = fwd3(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + lng(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                            lng(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
            acc_i = acc_i + lng(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                            lng(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
          end do
          v(row)%real = v(row)%real + acc_r
          v(row)%imag = v(row)%imag + acc_i
        end do
      end do
      do row = 1, 3
        dst(my_site)%c(row)%real = dst(my_site)%c(row)%real + v(row)%real
        dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag + v(row)%imag
      end do

      v(:)%real = 0.0_dp
      v(:)%imag = 0.0_dp
      do dir = 1, 4
        src_site = bck(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + fatbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                            fatbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
            acc_i = acc_i + fatbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                            fatbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
          end do
          v(row)%real = v(row)%real + acc_r
          v(row)%imag = v(row)%imag + acc_i
        end do
      end do
      do row = 1, 3
        dst(my_site)%c(row)%real = dst(my_site)%c(row)%real - v(row)%real
        dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag - v(row)%imag
      end do

      v(:)%real = 0.0_dp
      v(:)%imag = 0.0_dp
      do dir = 1, 4
        src_site = bck3(dir,my_site)
        do row = 1, 3
          acc_r = 0.0_dp
          acc_i = 0.0_dp
          do col = 1, 3
            acc_r = acc_r + lngbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%real - &
                            lngbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%imag
            acc_i = acc_i + lngbck(dir,my_site)%e(row,col)%real * src(src_site)%c(col)%imag + &
                            lngbck(dir,my_site)%e(row,col)%imag * src(src_site)%c(col)%real
          end do
          v(row)%real = v(row)%real + acc_r
          v(row)%imag = v(row)%imag + acc_i
        end do
      end do
      do row = 1, 3
        dst(my_site)%c(row)%real = dst(my_site)%c(row)%real - v(row)%real
        dst(my_site)%c(row)%imag = dst(my_site)%c(row)%imag - v(row)%imag
      end do
    end do
    !$omp end parallel do
  end subroutine dslash_host

  real(dp) function max_vector_error(lhs, rhs) result(max_error)
    type(su3_vector), intent(in) :: lhs(total_sites), rhs(total_sites)
    integer :: site, row

    max_error = 0.0_dp
    do site = 1, even_sites
      do row = 1, 3
        max_error = max(max_error, abs(lhs(site)%c(row)%real - rhs(site)%c(row)%real))
        max_error = max(max_error, abs(lhs(site)%c(row)%imag - rhs(site)%c(row)%imag))
      end do
    end do
  end function max_vector_error

end program dslash_main
