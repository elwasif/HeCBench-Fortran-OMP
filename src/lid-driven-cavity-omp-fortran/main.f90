! SPDX-License-Identifier: CC0-1.0
module lid_driven_cavity_mod
  use omp_lib
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: num = 512
  integer, parameter :: num_2 = num / 2
  integer, parameter :: block_size = 128
  integer, parameter :: full_stride = num + 2
  integer, parameter :: pres_stride = num_2 + 2
  integer, parameter :: size_grid = full_stride * full_stride
  integer, parameter :: size_pres = pres_stride * full_stride
  integer, parameter :: size_res = num / (2 * block_size) * num
  integer, parameter :: num_pressure_teams = (num * num / 2) / block_size

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: two = 2.0_dp
  real(dp), parameter :: four = 4.0_dp
  real(dp), parameter :: small = 1.0e-10_dp
  real(dp), parameter :: re_num = 1000.0_dp
  real(dp), parameter :: omega = 1.7_dp
  real(dp), parameter :: mix_param = 0.9_dp
  real(dp), parameter :: tau = 0.5_dp
  real(dp), parameter :: gx = 0.0_dp
  real(dp), parameter :: gy = 0.0_dp
  real(dp), parameter :: x_length = 1.0_dp
  real(dp), parameter :: y_length = 1.0_dp
  real(dp), parameter :: dx = x_length / real(num, dp)
  real(dp), parameter :: dy = y_length / real(num, dp)

contains

  pure integer function grid_idx(col, row) result(idx)
    integer, intent(in) :: col, row
    idx = col * full_stride + row
  end function grid_idx

  pure integer function pres_idx(col, row) result(idx)
    integer, intent(in) :: col, row
    idx = col * pres_stride + row
  end function pres_idx

  subroutine set_bcs_host(u, v)
    real(dp), intent(inout) :: u(0:), v(0:)
    integer :: ind

    do ind = 0, num + 1
      u(grid_idx(0, ind)) = zero
      v(grid_idx(0, ind)) = -v(grid_idx(1, ind))

      u(grid_idx(num, ind)) = zero
      v(grid_idx(num + 1, ind)) = -v(grid_idx(num, ind))

      u(grid_idx(ind, 0)) = -u(grid_idx(ind, 1))
      v(grid_idx(ind, 0)) = zero

      u(grid_idx(ind, num + 1)) = two - u(grid_idx(ind, num))
      v(grid_idx(ind, num)) = zero

      if (ind == num) then
        u(grid_idx(0, 0)) = zero
        v(grid_idx(0, 0)) = -v(grid_idx(1, 0))
        u(grid_idx(0, num + 1)) = zero
        v(grid_idx(0, num + 1)) = -v(grid_idx(1, num + 1))

        u(grid_idx(num, 0)) = zero
        v(grid_idx(num + 1, 0)) = -v(grid_idx(num, 0))
        u(grid_idx(num, num + 1)) = zero
        v(grid_idx(num + 1, num + 1)) = -v(grid_idx(num, num + 1))

        u(grid_idx(0, 0)) = -u(grid_idx(0, 1))
        v(grid_idx(0, 0)) = zero
        u(grid_idx(num + 1, 0)) = -u(grid_idx(num + 1, 1))
        v(grid_idx(num + 1, 0)) = zero

        u(grid_idx(0, num + 1)) = two - u(grid_idx(0, num))
        v(grid_idx(0, num)) = zero
        u(grid_idx(num + 1, num + 1)) = two - u(grid_idx(num + 1, num))
        v(grid_idx(ind, num + 1)) = zero
      end if
    end do
  end subroutine set_bcs_host

  subroutine calculate_f_kernel(dt, u, v, f)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u(0:), v(0:)
    real(dp), intent(inout) :: f(0:)
    integer :: col, row
    real(dp) :: u_ij, u_ip1j, u_ijp1, u_im1j, u_ijm1
    real(dp) :: v_ij, v_ip1j, v_ijm1, v_ip1jm1
    real(dp) :: du2dx, duvdy, d2udx2, d2udy2

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(col, row, u_ij, u_ip1j, u_ijp1, u_im1j, u_ijm1, v_ij, v_ip1j, &
    !$omp& v_ijm1, v_ip1jm1, du2dx, duvdy, d2udx2, d2udy2)
    do col = 1, num
      do row = 1, num
        if (col == num) then
          f(grid_idx(0, row)) = u(grid_idx(0, row))
          f(grid_idx(num, row)) = u(grid_idx(num, row))
        else
          u_ij = u(grid_idx(col, row))
          u_ip1j = u(grid_idx(col + 1, row))
          u_ijp1 = u(grid_idx(col, row + 1))
          u_im1j = u(grid_idx(col - 1, row))
          u_ijm1 = u(grid_idx(col, row - 1))

          v_ij = v(grid_idx(col, row))
          v_ip1j = v(grid_idx(col + 1, row))
          v_ijm1 = v(grid_idx(col, row - 1))
          v_ip1jm1 = v(grid_idx(col + 1, row - 1))

          du2dx = (((u_ij + u_ip1j) * (u_ij + u_ip1j) - (u_im1j + u_ij) * (u_im1j + u_ij)) + &
              mix_param * (abs(u_ij + u_ip1j) * (u_ij - u_ip1j) - &
              abs(u_im1j + u_ij) * (u_im1j - u_ij))) / (four * dx)
          duvdy = ((v_ij + v_ip1j) * (u_ij + u_ijp1) - (v_ijm1 + v_ip1jm1) * (u_ijm1 + u_ij) + &
              mix_param * (abs(v_ij + v_ip1j) * (u_ij - u_ijp1) - &
              abs(v_ijm1 + v_ip1jm1) * (u_ijm1 - u_ij))) / (four * dy)
          d2udx2 = (u_ip1j - (two * u_ij) + u_im1j) / (dx * dx)
          d2udy2 = (u_ijp1 - (two * u_ij) + u_ijm1) / (dy * dy)

          f(grid_idx(col, row)) = u_ij + dt * (((d2udx2 + d2udy2) / re_num) - du2dx - duvdy + gx)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine calculate_f_kernel

  subroutine calculate_g_kernel(dt, u, v, g)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u(0:), v(0:)
    real(dp), intent(inout) :: g(0:)
    integer :: col, row
    real(dp) :: u_ij, u_ijp1, u_im1j, u_im1jp1
    real(dp) :: v_ij, v_ijp1, v_ip1j, v_ijm1, v_im1j
    real(dp) :: dv2dy, duvdx, d2vdx2, d2vdy2

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(col, row, u_ij, u_ijp1, u_im1j, u_im1jp1, v_ij, v_ijp1, v_ip1j, &
    !$omp& v_ijm1, v_im1j, dv2dy, duvdx, d2vdx2, d2vdy2)
    do col = 1, num
      do row = 1, num
        if (row == num) then
          g(grid_idx(col, 0)) = v(grid_idx(col, 0))
          g(grid_idx(col, num)) = v(grid_idx(col, num))
        else
          u_ij = u(grid_idx(col, row))
          u_ijp1 = u(grid_idx(col, row + 1))
          u_im1j = u(grid_idx(col - 1, row))
          u_im1jp1 = u(grid_idx(col - 1, row + 1))

          v_ij = v(grid_idx(col, row))
          v_ijp1 = v(grid_idx(col, row + 1))
          v_ip1j = v(grid_idx(col + 1, row))
          v_ijm1 = v(grid_idx(col, row - 1))
          v_im1j = v(grid_idx(col - 1, row))

          dv2dy = ((v_ij + v_ijp1) * (v_ij + v_ijp1) - (v_ijm1 + v_ij) * (v_ijm1 + v_ij) + &
              mix_param * (abs(v_ij + v_ijp1) * (v_ij - v_ijp1) - &
              abs(v_ijm1 + v_ij) * (v_ijm1 - v_ij))) / (four * dy)
          duvdx = ((u_ij + u_ijp1) * (v_ij + v_ip1j) - (u_im1j + u_im1jp1) * (v_im1j + v_ij) + &
              mix_param * (abs(u_ij + u_ijp1) * (v_ij - v_ip1j) - &
              abs(u_im1j + u_im1jp1) * (v_im1j - v_ij))) / (four * dx)
          d2vdx2 = (v_ip1j - (two * v_ij) + v_im1j) / (dx * dx)
          d2vdy2 = (v_ijp1 - (two * v_ij) + v_ijm1) / (dy * dy)

          g(grid_idx(col, row)) = v_ij + dt * (((d2vdx2 + d2vdy2) / re_num) - dv2dy - duvdx + gy)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine calculate_g_kernel

  subroutine sum_pressure_kernel(pres_red, pres_black, pres_sum)
    real(dp), intent(in) :: pres_red(0:), pres_black(0:)
    real(dp), intent(inout) :: pres_sum(0:)
    integer :: lid, tid, gid, row, col, i
    real(dp) :: sum_cache(0:block_size - 1)
    real(dp) :: pres_r, pres_b

    !$omp target teams num_teams(num_pressure_teams) thread_limit(block_size) private(sum_cache)
      !$omp parallel private(lid, tid, gid, row, col, pres_r, pres_b, i) shared(sum_cache)
        lid = omp_get_thread_num()
        tid = omp_get_team_num()
        gid = tid * block_size + lid
        row = mod(gid, num_2) + 1
        col = gid / num_2 + 1

        pres_r = pres_red(pres_idx(col, row))
        pres_b = pres_black(pres_idx(col, row))
        sum_cache(lid) = (pres_r * pres_r) + (pres_b * pres_b)

        !$omp barrier
        i = block_size / 2
        do while (i /= 0)
          if (lid < i) sum_cache(lid) = sum_cache(lid) + sum_cache(lid + i)
          !$omp barrier
          i = i / 2
        end do

        if (lid == 0) pres_sum(tid) = sum_cache(0)
      !$omp end parallel
    !$omp end target teams
  end subroutine sum_pressure_kernel

  subroutine set_horz_pres_bcs_kernel(pres_red, pres_black)
    real(dp), intent(inout) :: pres_red(0:), pres_black(0:)
    integer :: col, pcol

    !$omp target teams distribute parallel do thread_limit(block_size) private(col, pcol)
    do col = 1, num_2
      pcol = (col * 2) - 1
      pres_black(pres_idx(pcol, 0)) = pres_red(pres_idx(pcol, 1))
      pres_red(pres_idx(pcol + 1, 0)) = pres_black(pres_idx(pcol + 1, 1))
      pres_red(pres_idx(pcol, num_2 + 1)) = pres_black(pres_idx(pcol, num_2))
      pres_black(pres_idx(pcol + 1, num_2 + 1)) = pres_red(pres_idx(pcol + 1, num_2))
    end do
    !$omp end target teams distribute parallel do
  end subroutine set_horz_pres_bcs_kernel

  subroutine set_vert_pres_bcs_kernel(pres_red, pres_black)
    real(dp), intent(inout) :: pres_red(0:), pres_black(0:)
    integer :: row

    !$omp target teams distribute parallel do thread_limit(block_size) private(row)
    do row = 1, num_2
      pres_black(pres_idx(0, row)) = pres_red(pres_idx(1, row))
      pres_red(pres_idx(0, row)) = pres_black(pres_idx(1, row))
      pres_black(pres_idx(num + 1, row)) = pres_red(pres_idx(num, row))
      pres_red(pres_idx(num + 1, row)) = pres_black(pres_idx(num, row))
    end do
    !$omp end target teams distribute parallel do
  end subroutine set_vert_pres_bcs_kernel

  subroutine red_kernel(dt, f, g, pres_black, pres_red)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: f(0:), g(0:), pres_black(0:)
    real(dp), intent(inout) :: pres_red(0:)
    integer :: col, row
    real(dp) :: p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(col, row, p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs)
    do col = 1, num
      do row = 1, num_2
        p_ij = pres_red(pres_idx(col, row))
        p_im1j = pres_black(pres_idx(col - 1, row))
        p_ip1j = pres_black(pres_idx(col + 1, row))
        p_ijm1 = pres_black(pres_idx(col, row - iand(col, 1)))
        p_ijp1 = pres_black(pres_idx(col, row + iand(col + 1, 1)))

        rhs = (((f(grid_idx(col, (2 * row) - iand(col, 1))) - &
            f(grid_idx(col - 1, (2 * row) - iand(col, 1)))) / dx) + &
            ((g(grid_idx(col, (2 * row) - iand(col, 1))) - &
            g(grid_idx(col, (2 * row) - iand(col, 1) - 1))) / dy)) / dt

        pres_red(pres_idx(col, row)) = p_ij * (one - omega) + omega * &
            (((p_ip1j + p_im1j) / (dx * dx)) + ((p_ijp1 + p_ijm1) / (dy * dy)) - rhs) / &
            ((two / (dx * dx)) + (two / (dy * dy)))
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine red_kernel

  subroutine black_kernel(dt, f, g, pres_red, pres_black)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: f(0:), g(0:), pres_red(0:)
    real(dp), intent(inout) :: pres_black(0:)
    integer :: col, row
    real(dp) :: p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(col, row, p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs)
    do col = 1, num
      do row = 1, num_2
        p_ij = pres_black(pres_idx(col, row))
        p_im1j = pres_red(pres_idx(col - 1, row))
        p_ip1j = pres_red(pres_idx(col + 1, row))
        p_ijm1 = pres_red(pres_idx(col, row - iand(col + 1, 1)))
        p_ijp1 = pres_red(pres_idx(col, row + iand(col, 1)))

        rhs = (((f(grid_idx(col, (2 * row) - iand(col + 1, 1))) - &
            f(grid_idx(col - 1, (2 * row) - iand(col + 1, 1)))) / dx) + &
            ((g(grid_idx(col, (2 * row) - iand(col + 1, 1))) - &
            g(grid_idx(col, (2 * row) - iand(col + 1, 1) - 1))) / dy)) / dt

        pres_black(pres_idx(col, row)) = p_ij * (one - omega) + omega * &
            (((p_ip1j + p_im1j) / (dx * dx)) + ((p_ijp1 + p_ijm1) / (dy * dy)) - rhs) / &
            ((two / (dx * dx)) + (two / (dy * dy)))
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine black_kernel

  subroutine calc_residual_kernel(dt, f, g, pres_red, pres_black, res_arr)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: f(0:), g(0:), pres_red(0:), pres_black(0:)
    real(dp), intent(inout) :: res_arr(0:)
    integer :: lid, tid, gid, row, col, i
    real(dp) :: sum_cache(0:block_size - 1)
    real(dp) :: p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs, res, res2

    !$omp target teams num_teams(num_pressure_teams) thread_limit(block_size) private(sum_cache)
      !$omp parallel private(lid, tid, gid, row, col, i, p_ij, p_im1j, p_ip1j, p_ijm1, p_ijp1, rhs, res, res2) &
      !$omp& shared(sum_cache)
        lid = omp_get_thread_num()
        tid = omp_get_team_num()
        gid = tid * block_size + lid
        row = mod(gid, num_2) + 1
        col = gid / num_2 + 1

        p_ij = pres_red(pres_idx(col, row))
        p_im1j = pres_black(pres_idx(col - 1, row))
        p_ip1j = pres_black(pres_idx(col + 1, row))
        p_ijm1 = pres_black(pres_idx(col, row - iand(col, 1)))
        p_ijp1 = pres_black(pres_idx(col, row + iand(col + 1, 1)))
        rhs = (((f(grid_idx(col, (2 * row) - iand(col, 1))) - &
            f(grid_idx(col - 1, (2 * row) - iand(col, 1)))) / dx) + &
            ((g(grid_idx(col, (2 * row) - iand(col, 1))) - &
            g(grid_idx(col, (2 * row) - iand(col, 1) - 1))) / dy)) / dt
        res = ((p_ip1j - (two * p_ij) + p_im1j) / (dx * dx)) + &
            ((p_ijp1 - (two * p_ij) + p_ijm1) / (dy * dy)) - rhs

        p_ij = pres_black(pres_idx(col, row))
        p_im1j = pres_red(pres_idx(col - 1, row))
        p_ip1j = pres_red(pres_idx(col + 1, row))
        p_ijm1 = pres_red(pres_idx(col, row - iand(col + 1, 1)))
        p_ijp1 = pres_red(pres_idx(col, row + iand(col, 1)))
        rhs = (((f(grid_idx(col, (2 * row) - iand(col + 1, 1))) - &
            f(grid_idx(col - 1, (2 * row) - iand(col + 1, 1)))) / dx) + &
            ((g(grid_idx(col, (2 * row) - iand(col + 1, 1))) - &
            g(grid_idx(col, (2 * row) - iand(col + 1, 1) - 1))) / dy)) / dt
        res2 = ((p_ip1j - (two * p_ij) + p_im1j) / (dx * dx)) + &
            ((p_ijp1 - (two * p_ij) + p_ijm1) / (dy * dy)) - rhs

        sum_cache(lid) = (res * res) + (res2 * res2)
        !$omp barrier
        i = block_size / 2
        do while (i /= 0)
          if (lid < i) sum_cache(lid) = sum_cache(lid) + sum_cache(lid + i)
          !$omp barrier
          i = i / 2
        end do

        if (lid == 0) res_arr(tid) = sum_cache(0)
      !$omp end parallel
    !$omp end target teams
  end subroutine calc_residual_kernel

  subroutine calculate_u_kernel(dt, f, pres_red, pres_black, u, max_u_arr)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: f(0:), pres_red(0:), pres_black(0:)
    real(dp), intent(inout) :: u(0:), max_u_arr(0:)
    integer :: lid, tid, gid, row, col, i
    real(dp) :: max_cache(0:block_size - 1)
    real(dp) :: p_ij, p_ip1j, new_u, new_u2

    !$omp target teams num_teams(num_pressure_teams) thread_limit(block_size) private(max_cache)
      !$omp parallel private(lid, tid, gid, row, col, i, p_ij, p_ip1j, new_u, new_u2) shared(max_cache)
        lid = omp_get_thread_num()
        tid = omp_get_team_num()
        gid = tid * block_size + lid
        row = mod(gid, num_2) + 1
        col = gid / num_2 + 1

        max_cache(lid) = zero
        new_u = zero

        if (col /= num) then
          p_ij = pres_red(pres_idx(col, row))
          p_ip1j = pres_black(pres_idx(col + 1, row))
          new_u = f(grid_idx(col, (2 * row) - iand(col, 1))) - (dt * (p_ip1j - p_ij) / dx)
          u(grid_idx(col, (2 * row) - iand(col, 1))) = new_u

          p_ij = pres_black(pres_idx(col, row))
          p_ip1j = pres_red(pres_idx(col + 1, row))
          new_u2 = f(grid_idx(col, (2 * row) - iand(col + 1, 1))) - (dt * (p_ip1j - p_ij) / dx)
          u(grid_idx(col, (2 * row) - iand(col + 1, 1))) = new_u2

          new_u = max(abs(new_u), abs(new_u2))
          if ((2 * row) == num) new_u = max(new_u, abs(u(grid_idx(col, num + 1))))
        else
          new_u = max(abs(u(grid_idx(num, 2 * row))), abs(u(grid_idx(0, 2 * row))))
          new_u = max(abs(u(grid_idx(num, (2 * row) - 1))), new_u)
          new_u = max(abs(u(grid_idx(0, (2 * row) - 1))), new_u)
          new_u = max(abs(u(grid_idx(num + 1, 2 * row))), new_u)
          new_u = max(abs(u(grid_idx(num + 1, (2 * row) - 1))), new_u)
        end if

        max_cache(lid) = new_u
        !$omp barrier
        i = block_size / 2
        do while (i /= 0)
          if (lid < i) max_cache(lid) = max(max_cache(lid), max_cache(lid + i))
          !$omp barrier
          i = i / 2
        end do

        if (lid == 0) max_u_arr(tid) = max_cache(0)
      !$omp end parallel
    !$omp end target teams
  end subroutine calculate_u_kernel

  subroutine calculate_v_kernel(dt, g, pres_red, pres_black, v, max_v_arr)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: g(0:), pres_red(0:), pres_black(0:)
    real(dp), intent(inout) :: v(0:), max_v_arr(0:)
    integer :: lid, tid, gid, row, col, i
    real(dp) :: max_cache(0:block_size - 1)
    real(dp) :: p_ij, p_ijp1, new_v, new_v2

    !$omp target teams num_teams(num_pressure_teams) thread_limit(block_size) private(max_cache)
      !$omp parallel private(lid, tid, gid, row, col, i, p_ij, p_ijp1, new_v, new_v2) shared(max_cache)
        lid = omp_get_thread_num()
        tid = omp_get_team_num()
        gid = tid * block_size + lid
        row = mod(gid, num_2) + 1
        col = gid / num_2 + 1

        max_cache(lid) = zero
        new_v = zero

        if (row /= num_2) then
          p_ij = pres_red(pres_idx(col, row))
          p_ijp1 = pres_black(pres_idx(col, row + iand(col + 1, 1)))
          new_v = g(grid_idx(col, (2 * row) - iand(col, 1))) - (dt * (p_ijp1 - p_ij) / dy)
          v(grid_idx(col, (2 * row) - iand(col, 1))) = new_v

          p_ij = pres_black(pres_idx(col, row))
          p_ijp1 = pres_red(pres_idx(col, row + iand(col, 1)))
          new_v2 = g(grid_idx(col, (2 * row) - iand(col + 1, 1))) - (dt * (p_ijp1 - p_ij) / dy)
          v(grid_idx(col, (2 * row) - iand(col + 1, 1))) = new_v2

          new_v = max(abs(new_v), abs(new_v2))
          if (col == num) new_v = max(new_v, abs(v(grid_idx(num + 1, 2 * row))))
        else
          if (iand(col, 1) == 1) then
            p_ij = pres_red(pres_idx(col, row))
            p_ijp1 = pres_black(pres_idx(col, row + iand(col + 1, 1)))
            new_v = g(grid_idx(col, (2 * row) - iand(col, 1))) - (dt * (p_ijp1 - p_ij) / dy)
            v(grid_idx(col, (2 * row) - iand(col, 1))) = new_v
          else
            p_ij = pres_black(pres_idx(col, row))
            p_ijp1 = pres_red(pres_idx(col, row + iand(col, 1)))
            new_v = g(grid_idx(col, (2 * row) - iand(col + 1, 1))) - (dt * (p_ijp1 - p_ij) / dy)
            v(grid_idx(col, (2 * row) - iand(col + 1, 1))) = new_v
          end if

          new_v = abs(new_v)
          new_v = max(abs(v(grid_idx(col, num))), new_v)
          new_v = max(abs(v(grid_idx(col, 0))), new_v)
          new_v = max(abs(v(grid_idx(col, num + 1))), new_v)
        end if

        max_cache(lid) = new_v
        !$omp barrier
        i = block_size / 2
        do while (i /= 0)
          if (lid < i) max_cache(lid) = max(max_cache(lid), max_cache(lid + i))
          !$omp barrier
          i = i / 2
        end do

        if (lid == 0) max_v_arr(tid) = max_cache(0)
      !$omp end parallel
    !$omp end target teams
  end subroutine calculate_v_kernel

  subroutine set_bcs_kernel(u, v)
    real(dp), intent(inout) :: u(0:), v(0:)
    integer :: ind

    !$omp target teams distribute parallel do thread_limit(block_size) private(ind)
    do ind = 1, num
      u(grid_idx(0, ind)) = zero
      v(grid_idx(0, ind)) = -v(grid_idx(1, ind))

      u(grid_idx(num, ind)) = zero
      v(grid_idx(num + 1, ind)) = -v(grid_idx(num, ind))

      u(grid_idx(ind, 0)) = -u(grid_idx(ind, 1))
      v(grid_idx(ind, 0)) = zero

      u(grid_idx(ind, num + 1)) = two - u(grid_idx(ind, num))
      v(grid_idx(ind, num)) = zero

      if (ind == num) then
        u(grid_idx(0, 0)) = zero
        v(grid_idx(0, 0)) = -v(grid_idx(1, 0))
        u(grid_idx(0, num + 1)) = zero
        v(grid_idx(0, num + 1)) = -v(grid_idx(1, num + 1))

        u(grid_idx(num, 0)) = zero
        v(grid_idx(num + 1, 0)) = -v(grid_idx(num, 0))
        u(grid_idx(num, num + 1)) = zero
        v(grid_idx(num + 1, num + 1)) = -v(grid_idx(num, num + 1))

        u(grid_idx(0, 0)) = -u(grid_idx(0, 1))
        v(grid_idx(0, 0)) = zero
        u(grid_idx(num + 1, 0)) = -u(grid_idx(num + 1, 1))
        v(grid_idx(num + 1, 0)) = zero

        u(grid_idx(0, num + 1)) = two - u(grid_idx(0, num))
        v(grid_idx(0, num)) = zero
        u(grid_idx(num + 1, num + 1)) = two - u(grid_idx(num + 1, num))
        v(grid_idx(ind, num + 1)) = zero
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine set_bcs_kernel

  function fmt_exp(value) result(text)
    real(dp), intent(in) :: value
    character(len=32) :: text
    integer :: pos

    write(text, '(es13.6e2)') value
    text = adjustl(text)
    do pos = 1, len(text)
      if (text(pos:pos) == 'E') text(pos:pos) = 'e'
    end do
  end function fmt_exp

  function fmt_fixed(value) result(text)
    real(dp), intent(in) :: value
    character(len=32) :: text

    write(text, '(f32.6)') value
    text = adjustl(text)
  end function fmt_fixed

  subroutine write_velocity_file(u, v)
    real(dp), intent(in) :: u(0:), v(0:)
    integer :: row, col, unit
    real(dp) :: u_ij, u_im1j, v_ij, v_ijm1
    character(len=1), parameter :: tab = char(9)

    open(newunit=unit, file='velocity_gpu.dat', status='replace', action='write')
    write(unit, '(a)') '#x' // tab // 'y' // tab // 'u' // tab // 'v'
    do row = 0, num - 1
      do col = 0, num - 1
        u_ij = u(col * num + row)
        if (col == 0) then
          u_im1j = zero
        else
          u_im1j = u((col - 1) * num + row)
        end if
        u_ij = (u_ij + u_im1j) / two

        v_ij = v(col * num + row)
        if (row == 0) then
          v_ijm1 = zero
        else
          v_ijm1 = v(col * num + row - 1)
        end if
        v_ij = (v_ij + v_ijm1) / two

        write(unit, '(a,a,a,a,a,a,a)') trim(fmt_fixed((real(col, dp) + 0.5_dp) * dx)), tab, &
            trim(fmt_fixed((real(row, dp) + 0.5_dp) * dy)), tab, trim(fmt_fixed(u_ij)), tab, trim(fmt_fixed(v_ij))
      end do
    end do
    close(unit)
  end subroutine write_velocity_file

end module lid_driven_cavity_mod

program main
  use lid_driven_cavity_mod
  implicit none

  integer, parameter :: it_max = 1000000
  real(dp), parameter :: tol = 0.001_dp
  real(dp), parameter :: time_start_value = 0.0_dp
  real(dp), parameter :: time_end_value = 0.001_dp
  real(dp) :: dt, time_value, dt_re, max_u, max_v, p0_norm, norm_l2
  real(dp) :: start_time, end_time
  integer :: iter, i, col, row
  real(dp), allocatable :: f(:), u(:), g(:), v(:)
  real(dp), allocatable :: pres_red(:), pres_black(:)
  real(dp), allocatable :: res_arr(:), max_u_arr(:), max_v_arr(:), pres_sum(:)

  allocate(f(0:size_grid - 1), u(0:size_grid - 1), g(0:size_grid - 1), v(0:size_grid - 1))
  allocate(pres_red(0:size_pres - 1), pres_black(0:size_pres - 1))
  allocate(res_arr(0:size_res - 1), max_u_arr(0:size_res - 1), max_v_arr(0:size_res - 1), pres_sum(0:size_res - 1))

  f = zero
  u = zero
  g = zero
  v = zero
  pres_red = zero
  pres_black = zero
  res_arr = zero
  max_u_arr = zero
  max_v_arr = zero
  pres_sum = zero

  write(*, '("Problem size: ", i0, " x ", i0, " ")') num, num

  call set_bcs_host(u, v)

  max_u = small
  max_v = small
  do col = 0, num + 1
    do row = 1, num + 1
      max_u = max(max_u, abs(u(grid_idx(col, row))))
    end do
  end do
  do col = 1, num + 1
    do row = 0, num + 1
      max_v = max(max_v, abs(v(grid_idx(col, row))))
    end do
  end do

  !$omp target data map(tofrom: u(0:size_grid - 1), v(0:size_grid - 1), &
  !$omp& pres_red(0:size_pres - 1), pres_black(0:size_pres - 1)) &
  !$omp& map(to: f(0:size_grid - 1), g(0:size_grid - 1)) &
  !$omp& map(alloc: pres_sum(0:size_res - 1), res_arr(0:size_res - 1), &
  !$omp& max_u_arr(0:size_res - 1), max_v_arr(0:size_res - 1))
    time_value = time_start_value
    dt = 0.02_dp
    dt_re = 0.5_dp * re_num / ((one / (dx * dx)) + (one / (dy * dy)))
    start_time = omp_get_wtime()

    do while (time_value < time_end_value)
      dt = min(dx / max_u, dy / max_v)
      dt = tau * min(dt_re, dt)

      if ((time_value + dt) >= time_end_value) dt = time_end_value - time_value

      call calculate_f_kernel(dt, u, v, f)
      call calculate_g_kernel(dt, u, v, g)
      call sum_pressure_kernel(pres_red, pres_black, pres_sum)

      !$omp target update from(pres_sum(0:size_res - 1))
      p0_norm = zero
      do i = 0, size_res - 1
        p0_norm = p0_norm + pres_sum(i)
      end do
      p0_norm = sqrt(p0_norm / real(num * num, dp))
      if (p0_norm < 0.0001_dp) p0_norm = one

      norm_l2 = zero
      do iter = 1, it_max
        call set_horz_pres_bcs_kernel(pres_red, pres_black)
        call set_vert_pres_bcs_kernel(pres_red, pres_black)
        call red_kernel(dt, f, g, pres_black, pres_red)
        call black_kernel(dt, f, g, pres_red, pres_black)
        call calc_residual_kernel(dt, f, g, pres_red, pres_black, res_arr)

        !$omp target update from(res_arr(0:size_res - 1))
        norm_l2 = zero
        do i = 0, size_res - 1
          norm_l2 = norm_l2 + res_arr(i)
        end do
        norm_l2 = sqrt(norm_l2 / real(num * num, dp)) / p0_norm
        if (norm_l2 < tol) exit
      end do

      write(*, '("Time = ", f8.6, ", delt = ", a, ", iter = ", i0, ", res = ", a)') &
          time_value + dt, trim(fmt_exp(dt)), iter, trim(fmt_exp(norm_l2))

      call calculate_u_kernel(dt, f, pres_red, pres_black, u, max_u_arr)
      !$omp target update from(max_u_arr(0:size_res - 1))
      call calculate_v_kernel(dt, g, pres_red, pres_black, v, max_v_arr)
      !$omp target update from(max_v_arr(0:size_res - 1))

      max_u = small
      max_v = small
      do i = 0, size_res - 1
        max_u = max(max_u, max_u_arr(i))
        max_v = max(max_v, max_v_arr(i))
      end do

      call set_bcs_kernel(u, v)
      time_value = time_value + dt
    end do

    end_time = omp_get_wtime()
    write(*, *)
    write(*, '("Total execution time of the iteration loop: ", f0.6, " (s)")') end_time - start_time
  !$omp end target data

  call write_velocity_file(u, v)

  deallocate(pres_red, pres_black, u, v, f, g, max_u_arr, max_v_arr, res_arr, pres_sum)
end program main
