! SPDX-License-Identifier: CC0-1.0
  subroutine pair32(pr, pi, a, b, idx, peer)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: idx, peer
    real(real32) :: tmp_real, tmp_imag
    tmp_real = pr(idx)
    tmp_imag = pi(idx)
    pr(idx) = a * tmp_real - b * pi(peer)
    pi(idx) = a * tmp_imag + b * pr(peer)
    pr(peer) = a * pr(peer) - b * tmp_imag
    pi(peer) = a * pi(peer) + b * tmp_real
  end subroutine pair32

  subroutine pair64(pr, pi, a, b, idx, peer)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: idx, peer
    real(real64) :: tmp_real, tmp_imag
    tmp_real = pr(idx)
    tmp_imag = pi(idx)
    pr(idx) = a * tmp_real - b * pi(peer)
    pi(idx) = a * tmp_imag + b * pr(peer)
    pr(peer) = a * pr(peer) - b * tmp_imag
    pi(peer) = a * pi(peer) + b * tmp_real
  end subroutine pair64

  subroutine kernel1_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 0, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel1_32

  subroutine kernel2_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel2_32

  subroutine kernel3_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 1, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel3_32

  subroutine kernel4_32(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel4_32

  subroutine kernel1_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 0, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel1_64

  subroutine kernel2_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel2_64

  subroutine kernel3_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    do x = 1, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel3_64

  subroutine kernel4_64(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
  end subroutine kernel4_64

  subroutine trotter_kernel32(p_real, p_imag, p2_real, p2_imag, a, b, width, height)
    integer, parameter :: STEPS = 1, BLOCK_X = 16, BLOCK_Y = 96, MARGIN_X = 3, MARGIN_Y = 4
    integer, parameter :: STRIDE_Y = 16, PARTS = BLOCK_Y / (STRIDE_Y * 2)
    real(real32), intent(in) :: p_real(0:), p_imag(0:)
    real(real32), intent(out) :: p2_real(0:), p2_imag(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: teamX, teamY, teams, threads

    teamX = (width + (BLOCK_X - 2 * STEPS * MARGIN_X) - 1) / (BLOCK_X - 2 * STEPS * MARGIN_X)
    teamY = (height + (BLOCK_Y - 2 * STEPS * MARGIN_Y) - 1) / (BLOCK_Y - 2 * STEPS * MARGIN_Y)
    teams = teamX * teamY
    threads = BLOCK_X * STRIDE_Y

    !$omp target teams num_teams(teams) thread_limit(threads)
    block
      real(real32) :: rl(0:BLOCK_Y-1, 0:BLOCK_X-1), im(0:BLOCK_Y-1, 0:BLOCK_X-1)
      !$omp parallel
      block
        integer :: blockIdx_x, blockIdx_y, threadIdx_x, threadIdx_y
        integer :: px, py, sx, sy, checkerboard_py, i, pidx, part
        integer :: ky, ky_peer, kx_peer, backwards, phase, write_sx, write_sy, write_px, write_py
        real(real32) :: cell_r(0:PARTS-1), cell_i(0:PARTS-1), peer_r, peer_i, old_cell_r, old_cell_i

        blockIdx_x = mod(omp_get_team_num(), teamX)
        blockIdx_y = omp_get_team_num() / teamX
        threadIdx_x = mod(omp_get_thread_num(), BLOCK_X)
        threadIdx_y = omp_get_thread_num() / BLOCK_X
        px = blockIdx_x * (BLOCK_X - 2 * STEPS * MARGIN_X) + threadIdx_x - STEPS * MARGIN_X
        py = blockIdx_y * (BLOCK_Y - 2 * STEPS * MARGIN_Y) + threadIdx_y - STEPS * MARGIN_Y

        if (px >= 0 .and. px < width) then
          do i = 0, BLOCK_Y / STRIDE_Y - 1
            if (py + i * STRIDE_Y >= 0 .and. py + i * STRIDE_Y < height) then
              pidx = (py + i * STRIDE_Y) * width + px
              rl(threadIdx_y + i * STRIDE_Y, threadIdx_x) = p_real(pidx)
              im(threadIdx_y + i * STRIDE_Y, threadIdx_x) = p_imag(pidx)
            end if
          end do
        end if

        !$omp barrier

        sx = threadIdx_x
        if (mod(STEPS * MARGIN_X, 2) == mod(STEPS * MARGIN_Y, 2)) then
          sy = 2 * threadIdx_y + mod(threadIdx_x, 2)
        else
          sy = 2 * threadIdx_y + 1 - mod(threadIdx_x, 2)
        end if
        checkerboard_py = blockIdx_y * (BLOCK_Y - 2 * STEPS * MARGIN_Y) + sy - STEPS * MARGIN_Y

        do part = 0, PARTS - 1
          cell_r(part) = rl(sy + part * 2 * STRIDE_Y, sx)
          cell_i(part) = im(sy + part * 2 * STRIDE_Y, sx)
        end do

        do phase = 1, 8
          select case (phase)
          case (1, 2, 7, 8)
            backwards = 0
          case default
            backwards = 1
          end select
          do part = 0, PARTS - 1
            ky = sy + part * 2 * STRIDE_Y
            if (phase == 1 .or. phase == 3 .or. phase == 6 .or. phase == 8) then
              ky_peer = ky + 1 - 2 * backwards
              if (checkerboard_py + part * 2 * STRIDE_Y >= backwards .and. &
                  checkerboard_py + part * 2 * STRIDE_Y < height - 1 + backwards .and. &
                  ky >= backwards .and. ky < BLOCK_Y - 1 + backwards) then
                peer_r = rl(ky_peer, sx)
                peer_i = im(ky_peer, sx)
                old_cell_r = cell_r(part)
                old_cell_i = cell_i(part)
                rl(ky_peer, sx) = a * peer_r - b * old_cell_i
                im(ky_peer, sx) = a * peer_i + b * old_cell_r
                cell_r(part) = a * old_cell_r - b * peer_i
                cell_i(part) = a * old_cell_i + b * peer_r
              end if
            else
              kx_peer = sx + 1 - 2 * backwards
              if (px >= backwards .and. px < width - 1 + backwards .and. &
                  sx >= backwards .and. sx < BLOCK_X - 1 + backwards) then
                peer_r = rl(ky, kx_peer)
                peer_i = im(ky, kx_peer)
                old_cell_r = cell_r(part)
                old_cell_i = cell_i(part)
                rl(ky, kx_peer) = a * peer_r - b * old_cell_i
                im(ky, kx_peer) = a * peer_i + b * old_cell_r
                cell_r(part) = a * old_cell_r - b * peer_i
                cell_i(part) = a * old_cell_i + b * peer_r
              end if
            end if
          end do
          !$omp barrier
        end do

        do part = 0, PARTS - 1
          rl(sy + part * 2 * STRIDE_Y, sx) = cell_r(part)
          im(sy + part * 2 * STRIDE_Y, sx) = cell_i(part)
        end do
        !$omp barrier

        write_sx = threadIdx_x + STEPS * MARGIN_X
        write_sy = threadIdx_y + STEPS * MARGIN_Y
        write_px = px + STEPS * MARGIN_X
        write_py = py + STEPS * MARGIN_Y
        if (write_sx < BLOCK_X - STEPS * MARGIN_X .and. write_px < width) then
          do i = 0, BLOCK_Y / STRIDE_Y - 1
            if (write_sy + i * STRIDE_Y < BLOCK_Y - STEPS * MARGIN_Y .and. &
                write_py + i * STRIDE_Y < height) then
              pidx = (write_py + i * STRIDE_Y) * width + write_px
              p2_real(pidx) = rl(write_sy + i * STRIDE_Y, write_sx)
              p2_imag(pidx) = im(write_sy + i * STRIDE_Y, write_sx)
            end if
          end do
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine trotter_kernel32

  subroutine trotter_kernel64(p_real, p_imag, p2_real, p2_imag, a, b, width, height)
    integer, parameter :: STEPS = 1, BLOCK_X = 16, BLOCK_Y = 32, MARGIN_X = 3, MARGIN_Y = 4
    integer, parameter :: STRIDE_Y = 16, PARTS = BLOCK_Y / (STRIDE_Y * 2)
    real(real64), intent(in) :: p_real(0:), p_imag(0:)
    real(real64), intent(out) :: p2_real(0:), p2_imag(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: teamX, teamY, teams, threads

    teamX = (width + (BLOCK_X - 2 * STEPS * MARGIN_X) - 1) / (BLOCK_X - 2 * STEPS * MARGIN_X)
    teamY = (height + (BLOCK_Y - 2 * STEPS * MARGIN_Y) - 1) / (BLOCK_Y - 2 * STEPS * MARGIN_Y)
    teams = teamX * teamY
    threads = BLOCK_X * STRIDE_Y

    !$omp target teams num_teams(teams) thread_limit(threads)
    block
      real(real64) :: rl(0:BLOCK_Y-1, 0:BLOCK_X-1), im(0:BLOCK_Y-1, 0:BLOCK_X-1)
      !$omp parallel
      block
        integer :: blockIdx_x, blockIdx_y, threadIdx_x, threadIdx_y
        integer :: px, py, sx, sy, checkerboard_py, i, pidx, part
        integer :: ky, ky_peer, kx_peer, backwards, phase, write_sx, write_sy, write_px, write_py
        real(real64) :: cell_r(0:PARTS-1), cell_i(0:PARTS-1), peer_r, peer_i, old_cell_r, old_cell_i

        blockIdx_x = mod(omp_get_team_num(), teamX)
        blockIdx_y = omp_get_team_num() / teamX
        threadIdx_x = mod(omp_get_thread_num(), BLOCK_X)
        threadIdx_y = omp_get_thread_num() / BLOCK_X
        px = blockIdx_x * (BLOCK_X - 2 * STEPS * MARGIN_X) + threadIdx_x - STEPS * MARGIN_X
        py = blockIdx_y * (BLOCK_Y - 2 * STEPS * MARGIN_Y) + threadIdx_y - STEPS * MARGIN_Y

        if (px >= 0 .and. px < width) then
          do i = 0, BLOCK_Y / STRIDE_Y - 1
            if (py + i * STRIDE_Y >= 0 .and. py + i * STRIDE_Y < height) then
              pidx = (py + i * STRIDE_Y) * width + px
              rl(threadIdx_y + i * STRIDE_Y, threadIdx_x) = p_real(pidx)
              im(threadIdx_y + i * STRIDE_Y, threadIdx_x) = p_imag(pidx)
            end if
          end do
        end if

        !$omp barrier

        sx = threadIdx_x
        if (mod(STEPS * MARGIN_X, 2) == mod(STEPS * MARGIN_Y, 2)) then
          sy = 2 * threadIdx_y + mod(threadIdx_x, 2)
        else
          sy = 2 * threadIdx_y + 1 - mod(threadIdx_x, 2)
        end if
        checkerboard_py = blockIdx_y * (BLOCK_Y - 2 * STEPS * MARGIN_Y) + sy - STEPS * MARGIN_Y

        do part = 0, PARTS - 1
          cell_r(part) = rl(sy + part * 2 * STRIDE_Y, sx)
          cell_i(part) = im(sy + part * 2 * STRIDE_Y, sx)
        end do

        do phase = 1, 8
          select case (phase)
          case (1, 2, 7, 8)
            backwards = 0
          case default
            backwards = 1
          end select
          do part = 0, PARTS - 1
            ky = sy + part * 2 * STRIDE_Y
            if (phase == 1 .or. phase == 3 .or. phase == 6 .or. phase == 8) then
              ky_peer = ky + 1 - 2 * backwards
              if (checkerboard_py + part * 2 * STRIDE_Y >= backwards .and. &
                  checkerboard_py + part * 2 * STRIDE_Y < height - 1 + backwards .and. &
                  ky >= backwards .and. ky < BLOCK_Y - 1 + backwards) then
                peer_r = rl(ky_peer, sx)
                peer_i = im(ky_peer, sx)
                old_cell_r = cell_r(part)
                old_cell_i = cell_i(part)
                rl(ky_peer, sx) = a * peer_r - b * old_cell_i
                im(ky_peer, sx) = a * peer_i + b * old_cell_r
                cell_r(part) = a * old_cell_r - b * peer_i
                cell_i(part) = a * old_cell_i + b * peer_r
              end if
            else
              kx_peer = sx + 1 - 2 * backwards
              if (px >= backwards .and. px < width - 1 + backwards .and. &
                  sx >= backwards .and. sx < BLOCK_X - 1 + backwards) then
                peer_r = rl(ky, kx_peer)
                peer_i = im(ky, kx_peer)
                old_cell_r = cell_r(part)
                old_cell_i = cell_i(part)
                rl(ky, kx_peer) = a * peer_r - b * old_cell_i
                im(ky, kx_peer) = a * peer_i + b * old_cell_r
                cell_r(part) = a * old_cell_r - b * peer_i
                cell_i(part) = a * old_cell_i + b * peer_r
              end if
            end if
          end do
          !$omp barrier
        end do

        do part = 0, PARTS - 1
          rl(sy + part * 2 * STRIDE_Y, sx) = cell_r(part)
          im(sy + part * 2 * STRIDE_Y, sx) = cell_i(part)
        end do
        !$omp barrier

        write_sx = threadIdx_x + STEPS * MARGIN_X
        write_sy = threadIdx_y + STEPS * MARGIN_Y
        write_px = px + STEPS * MARGIN_X
        write_py = py + STEPS * MARGIN_Y
        if (write_sx < BLOCK_X - STEPS * MARGIN_X .and. write_px < width) then
          do i = 0, BLOCK_Y / STRIDE_Y - 1
            if (write_sy + i * STRIDE_Y < BLOCK_Y - STEPS * MARGIN_Y .and. &
                write_py + i * STRIDE_Y < height) then
              pidx = (write_py + i * STRIDE_Y) * width + write_px
              p2_real(pidx) = rl(write_sy + i * STRIDE_Y, write_sx)
              p2_imag(pidx) = im(write_sy + i * STRIDE_Y, write_sx)
            end if
          end do
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine trotter_kernel64

  subroutine kernel1_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 0, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel1_32_device

  subroutine kernel2_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel2_32_device

  subroutine kernel3_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 1, width - 1, 2
      call pair32(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel3_32_device

  subroutine kernel4_32_device(pr, pi, a, b, width, height)
    real(real32), intent(inout) :: pr(0:), pi(0:)
    real(real32), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair32(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel4_32_device

  subroutine kernel1_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 0, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel1_64_device

  subroutine kernel2_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel2_64_device

  subroutine kernel3_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: x, y, idx, peer
    !$omp target teams distribute parallel do thread_limit(256)
    do x = 1, width - 1, 2
      call pair64(pr, pi, a, b, x, width + x)
    end do
    !$omp end target teams distribute parallel do
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 1, height - 2
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 1, 2
        peer = idx + width
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel3_64_device

  subroutine kernel4_64_device(pr, pi, a, b, width, height)
    real(real64), intent(inout) :: pr(0:), pi(0:)
    real(real64), intent(in) :: a, b
    integer, intent(in) :: width, height
    integer :: y, idx, peer
    !$omp target teams distribute parallel do private(idx, peer) thread_limit(256)
    do y = 0, height - 1
      do idx = y * width + 1 - mod(y, 2), (y + 1) * width - 2, 2
        peer = idx + 1
        call pair64(pr, pi, a, b, idx, peer)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine kernel4_64_device
