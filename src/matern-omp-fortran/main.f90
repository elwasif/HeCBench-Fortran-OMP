! SPDX-License-Identifier: CC0-1.0
program matern
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: real32, real64
  use omp_lib
  implicit none

  integer, parameter :: nsources = 50
  integer, parameter :: sx = 16
  integer, parameter :: sy = nsources
  integer(c_int), parameter :: c_rand_max = 2147483647_c_int
  real(real32), parameter :: sqrt5 = 2.2360679774997898_real32

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

  integer :: argc, npoints, repeat
  character(len=64) :: arg
  integer :: source_size, weight_size, ntargets, target_size, result_size
  real(real32), allocatable :: sources(:), targets(:), weights(:), result(:), result_ref(:)
  real(real32) :: l_values(7)
  character(len=7) :: l_text(7)
  integer :: i, li
  real(real64) :: start_time, end_time

  argc = command_argument_count()
  if (argc /= 2) then
    call get_command_argument(0, arg)
    print '(A,A,A)', 'Usage: ', trim(arg), ' <number of points> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) npoints
  call get_command_argument(2, arg)
  read(arg, *) repeat

  l_values = [0.1_real32, 1.0_real32, 10.0_real32, 100.0_real32, 1000.0_real32, &
              10000.0_real32, 100000.0_real32]
  l_text = ['1.0e-01', '1.0e+00', '1.0e+01', '1.0e+02', '1.0e+03', '1.0e+04', '1.0e+05']

  source_size = nsources * 3
  weight_size = nsources
  ntargets = npoints * npoints * npoints
  target_size = ntargets * 3
  result_size = ntargets

  allocate(sources(source_size), targets(target_size), weights(weight_size), result(result_size), result_ref(result_size))
  call initialize_data(sources, targets, weights)

  !$omp target data map(to: sources(1:source_size), weights(1:weight_size), targets(1:target_size)) &
  !$omp& map(alloc: result(1:result_size))
    print '(A)', '------------------------------------------------------------'
    print '(A)', 'Verifying the kernel results with the problem size (16 cube)'
    print '(A)', '------------------------------------------------------------'

    do li = 1, size(l_values)
      call matern_reference(16 * 16 * 16, l_values(li), sources, targets, weights, result_ref)
      call matern_kernel(16 * 16 * 16, l_values(li), sources, targets, weights, result)
      call check_result(16 * 16 * 16, result, result_ref, l_text(li))
    end do

    print '(A)', '--------------------------------------------------------------------'
    print '(A,I0,A)', 'Timing the kernel execution with the problem size (', npoints, ' cube)'
    print '(A)', '--------------------------------------------------------------------'

    do li = 1, size(l_values)
      print '(A)', 'Warmup..'
      do i = 1, repeat
        call matern_kernel(ntargets, l_values(li), sources, targets, weights, result)
      end do

      start_time = omp_get_wtime()
      do i = 1, repeat
        call matern_kernel(ntargets, l_values(li), sources, targets, weights, result)
      end do
      end_time = omp_get_wtime()
      print '(A,A,1X,A,F0.6,A)', 'Length scale = ', l_text(li), &
          'Average kernel execution time: ', (end_time - start_time) * 1.0e6_real64 / real(repeat, real64), ' (us)'
    end do
  !$omp end target data

  deallocate(sources, targets, weights, result, result_ref)

contains

  subroutine initialize_data(sources, targets, weights)
    real(real32), intent(out) :: sources(:), targets(:), weights(:)
    integer :: i

    call c_srand(123_c_int)
    do i = 1, size(sources)
      sources(i) = real(c_rand(), real32) / real(c_rand_max, real32)
    end do
    do i = 1, size(weights)
      weights(i) = real(c_rand(), real32) / real(c_rand_max, real32)
    end do
    do i = 1, size(targets)
      targets(i) = real(c_rand(), real32) / real(c_rand_max, real32)
    end do
  end subroutine initialize_data

  subroutine matern_kernel(num_targets, length_scale, sources, targets, weights, result)
    integer, intent(in) :: num_targets
    real(real32), intent(in) :: length_scale
    real(real32), intent(in) :: sources(:), targets(:), weights(:)
    real(real32), intent(inout) :: result(:)
    integer :: teams

    teams = (num_targets + sx - 1) / sx

    !$omp target teams num_teams(teams) thread_limit(sx * 64)
    block
      real(real32) :: local_result(sx * sy)
      real(real32) :: local_targets(sx * 3)
      real(real32) :: local_sources(sy * 3)
      real(real32) :: local_weights(sy)

      !$omp parallel
      block
        integer :: tx, ty, px, py, k
        real(real32) :: squared_diff, diff, res

        tx = mod(omp_get_thread_num(), sx)
        ty = omp_get_thread_num() / sx
        px = omp_get_team_num() * sx + tx
        py = ty

        if (px < num_targets .and. py < sy) then
          if (ty == 0) then
            do k = 1, 3
              local_targets(tx * 3 + k) = targets(px * 3 + k)
            end do
          end if

          if (tx == 0) then
            do k = 1, 3
              local_sources(ty * 3 + k) = sources(py * 3 + k)
            end do
            local_weights(ty + 1) = weights(ty + 1)
          end if
        end if
        !$omp barrier

        if (px < num_targets .and. py < sy) then
          squared_diff = 0.0_real32
          do k = 1, 3
            squared_diff = squared_diff + (local_targets(tx * 3 + k) - local_sources(ty * 3 + k)) * &
                (local_targets(tx * 3 + k) - local_sources(ty * 3 + k))
          end do
          diff = sqrt(squared_diff)
          local_result(tx * sy + ty + 1) = (1.0_real32 + sqrt5 * diff / length_scale + &
              5.0_real32 * squared_diff / (3.0_real32 * length_scale * length_scale)) * &
              exp(-sqrt5 * diff / length_scale) * local_weights(ty + 1)
        end if
        !$omp barrier

        if (px < num_targets .and. py < sy) then
          if (ty == 0) then
            res = 0.0_real32
            do k = 1, sy
              res = res + local_result(tx * sy + k)
            end do
            result(px + 1) = res
          end if
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
    !$omp target update from(result(1:num_targets))
  end subroutine matern_kernel

  subroutine matern_reference(num_targets, length_scale, sources, targets, weights, result)
    integer, intent(in) :: num_targets
    real(real32), intent(in) :: length_scale
    real(real32), intent(in) :: sources(:), targets(:), weights(:)
    real(real32), intent(out) :: result(:)
    integer :: t, s, k
    real(real32) :: squared_diff, diff, sum

    do t = 1, num_targets
      sum = 0.0_real32
      do s = 1, nsources
        squared_diff = 0.0_real32
        do k = 1, 3
          squared_diff = squared_diff + (sources((s - 1) * 3 + k) - targets((t - 1) * 3 + k)) * &
              (sources((s - 1) * 3 + k) - targets((t - 1) * 3 + k))
        end do
        diff = sqrt(squared_diff)
        sum = sum + (1.0_real32 + sqrt5 * diff / length_scale + &
            5.0_real32 * squared_diff / (3.0_real32 * length_scale * length_scale)) * &
            exp(-sqrt5 * diff / length_scale) * weights(s)
      end do
      result(t) = sum
    end do
  end subroutine matern_reference

  subroutine check_result(num_targets, result, result_ref, lstr)
    integer, intent(in) :: num_targets
    real(real32), intent(in) :: result(:), result_ref(:)
    character(len=*), intent(in) :: lstr
    integer :: i
    logical :: ok

    ok = .true.
    do i = 1, num_targets
      if (abs(result(i) - result_ref(i)) > 1.0e-3_real32) then
        print '(A,I0,A,F0.6,A,F0.6)', '@', i - 1, ' actual=', result(i), ' expected=', result_ref(i)
        ok = .false.
        exit
      end if
    end do
    if (ok) then
      print '(A,A,A)', 'Length scale = ', lstr, ' check = PASS'
    else
      print '(A,A,A)', 'Length scale = ', lstr, ' check = FAIL'
    end if
  end subroutine check_result

end program matern
