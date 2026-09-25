! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use, intrinsic :: iso_fortran_env, only : real64
  use omp_lib
  implicit none

  interface
    subroutine phmm_rng_seed(seed) bind(C, name="phmm_rng_seed")
      import :: c_int
      integer(c_int), value :: seed
    end subroutine phmm_rng_seed

    function phmm_rng_next() bind(C, name="phmm_rng_next") result(value)
      import :: c_double
      real(c_double) :: value
    end function phmm_rng_next
  end interface

  integer, parameter :: x_dim = 11
  integer, parameter :: y_dim = 40
  integer, parameter :: batch = 4
  integer, parameter :: states = 3
  integer, parameter :: nstate = states - 1

  character(len=256) :: arg0, arg
  integer :: repeat, count, i, j
  real(real64), allocatable, target :: cur_storage(:,:,:,:), next_storage(:,:,:,:)
  real(real64), allocatable, target :: cpu_cur_storage(:,:,:,:), cpu_next_storage(:,:,:,:)
  real(real64), pointer :: d_cur_forward(:,:,:,:), d_next_forward(:,:,:,:), tmp_forward(:,:,:,:)
  real(real64), pointer :: h_cur_forward(:,:,:,:), h_next_forward(:,:,:,:), h_tmp_forward(:,:,:,:)
  real(real64), allocatable :: emis(:,:,:,:), trans(:,:,:,:), like(:,:,:,:), start(:,:)
  real(real64) :: start_time, elapsed_ms, checksum
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 1) then
    write(*,'(A,A,A)') 'Usage: ', trim(arg0), ' <repeat>'
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) repeat

  allocate(cur_storage(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(next_storage(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(cpu_cur_storage(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(cpu_next_storage(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(emis(0:x_dim,0:y_dim,0:batch-1,0:nstate-1))
  allocate(trans(0:x_dim,0:batch-1,0:nstate-1,0:states-1))
  allocate(like(0:1,0:1,0:batch-1,0:nstate-1))
  allocate(start(0:batch-1,0:nstate-1))

  call initialize(cur_storage, emis, trans, like, start)
  next_storage = 0.0_real64
  cpu_cur_storage = cur_storage
  cpu_next_storage = 0.0_real64

  h_cur_forward => cpu_cur_storage
  h_next_forward => cpu_next_storage

  do count = 1, repeat
    do i = 1, x_dim
      do j = 1, y_dim
        call pair_hmm_forward_host(i, j, h_cur_forward, trans, emis, like, start, h_next_forward)
        h_tmp_forward => h_cur_forward
        h_cur_forward => h_next_forward
        h_next_forward => h_tmp_forward
      end do
    end do
  end do

  d_cur_forward => cur_storage
  d_next_forward => next_storage

  !$omp target data map(to: emis, trans, like, start) map(tofrom: cur_storage) map(alloc: next_storage)
  start_time = omp_get_wtime()
  do count = 1, repeat
    do i = 1, x_dim
      do j = 1, y_dim
        call pair_hmm_forward_device(i, j, d_cur_forward, trans, emis, like, start, d_next_forward)
        tmp_forward => d_cur_forward
        d_cur_forward => d_next_forward
        d_next_forward => tmp_forward
      end do
    end do
  end do
  elapsed_ms = (omp_get_wtime() - start_time) * 1.0e3_real64
  !$omp end target data

  ok = all(abs(h_cur_forward - d_cur_forward) <= 1.0e-9_real64)
  if (.not. ok) stop 1

  checksum = sum(d_cur_forward)
  write(*,'(A,F0.6,A)') 'Total execution time ', elapsed_ms, ' milliseconds'
  write(*,'(A,G0.6)') 'Checksum ', checksum

  deallocate(cur_storage, next_storage, cpu_cur_storage, cpu_next_storage, emis, trans, like, start)

contains

  subroutine initialize(cur, emis, trans, like, start)
    real(real64), intent(out) :: cur(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: emis(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: trans(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(out) :: like(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(out) :: start(0:batch-1,0:nstate-1)
    integer :: i, j, b, s, t, a, c

    call phmm_rng_seed(123_c_int)
    do i = 0, x_dim
      do j = 0, y_dim
        do b = 0, batch - 1
          do s = 0, nstate - 1
            cur(i,j,b,s) = next_random()
            emis(i,j,b,s) = next_random()
          end do
        end do
      end do
    end do
    do i = 0, x_dim
      do b = 0, batch - 1
        do s = 0, nstate - 1
          do t = 0, states - 1
            trans(i,b,s,t) = next_random()
          end do
        end do
      end do
    end do
    do b = 0, batch - 1
      do s = 0, nstate - 1
        start(b,s) = next_random()
      end do
    end do
    do a = 0, 1
      do c = 0, 1
        do b = 0, batch - 1
          do s = 0, nstate - 1
            like(a,c,b,s) = next_random()
          end do
        end do
      end do
    end do
  end subroutine initialize

  real(real64) function next_random()
    next_random = real(phmm_rng_next(), real64)
  end function next_random

  subroutine pair_hmm_forward_host(cur_i, cur_j, forward_in, transitions, emissions, likelihood, start_transitions, forward_out)
    integer, intent(in) :: cur_i, cur_j
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    real(real64), intent(inout) :: forward_out(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    integer :: b, s

    do b = 0, batch - 1
      do s = 0, nstate - 1
        forward_out(cur_i,cur_j,b,s) = pair_hmm_value(cur_i, cur_j, b, s, forward_in, transitions, emissions, likelihood, start_transitions)
      end do
    end do
  end subroutine pair_hmm_forward_host

  subroutine pair_hmm_forward_device(cur_i, cur_j, forward_in, transitions, emissions, likelihood, start_transitions, forward_out)
    integer, intent(in) :: cur_i, cur_j
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    real(real64), intent(inout) :: forward_out(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    !$omp target teams num_teams(batch) thread_limit(nstate)
    block
      real(real64) :: e(0:batch-1,0:nstate-1)
      real(real64) :: f01(0:0,0:batch-1,0:nstate-1)
      real(real64) :: mul_3d(0:0,0:batch-1,0:nstate-1)
      real(real64) :: mul_4d(0:3,0:batch-1,0:0,0:nstate-1)
      real(real64) :: t(0:1,0:1,0:batch-1,0:nstate-1,0:nstate-1)

      !$omp parallel
      block
        integer :: batch_id, states_id, k, l, m
        real(real64) :: t01(0:batch-1,0:nstate-1,0:nstate-1)
        real(real64) :: f(0:1,0:1,0:batch-1,0:0,0:nstate-1)
        real(real64) :: s_value, s0, s1, s2, s3, summation

        batch_id = omp_get_team_num()
        states_id = omp_get_thread_num()

        e(batch_id,states_id) = emissions(cur_i,cur_j,batch_id,states_id)

        do k = 0, nstate - 1
          do l = 0, nstate - 1
            t(0,0,batch_id,k,l) = transitions(cur_i - 1,batch_id,k,l)
            t(0,1,batch_id,k,l) = transitions(cur_i - 1,batch_id,k,l)
            t(1,0,batch_id,k,l) = transitions(cur_i,batch_id,k,l)
            t(1,1,batch_id,k,l) = transitions(cur_i,batch_id,k,l)
          end do
        end do
        !$omp barrier

        if (cur_i > 0 .and. cur_j == 0) then
          if (cur_i == 1) then
            forward_out(1,0,batch_id,states_id) = start_transitions(batch_id,states_id) * e(0,states_id)
          else
            do k = 0, nstate - 1
              do l = 0, nstate - 1
                t01(batch_id,k,l) = t(0,1,batch_id,k,l)
              end do
            end do

            f01(0,batch_id,states_id) = forward_in(cur_i - 1,cur_j,batch_id,states_id)
            !$omp barrier

            s_value = 0.0_real64
            do k = 0, nstate - 1
              s_value = s_value + f01(0,batch_id,k) * t01(batch_id,k,states_id)
            end do
            s_value = s_value * (e(batch_id,states_id) * likelihood(0,1,batch_id,states_id))
            mul_3d(0,batch_id,states_id) = s_value
            !$omp barrier

            forward_out(cur_i,0,batch_id,states_id) = mul_3d(0,batch_id,states_id)
          end if
        else if (cur_i > 0 .and. cur_j > 0) then
          do m = 0, nstate - 1
            f(0,0,batch_id,0,m) = forward_in(cur_i - 1,cur_j - 1,batch_id,m)
            f(0,1,batch_id,0,m) = forward_in(cur_i - 1,cur_j,batch_id,m)
            f(1,0,batch_id,0,m) = forward_in(cur_i,cur_j - 1,batch_id,m)
            f(1,1,batch_id,0,m) = forward_in(cur_i,cur_j,batch_id,m)
          end do
          !$omp barrier

          s0 = 0.0_real64
          s1 = 0.0_real64
          s2 = 0.0_real64
          s3 = 0.0_real64
          do k = 0, nstate - 1
            s0 = s0 + f(0,0,batch_id,0,k) * t(0,0,batch_id,k,states_id)
            s1 = s1 + f(0,1,batch_id,0,k) * t(0,1,batch_id,k,states_id)
            s2 = s2 + f(1,0,batch_id,0,k) * t(1,0,batch_id,k,states_id)
            s3 = s3 + f(1,1,batch_id,0,k) * t(1,1,batch_id,k,states_id)
          end do
          s0 = s0 * likelihood(0,0,batch_id,states_id)
          s1 = s1 * likelihood(0,1,batch_id,states_id)
          s2 = s2 * likelihood(1,0,batch_id,states_id)
          s3 = s3 * likelihood(1,1,batch_id,states_id)
          mul_4d(0,batch_id,0,states_id) = s0
          mul_4d(1,batch_id,0,states_id) = s1
          mul_4d(2,batch_id,0,states_id) = s2
          mul_4d(3,batch_id,0,states_id) = s3
          !$omp barrier

          do m = 0, nstate - 1
            summation = mul_4d(0,batch_id,0,m) + mul_4d(1,batch_id,0,m) + &
                        mul_4d(2,batch_id,0,m) + mul_4d(3,batch_id,0,m)
            summation = summation * e(batch_id,m)
            forward_out(cur_i,cur_j,batch_id,m) = summation
          end do
        end if
      end block
      !$omp end parallel
    end block
    !$omp end target teams
  end subroutine pair_hmm_forward_device

  real(real64) function pair_hmm_value(cur_i, cur_j, b, s, forward_in, transitions, emissions, likelihood, start_transitions)
    integer, intent(in) :: cur_i, cur_j, b, s
    real(real64), intent(in) :: forward_in(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: transitions(0:x_dim,0:batch-1,0:nstate-1,0:states-1)
    real(real64), intent(in) :: emissions(0:x_dim,0:y_dim,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: likelihood(0:1,0:1,0:batch-1,0:nstate-1)
    real(real64), intent(in) :: start_transitions(0:batch-1,0:nstate-1)
    integer :: k
    real(real64) :: total

    if (cur_i == 1 .and. cur_j == 0) then
      pair_hmm_value = start_transitions(b,s) * emissions(cur_i,cur_j,b,s)
    else if (cur_i > 0 .and. cur_j == 0) then
      total = 0.0_real64
      do k = 0, nstate - 1
        total = total + forward_in(cur_i - 1,cur_j,b,k) * transitions(cur_i - 1,b,k,s)
      end do
      pair_hmm_value = total * emissions(cur_i,cur_j,b,s) * likelihood(0,1,b,s)
    else
      total = 0.0_real64
      do k = 0, nstate - 1
        total = total + forward_in(cur_i - 1,cur_j - 1,b,k) * transitions(cur_i - 1,b,k,s) * likelihood(0,0,b,s)
        total = total + forward_in(cur_i - 1,cur_j,b,k) * transitions(cur_i - 1,b,k,s) * likelihood(0,1,b,s)
        total = total + forward_in(cur_i,cur_j - 1,b,k) * transitions(cur_i,b,k,s) * likelihood(1,0,b,s)
        total = total + forward_in(cur_i,cur_j,b,k) * transitions(cur_i,b,k,s) * likelihood(1,1,b,s)
      end do
      pair_hmm_value = total * emissions(cur_i,cur_j,b,s)
    end if
  end function pair_hmm_value

end program main
