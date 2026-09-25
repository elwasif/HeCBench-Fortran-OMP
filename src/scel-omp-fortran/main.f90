! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use, intrinsic :: iso_c_binding, only : c_int
  use omp_lib
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

  character(len=256) :: arg0, arg1, arg2, arg3
  integer :: outer_size, inner_size, repeat, input_size, output_size
  integer :: i, iter, log_d, log_d_trick, unjoined_lr_loss
  real(real32), allocatable :: logits(:), targets(:), output(:), ref_output(:)
  real(real32), parameter :: rand_max = 2147483647.0_real32
  real(real64) :: start_time, end_time, elapsed_us
  logical :: ok

  call get_command_argument(0, arg0)
  if (command_argument_count() /= 3) then
    print '(3A)', 'Usage: ', trim(arg0), ' <outer size> <inner_size> <repeat>'
    stop 1
  end if

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  read(arg1, *) outer_size
  read(arg2, *) inner_size
  read(arg3, *) repeat
  input_size = (outer_size + 1) * inner_size
  output_size = outer_size
  allocate(logits(input_size), targets(input_size), output(output_size), ref_output(output_size))

  call c_srand(123_c_int)
  do i = 1, input_size
    logits(i) = random_uniform_signed()
    targets(i) = random_uniform_signed() + 1.0_real32
  end do

  ok = .true.
  output = 0.0_real32
  ref_output = 0.0_real32

  !$omp target data map(to: logits(1:input_size), targets(1:input_size)) map(from: output(1:output_size))
  do unjoined_lr_loss = 0, 1
    if (unjoined_lr_loss == 0) then
      log_d = 1
    else
      log_d = 0
    end if

    do log_d_trick = 0, log_d
      start_time = omp_get_wtime()
      do iter = 1, repeat
        call sigmoid_cross_entropy_kernel(outer_size, inner_size, log_d_trick /= 0, &
            unjoined_lr_loss /= 0, logits, targets, output)
      end do
      end_time = omp_get_wtime()
      elapsed_us = (end_time - start_time) * 1.0e6_real64 / real(repeat, real64)
      print '(A,F0.6,A)', 'Average execution time of SigmoidCrossEntropyWithLogits kernel: ', elapsed_us, ' (us)'

      !$omp target update from(output(1:output_size))

      call reference(outer_size, inner_size, log_d_trick /= 0, unjoined_lr_loss /= 0, logits, targets, ref_output)
      do i = 1, output_size
        if (abs(ref_output(i) - output(i)) > 1.0e-3_real32) then
          ok = .false.
          exit
        end if
      end do
    end do
  end do
  !$omp end target data

  if (ok) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(logits, targets, output, ref_output)

contains

  function random_uniform_signed() result(value)
    real(real32) :: value

    value = real(c_rand(), real32) / rand_max
    value = value * 4.0_real32 - 2.0_real32
  end function random_uniform_signed

  real(real32) function sigmoid_xent_forward(lgt, tgt)
    real(real32), intent(in) :: lgt, tgt
    real(real32) :: gate

    gate = merge(1.0_real32, 0.0_real32, lgt >= 0.0_real32)
    sigmoid_xent_forward = lgt * (tgt - gate) - log(1.0_real32 + exp(lgt - 2.0_real32 * lgt * gate))
  end function sigmoid_xent_forward

  real(real32) function sigmoid_partition(lgt)
    real(real32), intent(in) :: lgt
    real(real32) :: gate

    gate = merge(1.0_real32, 0.0_real32, lgt >= 0.0_real32)
    sigmoid_partition = lgt * gate + log(1.0_real32 + exp(lgt - 2.0_real32 * lgt * gate))
  end function sigmoid_partition

  real(real32) function sigmoid_xent_forward_with_log_d_trick(lgt, tgt)
    real(real32), intent(in) :: lgt, tgt

    sigmoid_xent_forward_with_log_d_trick = (2.0_real32 * tgt - 1.0_real32) * (lgt - sigmoid_partition(lgt))
  end function sigmoid_xent_forward_with_log_d_trick

  real(real32) function unjoined_sigmoid_xent_forward(lgt, tgt)
    real(real32), intent(in) :: lgt, tgt
    real(real32) :: gate

    gate = merge(1.0_real32, 0.0_real32, lgt >= 0.0_real32)
    unjoined_sigmoid_xent_forward = lgt * tgt + (tgt - 1.0_real32) * lgt * gate - &
        (1.0_real32 - tgt) * log(1.0_real32 + exp(lgt - 2.0_real32 * lgt * gate))
  end function unjoined_sigmoid_xent_forward

  subroutine sigmoid_cross_entropy_kernel(outer_size, inner_size, log_d_trick, unjoined_lr_loss, logits, targets, output)
    integer, intent(in) :: outer_size, inner_size
    logical, intent(in) :: log_d_trick, unjoined_lr_loss
    real(real32), intent(in) :: logits(:), targets(:)
    real(real32), intent(out) :: output(:)
    integer :: outer, inner, idx
    real(real32) :: value, lgt, tgt

    !$omp target teams distribute num_teams(outer_size) private(inner, idx, value, lgt, tgt)
    do outer = 1, outer_size
      value = 0.0_real32
      !$omp parallel do reduction(+:value) num_threads(256) private(idx, lgt, tgt)
      do inner = 1, inner_size
        idx = (outer - 1) * inner_size + inner
        lgt = logits(idx)
        tgt = targets(idx)
        if (unjoined_lr_loss) then
          value = value + unjoined_sigmoid_xent_forward(lgt, tgt)
        else
          if (log_d_trick) then
            value = value + sigmoid_xent_forward_with_log_d_trick(lgt, tgt)
          else
            value = value + sigmoid_xent_forward(lgt, tgt)
          end if
        end if
      end do
      !$omp end parallel do
      output(outer) = -value / real(inner_size, real32)
    end do
    !$omp end target teams distribute
  end subroutine sigmoid_cross_entropy_kernel

  subroutine reference(outer_size, inner_size, log_d_trick, unjoined_lr_loss, logits, targets, output)
    integer, intent(in) :: outer_size, inner_size
    logical, intent(in) :: log_d_trick, unjoined_lr_loss
    real(real32), intent(in) :: logits(:), targets(:)
    real(real32), intent(out) :: output(:)
    integer :: outer, inner, idx
    real(real32) :: value, lgt, tgt

    do outer = 1, outer_size
      value = 0.0_real32
      do inner = 1, inner_size
        idx = (outer - 1) * inner_size + inner
        lgt = logits(idx)
        tgt = targets(idx)
        if (unjoined_lr_loss) then
          value = value + unjoined_sigmoid_xent_forward(lgt, tgt)
        else
          if (log_d_trick) then
            value = value + sigmoid_xent_forward_with_log_d_trick(lgt, tgt)
          else
            value = value + sigmoid_xent_forward(lgt, tgt)
          end if
        end if
      end do
      output(outer) = -value / real(inner_size, real32)
    end do
  end subroutine reference

end program main
