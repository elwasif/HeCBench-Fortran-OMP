! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  type, bind(C) :: float3
    real(c_float) :: x
    real(c_float) :: y
    real(c_float) :: z
    real(c_float) :: pad
  end type float3

  type, bind(C) :: float4
    real(c_float) :: x
    real(c_float) :: y
    real(c_float) :: z
    real(c_float) :: w
  end type float4

  type, bind(C) :: box
    integer(c_int) :: width
    integer(c_int) :: height
    integer(c_int) :: left
    integer(c_int) :: top
  end type box

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

  character(len=256) :: arg
  integer :: width, height, img_size, num_detections, i
  type(float3), allocatable :: input(:), output(:), ref_output(:)
  type(box), allocatable :: detections(:)
  type(float4) :: colors
  real(real64) :: start_time, end_time
  logical :: ok

  if (command_argument_count() /= 2) then
    call get_command_argument(0, arg)
    write(*,'(2A)', advance='no') 'Usage: ', trim(arg)
    write(*,'(A)') ' <width> <height>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *) width
  call get_command_argument(2, arg)
  read(arg, *) height
  if (width <= 64 .or. height <= 64) stop 1

  img_size = width * height
  num_detections = int(real(img_size, real32) * 0.8_real32)

  allocate(input(0:img_size - 1), output(0:img_size - 1), ref_output(0:img_size - 1), &
           detections(0:num_detections - 1))

  call initialize_inputs(input, output, ref_output, detections, img_size, width, height, num_detections)

  colors = float4(255.0_c_float, 204.0_c_float, 203.0_c_float, 1.0_c_float)

  !$omp target data map(to: input(0:img_size - 1)) map(tofrom: output(0:img_size - 1))
  start_time = omp_get_wtime()
  do i = 0, num_detections - 1
    call detection_overlay_box(input, output, width, height, detections(i)%left, detections(i)%top, &
                               detections(i)%width, detections(i)%height, colors)
  end do
  end_time = omp_get_wtime()
  !$omp end target data

  write(*,'(A,F0.6,A)') 'Total kernel execution time: ', end_time - start_time, ' (s)'

  call reference_overlay(input, ref_output, width, height, detections, num_detections, colors)

  ok = .true.
  do i = 0, img_size - 1
    if (abs(ref_output(i)%x - output(i)%x) > 1.0e-3_real32 .or. &
        abs(ref_output(i)%y - output(i)%y) > 1.0e-3_real32 .or. &
        abs(ref_output(i)%z - output(i)%z) > 1.0e-3_real32) then
      write(*,'(A,I0)') 'Error at index ', i
      ok = .false.
      exit
    end if
  end do

  if (ok) then
    write(*,'(A)') 'PASS'
  else
    write(*,'(A)') 'FAIL'
  end if

  deallocate(input, output, ref_output, detections)

contains

  subroutine initialize_inputs(input, output, ref_output, detections, img_size, width, height, num_detections)
    type(float3), intent(out) :: input(0:), output(0:), ref_output(0:)
    type(box), intent(out) :: detections(0:)
    integer, intent(in) :: img_size, width, height, num_detections
    integer :: i

    call c_srand(123_c_int)
    do i = 0, img_size - 1
      input(i)%x = real(rand_mod(256), c_float)
      input(i)%y = real(rand_mod(256), c_float)
      input(i)%z = real(rand_mod(256), c_float)
      input(i)%pad = 0.0_c_float
      output(i) = input(i)
      ref_output(i) = input(i)
    end do

    do i = 0, num_detections - 1
      detections(i)%width = 64 + rand_mod(128)
      detections(i)%height = 64 + rand_mod(128)
      detections(i)%left = rand_mod(width - 64)
      detections(i)%top = rand_mod(height - 64)
    end do
  end subroutine initialize_inputs

  integer function rand_mod(modulus)
    integer, intent(in) :: modulus

    rand_mod = int(modulo(c_rand(), int(modulus, c_int)))
  end function rand_mod

  subroutine detection_overlay_box(input, output, img_width, img_height, x0, y0, box_width, box_height, color)
    type(float3), intent(in) :: input(0:)
    type(float3), intent(inout) :: output(0:)
    integer, intent(in) :: img_width, img_height, x0, y0, box_width, box_height
    type(float4), intent(in) :: color
    integer :: box_x, box_y, x, y, idx
    type(float3) :: px
    real(real32) :: alpha, ialph

    !$omp target teams distribute parallel do collapse(2) thread_limit(64) &
    !$omp& private(box_x, box_y, x, y, idx, px, alpha, ialph)
    do box_y = 0, box_height - 1
      do box_x = 0, box_width - 1
        x = box_x + x0
        y = box_y + y0
        if (x < img_width .and. y < img_height) then
          idx = y * img_width + x
          px = input(idx)
          alpha = color%w / 255.0_real32
          ialph = 1.0_real32 - alpha
          px%x = alpha * color%x + ialph * px%x
          px%y = alpha * color%y + ialph * px%y
          px%z = alpha * color%z + ialph * px%z
          output(idx) = px
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine detection_overlay_box

  subroutine reference_overlay(input, output, img_width, img_height, detections, num_detections, color)
    type(float3), intent(in) :: input(0:)
    type(float3), intent(inout) :: output(0:)
    integer, intent(in) :: img_width, img_height, num_detections
    type(box), intent(in) :: detections(0:)
    type(float4), intent(in) :: color
    integer :: n, box_x, box_y, x, y, idx
    type(float3) :: px
    real(real32) :: alpha, ialph

    do n = 0, num_detections - 1
      do box_y = 0, detections(n)%height - 1
        do box_x = 0, detections(n)%width - 1
          x = box_x + detections(n)%left
          y = box_y + detections(n)%top
          if (x < img_width .and. y < img_height) then
            idx = y * img_width + x
            px = input(idx)
            alpha = color%w / 255.0_real32
            ialph = 1.0_real32 - alpha
            px%x = alpha * color%x + ialph * px%x
            px%y = alpha * color%y + ialph * px%y
            px%z = alpha * color%z + ialph * px%z
            output(idx) = px
          end if
        end do
      end do
    end do
  end subroutine reference_overlay

end program main
