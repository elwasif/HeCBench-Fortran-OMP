! SPDX-License-Identifier: CC0-1.0
module asmooth_mod
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: real32
  use omp_lib
  implicit none

  integer, parameter :: block_size = 256

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

contains

  subroutine fill_image(img)
    real(real32), intent(out) :: img(0:)
    integer :: i

    call c_srand(123_c_int)
    do i = 0, size(img) - 1
      img(i) = real(mod(c_rand(), 256_c_int), real32)
    end do
  end subroutine fill_image

  subroutine reference(lx, ly, threshold, max_rad, img, box, norm, out)
    integer, intent(in) :: lx, ly, threshold, max_rad
    real(real32), intent(inout) :: img(0:)
    integer, intent(out) :: box(0:)
    real(real32), intent(inout) :: norm(0:), out(0:)
    integer :: x, y, i, j, s, q, ksum
    real(real32) :: sum_value

    do x = 0, lx - 1
      do y = 0, ly - 1
        sum_value = 0.0_real32
        s = 1
        q = 1
        ksum = 0

        do while (sum_value < real(threshold, real32) .and. q < max_rad)
          s = q
          sum_value = 0.0_real32
          ksum = 0

          do i = -s, s
            do j = -s, s
              if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
                sum_value = sum_value + img((x + i) * ly + y + j)
                ksum = ksum + 1
              end if
            end do
          end do
          q = q + 1
        end do

        box(x * ly + y) = s

        do i = -s, s
          do j = -s, s
            if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
              if (ksum /= 0) norm((x + i) * ly + y + j) = norm((x + i) * ly + y + j) + &
                1.0_real32 / real(ksum, real32)
            end if
          end do
        end do
      end do
    end do

    do x = 0, lx - 1
      do y = 0, ly - 1
        if (norm(x * ly + y) /= 0.0_real32) img(x * ly + y) = img(x * ly + y) / norm(x * ly + y)
      end do
    end do

    do x = 0, lx - 1
      do y = 0, ly - 1
        s = box(x * ly + y)
        sum_value = 0.0_real32
        ksum = 0

        do i = -s, s
          do j = -s, s
            if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
              sum_value = sum_value + img((x + i) * ly + y + j)
              ksum = ksum + 1
            end if
          end do
        end do
        if (ksum /= 0) out(x * ly + y) = sum_value / real(ksum, real32)
      end do
    end do
  end subroutine reference

  subroutine verify(size_total, max_rad, norm, h_norm, out, h_out, box, h_box)
    integer, intent(in) :: size_total, max_rad
    real(real32), intent(in) :: norm(0:), h_norm(0:), out(0:), h_out(0:)
    integer, intent(in) :: box(0:), h_box(0:)
    integer, allocatable :: cnt(:)
    integer :: i, j
    logical :: ok

    allocate(cnt(0:max_rad - 1))
    cnt = 0
    ok = .true.

    do i = 0, size_total - 1
      if (abs(norm(i) - h_norm(i)) > 1.0e-3_real32) then
        write(*,'(A,I0,1X,F0.6,1X,F0.6)') 'norm: ', i, norm(i), h_norm(i)
        ok = .false.
        exit
      end if
      if (abs(out(i) - h_out(i)) > 1.0e-3_real32) then
        write(*,'(A,I0,1X,F0.6,1X,F0.6)') 'out: ', i, out(i), h_out(i)
        ok = .false.
        exit
      end if
      if (box(i) /= h_box(i)) then
        write(*,'(A,I0,1X,I0,1X,I0)') 'box: ', i, box(i), h_box(i)
        ok = .false.
        exit
      else
        do j = 0, max_rad - 1
          if (box(i) == j) then
            cnt(j) = cnt(j) + 1
            exit
          end if
        end do
      end if
    end do

    if (ok) then
      write(*,'(A)') 'PASS'
      write(*,'(A)') 'Distribution of box sizes:'
      do j = 1, max_rad - 1
        write(*,'(A,I0,A,F0.6)') 'size=', j, ': ', real(cnt(j), real32) / real(size_total, real32)
      end do
    else
      write(*,'(A)') 'FAIL'
    end if

    deallocate(cnt)
  end subroutine verify

end module asmooth_mod

program main
  use iso_fortran_env, only: real32, real64
  use omp_lib
  use asmooth_mod
  implicit none

  integer :: argc, lx, ly, size_total, threshold, max_rad, repeat, status
  integer :: x, y, i, j, s, q, ksum, iter
  real(real32), allocatable :: img(:), norm(:), h_norm(:), out(:), h_out(:), original_img(:)
  integer, allocatable :: box(:), h_box(:)
  real(real32) :: sum_value
  real(real64) :: time_total, start_time
  character(len=256) :: arg, prog

  argc = command_argument_count()
  if (argc /= 4) then
    call get_command_argument(0, prog)
    write(*,'(A,A,A)') './', trim(prog), ' <image dimension> <threshold> <max box size> <iterations>'
    stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=status) lx
  if (status /= 0) stop 1
  ly = lx
  size_total = lx * ly

  call get_command_argument(2, arg)
  read(arg, *, iostat=status) threshold
  if (status /= 0) stop 1
  call get_command_argument(3, arg)
  read(arg, *, iostat=status) max_rad
  if (status /= 0) stop 1
  call get_command_argument(4, arg)
  read(arg, *, iostat=status) repeat
  if (status /= 0) stop 1

  allocate(img(0:size_total - 1), original_img(0:size_total - 1))
  allocate(norm(0:size_total - 1), h_norm(0:size_total - 1))
  allocate(box(0:size_total - 1), h_box(0:size_total - 1))
  allocate(out(0:size_total - 1), h_out(0:size_total - 1))

  call fill_image(original_img)
  img = original_img
  norm = 0.0_real32
  h_norm = 0.0_real32
  box = 0
  h_box = 0
  out = 0.0_real32
  h_out = 0.0_real32
  time_total = 0.0_real64

  !$omp target data map(alloc: img(0:size_total - 1), norm(0:size_total - 1), box(0:size_total - 1)) &
  !$omp& map(to: out(0:size_total - 1))
  do iter = 1, repeat
    img = original_img
    norm = 0.0_real32
    !$omp target update to(img(0:size_total - 1))
    !$omp target update to(norm(0:size_total - 1))

    start_time = omp_get_wtime()

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(sum_value, s, q, ksum, i, j)
    do x = 0, lx - 1
      do y = 0, ly - 1
        sum_value = 0.0_real32
        s = 1
        q = 1
        ksum = 0

        do while (sum_value < real(threshold, real32) .and. q < max_rad)
          s = q
          sum_value = 0.0_real32
          ksum = 0

          do i = -s, s
            do j = -s, s
              if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
                sum_value = sum_value + img((x + i) * ly + y + j)
                ksum = ksum + 1
              end if
            end do
          end do
          q = q + 1
        end do

        box(x * ly + y) = s

        do i = -s, s
          do j = -s, s
            if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
              if (ksum /= 0) then
                !$omp atomic update
                norm((x + i) * ly + y + j) = norm((x + i) * ly + y + j) + 1.0_real32 / real(ksum, real32)
              end if
            end if
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size)
    do x = 0, lx - 1
      do y = 0, ly - 1
        if (norm(x * ly + y) /= 0.0_real32) img(x * ly + y) = img(x * ly + y) / norm(x * ly + y)
      end do
    end do
    !$omp end target teams distribute parallel do

    !$omp target teams distribute parallel do collapse(2) thread_limit(block_size) &
    !$omp& private(s, sum_value, ksum, i, j)
    do x = 0, lx - 1
      do y = 0, ly - 1
        s = box(x * ly + y)
        sum_value = 0.0_real32
        ksum = 0

        do i = -s, s
          do j = -s, s
            if (x - s >= 0 .and. x + s < lx .and. y - s >= 0 .and. y + s < ly) then
              sum_value = sum_value + img((x + i) * ly + y + j)
              ksum = ksum + 1
            end if
          end do
        end do
        if (ksum /= 0) out(x * ly + y) = sum_value / real(ksum, real32)
      end do
    end do
    !$omp end target teams distribute parallel do

    time_total = time_total + (omp_get_wtime() - start_time)
  end do

  write(*,'(A,F0.6,A)') 'Average filtering time ', time_total / real(repeat, real64), ' (s)'
  !$omp target update from(out(0:size_total - 1))
  !$omp target update from(box(0:size_total - 1))
  !$omp target update from(norm(0:size_total - 1))
  !$omp end target data

  img = original_img
  call reference(lx, ly, threshold, max_rad, img, h_box, h_norm, h_out)
  call verify(size_total, max_rad, norm, h_norm, out, h_out, box, h_box)

  deallocate(img, original_img, norm, h_norm, box, h_box, out, h_out)
end program main
