! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_float, c_int, c_long
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use omp_lib
  implicit none

  interface
    subroutine c_srand(seed) bind(C, name='srand')
      import :: c_int
      integer(c_int), value :: seed
    end subroutine c_srand

    integer(c_int) function c_rand() bind(C, name='rand')
      import :: c_int
    end function c_rand
  end interface

  integer, parameter :: number_par_per_box = 100
  integer, parameter :: number_threads = 128
  real(real32), parameter :: alpha = 0.5_real32
  real(real32), parameter :: tolerance = 2.0e-3_real32

  type, bind(C) :: three_vector
    real(c_float) :: x, y, z
  end type three_vector

  type, bind(C) :: four_vector
    real(c_float) :: v, x, y, z
  end type four_vector

  type, bind(C) :: nei_str
    integer(c_int) :: x, y, z
    integer(c_int) :: number
    integer(c_long) :: offset
  end type nei_str

  type, bind(C) :: box_str
    integer(c_int) :: x, y, z
    integer(c_int) :: number
    integer(c_long) :: offset
    integer(c_int) :: nn
    type(nei_str) :: nei(26)
  end type box_str

  integer :: boxes1d, number_boxes, space_elem, i
  type(box_str), allocatable :: box_cpu(:)
  type(four_vector), allocatable :: rv_cpu(:), fv_cpu(:), ref_cpu(:)
  real(real32), allocatable :: qv_cpu(:)
  real(real64) :: start_total, end_total, start_kernel, end_kernel
  logical :: ok
  character(len=256) :: flag

  write(*,'(A,I0,A)') 'WG size of kernel = ', number_threads, ' '

  if (command_argument_count() /= 2) then
    write(*,'(A)', advance='no') 'Provide boxes1d argument, example: -boxes1d 16'
    stop 0
  end if
  call get_command_argument(1, flag)
  if (trim(flag) /= '-boxes1d') then
    write(*,'(A)') 'ERROR: Unknown argument'
    stop 0
  end if
  boxes1d = read_int_arg(2)
  if (boxes1d <= 0) then
    write(*,'(A)') 'ERROR: Wrong value to -boxes1d argument, cannot be <=0'
    stop 0
  end if

  write(*,'(A,I0,A,I0,A,I0)') 'Configuration used: arch = ', 0, ', cores = ', 1, ', boxes1d = ', boxes1d

  number_boxes = boxes1d * boxes1d * boxes1d
  space_elem = number_boxes * number_par_per_box
  allocate(box_cpu(number_boxes))
  allocate(rv_cpu(space_elem), fv_cpu(space_elem), ref_cpu(space_elem))
  allocate(qv_cpu(space_elem))

  call initialize_boxes(boxes1d, number_boxes, box_cpu)
  call c_srand(2_c_int)
  do i = 1, space_elem
    rv_cpu(i)%v = lava_random()
    rv_cpu(i)%x = lava_random()
    rv_cpu(i)%y = lava_random()
    rv_cpu(i)%z = lava_random()
  end do
  do i = 1, space_elem
    qv_cpu(i) = lava_random()
  end do
  fv_cpu%v = 0.0_real32
  fv_cpu%x = 0.0_real32
  fv_cpu%y = 0.0_real32
  fv_cpu%z = 0.0_real32
  ref_cpu%v = 0.0_real32
  ref_cpu%x = 0.0_real32
  ref_cpu%y = 0.0_real32
  ref_cpu%z = 0.0_real32

  start_total = omp_get_wtime()
  !$omp target data map(to: box_cpu(1:number_boxes), rv_cpu(1:space_elem), qv_cpu(1:space_elem)) &
  !$omp& map(tofrom: fv_cpu(1:space_elem))
  start_kernel = omp_get_wtime()
  call compute_forces_device(number_boxes, box_cpu, rv_cpu, qv_cpu, fv_cpu)
  end_kernel = omp_get_wtime()
  !$omp end target data
  end_total = omp_get_wtime()

  call compute_forces_host(number_boxes, box_cpu, rv_cpu, qv_cpu, ref_cpu)
  ok = compare_forces(space_elem, fv_cpu, ref_cpu)

  write(*,'(A)') 'Device offloading time:'
  call print_seconds(real(end_total - start_total, real32))
  write(*,'(A)') 'Kernel execution time:'
  call print_seconds(real(end_kernel - start_kernel, real32))

  deallocate(box_cpu, rv_cpu, qv_cpu, fv_cpu, ref_cpu)

contains

  integer function read_int_arg(pos) result(value)
    integer, intent(in) :: pos
    character(len=256) :: buffer
    call get_command_argument(pos, buffer)
    read(buffer, *) value
  end function read_int_arg

  real(real32) function lava_random() result(value)
    integer(c_int) :: sample
    sample = c_rand()
    value = real(mod(sample, 10_c_int) + 1_c_int, real32) / 10.0_real32
  end function lava_random

  subroutine print_seconds(value)
    real(real32), intent(in) :: value
    character(len=32) :: buffer
    write(buffer, '(F20.12)') value
    buffer = adjustl(buffer)
    if (buffer(1:1) == '.') then
      write(*,'(A,A,A)') '0', trim(buffer), ' s'
    else
      write(*,'(A,A)') trim(buffer), ' s'
    end if
  end subroutine print_seconds

  subroutine initialize_boxes(boxes1d, number_boxes, box_cpu)
    integer, intent(in) :: boxes1d, number_boxes
    type(box_str), intent(out) :: box_cpu(:)
    integer :: i, j, k, l, m, n, nh, neighbor, nn
    do nh = 1, number_boxes
      do nn = 1, 26
        box_cpu(nh)%nei(nn)%x = 0_c_int
        box_cpu(nh)%nei(nn)%y = 0_c_int
        box_cpu(nh)%nei(nn)%z = 0_c_int
        box_cpu(nh)%nei(nn)%number = 0_c_int
        box_cpu(nh)%nei(nn)%offset = 0_c_long
      end do
    end do
    nh = 0
    do i = 0, boxes1d - 1
      do j = 0, boxes1d - 1
        do k = 0, boxes1d - 1
          nh = nh + 1
          box_cpu(nh)%x = k
          box_cpu(nh)%y = j
          box_cpu(nh)%z = i
          box_cpu(nh)%number = nh - 1
          box_cpu(nh)%offset = int(nh - 1, c_long) * int(number_par_per_box, c_long)
          box_cpu(nh)%nn = 0
          do l = -1, 1
            do m = -1, 1
              do n = -1, 1
                if ((i + l) >= 0 .and. (j + m) >= 0 .and. (k + n) >= 0 .and. &
                    (i + l) < boxes1d .and. (j + m) < boxes1d .and. (k + n) < boxes1d .and. &
                    .not. (l == 0 .and. m == 0 .and. n == 0)) then
                  nn = box_cpu(nh)%nn + 1
                  neighbor = ((i + l) * boxes1d * boxes1d) + ((j + m) * boxes1d) + (k + n)
                  box_cpu(nh)%nei(nn)%x = k + n
                  box_cpu(nh)%nei(nn)%y = j + m
                  box_cpu(nh)%nei(nn)%z = i + l
                  box_cpu(nh)%nei(nn)%number = neighbor
                  box_cpu(nh)%nei(nn)%offset = int(neighbor, c_long) * int(number_par_per_box, c_long)
                  box_cpu(nh)%nn = nn
                end if
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine initialize_boxes

  subroutine compute_forces_device(number_boxes, box_cpu, rv_cpu, qv_cpu, fv_cpu)
    integer, intent(in) :: number_boxes
    type(box_str), intent(in) :: box_cpu(:)
    type(four_vector), intent(in) :: rv_cpu(:)
    real(real32), intent(in) :: qv_cpu(:)
    type(four_vector), intent(inout) :: fv_cpu(:)
    integer :: bx, tx, wtx, k, j, first_i, first_j, pointer, home_idx
    real(real32) :: a2, r2, u2, vij, fs, dx, dy, dz
    real(real32) :: rA_shared_v(number_par_per_box), rA_shared_x(number_par_per_box)
    real(real32) :: rA_shared_y(number_par_per_box), rA_shared_z(number_par_per_box)
    real(real32) :: rB_shared_v(number_par_per_box), rB_shared_x(number_par_per_box)
    real(real32) :: rB_shared_y(number_par_per_box), rB_shared_z(number_par_per_box)
    real(real32) :: qB_shared(number_par_per_box)

    !$omp target teams num_teams(number_boxes) thread_limit(number_threads) &
    !$omp& private(rA_shared_v, rA_shared_x, rA_shared_y, rA_shared_z, rB_shared_v, rB_shared_x, &
    !$omp& rB_shared_y, rB_shared_z, qB_shared, a2)
    a2 = 2.0_real32 * alpha * alpha
    !$omp parallel private(bx, tx, wtx, k, j, first_i, first_j, pointer, home_idx, r2, u2, vij, fs, dx, dy, dz)
    bx = omp_get_team_num() + 1
    tx = omp_get_thread_num()
    wtx = tx
    if (bx <= number_boxes) then
      first_i = int(box_cpu(bx)%offset)

      do while (wtx < number_par_per_box)
        home_idx = first_i + wtx + 1
        rA_shared_v(wtx + 1) = rv_cpu(home_idx)%v
        rA_shared_x(wtx + 1) = rv_cpu(home_idx)%x
        rA_shared_y(wtx + 1) = rv_cpu(home_idx)%y
        rA_shared_z(wtx + 1) = rv_cpu(home_idx)%z
        wtx = wtx + number_threads
      end do
      wtx = tx
      !$omp barrier

      do k = 0, box_cpu(bx)%nn
        if (k == 0) then
          pointer = bx
        else
          pointer = box_cpu(bx)%nei(k)%number + 1
        end if
        first_j = int(box_cpu(pointer)%offset)

        do while (wtx < number_par_per_box)
          j = first_j + wtx + 1
          rB_shared_v(wtx + 1) = rv_cpu(j)%v
          rB_shared_x(wtx + 1) = rv_cpu(j)%x
          rB_shared_y(wtx + 1) = rv_cpu(j)%y
          rB_shared_z(wtx + 1) = rv_cpu(j)%z
          qB_shared(wtx + 1) = qv_cpu(j)
          wtx = wtx + number_threads
        end do
        wtx = tx
        !$omp barrier

        do while (wtx < number_par_per_box)
          home_idx = first_i + wtx + 1
          do j = 1, number_par_per_box
            r2 = rA_shared_v(wtx + 1) + rB_shared_v(j) - &
                (rA_shared_x(wtx + 1) * rB_shared_x(j) + &
                 rA_shared_y(wtx + 1) * rB_shared_y(j) + &
                 rA_shared_z(wtx + 1) * rB_shared_z(j))
            u2 = a2 * r2
            vij = exp(-u2)
            fs = 2.0_real32 * vij
            dx = rA_shared_x(wtx + 1) - rB_shared_x(j)
            dy = rA_shared_y(wtx + 1) - rB_shared_y(j)
            dz = rA_shared_z(wtx + 1) - rB_shared_z(j)
            fv_cpu(home_idx)%v = fv_cpu(home_idx)%v + qB_shared(j) * vij
            fv_cpu(home_idx)%x = fv_cpu(home_idx)%x + qB_shared(j) * fs * dx
            fv_cpu(home_idx)%y = fv_cpu(home_idx)%y + qB_shared(j) * fs * dy
            fv_cpu(home_idx)%z = fv_cpu(home_idx)%z + qB_shared(j) * fs * dz
          end do
          wtx = wtx + number_threads
        end do
        wtx = tx
        !$omp barrier
      end do
    end if
    !$omp end parallel
    !$omp end target teams
  end subroutine compute_forces_device

  subroutine compute_forces_host(number_boxes, box_cpu, rv_cpu, qv_cpu, fv_cpu)
    integer, intent(in) :: number_boxes
    type(box_str), intent(in) :: box_cpu(:)
    type(four_vector), intent(in) :: rv_cpu(:)
    real(real32), intent(in) :: qv_cpu(:)
    type(four_vector), intent(inout) :: fv_cpu(:)
    integer :: bx, pi, k, j, first_i, first_j, pointer, home_idx, nei_idx
    real(real32) :: a2, r2, u2, vij, fs, dx, dy, dz

    a2 = 2.0_real32 * alpha * alpha
    do bx = 1, number_boxes
      first_i = int(box_cpu(bx)%offset)
      do pi = 1, number_par_per_box
        home_idx = first_i + pi
        do k = 0, box_cpu(bx)%nn
          if (k == 0) then
            pointer = bx
          else
            pointer = box_cpu(bx)%nei(k)%number + 1
          end if
          first_j = int(box_cpu(pointer)%offset)
          do j = 1, number_par_per_box
            nei_idx = first_j + j
            r2 = rv_cpu(home_idx)%v + rv_cpu(nei_idx)%v - &
                (rv_cpu(home_idx)%x * rv_cpu(nei_idx)%x + rv_cpu(home_idx)%y * rv_cpu(nei_idx)%y + &
                 rv_cpu(home_idx)%z * rv_cpu(nei_idx)%z)
            u2 = a2 * r2
            vij = exp(-u2)
            fs = 2.0_real32 * vij
            dx = rv_cpu(home_idx)%x - rv_cpu(nei_idx)%x
            dy = rv_cpu(home_idx)%y - rv_cpu(nei_idx)%y
            dz = rv_cpu(home_idx)%z - rv_cpu(nei_idx)%z
            fv_cpu(home_idx)%v = fv_cpu(home_idx)%v + qv_cpu(nei_idx) * vij
            fv_cpu(home_idx)%x = fv_cpu(home_idx)%x + qv_cpu(nei_idx) * fs * dx
            fv_cpu(home_idx)%y = fv_cpu(home_idx)%y + qv_cpu(nei_idx) * fs * dy
            fv_cpu(home_idx)%z = fv_cpu(home_idx)%z + qv_cpu(nei_idx) * fs * dz
          end do
        end do
      end do
    end do
  end subroutine compute_forces_host

  logical function compare_forces(n, fv_cpu, ref_cpu) result(ok)
    integer, intent(in) :: n
    type(four_vector), intent(in) :: fv_cpu(:), ref_cpu(:)
    integer :: i
    ok = .true.
    do i = 1, n
      if (abs(fv_cpu(i)%v - ref_cpu(i)%v) > tolerance .or. &
          abs(fv_cpu(i)%x - ref_cpu(i)%x) > tolerance .or. &
          abs(fv_cpu(i)%y - ref_cpu(i)%y) > tolerance .or. &
          abs(fv_cpu(i)%z - ref_cpu(i)%z) > tolerance) then
        ok = .false.
        return
      end if
    end do
  end function compare_forces

end program main
