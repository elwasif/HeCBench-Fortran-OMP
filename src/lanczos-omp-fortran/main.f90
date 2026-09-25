! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int32, real32, real64
  use omp_lib
  implicit none

  integer, parameter :: threads_per_block = 256
  character(len=512) :: graph_file
  integer :: node_count, eigen_count
  logical :: double_precision

  call parse_options(graph_file, node_count, eigen_count, double_precision)
  if (double_precision) then
    call run_lanczos_real64(graph_file, node_count, eigen_count)
  else
    call run_lanczos_real32(graph_file, node_count, eigen_count)
  end if

contains

  subroutine usage(program_name)
    character(len=*), intent(in) :: program_name
    write(*,'(A,A,A)') 'usage: ', trim(program_name), ' [options]'
    write(*,'(A)') 'options:'
    write(*,'(A)') '  -g --graph <file>'
    write(*,'(A)') '  -n --nodes <n>'
    write(*,'(A)') '  -k --eigens <k>'
    write(*,'(A)') '  -d --double'
  end subroutine usage

  subroutine parse_options(graph_file, node_count, eigen_count, double_precision)
    character(len=*), intent(out) :: graph_file
    integer, intent(out) :: node_count, eigen_count
    logical, intent(out) :: double_precision
    character(len=512) :: arg, program_name
    integer :: argc, i

    graph_file = ''
    node_count = 0
    eigen_count = 0
    double_precision = .false.
    argc = command_argument_count()
    call get_command_argument(0, program_name)
    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('-h', '--help', '-?')
        call usage(program_name)
        stop 0
      case ('-g', '--graph')
        if (i == argc) then
          call usage(program_name)
          stop 1
        end if
        i = i + 1
        call get_command_argument(i, graph_file)
      case ('-n', '--nodes')
        if (i == argc) then
          call usage(program_name)
          stop 1
        end if
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *) node_count
      case ('-k', '--eigens')
        if (i == argc) then
          call usage(program_name)
          stop 1
        end if
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *) eigen_count
      case ('-d', '--double')
        double_precision = .true.
      case default
        ! getopt_long leaves the Makefile's trailing "float" operand ignored.
      end select
      i = i + 1
    end do

    if (len_trim(graph_file) == 0) then
      write(0,'(A,A)') trim(program_name), ': missing graph file'
      stop 1
    end if
    if (node_count <= 0) then
      write(0,'(A,A)') trim(program_name), ': invalid node count'
      stop 1
    end if
    if (eigen_count <= 0 .or. eigen_count > node_count) then
      write(0,'(A,A)') trim(program_name), ': invalid eigenvalue count'
      stop 1
    end if
  end subroutine parse_options


  subroutine run_lanczos_real64(graph_file, node_count, eigen_count)
    character(len=*), intent(in) :: graph_file
    integer, intent(in) :: node_count, eigen_count
    integer, allocatable :: rows(:), cols(:), row_ptr(:), col_ind(:)
    real(real64), allocatable :: values(:), eigen(:)
    integer :: nonzeros, steps
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    call read_graph_real64(node_count, graph_file, rows, cols, values, nonzeros)
    call coo_to_csr_real64(node_count, rows, cols, values, nonzeros, row_ptr, col_ind)
    end_time = omp_get_wtime()
    write(*,'("graph load time: ",G0.15," sec")') end_time - start_time

    write(*,'(A)') '*** running GPU Lanczos ***'
    steps = 2 * eigen_count + 1
    do while (steps < 64 * eigen_count)
      call gpu_lanczos_eigen_real64(node_count, row_ptr, col_ind, values, nonzeros, eigen_count, steps, eigen)
      call print_vector_real64(eigen)
      write(*,*)
      if (allocated(eigen)) deallocate(eigen)
      steps = steps + 16
    end do

    write(*,'(A)') '*** running CPU Lanczos ***'
    steps = 2 * eigen_count + 1
    do while (steps < 8 * eigen_count)
      call cpu_lanczos_eigen_real64(node_count, row_ptr, col_ind, values, eigen_count, steps, eigen)
      call print_vector_real64(eigen)
      write(*,*)
      if (allocated(eigen)) deallocate(eigen)
      steps = steps + 16
    end do

    deallocate(col_ind, row_ptr, values, cols, rows)
  end subroutine run_lanczos_real64

  subroutine read_graph_real64(node_count, path, rows, cols, values, nonzeros)
    integer, intent(in) :: node_count
    character(len=*), intent(in) :: path
    integer, allocatable, intent(out) :: rows(:), cols(:)
    real(real64), allocatable, intent(out) :: values(:)
    integer, intent(out) :: nonzeros
    integer :: unit, ios, i, j, count

    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0,'(A,A)') 'Error: failed to open the file: ', trim(path)
      stop 1
    end if

    count = 0
    do
      read(unit, *, iostat=ios) i, j
      if (ios /= 0) exit
      if (i < node_count .and. j < node_count) count = count + 1
    end do
    rewind(unit)

    allocate(rows(0:count-1), cols(0:count-1), values(0:count-1))
    nonzeros = 0
    do
      read(unit, *, iostat=ios) i, j
      if (ios /= 0) exit
      if (i < node_count .and. j < node_count) then
        rows(nonzeros) = i
        cols(nonzeros) = j
        values(nonzeros) = 1.0_real64
        nonzeros = nonzeros + 1
      end if
    end do
    close(unit)
  end subroutine read_graph_real64

  subroutine coo_to_csr_real64(n, rows, cols, values, nonzeros, row_ptr, col_ind)
    integer, intent(in) :: n, nonzeros
    integer, intent(in) :: rows(0:), cols(0:)
    real(real64), intent(inout) :: values(0:)
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
    real(real64), allocatable :: csr_values(:)
    integer :: i, r, pos

    allocate(row_ptr(0:n), col_ind(0:nonzeros-1), csr_values(0:nonzeros-1))
    row_ptr = 0
    do i = 0, nonzeros - 1
      row_ptr(rows(i)) = row_ptr(rows(i)) + 1
    end do
    do i = 1, n
      row_ptr(i) = row_ptr(i) + row_ptr(i - 1)
    end do
    do i = 0, nonzeros - 1
      r = rows(i)
      pos = row_ptr(r) - 1
      row_ptr(r) = pos
      col_ind(pos) = cols(i)
      csr_values(pos) = values(i)
    end do
    values = csr_values
    deallocate(csr_values)
  end subroutine coo_to_csr_real64

  subroutine cpu_lanczos_eigen_real64(n, row_ptr, col_ind, values, k, steps, eigen)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:), k, steps
    real(real64), intent(in) :: values(0:)
    real(real64), allocatable, intent(out) :: eigen(:)
    real(real64), allocatable :: alpha(:), beta_arr(:), basis(:, :), r(:), product(:)
    real(real64) :: beta_value, alpha_value, start_time, end_time
    integer :: t

    allocate(alpha(0:steps-1), beta_arr(0:steps-2), basis(0:n-1,0:steps-1), r(0:n-1), product(0:n-1))
    r = 0.0_real64
    r(0) = 1.0_real64
    beta_value = l2_norm_real64(r)
    start_time = omp_get_wtime()
    do t = 0, steps - 1
      if (t > 0) beta_arr(t - 1) = beta_value
      call multiply_inplace_host_real64(r, 1.0_real64 / beta_value)
      basis(:, t) = r
      call csr_multiply_host_real64(n, row_ptr, col_ind, values, r, product)
      alpha_value = dot_product_host_real64(basis(:, t), product)
      call saxpy_inplace_host_real64(product, -alpha_value, basis(:, t))
      if (t > 0) call saxpy_inplace_host_real64(product, -beta_value, basis(:, t - 1))
      alpha(t) = alpha_value
      beta_value = l2_norm_real64(product)
      r = product
    end do
    end_time = omp_get_wtime()
    write(*,'("CPU Lanczos iterations: ",I0)') steps
    write(*,'("CPU Lanczos time: ",G0.15," sec")') end_time - start_time
    call lanczos_no_spurious_real64(alpha, beta_arr, steps, k, eigen)
    deallocate(product, r, basis, beta_arr, alpha)
  end subroutine cpu_lanczos_eigen_real64

  subroutine gpu_lanczos_eigen_real64(n, row_ptr, col_ind, values, nonzeros, k, steps, eigen)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:), nonzeros, k, steps
    real(real64), intent(in) :: values(0:)
    real(real64), allocatable, intent(out) :: eigen(:)
    real(real64), allocatable, target :: x_storage(:), x_prev_storage(:), y_storage(:)
    real(real64), pointer :: x(:), x_prev(:), y(:), tmp(:)
    real(real64), allocatable :: alpha(:), beta_arr(:)
    real(real64) :: product, start_time, end_time
    integer :: iter

    allocate(alpha(0:steps), beta_arr(0:steps-1))
    allocate(x_storage(0:n-1), x_prev_storage(0:n-1), y_storage(0:n-1))
    x_storage = 0.0_real64
    x_storage(0) = 1.0_real64
    x => x_storage
    x_prev => x_prev_storage
    y => y_storage

    !$omp target data map(to: row_ptr(0:n), col_ind(0:nonzeros-1), values(0:nonzeros-1)) &
    !$omp& map(to: x_storage(0:n-1)) map(alloc: x_prev_storage(0:n-1), y_storage(0:n-1))
      start_time = omp_get_wtime()
      do iter = 0, steps - 1
        call gpu_csr_multiply_real64(n, row_ptr, col_ind, values, x, y)
        call gpu_dot_product_real64(n, x, y, product)
        alpha(iter) = product
        call gpu_saxpy_inplace_real64(n, y, x, -product)
        if (iter > 0) call gpu_saxpy_inplace_real64(n, y, x_prev, -beta_arr(iter - 1))
        tmp => x
        x => x_prev
        x_prev => tmp
        call gpu_dot_product_real64(n, y, y, product)
        beta_arr(iter) = sqrt(product)
        call gpu_multiply_inplace_real64(n, y, 1.0_real64 / beta_arr(iter))
        tmp => x
        x => y
        y => tmp
      end do
      end_time = omp_get_wtime()
    !$omp end target data

    write(*,'("GPU Lanczos iterations: ",I0)') steps
    write(*,'("GPU Lanczos time: ",G0.15," sec")') end_time - start_time
    call lanczos_no_spurious_real64(alpha(0:steps-1), beta_arr(0:steps-2), steps, k, eigen)
    deallocate(y_storage, x_prev_storage, x_storage, beta_arr, alpha)
  end subroutine gpu_lanczos_eigen_real64

  subroutine gpu_csr_multiply_real64(n, row_ptr, col_ind, values, x, y)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:)
    real(real64), intent(in) :: values(0:), x(0:)
    real(real64), intent(inout) :: y(0:)
    integer :: row_nonzeros, group_size, groups_per_block, multiply_blocks
    integer :: lid, index, r, lane, start_idx, end_idx, jj, half
    real(real64) :: scratch(0:threads_per_block-1)

    row_nonzeros = max(1, (row_ptr(n) / max(1, n)))
    if (row_nonzeros > 16) then
      group_size = 32
    else if (row_nonzeros > 8) then
      group_size = 16
    else if (row_nonzeros > 4) then
      group_size = 8
    else if (row_nonzeros > 2) then
      group_size = 4
    else
      group_size = 2
    end if
    groups_per_block = threads_per_block / group_size
    multiply_blocks = (n + groups_per_block - 1) / groups_per_block

    !$omp target teams num_teams(multiply_blocks) thread_limit(threads_per_block) private(scratch)
      !$omp parallel private(lid, index, r, lane, start_idx, end_idx, jj, half)
        lid = omp_get_thread_num()
        index = omp_get_team_num() * omp_get_num_threads() + lid
        r = index / group_size
        lane = mod(index, group_size)
        scratch(lid) = 0.0_real64
        if (r < n) then
          start_idx = row_ptr(r)
          end_idx = row_ptr(r + 1)
          do jj = start_idx + lane, end_idx - 1, group_size
            scratch(lid) = scratch(lid) + values(jj) * x(col_ind(jj))
          end do
          half = group_size / 2
          do while (half > 0)
            !$omp barrier
            if (lane < half) scratch(lid) = scratch(lid) + scratch(lid + half)
            half = half / 2
          end do
          if (lane == 0) y(r) = scratch(lid)
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine gpu_csr_multiply_real64

  subroutine gpu_dot_product_real64(n, x, y, result)
    integer, intent(in) :: n
    real(real64), intent(in) :: x(0:), y(0:)
    real(real64), intent(out) :: result
    integer :: idx, blocks

    blocks = (n + threads_per_block - 1) / threads_per_block
    result = 0.0_real64
    !$omp target teams distribute parallel do reduction(+:result) num_teams(blocks) thread_limit(threads_per_block) map(tofrom: result)
    do idx = 0, n - 1
      result = result + x(idx) * y(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_dot_product_real64

  subroutine gpu_saxpy_inplace_real64(n, y, x, a)
    integer, intent(in) :: n
    real(real64), intent(inout) :: y(0:)
    real(real64), intent(in) :: x(0:), a
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do idx = 0, n - 1
      y(idx) = y(idx) + a * x(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_saxpy_inplace_real64

  subroutine gpu_multiply_inplace_real64(n, x, a)
    integer, intent(in) :: n
    real(real64), intent(inout) :: x(0:)
    real(real64), intent(in) :: a
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do idx = 0, n - 1
      x(idx) = x(idx) * a
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_multiply_inplace_real64

  subroutine csr_multiply_host_real64(n, row_ptr, col_ind, values, x, product)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:)
    real(real64), intent(in) :: values(0:), x(0:)
    real(real64), intent(out) :: product(0:)
    integer :: r, jj

    product = 0.0_real64
    do r = 0, n - 1
      do jj = row_ptr(r), row_ptr(r + 1) - 1
        product(r) = product(r) + values(jj) * x(col_ind(jj))
      end do
    end do
  end subroutine csr_multiply_host_real64

  real(real64) function dot_product_host_real64(x, y) result(s)
    real(real64), intent(in) :: x(:), y(:)
    integer :: i

    s = 0.0_real64
    do i = 1, size(x)
      s = s + x(i) * y(i)
    end do
  end function dot_product_host_real64

  real(real64) function l2_norm_real64(x) result(norm)
    real(real64), intent(in) :: x(:)
    norm = sqrt(dot_product_host_real64(x, x))
  end function l2_norm_real64

  subroutine multiply_inplace_host_real64(x, a)
    real(real64), intent(inout) :: x(:)
    real(real64), intent(in) :: a
    integer :: i
    do i = 1, size(x)
      x(i) = x(i) * a
    end do
  end subroutine multiply_inplace_host_real64

  subroutine saxpy_inplace_host_real64(y, a, x)
    real(real64), intent(inout) :: y(:)
    real(real64), intent(in) :: a, x(:)
    integer :: i
    do i = 1, size(y)
      y(i) = y(i) + a * x(i)
    end do
  end subroutine saxpy_inplace_host_real64

  subroutine lanczos_no_spurious_real64(alpha, beta_arr, n, k, result)
    real(real64), intent(in) :: alpha(0:), beta_arr(0:)
    integer, intent(in) :: n, k
    real(real64), allocatable, intent(out) :: result(:)
    real(real64), allocatable :: eigen(:), test_eigen(:), reduced_alpha(:), reduced_beta(:), tmp(:)
    real(real64), parameter :: epsilon = 1.0e-3_real64
    real(real64) :: start_time, end_time
    integer :: i, j, count

    start_time = omp_get_wtime()
    call tqlrat_eigen_real64(alpha, beta_arr, n, eigen)
    if (n > 1) then
      allocate(reduced_alpha(0:n-2), reduced_beta(0:max(0, n-3)))
      reduced_alpha = alpha(1:n-1)
      if (n > 2) reduced_beta = beta_arr(1:n-2)
      call tqlrat_eigen_real64(reduced_alpha, reduced_beta, n - 1, test_eigen)
      deallocate(reduced_alpha, reduced_beta)
    else
      allocate(test_eigen(0:0))
      test_eigen = huge(1.0_real64)
    end if

    allocate(tmp(0:n-1))
    count = 0
    i = 0
    j = 0
    do while (j <= size(eigen))
      if (j < size(eigen) .and. abs(eigen(j) - eigen(i)) < epsilon) then
        j = j + 1
      else
        if ((j - i > 1) .or. (.not. approximate_contains_real64(test_eigen, eigen(i), epsilon))) then
          tmp(count) = eigen(i)
          count = count + 1
        end if
        i = j
        j = j + 1
      end if
    end do

    call sort_descending_real64(tmp, count)
    count = min(count, k)
    allocate(result(0:count-1))
    if (count > 0) result = tmp(0:count-1)
    end_time = omp_get_wtime()
    write(*,'("spurious removal time: ",G0.15," sec")') end_time - start_time
    deallocate(tmp, test_eigen, eigen)
  end subroutine lanczos_no_spurious_real64

  subroutine tqlrat_eigen_real64(alpha, beta_arr, n, eigen)
    real(real64), intent(in) :: alpha(0:), beta_arr(0:)
    integer, intent(in) :: n
    real(real64), allocatable, intent(out) :: eigen(:)
    real(real64), allocatable :: d(:), e2(:)
    real(real64), parameter :: epsilon = 1.0e-8_real64
    real(real64) :: b, b2, f, h, g, p2, r2, s2
    real(real64) :: start_time, end_time
    integer :: kk, m, ii, jj

    start_time = omp_get_wtime()
    allocate(d(0:n-1), e2(0:n-1), eigen(0:n-1))
    d = alpha(0:n-1)
    e2 = 0.0_real64
    if (n > 1) then
      do ii = 0, n - 2
        e2(ii) = beta_arr(ii) * beta_arr(ii)
      end do
    end if
    b = 0.0_real64
    b2 = 0.0_real64
    f = 0.0_real64
    do kk = 0, n - 1
      h = epsilon * epsilon * (d(kk) * d(kk) + e2(kk))
      if (b2 < h) then
        b = sqrt(h)
        b2 = h
      end if
      m = kk
      do while (m < n .and. e2(m) > b2)
        m = m + 1
      end do
      if (m == n) m = m - 1
      if (m > kk) then
        do
          g = d(kk)
          p2 = sqrt(e2(kk))
          h = (d(kk + 1) - g) / (2.0_real64 * p2)
          r2 = sqrt(h * h + 1.0_real64)
          if (h < 0.0_real64) then
            h = p2 / (h - r2)
          else
            h = p2 / (h + r2)
          end if
          d(kk) = h
          h = g - h
          f = f + h
          do ii = kk + 1, n - 1
            d(ii) = d(ii) - h
          end do
          if (abs(d(m)) < epsilon) then
            h = b
          else
            h = d(m)
          end if
          g = h
          s2 = 0.0_real64
          do ii = m - 1, kk, -1
            p2 = g * h
            r2 = p2 + e2(ii)
            e2(ii + 1) = s2 * r2
            s2 = e2(ii) / r2
            d(ii + 1) = h + s2 * (h + d(ii))
            g = d(ii) - e2(ii) / g
            if (abs(g) < epsilon) g = b
            h = g * p2 / r2
          end do
          e2(kk) = s2 * g * h
          d(kk) = h
          if (e2(kk) <= b2) exit
        end do
      end if
      h = d(kk) + f
      jj = kk
      do while (jj > 0)
        if (h >= d(jj - 1)) exit
        d(jj) = d(jj - 1)
        jj = jj - 1
      end do
      d(jj) = h
    end do
    eigen = d
    end_time = omp_get_wtime()
    write(*,'("TQLRAT time: ",G0.15," sec")') end_time - start_time
    deallocate(e2, d)
  end subroutine tqlrat_eigen_real64

  logical function approximate_contains_real64(values, target, eps) result(found)
    real(real64), intent(in) :: values(:), target, eps
    integer :: i

    found = .false.
    do i = 1, size(values)
      if (abs(values(i) - target) < eps) then
        found = .true.
        return
      end if
    end do
  end function approximate_contains_real64

  subroutine sort_descending_real64(values, count)
    real(real64), intent(inout) :: values(0:)
    integer, intent(in) :: count
    integer :: i, j
    real(real64) :: tmp

    do i = 0, count - 2
      do j = i + 1, count - 1
        if (values(j) > values(i)) then
          tmp = values(i)
          values(i) = values(j)
          values(j) = tmp
        end if
      end do
    end do
  end subroutine sort_descending_real64

  subroutine print_vector_real64(values)
    real(real64), intent(in) :: values(0:)
    integer :: i

    do i = 0, size(values) - 1
      write(*,'(G0.15,1X)', advance='no') values(i)
    end do
    write(*,*)
  end subroutine print_vector_real64


  subroutine run_lanczos_real32(graph_file, node_count, eigen_count)
    character(len=*), intent(in) :: graph_file
    integer, intent(in) :: node_count, eigen_count
    integer, allocatable :: rows(:), cols(:), row_ptr(:), col_ind(:)
    real(real32), allocatable :: values(:), eigen(:)
    integer :: nonzeros, steps
    real(real64) :: start_time, end_time

    start_time = omp_get_wtime()
    call read_graph_real32(node_count, graph_file, rows, cols, values, nonzeros)
    call coo_to_csr_real32(node_count, rows, cols, values, nonzeros, row_ptr, col_ind)
    end_time = omp_get_wtime()
    write(*,'("graph load time: ",G0.15," sec")') end_time - start_time

    write(*,'(A)') '*** running GPU Lanczos ***'
    steps = 2 * eigen_count + 1
    do while (steps < 64 * eigen_count)
      call gpu_lanczos_eigen_real32(node_count, row_ptr, col_ind, values, nonzeros, eigen_count, steps, eigen)
      call print_vector_real32(eigen)
      write(*,*)
      if (allocated(eigen)) deallocate(eigen)
      steps = steps + 16
    end do

    write(*,'(A)') '*** running CPU Lanczos ***'
    steps = 2 * eigen_count + 1
    do while (steps < 8 * eigen_count)
      call cpu_lanczos_eigen_real32(node_count, row_ptr, col_ind, values, eigen_count, steps, eigen)
      call print_vector_real32(eigen)
      write(*,*)
      if (allocated(eigen)) deallocate(eigen)
      steps = steps + 16
    end do

    deallocate(col_ind, row_ptr, values, cols, rows)
  end subroutine run_lanczos_real32

  subroutine read_graph_real32(node_count, path, rows, cols, values, nonzeros)
    integer, intent(in) :: node_count
    character(len=*), intent(in) :: path
    integer, allocatable, intent(out) :: rows(:), cols(:)
    real(real32), allocatable, intent(out) :: values(:)
    integer, intent(out) :: nonzeros
    integer :: unit, ios, i, j, count

    open(newunit=unit, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0,'(A,A)') 'Error: failed to open the file: ', trim(path)
      stop 1
    end if

    count = 0
    do
      read(unit, *, iostat=ios) i, j
      if (ios /= 0) exit
      if (i < node_count .and. j < node_count) count = count + 1
    end do
    rewind(unit)

    allocate(rows(0:count-1), cols(0:count-1), values(0:count-1))
    nonzeros = 0
    do
      read(unit, *, iostat=ios) i, j
      if (ios /= 0) exit
      if (i < node_count .and. j < node_count) then
        rows(nonzeros) = i
        cols(nonzeros) = j
        values(nonzeros) = 1.0_real32
        nonzeros = nonzeros + 1
      end if
    end do
    close(unit)
  end subroutine read_graph_real32

  subroutine coo_to_csr_real32(n, rows, cols, values, nonzeros, row_ptr, col_ind)
    integer, intent(in) :: n, nonzeros
    integer, intent(in) :: rows(0:), cols(0:)
    real(real32), intent(inout) :: values(0:)
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
    real(real32), allocatable :: csr_values(:)
    integer :: i, r, pos

    allocate(row_ptr(0:n), col_ind(0:nonzeros-1), csr_values(0:nonzeros-1))
    row_ptr = 0
    do i = 0, nonzeros - 1
      row_ptr(rows(i)) = row_ptr(rows(i)) + 1
    end do
    do i = 1, n
      row_ptr(i) = row_ptr(i) + row_ptr(i - 1)
    end do
    do i = 0, nonzeros - 1
      r = rows(i)
      pos = row_ptr(r) - 1
      row_ptr(r) = pos
      col_ind(pos) = cols(i)
      csr_values(pos) = values(i)
    end do
    values = csr_values
    deallocate(csr_values)
  end subroutine coo_to_csr_real32

  subroutine cpu_lanczos_eigen_real32(n, row_ptr, col_ind, values, k, steps, eigen)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:), k, steps
    real(real32), intent(in) :: values(0:)
    real(real32), allocatable, intent(out) :: eigen(:)
    real(real32), allocatable :: alpha(:), beta_arr(:), basis(:, :), r(:), product(:)
    real(real32) :: beta_value, alpha_value, start_time, end_time
    integer :: t

    allocate(alpha(0:steps-1), beta_arr(0:steps-2), basis(0:n-1,0:steps-1), r(0:n-1), product(0:n-1))
    r = 0.0_real32
    r(0) = 1.0_real32
    beta_value = l2_norm_real32(r)
    start_time = omp_get_wtime()
    do t = 0, steps - 1
      if (t > 0) beta_arr(t - 1) = beta_value
      call multiply_inplace_host_real32(r, 1.0_real32 / beta_value)
      basis(:, t) = r
      call csr_multiply_host_real32(n, row_ptr, col_ind, values, r, product)
      alpha_value = dot_product_host_real32(basis(:, t), product)
      call saxpy_inplace_host_real32(product, -alpha_value, basis(:, t))
      if (t > 0) call saxpy_inplace_host_real32(product, -beta_value, basis(:, t - 1))
      alpha(t) = alpha_value
      beta_value = l2_norm_real32(product)
      r = product
    end do
    end_time = omp_get_wtime()
    write(*,'("CPU Lanczos iterations: ",I0)') steps
    write(*,'("CPU Lanczos time: ",G0.15," sec")') end_time - start_time
    call lanczos_no_spurious_real32(alpha, beta_arr, steps, k, eigen)
    deallocate(product, r, basis, beta_arr, alpha)
  end subroutine cpu_lanczos_eigen_real32

  subroutine gpu_lanczos_eigen_real32(n, row_ptr, col_ind, values, nonzeros, k, steps, eigen)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:), nonzeros, k, steps
    real(real32), intent(in) :: values(0:)
    real(real32), allocatable, intent(out) :: eigen(:)
    real(real32), allocatable, target :: x_storage(:), x_prev_storage(:), y_storage(:)
    real(real32), pointer :: x(:), x_prev(:), y(:), tmp(:)
    real(real32), allocatable :: alpha(:), beta_arr(:)
    real(real32) :: product, start_time, end_time
    integer :: iter

    allocate(alpha(0:steps), beta_arr(0:steps-1))
    allocate(x_storage(0:n-1), x_prev_storage(0:n-1), y_storage(0:n-1))
    x_storage = 0.0_real32
    x_storage(0) = 1.0_real32
    x => x_storage
    x_prev => x_prev_storage
    y => y_storage

    !$omp target data map(to: row_ptr(0:n), col_ind(0:nonzeros-1), values(0:nonzeros-1)) &
    !$omp& map(to: x_storage(0:n-1)) map(alloc: x_prev_storage(0:n-1), y_storage(0:n-1))
      start_time = omp_get_wtime()
      do iter = 0, steps - 1
        call gpu_csr_multiply_real32(n, row_ptr, col_ind, values, x, y)
        call gpu_dot_product_real32(n, x, y, product)
        alpha(iter) = product
        call gpu_saxpy_inplace_real32(n, y, x, -product)
        if (iter > 0) call gpu_saxpy_inplace_real32(n, y, x_prev, -beta_arr(iter - 1))
        tmp => x
        x => x_prev
        x_prev => tmp
        call gpu_dot_product_real32(n, y, y, product)
        beta_arr(iter) = sqrt(product)
        call gpu_multiply_inplace_real32(n, y, 1.0_real32 / beta_arr(iter))
        tmp => x
        x => y
        y => tmp
      end do
      end_time = omp_get_wtime()
    !$omp end target data

    write(*,'("GPU Lanczos iterations: ",I0)') steps
    write(*,'("GPU Lanczos time: ",G0.15," sec")') end_time - start_time
    call lanczos_no_spurious_real32(alpha(0:steps-1), beta_arr(0:steps-2), steps, k, eigen)
    deallocate(y_storage, x_prev_storage, x_storage, beta_arr, alpha)
  end subroutine gpu_lanczos_eigen_real32

  subroutine gpu_csr_multiply_real32(n, row_ptr, col_ind, values, x, y)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:)
    real(real32), intent(in) :: values(0:), x(0:)
    real(real32), intent(inout) :: y(0:)
    integer :: row_nonzeros, group_size, groups_per_block, multiply_blocks
    integer :: lid, index, r, lane, start_idx, end_idx, jj, half
    real(real32) :: scratch(0:threads_per_block-1)

    row_nonzeros = max(1, (row_ptr(n) / max(1, n)))
    if (row_nonzeros > 16) then
      group_size = 32
    else if (row_nonzeros > 8) then
      group_size = 16
    else if (row_nonzeros > 4) then
      group_size = 8
    else if (row_nonzeros > 2) then
      group_size = 4
    else
      group_size = 2
    end if
    groups_per_block = threads_per_block / group_size
    multiply_blocks = (n + groups_per_block - 1) / groups_per_block

    !$omp target teams num_teams(multiply_blocks) thread_limit(threads_per_block) private(scratch)
      !$omp parallel private(lid, index, r, lane, start_idx, end_idx, jj, half)
        lid = omp_get_thread_num()
        index = omp_get_team_num() * omp_get_num_threads() + lid
        r = index / group_size
        lane = mod(index, group_size)
        scratch(lid) = 0.0_real32
        if (r < n) then
          start_idx = row_ptr(r)
          end_idx = row_ptr(r + 1)
          do jj = start_idx + lane, end_idx - 1, group_size
            scratch(lid) = scratch(lid) + values(jj) * x(col_ind(jj))
          end do
          half = group_size / 2
          do while (half > 0)
            !$omp barrier
            if (lane < half) scratch(lid) = scratch(lid) + scratch(lid + half)
            half = half / 2
          end do
          if (lane == 0) y(r) = scratch(lid)
        end if
      !$omp end parallel
    !$omp end target teams
  end subroutine gpu_csr_multiply_real32

  subroutine gpu_dot_product_real32(n, x, y, result)
    integer, intent(in) :: n
    real(real32), intent(in) :: x(0:), y(0:)
    real(real32), intent(out) :: result
    integer :: idx, blocks

    blocks = (n + threads_per_block - 1) / threads_per_block
    result = 0.0_real32
    !$omp target teams distribute parallel do reduction(+:result) num_teams(blocks) thread_limit(threads_per_block) map(tofrom: result)
    do idx = 0, n - 1
      result = result + x(idx) * y(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_dot_product_real32

  subroutine gpu_saxpy_inplace_real32(n, y, x, a)
    integer, intent(in) :: n
    real(real32), intent(inout) :: y(0:)
    real(real32), intent(in) :: x(0:), a
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do idx = 0, n - 1
      y(idx) = y(idx) + a * x(idx)
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_saxpy_inplace_real32

  subroutine gpu_multiply_inplace_real32(n, x, a)
    integer, intent(in) :: n
    real(real32), intent(inout) :: x(0:)
    real(real32), intent(in) :: a
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do idx = 0, n - 1
      x(idx) = x(idx) * a
    end do
    !$omp end target teams distribute parallel do
  end subroutine gpu_multiply_inplace_real32

  subroutine csr_multiply_host_real32(n, row_ptr, col_ind, values, x, product)
    integer, intent(in) :: n, row_ptr(0:), col_ind(0:)
    real(real32), intent(in) :: values(0:), x(0:)
    real(real32), intent(out) :: product(0:)
    integer :: r, jj

    product = 0.0_real32
    do r = 0, n - 1
      do jj = row_ptr(r), row_ptr(r + 1) - 1
        product(r) = product(r) + values(jj) * x(col_ind(jj))
      end do
    end do
  end subroutine csr_multiply_host_real32

  real(real32) function dot_product_host_real32(x, y) result(s)
    real(real32), intent(in) :: x(:), y(:)
    integer :: i

    s = 0.0_real32
    do i = 1, size(x)
      s = s + x(i) * y(i)
    end do
  end function dot_product_host_real32

  real(real32) function l2_norm_real32(x) result(norm)
    real(real32), intent(in) :: x(:)
    norm = sqrt(dot_product_host_real32(x, x))
  end function l2_norm_real32

  subroutine multiply_inplace_host_real32(x, a)
    real(real32), intent(inout) :: x(:)
    real(real32), intent(in) :: a
    integer :: i
    do i = 1, size(x)
      x(i) = x(i) * a
    end do
  end subroutine multiply_inplace_host_real32

  subroutine saxpy_inplace_host_real32(y, a, x)
    real(real32), intent(inout) :: y(:)
    real(real32), intent(in) :: a, x(:)
    integer :: i
    do i = 1, size(y)
      y(i) = y(i) + a * x(i)
    end do
  end subroutine saxpy_inplace_host_real32

  subroutine lanczos_no_spurious_real32(alpha, beta_arr, n, k, result)
    real(real32), intent(in) :: alpha(0:), beta_arr(0:)
    integer, intent(in) :: n, k
    real(real32), allocatable, intent(out) :: result(:)
    real(real32), allocatable :: eigen(:), test_eigen(:), reduced_alpha(:), reduced_beta(:), tmp(:)
    real(real32), parameter :: epsilon = 1.0e-3_real32
    real(real64) :: start_time, end_time
    integer :: i, j, count

    start_time = omp_get_wtime()
    call tqlrat_eigen_real32(alpha, beta_arr, n, eigen)
    if (n > 1) then
      allocate(reduced_alpha(0:n-2), reduced_beta(0:max(0, n-3)))
      reduced_alpha = alpha(1:n-1)
      if (n > 2) reduced_beta = beta_arr(1:n-2)
      call tqlrat_eigen_real32(reduced_alpha, reduced_beta, n - 1, test_eigen)
      deallocate(reduced_alpha, reduced_beta)
    else
      allocate(test_eigen(0:0))
      test_eigen = huge(1.0_real32)
    end if

    allocate(tmp(0:n-1))
    count = 0
    i = 0
    j = 0
    do while (j <= size(eigen))
      if (j < size(eigen) .and. abs(eigen(j) - eigen(i)) < epsilon) then
        j = j + 1
      else
        if ((j - i > 1) .or. (.not. approximate_contains_real32(test_eigen, eigen(i), epsilon))) then
          tmp(count) = eigen(i)
          count = count + 1
        end if
        i = j
        j = j + 1
      end if
    end do

    call sort_descending_real32(tmp, count)
    count = min(count, k)
    allocate(result(0:count-1))
    if (count > 0) result = tmp(0:count-1)
    end_time = omp_get_wtime()
    write(*,'("spurious removal time: ",G0.15," sec")') end_time - start_time
    deallocate(tmp, test_eigen, eigen)
  end subroutine lanczos_no_spurious_real32

  subroutine tqlrat_eigen_real32(alpha, beta_arr, n, eigen)
    real(real32), intent(in) :: alpha(0:), beta_arr(0:)
    integer, intent(in) :: n
    real(real32), allocatable, intent(out) :: eigen(:)
    real(real32), allocatable :: d(:), e2(:)
    real(real32), parameter :: epsilon = 1.0e-8_real32
    real(real32) :: b, b2, f, h, g, p2, r2, s2
    real(real64) :: start_time, end_time
    integer :: kk, m, ii, jj

    start_time = omp_get_wtime()
    allocate(d(0:n-1), e2(0:n-1), eigen(0:n-1))
    d = alpha(0:n-1)
    e2 = 0.0_real32
    if (n > 1) then
      do ii = 0, n - 2
        e2(ii) = beta_arr(ii) * beta_arr(ii)
      end do
    end if
    b = 0.0_real32
    b2 = 0.0_real32
    f = 0.0_real32
    do kk = 0, n - 1
      h = epsilon * epsilon * (d(kk) * d(kk) + e2(kk))
      if (b2 < h) then
        b = sqrt(h)
        b2 = h
      end if
      m = kk
      do while (m < n .and. e2(m) > b2)
        m = m + 1
      end do
      if (m == n) m = m - 1
      if (m > kk) then
        do
          g = d(kk)
          p2 = sqrt(e2(kk))
          h = (d(kk + 1) - g) / (2.0_real32 * p2)
          r2 = sqrt(h * h + 1.0_real32)
          if (h < 0.0_real32) then
            h = p2 / (h - r2)
          else
            h = p2 / (h + r2)
          end if
          d(kk) = h
          h = g - h
          f = f + h
          do ii = kk + 1, n - 1
            d(ii) = d(ii) - h
          end do
          if (abs(d(m)) < epsilon) then
            h = b
          else
            h = d(m)
          end if
          g = h
          s2 = 0.0_real32
          do ii = m - 1, kk, -1
            p2 = g * h
            r2 = p2 + e2(ii)
            e2(ii + 1) = s2 * r2
            s2 = e2(ii) / r2
            d(ii + 1) = h + s2 * (h + d(ii))
            g = d(ii) - e2(ii) / g
            if (abs(g) < epsilon) g = b
            h = g * p2 / r2
          end do
          e2(kk) = s2 * g * h
          d(kk) = h
          if (e2(kk) <= b2) exit
        end do
      end if
      h = d(kk) + f
      jj = kk
      do while (jj > 0)
        if (h >= d(jj - 1)) exit
        d(jj) = d(jj - 1)
        jj = jj - 1
      end do
      d(jj) = h
    end do
    eigen = d
    end_time = omp_get_wtime()
    write(*,'("TQLRAT time: ",G0.15," sec")') end_time - start_time
    deallocate(e2, d)
  end subroutine tqlrat_eigen_real32

  logical function approximate_contains_real32(values, target, eps) result(found)
    real(real32), intent(in) :: values(:), target, eps
    integer :: i

    found = .false.
    do i = 1, size(values)
      if (abs(values(i) - target) < eps) then
        found = .true.
        return
      end if
    end do
  end function approximate_contains_real32

  subroutine sort_descending_real32(values, count)
    real(real32), intent(inout) :: values(0:)
    integer, intent(in) :: count
    integer :: i, j
    real(real32) :: tmp

    do i = 0, count - 2
      do j = i + 1, count - 1
        if (values(j) > values(i)) then
          tmp = values(i)
          values(i) = values(j)
          values(j) = tmp
        end if
      end do
    end do
  end subroutine sort_descending_real32

  subroutine print_vector_real32(values)
    real(real32), intent(in) :: values(0:)
    integer :: i

    do i = 0, size(values) - 1
      write(*,'(G0.15,1X)', advance='no') values(i)
    end do
    write(*,*)
  end subroutine print_vector_real32


end program main
