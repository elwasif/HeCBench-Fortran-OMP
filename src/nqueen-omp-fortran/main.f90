! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int, c_int8_t
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use omp_lib
  implicit none

  integer, parameter :: queens_block_size = 128
  integer, parameter :: empty = -1
  integer, parameter :: max_depth = 20
  integer, parameter :: max_board = 12
  integer, parameter :: nMaxPrefixos = 75580635

  type, bind(C) :: queen_root
    integer(c_int) :: control
    integer(c_int8_t) :: board(0:max_board - 1)
  end type queen_root

  integer :: size, initial_depth, repeat
  integer(int64) :: tree_size, n_explorers, qtd_sols_global
  type(queen_root), allocatable :: root_prefixes(:)
  integer(int64), allocatable :: vector_of_tree_size(:), solutions(:)
  integer(int64) :: i

  if (command_argument_count() /= 3) then
    print '(A)', 'Usage: ./main <size> <initial depth> <repeat>'
    stop 1
  end if

  size = read_arg(1)
  initial_depth = read_arg(2)
  repeat = read_arg(3)

  write(*,'(/,A,I0,A,I0,A)', advance='no') '### Initial depth: ', initial_depth, ' - Size: ', size, ':'

  if (size > max_depth .or. initial_depth > max_board) then
    print '(A)', ''
    print '(A)', 'Error: unsupported board size or initial depth'
    stop 1
  end if

  allocate(root_prefixes(0:nMaxPrefixos - 1))
  allocate(vector_of_tree_size(0:nMaxPrefixos - 1), solutions(0:nMaxPrefixos - 1))
  root_prefixes%control = 0_c_int
  vector_of_tree_size = 0_int64
  solutions = 0_int64

  call bp_queens_prefixes(size, initial_depth, tree_size, root_prefixes, n_explorers)
  call nqueens(size, initial_depth, int(n_explorers), root_prefixes, vector_of_tree_size, solutions, repeat)

  print '(A,I0)', ''
  print '(A,I0)', 'Tree size: ', tree_size

  qtd_sols_global = 0_int64
  do i = 0, n_explorers - 1
    if (solutions(i) > 0_int64) qtd_sols_global = qtd_sols_global + solutions(i)
    if (vector_of_tree_size(i) > 0_int64) tree_size = tree_size + vector_of_tree_size(i)
  end do

  print '(A,I0,A,/,A,I0)', 'Number of solutions found: ', qtd_sols_global, ' ', 'Tree size: ', tree_size
  if (size == 15 .and. initial_depth == 7) then
    if (qtd_sols_global == 2279184_int64 .and. tree_size == 171129071_int64) then
      print '(A)', 'PASS'
    else
      print '(A)', 'FAIL'
    end if
  end if

  deallocate(root_prefixes, vector_of_tree_size, solutions)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=256) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  logical function host_still_legal(board, r)
    integer(c_int8_t), intent(in) :: board(0:)
    integer, intent(in) :: r
    integer :: i, ld, rd

    host_still_legal = .true.
    do i = 0, r - 1
      if (board(i) == board(r)) then
        host_still_legal = .false.
        return
      end if
    end do

    ld = int(board(r))
    rd = int(board(r))
    do i = r - 1, 0, -1
      ld = ld - 1
      rd = rd + 1
      if (board(i) == ld .or. board(i) == rd) then
        host_still_legal = .false.
        return
      end if
    end do
  end function host_still_legal

  subroutine prefixes_handle_sol(root_prefixes, flag, board, initial_depth, num_sol)
    type(queen_root), intent(inout) :: root_prefixes(0:)
    integer(c_int), intent(in) :: flag
    integer(c_int8_t), intent(in) :: board(0:)
    integer, intent(in) :: initial_depth
    integer(int64), intent(in) :: num_sol
    integer :: i

    root_prefixes(num_sol)%control = flag
    do i = 0, initial_depth - 1
      root_prefixes(num_sol)%board(i) = board(i)
    end do
  end subroutine prefixes_handle_sol

  subroutine bp_queens_prefixes(size, initial_depth, tree_size, root_prefixes, num_sol)
    integer, intent(in) :: size, initial_depth
    integer(int64), intent(out) :: tree_size, num_sol
    type(queen_root), intent(inout) :: root_prefixes(0:)
    integer(c_int) :: flag, bit_test
    integer(c_int8_t) :: vertice(0:max_depth - 1)
    integer :: nivel

    flag = 0_c_int
    tree_size = 0_int64
    num_sol = 0_int64
    vertice = int(empty, c_int8_t)
    nivel = 0

    do while (nivel >= 0)
      vertice(nivel) = vertice(nivel) + 1
      bit_test = shiftl(1_c_int, int(vertice(nivel)))
      if (vertice(nivel) == size) then
        vertice(nivel) = int(empty, c_int8_t)
      else if (host_still_legal(vertice, nivel) .and. (iand(flag, bit_test) == 0_c_int)) then
        flag = ior(flag, bit_test)
        nivel = nivel + 1
        tree_size = tree_size + 1_int64
        if (nivel == initial_depth) then
          call prefixes_handle_sol(root_prefixes, flag, vertice, initial_depth, num_sol)
          num_sol = num_sol + 1_int64
        else
          cycle
        end if
      else
        cycle
      end if

      nivel = nivel - 1
      if (nivel >= 0) flag = iand(flag, not(shiftl(1_c_int, int(vertice(nivel)))))
    end do
  end subroutine bp_queens_prefixes

  subroutine nqueens(size, initial_depth, n_explorers, root_prefixes, vector_of_tree_size, solutions, repeat)
    integer, intent(in) :: size, initial_depth, n_explorers, repeat
    type(queen_root), intent(in) :: root_prefixes(0:)
    integer(int64), intent(out) :: vector_of_tree_size(0:), solutions(0:)
    integer :: rep
    real(real64) :: start_time, end_time

    print '(A)', ''
    print '(A)', '### Regular BP-DFS search. ###'

    !$omp target data map(to: root_prefixes(0:n_explorers - 1)) &
    !$omp& map(from: vector_of_tree_size(0:n_explorers - 1), solutions(0:n_explorers - 1))
    start_time = omp_get_wtime()
    do rep = 1, repeat
      call bp_queens_root_dfs(size, n_explorers, initial_depth, root_prefixes, vector_of_tree_size, solutions)
    end do
    end_time = omp_get_wtime()
    print '(A,F8.6,A)', 'Average kernel execution time: ', (end_time - start_time) / real(repeat, real64), ' (s)'
    !$omp end target data
  end subroutine nqueens

  subroutine bp_queens_root_dfs(size, n_explorers, initial_depth, root_prefixes, vector_of_tree_size, solutions)
    integer, intent(in) :: size, n_explorers, initial_depth
    type(queen_root), intent(in) :: root_prefixes(0:)
    integer(int64), intent(out) :: vector_of_tree_size(0:), solutions(0:)
    integer :: idx, i, depth, depth_global, ld, rd
    integer(c_int) :: flag, bit_test
    integer(int64) :: qtd_solutions_thread, tree_size
    integer(c_int8_t) :: vertice(0:max_depth - 1)
    logical :: safe

    !$omp target teams distribute parallel do thread_limit(queens_block_size) &
    !$omp& private(idx, i, depth, depth_global, ld, rd, flag, bit_test, qtd_solutions_thread, tree_size, vertice, safe)
    do idx = 0, n_explorers - 1
      flag = root_prefixes(idx)%control
      vertice = int(empty, c_int8_t)
      do i = 0, initial_depth - 1
        vertice(i) = root_prefixes(idx)%board(i)
      end do

      depth_global = initial_depth
      depth = depth_global
      qtd_solutions_thread = 0_int64
      tree_size = 0_int64

      do while (depth >= depth_global)
        vertice(depth) = vertice(depth) + 1
        bit_test = shiftl(1_c_int, int(vertice(depth)))
        if (vertice(depth) == size) then
          vertice(depth) = int(empty, c_int8_t)
        else if (iand(flag, bit_test) == 0_c_int) then
          safe = .true.
          do i = 0, depth - 1
            if (vertice(i) == vertice(depth)) safe = .false.
          end do
          ld = int(vertice(depth))
          rd = int(vertice(depth))
          do i = depth - 1, 0, -1
            ld = ld - 1
            rd = rd + 1
            if (vertice(i) == ld .or. vertice(i) == rd) safe = .false.
          end do

          if (safe) then
            tree_size = tree_size + 1_int64
            flag = ior(flag, bit_test)
            depth = depth + 1
            if (depth == size) then
              qtd_solutions_thread = qtd_solutions_thread + 1_int64
            else
              cycle
            end if
          else
            cycle
          end if
        else
          cycle
        end if

        depth = depth - 1
        if (depth >= 0) flag = iand(flag, not(shiftl(1_c_int, int(vertice(depth)))))
      end do

      solutions(idx) = qtd_solutions_thread
      vector_of_tree_size(idx) = tree_size
    end do
    !$omp end target teams distribute parallel do
  end subroutine bp_queens_root_dfs

end program main
