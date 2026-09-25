! SPDX-License-Identifier: CC0-1.0
module grep_port
  use iso_fortran_env, only: int32, int64, real64
  use omp_lib
  implicit none

  integer, parameter :: any_tok = int(z'15')
  integer, parameter :: concatenate_tok = int(z'1b')
  integer, parameter :: alternate_tok = int(z'04')
  integer, parameter :: question_tok = int(z'02')
  integer, parameter :: star_tok = int(z'03')
  integer, parameter :: plus_tok = int(z'01')
  integer, parameter :: paren_open_tok = int(z'05')
  integer, parameter :: paren_close_tok = int(z'06')
  integer, parameter :: match_state = 256
  integer, parameter :: split_state = 257
  integer, parameter :: any_state = 258
  integer, parameter :: buffer_size = 8000
  integer, parameter :: max_states = 100
  integer, parameter :: max_stack = 1000
  integer, parameter :: max_teams = 512
  integer, parameter :: max_threads = 160

contains

  subroutine usage(prog)
    character(len=*), intent(in) :: prog
    print '(a)', 'Usage: '//trim(prog)//' [options] [pattern] '
    print '(a)', 'Program Options:'
    print '(a)', '  -v  Visualize the NFA then exit'
    print '(a)', '  -p  View postfix expression then exit'
    print '(a)', '  -s  View simplified expression then exit'
    print '(a)', '  -t  Print timing data'
    print '(a)', '  -f <FILE> --file Input file to be matched'
    print '(a)', '  -r <FILE> --regex Input file with regexs'
    print '(a)', '  -? This message'
    print '(a)', '[pattern] required only if -r or --regex is not used'
  end subroutine usage

  subroutine parse_args(file_name, pattern, timer_on, simplified, postfix, visualize, have_pattern)
    character(len=:), allocatable, intent(out) :: file_name, pattern
    logical, intent(out) :: timer_on, simplified, postfix, visualize, have_pattern
    character(len=4096) :: arg
    integer :: argc, i, n

    argc = command_argument_count()
    timer_on = .false.
    simplified = .false.
    postfix = .false.
    visualize = .false.
    have_pattern = .false.
    file_name = ''
    pattern = ''

    if (argc < 2) then
      call get_command_argument(0, arg, length=n)
      call usage(arg(:n))
      stop
    end if

    i = 1
    do while (i <= argc)
      call get_command_argument(i, arg, length=n)
      select case (arg(:n))
      case ('-v', '--visualize')
        visualize = .true.
      case ('-p', '--postfix')
        postfix = .true.
      case ('-s', '--simplified')
        simplified = .true.
      case ('-t', '--time')
        timer_on = .true.
      case ('-?', '--help')
        call get_command_argument(0, arg, length=n)
        call usage(arg(:n))
        stop
      case ('-f', '--file')
        i = i + 1
        if (i <= argc) then
          call get_command_argument(i, arg, length=n)
          file_name = arg(:n)
        end if
      case ('-r', '--regex')
        i = i + 1
      case default
        pattern = arg(:n)
        have_pattern = .true.
      end select
      i = i + 1
    end do
  end subroutine parse_args

  subroutine append_value(values, n, cap, value)
    integer, allocatable, intent(inout) :: values(:)
    integer, intent(inout) :: n, cap
    integer, intent(in) :: value
    integer, allocatable :: tmp(:)

    if (n >= cap) then
      allocate(tmp(max(16, cap * 2)))
      if (n > 0) tmp(1:n) = values(1:n)
      call move_alloc(tmp, values)
      cap = size(values)
    end if
    n = n + 1
    values(n) = value
  end subroutine append_value

  subroutine append_range(values, n, cap, first_char, last_char)
    integer, allocatable, intent(inout) :: values(:)
    integer, intent(inout) :: n, cap
    integer, intent(in) :: first_char, last_char
    integer :: c

    call append_value(values, n, cap, paren_open_tok)
    call append_value(values, n, cap, first_char)
    do c = first_char + 1, last_char
      call append_value(values, n, cap, alternate_tok)
      call append_value(values, n, cap, c)
    end do
    call append_value(values, n, cap, paren_close_tok)
  end subroutine append_range

  subroutine simplify_regex(pattern, regex)
    character(len=*), intent(in) :: pattern
    integer, allocatable, intent(out) :: regex(:)
    integer, allocatable :: values(:)
    integer :: i, n, cap, lenp, c1, c2

    lenp = len_trim(pattern)
    cap = max(16, lenp + 3)
    allocate(values(cap))
    n = 0
    i = 1
    do while (i <= lenp)
      select case (pattern(i:i))
      case ('\')
        if (i == lenp) then
          call append_value(values, n, cap, iachar('\'))
        else
          i = i + 1
          select case (pattern(i:i))
          case ('t')
            call append_value(values, n, cap, 9)
          case ('n')
            call append_value(values, n, cap, 10)
          case ('d')
            call append_range(values, n, cap, iachar('0'), iachar('9'))
          case ('w')
            call append_value(values, n, cap, paren_open_tok)
            call append_range(values, n, cap, iachar('a'), iachar('z'))
            call append_value(values, n, cap, alternate_tok)
            call append_range(values, n, cap, iachar('A'), iachar('Z'))
            call append_value(values, n, cap, alternate_tok)
            call append_value(values, n, cap, iachar('_'))
            call append_value(values, n, cap, paren_close_tok)
          case ('s')
            call append_value(values, n, cap, paren_open_tok)
            call append_value(values, n, cap, iachar(' '))
            call append_value(values, n, cap, alternate_tok)
            call append_value(values, n, cap, 9)
            call append_value(values, n, cap, alternate_tok)
            call append_value(values, n, cap, 10)
            call append_value(values, n, cap, paren_close_tok)
          case default
            call append_value(values, n, cap, iachar(pattern(i:i)))
          end select
        end if
      case ('.')
        call append_value(values, n, cap, any_tok)
      case ('+')
        call append_value(values, n, cap, plus_tok)
      case ('?')
        call append_value(values, n, cap, question_tok)
      case ('*')
        call append_value(values, n, cap, star_tok)
      case ('|')
        call append_value(values, n, cap, alternate_tok)
      case ('(')
        call append_value(values, n, cap, paren_open_tok)
      case (')')
        call append_value(values, n, cap, paren_close_tok)
      case ('[')
        if (i + 4 <= lenp .and. pattern(i+2:i+2) == '-' .and. pattern(i+4:i+4) == ']') then
          c1 = iachar(pattern(i+1:i+1))
          c2 = iachar(pattern(i+3:i+3))
          if (c1 <= c2 .and. c1 > 32) then
            call append_range(values, n, cap, c1, c2)
            i = i + 4
          else
            call append_value(values, n, cap, iachar(pattern(i:i)))
          end if
        else
          call append_value(values, n, cap, iachar(pattern(i:i)))
        end if
      case default
        call append_value(values, n, cap, iachar(pattern(i:i)))
      end select
      i = i + 1
    end do

    allocate(regex(n + 1))
    if (n > 0) regex(1:n) = values(1:n)
    regex(n + 1) = 0
  end subroutine simplify_regex

  subroutine print_regex(label, regex)
    character(len=*), intent(in) :: label
    integer, intent(in) :: regex(:)
    character(len=:), allocatable :: out
    integer :: i, n

    n = 0
    do while (n < size(regex) .and. regex(n + 1) /= 0)
      n = n + 1
    end do
    allocate(character(len=n) :: out)
    do i = 1, n
      select case (regex(i))
      case (any_tok)
        out(i:i) = '.'
      case (concatenate_tok)
        out(i:i) = '`'
      case (alternate_tok)
        out(i:i) = '|'
      case (question_tok)
        out(i:i) = '?'
      case (star_tok)
        out(i:i) = '*'
      case (plus_tok)
        out(i:i) = '+'
      case (paren_open_tok)
        out(i:i) = '('
      case (paren_close_tok)
        out(i:i) = ')'
      case default
        out(i:i) = achar(regex(i))
      end select
    end do
    print '(a)'
    print '(a)', trim(label)//out
  end subroutine print_regex

  logical function pre2post_host(re, post) result(ok)
    integer, intent(in) :: re(:)
    integer, allocatable, intent(out) :: post(:)
    integer :: tmp(buffer_size)
    integer :: paren_alt(100), paren_atom(100)
    integer :: nalt, natom, p, dst, i, token, lenr

    lenr = 0
    do while (lenr < size(re) .and. re(lenr + 1) /= 0)
      lenr = lenr + 1
    end do
    ok = .false.
    if (lenr >= buffer_size / 2) return
    nalt = 0
    natom = 0
    p = 1
    dst = 1
    do i = 1, lenr
      token = re(i)
      select case (token)
      case (paren_open_tok)
        if (natom > 1) then
          natom = natom - 1
          tmp(dst) = concatenate_tok
          dst = dst + 1
        end if
        if (p > 100) return
        paren_alt(p) = nalt
        paren_atom(p) = natom
        p = p + 1
        nalt = 0
        natom = 0
      case (alternate_tok)
        if (natom == 0) return
        do
          natom = natom - 1
          if (natom <= 0) exit
          tmp(dst) = concatenate_tok
          dst = dst + 1
        end do
        nalt = nalt + 1
      case (paren_close_tok)
        if (p == 1 .or. natom == 0) return
        do
          natom = natom - 1
          if (natom <= 0) exit
          tmp(dst) = concatenate_tok
          dst = dst + 1
        end do
        do while (nalt > 0)
          tmp(dst) = alternate_tok
          dst = dst + 1
          nalt = nalt - 1
        end do
        p = p - 1
        nalt = paren_alt(p)
        natom = paren_atom(p) + 1
      case (star_tok, plus_tok, question_tok)
        if (natom == 0) return
        tmp(dst) = token
        dst = dst + 1
      case default
        if (natom > 1) then
          natom = natom - 1
          tmp(dst) = concatenate_tok
          dst = dst + 1
        end if
        tmp(dst) = token
        dst = dst + 1
        natom = natom + 1
      end select
    end do
    if (p /= 1) return
    do
      natom = natom - 1
      if (natom <= 0) exit
      tmp(dst) = concatenate_tok
      dst = dst + 1
    end do
    do while (nalt > 0)
      tmp(dst) = alternate_tok
      dst = dst + 1
      nalt = nalt - 1
    end do
    tmp(dst) = 0
    allocate(post(dst))
    post(1:dst) = tmp(1:dst)
    ok = .true.
  end function pre2post_host

  subroutine read_file_bytes(file_name, bytes, table, num_lines, file_len, end_read_file)
    character(len=*), intent(in) :: file_name
    integer, allocatable, intent(out) :: bytes(:), table(:)
    integer, intent(out) :: num_lines, file_len
    real(real64), intent(out) :: end_read_file
    character(len=:), allocatable :: text
    integer :: unit, stat, i
    integer(int64) :: fsize

    inquire(file=trim(file_name), size=fsize)
    if (fsize < 0_int64) then
      allocate(bytes(0:0), table(0:0))
      bytes = 0
      table = 0
      num_lines = 0
      file_len = 0
      end_read_file = omp_get_wtime()
      return
    end if
    file_len = int(fsize)
    allocate(character(len=file_len) :: text)
    open(newunit=unit, file=trim(file_name), access='stream', form='unformatted', status='old', action='read', iostat=stat)
    if (stat /= 0) stop 'Error opening file'
    if (file_len > 0) read(unit) text
    close(unit)
    end_read_file = omp_get_wtime()

    allocate(bytes(0:file_len))
    allocate(table(0:max(0, file_len)))
    do i = 1, file_len
      bytes(i - 1) = iachar(text(i:i))
    end do
    bytes(file_len) = 0

    table(0) = 0
    num_lines = 0
    do i = 0, file_len - 1
      if (bytes(i) == 10) then
        num_lines = num_lines + 1
        table(num_lines) = i + 1
        bytes(i) = 0
      end if
    end do
    if (file_len > 0) then
      if (bytes(file_len - 1) == 0) num_lines = num_lines - 1
    end if
  end subroutine read_file_bytes

  subroutine run_kernel(regex, postsize, table, line_data, file_len, num_lines, result)
    integer, intent(in) :: regex(:)
    integer, intent(in) :: postsize, file_len, num_lines
    integer, intent(in) :: table(0:), line_data(0:)
    integer, intent(out) :: result(0:)
    integer, allocatable :: buf(:, :), state_c(:, :), state_out(:, :), state_out1(:, :)
    integer, allocatable :: list_a(:, :, :), list_b(:, :, :)
    integer, allocatable :: starts(:), nstates(:)
    integer :: team, tid, stride, idx, start_state, n1, n2
    logical :: ok

    if (num_lines <= 0) return
    allocate(buf(max_teams, buffer_size))
    allocate(state_c(max_teams, 0:max_states))
    allocate(state_out(max_teams, 0:max_states))
    allocate(state_out1(max_teams, 0:max_states))
    allocate(list_a(max_teams, max_threads, max_states))
    allocate(list_b(max_teams, max_threads, max_states))
    allocate(starts(max_teams), nstates(max_teams))

    !$omp target data map(to: regex(1:postsize), table(0:file_len), line_data(0:file_len)) &
    !$omp& map(from: result(0:num_lines - 1)) &
    !$omp& map(alloc: buf(1:max_teams, 1:buffer_size), state_c(1:max_teams, 0:max_states), &
    !$omp& state_out(1:max_teams, 0:max_states), state_out1(1:max_teams, 0:max_states), &
    !$omp& list_a(1:max_teams, 1:max_threads, 1:max_states), &
    !$omp& list_b(1:max_teams, 1:max_threads, 1:max_states), &
    !$omp& starts(1:max_teams), nstates(1:max_teams))
    !$omp target teams num_teams(max_teams) thread_limit(max_threads) private(team)
    team = omp_get_team_num() + 1
    !$omp parallel private(tid, stride, idx, start_state, n1, n2, ok)
    tid = omp_get_thread_num()
    if (tid == 0) then
      call pre2post_device(regex, buf(team, :), ok)
      nstates(team) = 0
      starts(team) = 0
      if (ok) call post2nfa_device(buf(team, :), state_c(team, :), state_out(team, :), &
          state_out1(team, :), nstates(team), starts(team))
    end if
    !$omp barrier
    start_state = starts(team)
    stride = omp_get_num_threads() * omp_get_num_teams()
    idx = omp_get_team_num() * omp_get_num_threads() + tid
    do while (idx < num_lines)
      if (start_state > 0) then
        if (match_line_device(start_state, table(idx), line_data, state_c(team, :), &
            state_out(team, :), state_out1(team, :), list_a(team, tid + 1, :), &
            list_b(team, tid + 1, :), n1, n2)) then
          result(idx) = 1
        else
          result(idx) = 0
        end if
      else
        result(idx) = 0
      end if
      idx = idx + stride
    end do
      !$omp end parallel
    !$omp end target teams
    !$omp end target data

    deallocate(buf, state_c, state_out, state_out1, list_a, list_b, starts, nstates)
  end subroutine run_kernel

  subroutine pre2post_device(re, dst, ok)
    !$omp declare target
    integer, intent(in) :: re(:)
    integer, intent(out) :: dst(:)
    logical, intent(out) :: ok
    integer :: paren_alt(100), paren_atom(100)
    integer :: nalt, natom, p, outp, i, token, lenr

    lenr = 0
    do while (lenr < size(re) .and. re(lenr + 1) /= 0)
      lenr = lenr + 1
    end do
    ok = .false.
    if (lenr >= buffer_size / 2) return
    nalt = 0
    natom = 0
    p = 1
    outp = 1
    do i = 1, lenr
      token = re(i)
      select case (token)
      case (paren_open_tok)
        if (natom > 1) then
          natom = natom - 1
          dst(outp) = concatenate_tok
          outp = outp + 1
        end if
        if (p > 100) return
        paren_alt(p) = nalt
        paren_atom(p) = natom
        p = p + 1
        nalt = 0
        natom = 0
      case (alternate_tok)
        if (natom == 0) return
        do
          natom = natom - 1
          if (natom <= 0) exit
          dst(outp) = concatenate_tok
          outp = outp + 1
        end do
        nalt = nalt + 1
      case (paren_close_tok)
        if (p == 1 .or. natom == 0) return
        do
          natom = natom - 1
          if (natom <= 0) exit
          dst(outp) = concatenate_tok
          outp = outp + 1
        end do
        do while (nalt > 0)
          dst(outp) = alternate_tok
          outp = outp + 1
          nalt = nalt - 1
        end do
        p = p - 1
        nalt = paren_alt(p)
        natom = paren_atom(p) + 1
      case (star_tok, plus_tok, question_tok)
        if (natom == 0) return
        dst(outp) = token
        outp = outp + 1
      case default
        if (natom > 1) then
          natom = natom - 1
          dst(outp) = concatenate_tok
          outp = outp + 1
        end if
        dst(outp) = token
        outp = outp + 1
        natom = natom + 1
      end select
    end do
    if (p /= 1) return
    do
      natom = natom - 1
      if (natom <= 0) exit
      dst(outp) = concatenate_tok
      outp = outp + 1
    end do
    do while (nalt > 0)
      dst(outp) = alternate_tok
      outp = outp + 1
      nalt = nalt - 1
    end do
    dst(outp) = 0
    ok = .true.
  end subroutine pre2post_device

  subroutine post2nfa_device(post, state_c, state_out, state_out1, nstate, start_state)
    !$omp declare target
    integer, intent(in) :: post(:)
    integer, intent(inout) :: state_c(0:), state_out(0:), state_out1(0:)
    integer, intent(out) :: nstate, start_state
    integer :: frag_start(max_stack), frag_head(max_stack)
    integer :: pend_state(max_stack), pend_slot(max_stack), pend_next(max_stack)
    integer :: stackp, pendp, i, token, s, e1_start, e1_head, e2_start, e2_head, e_start, e_head

    nstate = 0
    stackp = 0
    pendp = 0
    i = 1
    do while (i <= size(post) .and. post(i) /= 0)
      token = post(i)
      select case (token)
      case (any_tok)
        call new_state(any_state, 0, 0, state_c, state_out, state_out1, nstate, s)
        call list1_pending(s, 1, pend_state, pend_slot, pend_next, pendp, e_head)
        call push_frag(s, e_head, frag_start, frag_head, stackp)
      case (concatenate_tok)
        call pop_frag(e2_start, e2_head, frag_start, frag_head, stackp)
        call pop_frag(e1_start, e1_head, frag_start, frag_head, stackp)
        call patch_pending(e1_head, e2_start, pend_state, pend_slot, pend_next, state_out, state_out1)
        call push_frag(e1_start, e2_head, frag_start, frag_head, stackp)
      case (alternate_tok)
        call pop_frag(e2_start, e2_head, frag_start, frag_head, stackp)
        call pop_frag(e1_start, e1_head, frag_start, frag_head, stackp)
        call new_state(split_state, e1_start, e2_start, state_c, state_out, state_out1, nstate, s)
        e_head = append_pending(e1_head, e2_head, pend_next)
        call push_frag(s, e_head, frag_start, frag_head, stackp)
      case (question_tok)
        call pop_frag(e_start, e_head, frag_start, frag_head, stackp)
        call new_state(split_state, e_start, 0, state_c, state_out, state_out1, nstate, s)
        call list1_pending(s, 2, pend_state, pend_slot, pend_next, pendp, e2_head)
        e_head = append_pending(e_head, e2_head, pend_next)
        call push_frag(s, e_head, frag_start, frag_head, stackp)
      case (star_tok)
        call pop_frag(e_start, e_head, frag_start, frag_head, stackp)
        call new_state(split_state, e_start, 0, state_c, state_out, state_out1, nstate, s)
        call patch_pending(e_head, s, pend_state, pend_slot, pend_next, state_out, state_out1)
        call list1_pending(s, 2, pend_state, pend_slot, pend_next, pendp, e2_head)
        call push_frag(s, e2_head, frag_start, frag_head, stackp)
      case (plus_tok)
        call pop_frag(e_start, e_head, frag_start, frag_head, stackp)
        call new_state(split_state, e_start, 0, state_c, state_out, state_out1, nstate, s)
        call patch_pending(e_head, s, pend_state, pend_slot, pend_next, state_out, state_out1)
        call list1_pending(s, 2, pend_state, pend_slot, pend_next, pendp, e2_head)
        call push_frag(e_start, e2_head, frag_start, frag_head, stackp)
      case default
        call new_state(token, 0, 0, state_c, state_out, state_out1, nstate, s)
        call list1_pending(s, 1, pend_state, pend_slot, pend_next, pendp, e_head)
        call push_frag(s, e_head, frag_start, frag_head, stackp)
      end select
      i = i + 1
    end do

    call pop_frag(e_start, e_head, frag_start, frag_head, stackp)
    state_c(0) = match_state
    state_out(0) = 0
    state_out1(0) = 0
    call patch_pending(e_head, 0, pend_state, pend_slot, pend_next, state_out, state_out1)
    start_state = e_start
  end subroutine post2nfa_device

  subroutine new_state(c, out, out1, state_c, state_out, state_out1, nstate, s)
    !$omp declare target
    integer, intent(in) :: c, out, out1
    integer, intent(inout) :: state_c(0:), state_out(0:), state_out1(0:)
    integer, intent(inout) :: nstate
    integer, intent(out) :: s

    nstate = nstate + 1
    if (nstate > max_states) then
      s = 0
      return
    end if
    s = nstate
    state_c(s) = c
    state_out(s) = out
    state_out1(s) = out1
  end subroutine new_state

  subroutine push_frag(start, head, frag_start, frag_head, stackp)
    !$omp declare target
    integer, intent(in) :: start, head
    integer, intent(inout) :: frag_start(:), frag_head(:), stackp
    stackp = stackp + 1
    frag_start(stackp) = start
    frag_head(stackp) = head
  end subroutine push_frag

  subroutine pop_frag(start, head, frag_start, frag_head, stackp)
    !$omp declare target
    integer, intent(out) :: start, head
    integer, intent(in) :: frag_start(:), frag_head(:)
    integer, intent(inout) :: stackp
    start = frag_start(stackp)
    head = frag_head(stackp)
    stackp = stackp - 1
  end subroutine pop_frag

  subroutine list1_pending(state, slot, pend_state, pend_slot, pend_next, pendp, head)
    !$omp declare target
    integer, intent(in) :: state, slot
    integer, intent(inout) :: pend_state(:), pend_slot(:), pend_next(:), pendp
    integer, intent(out) :: head
    pendp = pendp + 1
    head = pendp
    pend_state(head) = state
    pend_slot(head) = slot
    pend_next(head) = 0
  end subroutine list1_pending

  integer function append_pending(head1, head2, pend_next) result(head)
    !$omp declare target
    integer, intent(in) :: head1, head2
    integer, intent(inout) :: pend_next(:)
    integer :: p
    if (head1 == 0) then
      head = head2
      return
    end if
    p = head1
    do while (pend_next(p) /= 0)
      p = pend_next(p)
    end do
    pend_next(p) = head2
    head = head1
  end function append_pending

  subroutine patch_pending(head, target, pend_state, pend_slot, pend_next, state_out, state_out1)
    !$omp declare target
    integer, intent(in) :: head, target
    integer, intent(in) :: pend_state(:), pend_slot(:), pend_next(:)
    integer, intent(inout) :: state_out(0:), state_out1(0:)
    integer :: p, s

    p = head
    do while (p /= 0)
      s = pend_state(p)
      if (pend_slot(p) == 1) then
        state_out(s) = target
      else
        state_out1(s) = target
      end if
      p = pend_next(p)
    end do
  end subroutine patch_pending

  logical function match_line_device(start_state, offset, line_data, state_c, state_out, state_out1, &
      list1, list2, n1, n2) result(matched)
    !$omp declare target
    integer, intent(in) :: start_state, offset
    integer, intent(in) :: line_data(0:), state_c(0:), state_out(0:), state_out1(0:)
    integer, intent(inout) :: list1(:), list2(:), n1, n2
    integer :: pos, c, tmpn, i, cur

    n1 = 0
    call add_state_device(list1, n1, start_state, state_c, state_out, state_out1)
    n2 = 0
    pos = offset
    do while (line_data(pos) /= 0)
      c = iand(line_data(pos), 255)
      n2 = 0
      do i = 1, n1
        cur = list1(i)
        if (state_c(cur) == c .or. state_c(cur) == any_state) then
          call add_state_device(list2, n2, state_out(cur), state_c, state_out, state_out1)
        end if
      end do
      do i = 1, n2
        list1(i) = list2(i)
      end do
      tmpn = n1
      n1 = n2
      n2 = tmpn
      pos = pos + 1
    end do
    matched = .false.
    do i = 1, n1
      if (state_c(list1(i)) == match_state) then
        matched = .true.
        return
      end if
    end do
  end function match_line_device

  subroutine add_state_device(list, n, start, state_c, state_out, state_out1)
    !$omp declare target
    integer, intent(inout) :: list(:), n
    integer, intent(in) :: start
    integer, intent(in) :: state_c(0:), state_out(0:), state_out1(0:)
    integer :: stack(max_states), sp, s

    sp = 0
    if (start >= 0) then
      sp = sp + 1
      stack(sp) = start
    end if
    do while (sp > 0)
      s = stack(sp)
      sp = sp - 1
      if (s < 0) cycle
      if (state_c(s) == split_state) then
        if (state_out(s) >= 0 .and. sp < max_states) then
          sp = sp + 1
          stack(sp) = state_out(s)
        end if
        if (state_out1(s) >= 0 .and. sp < max_states) then
          sp = sp + 1
          stack(sp) = state_out1(s)
        end if
      else
        if (n < max_states) then
          n = n + 1
          list(n) = s
        end if
      end if
    end do
  end subroutine add_state_device

  subroutine visualize_nfa(post)
    integer, intent(in) :: post(:)
    integer :: state_c(0:max_states), state_out(0:max_states), state_out1(0:max_states)
    integer :: count(0:max_states), nstate, start_state

    state_c = 0
    state_out = 0
    state_out1 = 0
    count = 0
    call post2nfa_device(post, state_c, state_out, state_out1, nstate, start_state)
    write(*, '(a)', advance='no') '['
    call visualize_nfa_help(start_state, state_c, state_out, state_out1, count)
    write(*, '(a)') ']'
  end subroutine visualize_nfa

  recursive subroutine visualize_nfa_help(start, state_c, state_out, state_out1, count)
    integer, intent(in) :: start
    integer, intent(in) :: state_c(0:), state_out(0:), state_out1(0:)
    integer, intent(inout) :: count(0:)
    integer :: out_id, out1_id
    character(len=16) :: data

    if (start < 0 .or. start > max_states) return
    if (count(start) > 0) return
    count(start) = count(start) + 1

    select case (state_c(start))
    case (match_state)
      data = 'Match'
      out_id = -1
      out1_id = -1
    case (split_state)
      data = 'Split'
      out_id = state_out(start)
      out1_id = state_out1(start)
    case (any_state)
      data = 'Any'
      out_id = state_out(start)
      out1_id = -1
    case default
      data = 'Char '//achar(iand(state_c(start), 255))
      out_id = state_out(start)
      out1_id = -1
    end select

    write(*, '("{ ""id"": """, i0, """, ""data"":""", a, """, ""out"":""", i0, """, ""out1"":""", i0, """ ")') &
        start, trim(data), out_id, out1_id
    write(*, '(a)', advance='no') '},'

    call visualize_nfa_help(state_out(start), state_c, state_out, state_out1, count)
    call visualize_nfa_help(state_out1(start), state_c, state_out, state_out1, count)
  end subroutine visualize_nfa_help

end module grep_port

program main
  use grep_port
  implicit none

  character(len=:), allocatable :: file_name, pattern
  integer, allocatable :: regex(:), post(:), line_data(:), table(:), result(:)
  logical :: timer_on, simplified, postfix, visualize, have_pattern, ok
  integer :: num_lines, file_len, i
  real(real64) :: start_time, end_read_file, end_setup, end_kernel, end_time

  call parse_args(file_name, pattern, timer_on, simplified, postfix, visualize, have_pattern)
  if (.not. have_pattern) then
    call usage('./main')
    stop
  end if

  call simplify_regex(pattern, regex)
  ok = pre2post_host(regex, post)
  if (.not. ok) then
    write(*, '(a)') 'bad regexp '//pattern
    stop 1
  end if
  if (simplified) then
    call print_regex('Simplified Regex: ', regex)
    stop
  end if
  if (postfix) then
    call print_regex('Postfix buffer: ', post)
    stop
  end if
  if (visualize) then
    call visualize_nfa(post)
    stop
  end if

  if (len_trim(file_name) == 0) then
    print '(a)', 'Enter a file '
    stop
  end if

  start_time = omp_get_wtime()
  call read_file_bytes(file_name, line_data, table, num_lines, file_len, end_read_file)

  allocate(result(0:max(0, num_lines - 1)))
  result = 0
  end_setup = omp_get_wtime()
  call run_kernel(regex, size(regex), table, line_data, file_len, num_lines, result)
  end_kernel = omp_get_wtime()

  if (.not. timer_on) then
    do i = 0, num_lines - 1
      if (result(i) == 1) call print_line(line_data, table(i))
    end do
  end if

  end_time = omp_get_wtime()
  if (timer_on) then
    print '(a,f0.4,1x)', ''
    call print_timing('ReadFile time ', end_read_file - start_time)
    print '(a)'
    call print_timing('Device setup time ', end_setup - end_read_file)
    print '(a)'
    call print_timing('Kernel execution Time ', end_kernel - end_setup)
    print '(a)'
    call print_timing('Total time ', end_time - start_time)
    print '(a)'
  end if

contains

  subroutine print_timing(label, value)
    character(len=*), intent(in) :: label
    real(real64), intent(in) :: value
    character(len=64) :: text

    write(text, '(f0.4)') value
    text = adjustl(text)
    if (text(1:1) == '.') text = '0'//trim(text)
    if (len_trim(text) >= 2 .and. text(1:2) == '-.') text = '-0'//trim(text(2:))
    write(*, '(a,a," ")') label, trim(text)
  end subroutine print_timing

  subroutine print_line(bytes, offset)
    integer, intent(in) :: bytes(0:), offset
    integer :: pos, n
    character(len=:), allocatable :: line

    n = 0
    do while (bytes(offset + n) /= 0)
      n = n + 1
    end do
    allocate(character(len=n) :: line)
    do pos = 1, n
      line(pos:pos) = achar(bytes(offset + pos - 1))
    end do
    print '(a)', line
  end subroutine print_line

end program main
