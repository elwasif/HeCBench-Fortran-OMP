! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  integer :: seq_length, num_layers, hidden_size, mini_batch, num_runs
  integer :: num_elements, i, j, m, run, errors
  real(real32), allocatable :: out_i(:), out_h(:), out_c(:)
  real(real32), allocatable :: ref_i(:), ref_h(:), ref_c(:)
  real(real64) :: total_time

  if (command_argument_count() == 5) then
    seq_length = read_arg(1)
    num_layers = read_arg(2)
    hidden_size = read_arg(3)
    mini_batch = read_arg(4)
    num_runs = read_arg(5)
  else if (command_argument_count() == 0) then
    print '(A)', 'Running with default settings'
    seq_length = 100
    num_layers = 4
    hidden_size = 512
    mini_batch = 64
    num_runs = 1
  else
    print '(A)', 'Usage: ./main <seqLength> <numLayers> <hiddenSize> <miniBatch> <repeat>'
    stop 1
  end if

  write(*, '(A,I0,A,I0,A,I0,A,I0)') 'seqLength ', seq_length, ', numLayers ', num_layers, &
      ', hiddenSize ', hidden_size, ', miniBatch ', mini_batch

  num_elements = hidden_size * mini_batch
  allocate(out_i(num_elements * seq_length), out_h(num_elements * num_layers), out_c(num_elements * num_layers))
  allocate(ref_i(num_elements * seq_length), ref_h(num_elements * num_layers), ref_c(num_elements * num_layers))

  total_time = 0.0_real64
  do run = 1, num_runs
    call test_device(hidden_size, mini_batch, seq_length, num_layers, out_i, out_h, out_c, total_time)
    call test_ref(hidden_size, mini_batch, seq_length, num_layers, ref_i, ref_h, ref_c)
  end do

  write(*, '(A,F0.6,A)') 'Average kernel execution time: ', (total_time / real(num_runs, real64)), ' (s)'

  errors = 0
  do m = 0, mini_batch - 1
    do j = 0, seq_length - 1
      do i = 0, hidden_size - 1
        if (abs(out_i(j * num_elements + m * hidden_size + i + 1) - ref_i(j * num_elements + m * hidden_size + i + 1)) > 1.0e-4_real32) errors = errors + 1
      end do
    end do
    do j = 0, num_layers - 1
      do i = 0, hidden_size - 1
        if (abs(out_h(j * num_elements + m * hidden_size + i + 1) - ref_h(j * num_elements + m * hidden_size + i + 1)) > 1.0e-4_real32) errors = errors + 1
        if (abs(out_c(j * num_elements + m * hidden_size + i + 1) - ref_c(j * num_elements + m * hidden_size + i + 1)) > 1.0e-4_real32) errors = errors + 1
      end do
    end do
  end do

  if (errors == 0) then
    print '(A)', 'PASS'
  else
    print '(A)', 'FAIL'
  end if

  deallocate(out_i, out_h, out_c, ref_i, ref_h, ref_c)

contains

  integer function read_arg(position)
    integer, intent(in) :: position
    character(len=128) :: buffer
    call get_command_argument(position, buffer)
    read(buffer, *) read_arg
  end function read_arg

  real(real32) function lcg_value(index, size)
    integer, intent(in) :: index, size
    integer(int64) :: state
    state = int(ieor(index, size), int64)
    state = mod(26757677_int64 * state + 1_int64, 2147483648_int64)
    lcg_value = real(state, real32) / 2147483648.0_real32
  end function lcg_value

  subroutine init_array(data)
    real(real32), intent(out) :: data(:)
    integer :: idx, n
    n = size(data)
    do idx = 1, n
      data(idx) = lcg_value(idx - 1, n)
    end do
  end subroutine init_array

  subroutine init_dev(data, n)
    real(real32), intent(out) :: data(:)
    integer, intent(in) :: n
    integer :: idx

    !$omp target teams distribute parallel do thread_limit(256)
    do idx = 1, n
      data(idx) = lcg_value(idx - 1, n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine init_dev

  subroutine test_ref(hidden_size, mini_batch, seq_length, num_layers, test_i, test_h, test_c)
    integer, intent(in) :: hidden_size, mini_batch, seq_length, num_layers
    real(real32), intent(out) :: test_i(:), test_h(:), test_c(:)
    integer :: num_elements, hc_size, i_size, bias_size, tmp_h_size, tmp_i_size, act_size
    integer :: l_start, l_end, r_start, r_end, recur_batch_size, layer, step
    real(real32), allocatable :: h_data(:), i_data(:), c_data(:), bias(:), tmp_h(:), tmp_i(:), gates(:)

    num_elements = hidden_size * mini_batch
    hc_size = (seq_length + 1) * num_layers * num_elements
    i_size = seq_length * (num_layers + 1) * num_elements
    bias_size = num_layers * hidden_size * 8
    tmp_h_size = 4 * num_layers * num_elements
    tmp_i_size = 4 * seq_length * num_elements
    act_size = 4 * seq_length * num_layers * num_elements
    allocate(h_data(hc_size), i_data(i_size), c_data(hc_size), bias(bias_size), tmp_h(tmp_h_size), tmp_i(tmp_i_size), gates(act_size))
    h_data = 0.0_real32
    i_data = 0.0_real32
    gates = 0.0_real32
    call init_array(tmp_h)
    call init_array(tmp_i)
    call init_array(c_data)
    call init_array(bias)

    l_start = 0
    l_end = 0
    r_start = 0
    recur_batch_size = 2
    do
      call next_tile(l_start, l_end, r_start, r_end, recur_batch_size, seq_length, num_layers)
      if (l_end < 0) exit
      do layer = l_start, l_end - 1
        do step = r_start, r_end - 1
          call elementwise_ref(hidden_size, mini_batch, tmp_h, 4 * layer * num_elements, tmp_i, 4 * step * num_elements, &
              bias, 8 * layer * hidden_size, gates, 4 * (step * num_elements + layer * seq_length * num_elements), &
              h_data, (step + 1) * num_elements + layer * (seq_length + 1) * num_elements, &
              i_data, step * num_elements + (layer + 1) * seq_length * num_elements, &
              c_data, step * num_elements + layer * (seq_length + 1) * num_elements, &
              c_data, (step + 1) * num_elements + layer * (seq_length + 1) * num_elements)
        end do
      end do
    end do

    test_i = i_data(num_layers * seq_length * num_elements + 1 : (num_layers + 1) * seq_length * num_elements)
    do layer = 0, num_layers - 1
      test_h(layer * num_elements + 1 : (layer + 1) * num_elements) = &
          h_data(seq_length * num_elements + layer * (seq_length + 1) * num_elements + 1 : &
                 seq_length * num_elements + layer * (seq_length + 1) * num_elements + num_elements)
      test_c(layer * num_elements + 1 : (layer + 1) * num_elements) = &
          c_data(seq_length * num_elements + layer * (seq_length + 1) * num_elements + 1 : &
                 seq_length * num_elements + layer * (seq_length + 1) * num_elements + num_elements)
    end do
    deallocate(h_data, i_data, c_data, bias, tmp_h, tmp_i, gates)
  end subroutine test_ref

  subroutine test_device(hidden_size, mini_batch, seq_length, num_layers, test_i, test_h, test_c, time)
    integer, intent(in) :: hidden_size, mini_batch, seq_length, num_layers
    real(real32), intent(out) :: test_i(:), test_h(:), test_c(:)
    real(real64), intent(inout) :: time
    integer :: num_elements, hc_size, i_size, bias_size, tmp_h_size, tmp_i_size, act_size
    integer :: l_start, l_end, r_start, r_end, recur_batch_size, layer, step
    real(real64) :: start_time, end_time
    real(real32), allocatable :: h_data(:), i_data(:), c_data(:), bias(:), tmp_h(:), tmp_i(:), gates(:)

    num_elements = hidden_size * mini_batch
    hc_size = (seq_length + 1) * num_layers * num_elements
    i_size = seq_length * (num_layers + 1) * num_elements
    bias_size = num_layers * hidden_size * 8
    tmp_h_size = 4 * num_layers * num_elements
    tmp_i_size = 4 * seq_length * num_elements
    act_size = 4 * seq_length * num_layers * num_elements
    allocate(h_data(hc_size), i_data(i_size), c_data(hc_size), bias(bias_size), tmp_h(tmp_h_size), tmp_i(tmp_i_size), gates(act_size))

    l_start = 0
    l_end = 0
    r_start = 0
    recur_batch_size = 2
    !$omp target data map(alloc: h_data(1:hc_size), i_data(1:i_size), c_data(1:hc_size), bias(1:bias_size), &
    !$omp& tmp_h(1:tmp_h_size), tmp_i(1:tmp_i_size), gates(1:act_size))
    call init_dev(tmp_h, tmp_h_size)
    call init_dev(tmp_i, tmp_i_size)
    call init_dev(c_data, hc_size)
    call init_dev(bias, bias_size)
    do
      call next_tile(l_start, l_end, r_start, r_end, recur_batch_size, seq_length, num_layers)
      if (l_end < 0) exit
      start_time = omp_get_wtime()
      do layer = l_start, l_end - 1
        do step = r_start, r_end - 1
          call elementwise_dev(hidden_size, mini_batch, tmp_h, 4 * layer * num_elements, tmp_i, 4 * step * num_elements, &
              bias, 8 * layer * hidden_size, gates, 4 * (step * num_elements + layer * seq_length * num_elements), &
              h_data, (step + 1) * num_elements + layer * (seq_length + 1) * num_elements, &
              i_data, step * num_elements + (layer + 1) * seq_length * num_elements, &
              c_data, step * num_elements + layer * (seq_length + 1) * num_elements, &
              c_data, (step + 1) * num_elements + layer * (seq_length + 1) * num_elements)
        end do
      end do
      end_time = omp_get_wtime()
      time = time + (end_time - start_time)
    end do
    !$omp target update from(i_data(1:i_size), h_data(1:hc_size), c_data(1:hc_size))
    !$omp end target data

    test_i = i_data(num_layers * seq_length * num_elements + 1 : (num_layers + 1) * seq_length * num_elements)
    do layer = 0, num_layers - 1
      test_h(layer * num_elements + 1 : (layer + 1) * num_elements) = &
          h_data(seq_length * num_elements + layer * (seq_length + 1) * num_elements + 1 : &
                 seq_length * num_elements + layer * (seq_length + 1) * num_elements + num_elements)
      test_c(layer * num_elements + 1 : (layer + 1) * num_elements) = &
          c_data(seq_length * num_elements + layer * (seq_length + 1) * num_elements + 1 : &
                 seq_length * num_elements + layer * (seq_length + 1) * num_elements + num_elements)
    end do
    deallocate(h_data, i_data, c_data, bias, tmp_h, tmp_i, gates)
  end subroutine test_device

  subroutine next_tile(l_start, l_end, r_start, r_end, recur_batch_size, seq_length, num_layers)
    integer, intent(inout) :: l_start, l_end, r_start
    integer, intent(out) :: r_end
    integer, intent(in) :: recur_batch_size, seq_length, num_layers
    if (l_end == 0) then
      l_start = 0
      l_end = 1
      r_start = 0
    else
      l_start = l_start + 1
      l_end = l_end + 1
      r_start = r_start - recur_batch_size
      if (l_end > num_layers .or. r_start < 0) then
        r_start = r_start + (l_start + 1) * recur_batch_size
        l_start = 0
        l_end = 1
      end if
      do while (r_start >= seq_length .and. l_end <= num_layers)
        l_start = l_start + 1
        l_end = l_end + 1
        r_start = r_start - recur_batch_size
      end do
      if (l_end > num_layers .or. r_start < 0) then
        l_end = -1
        r_end = -1
        return
      end if
    end if
    r_end = min(r_start + recur_batch_size, seq_length)
  end subroutine next_tile

  subroutine elementwise_ref(hidden_size, mini_batch, tmp_h, tmp_h_base, tmp_i, tmp_i_base, bias, bias_base, gates, gates_base, &
      h_out, h_base, i_out, i_base, c_in, c_in_base, c_out, c_out_base)
    integer, intent(in) :: hidden_size, mini_batch, tmp_h_base, tmp_i_base, bias_base, gates_base, h_base, i_base, c_in_base, c_out_base
    real(real32), intent(in) :: tmp_h(:), tmp_i(:), bias(:), c_in(:)
    real(real32), intent(inout) :: gates(:), h_out(:), i_out(:), c_out(:)
    integer :: idx, batch, hidden_idx, gate_index, gidx
    real(real32) :: g0, g1, g2, g3, in_gate, forget_gate, in_gate2, out_gate, val
    do idx = 0, mini_batch * hidden_size - 1
      batch = idx / hidden_size
      hidden_idx = mod(idx, hidden_size)
      gate_index = hidden_idx + 4 * batch * hidden_size
      gidx = gate_index + 1
      g0 = tmp_i(tmp_i_base + gidx) + tmp_h(tmp_h_base + gidx) + bias(bias_base + hidden_idx + 1) + bias(bias_base + 4 * hidden_size + hidden_idx + 1)
      g1 = tmp_i(tmp_i_base + hidden_size + gidx) + tmp_h(tmp_h_base + hidden_size + gidx) + bias(bias_base + hidden_size + hidden_idx + 1) + bias(bias_base + 5 * hidden_size + hidden_idx + 1)
      g2 = tmp_i(tmp_i_base + 2 * hidden_size + gidx) + tmp_h(tmp_h_base + 2 * hidden_size + gidx) + bias(bias_base + 2 * hidden_size + hidden_idx + 1) + bias(bias_base + 6 * hidden_size + hidden_idx + 1)
      g3 = tmp_i(tmp_i_base + 3 * hidden_size + gidx) + tmp_h(tmp_h_base + 3 * hidden_size + gidx) + bias(bias_base + 3 * hidden_size + hidden_idx + 1) + bias(bias_base + 7 * hidden_size + hidden_idx + 1)
      gates(gates_base + gidx) = g0
      gates(gates_base + hidden_size + gidx) = g1
      gates(gates_base + 2 * hidden_size + gidx) = g2
      gates(gates_base + 3 * hidden_size + gidx) = g3
      in_gate = 1.0_real32 / (1.0_real32 + exp(-g0))
      forget_gate = 1.0_real32 / (1.0_real32 + exp(-g1))
      in_gate2 = tanh(g2)
      out_gate = 1.0_real32 / (1.0_real32 + exp(-g3))
      val = forget_gate * c_in(c_in_base + idx + 1) + in_gate * in_gate2
      c_out(c_out_base + idx + 1) = val
      val = out_gate * tanh(val)
      h_out(h_base + idx + 1) = val
      i_out(i_base + idx + 1) = val
    end do
  end subroutine elementwise_ref

  subroutine elementwise_dev(hidden_size, mini_batch, tmp_h, tmp_h_base, tmp_i, tmp_i_base, bias, bias_base, gates, gates_base, &
      h_out, h_base, i_out, i_base, c_in, c_in_base, c_out, c_out_base)
    integer, intent(in) :: hidden_size, mini_batch, tmp_h_base, tmp_i_base, bias_base, gates_base, h_base, i_base, c_in_base, c_out_base
    real(real32), intent(in) :: tmp_h(:), tmp_i(:), bias(:), c_in(:)
    real(real32), intent(inout) :: gates(:), h_out(:), i_out(:), c_out(:)
    integer :: idx, batch, hidden_idx, gate_index, gidx
    real(real32) :: g0, g1, g2, g3, in_gate, forget_gate, in_gate2, out_gate, val
    !$omp target teams distribute parallel do thread_limit(256) &
    !$omp& private(batch, hidden_idx, gate_index, gidx, g0, g1, g2, g3, in_gate, forget_gate, in_gate2, out_gate, val)
    do idx = 0, mini_batch * hidden_size - 1
      batch = idx / hidden_size
      hidden_idx = mod(idx, hidden_size)
      gate_index = hidden_idx + 4 * batch * hidden_size
      gidx = gate_index + 1
      g0 = tmp_i(tmp_i_base + gidx) + tmp_h(tmp_h_base + gidx) + bias(bias_base + hidden_idx + 1) + bias(bias_base + 4 * hidden_size + hidden_idx + 1)
      g1 = tmp_i(tmp_i_base + hidden_size + gidx) + tmp_h(tmp_h_base + hidden_size + gidx) + bias(bias_base + hidden_size + hidden_idx + 1) + bias(bias_base + 5 * hidden_size + hidden_idx + 1)
      g2 = tmp_i(tmp_i_base + 2 * hidden_size + gidx) + tmp_h(tmp_h_base + 2 * hidden_size + gidx) + bias(bias_base + 2 * hidden_size + hidden_idx + 1) + bias(bias_base + 6 * hidden_size + hidden_idx + 1)
      g3 = tmp_i(tmp_i_base + 3 * hidden_size + gidx) + tmp_h(tmp_h_base + 3 * hidden_size + gidx) + bias(bias_base + 3 * hidden_size + hidden_idx + 1) + bias(bias_base + 7 * hidden_size + hidden_idx + 1)
      gates(gates_base + gidx) = g0
      gates(gates_base + hidden_size + gidx) = g1
      gates(gates_base + 2 * hidden_size + gidx) = g2
      gates(gates_base + 3 * hidden_size + gidx) = g3
      in_gate = 1.0_real32 / (1.0_real32 + exp(-g0))
      forget_gate = 1.0_real32 / (1.0_real32 + exp(-g1))
      in_gate2 = tanh(g2)
      out_gate = 1.0_real32 / (1.0_real32 + exp(-g3))
      val = forget_gate * c_in(c_in_base + idx + 1) + in_gate * in_gate2
      c_out(c_out_base + idx + 1) = val
      val = out_gate * tanh(val)
      h_out(h_base + idx + 1) = val
      i_out(i_base + idx + 1) = val
    end do
    !$omp end target teams distribute parallel do
  end subroutine elementwise_dev

end program main
