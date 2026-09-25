! SPDX-License-Identifier: CC0-1.0
program main
  use, intrinsic :: iso_c_binding, only : c_int
  use, intrinsic :: iso_fortran_env, only : int64, real32, real64
  use omp_lib
  implicit none

  real(real64), parameter :: eps_0 = 8.85418782e-12_real64
  real(real64), parameter :: k_b = 1.38065e-23_real64
  real(real64), parameter :: me = 9.10938215e-31_real64
  real(real64), parameter :: qe = 1.602176565e-19_real64
  real(real64), parameter :: amu = 1.660538921e-27_real64
  real(real64), parameter :: ev_to_k = 11604.52_real64
  real(real64), parameter :: plasma_den = 1.0e16_real64
  integer, parameter :: num_ions = 500000
  integer, parameter :: num_electrons = 500000
  real(real64), parameter :: dx = 1.0e-4_real64
  integer, parameter :: nc = 100
  integer, parameter :: num_ts = 1000
  real(real64), parameter :: dt = 1.0e-11_real64
  real(real64), parameter :: electron_temp = 3.0_real64
  real(real64), parameter :: ion_temp = 1.0_real64
  real(real64), parameter :: x0 = 0.0_real64
  real(real64), parameter :: xl = real(nc, real64) * dx
  real(real64), parameter :: xmax = x0 + xl
  integer, parameter :: ni = nc + 1
  integer, parameter :: threads_per_block = 256
  integer, parameter :: particle_x = 0
  integer, parameter :: particle_v = 1
  integer, parameter :: particle_alive = 2
  integer, parameter :: particle_stride = 3

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

  real(real64), allocatable :: phi(:), rho(:), ef(:), phi_direct(:)
  real(real32), allocatable :: ndi(:), nde(:), ndi_ref(:), nde_ref(:)
  real(real64), allocatable :: ions_part(:), electrons_part(:)
  real(real64) :: ions_mass, ions_charge, ions_spwt
  real(real64) :: electrons_mass, electrons_charge, electrons_spwt
  real(real64) :: delta_ions, delta_electrons, v_thi, v_the
  real(real64) :: sp_time, start_time, end_time, total_time, max_phi
  integer :: p, ts, ios

  allocate(phi(0:ni - 1), rho(0:ni - 1), ef(0:ni - 1), phi_direct(0:ni - 1))
  allocate(ndi(0:ni - 1), nde(0:ni - 1), ndi_ref(0:ni - 1), nde_ref(0:ni - 1))
  allocate(ions_part(0:particle_stride * num_ions - 1), &
    electrons_part(0:particle_stride * num_electrons - 1))

  phi = 0.0_real64
  rho = 0.0_real64
  ef = 0.0_real64
  ndi = 0.0_real32
  nde = 0.0_real32
  sp_time = 0.0_real64

  ions_mass = 16.0_real64 * amu
  ions_charge = qe
  ions_spwt = plasma_den * xl / real(num_ions, real64)
  electrons_mass = me
  electrons_charge = -qe
  electrons_spwt = plasma_den * xl / real(num_electrons, real64)

  call c_srand(123_c_int)
  delta_ions = xl / real(num_ions, real64)
  v_thi = sqrt(2.0_real64 * k_b * ion_temp * ev_to_k / ions_mass)
  do p = 0, num_ions - 1
    ions_part(particle_stride * p + particle_x) = x0 + real(p, real64) * delta_ions
    ions_part(particle_stride * p + particle_v) = sample_vel(v_thi)
    ions_part(particle_stride * p + particle_alive) = 1.0_real64
  end do

  delta_electrons = xl / real(num_electrons, real64)
  v_the = sqrt(2.0_real64 * k_b * electron_temp * ev_to_k / electrons_mass)
  do p = 0, num_electrons - 1
    electrons_part(particle_stride * p + particle_x) = x0 + real(p, real64) * delta_electrons
    electrons_part(particle_stride * p + particle_v) = sample_vel(v_the)
    electrons_part(particle_stride * p + particle_alive) = 1.0_real64
  end do

  !$omp target data map(to: ions_part(0:particle_stride * num_ions - 1), &
  !$omp& electrons_part(0:particle_stride * num_electrons - 1)) &
  !$omp& map(alloc: ndi(0:ni - 1), nde(0:ni - 1), ef(0:ni - 1))

  call scatter_species(num_ions, ions_part, ndi, ions_spwt, sp_time)
  call scatter_species(num_electrons, electrons_part, nde, electrons_spwt, sp_time)
  call scatter_species_host(num_ions, ions_part, ndi_ref, ions_spwt)
  call scatter_species_host(num_electrons, electrons_part, nde_ref, electrons_spwt)
  call validate_density(ndi, ndi_ref, 'ion')
  call validate_density(nde, nde_ref, 'electron')

  call compute_rho(ions_charge, electrons_charge, ndi, nde, rho)
  call solve_potential(phi, rho)
  call compute_ef(phi, ef)
  call rewind_species(num_electrons, electrons_part, electrons_charge / electrons_mass, ef)
  call rewind_species(num_ions, ions_part, ions_charge / ions_mass, ef)

  open(unit=10, file='result.dat', status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open result.dat'
  write(10,'(A)') 'VARIABLES = x nde ndi rho phi ef'
  call write_results(10, 0, nde, ndi, rho, phi, ef)

  start_time = omp_get_wtime()
  do ts = 1, num_ts
    call scatter_species(num_ions, ions_part, ndi, ions_spwt, sp_time)
    call scatter_species(num_electrons, electrons_part, nde, electrons_spwt, sp_time)
    call compute_rho(ions_charge, electrons_charge, ndi, nde, rho)
    call solve_potential(phi, rho)
    call compute_ef(phi, ef)
    call push_species(num_electrons, electrons_part, electrons_charge / electrons_mass, ef)
    call push_species(num_ions, ions_part, ions_charge / ions_mass, ef)

    if (mod(ts, 25) == 0) then
      max_phi = maxval(abs(phi))
      write(*,'(A,I0,A,I0,A,I0,A,A)') 'TS:', ts, char(9)//'np_i:', num_ions, &
        char(9)//'np_e:', num_electrons, char(9)//'dphi:', trim(format_g3(max_phi - phi(0)))
    end if

    if (mod(ts, 1000) == 0) call write_results(10, ts, nde, ndi, rho, phi, ef)
  end do
  end_time = omp_get_wtime()
  total_time = end_time - start_time
  close(10)

  !$omp end target data

  call solve_potential_direct(phi_direct, rho)
  call validate_phi(phi, phi_direct)

  write(*,'(A,A,A)') 'Total kernel execution time (scatter particles) : ', trim(format_g3(sp_time)), ' (s)'
  write(*,'(A,I0,A,A,A)') 'Total time for ', num_ts, ' time steps: ', trim(format_g3(total_time)), ' (s)'
  write(*,'(A,A,A)') 'Time per time step: ', &
    trim(format_g3((total_time * 1.0e3_real64) / real(num_ts, real64))), ' (ms)'

  deallocate(electrons_part, ions_part)
  deallocate(nde_ref, ndi_ref, nde, ndi, phi_direct, ef, rho, phi)

contains

  real(real64) function rnd() result(value)
    value = real(c_rand(), real64) / 2147483647.0_real64
  end function rnd

  real(real64) function sample_vel(v_th) result(value)
    real(real64), intent(in) :: v_th
    integer, parameter :: m = 12
    integer :: i
    real(real64) :: sum_value

    sum_value = 0.0_real64
    do i = 1, m
      sum_value = sum_value + rnd()
    end do
    value = sqrt(0.5_real64) * v_th * (sum_value - real(m, real64) / 2.0_real64) / &
      sqrt(real(m, real64) / 12.0_real64)
  end function sample_vel

  character(len=32) function format_g3(value) result(text)
    real(real64), intent(in) :: value
    character(len=32) :: buffer
    integer :: dot_pos, exp_pos, last

    write(buffer,'(G0.3)') value
    buffer = adjustl(buffer)
    exp_pos = index(buffer, 'E')
    if (exp_pos == 0) exp_pos = index(buffer, 'D')
    if (exp_pos > 0) then
      last = exp_pos - 1
    else
      last = len_trim(buffer)
    end if

    dot_pos = index(buffer(:last), '.')
    if (dot_pos > 0) then
      do while (last > dot_pos .and. buffer(last:last) == '0')
        last = last - 1
      end do
      if (last == dot_pos) last = last - 1
    end if

    if (exp_pos > 0) then
      text = buffer(:last) // buffer(exp_pos:len_trim(buffer))
    else
      text = buffer(:last)
    end if
  end function format_g3

  subroutine scatter_species(npart, particles, den, spwt, time_accum)
    integer, intent(in) :: npart
    real(real64), intent(in) :: particles(0:)
    real(real32), intent(inout) :: den(0:)
    real(real64), intent(in) :: spwt
    real(real64), intent(inout) :: time_accum
    integer :: i
    real(real64) :: t0, t1

    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do i = 0, ni - 1
      den(i) = 0.0_real32
    end do

    t0 = omp_get_wtime()
    !$omp target teams distribute parallel do thread_limit(threads_per_block)
    do i = 0, npart - 1
      if (particles(particle_stride * i + particle_alive) /= 0.0_real64) &
        call scatter_value(particles(particle_stride * i + particle_x) / dx, 1.0_real32, den)
    end do
    t1 = omp_get_wtime()
    time_accum = time_accum + (t1 - t0)

    !$omp target update from(den(0:ni - 1))
    do i = 0, ni - 1
      den(i) = den(i) * real(spwt / dx, real32)
    end do
    den(0) = den(0) * 2.0_real32
    den(ni - 1) = den(ni - 1) * 2.0_real32
  end subroutine scatter_species

  subroutine scatter_species_host(npart, particles, den, spwt)
    integer, intent(in) :: npart
    real(real64), intent(in) :: particles(0:)
    real(real32), intent(out) :: den(0:)
    real(real64), intent(in) :: spwt
    integer :: i

    den = 0.0_real32
    do i = 0, npart - 1
      if (particles(particle_stride * i + particle_alive) /= 0.0_real64) &
        call scatter_value_host(particles(particle_stride * i + particle_x) / dx, 1.0_real32, den)
    end do
    do i = 0, ni - 1
      den(i) = den(i) * real(spwt / dx, real32)
    end do
    den(0) = den(0) * 2.0_real32
    den(ni - 1) = den(ni - 1) * 2.0_real32
  end subroutine scatter_species_host

  subroutine compute_rho(ion_charge, electron_charge, ion_den, electron_den, charge_density)
    real(real64), intent(in) :: ion_charge, electron_charge
    real(real32), intent(in) :: ion_den(0:), electron_den(0:)
    real(real64), intent(out) :: charge_density(0:)
    integer :: i

    do i = 0, ni - 1
      charge_density(i) = ion_charge * real(ion_den(i), real64) + &
        electron_charge * real(electron_den(i), real64)
    end do
  end subroutine compute_rho

  subroutine solve_potential(potential, charge_density)
    real(real64), intent(inout) :: potential(0:)
    real(real64), intent(in) :: charge_density(0:)
    integer :: solver_it, i
    real(real64) :: dx2, g, residual, sum_value, l2

    dx2 = dx * dx
    potential(0) = 0.0_real64
    potential(ni - 1) = 0.0_real64
    l2 = huge(1.0_real64)

    do solver_it = 0, 39999
      do i = 1, ni - 2
        g = 0.5_real64 * (potential(i - 1) + potential(i + 1) + dx2 * charge_density(i) / eps_0)
        potential(i) = potential(i) + 1.4_real64 * (g - potential(i))
      end do

      if (mod(solver_it, 25) == 0) then
        sum_value = 0.0_real64
        do i = 1, ni - 2
          residual = -charge_density(i) / eps_0 - &
            (potential(i - 1) - 2.0_real64 * potential(i) + potential(i + 1)) / dx2
          sum_value = sum_value + residual * residual
        end do
        l2 = sqrt(sum_value) / real(ni, real64)
        if (l2 < 1.0e-4_real64) return
      end if
    end do
    write(*,'(A,ES9.3,A)') 'Gauss-Seidel solver failed to converge, L2=', l2, '!'
  end subroutine solve_potential

  subroutine solve_potential_direct(solution, charge_density)
    real(real64), intent(out) :: solution(0:)
    real(real64), intent(in) :: charge_density(0:)
    real(real64) :: a(0:ni - 1), b(0:ni - 1), c(0:ni - 1)
    real(real64) :: dx2, id
    integer :: i

    dx2 = dx * dx
    do i = 1, ni - 2
      a(i) = 1.0_real64
      b(i) = -2.0_real64
      c(i) = 1.0_real64
      solution(i) = -charge_density(i) * dx2 / eps_0
    end do
    a(0) = 0.0_real64
    b(0) = 1.0_real64
    c(0) = 0.0_real64
    solution(0) = 0.0_real64
    a(ni - 1) = 0.0_real64
    b(ni - 1) = 1.0_real64
    c(ni - 1) = 0.0_real64
    solution(ni - 1) = 0.0_real64

    c(0) = c(0) / b(0)
    solution(0) = solution(0) / b(0)
    do i = 1, ni - 1
      id = b(i) - c(i - 1) * a(i)
      c(i) = c(i) / id
      solution(i) = (solution(i) - solution(i - 1) * a(i)) / id
    end do
    do i = ni - 2, 0, -1
      solution(i) = solution(i) - c(i) * solution(i + 1)
    end do
  end subroutine solve_potential_direct

  subroutine compute_ef(potential, field)
    real(real64), intent(in) :: potential(0:)
    real(real64), intent(inout) :: field(0:)
    integer :: i

    do i = 1, ni - 2
      field(i) = -(potential(i + 1) - potential(i - 1)) / (2.0_real64 * dx)
    end do
    field(0) = -(potential(1) - potential(0)) / dx
    field(ni - 1) = -(potential(ni - 1) - potential(ni - 2)) / dx
    !$omp target update to(field(0:ni - 1))
  end subroutine compute_ef

  subroutine push_species(npart, particles, qm, field)
    integer, intent(in) :: npart
    real(real64), intent(inout) :: particles(0:)
    real(real64), intent(in) :: qm
    real(real64), intent(in) :: field(0:)
    integer :: i
    real(real64) :: lc, part_ef

    !$omp target teams distribute parallel do thread_limit(threads_per_block) private(lc, part_ef)
    do i = 0, npart - 1
      if (particles(particle_stride * i + particle_alive) /= 0.0_real64) then
        lc = particles(particle_stride * i + particle_x) / dx
        part_ef = gather_value(lc, field)
        particles(particle_stride * i + particle_v) = &
          particles(particle_stride * i + particle_v) + dt * qm * part_ef
        particles(particle_stride * i + particle_x) = &
          particles(particle_stride * i + particle_x) + &
          dt * particles(particle_stride * i + particle_v)
        if (particles(particle_stride * i + particle_x) < x0 .or. &
            particles(particle_stride * i + particle_x) >= xmax) &
          particles(particle_stride * i + particle_alive) = 0.0_real64
      end if
    end do
  end subroutine push_species

  subroutine rewind_species(npart, particles, qm, field)
    integer, intent(in) :: npart
    real(real64), intent(inout) :: particles(0:)
    real(real64), intent(in) :: qm
    real(real64), intent(in) :: field(0:)
    integer :: i
    real(real64) :: lc, part_ef

    !$omp target teams distribute parallel do thread_limit(threads_per_block) private(lc, part_ef)
    do i = 0, npart - 1
      if (particles(particle_stride * i + particle_alive) /= 0.0_real64) then
        lc = particles(particle_stride * i + particle_x) / dx
        part_ef = gather_value(lc, field)
        particles(particle_stride * i + particle_v) = &
          particles(particle_stride * i + particle_v) - 0.5_real64 * dt * qm * part_ef
      end if
    end do
  end subroutine rewind_species

  subroutine write_results(unit_id, step, electron_den, ion_den, charge_density, potential, field)
    integer, intent(in) :: unit_id, step
    real(real32), intent(in) :: electron_den(0:), ion_den(0:)
    real(real64), intent(in) :: charge_density(0:), potential(0:), field(0:)
    integer :: i

    write(unit_id,'(A,I0,A,I6.6)') 'ZONE I=', ni, ' T=ZONE_', step
    do i = 0, ni - 1
      write(unit_id,'(6(ES14.6,1X))') real(i, real64) * dx, real(electron_den(i), real64), &
        real(ion_den(i), real64), charge_density(i), potential(i), field(i)
    end do
    flush(unit_id)
  end subroutine write_results

  subroutine validate_density(actual, expected, label)
    real(real32), intent(in) :: actual(0:), expected(0:)
    character(len=*), intent(in) :: label
    real(real64) :: diff, scale

    diff = maxval(abs(real(actual - expected, real64)))
    scale = max(1.0_real64, maxval(abs(real(expected, real64))))
    if (diff / scale > 5.0e-4_real64) then
      write(*,'(A,A,A,ES12.4)') 'Density validation failed for ', trim(label), ', relative error ', diff / scale
      error stop 2
    end if
  end subroutine validate_density

  subroutine validate_phi(actual, expected)
    real(real64), intent(in) :: actual(0:), expected(0:)
    real(real64) :: diff, scale

    diff = maxval(abs(actual - expected))
    scale = max(1.0_real64, maxval(abs(expected)))
    if (diff / scale > 1.0e-3_real64) then
      write(*,'(A,ES12.4)') 'Potential validation failed, relative error ', diff / scale
      error stop 3
    end if
  end subroutine validate_phi

  subroutine scatter_value(lc, value, field)
    !$omp declare target
    real(real64), intent(in) :: lc
    real(real32), intent(in) :: value
    real(real32), intent(inout) :: field(0:)
    integer :: i
    real(real32) :: di

    i = int(lc)
    di = real(lc - real(i, real64), real32)
    !$omp atomic update
    field(i) = field(i) + value * (1.0_real32 - di)
    !$omp atomic update
    field(i + 1) = field(i + 1) + value * di
  end subroutine scatter_value

  subroutine scatter_value_host(lc, value, field)
    real(real64), intent(in) :: lc
    real(real32), intent(in) :: value
    real(real32), intent(inout) :: field(0:)
    integer :: i
    real(real32) :: di

    i = int(lc)
    di = real(lc - real(i, real64), real32)
    field(i) = field(i) + value * (1.0_real32 - di)
    field(i + 1) = field(i + 1) + value * di
  end subroutine scatter_value_host

  real(real64) function gather_value(lc, field) result(value)
    !$omp declare target
    real(real64), intent(in) :: lc
    real(real64), intent(in) :: field(0:)
    integer :: i
    real(real64) :: di

    i = int(lc)
    di = lc - real(i, real64)
    value = field(i) * (1.0_real64 - di) + field(i + 1) * di
  end function gather_value

end program main
