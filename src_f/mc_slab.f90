program mc_slab
    use mcnp_random
    use self_library
    use mc_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+7, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4)

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_cells), &
        !tot_const = 1.0d+0, &  ! 1/cm
        !first_lambda = dble(99) / dble(10), &
        !second_lambda = dble(11) / dble(10), &
        first_lambda = dble(101) / dble(20), &
        second_lambda = dble(101) / dble(20), &
        !first_xs = dble(10) / dble(99), &
        !second_xs = dble(100) / dble(11), &
        first_xs = dble(2) / dble(101), &
        second_xs = dble(200) / dble(101), &
        first_prob = first_lambda / (first_lambda + second_lambda), &
        second_prob = 1.0d+0 - first_prob, &
        tot_const = first_xs * first_prob + second_xs * second_prob, &  ! 1/cm
        scat_const = 0.9d+0 * first_xs * first_prob + 0.9d+0 * second_xs * second_prob  ! 1/cm

    ! Material variables
    real(8), dimension(num_cells) :: &
        delta_x, macro_scat, macro_tot

    ! Calculation variables
    integer :: &
        c, cell_index
    integer(8) :: &
        particle_counter
    real(8) :: &
        leakage_l, leakage_r, mu, mu_0, azimuth, psi_bound_l, psi_bound_r, distance, absorbed, &
        flight_distance, dist_in_cell
    real(8), dimension(num_cells) :: &
        phi
    logical :: &
        exists

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        point_vector

    ! Assignment of material variables
    delta_x(:) = thickness / dble(num_cells)  ! cm
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_scat(:) = scat_const  ! 1/cm
    macro_tot(:) = tot_const  ! 1/cm

    ! Boundary conditions: 1/cm^2-s-MeV
    psi_bound_l = 1.0d+0  ! Isotropic
    psi_bound_r = 0.0d+0  ! Vacuum

    ! Tallies
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    absorbed = 0.0d+0

    call RN_init_problem(int(123456, 8), int(1, 4))

    !$omp parallel do default(private) shared(delta_x,macro_scat,macro_tot) reduction(+:leakage_l,leakage_r,absorbed,phi)
    do particle_counter = int(1, 8), num_particles, int(1, 8)
        ! Spawn the particle
        exists = .true.
        mu = dsqrt(rang())
        flight_distance = distance_sample(tot_const, rang())  ! cm
        distance = flight_distance * mu  ! cm
        call check_distance(thickness, scat_const, tot_const, distance, exists, leakage_l, leakage_r, absorbed)
        ! Start of geometry
        cell_index = 1
        dist_in_cell = 0.0d+0  ! cm
        call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
        ! Continue distance sampling and checking after a scattering collision
        do while (exists)
            mu_0 = 2.0d+0 * rang() - 1.0d+0
            azimuth = 2.0d+0 * PI * rang()
            mu = mu * mu_0 + dsqrt(1.0d+0 - mu * mu) * dsqrt(1.0d+0 - mu_0 * mu_0) * dcos(azimuth)
            flight_distance = distance_sample(tot_const, rang())  ! cm
            distance = distance + flight_distance * mu  ! cm
            call check_distance(thickness, scat_const, tot_const, distance, exists, leakage_l, leakage_r, absorbed)
            call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
        end do
    end do
    !$omp end parallel do

    print *, "Leakage Left: ", leakage_l / num_particles
    print *, "Leakage Right: ", leakage_r / num_particles
    print *, "Absorbed: ", absorbed / num_particles
    print *, "Total: ", (leakage_l + leakage_r + absorbed) / num_particles

    ! Average the values for flux
    do c = 1, num_cells
        phi(c) = phi(c) / (delta_x(c) * num_particles)
    end do

    ! Create plot
    call linspace(point_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/mc_slab.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(7,*) point_vector(c), phi(c)
    end do
    close(7)
end program mc_slab