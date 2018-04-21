program mc_slab
    use mcnp_random
    use self_library
    use mc_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+7, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4), &
        num_materials = 2

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
        first_scat_rat = 0.9d+0, &
        second_scat_rat = 0.9d+0!, &
        !tot_const = first_xs * first_prob + second_xs * second_prob, &  ! 1/cm
        !scat_const = first_scat_rat * first_xs * first_prob + second_scat_rat * second_xs * second_prob  ! 1/cm

    ! Material variables
    real(8), dimension(num_cells) :: &
        delta_x
    real(8), dimension(num_materials) :: &
        macro_scat, macro_tot, prob, &
        chords = (/first_lambda, second_lambda/)

    ! Calculation variables
    integer :: &
        c, cell_index, mat_num, flight_mat_num
    integer(8) :: &
        particle_counter
    real(8) :: &
        leakage_l, leakage_r, leakage_l_1, leakage_l_2, leakage_r_1, leakage_r_2, mu, mu_0, azimuth, distance, absorbed, &
        collision_distance, interface_distance, flight_distance, dist_in_cell!, psi_bound_l, psi_bound_r
    real(8), dimension(num_cells) :: &
        phi, phi_1, phi_2
    logical :: &
        exists, mat_change

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        point_vector

    ! Assignment of material variables
    delta_x(:) = thickness / dble(num_cells)  ! cm
    ! Initialization of fluxes
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    phi_1(:) = 0.0d+0  ! 1/cm^2-s-MeV
    phi_2(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_tot(1) = first_xs  ! 1/cm
    macro_tot(2) = second_xs  ! 1/cm
    macro_scat(1) = macro_tot(1) * first_scat_rat  ! 1/cm
    macro_scat(2) = macro_tot(2) * second_scat_rat  ! 1/cm
    prob(1) = first_prob
    prob(2) = second_prob

    ! Boundary conditions: 1/cm^2-s-MeV
    !psi_bound_l = 1.0d+0  ! Isotropic
    !psi_bound_r = 0.0d+0  ! Vacuum

    ! Tallies
    leakage_l_1 = 0.0d+0
    leakage_r_1 = 0.0d+0
    leakage_l_2 = 0.0d+0
    leakage_r_2 = 0.0d+0
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    absorbed = 0.0d+0

    call RN_init_problem(int(123456, 8), int(1, 4))

    !$omp parallel do default(private) shared(delta_x,macro_scat,macro_tot) reduction(+:leakage_l,leakage_r,absorbed,phi)
    do particle_counter = int(1, 8), num_particles, int(1, 8)
        ! Sample the material number
        mat_num = material_sample(first_prob, rang())
        flight_mat_num = mat_num
        ! Spawn the particle
        exists = .true.
        mu = dsqrt(rang())
        collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
        interface_distance = interface_distance_sample(chords(mat_num), mu, rang())  ! cm
        flight_distance = 0.0d+0  ! cm
        if (mat_num == 1) then
            call check_distance_lp(thickness, macro_scat(mat_num), macro_tot(mat_num), mu, &
            distance, interface_distance, collision_distance, exists, leakage_l_1, leakage_r_1, &
            absorbed, prob(mat_num), flight_distance, mat_num, mat_change)
        else
            call check_distance_lp(thickness, macro_scat(mat_num), macro_tot(mat_num), mu, &
            distance, interface_distance, collision_distance, exists, leakage_l_2, leakage_r_2, &
            absorbed, prob(mat_num), flight_distance, mat_num, mat_change)
        end if
        distance = flight_distance * mu  ! cm
        ! Start of geometry
        cell_index = 1
        dist_in_cell = 0.0d+0  ! cm
        if (flight_mat_num == 1) then
            call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
        else
            call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
        end if
        ! Continue distance sampling and checking after a scattering collision or material change
        do while (exists)
            ! Do not recompute with a material transfer
            if (.not. mat_change) then
                mu_0 = 2.0d+0 * rang() - 1.0d+0
                azimuth = 2.0d+0 * PI * rang()
                mu = mu * mu_0 + dsqrt(1.0d+0 - mu * mu) * dsqrt(1.0d+0 - mu_0 * mu_0) * dcos(azimuth)
            end if
            collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
            interface_distance = interface_distance_sample(chords(mat_num), mu, rang())  ! cm
            flight_mat_num = mat_num
            if (mat_num == 1) then
                call check_distance_lp(thickness, macro_scat(mat_num), macro_tot(mat_num), mu, &
                distance, interface_distance, collision_distance, exists, leakage_l_1, leakage_r_1, &
                absorbed, prob(mat_num), flight_distance, mat_num, mat_change)
            else
                call check_distance_lp(thickness, macro_scat(mat_num), macro_tot(mat_num), mu, &
                distance, interface_distance, collision_distance, exists, leakage_l_2, leakage_r_2, &
                absorbed, prob(mat_num), flight_distance, mat_num, mat_change)
            end if
            distance = distance + flight_distance * mu  ! cm
            if (flight_mat_num == 1) then
                call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
            else
                call tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
            end if
            ! TODO: Bring this print statement to non-existance
            if (((cell_index - 1) * delta_x(1) > distance) .or. (cell_index * delta_x(1) < distance)) then
                print *, (cell_index - 1) * delta_x(1), distance, cell_index * delta_x(1)
            end if
        end do
    end do
    !$omp end parallel do

    leakage_l = leakage_l_1 + leakage_l_2
    leakage_r = leakage_r_1 + leakage_r_2
    print *, "Leakage Left: ", leakage_l / num_particles
    print *, "Leakage Right: ", leakage_r / num_particles
    print *, "Absorbed: ", absorbed / num_particles
    print *, "Total: ", (leakage_l + leakage_r + absorbed) / num_particles

    ! Average the values for flux
    do c = 1, num_cells
        !phi_1(c) = phi_1(c) / (delta_x(c) * num_particles * first_prob)
        !phi_2(c) = phi_2(c) / (delta_x(c) * num_particles * second_prob)
        !phi(c) = first_prob * phi_1(c) + second_prob * phi_2(c)
        phi(c) = phi(c) / (delta_x(c) * num_particles)
    end do

    ! Create plot
    call linspace(point_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/mc_closure_slab.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(7,*) point_vector(c), phi(c)
    end do
    close(7)
    open(unit=8, file="./out/mc_closure_slab_1.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(8,*) point_vector(c), phi_1(c)
    end do
    close(8)
    open(unit=9, file="./out/mc_closure_slab_2.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(9,*) point_vector(c), phi_2(c)
    end do
    close(9)
end program mc_slab
