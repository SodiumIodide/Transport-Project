program mc_slab
    use mcnp_random
    use self_library
    use mc_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+9, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4)

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_cells), &
        !first_lambda = dble(99) / dble(100), &
        !second_lambda = dble(11) / dble(100), &
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
        delta_x

    ! Calculation variables
    integer :: &
        c, cell_index
    integer(8) :: &
        particle_counter
    real(8) :: &
        leakage_l, leakage_r, mu, mu_0, azimuth, distance, absorbed, &
        collision_distance, dist_in_cell, macro_scat, macro_tot, &
        leakage_l_var, leakage_r_var, absorbed_var
    real(8), dimension(num_cells) :: &
        phi
    logical :: &
        exists, scattered

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        point_vector

    ! Assignment of material variables
    delta_x(:) = struct_thickness  ! cm
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_scat = scat_const  ! 1/cm
    macro_tot = tot_const  ! 1/cm

    ! Tallies
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    absorbed = 0.0d+0
    ! Variances
    leakage_l_var = 0.0d+0
    leakage_r_var = 0.0d+0
    absorbed_var = 0.0d+0

    call RN_init_problem(int(123456, 8), int(1, 4))

    !$omp parallel do default(private) shared(delta_x,macro_scat,macro_tot) reduction(+:leakage_l,leakage_r,absorbed,phi)
    do particle_counter = int(1, 8), num_particles, int(1, 8)
        ! Spawn the particle
        exists = .true.
        scattered = .false.
        mu = dsqrt(rang())
        collision_distance = collision_distance_sample(macro_tot, rang())  ! cm
        ! Start of geometry
        cell_index = 1
        dist_in_cell = 0.0d+0  ! cm
        distance = 0.0d+0  ! cm
        call move_particle(struct_thickness, struct_thickness, dist_in_cell, mu, macro_scat, macro_tot, &
            collision_distance, distance, exists, scattered, leakage_l, leakage_r, &
            absorbed, cell_index, num_cells, phi(cell_index))
        ! Continue distance sampling and checking as long as particle exists (scatter or in-geometry transfer)
        do while (exists)
            if (scattered) then
                mu_0 = 2.0d+0 * rang() - 1.0d+0
                azimuth = 2.0d+0 * PI * rang()
                mu = mu * mu_0 + dsqrt(1.0d+0 - mu * mu) * dsqrt(1.0d+0 - mu_0 * mu_0) * dcos(azimuth)
            end if
            collision_distance = collision_distance_sample(macro_tot, rang())  ! cm
            call move_particle(struct_thickness, struct_thickness, dist_in_cell, mu, macro_scat, macro_tot, &
                collision_distance, distance, exists, scattered, leakage_l, leakage_r, &
                absorbed, cell_index, num_cells, phi(cell_index))
        end do
    end do
    !$omp end parallel do

    leakage_l_var = dsqrt((leakage_l / dble(num_particles) - (leakage_l / dble(num_particles))**2) &
        * dble(num_particles) / dble(num_particles - 1)) / dsqrt(dble(num_particles))
    leakage_r_var = dsqrt((leakage_r / dble(num_particles) - (leakage_r / dble(num_particles))**2) &
        * dble(num_particles) / dble(num_particles - 1)) / dsqrt(dble(num_particles))
    absorbed_var = dsqrt((absorbed / dble(num_particles) - (absorbed / dble(num_particles))**2) &
        * dble(num_particles) / dble(num_particles - 1)) / dsqrt(dble(num_particles))
    print *, "Leakage Left: ", leakage_l / dble(num_particles), " +/- ", leakage_l_var
    print *, "Leakage Right: ", leakage_r / dble(num_particles), " +/- ", leakage_r_var
    print *, "Absorbed: ", absorbed / dble(num_particles), " +/- ", absorbed_var
    print *, "Balance (should be 1.0): ", (leakage_l + leakage_r + absorbed) / dble(num_particles)

    ! Average the values for flux
    do c = 1, num_cells
        phi(c) = phi(c) / (delta_x(c) * dble(num_particles))
    end do

    ! Create plot
    call linspace(point_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/mc_slab_non_mark.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(7,*) point_vector(c), phi(c)
    end do
    close(7)
end program mc_slab
