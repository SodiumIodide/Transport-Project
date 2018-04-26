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
    logical, parameter :: &
        use_alpha = .true.

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_cells), &
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
        second_scat_rat = 0.9d+0

    ! Material variables
    real(8), dimension(num_cells) :: &
        delta_x
    real(8), dimension(num_materials) :: &
        macro_scat, macro_tot, prob, &
        chords = (/first_lambda, second_lambda/)

    ! Calculation variables
    integer :: &
        c, cell_index, mat_num
    integer(8) :: &
        particle_counter
    real(8) :: &
        leakage_l, leakage_r, mu, mu_0, azimuth, distance, absorbed, &
        collision_distance, interface_distance, dist_in_cell, &
        alpha, sig_a_av, sig_t_av, first_abs, second_abs, &
        leakage_l_var, leakage_r_var, absorbed_var
    real(8), dimension(num_cells, num_materials) :: &
        phi_mat
    real(8), dimension(num_cells) :: &
        phi
    logical :: &
        exists, scattered

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        point_vector
    character(:), allocatable :: &
        filename

    ! Assignment of material variables
    delta_x(:) = struct_thickness  ! cm
    ! Initialization of fluxes
    phi_mat(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    !phi_1(:) = 0.0d+0  ! 1/cm^2-s-MeV
    !phi_2(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_tot(1) = first_xs  ! 1/cm
    macro_tot(2) = second_xs  ! 1/cm
    macro_scat(1) = macro_tot(1) * first_scat_rat  ! 1/cm
    macro_scat(2) = macro_tot(2) * second_scat_rat  ! 1/cm
    prob(1) = first_prob
    prob(2) = second_prob

    if (use_alpha) then
        first_abs = (1.0d+0 - first_scat_rat) * first_xs
        second_abs = (1.0d+0 - second_scat_rat) * second_xs
        sig_a_av = first_prob * first_abs + second_prob * second_abs
        sig_t_av = first_prob * first_xs + second_prob * second_xs
        alpha = dsqrt(sig_a_av * sig_t_av) &
            * (sig_a_av * (second_xs - first_xs)**2 + sig_t_av * (second_abs - first_abs)**2) &
            / ((sig_a_av * (second_xs - first_xs))**2 + (sig_t_av * (second_abs - first_abs))**2)
    else
        alpha = 1.0d+0
    end if
    print *, "Alpha: ", alpha

    ! Tallies
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    absorbed = 0.0d+0
    ! Variances
    leakage_l_var = 0.0d+0
    leakage_r_var = 0.0d+0
    absorbed_var = 0.0d+0

    call RN_init_problem(int(123456, 8), int(1, 4))

    !$omp parallel do default(private) shared(delta_x,macro_scat,macro_tot,chords,alpha) &
    !$omp& reduction(+:leakage_l,leakage_r,absorbed,phi_mat)
    do particle_counter = int(1, 8), num_particles, int(1, 8)
        ! Sample the material number
        mat_num = material_sample(first_prob, rang())
        ! Spawn the particle
        exists = .true.
        scattered = .false.
        mu = dsqrt(rang())
        collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
        interface_distance = interface_distance_sample(chords(mat_num), mu, rang(), alpha)  ! cm
        ! Start of geometry
        cell_index = 1
        dist_in_cell = 0.0d+0  ! cm
        distance = 0.0d+0  ! cm
        call move_particle_lp(delta_x(cell_index), dist_in_cell, mu, macro_scat(mat_num), macro_tot(mat_num), &
            collision_distance, interface_distance, distance, exists, scattered, leakage_l, leakage_r, &
            absorbed, cell_index, num_cells, phi_mat(cell_index, mat_num), mat_num)
        ! Continue distance sampling and checking after a scattering collision or material change
        do while (exists)
            ! Do not recompute with a material transfer
            if (scattered) then
                mu_0 = 2.0d+0 * rang() - 1.0d+0
                azimuth = 2.0d+0 * PI * rang()
                mu = mu * mu_0 + dsqrt(1.0d+0 - mu * mu) * dsqrt(1.0d+0 - mu_0 * mu_0) * dcos(azimuth)
            end if
            collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
            interface_distance = interface_distance_sample(chords(mat_num), mu, rang(), alpha)  ! cm
            call move_particle_lp(delta_x(cell_index), dist_in_cell, mu, macro_scat(mat_num), macro_tot(mat_num), &
                collision_distance, interface_distance, distance, exists, scattered, leakage_l, leakage_r, &
                absorbed, cell_index, num_cells, phi_mat(cell_index, mat_num), mat_num)
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
        phi_mat(c, 1) = phi_mat(c, 1) / (delta_x(c) * dble(num_particles) * first_prob)
        phi_mat(c, 2) = phi_mat(c, 2) / (delta_x(c) * dble(num_particles) * second_prob)
        phi(c) = first_prob * phi_mat(c, 1) + second_prob * phi_mat(c, 2)
    end do

    ! Create plot
    call linspace(point_vector, 0.0d+0, thickness, num_cells)
    if (use_alpha) then
        filename = "./out/mc_slab_closure_alpha.out"
    else
        filename = "./out/mc_slab_closure.out"
    end if
    open(unit=7, file=filename, form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(7,*) point_vector(c), phi(c)
    end do
    close(7)
    if (use_alpha) then
        filename = "./out/mc_slab_closure_alpha_1.out"
    else
        filename = "./out/mc_slab_closure_1.out"
    end if
    open(unit=8, file=filename, form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(8,*) point_vector(c), phi_mat(c, 1)
    end do
    close(8)
    if (use_alpha) then
        filename = "./out/mc_slab_closure_alpha_2.out"
    else
        filename = "./out/mc_slab_closure_2.out"
    end if
    open(unit=9, file=filename, form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(9,*) point_vector(c), phi_mat(c, 2)
    end do
    close(9)

    if (allocated(filename)) then
        deallocate(filename)
    end if
end program mc_slab
