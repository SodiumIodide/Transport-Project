program mc_slab
    use mcnp_random
    use self_library
    use mc_library
    use geometry_gen
    use mesh_map

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+4, 8), &
        num_realizations = int(1.0d+4, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4), &
        num_materials = 2, &
        num_geom_divs = 10

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
        second_scat_rat = 0.9d+0

    ! Material variables
    real(8), dimension(num_materials) :: &
        macro_scat_mat, macro_tot_mat, prob, &
        chords = (/first_lambda, second_lambda/)

    ! Allocated as (num_ind_cells) (by subroutine)
    real(8), dimension(:), allocatable :: &
        macro_scat, macro_tot, delta_x
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials

    ! Calculation variables
    integer :: &
        c, k, cell_index, num_ind_cells
    integer(8) :: &
        particle_counter, realization_counter
    real(8) :: &
        leakage_l, leakage_r, mu, mu_0, azimuth, distance, absorbed, &
        collision_distance, interface_distance, dist_in_cell, last_delta_x
    real(8), dimension(num_cells, num_materials) :: &
        phi_mat
    real(8), dimension(num_cells) :: &
        phi
    logical :: &
        exists, scattered
    real(8), dimension(:), allocatable :: &
        phi_morph

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        point_vector

    ! Assignment of material variables
    ! Initialization of fluxes
    phi_mat(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    !phi_1(:) = 0.0d+0  ! 1/cm^2-s-MeV
    !phi_2(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_tot_mat(1) = first_xs  ! 1/cm
    macro_tot_mat(2) = second_xs  ! 1/cm
    macro_scat_mat(1) = macro_tot_mat(1) * first_scat_rat  ! 1/cm
    macro_scat_mat(2) = macro_tot_mat(2) * second_scat_rat  ! 1/cm
    prob(1) = first_prob
    prob(2) = second_prob

    ! Tallies
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    absorbed = 0.0d+0

    call RN_init_problem(int(123456, 8), int(1, 4))

    do realization_counter = int(1, 8), num_realizations, int(1, 8)
        ! Fill Markovian geometry
        ! Allocations and deallocations of delta_x and materials should be
        ! Handled by the subroutine
        call get_geometry(delta_x, x_points, materials, chords(1), chords(2), &
            thickness, num_ind_cells, num_geom_divs)

        ! Allocate and define properties
        ! Anything that depends on the number of individual cells must be reallocated

        ! Allocate and define calculational arrays
        if (allocated(macro_tot)) then
            deallocate(macro_tot)
        end if
        if (allocated(macro_scat)) then
            deallocate(macro_scat)
        end if
        allocate(macro_tot(num_ind_cells))
        allocate(macro_scat(num_ind_cells))
        do k = 1, num_ind_cells
            if (materials(k) == 1) then
                macro_tot(k) = macro_tot_mat(1)
                macro_scat(k) = macro_scat_mat(1)
            else
                macro_tot(k) = macro_tot_mat(2)
                macro_scat(k) = macro_scat_mat(2)
            end if
        end do
        if (allocated(phi_morph)) then
            deallocate(phi_morph)
        end if
        allocate(phi_morph(num_ind_cells))
        phi_morph(:) = 0.0d+0  ! 1/cm^2-s-MeV

        !$omp parallel do default(private) shared(delta_x,macro_scat,macro_tot,chords) reduction(+:leakage_l,leakage_r,absorbed,phi_mat)
        do particle_counter = int(1, 8), num_particles, int(1, 8)
            ! Sample the material number
            mat_num = material_sample(first_prob, rang())
            flight_mat_num = mat_num
            ! Spawn the particle
            exists = .true.
            scattered = .false.
            mu = dsqrt(rang())
            collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
            ! Start of geometry
            cell_index = 1
            dist_in_cell = 0.0d+0  ! cm
            distance = 0.0d+0  ! cm
            last_delta_x = 0.0d+0  ! cm
            call move_particle(delta_x(cell_index), last_delta_x, dist_in_cell, mu, macro_scat(mat_num), macro_tot(mat_num), &
                collision_distance, distance, exists, scattered, leakage_l, leakage_r, &
                absorbed, cell_index, num_ind_cells, phi_morph(cell_index))
            ! Continue distance sampling and checking after a scattering collision or material change
            do while (exists)
                ! Do not recompute with a material transfer
                if (scattered) then
                    mu_0 = 2.0d+0 * rang() - 1.0d+0
                    azimuth = 2.0d+0 * PI * rang()
                    mu = mu * mu_0 + dsqrt(1.0d+0 - mu * mu) * dsqrt(1.0d+0 - mu_0 * mu_0) * dcos(azimuth)
                end if
                collision_distance = collision_distance_sample(macro_tot(mat_num), rang())  ! cm
                if (cell_index > 1) then
                    last_delta_x = delta_x(cell_index - 1)  ! cm
                else
                    last_delta_x = 0.0d+0  ! cm
                end if
                call move_particle(delta_x(cell_index), last_delta_x, dist_in_cell, mu, macro_scat(mat_num), macro_tot(mat_num), &
                    collision_distance, distance, exists, scattered, leakage_l, leakage_r, &
                    absorbed, cell_index, num_ind_cells, phi_morph(cell_index))
            end do
        end do
        !$omp end parallel do

        call material_calc(phi_morph, delta_x, num_ind_cells, materials, phi_mat, &
            struct_thickness, num_cells, num_materials)
    end do

    !leakage_l = leakage_l_1 + leakage_l_2
    !leakage_r = leakage_r_1 + leakage_r_2
    print *, "Leakage Left: ", leakage_l / num_particles
    print *, "Leakage Right: ", leakage_r / num_particles
    print *, "Absorbed: ", absorbed / num_particles
    print *, "Total: ", (leakage_l + leakage_r + absorbed) / num_particles

    ! Average the values for flux
    do c = 1, num_cells
        !phi_1(c) = phi_1(c) / (delta_x(c) * num_particles * first_prob)
        !phi_2(c) = phi_2(c) / (delta_x(c) * num_particles * second_prob)
        !phi(c) = first_prob * phi_1(c) + second_prob * phi_2(c)
        phi_mat(c, 1) = phi_mat(c, 1) / (delta_x(c) * num_particles * first_prob)
        phi_mat(c, 2) = phi_mat(c, 2) / (delta_x(c) * num_particles * second_prob)
        phi(c) = first_prob * phi_mat(c, 1) + second_prob * phi_mat(c, 2)
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
        write(8,*) point_vector(c), phi_mat(c, 1)
    end do
    close(8)
    open(unit=9, file="./out/mc_closure_slab_2.out", form="formatted", status="replace", action="write")
    do c = 1, num_cells
        write(9,*) point_vector(c), phi_mat(c, 2)
    end do
    close(9)

    if (allocated(macro_tot)) then
        deallocate(macro_tot)
    end if
    if (allocated(macro_scat)) then
        deallocate(macro_scat)
    end if
    if (allocated(materials)) then
        deallocate(materials)
    end if
    if (allocated(delta_x)) then
        deallocate(delta_x)
    end if
end program mc_slab
