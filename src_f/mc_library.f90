module mc_library
    private
    public :: &
        collision_distance_sample, &
        interface_distance_sample, &
        material_sample, &
        move_particle, &
        move_particle_lp

    contains

    ! Markovian sample of a distance
    function collision_distance_sample(sigma_tot, rand_num) result(dist)
        implicit none

        real(8), intent(in) :: &
            sigma_tot, rand_num
        real(8) :: &
            dist

        dist = -dlog(rand_num) / sigma_tot
    end function collision_distance_sample

    function interface_distance_sample(chord, mu, rand_num) result(dist)
        implicit none

        real(8), intent(in) :: &
            chord, mu, rand_num
        real(8) :: &
            dist

        dist = -dlog(rand_num) * chord / dabs(mu)
    end function interface_distance_sample

    function material_sample(first_prob, rand_num) result(mat_num)
        implicit none

        real(8), intent(in) :: &
            first_prob, rand_num
        integer :: &
            mat_num

        if (rand_num < first_prob) then
            mat_num = 1
        else
            mat_num = 2
        end if
    end function material_sample

    subroutine move_particle(delta_x, last_delta_x, dist_in_cell, mu, sigma_scat, sigma_tot, distance_to_collision, &
        distance, exists, scattered, leakage_l, leakage_r, absorbed, cell_index, num_cells, phi)
        use mcnp_random

        implicit none

        real(8), intent(in) :: &
            sigma_scat, sigma_tot, delta_x, last_delta_x, mu, distance_to_collision
        integer, intent(in) :: &
            num_cells
        integer, intent(inout) :: &
            cell_index
        real(8), intent(inout) :: &
            leakage_l, leakage_r, absorbed, phi, distance, dist_in_cell
        logical, intent(out) :: &
            exists, scattered
        real(8) :: &
            distance_to_boundary, flight_distance

        if (mu > 0.0d+0) then
            distance_to_boundary = (delta_x - dist_in_cell) / dabs(mu)
        else
            distance_to_boundary = dist_in_cell / dabs(mu)
        end if

        ! Take the minimum of compared event distances
        if (distance_to_boundary < distance_to_collision) then
            ! Particle escapes region
            scattered = .false.
            flight_distance = distance_to_boundary
            if (mu > 0.0d+0) then
                cell_index = cell_index + 1
                if (cell_index > num_cells) then
                    leakage_r = leakage_r + 1.0d+0
                    exists = .false.
                else
                    dist_in_cell = 0.0d+0  ! cm
                end if
            else
                cell_index = cell_index - 1
                if (cell_index < 1) then
                    leakage_l = leakage_l + 1.0d+0
                    exists = .false.
                else
                    ! This may need to be changed for the use of unstructured meshes
                    dist_in_cell = last_delta_x  ! cm
                end if
            end if
        else
            ! Collision occurs
            flight_distance = distance_to_collision
            if (rang() > sigma_scat / sigma_tot) then
                ! Absorption
                absorbed = absorbed + 1.0d+0
                exists = .false.
                scattered = .false.
            else
                ! Scatter
                scattered = .true.
                dist_in_cell = dist_in_cell + distance_to_collision * mu  ! cm
            end if
        end if

        ! Apply distance tally for flux tracking
        distance = distance + flight_distance * mu  ! cm
        phi = phi + flight_distance
    end subroutine move_particle

    subroutine move_particle_lp(delta_x, dist_in_cell, mu, sigma_scat, sigma_tot, &
        distance_to_collision, distance_to_interface, distance, exists, scattered, &
        leakage_l, leakage_r, absorbed, cell_index, num_cells, phi, mat_num)
        use mcnp_random
    
        implicit none
    
        real(8), intent(in) :: &
            delta_x, mu, sigma_scat, sigma_tot, distance_to_interface, distance_to_collision
        integer, intent(in) :: &
            num_cells
        real(8), intent(inout) :: &
            leakage_l, leakage_r, absorbed, distance, dist_in_cell, phi
        integer, intent(inout) :: &
            mat_num, cell_index
        logical, intent(out) :: &
            exists, scattered
        real(8) :: &
            distance_to_boundary, flight_distance
    
        ! Compute the distance to cell boundary
        if (mu > 0.0d+0) then
            distance_to_boundary = (delta_x - dist_in_cell) / dabs(mu)  ! cm
        else
            distance_to_boundary = dist_in_cell / dabs(mu)  ! cm
        end if
    
        ! Take the minimum of compared event distances
        if ((distance_to_boundary < distance_to_interface) .and. (distance_to_boundary < distance_to_collision)) then
            ! Particle escapes region
            scattered = .false.
            flight_distance = distance_to_boundary  ! cm
            if (mu > 0.0d+0) then
                cell_index = cell_index + 1
                if (cell_index > num_cells) then
                    leakage_r = leakage_r + 1.0d+0
                    exists = .false.
                else
                    dist_in_cell = 0.0d+0  ! cm
                end if
            else
                cell_index = cell_index - 1
                if (cell_index < 1) then
                    leakage_l = leakage_l + 1.0d+0
                    exists = .false.
                else
                    ! This may need to be changed for the use of unstructured meshes
                    dist_in_cell = delta_x  ! cm
                end if
            end if
        else if ((distance_to_collision < distance_to_interface) .and. (distance_to_collision < distance_to_boundary)) then
            ! Collision occurs
            flight_distance = distance_to_collision  ! cm
            if (rang() > sigma_scat / sigma_tot) then
                ! Absorption
                absorbed = absorbed + 1.0d+0
                exists = .false.
                scattered = .false.
            else
                ! Scatter
                scattered = .true.
                dist_in_cell = dist_in_cell + distance_to_collision * mu  ! cm
            end if
        else
            ! Minimum is distance_to_interface
            flight_distance = distance_to_interface  ! cm
            scattered = .false.
            ! Two material equations
            if (mat_num /= 1) then
                mat_num = 1
            else
                mat_num = 2
            end if
            dist_in_cell = dist_in_cell + distance_to_interface * mu  ! cm
        end if

        ! Apply distance tally for flux tracking
        distance = distance + flight_distance * mu  ! cm
        phi = phi + flight_distance
    end subroutine move_particle_lp
end module mc_library
