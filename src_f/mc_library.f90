module mc_library
    private
    public :: &
        collision_distance_sample, &
        interface_distance_sample, &
        material_sample, &
        move_particle, &
        check_distance_lp, &
        tally_cells

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

    subroutine move_particle(delta_x, dist_in_cell, mu, sigma_scat, sigma_tot, &
        distance_to_collision, distance, exists, scattered, leakage_l, leakage_r, absorbed, cell_index, num_cells, phi)
        use mcnp_random

        implicit none

        real(8), intent(in) :: &
            sigma_scat, sigma_tot, delta_x, mu, distance_to_collision
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
                    dist_in_cell = delta_x  ! cm
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
                dist_in_cell = dist_in_cell + distance_to_collision * mu
            end if
        end if

        ! Apply distance tally for flux tracking
        distance = distance + flight_distance * mu
        phi = phi + flight_distance
    end subroutine move_particle

    subroutine check_distance_x(thickness, sigma_scat, sigma_tot, distance, exists, leakage_l, leakage_r, absorbed)
        use mcnp_random

        implicit none

        real(8), intent(in) :: &
            thickness, sigma_scat, sigma_tot, distance
        real(8), intent(inout) :: &
            leakage_l, leakage_r, absorbed
        logical, intent(out) :: &
            exists

        if ((distance > 0.0d+0) .and. (distance < thickness)) then
            if (rang() > sigma_scat / sigma_tot) then
                exists = .false.
                absorbed = absorbed + 1.0d+0
            end if
        else if (distance < 0.0d+0) then
            exists = .false.
            leakage_l = leakage_l + 1.0d+0
        else if (distance > thickness) then
            exists = .false.
            leakage_r = leakage_r + 1.0d+0
        end if
    end subroutine check_distance_x

    subroutine check_distance_lp(thickness, sigma_scat, sigma_tot, mu, distance, distance_to_interface, &
        distance_to_collision, exists, leakage_l, leakage_r, absorbed, prob, flight_distance, mat_num, mat_change)
        use mcnp_random

        implicit none

        real(8), intent(in) :: &
            thickness, mu, sigma_scat, sigma_tot, distance_to_interface, distance_to_collision, prob
        real(8), intent(inout) :: &
            leakage_l, leakage_r, absorbed, distance
        integer, intent(inout) :: &
            mat_num
        real(8), intent(out) :: &
            flight_distance
        logical, intent(out) :: &
            exists, mat_change
        real(8) :: &
            distance_to_boundary

        ! "distance" should be the current value of the particle at its location prior to flight
        if (mu < 0.0d+0) then
            distance_to_boundary = distance / dabs(mu)  ! cm
        else
            distance_to_boundary = (thickness - distance) / dabs(mu)  ! cm
        end if

        ! Determine what the minimum value is
        if ((distance_to_boundary < distance_to_interface) .and. (distance_to_boundary < distance_to_collision)) then
            ! Minimum is distance_to_boundary
            exists = .false.
            mat_change = .false.
            flight_distance = distance_to_boundary  ! cm
            if (mu < 0.0d+0) then
                leakage_l = leakage_l + 1.0d+0
            else
                leakage_r = leakage_r + 1.0d+0
            end if
        else if ((distance_to_collision < distance_to_interface) .and. (distance_to_collision < distance_to_boundary)) then
            ! Minimum is distance_to_collision
            flight_distance = distance_to_collision  ! cm
            mat_change = .false.
            if (rang() > sigma_scat / sigma_tot) then
                exists = .false.
                absorbed = absorbed + 1.0d+0
            end if
        else
            ! Minimum is distance_to_interface
            flight_distance = distance_to_interface  ! cm
            mat_change = .true.
            ! Two material equations
            if (mat_num /= 1) then
                mat_num = 1
            else
                mat_num = 2
            end if
        end if
    end subroutine check_distance_lp

    subroutine tally_cells(cell_index, dist_in_cell, flight_distance, mu, delta_x, num_cells, phi)
        implicit none

        integer, intent(in) :: &
            num_cells
        real(8), intent(in) :: &
            mu
        real(8), dimension(num_cells), intent(in) :: &
            delta_x
        real(8), dimension(num_cells), intent(inout) :: &
            phi
        real(8), intent(inout) :: &
            flight_distance, dist_in_cell
        integer, intent(inout) :: &
            cell_index
        real(8) :: &
            in_cell, dist_to_boundary
        logical :: &
            boundary_cross, cont_tally

        ! First distance to boundary
        if (mu < 0.0d+0) then
            dist_to_boundary = dist_in_cell  ! cm
        else
            dist_to_boundary = delta_x(cell_index) - dist_in_cell  ! cm
        end if
        in_cell = 0.0d+0  ! cm
        boundary_cross = .false.
        cont_tally = .true.
        ! Loop for tracing the particle along its entire flight distance
        do while (cont_tally)
            if (dabs(flight_distance * mu) < dist_to_boundary) then
                in_cell = dabs(flight_distance)
                cont_tally = .false.
            else
                in_cell = dabs(dist_to_boundary / mu)
                boundary_cross = .true.
            end if
            ! Tally the "flux"
            phi(cell_index) = phi(cell_index) + in_cell
            ! Reduce the flight distance by the appropriate amount
            flight_distance = flight_distance - in_cell  ! cm
            ! Adjust the cell index if the particle crosses a boundary
            if (boundary_cross) then
                boundary_cross = .false.
                if (mu < 0.0d+0) then
                    cell_index = cell_index - 1
                    if (cell_index < 1) then
                        ! Particle escaped, stop tallies
                        flight_distance = 0.0d+0  ! cm
                        cont_tally = .false.
                    end if
                else
                    cell_index = cell_index + 1
                    if (cell_index > num_cells) then
                        ! Particle escaped, stop tallies
                        flight_distance = 0.0d+0  ! cm
                        cont_tally = .false.
                    end if
                end if
                ! Set the distance to boundary
                dist_to_boundary = delta_x(cell_index)
            end if
        end do
        ! Set the particle's location in the current cell (measured from left face)
        if (mu < 0.0d+0) then
            dist_in_cell = delta_x(cell_index) + in_cell * mu
        else
            dist_in_cell = in_cell * dabs(mu)
        end if
    end subroutine tally_cells
end module mc_library
