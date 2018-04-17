module mc_library
    private
    public :: &
        distance_sample, &
        check_distance, &
        tally_cells

    contains

    ! Markovian sample of a distance
    function distance_sample(sigma_tot, rand_num) result(dist)
        implicit none

        real(8), intent(in) :: &
            sigma_tot, rand_num
        real(8) :: &
            dist

        dist = -dlog(rand_num) / sigma_tot
    end function distance_sample

    subroutine check_distance(thickness, sigma_scat, sigma_tot, distance, exists, leakage_l, leakage_r, absorbed)
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
    end subroutine check_distance

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
