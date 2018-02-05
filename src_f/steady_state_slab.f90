program steady_state_slab
    use self_library
    use geometry_gen
    use mesh_map

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter = 5000
    integer, parameter :: &
        num_cells = 100, &
        num_groups = 1, &
        num_ords = 16, &
        num_materials = 2

    ! Material properties
    double precision, parameter :: &
        scat_const = 0.2d+0, &  ! 1/cm
        fis_const = 0.0d+0, &  ! 1/cm
        tot_const = 1.0d+0, &  ! 1/cm
        thickness = 1.0d+0, &  ! cm
        chord_a = 0.05d+0, &  ! cm
        chord_b = 0.05d+0, &  ! cm
        alpha_l = 0.0d+0, &
        alpha_r = 0.0d+0, &
        struct_thickness = thickness / dble(num_cells)  ! cm
    ! For alpha (albedo boundary), 0.0 = no refl., 1.0 = total refl.

    ! Material variables
    ! Allocated as (num_ind_cells) (by subroutine)
    double precision, dimension(:), allocatable :: &
        delta_x, x_points
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials
    ! Allocated as (num_groups, num_groups, num_cells)
    double precision, dimension(:, :, :), allocatable :: &
        macro_scat
    ! Allocated as (num_groups, num_cells)
    double precision, dimension(:, :), allocatable :: &
        macro_fis, macro_tot, phi_morph
    double precision, dimension(num_groups) :: &
        chi, nu

    ! Calculation variables
    integer(8) :: &
        iterations
    integer :: &
        i, c, g, g_prime, m, num_ind_cells
    logical :: &
        cont_calc, point_found
    double precision :: &
        tolerance, scatter_into, fission_into, weighted_sum, err, cons_distance
    double precision, dimension(num_groups, num_cells) :: &
        phi_new, phi_old, scat_source, fis_source, spont_source, tot_source
    double precision, dimension(num_cells, num_materials) :: &
        phi_mat
    ! Allocated as (num_groups, num_cells, num_ords)
    double precision, dimension(:, :, :), allocatable :: &
        psi, psi_refl, psi_i_p, psi_i_m
    double precision, dimension(num_groups, num_ords) :: &
        psi_bound_l, psi_bound_r
    double precision, dimension(num_ords) :: &
        ordinates, weights, mu

    ! Additional variables (plotting, etc.)
    double precision, dimension(num_cells) :: &
        cell_vector

    ! Assigment of material variables
    chi(:) = 1.0d+0 / dble(num_groups)  ! fission energy spectrum
    nu(:) = 2.4144d+0  ! average number of neutrons emitted per fission
    ! Assignment of initial calculation variables
    phi_new(:, :) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    ! Can adjust boundaries by group and by ordinate
    psi_bound_l(:, :) = 1.0d+0  ! 1/cm^2-s-MeV-strad
    psi_bound_r(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    scat_source(:, :) = 0.0d+0  ! 1/cm^3-s
    fis_source(:, :) = 0.0d+0  ! 1/cm^3-s
    spont_source(:, :) = 0.0d+0  ! 1/cm^3-s
    tot_source(:, :) = 0.0d+0  ! 1/cm^3-s

    ! Initial value (assigned via get_geometry subroutine)
    num_ind_cells = 0

    ! Legendre Gauss Quadrature values over chosen ordinates
    call legendre_gauss_quad(num_ords, -1.0d+0, 1.0d+0, ordinates, weights)
    mu = ordinates(num_ords:1:-1)
    weights = weights(num_ords:1:-1)
    tolerance = 1.0e-8

    ! Calculation: iterations
    cont_calc = .true.
    ! Start counter at zero
    iterations = int(0, 8)
    do while (cont_calc)
        phi_old = phi_new  ! 1/cm^2-s-MeV

        ! Fill Markovian geometry
        ! Allocations and deallocations of delta_x and materials should be
        ! Handled by the subroutine
        call get_geometry(delta_x, x_points, materials, chord_a, chord_b, &
                          thickness, num_ind_cells)

        ! Allocate and define properties
        allocate(macro_scat(num_groups, num_groups, num_ind_cells))
        macro_scat(:, :, :) = scat_const  ! 1/cm
        allocate(macro_fis(num_groups, num_cells))
        macro_fis(:, :) = fis_const  ! 1/cm
        allocate(macro_tot(num_groups, num_cells))
        macro_tot(:, :) = tot_const  ! 1/cm

        ! Allocate and define calculational arrays
        allocate(psi(num_groups, num_ind_cells, num_ords))
        psi(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_refl(num_groups, num_ind_cells, num_ords))
        psi_refl(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_p(num_groups, num_ind_cells, num_ords))
        psi_i_p = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_m(num_groups, num_ind_cells, num_ords))
        psi_i_m = 0.0d+0  ! 1/cm^2-s-MeV-strad

        ! Allocate phi calculations
        allocate(phi_morph(num_groups, num_ind_cells))
        call struct_to_unstruct(phi_new, )

        ! Determine sources for each cell and group
        do c = 1, num_ind_cells
            do g = 1, num_groups
                ! Scatter into, fission into, in-group scattering
                scatter_into = 0.0d+0  ! neutrons
                fission_into = 0.0d+0  ! neutrons
                cons_distance = 0.0d+0  ! cm
                ! For generating sources, go until the point is found (mapping)
                point_found = .false.
                do g_prime = 1, num_groups
                    scatter_into = scatter_into + macro_scat(g_prime, g, c) &
                                   * phi_new(g_prime, c)
                    fission_into = fission_into + chi(g) * nu(g_prime) &
                                   * macro_fis(g_prime, c) * phi_new(g_prime, c)
                end do
                scat_source(g, c) = scatter_into
                fis_source(g, c) = fission_into
                tot_source(g, c) = scat_source(g, c) + fis_source(g, c) &
                                   + spont_source(g, c)
            end do
        end do

        ! Forward sweep (left to right)
        ! First cell (left boundary)
        do g = 1, num_groups
            ! Ordinate loop, only consider the pos. ords for forward motion
            do m = (num_ords / 2 + 1), num_ords
                ! Lewis and Miller Eq. 3-40
                psi(g, 1, m) = (1.0d+0 + (macro_tot(g, 1) * delta_x(1)) &
                                / (2.0d+0 * dabs(mu(m))))**(-1) &
                               * (psi_bound_l(g, m) + psi_refl(g, 1, m) &
                                  + (tot_source(g, 1) * delta_x(1)) &
                                  / (2.d0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-41
                psi_i_p(g, 1, m) = 2.0d+0 * psi(g, 1, m) - psi_bound_l(g, m) &
                                   - psi_refl(g, 1, m)
            end do
        end do
        ! Rest of the cells (sans left bounding cell)
        do c = 2, num_cells
            do g = 1, num_groups
                do m = (num_ords / 2 + 1), num_ords
                    ! Continuity of boundaries
                    psi_i_m(g, c, m) = psi_i_p(g, c-1, m)
                    ! Lewis and Miller Eq. 3-40
                    psi(g, c, m) = (1.0d+0 + (macro_tot(g, c) * delta_x(c)) &
                                    / (2.0d+0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_m(g, c, m) + (tot_source(g, c) &
                                                          * delta_x(c)) &
                                      / (2.0d+0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-41
                    psi_i_p(g, c, m) = 2.0d+0 * psi(g, c, m) - psi_i_m(g, c, m)
                end do
            end do
        end do

        ! Backward sweep (right to left)
        ! First cell (right boundary)
        do g = 1, num_groups
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                ! Lewis and Miller Eq. 3-42
                psi(g, num_cells, m) = (1.0d+0 + (macro_tot(g, num_cells) &
                                                * delta_x(num_cells)) &
                                        / (2.0d+0 * dabs(mu(m))))**(-1) &
                                       * (psi_bound_r(g, m) &
                                          + psi_refl(g, num_cells, m) &
                                          + (tot_source(g, num_cells) &
                                             * delta_x(num_cells)) &
                                          / (2.0d+0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-43
                psi_i_m(g, num_cells, m) = 2.0d+0 * psi(g, num_cells, m) &
                                           - psi_bound_r(g, m) &
                                           - psi_refl(g, num_cells, m)
            end do
        end do
        ! Rest of the cells (sans right bounding cell)
        do c = (num_cells - 1), 1, -1
            do g = 1, num_groups
                do m = 1, (num_ords / 2)
                    ! Continuation of boundaries
                    psi_i_p(g, c, m) = psi_i_m(g, c+1, m)
                    ! Lewis and Miller Eq. 3-42
                    psi(g, c, m) = (1.0d+0 + (macro_tot(g, c) * delta_x(c)) &
                                    / (2.0d+0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_p(g, c, m) + (tot_source(g, c) &
                                                          * delta_x(c)) &
                                      / (2.0d+0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-43
                    psi_i_m(g, c, m) = 2.0d+0 * psi(g, c, m) - psi_i_p(g, c, m)
                end do
            end do
        end do

        ! Reflected angular fluxes at the left and right boundaries
        do g = 1, num_groups
            do m = 1, num_ords
                if (m >= num_ords / 2 + 1) then
                    ! Last cell, pos. angular flux reflected to neg. direction
                    psi_refl(g, num_cells, num_ords-m+1) = alpha_r &
                                                           * psi_i_p(g, num_cells, m)
                else
                    ! First cell, neg. angular flux reflected to pos. direction
                    psi_refl(g, 1, num_ords-m+1) = alpha_l &
                                                   * psi_i_m(g, 1, m)
                end if
            end do
        end do

        ! Calculate phi from psi
        ! Lewis and Miller Eq. 3-5
        do c=1,num_cells
            do g = 1,num_groups
                weighted_sum = 0.0d+0
                do m=1,num_ords
                    weighted_sum = weighted_sum + weights(m) * psi(g, c, m)
                end do
                phi_new(g, c) = 0.50d+0 * weighted_sum
            end do
        end do

        ! Relative error
        iterations = iterations + int(1, 8)
        err = maxval(dabs((phi_new - phi_old)) / phi_new)
        if (err <= tolerance) then
            cont_calc = .false.
            print *, "Quit after ", iterations , " iterations"
        else if (iterations > num_iter) then
            cont_calc = .false.
            print *, "Quit after maximum ", iterations, " iterations"
        end if

        ! Deallocations
        deallocate(macro_scat)
        deallocate(macro_fis)
        deallocate(macro_tot)
        deallocate(psi)
        deallocate(psi_refl)
        deallocate(psi_i_p)
        deallocate(psi_i_m)
    end do

    ! Create plot
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/steady_state_slab.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7,*) cell_vector(i), phi_new(1,i)
    end do
    close(7)
end program
