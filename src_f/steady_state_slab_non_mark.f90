program steady_state_slab
    use self_library
    use mcnp_random

    implicit none

    ! Constant parameters
    integer, parameter :: &
        num_cells = 100, num_groups = 1, num_ords = 16, num_iter = 5000

    ! Material parameters
    double precision :: &
        thickness, alpha_l, alpha_r
    double precision, dimension(num_cells) :: &
        delta_x
    double precision, dimension(num_groups, num_groups, num_cells) :: &
        macro_scat
    double precision, dimension(num_groups, num_cells) :: &
        macro_fis, macro_tot
    double precision, dimension(num_groups) :: &
        chi, nu

    ! Calculation variables
    integer :: &
        i, c, g, g_prime, m, iterations
    logical :: &
        cont_calc
    double precision :: &
        tolerance, scatter_into, fission_into, weighted_sum, err
    double precision, dimension(num_groups, num_cells) :: &
        phi_new, phi_old, scat_source, fis_source, spont_source, tot_source
    double precision, dimension(num_groups, num_cells, num_ords) :: &
        psi, psi_refl, psi_i_p, psi_i_m
    double precision, dimension(num_groups, num_ords) :: &
        psi_bound_l, psi_bound_r
    double precision, dimension(num_ords) :: &
        ordinates, weights, mu

    ! Additional variables (plotting, etc.)
    double precision, dimension(num_cells) :: &
        cell_vector

    ! Assigment of material parameters
    thickness = 1.d0  ! cm
    delta_x(:) = thickness / dble(num_cells)  ! cm
    macro_scat(:, :, :) = 0.2d0 / dble(num_groups)  ! cm^-1
    macro_fis(:, :) = 0.d0  ! cm^-1
    macro_tot(:, :) = 1.0d0  ! cm^-1
    chi(:) = 1.d0 / dble(num_groups)  ! fission energy spectrum
    nu(:) = 2.4144  ! average number of neutrons emitted per fission
    ! Assignment of initial calculation variables
    phi_new(:, :) = 1.d0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old(:, :) = 0.d0  ! 1/cm^2-s-MeV
    psi(:, :, :) = 0.d0  ! 1/cm^2-s-MeV-strad, angular neutron flux
    psi_refl(:, :, :) = 0.d0  ! 1/cm^2-s-MeV-strad
    psi_i_p(:, :, :) = 0.d0  ! 1/cm^2-s-MeV-strad
    psi_i_m(:, :, :) = 0.d0  ! 1/cm^2-s-MeV-strad
    psi_bound_l(:, :) = 1.d0  ! 1/cm^2-s-MeV-strad
    psi_bound_r(:, :) = 0.d0  ! 1/cm^2-s-MeV-strad
    scat_source(:, :) = 0.0d0  ! 1/cm^3-s
    fis_source(:, :) = 0.d0  ! 1/cm^3-s
    spont_source(:, :) = 0.d0  ! 1/cm^3-s
    tot_source(:, :) = 0.d0  ! 1/cm^3-s
    alpha_l = 0.d0  ! Albedo boundary condition, left
    alpha_r = 0.d0  ! Albedo boundary condition, right
    ! alpha = 1.d0 for perfect reflection, alpha = 0.d0 for zero reflection

    ! Legendre Gauss Quadrature values
    call legendre_gauss_quad(num_ords, -1.d0, 1.d0, ordinates, weights)
    mu = ordinates(num_ords:1:-1)
    weights = weights(num_ords:1:-1)
    tolerance = 1e-8

    ! Calculation: iterations
    cont_calc = .true.
    iterations = 0
    do while (cont_calc)
        phi_old = phi_new  ! 1/cm^2-s-MeV
        ! Determine sources for each cell and group
        do c = 1, num_cells
            do g = 1, num_groups
                ! Scatter into, fission into, in-group scattering
                scatter_into = 0.d0  ! neutrons
                fission_into = 0.d0  ! neutrons
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
                psi(g, 1, m) = (1.d0 + (macro_tot(g, 1) * delta_x(1)) &
                                / (2.d0 * dabs(mu(m))))**(-1) &
                               * (psi_bound_l(g, m) + psi_refl(g, 1, m) &
                                  + (tot_source(g, 1) * delta_x(1)) &
                                  / (2.d0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-41
                psi_i_p(g, 1, m) = 2.d0 * psi(g, 1, m) - psi_bound_l(g, m) &
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
                    psi(g, c, m) = (1.d0 + (macro_tot(g, c) * delta_x(c)) &
                                    / (2.d0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_m(g, c, m) + (tot_source(g, c) &
                                                          * delta_x(c)) &
                                      / (2.d0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-41
                    psi_i_p(g, c, m) = 2.d0 * psi(g, c, m) - psi_i_m(g, c, m)
                end do
            end do
        end do
        ! Backward sweep (right to left)
        ! First cell (right boundary)
        do g = 1, num_groups
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                ! Lewis and Miller Eq. 3-42
                psi(g, num_cells, m) = (1.d0 + (macro_tot(g, num_cells) &
                                                * delta_x(num_cells)) &
                                        / (2.d0 * dabs(mu(m))))**(-1) &
                                       * (psi_bound_r(g, m) &
                                          + psi_refl(g, num_cells, m) &
                                          + (tot_source(g, num_cells) &
                                             * delta_x(num_cells)) &
                                          / (2.d0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-43
                psi_i_m(g, num_cells, m) = 2.d0 * psi(g, num_cells, m) &
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
                    psi(g, c, m) = (1.d0 + (macro_tot(g, c) * delta_x(c)) &
                                    / (2.d0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_p(g, c, m) + (tot_source(g, c) &
                                                          * delta_x(c)) &
                                      / (2.d0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-43
                    psi_i_m(g, c, m) = 2.d0 * psi(g, c, m) - psi_i_p(g, c, m)
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
                weighted_sum = 0.d0
                do m=1,num_ords
                    weighted_sum = weighted_sum + weights(m) * psi(g, c, m)
                end do
                phi_new(g, c) = 0.5d0 * weighted_sum
            end do
        end do
        ! Relative error
        iterations = iterations + 1
        err = maxval(dabs((phi_new - phi_old)) / phi_new)
        if (err <= tolerance) then
            cont_calc = .false.
            print *, "Quit after ", iterations , " iterations"
        else if (iterations > num_iter) then
            cont_calc = .false.
            print *, "Quit after ", iterations, " iterations"
        end if
    end do

    ! Create plot
    call linspace(cell_vector, 0.d0, thickness, num_cells)
    open(unit=7, file="./out/steady_state_slab.out", form="formatted", &
         status="replace", action="write")
    do i=1,num_cells
        write(7,*) cell_vector(i), phi_new(1,i)
    end do
    close(7)
end program
