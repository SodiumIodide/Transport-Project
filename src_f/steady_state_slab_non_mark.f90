program steady_state_slab
    use self_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter = 5000
    integer, parameter :: &
        num_ords = 16, &
        num_cells = int(1e4, 4)

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+0, &  ! cm
        scat_const = 0.2d+0, &  ! 1/cm
        tot_const = 1.0d+0  ! 1/cm

    ! Material parameters
    real(8), dimension(num_cells) :: &
        delta_x, macro_scat, macro_tot

    ! Calculation variables
    integer(8) :: &
        iterations
    integer :: &
        i, c, m
    logical :: &
        cont_calc
    real(8) :: &
        tolerance, weighted_sum, err, total_abs, &
        leakage_l, leakage_r, balance_source_l, balance_source_r, source_dist, balance
    real(8), dimension(num_cells) :: &
        phi_new, phi_old, scat_source, spont_source, tot_source
    real(8), dimension(num_cells, num_ords) :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu, psi_bound_l, psi_bound_r

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector

    ! Assigment of material parameters
    delta_x(:) = thickness / dble(num_cells)  ! cm
    ! Assignment of initial calculation variables
    phi_new(:) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old(:) = 0.0d+0  ! 1/cm^2-s-MeV
    psi(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad, angular neutron flux
    psi_i_p(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    psi_i_m(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    scat_source(:) = 0.0d+0  ! 1/cm^3-s
    spont_source(:) = 0.0d+0  ! 1/cm^3-s
    tot_source(:) = 0.0d+0  ! 1/cm^3-s
    macro_scat(:) = scat_const  ! 1/cm
    macro_tot(:) = tot_const  ! 1/cm

    ! Legendre Gauss Quadrature values
    call legendre_gauss_quad(num_ords, -1.0d+0, 1.0d+0, ordinates, weights)
    !call gauleg(-1.0d+0, 1.0d+0, ordinates, weights, num_ords)
    !mu = ordinates
    mu = ordinates(num_ords:1:-1)
    weights = weights(num_ords:1:-1)

    do i = 1, num_ords
        print *, mu(i), weights(i)
    end do

    ! Boundary conditions: 1/cm^2-s-MeV-strad
    ! Left boundary
    !psi_bound_l(:) = 1.0d+0  ! Isotropic source
    psi_bound_l(:) = 0.0d+0  ! Vacuum
    ! Beam source (in conj. with vacuum):
    psi_bound_l(num_ords) = 1.0d+0 / (mu(num_ords) * weights(num_ords))

    ! Right boundary
    !psi_bound_r(:) = 1.0d+0  ! Isotropic source
    psi_bound_r(:) = 0.0d+0  ! Vacuum
    ! Beam source (in conj. with vacuum):
    !psi_bound_r(1) = 1.0d+0 / (mu(1) * weights(1))

    ! Tolerance for ending calculation
    tolerance = epsilon(1.0d+0)

    ! Calculation: iterations
    cont_calc = .true.
    ! Start counter at zero
    iterations = int(0, 8)
    do while (cont_calc)
        phi_old = phi_new  ! 1/cm^2-s-MeV

        ! Determine sources for each cell and group
        do c = 1, num_cells
            ! Scatter into, one-group scattering
            scat_source(c) = macro_scat(c) / 2.0d+0 * phi_new(c)
            tot_source(c) = scat_source(c) + spont_source(c) / 2.0d+0
        end do

        ! Forward sweep (left to right)
        ! First cell (left boundary)
        ! Ordinate loop, only consider the pos. ords for forward motion
        do m = (num_ords / 2 + 1), num_ords
            psi(1, m) = (1.0d+0 + (macro_tot(1) * delta_x(1)) &
                         / (2.0d+0 * dabs(mu(m))))**(-1) &
                        * (psi_bound_l(m) &
                           + (tot_source(1) * delta_x(1)) &
                             / (2.0d+0 * dabs(mu(m))))
            psi_i_p(1, m) = 2.0d+0 * psi(1, m) - psi_bound_l(m)
        end do
        ! Rest of the cells (sans left bounding cell)
        do c = 2, num_cells
            do m = (num_ords / 2 + 1), num_ords
                ! Continuity of boundaries
                psi_i_m(c, m) = psi_i_p(c - 1, m)
                psi(c, m) = (1.0d+0 + (macro_tot(c) * delta_x(c)) &
                             / (2.0d+0 * dabs(mu(m))))**(-1) &
                            * (psi_i_m(c, m) &
                               + (tot_source(c) * delta_x(c)) &
                                 / (2.0d+0 * dabs(mu(m))))
                psi_i_p(c, m) = 2.0d+0 * psi(c, m) - psi_i_m(c, m)
            end do
        end do

        ! Backward sweep (right to left)
        ! First cell (right boundary)
        ! Ordinate loop, only consider neg. ords for backwards motion
        do m = 1, (num_ords / 2)
            ! Lewis and Miller Eq. 3-42
            psi(num_cells, m) = (1.0d+0 + (macro_tot(num_cells) * delta_x(num_cells)) &
                                 / (2.0d+0 * dabs(mu(m))))**(-1) &
                                * (psi_bound_r(m) + (tot_source(num_cells) &
                                                     * delta_x(num_cells)) &
                                   / (2.0d+0 * dabs(mu(m))))
            ! Lewis and Miller Eq. 3-43
            psi_i_m(num_cells, m) = 2.0d+0 * psi(num_cells, m) - psi_bound_r(m)
        end do
        ! Rest of the cells (sans right bounding cell)
        do c = (num_cells - 1), 1, -1
            do m = 1, (num_ords / 2)
                ! Continuation of boundaries
                psi_i_p(c, m) = psi_i_m(c + 1, m)
                psi(c, m) = (1.0d+0 + (macro_tot(c) * delta_x(c)) &
                             / (2.0d+0 * dabs(mu(m))))**(-1) &
                            * (psi_i_p(c, m) + (tot_source(c) &
                                                * delta_x(c)) &
                               / (2.0d+0 * dabs(mu(m))))
                psi_i_m(c, m) = 2.0d+0 * psi(c, m) - psi_i_p(c, m)
            end do
        end do

        ! Calculate phi from psi
        do c = 1, num_cells
            weighted_sum = 0.0d+0
            do m = 1, num_ords
                weighted_sum = weighted_sum + weights(m) * psi(c, m)
            end do
            phi_new(c) = weighted_sum  ! 1/cm^2-s-MeV
        end do

        ! Relative error
        iterations = iterations + int(1, 8)
        err = maxval(dabs((phi_new - phi_old)) / phi_new)
        if (err <= tolerance) then
            cont_calc = .false.
            print *, "Quit after ", iterations , " iterations"
        else if (iterations > num_iter) then
            cont_calc = .false.
            print *, "Quit after ", iterations, " iterations"
        end if
    end do

    ! Check balances
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    total_abs = 0.0d+0
    balance_source_l = 0.0d+0
    balance_source_r = 0.0d+0
    source_dist = 0.0d+0
    ! Tally the boundary sources (currents)
    do m = 1, (num_ords / 2)
        balance_source_r = balance_source_r + weights(m) * dabs(mu(m)) * psi_bound_r(m)
    end do
    do m = (num_ords / 2 + 1), num_ords
        balance_source_l = balance_source_l + weights(m) * dabs(mu(m)) * psi_bound_l(m)
    end do
    ! Tally the overall losses due to leakage and absorption, as well as the distributed source
    do c = 1, num_cells
        ! Leakage in neg. direction from left face
        if (c == 1) then
            do m = 1, (num_ords / 2)
                leakage_l = leakage_l + dabs(mu(m)) * weights(m) * psi_i_m(c, m)
            end do
        ! Leakage in pos. direction from right face
        else if (c == num_cells) then
            do m = (num_ords / 2 + 1), num_ords
                leakage_r = leakage_r + dabs(mu(m)) * weights(m) * psi_i_p(c, m)
            end do
        end if
        ! Total absorption in system
        total_abs = total_abs + (macro_tot(c) - macro_scat(c)) * delta_x(c) * phi_new(c)
        ! Distributed source
        source_dist = source_dist + spont_source(c) * delta_x(c)
    end do
    print *, "Leakage left: ", leakage_l
    print *, "Leakage right: ", leakage_r
    print *, "Source left: ", balance_source_l
    print *, "Source right: ", balance_source_r
    print *, "Distributed source: ", source_dist
    print *, "Absorption loss: ", total_abs
    print *, "Source is ", balance_source_l + balance_source_r + source_dist
    print *, "Loss is ", leakage_l + leakage_r + total_abs
    balance = balance_source_l + balance_source_r + source_dist - leakage_l - leakage_r - total_abs
    print *, "Balance (source - loss) is ", balance

    ! Create plot
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/steady_state_slab.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7,*) cell_vector(i), phi_new(i)
    end do
    close(7)
end program
