program steady_state_slab_closure_simul
    use self_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter = 5000
    integer, parameter :: &
        num_ords = 16, &
        num_cells = int(5.0d+2, 4), &
        num_materials = 2
    logical, parameter :: &
        use_alpha = .true.

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1  ! cm

    ! Individual material properties
    real(8), dimension(num_materials), parameter :: &
        tot_const = (/dble(2)/dble(101), dble(200)/dble(101)/), &  ! 1/cm
        scat_const = (/dble(2)/dble(101)*0.9d+0, dble(200)/dble(101)*0.9d+0/), &  ! 1/cm
        abs_const = tot_const - scat_const, &  ! 1/cm
        chord = (/dble(101)/dble(20), dble(101)/dble(20)/)  ! cm
    real(8), dimension(num_materials) :: &
        prob

    ! Material parameters
    real(8), dimension(num_cells) :: &
        delta_x
    real(8), dimension(num_cells, num_materials) :: &
        macro_scat, macro_tot

    ! Calculation variables
    integer(8) :: &
        iterations
    integer :: &
        i, c, m, k, ok
    logical :: &
        cont_calc
    real(8) :: &
        tolerance, weighted_sum, err, total_chord, &
        leakage_l, leakage_r, alpha, sig_a_av, sig_t_av
    real(8), dimension(num_cells, num_materials) :: &
        phi_new, phi_old, spont_source, scat_source, tot_source
    real(8), dimension(num_cells, num_ords, num_materials) :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_cells, num_ords) :: &
        psi_overall
    real(8), dimension(num_cells) :: &
        phi_overall
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu, psi_bound_l, psi_bound_r
    real(8), dimension(num_materials) :: &
        flux_solution, internal_source, streaming_source
    real(8), dimension(num_materials, num_materials) :: &
        multiplier_mat
    integer, dimension(num_materials) :: &
        pivot

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector
    character(:), allocatable :: &
        filename

    ! Assigment of material parameters
    delta_x(:) = thickness / dble(num_cells)  ! cm
    ! Assignment of initial calculation variables
    phi_new(:, :) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    psi(:, :, :) = 1.0d+0  ! 1/cm^2-s-MeV-strad, angular neutron flux
    psi_i_p(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    psi_i_m(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    scat_source(:, :) = 0.0d+0  ! 1/cm^3-s
    tot_source(:, :) = 0.0d+0  ! 1/cm^3-s
    internal_source(:) = 0.0d+0
    streaming_source(:) = 0.0d+0
    multiplier_mat(:, :) = 0.0d+0
    flux_solution(:) = 0.0d+0
    total_chord = 0.0d+0  ! cm
    do k = 1, num_materials
        macro_scat(:, k) = scat_const(k)  ! 1/cm
        macro_tot(:, k) = tot_const(k)  ! 1/cm
        spont_source(:, k) = 0.0d+0  ! 1/cm^3-s
        total_chord = total_chord + chord(k)  ! cm
    end do
    prob(:) = chord(:) / total_chord
    sig_a_av = prob(1) * abs_const(1) + prob(2) * abs_const(2)
    sig_t_av = prob(1) * tot_const(1) + prob(2) * tot_const(2)
    if (use_alpha) then
        alpha = dsqrt(sig_a_av * sig_t_av) &
            * (sig_a_av * (tot_const(2) - tot_const(1))**2 + sig_t_av * (abs_const(2) - abs_const(1))**2) &
            / ((sig_a_av * (tot_const(2) - tot_const(1)))**2 + (sig_t_av * (abs_const(2) - abs_const(1)))**2)
    else
        alpha = 1.0d+0
    end if
    print *, "Alpha: ", alpha

    ! Legendre Gauss Quadrature values
    call legendre_gauss_quad(num_ords, -1.0d+0, 1.0d+0, ordinates, weights)
    mu = ordinates(num_ords:1:-1)
    weights = weights(num_ords:1:-1)

    ! Boundary conditions: 1/cm^2-s-MeV-strad
    ! Left boundary
    psi_bound_l(:) = 2.0d+0  ! Isotropic source
    !psi_bound_l(:) = 0.0d+0  ! Vacuum
    ! Beam source (in conj. with vacuum):
    !psi_bound_l(num_ords) = 1.0d+0 / (mu(num_ords) * weights(num_ords))

    ! Right boundary
    !psi_bound_r(:) = 2.0d+0  ! Isotropic source
    psi_bound_r(:) = 0.0d+0  ! Vacuum
    ! Beam source (in conj. with vacuum):
    !psi_bound_r(1) = 1.0d+0 / (mu(1) * weights(1))

    ! Tolerance for ending calculation
    tolerance = 1.0d-8

    psi_overall(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    phi_overall(:) = 0.0d+0  ! 1/cm^2-s-MeV

    iterations = int(0, 8)
    cont_calc = .true.
    ! Inner loop
    do while (cont_calc)
        phi_old = phi_new  ! 1/cm^2-s-MeV

        ! Determine sources for each cell
        do c = 1, num_cells
            do k = 1, num_materials
                ! Scatter into, one-group scattering
                scat_source(c, k) = macro_scat(c, k) / 2.0d+0 * phi_new(c, k)
                tot_source(c, k) = scat_source(c, k) + spont_source(c, k) / 2.0d+0
            end do
        end do

        ! Forward sweep (left to right)
        ! First cell (left boundary)
        ! Ordinate loop, only consider the pos. ords for forward motion
        do m = (num_ords / 2 + 1), num_ords
            do k = 1, num_materials
                internal_source(k) = tot_source(1, k) * delta_x(1) / (2.0d+0 * dabs(mu(m)))
                streaming_source(k) = psi_bound_l(m)
                flux_solution(k) = internal_source(k) + streaming_source(k)
            end do
            multiplier_mat(1, 1) = 1.0d+0 + delta_x(1) * macro_tot(1, 1) / (2.0d+0 * dabs(mu(m))) &
                + alpha * delta_x(1) / (2.0d+0 * chord(1))
            multiplier_mat(1, 2) = -delta_x(1) * alpha * prob(2) / (2.0d+0 * prob(1) * chord(2))
            multiplier_mat(2, 2) = 1.0d+0 + delta_x(1) * macro_tot(1, 2) / (2.0d+0 * dabs(mu(m))) &
                + alpha * delta_x(1) / (2.0d+0 * chord(2))
            multiplier_mat(2, 1) = -delta_x(1) * alpha * prob(1) / (2.0d+0 * prob(2) * chord(1))
            call dgesv(2, 1, multiplier_mat, 2, pivot, flux_solution, 2, ok)
            do k = 1, num_materials
                psi(1, m, k) = flux_solution(k)
                psi_i_p(1, m, k) = 2.0d+0 * psi(1, m, k) - psi_bound_l(m)
            end do
        end do
        ! Rest of the cells (sans left bounding cell)
        do c = 2, num_cells
            do m = (num_ords / 2 + 1), num_ords
                do k = 1, num_materials
                    ! Continuity of boundaries
                    psi_i_m(c, m, k) = psi_i_p(c - 1, m, k)
                    internal_source(k) = tot_source(c, k) * delta_x(c) / (2.0d+0 * dabs(mu(m)))
                    streaming_source(k) = psi_i_m(c, m, k)
                    flux_solution(k) = internal_source(k) + streaming_source(k)
                end do
                multiplier_mat(1, 1) = 1.0d+0 + delta_x(c) * macro_tot(c, 1) / (2.0d+0 * dabs(mu(m))) &
                    + alpha * delta_x(c) / (2.0d+0 * chord(1))
                multiplier_mat(1, 2) = -delta_x(c) * alpha * prob(2) / (2.0d+0 * prob(1) * chord(2))
                multiplier_mat(2, 2) = 1.0d+0 + delta_x(c) * macro_tot(c, 2) / (2.0d+0 * dabs(mu(m))) &
                    + alpha * delta_x(c) / (2.0d+0 * chord(2))
                multiplier_mat(2, 1) = -delta_x(c) * alpha * prob(1) / (2.0d+0 * prob(2) * chord(1))
                call dgesv(2, 1, multiplier_mat, 2, pivot, flux_solution, 2, ok)
                do k = 1, num_materials
                    psi(c, m, k) = flux_solution(k)
                    psi_i_p(c, m, k) = 2.0d+0 * psi(c, m, k) - psi_i_m(c, m, k)
                end do
            end do
        end do

        ! Backward sweep (right to left)
        ! First cell (right boundary)
        ! Ordinate loop, only consider neg. ords for backwards motion
        do m = 1, (num_ords / 2)
            do k = 1, num_materials
                internal_source(k) = tot_source(num_cells, k) * delta_x(num_cells) / (2.0d+0 * dabs(mu(m)))
                streaming_source(k) = psi_bound_r(m)
                flux_solution(k) = internal_source(k) + streaming_source(k)
            end do
            multiplier_mat(1, 1) = 1.0d+0 + delta_x(num_cells) * macro_tot(num_cells, 1) / (2.0d+0 * dabs(mu(m))) &
                + alpha * delta_x(num_cells) / (2.0d+0 * chord(1))
            multiplier_mat(1, 2) = -delta_x(num_cells) * alpha * prob(2) / (2.0d+0 * prob(1) * chord(2))
            multiplier_mat(2, 2) = 1.0d+0 + delta_x(num_cells) * macro_tot(num_cells, 2) / (2.0d+0 * dabs(mu(m))) &
                + alpha * delta_x(num_cells) / (2.0d+0 * chord(2))
            multiplier_mat(2, 1) = -delta_x(num_cells) * alpha * prob(1) / (2.0d+0 * prob(2) * chord(1))
            call dgesv(2, 1, multiplier_mat, 2, pivot, flux_solution, 2, ok)
            do k = 1, num_materials
                psi(num_cells, m, k) = flux_solution(k)
                psi_i_m(num_cells, m, k) = 2.0d+0 * psi(num_cells, m, k) - psi_bound_r(m)
            end do
        end do
        ! Rest of the cells (sans right bounding cell)
        do c = (num_cells - 1), 1, -1
            do m = 1, (num_ords / 2)
                do k = 1, num_materials
                    ! Continuation of boundaries
                    psi_i_p(c, m, k) = psi_i_m(c + 1, m, k)
                    internal_source(k) = tot_source(c, k) * delta_x(c) / (2.0d+0 * dabs(mu(m)))
                    streaming_source(k) = psi_i_p(c, m, k)
                    flux_solution(k) = internal_source(k) + streaming_source(k)
                end do
                multiplier_mat(1, 1) = 1.0d+0 + delta_x(c) * macro_tot(c, 1) / (2.0d+0 * dabs(mu(m))) &
                    + alpha * delta_x(c) / (2.0d+0 * chord(1))
                multiplier_mat(1, 2) = -delta_x(c) * alpha * prob(2) / (2.0d+0 * prob(1) * chord(2))
                multiplier_mat(2, 2) = 1.0d+0 + delta_x(c) * macro_tot(c, 2) / (2.0d+0 * dabs(mu(m))) &
                    + alpha * delta_x(c) / (2.0d+0 * chord(2))
                multiplier_mat(2, 1) = -delta_x(c) * alpha * prob(1) / (2.0d+0 * prob(2) * chord(1))
                call dgesv(2, 1, multiplier_mat, 2, pivot, flux_solution, 2, ok)
                do k = 1, num_materials
                    psi(c, m, k) = flux_solution(k)
                    psi_i_m(c, m, k) = 2.0d+0 * psi(c, m, k) - psi_i_p(c, m, k)
                end do
            end do
        end do

        ! Calculate phi from the individual psi calculations
        do c = 1, num_cells
            do k = 1, num_materials
                weighted_sum = 0.0d+0
                do m = 1, num_ords
                    weighted_sum = weighted_sum + weights(m) * psi(c, m, k)
                end do
                phi_new(c, k) = weighted_sum  ! 1/cm^2-s-MeV
            end do
        end do

        ! Relative error for the inner loop
        iterations = iterations + int(1, 8)
        err = maxval(dabs((phi_new - phi_old)) / phi_new)
        if (err <= tolerance) then
            cont_calc = .false.
            print *, "Converged after ", iterations, " iterations"
        else if (iterations > num_iter) then
            cont_calc = .false.
            print *, "No convergence after ", iterations
        end if
    end do

    ! Calculate the overall psi
    do c = 1, num_cells
        do m = 1, num_ords
            do k = 1, num_materials
                psi_overall(c, m) = psi_overall(c, m) + prob(k) * psi(c, m, k)
            end do
        end do
    end do

    ! Calculate phi from the overall psi
    do c = 1, num_cells
        weighted_sum = 0.0d+0
        do m = 1, num_ords
            weighted_sum = weighted_sum + weights(m) * psi_overall(c, m)
        end do
        phi_overall(c) = weighted_sum  ! 1/cm^2-s-MeV
    end do

    ! Check balances
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    ! Tally the overall losses due to leakage and absorption, as well as the distributed source
    do k = 1, num_materials
        ! Leakage in neg. direction from left face
        do m = 1, (num_ords / 2)
            leakage_l = leakage_l + dabs(mu(m)) * weights(m) * prob(k) * psi_i_m(1, m, k)
        end do
        ! Leakage in pos. direction from right face
        do m = (num_ords / 2 + 1), num_ords
            leakage_r = leakage_r + dabs(mu(m)) * weights(m) * prob(k) * psi_i_p(num_cells, m, k)
        end do
    end do
    print *, "Reflection from left: ", leakage_l
    print *, "Transmission on right: ", leakage_r

    ! Create plot
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)
    if (use_alpha) then
        filename = "./out/steady_state_slab_closure_alpha.out"
    else
        filename = "./out/steady_state_slab_closure.out"
    end if
    open(unit=7, file=filename, form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7, *) cell_vector(i), phi_overall(i)
    end do
    close(7)
    if (use_alpha) then
        filename = "./out/steady_state_slab_closure_alpha_1.out"
    else
        filename = "./out/steady_state_slab_closure_1.out"
    end if
    open(unit=8, file=filename, form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(8, *) cell_vector(i), phi_new(i, 1)
    end do
    close(8)
    if (use_alpha) then
        filename = "./out/steady_state_slab_closure_alpha_2.out"
    else
        filename = "./out/steady_state_slab_closure_2.out"
    end if
    open(unit=9, file=filename, form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(9, *) cell_vector(i), phi_new(i, 2)
    end do
    close(9)

    if (allocated(filename)) then
        deallocate(filename)
    end if
end program steady_state_slab_closure_simul
