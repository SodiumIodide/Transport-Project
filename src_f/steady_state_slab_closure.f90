program steady_state_slab_closure
    use self_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter_inner = 5000, &
        num_iter_outer = 100000
    integer, parameter :: &
        num_ords = 16, &
        num_cells = int(1e2, 4), &
        num_materials = 2

    ! Material properties
    real(8), parameter :: &
        thickness = 10.0d+0  ! cm

    ! Individual material properties
    real(8), dimension(num_materials), parameter :: &
        scat_const = (/0.2d+0, 0.3d+0/), &  ! 1/cm
        tot_const = (/1.0d+0, 1.0d+0/), &  ! 1/cm
        chord = (/0.5, 0.4/)  ! cm
    real(8), dimension(num_materials) :: &
        prob

    ! Material parameters
    real(8), dimension(num_cells) :: &
        delta_x
    real(8), dimension(num_cells, num_materials) :: &
        macro_scat, macro_tot

    ! Calculation variables
    integer(8) :: &
        iterations_inner, iterations_outer
    integer :: &
        i, c, m, k, material_num, material_off_num
    logical :: &
        cont_calc_inner, cont_calc_outer
    real(8) :: &
        tolerance_inner, tolerance_outer, weighted_sum, err, total_abs, total_chord, &
        leakage_l, leakage_r, balance_source_l, balance_source_r, source_dist, balance
    real(8), dimension(num_cells) :: &
        phi_new_outer, phi_old_outer, scat_source, tot_source
    real(8), dimension(num_cells, num_materials) :: &
        phi_new_inner, phi_old_inner, spont_source
    real(8), dimension(num_cells, num_ords, num_materials) :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_cells, num_ords) :: &
        psi_overall
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu, psi_bound_l, psi_bound_r

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector

    ! Assigment of material parameters
    delta_x(:) = thickness / dble(num_cells)  ! cm
    ! Assignment of initial calculation variables
    phi_new_outer(:) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old_outer(:) = 0.0d+0  ! 1/cm^2-s-MeV
    psi(:, :, :) = 1.0d+0  ! 1/cm^2-s-MeV-strad, angular neutron flux
    psi_i_p(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    psi_i_m(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    scat_source(:) = 0.0d+0  ! 1/cm^3-s
    tot_source(:) = 0.0d+0  ! 1/cm^3-s
    total_chord = 0.0d+0  ! cm
    do k = 1, num_materials
        macro_scat(:, k) = scat_const(k)  ! 1/cm
        macro_tot(:, k) = tot_const(k)  ! 1/cm
        spont_source(:, k) = 0.0d+0  ! 1/cm^3-s
        total_chord = total_chord + chord(k)  ! cm
    end do
    prob(:) = chord(:) / total_chord

    ! Legendre Gauss Quadrature values
    call legendre_gauss_quad(num_ords, -1.0d+0, 1.0d+0, ordinates, weights)
    mu = ordinates(num_ords:1:-1)
    weights = weights(num_ords:1:-1)

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
    tolerance_inner = epsilon(1.0d+0)
    tolerance_outer = epsilon(1.0d+0)

    ! Initial material number
    material_num = 1
    material_off_num = 2

    ! Calculation: iterations
    cont_calc_outer = .true.
    ! Start counter at zero
    iterations_outer = int(0, 8)
    do while (cont_calc_outer)
        phi_old_outer = phi_new_outer  ! 1/cm^2-s-MeV

        ! Switch value of two variables
        material_num = material_num + material_off_num
        material_off_num = material_num - material_off_num
        material_num = material_num - material_off_num

        phi_new_inner(:, :) = 1.0d+0  ! 1/cm^2-s-MeV
        phi_old_inner(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
        psi_overall(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

        iterations_inner = int(0, 8)
        cont_calc_inner = .true.
        ! Inner loop
        do while (cont_calc_inner)
            phi_old_inner = phi_new_inner  ! 1/cm^2-s-MeV

            ! Determine sources for each cell and group
            do c = 1, num_cells
                ! Scatter into, one-group scattering
                scat_source(c) = macro_scat(c, material_num) / 2.0d+0 * phi_new_inner(c, material_num)
                tot_source(c) = scat_source(c) + spont_source(c, material_num) / 2.0d+0
            end do

            ! Forward sweep (left to right)
            ! First cell (left boundary)
            ! Ordinate loop, only consider the pos. ords for forward motion
            do m = (num_ords / 2 + 1), num_ords
                psi(1, m, material_num) = (1.0d+0 + (macro_tot(1, material_num) * delta_x(1)) &
                    / (2.0d+0 * dabs(mu(m))) + delta_x(1) &
                    / (2.0d+0 * chord(material_num)))**(-1) &
                    * (psi_bound_l(m) + (tot_source(1) * delta_x(1)) &
                    / (2.0d+0 * dabs(mu(m))) &
                    + (delta_x(1) * prob(material_off_num) * psi(1, m, material_off_num)) &
                    / (2.0d+0 * prob(material_num) * chord(material_off_num)))
                psi_i_p(1, m, material_num) = 2.0d+0 * psi(1, m, material_num) - psi_bound_l(m)
            end do
            ! Rest of the cells (sans left bounding cell)
            do c = 2, num_cells
                do m = (num_ords / 2 + 1), num_ords
                    ! Continuity of boundaries
                    psi_i_m(c, m, material_num) = psi_i_p(c - 1, m, material_num)
                    psi(c, m, material_num) = (1.0d+0 + (macro_tot(c, material_num) * delta_x(c)) &
                        / (2.0d+0 * dabs(mu(m))) + delta_x(c) &
                        / (2.0d+0 * chord(material_num)))**(-1) &
                        * (psi_i_m(c, m, material_num) + (tot_source(c) * delta_x(c)) &
                        / (2.0d+0 * dabs(mu(m))) &
                        + (delta_x(c) * prob(material_off_num) * psi(c, m, material_off_num)) &
                        / (2.0d+0 * prob(material_num) * chord(material_off_num)))
                    psi_i_p(c, m, material_num) = 2.0d+0 * psi(c, m, material_num) - psi_i_m(c, m, material_num)
                end do
            end do

            ! Backward sweep (right to left)
            ! First cell (right boundary)
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                psi(num_cells, m, material_num) = (1.0d+0 + (macro_tot(num_cells, material_num) * delta_x(num_cells)) &
                    / (2.0d+0 * dabs(mu(m))) + delta_x(num_cells) &
                    / (2.0d+0 * chord(material_num)))**(-1) &
                    * (psi_bound_r(m) + (tot_source(num_cells) * delta_x(num_cells)) &
                    / (2.0d+0 * dabs(mu(m))) &
                    + (delta_x(num_cells) * prob(material_off_num) * psi(num_cells, m, material_off_num)) &
                    / (2.0d+0 * prob(material_num) * chord(material_off_num)))
                psi_i_m(num_cells, m, material_num) = 2.0d+0 * psi(num_cells, m, material_num) - psi_bound_r(m)
            end do
            ! Rest of the cells (sans right bounding cell)
            do c = (num_cells - 1), 1, -1
                do m = 1, (num_ords / 2)
                    ! Continuation of boundaries
                    psi_i_p(c, m, material_num) = psi_i_m(c + 1, m, material_num)
                    psi(c, m, material_num) = (1.0d+0 + (macro_tot(c, material_num) * delta_x(c)) &
                        / (2.0d+0 * dabs(mu(m))) + delta_x(c) &
                        / (2.0d+0 * chord(material_num)))**(-1) &
                        * (psi_i_p(c, m, material_num) + (tot_source(c) * delta_x(c)) &
                        / (2.0d+0 * dabs(mu(m))) &
                        + (delta_x(c) * prob(material_off_num) * psi(c, m, material_off_num)) &
                        / (2.0d+0 * prob(material_num) * chord(material_off_num)))
                    psi_i_m(c, m, material_num) = 2.0d+0 * psi(c, m, material_num) - psi_i_p(c, m, material_num)
                end do
            end do

            ! Calculate phi from the individual psi calculations
            do c = 1, num_cells
                do k = 1, num_materials
                    weighted_sum = 0.0d+0
                    do m = 1, num_ords
                        weighted_sum = weighted_sum + weights(m) * psi(c, m, k)
                    end do
                    phi_new_inner(c, k) = weighted_sum  ! 1/cm^2-s-MeV
                end do
            end do

            ! Relative error for the inner loop
            iterations_inner = iterations_inner + int(1, 8)
            err = maxval(dabs((phi_new_inner - phi_old_inner)) / phi_new_inner)
            if (err <= tolerance_inner) then
                cont_calc_inner = .false.
                !print *, "Converged after ", iterations_inner, " inner loops"
            else if (iterations_inner > num_iter_inner) then
                cont_calc_inner = .false.
                print *, "No convergence on outer iteration number ", iterations_outer
            end if
        end do  ! Inner

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
            phi_new_outer(c) = weighted_sum  ! 1/cm^2-s-MeV
        end do

        ! Relative error for the outer loop
        iterations_outer = iterations_outer + int(1, 8)
        err = maxval(dabs((phi_new_outer - phi_old_outer)) / phi_new_outer)
        if (err <= tolerance_outer) then
            print *, "Converged after ", iterations_outer, " outer iterations"
            cont_calc_outer = .false.
        else if (iterations_outer > num_iter_outer) then
            cont_calc_outer = .false.
            print *, "No convergence on outer loop; quit after ", iterations_outer, " iterations"
        end if
    end do  ! Outer

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
        do k = 1, num_materials
            ! Leakage in neg. direction from left face
            if (c == 1) then
                do m = 1, (num_ords / 2)
                    leakage_l = leakage_l + dabs(mu(m)) * weights(m) * prob(k) * psi_i_m(c, m, k)
                end do
            ! Leakage in pos. direction from right face
            else if (c == num_cells) then
                do m = (num_ords / 2 + 1), num_ords
                    leakage_r = leakage_r + dabs(mu(m)) * weights(m) * prob(k) * psi_i_p(c, m, k)
                end do
            end if
            ! Total absorption in system
            total_abs = total_abs + (macro_tot(c, k) - macro_scat(c, k)) * delta_x(c) * prob(k) * phi_new_inner(c, k)
            ! Distributed source
            source_dist = source_dist + spont_source(c, k) * delta_x(c) * prob(k)
        end do
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
    open(unit=7, file="./out/steady_state_slab_closure.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7,*) cell_vector(i), phi_new_outer(i)
    end do
    close(7)
end program steady_state_slab_closure
