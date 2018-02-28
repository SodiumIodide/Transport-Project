program steady_state_slab
    use self_library
    use geometry_gen
    use mesh_map
    use mcnp_random

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter_outer = int(1.0d+7, 8), &
        num_iter_inner = int(1.0d+6, 8), &
        seed = int(123456, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4), &
        num_ords = int(16, 4), &
        num_materials = int(2, 4)
    logical, parameter :: &
        tally = .false.

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_cells), &  ! cm
        inner_tolerance = 1.0d-7

    real(8), dimension(num_materials), parameter :: &
        tot_const = (/dble(10)/dble(99), dble(100)/dble(11)/), &  ! 1/cm
        scat_const = (/dble(10)/dble(99)*0.9d+0, dble(100)/dble(11)*0.9d+0/), &  ! 1/cm
        chord = (/dble(99)/dble(10), dble(11)/dble(10)/), &  ! cm
        spont_source_const = (/0.0d+0, 0.0d+0/)  ! 1/cm^3

    ! Material variables
    ! Allocated as (num_ind_cells) (by subroutine)
    real(8), dimension(:), allocatable :: &
        delta_x, x_points, scat_source, spont_source, tot_source
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials
    real(8), dimension(num_materials) :: &
        macro_scat, macro_tot, prob

    ! Calculation variables
    integer(8) :: &
        iterations_inner, iterations_outer
    integer :: &
        i, c, m, k, num_ind_cells
    logical :: &
        cont_calc_outer, cont_calc_inner, left_switch, right_switch, converged
    real(8) :: &
        weighted_sum, err_inner, total_chord, &
        leakage_l, leakage_r, total_abs, balance_source_l, balance_source_r, &
        source_dist, balance
    real(8), dimension(num_cells) :: &
        phi_outer, phi_real
    real(8), dimension(num_cells, num_materials) :: &
        phi_mat_old, phi_mat_new
    real(8), dimension(num_ords, num_materials) :: &
        psi_mat_leak_l, psi_mat_leak_r
    ! Allocated as (num_ind_cells)
    real(8), dimension(:), allocatable :: &
        phi_morph_new, phi_morph_old
    ! Allocated as (num_ind_cells, num_ords)
    real(8), dimension(:, :), allocatable :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu, psi_bound_l, psi_bound_r

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector
    character(len=1024) :: &
        filename

    ! Assigment of material variables
    macro_tot = tot_const  ! 1/cm
    macro_scat = scat_const  ! 1/cm
    phi_mat_new(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    phi_mat_old(:, :) = 1.0d+0  ! 1/cm^2-s-MeV
    ! Assignment of initial calculation variables
    phi_outer(:) = 1.0d+0  ! 1/cm^2-s-MeV
    do k = 1, num_materials
        total_chord = total_chord + chord(k)
    end do
    prob = chord / total_chord
    psi_mat_leak_l(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    psi_mat_leak_r(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    ! Tally averages for boundary leakage
    left_switch = .false.
    right_switch = .false.

    ! Initial value (assigned via get_geometry subroutine)
    num_ind_cells = 0

    ! X values of flux for plotting
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)

    ! Legendre Gauss Quadrature values over chosen ordinates
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
    !psi_bound_r(:) = 1.0d+0  ! Isotropic source
    psi_bound_r(:) = 0.0d+0  ! Vacuum
    ! Beam source (in conj. with vacuum):
    !psi_bound_r(1) = 1.0d+0 / (mu(1) * weights(1))

    ! Start random generation seed
    call RN_init_problem(seed, 1)

    ! Calculation: iterations
    cont_calc_outer = .true.
    ! Start counter at zero
    iterations_outer = int(0, 8)

    ! Outer loop over overall mixed system
    do while (cont_calc_outer)
        cont_calc_inner = .true.
        ! Start inner iterations at zero
        iterations_inner = int(0, 8)

        ! Fill Markovian geometry
        ! Allocations and deallocations of delta_x and materials should be
        ! Handled by the subroutine
        call get_geometry(delta_x, x_points, materials, chord(1), chord(2), &
                          thickness, num_ind_cells)

        ! Allocate and define properties
        ! Anything that depends on the number of individual cells must be reallocated

        ! Allocate and define calculational arrays
        if (allocated(psi)) then
            deallocate(psi)
        end if
        allocate(psi(num_ind_cells, num_ords))
        psi(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        if (allocated(psi_i_p)) then
            deallocate(psi_i_p)
        end if
        allocate(psi_i_p(num_ind_cells, num_ords))
        psi_i_p(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        if (allocated(psi_i_m)) then
            deallocate(psi_i_m)
        end if
        allocate(psi_i_m(num_ind_cells, num_ords))
        psi_i_m(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

        ! Allocate and define initial source terms
        if (allocated(scat_source)) then
            deallocate(scat_source)
        end if
        allocate(scat_source(num_ind_cells))
        scat_source(:) = 0.0d+0  ! 1/cm^3-s
        if (allocated(spont_source)) then
            deallocate(spont_source)
        end if
        allocate(spont_source(num_ind_cells))
        spont_source(:) = 0.0d+0  ! 1/cm^3-s
        if (allocated(tot_source)) then
            deallocate(tot_source)
        end if
        allocate(tot_source(num_ind_cells))
        tot_source(:) = 0.0d+0  ! 1/cm^3-s

        ! Allocate phi calculations
        if (allocated(phi_morph_new)) then
            deallocate(phi_morph_new)
        end if
        allocate(phi_morph_new(num_ind_cells))
        phi_morph_new(:) = 1.0d+0
        if (allocated(phi_morph_old)) then
            deallocate(phi_morph_old)
        end if
        allocate(phi_morph_old(num_ind_cells))
        phi_morph_old(:) = 0.0d+0  ! 1/cm^2-s-MeV

        ! Flag to tell whether to use the inner loop data in averaging or not
        converged = .false.
        ! Inner loop over generated geometry
        do while (cont_calc_inner)
            phi_morph_old = phi_morph_new  ! 1/cm^2-s-MeV
            ! Determine sources for each cell and group
            do c = 1, num_ind_cells
                scat_source(c) = macro_scat(materials(c)) / 2.0d+0 * phi_morph_new(c)
                spont_source(c) = spont_source_const(materials(c)) / 2.0d+0
                tot_source(c) = scat_source(c) + spont_source(c)
            end do

            ! Forward sweep (left to right)
            ! First cell (left boundary)
            ! Ordinate loop, only consider the pos. ords for forward motion
            do m = (num_ords / 2 + 1), num_ords
                psi(1, m) = (1.0d+0 + (macro_tot(materials(1)) * delta_x(1)) &
                             / (2.0d+0 * dabs(mu(m))))**(-1) &
                            * (psi_bound_l(m) &
                               + (tot_source(1) * delta_x(1)) &
                               / (2.0d+0 * dabs(mu(m))))
                psi_i_p(1, m) = 2.0d+0 * psi(1, m) - psi_bound_l(m)
            end do
            ! Rest of the cells (sans left bounding cell)
            do c = 2, num_ind_cells
                do m = (num_ords / 2 + 1), num_ords
                    ! Continuity of boundaries
                    psi_i_m(c, m) = psi_i_p(c - 1, m)
                    psi(c, m) = (1.0d+0 + (macro_tot(materials(c)) * delta_x(c)) &
                                 / (2.0d+0 * dabs(mu(m))))**(-1) &
                                * (psi_i_m(c, m) + (tot_source(c) &
                                                    * delta_x(c)) / (2.0d+0 * dabs(mu(m))))
                    psi_i_p(c, m) = 2.0d+0 * psi(c, m) - psi_i_m(c, m)
                end do
            end do

            ! Backward sweep (right to left)
            ! First cell (right boundary)
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                psi(num_ind_cells, m) = (1.0d+0 + (macro_tot(materials(num_ind_cells)) &
                                                   * delta_x(num_ind_cells)) &
                                         / (2.0d+0 * dabs(mu(m))))**(-1) &
                                        * (psi_bound_r(m) &
                                           + (tot_source(num_ind_cells) &
                                              * delta_x(num_ind_cells)) &
                                           / (2.0d+0 * dabs(mu(m))))
                psi_i_m(num_ind_cells, m) = 2.0d+0 * psi(num_ind_cells, m) - psi_bound_r(m)
            end do
            ! Rest of the cells (sans right bounding cell)
            do c = (num_ind_cells - 1), 1, -1
                do m = 1, (num_ords / 2)
                    ! Continuation of boundaries
                    psi_i_p(c, m) = psi_i_m(c + 1, m)
                    psi(c, m) = (1.0d+0 + (macro_tot(materials(c)) * delta_x(c)) &
                                 / (2.0d+0 * dabs(mu(m))))**(-1) &
                                * (psi_i_p(c, m) + (tot_source(c) &
                                                    * delta_x(c)) / (2.0d+0 * dabs(mu(m))))
                    psi_i_m(c, m) = 2.0d+0 * psi(c, m) - psi_i_p(c, m)
                end do
            end do

            ! Calculate phi from psi
            do c = 1, num_ind_cells
                weighted_sum = 0.0d+0
                do m = 1, num_ords
                    weighted_sum = weighted_sum + weights(m) * psi(c, m)
                end do
                phi_morph_new(c) = weighted_sum  ! 1/cm^2-s-MeV
            end do

            ! Relative error for inner loop
            iterations_inner = iterations_inner + int(1, 8)
            err_inner = maxval(dabs((phi_morph_old - phi_morph_new)) / phi_morph_new)
            if (err_inner <= inner_tolerance) then
                cont_calc_inner = .false.
                converged = .true.
            else if (iterations_inner > num_iter_inner) then
                cont_calc_inner = .false.
                print *, "No convergence on inner loop number ", iterations_outer, &
                         ": quit after maximum ", iterations_inner, &
                         " iterations; data will not be used"
            end if
        end do  ! Inner loop

        ! Print balance checks for each realization
        if (tally) then
            ! Check balances inner loop
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
            do c = 1, num_ind_cells
                ! Leakage in neg. direction from left face
                if (c == 1) then
                    do m = 1, (num_ords / 2)
                        leakage_l = leakage_l + dabs(mu(m)) * weights(m) * psi_i_m(c, m)
                    end do
                ! Leakage in pos. direction from right face
                else if (c == num_ind_cells) then
                    do m = (num_ords / 2 + 1), num_ords
                        leakage_r = leakage_r + dabs(mu(m)) * weights(m) * psi_i_p(c, m)
                    end do
                end if
                ! Total absorption in system
                total_abs = total_abs + (macro_tot(materials(c)) - macro_scat(materials(c))) * delta_x(c) * phi_morph_new(c)
                ! Distributed source
                source_dist = source_dist + spont_source(materials(c)) * delta_x(c)
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
        end if

        if (converged) then
            iterations_outer = iterations_outer + int(1, 8)

            ! Obtain material balances
            ! Average flux in structured cells
            call material_calc(phi_morph_new, delta_x, num_ind_cells, materials, phi_mat_new, &
                            struct_thickness, num_cells, num_materials, iterations_outer)
            ! Average leakage from material boundaries (for balance/tally purposes)
            do k = 1, num_materials
                ! Adjust the left switch
                if (materials(1) == k) then
                    left_switch = .true.
                else
                    left_switch = .false.
                end if
                ! Adjust the right switch
                if (materials(num_ind_cells) == k) then
                    right_switch = .true.
                else
                    right_switch = .false.
                end if
                ! Calculate the average negative leakage from left boundary
                if (left_switch) then
                    do m = 1, (num_ords / 2)
                        ! Don't weight the initial zero
                        if (psi_mat_leak_l(m, k) == 0.0d+0) then
                            psi_mat_leak_l(m, k) = psi_i_m(1, m)
                        else
                            psi_mat_leak_l(m, k) = psi_mat_leak_l(m, k) &
                                + (psi_i_m(1, m) - psi_mat_leak_l(m, k)) / dble(iterations_outer)
                        end if
                    end do
                end if
                ! Calculate the average positive leakage from right boundary
                if (right_switch) then
                    do m = (num_ords / 2 + 1), num_ords
                        ! Don't weight the initial zero
                        if (psi_mat_leak_r(m, k) == 0.0d+0) then
                            psi_mat_leak_r(m, k) = psi_i_p(num_ind_cells, m)
                        else
                            psi_mat_leak_r(m, k) = psi_mat_leak_r(m, k) &
                                + (psi_i_p(num_ind_cells, m) - psi_mat_leak_r(m, k)) / dble(iterations_outer)
                        end if
                    end do
                end if
            end do

            if (mod(iterations_outer, 1000) == 0) then
                print *, "Realization number ", iterations_outer
            end if

        ! Save the first few realizations
            if (iterations_outer < 4) then
                phi_real(:) = 0.0d+0  ! 1/cm^2-s-MeV
                ! Map the final realization onto phi_real
                call unstruct_to_struct(phi_morph_new, delta_x, num_ind_cells, phi_real, struct_thickness, num_cells)
                ! Write the realization to an output file
                write(filename, "(A18,I1,A4)") "./out/realization_", iterations_outer, ".out"
                open(unit=7, file=filename, form="formatted", status="replace", action="write")
                do i = 1, num_cells
                    write(7,*) cell_vector(i), phi_real(i)
                end do
                close(7)
            end if

            if (iterations_outer > num_iter_outer) then
                cont_calc_outer = .false.
            end if
        end if  ! Logical test for inner convergence: don't average nonconverged samples
    end do  ! Outer loop

    ! Print the final realization onto phi_real for plotting purposes
    phi_real(:) = 0.0d+0  ! 1/cm^2-s-MeV
    call unstruct_to_struct(phi_morph_new, delta_x, num_ind_cells, phi_real, struct_thickness, num_cells)

    ! Normalize the resulting arrays
    !do k = 1, num_materials
    !    do c = 1, num_cells
    !        phi_mat_new(c, k) = phi_mat_new(c, k) / dble(num_iter_outer)
    !    end do
    !end do

    ! Calculate outer loop reflection and transmission
    leakage_l = 0.0d+0
    leakage_r = 0.0d+0
    ! Tally the overall losses due to refleciton and transmission
    do k = 1, num_materials
        ! Leakage in neg. direction from left face
        do m = 1, (num_ords / 2)
            leakage_l = leakage_l + dabs(mu(m)) * weights(m) * prob(k) * psi_mat_leak_l(m, k)
        end do
        ! Leakage in pos. direction from right face
        do m = (num_ords / 2 + 1), num_ords
            leakage_r = leakage_r + dabs(mu(m)) * weights(m) * prob(k) * psi_mat_leak_r(m, k)
        end do
    end do
    print *, "Reflection on left: ", leakage_l
    print *, "Transmission on right: ", leakage_r

    ! Create plot
    open(unit=7, file="./out/steady_state_slab.out", form="formatted", &
         status="replace", action="write")
    open(unit=8, file="./out/steady_state_slab_1.out", form="formatted", &
         status="replace", action="write")
    open(unit=9, file="./out/steady_state_slab_2.out", form="formatted", &
         status="replace", action="write")
    open(unit=10, file="./out/last_realization.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        phi_outer(i) = prob(1) * phi_mat_new(i, 1) + prob(2) * phi_mat_new(i, 2)
        write(7,*) cell_vector(i), phi_outer(i)
        write(8,*) cell_vector(i), phi_mat_new(i, 1)
        write(9,*) cell_vector(i), phi_mat_new(i, 2)
        write(10,*) cell_vector(i), phi_real(i)
    end do

    ! Close all open file streams
    close(7)
    close(8)
    close(9)
    close(10)

    ! Deallocate all variable-width unstructured arrays
    if (allocated(psi)) then
        deallocate(psi)
    end if
    if (allocated(psi_i_p)) then
        deallocate(psi_i_p)
    end if
    if (allocated(psi_i_m)) then
        deallocate(psi_i_m)
    end if
    if (allocated(scat_source)) then
        deallocate(scat_source)
    end if
    if (allocated(spont_source)) then
        deallocate(spont_source)
    end if
    if (allocated(tot_source)) then
        deallocate(tot_source)
    end if
    if (allocated(phi_morph_old)) then
        deallocate(phi_morph_old)
    end if
    if (allocated(phi_morph_new)) then
        deallocate(phi_morph_new)
    end if
end program steady_state_slab
