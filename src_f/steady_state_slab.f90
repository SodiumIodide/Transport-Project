program steady_state_slab
    use self_library
    use geometry_gen
    use mesh_map
    use mcnp_random

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter_outer = int(100000, 8), &
        num_iter_inner = int(1000, 8)
    integer, parameter :: &
        num_cells = int(100, 4), &
        num_ords = int(16, 4), &
        num_materials = int(2, 4)
    logical, parameter :: &
        tally = .false.

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+0, &  ! cm
        struct_thickness = thickness / dble(num_cells), &  ! cm
        inner_tolerance = epsilon(1.0d+0), &
        outer_tolerance = 1.0d-5
    ! For alpha (albedo boundary), 0.0 = no refl., 1.0 = total refl.

    real(8), dimension(num_materials), parameter :: &
        scat_const = (/0.2d+0, 0.3d+0/), &  ! 1/cm
        tot_const = (/1.0d+0, 1.0d+0/), &  ! 1/cm
        spont_source_const = (/0.0d+0, 0.0d+0/),&  ! 1/cm^3
        chord = (/0.05d+0, 0.05d+0/)  ! cm

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
        i, c, m, k, num_ind_cells, material_num
    logical :: &
        cont_calc_outer, cont_calc_inner
    real(8) :: &
        weighted_sum, err_inner, err_outer, total_chord, &
        leakage_l, leakage_r, total_abs, balance_source_l, balance_source_r, source_dist, balance
    real(8), dimension(num_cells) :: &
        phi_outer, phi_real
    real(8), dimension(num_cells, num_materials) :: &
        phi_mat_old, phi_mat_new
    ! Allocated as (num_ind_cells)
    real(8), dimension(:), allocatable :: &
        phi_morph_new, phi_morph_old
    ! Allocated as (num_ind_cells, num_ords)
    real(8), dimension(:, :), allocatable :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu, psi_bound_l, psi_bound_r

    ! Allocation variables
    !real(8), dimension(:, :), allocatable :: &
    !    two_d_arr
    !real(8), dimension(:), allocatable :: &
    !    one_d_arr

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector

    ! Assigment of material variables
    macro_tot = tot_const  ! 1/cm
    macro_scat = scat_const  ! 1/cm
    phi_real(:) = 0.0d+0  ! 1/cm^2-s-MeV
    phi_mat_new(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    phi_mat_old(:, :) = 1.0d+0  ! 1/cm^2-s-MeV
    ! Assignment of initial calculation variables
    phi_outer(:) = 1.0d+0  ! 1/cm^2-s-MeV
    do k = 1, num_materials
        total_chord = total_chord + chord(k)
    end do
    prob = chord / total_chord
    ! Material dependent fluxes
    ! Can adjust boundaries by group and by ordinate

    ! Initial value (assigned via get_geometry subroutine)
    num_ind_cells = 0

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
    call RN_init_problem(int(12345678, 8), 1)

    ! Calculation: iterations
    cont_calc_outer = .true.
    ! Start counter at zero
    iterations_outer = int(0, 8)

    ! Outer loop over overall mixed system
    do while (cont_calc_outer)
        phi_mat_old = phi_mat_new  ! 1/cm^2-s-MeV

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
!        call struct_to_unstruct(phi_new_outer, struct_thickness, num_cells, &
!                                phi_morph_new, delta_x, num_ind_cells)
        if (allocated(phi_morph_old)) then
            deallocate(phi_morph_old)
        end if
        allocate(phi_morph_old(num_ind_cells))
        phi_morph_old(:) = 0.0d+0  ! 1/cm^2-s-MeV

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
                material_num = materials(1)
                psi(1, m) = (1.0d+0 + (macro_tot(material_num) * delta_x(1)) &
                             / (2.0d+0 * dabs(mu(m))))**(-1) &
                            * (psi_bound_l(m) &
                               + (tot_source(1) * delta_x(1)) &
                               / (2.0d+0 * dabs(mu(m))))
                psi_i_p(1, m) = 2.0d+0 * psi(1, m) - psi_bound_l(m)
            end do
            ! Rest of the cells (sans left bounding cell)
            do c = 2, num_ind_cells
                do m = (num_ords / 2 + 1), num_ords
                    material_num = materials(c)
                    ! Continuity of boundaries
                    psi_i_m(c, m) = psi_i_p(c - 1, m)
                    psi(c, m) = (1.0d+0 + (macro_tot(material_num) * delta_x(c)) &
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
                material_num = materials(num_ind_cells)
                psi(num_ind_cells, m) = (1.0d+0 + (macro_tot(material_num) &
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
                    material_num = materials(c)
                    ! Continuation of boundaries
                    psi_i_p(c, m) = psi_i_m(c + 1, m)
                    psi(c, m) = (1.0d+0 + (macro_tot(material_num) * delta_x(c)) &
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
            else if (iterations_inner > num_iter_inner) then
                cont_calc_inner = .false.
                print *, "No convergence on inner loop number ", iterations_outer, &
                         ": quit after maximum ", iterations_inner, " iterations"
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

        iterations_outer = iterations_outer + int(1, 8)

        ! Obtain material balances
        call material_calc(phi_morph_new, delta_x, num_ind_cells, materials, phi_mat_new, &
                           struct_thickness, num_cells, num_materials, iterations_outer)
        call material_boundaries()

        ! Map the unstructured phi onto the structured phi
        !call unstruct_to_struct(phi_morph_new, delta_x, phi_real, struct_thickness, num_cells)

        ! Relative error for outer loop
        err_outer = maxval(dabs((phi_mat_new - phi_mat_old)) / phi_mat_new)
        if (err_outer <= outer_tolerance) then
            cont_calc_outer = .false.
            print *, "Quit after ", iterations_outer , " outer iterations"
        else if (iterations_outer > num_iter_outer) then
            cont_calc_outer = .false.
            print *, "No convergence on outer loop: quit after maximum ", iterations_outer, " outer iterations"
        end if
    end do  ! Outer loop

    call unstruct_to_struct(phi_morph_new, delta_x, num_ind_cells, phi_real, struct_thickness, num_cells)

    ! Check balances outer loop
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

    ! Create plot
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)
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
    deallocate(psi)
    deallocate(psi_i_p)
    deallocate(psi_i_m)
    deallocate(scat_source)
    deallocate(spont_source)
    deallocate(tot_source)
    deallocate(phi_morph_old)
    deallocate(phi_morph_new)
end program steady_state_slab
