program steady_state_slab
    use self_library
    use geometry_gen
    use mesh_map
    use mcnp_random

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_iter_outer = int(20000, 8), &
        num_iter_inner = int(5000, 8)
    integer, parameter :: &
        num_cells = int(100, 4), &
        num_ords = int(16, 4), &
        num_materials = int(2, 4)

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+0, &  ! cm
        struct_thickness = thickness / dble(num_cells)  ! cm
    ! For alpha (albedo boundary), 0.0 = no refl., 1.0 = total refl.

    real(8), dimension(num_materials), parameter :: &
        scat_const = (/0.2d+0, 0.3d+0/), &  ! 1/cm
        tot_const = (/1.0d+0, 1.0d+0/), &  ! 1/cm
        spont_source_const = (/0.0d+0, 0.0d+0/),&  ! 1/cm^3
        chord = (/0.05d+0, 0.05d+0/)  ! cm

    ! Material variables
    ! Allocated as (num_ind_cells) (by subroutine)
    real(8), dimension(:), allocatable :: &
        delta_x, x_points
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials
    ! Allocated as (num_ind_cells, num_materials)
    real(8), dimension(:, :), allocatable :: &
        macro_scat, macro_tot, scat_source, spont_source, tot_source, phi_div

    ! Calculation variables
    integer(8) :: &
        iterations_inner, iterations_outer
    integer :: &
        i, c, m, k, num_ind_cells, material_num
    logical :: &
        cont_calc_outer, cont_calc_inner
    real(8) :: &
        tolerance, weighted_sum, err_inner, err_outer, cons_distance
    real(8), dimension(num_cells) :: &
        phi_new_outer, phi_old_outer
    real(8), dimension(num_cells, num_materials) :: &
        psi_mat, phi_mat
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
    phi_mat(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    ! Assignment of initial calculation variables
    phi_new_outer(:) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old_outer(:) = 0.0d+0  ! 1/cm^2-s-MeV
    ! Material dependent fluxes
    phi_div(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
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

    ! Tolerance for ending each instance
    tolerance = 1.0e-8

    ! Calculation: iterations
    cont_calc_outer = .true.
    ! Start counter at zero
    iterations_outer = int(0, 8)

    ! Outer loop over overall mixed system
    do while (cont_calc_outer)
        phi_old_outer = phi_new_outer  ! 1/cm^2-s-MeV

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
        allocate(macro_scat(num_ind_cells, num_materials))
        allocate(macro_tot(num_ind_cells, num_materials))
        do k = 1, num_materials
            macro_scat(:, k) = scat_const(k)  ! 1/cm
            macro_tot(:, k) = tot_const(k)  ! 1/cm
        end do

        ! Allocate and define calculational arrays
        allocate(psi(num_ind_cells, num_ords))
        psi(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_p(num_ind_cells, num_ords))
        psi_i_p(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_m(num_ind_cells, num_ords))
        psi_i_m(:, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

        ! Allocate and define initial source terms
        allocate(scat_source(num_ind_cells, num_materials))
        scat_source(:, :) = 0.0d+0  ! 1/cm^3-s
        allocate(spont_source(num_ind_cells, num_materials))
        do k = 1, num_materials
            spont_source(:, k) = spont_source_const(k)  ! 1/cm^3-s
        end do
        allocate(tot_source(num_ind_cells, num_materials))
        tot_source(:, :) = 0.0d+0  ! 1/cm^3-s

        ! Allocate phi calculations
        allocate(phi_morph_new(num_ind_cells))
        phi_morph_new(:) = 1.0d+0
!        call struct_to_unstruct(phi_new_outer, struct_thickness, num_cells, &
!                                phi_morph_new, delta_x, num_ind_cells)
        allocate(phi_morph_old(num_ind_cells))
        phi_morph_old(:) = 0.0d+0  ! 1/cm^2-s-MeV

        ! Inner loop over generated geometry
        do while (cont_calc_inner)
            phi_morph_old = phi_morph_new  ! 1/cm^2-s-MeV
            ! Determine sources for each cell and group
            do c = 1, num_ind_cells
                do k = 1, num_materials
                    ! Scatter into, fission into, in-group scattering
                    scat_source(c, k) = macro_scat(c, k) / 2.0d+0 * phi_morph_new(c)
                    tot_source(c, k) = scat_source(c, k) + spont_source(c, k) / 2.0d+0
                end do
            end do

            ! Forward sweep (left to right)
            ! First cell (left boundary)
            ! Ordinate loop, only consider the pos. ords for forward motion
            do m = (num_ords / 2 + 1), num_ords
                material_num = materials(1)
                psi(1, m) = (1.0d+0 + (macro_tot(1, material_num) * delta_x(1)) &
                             / (2.0d+0 * dabs(mu(m))))**(-1) &
                            * (psi_bound_l(m) &
                               + (tot_source(1, material_num) * delta_x(1)) &
                               / (2.0d+0 * dabs(mu(m))))
                psi_i_p(1, m) = 2.0d+0 * psi(1, m) - psi_bound_l(m)
            end do
            ! Rest of the cells (sans left bounding cell)
            do c = 2, num_ind_cells
                do m = (num_ords / 2 + 1), num_ords
                    ! Continuity of boundaries
                    psi_i_m(c, m) = psi_i_p(c - 1, m)
                    psi(c, m) = (1.0d+0 + (macro_tot(c, material_num) * delta_x(c)) &
                                 / (2.0d+0 * dabs(mu(m))))**(-1) &
                                * (psi_i_m(c, m) + (tot_source(c, material_num) &
                                                    * delta_x(c)) / (2.0d+0 * dabs(mu(m))))
                    psi_i_p(c, m) = 2.0d+0 * psi(c, m) - psi_i_m(c, m)
                end do
            end do

            ! Backward sweep (right to left)
            ! First cell (right boundary)
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                material_num = materials(num_ind_cells)
                psi(num_ind_cells, m) = (1.0d+0 + (macro_tot(num_ind_cells, material_num) &
                                                   * delta_x(num_ind_cells)) &
                                         / (2.0d+0 * dabs(mu(m))))**(-1) &
                                        * (psi_bound_r(m) &
                                           + (tot_source(num_ind_cells, material_num) &
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
                    psi(c, m) = (1.0d+0 + (macro_tot(c, material_num) * delta_x(c)) &
                                 / (2.0d+0 * dabs(mu(m))))**(-1) &
                                * (psi_i_p(c, m) + (tot_source(c, material_num) &
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
                phi_morph_new(c) = 0.50d+0 * weighted_sum  ! 1/cm^2-s-MeV
            end do

            ! Relative error for inner loop
            iterations_inner = iterations_inner + int(1, 8)
            err_inner = maxval(dabs((phi_morph_old - phi_morph_new)) / phi_morph_new)
            if (err_inner <= tolerance) then
                cont_calc_inner = .false.
            else if (iterations_inner > num_iter_inner) then
                cont_calc_inner = .false.
                print *, "No convergence on inner loop number ", iterations_outer, &
                         ": quit after maximum ", iterations_inner, " iterations"
            end if
        end do  ! Inner loop

        ! Obtain material balances
        !do k = 1, num_materials
        !    do c = 1, num_cells
        !        do m = 1, num_ords
        !            if (phi_mat(c, k) == 0.0d+0) then
        !                phi_mat(c, k) = 
        !            end if
        !        end do
        !    end do
        !end do

        ! Map the unstructured phi onto the structured phi
        call unstruct_to_struct(phi_morph_new, delta_x, phi_new_outer, struct_thickness, num_cells)
        ! NOTE: Compare these two, need to average in some way

        ! Relative error for outer loop
        iterations_outer = iterations_outer + int(1, 8)
        err_outer = maxval(dabs((phi_new_outer - phi_old_outer)) / phi_new_outer)
        if (err_outer <= tolerance) then
            cont_calc_outer = .false.
            print *, "Quit after ", iterations_outer , " outer iterations"
        else if (iterations_outer > num_iter_outer) then
            cont_calc_outer = .false.
            print *, "No convergence on outer loop: quit after maximum ", iterations_outer, " outer iterations"
        end if

        ! Deallocations of variable-width unstructured arrays
        deallocate(macro_scat)
        deallocate(macro_tot)
        deallocate(psi)
        deallocate(psi_i_p)
        deallocate(psi_i_m)
        deallocate(scat_source)
        deallocate(spont_source)
        deallocate(tot_source)
        deallocate(phi_morph_old)
        deallocate(phi_morph_new)
    end do  ! Outer loop

    ! Create plot
    call linspace(cell_vector, 0.0d+0, thickness, num_cells)
    open(unit=7, file="./out/steady_state_slab.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7,*) cell_vector(i), phi_new_outer(i)
    end do
    close(7)
end program steady_state_slab
