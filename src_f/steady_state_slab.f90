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
    real(8), parameter :: &
        thickness = 1.0d+0, &  ! cm
        chord_a = 0.05d+0, &  ! cm
        chord_b = 0.05d+0, &  ! cm
        alpha_l = 0.0d+0, &
        alpha_r = 0.0d+0, &
        struct_thickness = thickness / dble(num_cells), &  ! cm
        spont_source_const = 0.0d+0  ! 1/cm^3
    ! For alpha (albedo boundary), 0.0 = no refl., 1.0 = total refl.

    real(8), dimension(2), parameter :: &
        scat_const = (/0.2d+0, 0.3d+0/), &  ! 1/cm
        fis_const = (/0.0d+0, 0.0d+0/), &  ! 1/cm
        tot_const = (/1.0d+0, 1.0d+0/)  ! 1/cm

    ! Material variables
    ! Allocated as (num_ind_cells) (by subroutine)
    real(8), dimension(:), allocatable :: &
        delta_x, x_points
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials
    ! Allocated as (num_groups, num_groups, num_cells, num_materials)
    real(8), dimension(:, :, :, :), allocatable :: &
        macro_scat
    ! Allocated as (num_groups, num_cells, num_materials)
    real(8), dimension(:, :, :), allocatable :: &
        scat_source, fis_source, spont_source, tot_source
    ! Allocated as (num_groups, num_cells, num_materials)
    real(8), dimension(:, :, :), allocatable :: &
        macro_fis, macro_tot
    real(8), dimension(num_groups) :: &
        chi, nu

    ! Calculation variables
    integer(8) :: &
        iterations
    integer :: &
        i, c, g, g_prime, m, k, num_ind_cells, material_num
    logical :: &
        cont_calc, point_found
    real(8) :: &
        tolerance, scatter_into, fission_into, weighted_sum, err, cons_distance
    real(8), dimension(num_groups, num_cells) :: &
        phi_new, phi_old
    real(8), dimension(num_groups, num_cells, num_materials) :: &
        phi_div
    ! Allocated as (num_groups, num_cells)
    real(8), dimension(:, :), allocatable :: &
        phi_morph_1, phi_morph_2, phi_morph
    ! Allocated as (num_groups, num_cells, num_ords, num_materials)
    real(8), dimension(:, :, :, :), allocatable :: &
        psi, psi_i_p, psi_i_m
    real(8), dimension(num_groups, num_ords, num_materials) :: &
        psi_bound_l, psi_bound_r, psi_refl_l, psi_refl_r
    real(8), dimension(num_ords) :: &
        ordinates, weights, mu

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector

    ! Assigment of material variables
    chi(:) = 1.0d+0 / dble(num_groups)  ! fission energy spectrum
    nu(:) = 2.4144d+0  ! average number of neutrons emitted per fission
    ! Assignment of initial calculation variables
    phi_new(:, :) = 1.0d+0  ! 1/cm^2-s-MeV, assume scalar flux is init. const.
    phi_old(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
    ! Material dependent fluxes
    phi_div(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV
    ! Can adjust boundaries by group and by ordinate
    psi_bound_l(:, :, :) = 1.0d+0  ! 1/cm^2-s-MeV-strad
    psi_bound_r(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

    ! Initial reflection values
    psi_refl_l(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
    psi_refl_r(:, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

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
        allocate(macro_scat(num_groups, num_groups, num_ind_cells, num_materials))
        allocate(macro_fis(num_groups, num_cells, num_materials))
        allocate(macro_tot(num_groups, num_cells, num_materials))
        do k = 1, num_materials
            macro_scat(:, :, :, k) = scat_const(k)  ! 1/cm
            macro_fis(:, :, k) = fis_const(k)  ! 1/cm
            macro_tot(:, :, k) = tot_const(k)  ! 1/cm
        end do

        ! Allocate and define calculational arrays
        allocate(psi(num_groups, num_ind_cells, num_ords, num_materials))
        psi(:, :, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_p(num_groups, num_ind_cells, num_ords, num_materials))
        psi_i_p(:, :, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad
        allocate(psi_i_m(num_groups, num_ind_cells, num_ords, num_materials))
        psi_i_m(:, :, :, :) = 0.0d+0  ! 1/cm^2-s-MeV-strad

        ! Allocate and define initial source terms
        allocate(scat_source(num_groups, num_ind_cells, num_materials))
        scat_source(:, :, :) = 0.0d+0  ! 1/cm^3-s
        allocate(fis_source(num_groups, num_ind_cells, num_materials))
        fis_source(:, :, :) = 0.0d+0  ! 1/cm^3-s
        allocate(spont_source(num_groups, num_ind_cells, num_materials))
        spont_source(:, :, :) = spont_source_const  ! 1/cm^3-s
        allocate(tot_source(num_groups, num_ind_cells, num_materials))
        tot_source(:, :, :) = 0.0d+0  ! 1/cm^3-s

        ! Allocate phi calculations
        allocate(phi_morph(num_groups, num_ind_cells))
        call struct_to_unstruct(phi_new, struct_thickness, num_cells, &
                                phi_morph, delta_x, num_ind_cells)
        allocate(phi_morph_1(num_groups, num_ind_cells))
        phi_morph_1(:, :) = 0.0d+0  ! 1/cm^2-s-MeV
        allocate(phi_morph_2(num_groups, num_ind_cells))
        phi_morph_2(:, :) = 0.0d+0  ! 1/cm^2-s-MeV

        ! Determine sources for each cell and group
        do c = 1, num_ind_cells
            do g = 1, num_groups
                do k = 1, num_materials
                    ! Scatter into, fission into, in-group scattering
                    scatter_into = 0.0d+0  ! neutrons
                    fission_into = 0.0d+0  ! neutrons
                    cons_distance = 0.0d+0  ! cm
                    ! For generating sources, go until the point is found (mapping)
                    point_found = .false.
                    do g_prime = 1, num_groups
                        scatter_into = scatter_into + macro_scat(g_prime, g, c, k) &
                                       * phi_morph(g_prime, c)
                        fission_into = fission_into + chi(g) * nu(g_prime) &
                                       * macro_fis(g_prime, c, k) * phi_morph(g_prime, c)
                    end do
                    scat_source(g, c, k) = scatter_into
                    fis_source(g, c, k) = fission_into
                    tot_source(g, c, k) = scat_source(g, c, k) + fis_source(g, c, k) &
                                       + spont_source(g, c, k)
                end do
            end do
        end do

        ! Forward sweep (left to right)
        ! First cell (left boundary)
        do g = 1, num_groups
            ! Ordinate loop, only consider the pos. ords for forward motion
            do m = (num_ords / 2 + 1), num_ords
                material_num = materials(1)
                ! Lewis and Miller Eq. 3-40
                psi(g, 1, m, material_num) = (1.0d+0 + (macro_tot(g, 1, material_num) * delta_x(1)) &
                                / (2.0d+0 * dabs(mu(m))))**(-1) &
                               * (psi_bound_l(g, m, material_num) + psi_refl_l(g, m, material_num) &
                                  + (tot_source(g, 1, material_num) * delta_x(1)) &
                                  / (2.d0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-41
                psi_i_p(g, 1, m, material_num) = 2.0d+0 * psi(g, 1, m, material_num) - psi_bound_l(g, m, material_num) &
                                   - psi_refl_l(g, m, material_num)
            end do
        end do
        ! Rest of the cells (sans left bounding cell)
        do c = 2, num_ind_cells
            do g = 1, num_groups
                do m = (num_ords / 2 + 1), num_ords
                    material_num = materials(c)
                    ! Continuity of boundaries
                    psi_i_m(g, c, m, material_num) = psi_i_p(g, c-1, m, material_num)
                    ! Lewis and Miller Eq. 3-40
                    psi(g, c, m, material_num) = (1.0d+0 + (macro_tot(g, c, material_num) * delta_x(c)) &
                                    / (2.0d+0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_m(g, c, m, material_num) + (tot_source(g, c, material_num) &
                                                          * delta_x(c)) &
                                      / (2.0d+0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-41
                    psi_i_p(g, c, m, material_num) = 2.0d+0 * psi(g, c, m, material_num) - psi_i_m(g, c, m, material_num)
                end do
            end do
        end do

        ! Backward sweep (right to left)
        ! First cell (right boundary)
        do g = 1, num_groups
            ! Ordinate loop, only consider neg. ords for backwards motion
            do m = 1, (num_ords / 2)
                material_num = materials(num_ind_cells)
                ! Lewis and Miller Eq. 3-42
                psi(g, num_ind_cells, m, material_num) = (1.0d+0 + (macro_tot(g, num_ind_cells, material_num) &
                                                * delta_x(num_ind_cells)) &
                                        / (2.0d+0 * dabs(mu(m))))**(-1) &
                                       * (psi_bound_r(g, m, material_num) &
                                          + psi_refl_r(g, m, material_num) &
                                          + (tot_source(g, num_ind_cells, material_num) &
                                             * delta_x(num_ind_cells)) &
                                          / (2.0d+0 * dabs(mu(m))))
                ! Lewis and Miller Eq. 3-43
                psi_i_m(g, num_ind_cells, m, material_num) = 2.0d+0 * psi(g, num_ind_cells, m, material_num) &
                                           - psi_bound_r(g, m, material_num) &
                                           - psi_refl_r(g, m, material_num)
            end do
        end do
        ! Rest of the cells (sans right bounding cell)
        do c = (num_ind_cells - 1), 1, -1
            do g = 1, num_groups
                do m = 1, (num_ords / 2)
                    material_num = materials(c)
                    ! Continuation of boundaries
                    psi_i_p(g, c, m, material_num) = psi_i_m(g, c+1, m, material_num)
                    ! Lewis and Miller Eq. 3-42
                    psi(g, c, m, material_num) = (1.0d+0 + (macro_tot(g, c, material_num) * delta_x(c)) &
                                    / (2.0d+0 * dabs(mu(m))))**(-1) &
                                   * (psi_i_p(g, c, m, material_num) + (tot_source(g, c, material_num) &
                                                          * delta_x(c)) &
                                      / (2.0d+0 * dabs(mu(m))))
                    ! Lewis and Miller Eq. 3-43
                    psi_i_m(g, c, m, material_num) = 2.0d+0 * psi(g, c, m, material_num) - psi_i_p(g, c, m, material_num)
                end do
            end do
        end do

        ! Reflected angular fluxes at the left and right boundaries
        do g = 1, num_groups
            do m = 1, num_ords
                if (m >= num_ords / 2 + 1) then
                    ! Last cell, pos. angular flux reflected to neg. direction
                    psi_refl_r(g, num_ords - m + 1, material_num) = alpha_r &
                                                  * psi_i_p(g, num_cells, m, material_num)
                else
                    ! First cell, neg. angular flux reflected to pos. direction
                    psi_refl_l(g, num_ords - m + 1, material_num) = alpha_l &
                                                  * psi_i_m(g, 1, m, material_num)
                end if
            end do
        end do

        ! Calculate phi from psi
        ! Lewis and Miller Eq. 3-5
        do c = 1, num_ind_cells
            do g = 1, num_groups
                weighted_sum = 0.0d+0
                do m = 1, num_ords
                    weighted_sum = weighted_sum + weights(m) * psi(g, c, m, material_num)
                end do
                phi_morph(g, c) = 0.50d+0 * weighted_sum  ! 1/cm^2-s-MeV
            end do
        end do

        ! Map the unstructured phi onto the structured phi
        call unstruct_to_struct(phi_morph, delta_x, phi_new, struct_thickness, num_cells)

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

        ! Deallocations of variable-width unstructured arrays
        deallocate(macro_scat)
        deallocate(macro_fis)
        deallocate(macro_tot)
        deallocate(psi)
        deallocate(psi_i_p)
        deallocate(psi_i_m)
        deallocate(scat_source)
        deallocate(fis_source)
        deallocate(spont_source)
        deallocate(tot_source)
        deallocate(phi_morph)
        deallocate(phi_morph_1)
        deallocate(phi_morph_2)
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
