program mc_slab
    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+4, 8)
    integer, parameter :: &
        num_cells = int(2.0d+2, 4), &

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_cells), &
        tot_const = 1.0d+0, &  ! 1/cm
        scat_const = 0.9d+0 * tot_const  ! 1/cm

    ! Material variables
    ! Allocated as (num_ind_cells) (by subroutine)
    real(8), dimension(:), allocatable :: &
        delta_x, x_points, macro_scat, macro_tot
    ! Allocated as (num_ind_cells) (by subroutine)
    integer, dimension(:), allocatable :: &
        materials
    real(8), dimension(num_materials) :: &
        prob

    ! Calculation variables
    integer :: &
        num_ind_cells, k
    real(8) :: &
        total_chord, leakage_l, leakage_r, mu, psi_bound_l, psi_bound_r
    real(8), dimension(num_cells) :: &
        phi_real
    real(8), dimension(num_cells, num_materials) :: &
        phi_mat
    ! Allocated as (num_ind_cells)
    real(8), dimension(:), allocatable :: &
        phi_morph

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_cells) :: &
        cell_vector
    character(len=1024) :: &
        filename

    ! Assignment of material variables
    total_chord = 0.0d+0
    do k = 1, num_materials
        total_chord = total_chord + chord(k)
    end do
    prob = chord / total_chord
end program mc_slab
