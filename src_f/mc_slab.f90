program mc_slab
    use mcnp_random
    use self_library

    implicit none

    ! Constant parameters
    integer(8), parameter :: &
        num_particles = int(1.0d+8, 8)
    integer, parameter :: &
        num_points = int(2.0d+1, 4)

    ! Material properties
    real(8), parameter :: &
        thickness = 1.0d+1, &  ! cm
        struct_thickness = thickness / dble(num_points), &
        tot_const = 1.0d+0, &  ! 1/cm
        scat_const = 0.9d+0 * tot_const  ! 1/cm

    ! Material variables
    real(8), dimension(num_points) :: &
        delta_x, macro_scat, macro_tot

    ! Calculation variables
    integer :: &
        particle_counter, c
    real(8) :: &
        leakage_l, leakage_r, mu, psi_bound_l, psi_bound_r
    real(8), dimension(num_points) :: &
        phi

    ! Additional variables (plotting, etc.)
    real(8), dimension(num_points) :: &
        point_vector

    ! Assignment of material variables
    delta_x(:) = thickness / dble(num_points)  ! cm
    phi(:) = 0.0d+0  ! 1/cm^2-s-MeV
    macro_scat(:) = scat_const  ! 1/cm
    macro_tot(:) = tot_const  ! 1/cm

    ! Boundary conditions: 1/cm^2-s-MeV
    psi_bound_l = 1.0d+0  ! Isotropic
    psi_bound_r = 0.0d+0  ! Vacuum

    !$omp parallel do default(private) reduction(+:phi)
    do particle_counter = 1, num_particles
        do c = 1, num_points
            mu = rang()
            phi(c) = phi(c) + mu
        end do
    end do
    !$omp end parallel do

    do c = 1, num_points
        print *, phi(c) / dble(num_particles)
    end do

    ! Create plot
    call linspace(point_vector, 0.0d+0, thickness, num_points)
    open(unit=7, file="./out/mc_slab.out", form="formatted", status="replace", action="write")
    do c = 1, num_points
        write(7,*) point_vector(c), phi(c)
    end do
    close(7)
end program mc_slab
