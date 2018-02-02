program geometry_test
    !use dynamic_arrays
    !use mcnp_random
    use geometry_gen
    implicit none

    double precision, dimension(:), allocatable :: &
        x_dist
    integer, dimension(:), allocatable :: &
        materials
    double precision :: &
        thickness, chord_a, chord_b, cons_thickness, rand_num, distance
    integer :: &
        num_cells

    thickness = 10.d+0  ! cm
    cons_thickness = 0.d+0  ! cm
    chord_a = 0.05d+0  ! cm
    chord_b = 0.05d+0  ! cm
    num_cells = 0

    allocate(x_dist(0))
    allocate(materials(0))

    call get_geometry(x_dist, materials, chord_a, chord_b, thickness, num_cells)

    !call RN_init_problem(int(1234567, 8), 1)
    !do while (cons_thickness < thickness)
    !    rand_num = rang()
    !    distance = chord_a * dlog(1.0d+0 / (1.0d+0 - rand_num))  ! cm
    !    cons_thickness = cons_thickness + distance  ! cm
    !    num_cells = num_cells + 1
    !    call push_double(x_dist, distance)
    !end do

    print *, x_dist
    print *, materials

    deallocate(x_dist)

    !print *, cons_thickness
    !print *, num_cells
end program geometry_test
