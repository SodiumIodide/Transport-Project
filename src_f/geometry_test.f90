program geometry_test
    use dynamic_arrays
    use mcnp_random
    implicit none

    double precision, dimension(:), allocatable :: &
        x_dist
    double precision :: &
        thickness, chord_a, chord_b, cons_thickness, rand_num, distance
    integer :: &
        num_cells

    thickness = 10.d0  ! cm
    cons_thickness = 0.d0  ! cm
    chord_a = 0.05d0  ! cm
    chord_b = 0.05d0  ! cm
    num_cells = 0

    allocate(x_dist(0))

    call RN_init_problem(int(1234567, 8), 1)
    do while (cons_thickness < thickness)
        rand_num = rang()
        distance = chord_a * dlog(1.d0 / (1.d0 - rand_num))
        cons_thickness = cons_thickness + distance
        num_cells = num_cells + 1
        call push(x_dist, distance)
    end do

    print *, x_dist

    deallocate(x_dist)

    !print *, cons_thickness
    !print *, num_cells
end program geometry_test
