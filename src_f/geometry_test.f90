program geometry_test
    !use dynamic_arrays
    use mcnp_random
    use geometry_gen

    implicit none

    real(8), dimension(:), allocatable :: &
        x_dist, x_arr
    integer, dimension(:), allocatable :: &
        materials
    real(8) :: &
        thickness, chord_a, chord_b, cons_thickness!, rand_num, distance
    integer :: &
        num_cells, i
    integer(8), parameter :: &
        seed = int(123456, 8)

    thickness = 1.0d+1  ! cm
    cons_thickness = 0.0d+0  ! cm
    chord_a =  dble(101)/dble(20) ! cm
    chord_b = dble(101)/dble(20)  ! cm
    num_cells = 0

    allocate(x_dist(0))
    allocate(materials(0))

    call RN_init_problem(seed, 1)

    call get_geometry(x_dist, x_arr, materials, chord_a, chord_b, thickness, num_cells)

    !call RN_init_problem(int(1234567, 8), 1)
    !do while (cons_thickness < thickness)
    !    rand_num = rang()
    !    distance = chord_a * dlog(1.0d+0 / (1.0d+0 - rand_num))  ! cm
    !    cons_thickness = cons_thickness + distance  ! cm
    !    num_cells = num_cells + 1
    !    call push_double(x_dist, distance)
    !end do

    print *, x_dist
    print *, x_arr
    print *, materials

    open(unit=7, file="./out/geometry_test.out", form="formatted", &
         status="replace", action="write")
    do i = 1, num_cells
        write(7,*) x_arr(i), materials(i)
    end do
    close(7)

    deallocate(x_dist)
    deallocate(x_arr)
    deallocate(materials)
end program geometry_test
