module geometry_gen
    use mcnp_random
    use dynamic_arrays

    private
    public :: &
        get_geometry

    contains

    ! Currently written for 2 materials
    subroutine get_geometry(x_dist, x_arr, materials, chord_a, chord_b, thickness, num_cells)
        implicit none

        ! x_dist measures the delta_x values
        ! x_arr measures the overall distance
        real(8), dimension(:), allocatable, intent(out) :: &
            x_dist, x_arr
        integer, dimension(:), allocatable, intent(out) :: &
            materials
        real(8), intent(in) :: &
            chord_a, chord_b, thickness
        integer, intent(out) :: &
            num_cells
        integer :: &
            material_num, i
        real(8) :: &
            cons_thickness, rand_num, distance, chord, prob_a
        ! The total number of cells to utilize for each geometry segment
        integer, parameter :: &
            num_divs = 10

        !if (thickness > 5.0d+0) then
        !    num_divs = 100
        !else
        !    num_divs = 1000
        !end if

        ! Determine first material to use
        prob_a = chord_a / (chord_a + chord_b)
        rand_num = rang()
        if (rand_num < prob_a) then
            material_num = 1
        else
            material_num = 2
        end if

        ! Updateable values
        cons_thickness = 0.0d+0  ! cm
        distance = 0.0d+0  ! cm

        ! Checks and (de)allocations
        if (allocated(x_dist)) then
            deallocate(x_dist)
        end if
        allocate(x_dist(0))

        if (allocated(materials)) then
            deallocate(materials)
        end if
        allocate(materials(0))

        if (allocated(x_arr)) then
            deallocate(x_arr)
        end if
        allocate(x_arr(0))

        if (num_cells > 0) then
            num_cells = 0
        end if

        ! Assign geometry via Markovian generation
        do while (cons_thickness < thickness)
            ! Generate random number
            rand_num = rang()

            ! Assign a chord length based on material number
            if (material_num == 1) then
                chord = chord_a  ! cm
            else
                chord = chord_b  ! cm
            end if

            ! Calculate and append the material length
            distance = chord * dlog(1.0d+0 / (1.0d+0 - rand_num))  ! cm
            cons_thickness = cons_thickness + distance  ! cm

            ! Check on thickness to not overshoot boundary
            if (cons_thickness > thickness) then
                distance = distance + thickness - cons_thickness  ! cm
                cons_thickness = thickness  ! cm
            end if

            ! Further discretize geometry
            do i = 1, num_divs
                call push_double(x_dist, (distance / dble(num_divs)))
                call push_double(x_arr, (cons_thickness - distance + (distance / dble(num_divs) * dble(i))))
                call push_int(materials, material_num)
                num_cells = num_cells + 1
            end do

            ! Update material number
            if (material_num == 1) then
                material_num = 2
            else
                material_num = 1
            end if
        end do
    end subroutine get_geometry

end module geometry_gen
