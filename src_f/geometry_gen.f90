module geometry_gen
    use mcnp_random
    use dynamic_arrays

    private
    public :: &
        get_geometry

    integer :: &
        num_materials = 2
    integer :: &
        init_material = 1

    contains

    subroutine get_geometry(x_dist, x_arr, materials, chord_a, chord_b, thickness, num_cells)
        implicit none

        ! x_dist measures the delta_x values
        ! x_arr measures the overall distance
        double precision, dimension(:), allocatable, intent(inout) :: &
            x_dist, x_arr
        integer, dimension(:), allocatable, intent(inout) :: &
            materials
        double precision, intent(in) :: &
            chord_a, chord_b, thickness
        integer, intent(inout) :: &
            num_cells
        integer :: &
            material_num
        double precision :: &
            cons_thickness, rand_num, distance, chord

        ! Updateable values
        cons_thickness = 0.0d+0  ! cm
        material_num = init_material
        distance = 0.0d+0  ! cm

        ! Overwrite counter variable
        if (num_cells > 0) then
            num_cells = 0
        end if

        ! Checks and (de)allocations
        ! Reallocate to 1 for initial assignment
        ! (know that arrays will be at least size 1)
        if (allocated(x_dist)) then
            deallocate(x_dist)
        end if
        allocate(x_dist(1))

        if (allocated(materials)) then
            deallocate(materials)
        end if
        allocate(materials(1))

        if (allocated(x_arr)) then
            deallocate(x_arr)
        end if
        allocate(x_arr(1))

        ! Assign first unit onto newly-allocated arrays
        x_dist(1) = distance  ! cm
        x_arr(1) = cons_thickness  ! cm
        materials(1) = material_num

        ! Start random generation seed
        call RN_init_problem(int(12345678, 8), 1)

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

            ! Update material number
            if (material_num == init_material) then
                material_num = num_materials
            else
                material_num = init_material
            end if

            ! Calculate and append the material length
            distance = chord * dlog(1.0d+0 / (1.0d+0 - rand_num))  ! cm
            cons_thickness = cons_thickness + distance  ! cm

            ! Check on thickness to not overshoot boundary
            if (cons_thickness > thickness) then
                distance = cons_thickness - thickness  ! cm
                cons_thickness = thickness  ! cm
            end if

            ! Push results onto arrays
            call push_double(x_dist, distance)
            call push_double(x_arr, cons_thickness)
            call push_int(materials, material_num)

            ! Tally cells
            num_cells = num_cells + 1
        end do
    end subroutine get_geometry

end module geometry_gen
