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

    subroutine get_geometry(x_dist, materials, chord_a, chord_b, thickness, num_cells)
        implicit none

        double precision, dimension(:), allocatable, intent(inout) :: &
            x_dist
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

        if (allocated(x_dist)) then
            deallocate(x_dist)
        end if
        allocate(x_dist(0))

        if (num_cells > 0) then
            num_cells = 0
        end if

        cons_thickness = 0.0d+0  ! cm

        material_num = init_material

        call RN_init_problem(int(1234567, 8), 1)
        do while (cons_thickness < thickness)
            rand_num = rang()
            if (material_num == 1) then
                chord = chord_a
            else
                chord = chord_b
            end if
            distance = chord * dlog(1.0d+0 / (1.0d+0 - rand_num))  ! cm
            cons_thickness = cons_thickness + distance  ! cm
            call push_double(x_dist, distance)
            call push_int(materials, material_num)
            if (material_num == init_material) then
                material_num = num_materials
            else
                material_num = init_material
            end if
            num_cells = num_cells + 1
        end do

    end subroutine get_geometry

end module geometry_gen
