module mesh_map
    contains

    subroutine struct_to_unstruct(struct, struct_delta, struct_size, unstruct, unstruct_delta, unstruct_size, new_unstruct)
        integer, intent(in) :: &
            struct_size, unstruct_size
        double precision, dimension(:), allocatable, intent(in) :: &
            unstruct, unstruct_delta
        double precision, dimension(struct_size), intent(in) :: &
            struct
        double precision, dimension(:), allocatable, intent(inout) :: &
            new_unstruct
        double precision, intent(in) :: &
            struct_delta
        integer :: &
            i

        ! Allocate the new unstructured array map
        allocate(new_unstruct(unstruct_size))
        do i = 1, num_cells
            
        end do
    end subroutine struct_to_unstruct
end module mesh_map
