module mesh_map
    private
    public :: &
        struct_to_unstruct, &
        unstruct_to_struct, &
        material_calc

    contains

    subroutine struct_to_unstruct(struct, struct_delta, struct_size, unstruct, unstruct_delta, unstruct_size)
        implicit none

        integer, intent(in) :: &
            struct_size, unstruct_size
        ! Allocated as (num_ind_cells)
        real(8), dimension(:), allocatable, intent(in) :: &
            unstruct_delta
        real(8), dimension(struct_size), intent(in) :: &
            struct
        real(8), dimension(:), allocatable, intent(out) :: &
            unstruct
        real(8), intent(in) :: &
            struct_delta
        integer :: &
            i, counter
        real(8) :: &
            distance_tally, unstruct_distance_tally, struct_distance_tally, weight_tally, delta, leftover_distance
        logical :: &
            distance_overlap

        ! Assign counter to zero due to pre-fetch increment
        counter = 0
        ! Assign tallies
        distance_tally = 0.0d+0  ! cm
        unstruct_distance_tally = 0.0d+0  ! cm
        struct_distance_tally = 0.0d+0  ! cm
        leftover_distance = 0.0d+0  ! cm

        ! Allocate the new unstructured array map
        if (allocated(unstruct)) then
            deallocate(unstruct)
        end if
        allocate(unstruct(unstruct_size))

        ! Loop for mapping the structured to the unstructured mesh
        do i = 1, unstruct_size
            ! Start from the border with a zero weight for the current cell
            distance_overlap = .false.
            weight_tally = 0.0d+0  ! unit*cm

            ! Increment unstructured distance tally
            unstruct_distance_tally = unstruct_distance_tally + unstruct_delta(i)  ! cm

            ! Carry over leftover distance
            if (leftover_distance /= 0.0d+0) then
                ! If leftover_distance is still over-reaching the tally boundaries
                if ((distance_tally + leftover_distance) >= unstruct_distance_tally) then
                    delta = unstruct_delta(i)  ! cm
                    leftover_distance = struct_distance_tally - unstruct_distance_tally  ! cm
                    distance_overlap = .true.
                else
                    delta = leftover_distance  ! cm
                    leftover_distance = 0.0d+0  ! cm
                end if
                weight_tally = weight_tally + delta * struct(counter)  ! unit*cm
                distance_tally = distance_tally + delta  ! cm
            end if

            do while ((.not. distance_overlap) .and. (counter < struct_size))
                ! Increment counter (structured index)
                counter = counter + 1

                ! Increment structured distance tally
                struct_distance_tally = struct_distance_tally + struct_delta  ! cm

                ! Check for boundary overlap
                if (struct_distance_tally >= unstruct_distance_tally) then
                    delta = unstruct_distance_tally - distance_tally  ! cm
                    leftover_distance = struct_distance_tally - unstruct_distance_tally  ! cm
                    distance_overlap = .true.
                else
                    delta = struct_delta  ! cm
                end if

                ! Increment the known distance tally
                distance_tally = distance_tally + delta  ! cm

                ! Apply linear weighted tally
                weight_tally = weight_tally + delta * struct(counter)  ! unit*cm

                ! Reset counter to use leftover distance and re-use value
                if (distance_overlap) then
                    counter = counter - 1
                end if
            end do  ! Structured loop
            unstruct(i) = weight_tally / unstruct_delta(i)  ! unit
        end do  ! Unstructured loop
    end subroutine struct_to_unstruct

    subroutine unstruct_to_struct(unstruct, unstruct_delta, unstruct_size, struct, struct_delta, struct_size)
        implicit none

        integer, intent(in) :: &
            struct_size, unstruct_size
        ! Allocated as (num_ind_cells)
        real(8), dimension(:), allocatable, intent(in) :: &
            unstruct_delta, unstruct
        real(8), dimension(struct_size), intent(out) :: &
            struct
        real(8), intent(in) :: &
            struct_delta
        integer :: &
            i, counter
        real(8) :: &
            distance_tally, unstruct_distance_tally, struct_distance_tally, weight_tally, delta, leftover_distance
        logical :: &
            distance_overlap

        ! Assign counter to zero due to pre-fetch increment
        counter = 0
        ! Assign tallies
        distance_tally = 0.0d+0  ! cm
        unstruct_distance_tally = 0.0d+0  ! cm
        struct_distance_tally = 0.0d+0  ! cm
        leftover_distance = 0.0d+0  ! cm

        ! Structured array does not require an allocatable assignment

        ! Loop for mapping the unstructured to the structured mesh
        do i = 1, struct_size
            ! Start from the border with a zero weight for the current cell
            distance_overlap = .false.
            weight_tally = 0.0d+0  ! unit*cm

            ! Increment structured distance tally
            struct_distance_tally = struct_distance_tally + struct_delta  ! cm

            ! Carry over leftover distance
            if (leftover_distance /= 0.0d+0) then
                ! If leftover_distance is still over-reaching the tally boundaries
                if ((distance_tally + leftover_distance) >= struct_distance_tally) then
                    delta = struct_delta  ! cm
                    leftover_distance = unstruct_distance_tally - struct_distance_tally  ! cm
                    distance_overlap = .true.
                else
                    delta = leftover_distance  ! cm
                    leftover_distance = 0.0d+0  ! cm
                end if
                weight_tally = weight_tally + delta * unstruct(counter)  ! unit*cm
                distance_tally = distance_tally + delta  ! cm
            end if

            do while ((.not. distance_overlap) .and. (counter < unstruct_size))
                ! Increment counter (unstructured index)
                counter = counter + 1

                ! Increment structured distance tally
                unstruct_distance_tally = unstruct_distance_tally + unstruct_delta(counter)  ! cm

                ! Check for boundary overlap
                if ((unstruct_distance_tally >= struct_distance_tally) .or. (counter == unstruct_size)) then
                    delta = struct_distance_tally - distance_tally  ! cm
                    leftover_distance = unstruct_distance_tally - struct_distance_tally  ! cm
                    distance_overlap = .true.
                else
                    delta = unstruct_delta(counter)  ! cm
                    leftover_distance = 0.0d+0  ! cm
                end if

                ! Increment the known distance tally
                distance_tally = distance_tally + delta  ! cm

                ! Apply linear weighted tally
                weight_tally = weight_tally + delta * unstruct(counter)  ! unit*cm
            end do  ! Unstructured loop
            struct(i) = weight_tally / struct_delta  ! unit
        end do  ! Structured loop
    end subroutine unstruct_to_struct

    subroutine material_calc(unstruct, unstruct_delta, unstruct_size, materials, &
                             material_struct, struct_delta, struct_size, &
                             & num_materials)
        implicit none

        integer, intent(in) :: &
            struct_size, num_materials, unstruct_size
        real(8), dimension(:), allocatable, intent(in) :: &
            unstruct_delta
        ! Allocated as (num_ind_cells)
        real(8), dimension(:), allocatable, intent(in) :: &
            unstruct
        ! Allocated as (num_ind_cells)
        integer, dimension(:), allocatable, intent(in) :: &
            materials
        real(8), intent(in) :: &
            struct_delta
        real(8), dimension(struct_size, num_materials), intent(out) :: &
            material_struct
        integer :: &
            i, k, counter
        real(8) :: &
            distance_tally, unstruct_distance_tally, struct_distance_tally, weight_tally, delta, leftover_distance, switch
        logical :: &
            distance_overlap

        do k = 1, num_materials
            ! Assign counter to zero due to pre-fetch increment
            counter = 0
            ! Assign tallies
            distance_tally = 0.0d+0  ! cm
            unstruct_distance_tally = 0.0d+0  ! cm
            struct_distance_tally = 0.0d+0  ! cm
            leftover_distance = 0.0d+0  ! cm
            switch = 0.0d+0

            ! Structured array does not require an allocatable assignment

            ! Loop for mapping the unstructured to the structured mesh
            do i = 1, struct_size
                ! Start from the border with a zero weight for the current cell
                distance_overlap = .false.
                weight_tally = 0.0d+0  ! unit*cm

                ! Increment unstructured distance tally
                struct_distance_tally = struct_distance_tally + struct_delta  ! cm

                ! Carry over leftover distance
                if (leftover_distance > 0.0d+0) then
                    ! If leftover distance is still over-reaching tally boundaries
                    if ((distance_tally + leftover_distance) >= struct_distance_tally) then
                        delta = struct_delta  ! cm
                        leftover_distance = unstruct_distance_tally - struct_distance_tally  ! cm
                        distance_overlap = .true.
                    else
                        delta = leftover_distance  ! cm
                        leftover_distance = 0.0d+0  ! cm
                    end if
                    weight_tally = weight_tally + delta * unstruct(counter) * switch  ! unit*cm
                    distance_tally = distance_tally + delta  ! cm
                end if

                do while ((.not. distance_overlap) .and. (counter < unstruct_size))
                    ! Increment counter (unstructured index)
                    counter = counter + 1
                    ! Increment structured distance tally
                    unstruct_distance_tally = unstruct_distance_tally + unstruct_delta(counter)  ! cm

                    ! Material number for calculations
                    ! Tally switch for each material
                    if (materials(counter) == k) then
                        switch = 1.0d+0
                    else
                        switch = 0.0d+0
                    end if

                    ! Check for boundary overlap
                    if ((unstruct_distance_tally >= struct_distance_tally) .or. (counter == unstruct_size)) then
                        delta = struct_distance_tally - distance_tally  ! cm
                        leftover_distance = unstruct_distance_tally - struct_distance_tally  ! cm
                        distance_overlap = .true.
                    else
                        delta = unstruct_delta(counter)  ! cm
                    end if

                    ! Increment the known distance tally
                    distance_tally = distance_tally + delta  ! cm

                    ! Apply linear weighted tally
                    weight_tally = weight_tally + delta * unstruct(counter) * switch  ! unit*cm
                end do  ! Unstructured loop

                ! Average the results, or just append if no results previously
                material_struct(i, k) = material_struct(i, k) + (weight_tally / struct_delta)  ! unit
            end do  ! Structured loop
        end do  ! Material loop
    end subroutine material_calc
end module mesh_map
