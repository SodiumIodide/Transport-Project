module mesh_map
    private
    public :: &
        struct_to_unstruct, &
        unstruct_to_struct

    contains

    subroutine struct_to_unstruct(struct, struct_delta, struct_size, unstruct, unstruct_delta, unstruct_size)
        integer, intent(in) :: &
            struct_size, unstruct_size
        double precision, dimension(:), allocatable, intent(in) :: &
            unstruct_delta
        double precision, dimension(struct_size), intent(in) :: &
            struct
        double precision, dimension(:, :), allocatable, intent(inout) :: &
            unstruct
        double precision, intent(in) :: &
            struct_delta
        integer :: &
            i, counter
        double precision :: &
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
        allocate(unstruct(1, unstruct_size))

        ! Loop for mapping the structured to the unstructured mesh
        do i = 1, unstruct_size
            ! Start from the border with a zero weight for the current cell
            distance_overlap = .false.
            weight_tally = 0.0d+0  ! unit*cm

            ! Increment unstructured distance tally
            unstruct_distance_tally = unstruct_distance_tally + unstruct_delta(i)  ! cm

            ! Carry over leftover distance
            if (leftover_distance /= 0.0d+0) then
                weight_tally = weight_tally + leftover_distance * struct(counter)  ! unit*cm
                distance_tally = distance_tally + leftover_distance  ! cm
                leftover_distance = 0.0d+0  ! cm
            end if

            do while (.not. distance_overlap)
                ! Increment counter (structured index)
                counter = counter + 1

                ! Increment structured distance tally
                struct_distance_tally = struct_distance_tally + struct_delta  ! cm

                ! Check for boundary overlap
                if ((distance_tally + struct_delta) > unstruct_distance_tally) then
                    delta = unstruct_distance_tally - (distance_tally + struct_delta)  ! cm
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
            end do
            unstruct(1, i) = weight_tally / unstruct_delta(i)  ! unit
        end do
    end subroutine struct_to_unstruct

    subroutine unstruct_to_struct(unstruct, unstruct_delta, struct, struct_delta, struct_size)
        integer, intent(in) :: &
            struct_size
        double precision, dimension(:), allocatable, intent(in) :: &
            unstruct_delta
        double precision, dimension(:, :), allocatable, intent(in) :: &
            unstruct
        double precision, dimension(struct_size), intent(inout) :: &
            struct
        double precision, intent(in) :: &
            struct_delta
        integer :: &
            i, counter
        double precision :: &
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

            ! Increment unstructured distance tally
            struct_distance_tally = struct_distance_tally + struct_delta  ! cm

            ! Carry over leftover distance
            if (leftover_distance /= 0.0d+0) then
                weight_tally = weight_tally + leftover_distance * unstruct(1, counter)  ! unit*cm
                distance_tally = distance_tally + leftover_distance  ! cm
                leftover_distance = 0.0d+0  ! cm
            end if

            do while (.not. distance_overlap)
                ! Increment counter (structured index)
                counter = counter + 1

                ! Increment structured distance tally
                unstruct_distance_tally = unstruct_distance_tally + unstruct_delta(counter)  ! cm

                ! Check for boundary overlap
                if ((distance_tally + unstruct_delta(counter)) > struct_distance_tally) then
                    delta = struct_distance_tally - (distance_tally + unstruct_delta(counter))  ! cm
                    leftover_distance = unstruct_distance_tally - struct_distance_tally  ! cm
                    distance_overlap = .true.
                else
                    delta = unstruct_delta(counter)  ! cm
                end if

                ! Increment the known distance tally
                distance_tally = distance_tally + delta  ! cm

                ! Apply linear weighted tally
                weight_tally = weight_tally + delta * unstruct(1, counter)  ! unit*cm

                ! Reset counter to use leftover distance and re-use value
                if (distance_overlap) then
                    counter = counter - 1
                end if
            end do
            struct(i) = weight_tally / struct_delta  ! unit
        end do
    end subroutine unstruct_to_struct
end module mesh_map
