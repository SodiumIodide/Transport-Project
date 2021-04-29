module dynamic_arrays
    private
    public :: &
        push_double, &
        push_int

    contains

    subroutine push_double(list, element)
        implicit none

        integer :: &
            i, isize
        real(8), intent(in) :: &
            element
        real(8), dimension(:), allocatable, intent(inout) :: &
            list
        real(8), dimension(:), allocatable :: &
            clist

        ! Match the size of the list
        if (allocated(list)) then
            isize = size(list)
            allocate(clist(isize + 1))

            ! Copy the list definitions
            do i = 1, isize
                clist(i) = list(i)
            end do

            ! Append the element
            clist(isize + 1) = element

            ! Switch memory indexing
            deallocate(list)
            call move_alloc(clist, list)
        else
            ! Allocate a new list
            allocate(list(1))
            list(1) = element
        end if
    end subroutine push_double

    subroutine push_int(list, element)
        implicit none

        integer :: &
            i, isize
        integer, intent(in) :: &
            element
        integer, dimension(:), allocatable, intent(inout) :: &
            list
        integer, dimension(:), allocatable :: &
            clist

        ! Match the size of the list
        if (allocated(list)) then
            isize = size(list)
            allocate(clist(isize + 1))

            ! Copy the list definitions
            do i = 1, isize
                clist(i) = list(i)
            end do

            ! Append the element
            clist(isize + 1) = element

            ! Switch memory indexing
            deallocate(list)
            call move_alloc(clist, list)
        else
            ! Allocate a new list
            allocate(list(1))
            list(1) = element
        end if
    end subroutine push_int
end module dynamic_arrays
