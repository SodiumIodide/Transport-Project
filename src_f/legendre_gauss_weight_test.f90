program legendre_gauss_weight
    use self_library

    implicit none

    double precision, parameter :: intv_a = 1.d0
    double precision, parameter :: intv_b = 2.d0
    integer, parameter :: order = 5
    double precision, dimension(order) :: values, weights
    integer :: i
    call legendre_gauss_quad(order, intv_a, intv_b, values, weights)
    print *, "Values:"
    do i = 1,order
        print *, values(i)
    end do
    print *, "Weights:"
    do i = 1,order
        print *, weights(i)
    end do
end program
