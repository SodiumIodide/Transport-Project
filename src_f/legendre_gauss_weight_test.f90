program legendre_gauss_weight
    use self_library

    implicit none

    real(8), parameter :: intv_a = 1.d0
    real(8), parameter :: intv_b = 2.d0
    integer, parameter :: order = 5
    real(8), dimension(order) :: values, weights
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
end program legendre_gauss_weight
