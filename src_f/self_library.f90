module self_library
    implicit none
    public :: legendre_gauss_quad, linspace, PI

    double precision, parameter :: PI = 4.d0 * datan(1.d0)

    contains

    subroutine linspace(arr_x, x_start, x_end, x_len)
        double precision, dimension(:), intent(out) :: arr_x
        double precision, intent(in) :: x_start, x_end
        double precision :: delta_x
        integer :: x_len, i
        delta_x = (x_end - x_start) / (x_len - 1)
        arr_x(1:x_len) = [(x_start + (i-1)*delta_x, i=1,x_len)]
    end subroutine

    subroutine legendre_gauss_quad(order, intv_a, intv_b, values, weights)
        ! Compute Legendre-Gauss values and weights on an interval
        ! [intv_a, intv_b] with an input truncation order
        ! Inputs and outputs
        integer, intent(in) :: order
        double precision, intent(in) :: intv_a, intv_b
        double precision, dimension(order), intent(out) :: values, weights
        ! Vector variables
        double precision, dimension(order) :: x_space, y_space, y_hold, prime
        ! Scalar variables
        integer :: order_place, order_one, order_two, i
        ! Initial parameter for constant
        double precision, parameter :: y_const = 2.d0
        ! Legendre-Gauss Vendermonde matrix and its derivative
        double precision, dimension(order, order+1) :: legendre
        ! Test for truncation order
        if (order < 1) then
            print *, "Truncation order must be an integer greater than 0"
            call exit(1)
        end if
        ! Used for indexing purposes
        order_place = order - 1
        order_one = order_place + 1
        order_two = order_place + 2
        call linspace(x_space, -1.d0, 1.d0, order_one)
        do i = 0,order_place
            y_space(i+1) = dcos((2.d0*i + 1.d0)/(2.d0*order_place + 2.d0)*PI)+&
                     (0.27d0/order_one)*dsin((PI*x_space(i+1)*order_place)/&
                                             order_two)
        end do
        legendre(:,:) = 0.d0
        ! Compute the zeros of the N+1 Legendre Polynomial using the recursion
        ! relation and the Newton-Raphson method
        y_hold(:) = y_const
        do while (maxval(dabs(y_space - y_hold)) > epsilon(y_const))
            ! First Legendre Polynomial
            legendre(:,1) = 1.d0
            ! Second Legendre Polynomial
            if (order_one > 1) then
                legendre(:,2) = y_space
            end if
            ! Derivative calculations (only necessary for uppermost order derivative)
            if (order_one == 1) then
                prime(:) = 0.d0
            else if (order_one == 2) then
                prime(:) = 1.d0
            else
                ! Remaining Legendre Polynomials
                do i = 2,order_one
                    legendre(:,i+1) = ((2.d0*i - 1.d0)*y_space*legendre(:,i)-&
                                       (i-1.d0)*legendre(:,i-1))/i
                end do
                prime = order_two*(legendre(:,order_one)-&
                                   y_space*legendre(:,order_two))/(1-y_space**2)
            end if
            y_hold = y_space
            y_space(:) = y_hold - legendre(:,order_two) / prime
        end do
        ! Linear map from [-1, 1] to [intv_a, intv_b]
        values = intv_a * (1.d0 - y_space)/2.d0 + intv_b*(1.d0 + y_space)/2.d0
        ! Compute the weights
        weights = dble(intv_b - intv_a) / ((1.d0-y_space**2)*prime**2) * &
                  (dble(order_two)/dble(order_one))**2
    end subroutine
end module
