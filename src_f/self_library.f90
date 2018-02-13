module self_library
    private
    public :: &
        legendre_gauss_quad, &
        linspace, &
        gauleg, &
        PI

    real(8), parameter :: &
        PI = 4.d0 * datan(1.0d+0)

    contains

    subroutine linspace(arr_x, x_start, x_fin, x_len)
        implicit none

        real(8), dimension(:), intent(out) :: &
            arr_x
        real(8), intent(in) :: &
            x_start, x_fin
        real(8) :: &
            delta_x
        integer :: &
            x_len, i

        delta_x = (x_fin - x_start) / (x_len - 1)
        arr_x(1:x_len) = [(x_start + (i - 1) * delta_x, i = 1, x_len)]
    end subroutine linspace

    subroutine legendre_gauss_quad(order, intv_a, intv_b, values, weights)
        implicit none

        ! Compute Legendre-Gauss values and weights on an interval
        ! [intv_a, intv_b] with an input truncation order
        ! Inputs and outputs
        integer, intent(in) :: &
            order
        real(8), intent(in) :: &
            intv_a, intv_b
        real(8), dimension(order), intent(out) :: &
            values, weights
        ! Vector variables
        real(8), dimension(order) :: &
            x_space, y_space, y_hold, prime
        ! Scalar variables
        integer :: &
            order_place, order_one, order_two, i
        ! Initial parameter for constant
        real(8), parameter :: &
            y_const = 2.d0
        ! Legendre-Gauss Vendermonde matrix and its derivative
        real(8), dimension(order, order + 1) :: &
            legendre

        ! Test for truncation order
        if (order < 1) then
            print *, "Truncation order must be an integer greater than 0"
            call exit(1)
        end if

        ! Used for indexing purposes
        order_place = order - 1
        order_one = order_place + 1
        order_two = order_place + 2

        ! Assign the arrays
        call linspace(x_space, -1.0d+0, 1.0d+0, order_one)
        do i = 0, order_place
            y_space(i + 1) = dcos((2.0d+0 * i + 1.0d+0) / (2.0d+0 * order_place + 2.0d+0) * PI) + &
                     (0.27d+0 / order_one) * dsin((PI * x_space(i + 1) * order_place) / &
                                             order_two)
        end do
        legendre(:, :) = 0.0d+0

        ! Compute the zeros of the N+1 Legendre Polynomial using the recursion
        ! relation and the Newton-Raphson method
        y_hold(:) = y_const
        do while (maxval(dabs(y_space - y_hold)) > epsilon(y_const))
            ! First Legendre Polynomial
            legendre(:, 1) = 1.0d+0
            ! Second Legendre Polynomial
            legendre(:, 2) = y_space

            do i = 2, order_one
                legendre(:, i + 1) = ((2.0d+0 * i - 1.0d+0) * y_space * legendre(:, i) &
                                      - (i - 1.0d+0) * legendre(:, i - 1)) / dble(i)
            end do
            
            ! Derivative of the Legendre polynomial
            prime = order_two * (legendre(:, order_one) - &
                                 y_space * legendre(:, order_two)) / (1 - y_space**2)

            y_hold = y_space
            y_space(:) = y_hold - legendre(:, order_two) / prime
        end do

        ! Linear map from [-1, 1] to [intv_a, intv_b]
        values = intv_a * (1.0d+0 - y_space) / 2.0d+0 + intv_b * (1.0d+0 + y_space) / 2.0d+0

        ! Compute the weights
        weights = dble(intv_b - intv_a) / ((1.0d+0 - y_space**2) * prime**2) * &
                  (dble(order_two) / dble(order_one))**2
    end subroutine legendre_gauss_quad

    subroutine mat_diff(mat1, mat2, num_x, num_y, err)
        implicit none

        integer, intent(in) :: &
            num_x, num_y
        real(8), dimension(num_x, num_y), intent(in) :: &
            mat1, mat2
        real(8), intent(out) :: &
            err
        integer :: &
            i, j
        real(8) :: &
            difference

        err = 0.0d+0

        do i = 1, num_x
            do j = 1, num_y
                if ((mat1(i, j) /= 0.0d+0) .and. (mat2(i, j) /= 0.0d+0)) then
                    difference = mat1(i, j) - mat2(i, j)
                    if (difference > err) then
                        err = difference
                    end if
                end if
            end do
        end do
    end subroutine mat_diff
end module self_library
