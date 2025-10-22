program wallis_pi_parallel

    use omp_lib
    use iso_fortran_env, only: real64

    implicit none

    integer, parameter :: dp = real64

    ! Program variables
    integer :: n_terms           ! Number of terms to use in the product
    integer, parameter :: max_terms = 200000000 ! Limit to avoid infinite loop
    integer, parameter :: step_size = 10000     ! Increment used to check convergence
    integer :: i

    real(kind=dp) :: wallis_product  ! Variable used for reduction
    real(kind=dp) :: pi_approx
    real(kind=dp) :: term

    ! Reference value of Pi and tolerance for 5 correct digits
    real(kind=dp), parameter :: pi_true = 3.141592653589793_dp
    real(kind=dp), parameter :: tolerance = 0.000005_dp ! 0.5 * 10^-5

    ! Outer loop to find the number of terms 'n_terms'
    do n_terms = step_size, max_terms, step_size
        
        ! Reset the product for each new calculation
        wallis_product = 1.0_dp

        ! Parallel loop that computes the Wallis product up to 'n_terms'
        ! The 'reduction(*:wallis_product)' clause multiplies partial results
        ! from each thread. 'term' is a private working variable.
        !$omp parallel do reduction(*:wallis_product) private(term)
        do i = 1, n_terms
            term = (2.0_dp * i / (2.0_dp * i - 1.0_dp)) * &
                   (2.0_dp * i / (2.0_dp * i + 1.0_dp))
            wallis_product = wallis_product * term
        end do
        !$omp end parallel do

        ! Compute the Pi approximation from the formula
        pi_approx = 2.0_dp * wallis_product

        ! Check if desired precision has been reached
        if (abs(pi_approx - pi_true) < tolerance) then
            print '("Convergence reached.")'
            print *
            print '("Number of terms required for covergence: ", I10)', n_terms
            print '("Approximated Pi value:    ", F20.15)', pi_approx
            print '("True Pi value:            ", F20.15)', pi_true
            exit ! Exit the loop once convergence is achieved
        end if
        
        ! Print progress (optional)
        if (mod(n_terms, 1000000) == 0) then
            print '("Testing with ", I10, " terms. Current Pi: ", F20.15)', n_terms, pi_approx
        end if

    end do

    ! Message if convergence is not reached
    if (n_terms >= max_terms) then
        print *, "Desired precision not reached with the maximum number of terms."
        print '("Last approximation: ", F20.15, " with ", I12, " terms.")', pi_approx, n_terms
    end if

end program wallis_pi_parallel

