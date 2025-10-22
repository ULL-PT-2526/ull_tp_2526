program integral_parallel

    use omp_lib
    use iso_fortran_env, only: real64 

    implicit none
    integer, parameter :: dp = real64 
    
    ! Integration limits and number of trapezoids
    real(kind=dp), parameter :: a = 0.0_dp
    real(kind=dp), parameter :: b = 3.0_dp
    integer, parameter :: num_steps = 10000000

    real(kind=dp) :: step_size
    real(kind=dp) :: integral_sum
    integer :: i
    real(kind=dp) :: x

    step_size = (b - a) / real(num_steps, kind=dp)

    integral_sum = 0.0_dp

    !$omp parallel do reduction(+:integral_sum) private(x)
    do i = 1, num_steps - 1
        x = a + real(i, kind=dp) * step_size
        integral_sum = integral_sum + f(x)
    end do
    !$omp end parallel do

    integral_sum = step_size * ( (f(a) + f(b))/2.0_dp + integral_sum )

    print '("The approximate value of the integral is: ", F12.8)', integral_sum

contains

    !-----------------------------------------------------------------------
    ! Function f(x). It must also use the same precision 'dp'.
    !-----------------------------------------------------------------------
    function f(x_val) result(f_val)
        ! 'dp' is inherited from the main program
        implicit none
        real(kind=dp), intent(in) :: x_val
        real(kind=dp) :: f_val

        f_val = sin(x_val) * x_val
    end function f

end program integral_parallel