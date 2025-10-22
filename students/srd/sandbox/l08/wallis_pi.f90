program wallis_pi
  implicit none
  integer :: i, n
  real(8) :: pi, term
  real(8), parameter :: pi_real = 3.141592653589793d0

  ! Number of steps
  n = 10000

  pi = 1.0d0

  !$omp parallel do private(term) reduction(*:pi)
  do i = 1, n
     term = (4.0d0 * i**2) / (4.0d0 * i**2 - 1.0d0)
     pi = pi * term
  end do
  !$omp end parallel do

  pi = 2.0d0 * pi

  print '(A, F20.15)', 'Approximation to pi = ', pi
  print '(A, E15.5)',  'Absolute error       = ', abs(pi - pi_real)

end program wallis_pi

