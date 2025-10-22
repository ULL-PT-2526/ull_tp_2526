! Another version of ex1.f90 with structural improvements
program ex1v2

    use iso_fortran_env, only: real64
    use geometry
    use particle
    implicit none

    ! Define a parameter for double-precision floating point
    ! Not strictly necessary here, but useful for consistency
    integer, parameter :: dp = real64

    integer :: i, n
    real(dp) :: dt, t_end, dt_out, t_out
    real(dp) :: t

    type(particle3d), allocatable :: particles(:)
    type(vector3d), allocatable :: a(:) ! Array for accelerations

    character(len=100) :: input_filename
    integer :: unit_in, unit_out

    real(dp), parameter :: epsilon2 = 1.0e-6_dp

    !====================================================================
    ! MAIN PROGRAM
    !====================================================================

    ! Read input filename from the command line
    call get_command_argument(1, input_filename)

    ! Open files using modern syntax with automatic unit assignment
    open(newunit=unit_in, file=trim(input_filename), status='old', action='read')
    open(newunit=unit_out, file='output.dat', status='replace', action='write')

    ! Read simulation parameters
    read(unit_in, *) dt
    read(unit_in, *) dt_out
    read(unit_in, *) t_end
    read(unit_in, *) n

    ! Allocate memory
    allocate(particles(n))
    allocate(a(n))

    ! Read particle data
    do i = 1, n
      read(unit_in, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                       particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    end do

    ! Compute initial accelerations (calling the subroutine)
    call compute_accelerations(particles, a, n)

    ! Initialize time variables
    t_out = 0.0_dp
    t = 0.0_dp

    ! Main loop
    do while (t <= t_end)

      ! Leapfrog integration (first part: v(t+dt/2), p(t+dt))
      do i = 1, n
        particles(i)%v = particles(i)%v + a(i) * (dt / 2.0_dp)
        particles(i)%p = particles(i)%p + particles(i)%v * dt ! Uses the new v(t+dt/2)
      end do

      ! Reset and calculate new accelerations a(t+dt)
      call compute_accelerations(particles, a, n)

      ! Complete the velocity v(t+dt)
      do i = 1, n
        particles(i)%v = particles(i)%v + a(i) * (dt / 2.0_dp)
      end do

      ! Write to the output file
      t_out = t_out + dt
      if (t_out >= dt_out) then
        write(unit_out, '(E18.8, 1X)', advance='no') t
        do i = 1, n
          write(unit_out, '(3(E18.8,1X))', advance='no') &
            particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
        end do
        write(unit_out, *)
        t_out = 0.0_dp
      end if

      t = t + dt
    end do

    ! Close files
    close(unit_in)
    close(unit_out)

    ! Deallocate memory for particles and accelerations
    deallocate(particles)
    deallocate(a)

  contains

    !====================================================================
    ! SUBROUTINE to compute accelerations for all particles
    !====================================================================

    subroutine compute_accelerations(particles, a, n, use_softening, epsilon2)
      implicit none
      integer, intent(in) :: n
      type(particle3d), intent(in) :: particles(n)
      type(vector3d), intent(out) :: a(n)
      logical, intent(in), optional :: use_softening
      real(dp), intent(in), optional :: epsilon2

      integer :: i, j
      type(vector3d) :: rij
      real(dp) :: r2, r3
      logical :: soften
      real(dp) :: eps2_val

      ! We configure the default softening values. 
      ! For the simulation with the input.dat provided, 
      ! this is unnecessary, so we leave it set to false
      soften = .false.
      eps2_val = 0.0_dp

      if (present(use_softening)) soften = use_softening
      if (present(epsilon2)) eps2_val = epsilon2

      ! Reset accelerations to zero
      do i = 1, n
        a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
      end do

      ! Compute acceleration of each particle due to all others
      do i = 1, n
        do j = i + 1, n
          rij = particles(j)%p - particles(i)%p

          ! Calculate squared distance
          r2 = distance2(particles(j)%p, particles(i)%p)

          if (soften) then
            r3 = (r2 + eps2_val)**(1.5_dp)
          else
            r3 = r2**(1.5_dp)
          end if

          a(i) = a(i) + (particles(j)%m / r3) * rij
          a(j) = a(j) - (particles(i)%m / r3) * rij
        end do
      end do

  end subroutine compute_accelerations

end program ex1v2