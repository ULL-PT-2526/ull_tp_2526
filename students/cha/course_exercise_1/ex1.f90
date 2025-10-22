program ex1

  use iso_fortran_env, only: real64
  use geometry
  use particle
  implicit none

  integer :: i, j, n
  real(real64) :: dt, t_end, dt_out, t_out
  real(real64) :: t

  type(particle3d), allocatable :: particles(:)
  type(vector3d), allocatable :: a(:) ! Array for accelerations

  ! Variable to manage input and output files
  character(len=100) :: input_filename

  ! Variables for the calculation
  type(vector3d) :: rij
  real(real64) :: r2, r3
  real(real64), parameter :: epsilon = 1.0e-6

  ! Read the filename from the command line
  call get_command_argument(1, input_filename)

  ! Open input file 
  open(10, file=trim(input_filename), status='old', action='read')

  ! Open output file for writing
  open(20, file='output.dat', status='replace', action='write')

  ! Read the simulation parameters from the file
  read(10,*) dt
  read(10,*) dt_out
  read(10,*) t_end
  read(10,*) n

  ! Allocate memory for particles and accelerations
  allocate(particles(n))
  allocate(a(n))

  ! Read data for each particle from the file
  do i = 1, n
    read(10,*) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
               particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
  end do

  ! Initialize accelerations to zero
  do i = 1, n
      a(i) = vector3d(0.0, 0.0, 0.0)
  end do

  ! Calculate initial accelerations
  do i = 1, n
    do j = i + 1, n
      rij = particles(j)%p - particles(i)%p
      r2 = rij%x**2 + rij%y**2 + rij%z**2
      r3 = (r2 + epsilon)**(1.5)  ! epsilon avoids possible singularities
      a(i) = a(i) + (particles(j)%m / r3) * rij
      a(j) = a(j) - (particles(i)%m / r3) * rij
    end do
  end do

  t_out = 0.0
  t = 0.0

  ! Main loop
  do while (t <= t_end)

    ! Leapfrog integration (first part: v(t+dt/2), p(t+dt))
    do i = 1, n
      particles(i)%v = particles(i)%v + a(i) * (dt / 2.0)
      particles(i)%p = particles(i)%p + particles(i)%v * dt ! Uses the new v(t+dt/2)
    end do

    ! Reset and calculate new accelerations a(t+dt)
    do i = 1, n
      a(i) = vector3d(0.0, 0.0, 0.0)
    end do
    
    do i = 1, n
      do j = i + 1, n
        rij = particles(j)%p - particles(i)%p
        r2 = rij%x**2 + rij%y**2 + rij%z**2
        r3 = (r2 + epsilon)**(1.5)
        a(i) = a(i) + (particles(j)%m / r3) * rij
        a(j) = a(j) - (particles(i)%m / r3) * rij
      end do
    end do
      
    ! Complete the velocity v(t+dt)
    do i = 1, n
      particles(i)%v = particles(i)%v + a(i) * (dt / 2.0)
    end do

    ! Write to the output file
    t_out = t_out + dt
    if (t_out >= dt_out) then
      write(20, '(E18.8, 1X)', advance='no') t
      do i = 1, n
        write(20, '(3(E18.8,1X))', advance='no') particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
      end do
      write(20,*) 
      t_out = 0.0
    end if

    t = t + dt
  end do
  
  close(10)
  close(20)

  ! Deallocate memory for particles and accelerations
  deallocate(particles)
  deallocate(a)

end program ex1
