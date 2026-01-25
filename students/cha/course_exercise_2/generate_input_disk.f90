program generate_input_disk
  !============================================================
  ! Generate initial conditions for a thin disk around a central mass
  ! The disk radius scales with sqrt(N) to keep constant surface density.
  !============================================================

  use iso_fortran_env, only: real64
  implicit none

  ! Precision
  integer, parameter :: dp = real64

  ! Simulation parameters
  integer, parameter :: n_particles = 10000
  
  real(dp) :: disk_radius 
  real(dp), parameter :: density_scale_factor = 0.16_dp  ! we set a constant density for each generated disk


  character(len=64) :: filename 

  ! Input file parameters
  real(dp), parameter :: dt = 0.001_dp
  real(dp), parameter :: dt_out = 0.01_dp
  real(dp), parameter :: t_end = 10.0_dp

  real(dp) :: particle_mass  ! we don't set it as parameter because we are going to recalculate it depending on the
                             ! number of particles, in order to have a normalized-mass disk
  real(dp), parameter :: M_disk_total = 0.1_dp  ! we set the total disk mass to 5% of the BH mass
  real(dp), parameter :: M_central = 1.0_dp
  real(dp), parameter :: G = 1.0_dp

  real(dp), parameter :: z_dispersion = 0.05_dp
  real(dp), parameter :: vz_dispersion = 0.05_dp

  real(dp), allocatable :: pos_x(:), pos_y(:), pos_z(:)
  real(dp), allocatable :: vel_x(:), vel_y(:), vel_z(:)
  real(dp), allocatable :: radius(:), angle(:)
  integer :: i
  real(dp) :: velocity_mag

  call random_seed()

  ! R = scale_factor * sqrt(N)
  disk_radius = density_scale_factor * sqrt(real(n_particles, dp))
  if (disk_radius < 2.0_dp) disk_radius = 2.0_dp

  ! we calculate the mass of each individual particle
  particle_mass = M_disk_total / real(n_particles, dp)

  print *, "Generating system with ", n_particles, " particles"
  write(*, '(A, ES9.2)') "Particle Mass: ", particle_mass

  write(filename, '(a,i0,a)') "input_disk_", n_particles, ".dat"

  allocate(pos_x(n_particles), pos_y(n_particles), pos_z(n_particles))
  allocate(vel_x(n_particles), vel_y(n_particles), vel_z(n_particles))
  allocate(radius(n_particles), angle(n_particles))

  call random_number(radius)
  radius = sqrt(radius) * disk_radius

  call random_number(angle)
  angle = angle * 2.0_dp * acos(-1.0_dp) 

  call random_number(pos_z)
  pos_z = (pos_z - 0.5_dp) * 2.0_dp * z_dispersion 

  pos_x = radius * cos(angle)
  pos_y = radius * sin(angle)

  velocity_mag = 0.0_dp
  do i = 1, n_particles
      velocity_mag = sqrt(G * M_central / max(radius(i), 0.05_dp))
      
      vel_x(i) = -velocity_mag * sin(angle(i))
      vel_y(i) =  velocity_mag * cos(angle(i))
      
      call random_number(vel_z(i))
      vel_z(i) = (vel_z(i) - 0.5_dp) * 2.0_dp * vz_dispersion
  end do

  ! Write file
  open(unit=10, file=trim(filename), status='replace', action='write')

  write(10, '(ES12.3)') dt
  write(10, '(ES12.3)') dt_out
  write(10, '(F10.3)') t_end
  write(10, '(F10.3)') disk_radius
  write(10, *) n_particles+1

  write(10, '(F12.6,6F12.6)') M_central, 0.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,0.0_dp

  do i = 1, n_particles
      write(10, '(F12.6,6F12.6)') particle_mass, pos_x(i), pos_y(i), pos_z(i), &
                                  vel_x(i), vel_y(i), vel_z(i)
  end do

  close(10)

  print *, "File '", trim(filename), "' successfully generated."

  deallocate(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, radius, angle)

end program generate_input_disk