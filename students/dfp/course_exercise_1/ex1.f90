program ex1
  use geometry
  use particle
  implicit none

  integer :: k, n, nstep
  real(kind=dp) :: dt, t_end, dt_out, t_out, t
  type(particle3d), allocatable :: particles(:)
  type(vector3d), allocatable :: accelerations(:)
  integer :: outunit

  ! Open output file
  open(newunit=outunit, file="output.dat", status="replace", action="write")

  ! Read input parameters and particle data
  call read_input(dt, dt_out, t_end, n, particles)
  
  ! Allocate acceleration array
  allocate(accelerations(n))

  ! Compute initial accelerations and half-step velocities
  call compute_accelerations(particles, n, accelerations, use_smoothing=.false.)
  call update_velocities(particles, n, accelerations, dt/2.0_dp)

  ! Initialize simulation variables
  nstep = int(t_end/dt)
  t_out = 0.0_dp
  t = 0.0_dp

  ! Write initial state
  call write_output(outunit, t, particles, n)

  ! Main simulation loop
  do k = 1, nstep
     ! Leapfrog integration step
     call update_positions(particles, n, dt)
     call compute_accelerations(particles, n, accelerations, use_smoothing=.false.)
     call update_velocities(particles, n, accelerations, dt)

     ! Update time and output if needed
     t = t + dt
     t_out = t_out + dt
     if (t_out >= dt_out .or. k == nstep) then
        call write_output(outunit, t, particles, n)
        t_out = 0.0_dp
     end if
  end do

  close(outunit)
  
  ! Deallocate arrays for good memory management practice
  deallocate(particles)
  deallocate(accelerations)

contains

  subroutine read_input(dt, dt_out, t_end, n, particles)
    implicit none
    real(kind=dp), intent(out) :: dt, dt_out, t_end
    integer, intent(out) :: n
    type(particle3d), allocatable, intent(out) :: particles(:)
    integer :: i

    read *, dt
    read *, dt_out
    read *, t_end
    read *, n

    allocate(particles(n))

    do i = 1, n
       read *, particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    end do
  end subroutine read_input

  subroutine compute_accelerations(particles, n, accelerations, use_smoothing)
    implicit none
    type(particle3d), intent(in) :: particles(:)
    integer, intent(in) :: n
    type(vector3d), intent(out) :: accelerations(:)
    logical, intent(in), optional :: use_smoothing
    real(kind=dp), parameter :: eps = 0.2_dp
    type(vector3d) :: rji, a_i
    real(kind=dp) :: r, r3
    integer :: i, j
    logical :: smooth

    smooth = .false.
    if (present(use_smoothing)) smooth = use_smoothing

    do i = 1, n
       a_i = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
       do j = 1, n
          if (i /= j) then
             rji = particles(j)%p - particles(i)%p
             ! Use distance function from geometry module
             r = distance(particles(i)%p, particles(j)%p)
             if (smooth) then
                r = sqrt(r**2 + eps**2) ! Adding smoothing if needed
             end if
             r3 = r**3
             if (r > 0.0_dp) then
                a_i = a_i + (particles(j)%m / r3) * rji
             end if
          end if
       end do
       accelerations(i) = a_i
    end do
  end subroutine compute_accelerations

  subroutine update_positions(particles, n, dt)
    implicit none
    type(particle3d), intent(inout) :: particles(:)
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: dt
    integer :: i

    do i = 1, n
       particles(i)%p = particles(i)%p + particles(i)%v * dt
    end do
  end subroutine update_positions

  subroutine update_velocities(particles, n, accelerations, dt)
    implicit none
    type(particle3d), intent(inout) :: particles(:)
    integer, intent(in) :: n
    type(vector3d), intent(in) :: accelerations(:)
    real(kind=dp), intent(in) :: dt
    integer :: i

    do i = 1, n
       particles(i)%v = particles(i)%v + accelerations(i) * dt
    end do
  end subroutine update_velocities

  subroutine write_output(outunit, t, particles, n)
    implicit none
    integer, intent(in) :: outunit
    real(kind=dp), intent(in) :: t
    type(particle3d), intent(in) :: particles(:)
    integer, intent(in) :: n
    integer :: i

    write(outunit, '(F10.5)', advance='no') t
    do i = 1, n
       write(outunit, '(3F15.8)', advance='no') particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
    end do
    write(outunit, *)
  end subroutine write_output

end program ex1