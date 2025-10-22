program leapfrog
  use geometry
  use particle
  implicit none

  integer :: i, j
  integer :: n

  real(kind=bit64) :: dt, t_end, t, dt_out, t_out ! With our defined precision of 64 bits (in geometry module)
  real(kind=bit64) :: r2, r3
  real(kind=bit64), parameter :: epsilon = 1.0e-5 ! Parameter to avoid divisions by zero
  
  type(particle3d), allocatable :: p(:) ! We re-define those using the definitions from the modules
  type(vector3d), allocatable :: a(:)
  type(vector3d) :: rji
  character(len=100) :: infile
  integer :: ios, iu

  ! The name of the input file for the initial conditions
  infile = "input_ex.dat"

  iu = 10
  open(unit=iu, file=trim(infile), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     print *, "Error: can't open file ", trim(infile)
     stop
  end if

  read(iu, *) dt 		! In the example input we have the parameters ordered as here
  read(iu, *) dt_out
  read(iu, *) t_end
  read(iu, *) n

  allocate(p(n))
  allocate(a(n))

  ! In this code we just made the appropiate changes of the old leapfrog code to make it work with the new modules "geometry" and "particles"

  ! Read particles (masses, positions and velocities) from the input file 
  do i = 1, n
     read(iu, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
                  p(i)%v%x, p(i)%v%y, p(i)%v%z
  end do

  close(iu)

  ! Calculating initial acelerations
  a = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64)
  do i = 1, n
     do j = i + 1, n
        rji%x = p(j)%p%x - p(i)%p%x ! Also possible to define a function "subpp" in geometry and use it here
        rji%y = p(j)%p%y - p(i)%p%y
        rji%z = p(j)%p%z - p(i)%p%z

        r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2 ! Classical gravitational softening form
        r3 = r2 * sqrt(r2)

        a(i) = a(i) + (p(j)%m / r3) * rji

        a(j) = a(j) - (p(i)%m / r3) * rji

     end do
  end do

  ! Creating the output file
  open(unit=20, file="output.dat", status="replace", action="write", iostat=ios)
  if (ios /= 0) then
     print *, "Error: can't create output.dat"
     stop
  end if

  ! Integration using Leapfrog method
  t_out = 0.0_bit64
  t = 0.0_bit64
  do while (t <= t_end) ! With a do while loop there's no problem on "t" and "t_end" being real
  ! do t = 0.0_bit64, t_end, dt ! A do loop with real numbers (and not integers) it's risky a gives a warning

     do i = 1, n
        p(i)%v = p(i)%v + (0.5_bit64 * dt) * a(i)
     end do

     do i = 1, n
        p(i)%p = p(i)%p + dt * p(i)%v
     end do

     a = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64)
     do i = 1, n
        do j = i + 1, n
        
           rji%x = p(j)%p%x - p(i)%p%x
           rji%y = p(j)%p%y - p(i)%p%y
           rji%z = p(j)%p%z - p(i)%p%z
           r2  = mulvv(rji, rji) + epsilon**2
           r3  = r2 * sqrt(r2)

           a(i) = a(i) + (p(j)%m / r3) * rji
           a(j) = a(j) - (p(i)%m / r3) * rji

        end do
     end do

     do i = 1, n
        p(i)%v = p(i)%v + (0.5_bit64 * dt) * a(i)
     end do

     t_out = t_out + dt
     if (t_out >= dt_out) then
        write(20,'(F12.6, 999(F12.6))') t, (p(i)%p%x, p(i)%p%y, p(i)%p%z, i=1,n)
        t_out = 0.0_bit64
     end if
     t = t + dt
  end do

  close(20)

  deallocate(p)
  deallocate(a)

end program leapfrog
