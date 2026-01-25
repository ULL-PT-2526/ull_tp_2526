module particle
  use geometry
  implicit none

  type :: particle3d
    type(point3d) :: p   ! Particle position
    type(vector3d) :: v  ! Particle velocity
    real(kind=dp) :: m   ! Particle mass
  end type particle3d

end module particle