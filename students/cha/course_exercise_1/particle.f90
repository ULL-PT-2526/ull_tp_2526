module particle

  use iso_fortran_env, only: real64
  use geometry
  implicit none

  type :: particle3d
    real(real64) :: m     ! Particle mass
    type(point3d) :: p    ! Particle position
    type(vector3d) :: v   ! Particle velocity
  end type particle3d

end module particle