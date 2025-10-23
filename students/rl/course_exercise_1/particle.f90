! Create a module particle stored in particle.f90 that uses the geometry module and
! contains a type particle3d. This type should have components: a point3d variable p
! storing the particle’s position, a vector3d variable v storing the particle’s velocity, and a
! real variable m storing the particle’s mass.

module particle
  use geometry
  implicit none

  type :: particle3d
     type(point3d) :: p  ! position
     type(vector3d) :: v  ! velocity
     real(dp) :: m  ! mass
  end type particle3d

end module particle