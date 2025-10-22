module particle
  use geometry
  implicit none
  
  type :: particle3d
     type(point3d) :: p     ! position with a point3d
     type(vector3d) :: v    ! velocity with a vector3d
     real(kind=bit64) :: m  ! mass with a real number 
  end type particle3d

end module particle

