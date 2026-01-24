module particle
use iso_fortran_env, only: real64
!Creating the particle module that uses the geometry module
use geometry, only: point3d, vector3d
    implicit none
    type :: particle3d
        type(point3d) :: p  ! position
        type(vector3d) :: v   ! velocity
        real(real64) :: m   ! mass
    end type particle3d

contains
    
end module particle