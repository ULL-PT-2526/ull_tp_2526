module particle
    use iso_fortran_env, only: real64
    use geometry
    implicit none

    type particle3d
        real(real64) :: m       ! particle's mass
        type(point3d) :: p      ! particle's position
        type(vector3d) :: v     ! particle's velocity
    end type particle3d

end module particle