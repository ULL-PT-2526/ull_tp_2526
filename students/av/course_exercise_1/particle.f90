module particle

    ! Using the geometry module created for vector/real operations
    use geometry 
    implicit none

    ! Making private internal operators and functions # Lecture 4
    private

    ! Needed types are public # Lecture 4
    public :: particle3d

        ! Defining particle3d type # Lecture 5
        type particle3d
            type(point3d)  :: p ! Particle's position
            type(vector3d) :: v ! Particle's velocity
            real(kind=dp)  :: m ! Particle's mass
        end type particle3d 

end module particle