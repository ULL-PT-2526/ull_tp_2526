! This module contains:
! - vector3d and point3d data types
! - functions for adding and substracting points and vector, and for multiplying
!   and dividing vectors with reals, as well as among themselves
! - operators matching these functions
! - function distance: calculates the dist. between two 3d points
! - function angle: calculates the angle between two vectors
! - function normalize: takes a vector and divides it by its length
! - function cross_product: takes two 3d vectors and gives their cross product

module geometry
    implicit none
    integer, parameter :: dp = selected_real_kind(15)

    type :: vector3d
        real(dp) :: x, y, z
    end type vector3d

    type :: point3d
        real(dp) :: x, y, z
    end type point3d

    interface operator(+)
        module procedure sumvp, sumpv, sumvv, sumpp
    end interface

    interface operator(-)
        module procedure subvp, subpv, subvv, subpp 
    end interface

    interface operator(*)
        module procedure mulrv, mulvr 
    end interface

    interface operator(/)
        module procedure divvr
    end interface

contains

    ! sumvp: adds a 3D vector to a 3D point 
    pure type(point3d) function sumvp(v, p)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
    end function sumvp

    ! sumpv: adds a 3D point to a 3D vector
    pure type(point3d) function sumpv(p, v)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
    end function sumpv

    ! sumvv: adds a 3D vector to another 3D vector
    pure type(vector3d) function sumvv(v1,v2)
        type(vector3d), intent(in) :: v1
        type(vector3d), intent(in) :: v2
        sumvv = vector3d(v1%x + v2%x, v1%y + v2%y, v1%z + v2%z)
    end function sumvv

    ! sumpp: adds two 3D points, returning a vector3d
    pure type(vector3d) function sumpp(p1, p2)
        type(point3d), intent(in) :: p1, p2
        sumpp = vector3d(p1%x + p2%x, p1%y + p2%y, p1%z + p2%z)
    end function sumpp

    ! subpv: subtracts a 3D vector from a 3D point
    pure type(point3d) function subpv(p, v)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        subpv = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
    end function subpv

    ! subvp: subtracts a 3D point from a 3D vector
    pure type(point3d) function subvp(v, p)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z)
    end function subvp

    ! subpp: subtracts two 3D points, returning a vector3d
    pure type(vector3d) function subpp(p1, p2)
        type(point3d), intent(in) :: p1, p2
        subpp = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
    end function subpp

    ! subvv: substracts a 3D vector from another 3D vector
    pure type(vector3d) function subvv(v1,v2)
        type(vector3d), intent(in) :: v1
        type(vector3d), intent(in) :: v2
        subvv = vector3d(v1%x - v2%x, v1%y - v2%y, v1%z - v2%z)
    end function subvv

    ! mulrv: multiplies a real number by a 3D vector
    pure type(vector3d) function mulrv(r, v)
        real(dp), intent(in) :: r
        type(vector3d), intent(in) :: v
        mulrv = vector3d(r*v%x, r*v%y, r*v%z)
    end function mulrv

    ! mulvr: multiplies a 3D vector by a real number
    pure type(vector3d) function mulvr(v, r)
        type(vector3d), intent(in) :: v
        real(dp), intent(in) :: r
        mulvr = vector3d(v%x*r, v%y*r, v%z*r)
    end function mulvr

    ! divvr: divides a 3D vector by a real number
    pure type(vector3d) function divvr(v, r) 
        type(vector3d), intent(in) :: v
        real(dp), intent(in) :: r
            divvr = vector3d(v%x/r, v%y/r, v%z/r)
    end function divvr

    ! distance: returns distance between two 3D points
    pure real(dp) function distance(p1, p2)
        type(point3d), intent(in) :: p1, p2
        distance = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z))
    end function distance
        
    ! norm: calculates the norm of a 3D vector
    pure real(dp) function norm(v)
        type(vector3d), intent(in) :: v
        norm = sqrt(v%x**2 + v%y**2 + v%z**2)
    end function norm

    ! normalize: normalizes a 3D vector into a unit vector
    pure type(vector3d) function normalize(v) 
        type(vector3d), intent(in)  :: v
        normalize = v/norm(v)
    end function normalize

    ! dot_prod: calculates the dot product between two 3D vectors
    pure real(dp) function dot_prod(v1, v2)
        type(vector3d), intent(in) :: v1, v2
        dot_prod = v1%x*v2%x + v1%y*v2%y + v1%z*v2%z
    end function dot_prod

    ! angle: calculates the angle between two 3D vectors (result in radians)
    pure real(dp) function angle(a,b)
        type(vector3d), intent(in)  :: a, b
        angle = acos(dot_prod(a,b)/(norm(a)*norm(b)))
    end function angle

    
    ! cross_product: calculates the cross product between two 3D vectors
    pure type(vector3d) function cross_product(v1, v2)
        type(vector3d), intent(in) :: v1, v2
        cross_product = vector3d(v1%y*v2%z - v1%z*v2%y, v1%z*v2%x - v1%x*v2%z, v1%x*v2%y - v1%y*v2%x)
    end function cross_product

end module geometry