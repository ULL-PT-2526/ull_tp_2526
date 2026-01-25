module geometry
    implicit none

    ! We use double precision for real numbers
    integer, parameter :: dp = selected_real_kind(15)


    type :: vector3d
        real(dp) :: x,y,z
    end type vector3d

    type :: point3d
        real(dp) :: x,y,z
    end type point3d

   ! We link the functions with operators

    interface operator(+)
        module procedure sumvp, sumpv, sumvv, sumpp
    end interface

    interface operator(-)
        module procedure subpv,subvp, subvv, subpp
    end interface

    interface operator(*)
        module procedure mulrv, mulvr, mulpr, mulrp
    end interface

    interface operator(/)
        module procedure divvr, divpr
    end interface 

 contains


    
    pure type(point3d) function sumvp(v, p) ! The result is a point
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z) 
    end function sumvp

! The sum of a point and a vector is conmutative sumvp=sumpv

    pure type(point3d) function sumpv(p, v) 
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        sumpv = sumvp(v,p) ! We reuse the function
    end function sumpv

    pure type(point3d) function sumpp(p1, p2)
        type(point3d), intent(in) :: p1
        type(point3d), intent(in) :: p2
        sumpp = point3d(p1%x + p2%x, p1%y + p2%y, p1%z + p2%z) 
    end function sumpp

    pure type(point3d) function subvp(v,p)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z) 
    end function subvp

    pure type(vector3d) function sumvv(v1,v2)
        type(vector3d), intent(in) :: v1
        type(vector3d), intent(in) :: v2
        sumvv = vector3d(v1%x + v2%x, v1%y + v2%y, v1%z + v2%z)
    end function sumvv


    pure type(vector3d) function subpp(p1, p2) 
        type(point3d), intent(in) :: p1
        type(point3d), intent(in) :: p2
        subpp = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z) 
    end function subpp ! to compute rji in main


    pure type(vector3d) function subvv(v1,v2)
        type(vector3d), intent(in) :: v1
        type(vector3d), intent(in) :: v2
        subvv = vector3d(v1%x - v2%x, v1%y - v2%y, v1%z - v2%z)
    end function subvv

    pure type(point3d) function subpv(p,v)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        subpv = subvp(v,p) ! Here we can reuse functions too
        subpv = point3d(-subpv%x, -subpv%y, -subpv%z)
    end function subpv



    pure type(vector3d) function mulrv(a, v)
        real(dp), intent(in) :: a
        type(vector3d), intent(in) :: v
        mulrv = vector3d(a*v%x, a*v%y, a*v%z) ! The output is a 3D vector
    end function mulrv


    ! The factor order do not affect the product in this case

    pure type(vector3d) function mulvr(v, a)
        type(vector3d), intent(in) :: v
        real(dp), intent(in) :: a
        mulvr = mulrv(a,v) 
    end function mulvr

    pure type(point3d) function mulrp(a, p)
        real(dp), intent(in) :: a
        type(point3d), intent(in) :: p
        mulrp = point3d(a*p%x, a*p%y, a*p%z) 
    end function mulrp

    pure type(point3d) function mulpr(p, a)
        type(point3d), intent(in) :: p
        real(dp), intent(in) :: a
        mulpr = mulrp(a,p) 
    end function mulpr

    pure type(vector3d) function divvr(v, a)
        type(vector3d), intent(in) :: v
        real(dp), intent(in) :: a
        divvr = vector3d(v%x/a, v%y/a, v%z/a) 
    end function divvr

    pure type(point3d) function divpr(p, a)
        type(point3d), intent(in) :: p
        real(dp), intent(in) :: a
        divpr = point3d(p%x/a, p%y/a, p%z/a) 
    end function divpr

    pure real(dp) function distance(p1,p2)
        type(point3d), intent(in) :: p1, p2
        distance = real(sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z)**2), dp)
    end function distance


    ! Fortran trigonometrics like sin(), works in radians

    pure real(dp)function angle(v1,v2)
        type(vector3d), intent(in) :: v1, v2
        ! We make use of the dot product
        ! We calculate the magnitudes of the vectors
        real(dp) :: mag1, mag2, prod
        mag1 = sqrt(v1%x**2 + v1%y**2 + v1%z**2)
        mag2 = sqrt(v2%x**2 + v2%y**2 + v2%z**2)
        prod = v1%x*v2%x + v1%y*v2%y + v1%z*v2%z
        angle = real(acos(prod/(mag1*mag2)), dp)
    end function angle


    pure type(vector3d) function normalize(v)
        type(vector3d), intent(in) :: v
        real(dp) :: mag
        mag = sqrt(v%x**2 + v%y**2 + v%z**2)
        normalize = divvr(v, mag) ! We reuse the division function
    end function normalize

    pure type(vector3d) function cross_product(v1,v2)
        type(vector3d), intent(in) :: v1, v2
        cross_product = vector3d(v1%y*v2%z - v1%z*v2%y, v1%z*v2%x - v1%x*v2%z, v1%x*v2%y - v1%y*v2%x)
    end function cross_product

end module geometry