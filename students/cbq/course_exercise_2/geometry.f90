module geometry
    use iso_fortran_env, only: real64
    implicit none
    type :: vector3d
        real (real64) :: x , y, z  ! 64-bit real
    end type vector3d

    type :: point3d
    real(real64) :: x, y, z  ! 64-bit real
  end type point3d

    interface operator (+)
        module procedure sumvp , sumpv, sumvv
    end interface

    interface operator (-)
        module procedure subvp, subpv, subvv, subpp
    end interface

    interface operator (*)
        module procedure mulrv, mulvr
    end interface

    interface operator (/)
        module procedure divvr
    end interface

    contains

    ! Vector + point = point
    pure type(point3d) function sumvp(v, p)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
    end function sumvp

    ! Point + vector = point
    pure type(point3d) function sumpv(p, v)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
    end function sumpv

    ! Vector + vector = vector
    pure type(vector3d) function sumvv(a, b)
        type(vector3d), intent(in) :: a, b
        sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
    end function sumvv

    ! Vector - point = point
    pure type(point3d) function subvp(v, p)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z)
    end function subvp

    ! Point - vector = point
    pure type(point3d) function subpv(p, v)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        subpv = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
    end function subpv

    ! Point - point = vector
    pure type(vector3d) function subpp(p1, p2)
        type(point3d), intent(in) :: p1, p2
        subpp = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
    end function subpp

    ! Vector - vector = vector
    pure type(vector3d) function subvv(a, b)
        type(vector3d), intent(in) :: a, b
        subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
    end function subvv

    ! Real * vector = vector
    pure type(vector3d) function mulrv(r, v)
        real(real64), intent(in) :: r
        type(vector3d), intent(in) :: v
        mulrv = vector3d(r * v%x, r * v%y, r * v%z)
    end function mulrv

    ! Vector * real = vector
    pure type(vector3d) function mulvr(v, r)
        type(vector3d), intent(in) :: v
        real(real64), intent(in) :: r
        mulvr = vector3d(v%x * r, v%y * r, v%z * r)
    end function mulvr

    ! Vector / real = Vector
    pure type(vector3d) function divvr(v, r)
        type(vector3d), intent(in) :: v
        real(real64), intent(in) :: r
        divvr = vector3d(v%x / r, v%y / r, v%z / r)
    end function divvr


    ! Distance between two points
    pure real(real64) function distance(p1, p2)
        type(point3d), intent(in) :: p1, p2
        real(real64) :: dx, dy, dz
        dx = p1%x - p2%x
        dy = p1%y - p2%y
        dz = p1%z - p2%z
        distance = sqrt(dx*dx + dy*dy + dz*dz)
    end function distance

    ! Angle between two vectors (radians)
    pure real(real64) function angle(a, b)
        type(vector3d), intent(in) :: a, b
        real(real64) :: dot_product, mag_a, mag_b, cos_angle
        
        dot_product = a%x * b%x + a%y * b%y + a%z * b%z
        mag_a = sqrt(a%x*a%x + a%y*a%y + a%z*a%z)
        mag_b = sqrt(b%x*b%x + b%y*b%y + b%z*b%z)
        cos_angle = dot_product / (mag_a * mag_b)
        angle = acos(cos_angle)
    end function angle

    ! Normalize a vector (unit vector)
    pure type(vector3d) function normalize(a)
        type(vector3d), intent(in) :: a
        real(real64) :: mag
        
        mag = sqrt(a%x*a%x + a%y*a%y + a%z*a%z)
        if (mag == 0.0_real64) then
            normalize = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
        else
            normalize = vector3d(a%x/mag, a%y/mag, a%z/mag)
        end if
    end function normalize


    ! Cross-product
    pure type(vector3d) function xproduct(a, b)
        type(vector3d), intent(in) :: a, b
        xproduct= vector3d(a%y*b%z - a%z*b%y,  a%z* b%x - a%x*b%z, a%x*b%y - a%y*b%x )
    end function xproduct

end module geometry