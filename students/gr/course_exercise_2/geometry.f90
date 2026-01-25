module geometry
    use iso_fortran_env, only: real64
    implicit none

    ! Defining types
    type vector3d
        real(real64) :: x, y, z
    end type vector3d

    type point3d
        real(real64) :: x, y, z
    end type point3d

    ! Custom operators
    interface operator (+)
        module procedure sumpr, sumrp, sumpp, sumvp , sumpv, sumvv
    end interface

    interface operator (-)
        module procedure subpr, subrp, subpp, subvp , subpv, subvv
    end interface

    interface operator (*)
        module procedure mulrp, mulpr, mulrv , mulvr, mulvv
    end interface

    interface operator (/)
        module procedure divpr, divvr
    end interface

    contains

        ! sum functions
        elemental type(point3d) function sumpr(a, b)
            type(point3d), intent(in) :: a
            real(real64), intent(in) :: b
            sumpr = point3d(a%x + b, a%y + b, a%z + b)
        end function sumpr

        elemental type(point3d) function sumrp(a, b)
            real(real64), intent(in) :: a
            type(point3d), intent(in) :: b
            sumrp = point3d(a + b%x, a + b%y, a + b%z)
        end function sumrp
        
        elemental type(point3d) function sumpp(a, b)
            type(point3d), intent(in) :: a
            type(point3d), intent(in) :: b
            sumpp = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumpp

        elemental type(point3d) function sumvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            sumvp = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvp

        elemental type(point3d) function sumpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            sumpv = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumpv

        elemental type(vector3d) function sumvv(a, b)
            type(vector3d), intent(in) :: a, b
            sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvv

        ! subtract functions
        elemental type(point3d) function subpr(a, b)
            type(point3d), intent(in) :: a
            real(real64), intent(in) :: b
            subpr = point3d(a%x - b, a%y - b, a%z - b)
        end function subpr

        elemental type(point3d) function subrp(a, b)
            real(real64), intent(in) :: a
            type(point3d), intent(in) :: b
            subrp = point3d(a - b%x, a - b%y, a - b%z)
        end function subrp

        elemental type(point3d) function subvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            subvp = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvp

        elemental type(point3d) function subpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            subpv = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpv

        elemental type(vector3d) function subvv(a, b)
            type(vector3d), intent(in) :: a, b
            subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvv

        elemental type(vector3d) function subpp(a, b)
            type(point3d), intent(in) :: a, b
            subpp = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpp

        ! multiplication functions
        elemental type(point3d) function mulrp(a, b)
            real(real64), intent(in) :: a
            type(point3d), intent(in) :: b
            mulrp = point3d ( a * b%x, a * b%y, a * b%z)
        end function mulrp

        elemental type(point3d) function mulpr(a, b)
            type(point3d), intent(in) :: a
            real(real64), intent(in) :: b
            mulpr = point3d ( a%x * b, a%y * b, a%z * b)
        end function mulpr

        elemental type(vector3d) function mulrv(a, b)
            real(real64), intent(in) :: a
            type(vector3d), intent(in) :: b
            mulrv = vector3d ( a * b%x, a * b%y, a * b%z)
        end function mulrv

        elemental type(vector3d) function mulvr(a, b)
            type(vector3d), intent(in) :: a
            real(real64), intent(in) :: b
            mulvr = vector3d ( a%x * b, a%y * b, a%z * b)
        end function mulvr

        elemental type(vector3d) function mulvv(a, b)
            type(vector3d), intent(in) :: a, b
            mulvv = vector3d ( a%x * b%x, a%y * b%y, a%z * b%z)
        end function mulvv

        ! division function
        elemental type(point3d) function divpr(a, b)
            type(point3d), intent(in) :: a
            real(real64), intent(in) :: b
            divpr = point3d ( a%x / b, a%y / b, a%z / b)
        end function divpr

        elemental type(vector3d) function divvr(a, b)
            type(vector3d), intent(in) :: a
            real(real64), intent(in) :: b
            divvr = vector3d ( a%x / b, a%y / b, a%z / b)
        end function divvr

        ! function distance that calculates the distance between two 3d points
        pure real(real64) function distance(a, b)
            type(point3d), intent(in) :: a, b
            distance = sqrt( (b%x - a%x)**2 + (b%y - a%y)**2 + (b%z - a%z)**2 )
        end function distance

        ! function angle that calculates the angle between two vectors a and b (in radians)
        pure real(real64) function angle(a, b)
            type(vector3d), intent(in) :: a, b
            real(real64) :: norm_a, norm_b
            norm_a = sqrt( a%x**2 + a%y**2 + a%z**2 )
            norm_b = sqrt( b%x**2 + b%y**2 + b%z**2 )
            angle = acos( (a%x * b%x + a%y * b%y + a%z * b%z) / (norm_a * norm_b) )
        end function angle

        ! function normalize that takes a vector a and returns it divided by its length
        pure type(vector3d) function normalize(a)
            type(vector3d), intent(in) :: a
            real(real64) :: norm_a
            norm_a = sqrt( a%x**2 + a%y**2 + a%z**2 )
            normalize = a / norm_a
        end function normalize

        ! function cross product that takes two 3d vectors (a and b) and returns their cross product
        pure type(vector3d) function cross_product(a, b)
            type(vector3d), intent(in) :: a, b
            cross_product%x = a%y*b%z - a%z*b%y
            cross_product%y = a%z*b%x - a%x*b%z
            cross_product%z = a%x*b%y - a%y*b%x
        end function cross_product

end module geometry


