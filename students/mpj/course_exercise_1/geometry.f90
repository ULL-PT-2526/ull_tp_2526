! Exercise: Direct numerical integration of a system consisting of gravitationally interacting
! particles using the leapfrog integration method.

! Module with all the rutines that I will import after in other files with the command "use"
module geometry
    implicit none

    integer, parameter :: dp = SELECTED_REAL_KIND(p=15, r=307) ! Defining our desired precision (64 bits)

    ! we create vector3d, that is a type of data personalised
    type :: vector3d
        real(dp) :: x = 0.0_dp   ! 0 is the value for defect
        real(dp) :: y = 0.0_dp
        real(dp) :: z = 0.0_dp
    end type vector3d

    type :: point3d
        real(dp) :: x = 0.0_dp
        real(dp) :: y = 0.0_dp
        real(dp) :: z = 0.0_dp
    end type point3d

    ! Operators matching functions defined below 
    interface operator(+)
        module procedure sumvp, sumpv, sumvv
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

   ! Functions to add and substract points and vectors:
   ! vector + point
    function sumvp(vector, point) result(res)
        implicit none
        type(vector3d), intent(in) :: vector
        type(point3d),  intent(in) :: point

        type(vector3d) :: res
        res%x = vector%x + point%x
        res%y = vector%y + point%y
        res%z = vector%z + point%z
    end function sumvp

    ! point + vector
    function sumpv(point, vector) result(res)
        implicit none
        type(point3d),  intent(in) :: point
        type(vector3d), intent(in) :: vector

        type(point3d) :: res
        res%x = point%x + vector%x
        res%y = point%y + vector%y
        res%z = point%z + vector%z
    end function sumpv

    ! vector + vector = vector
    function sumvv(vector1, vector2) result(res)
        implicit none
        type(vector3d), intent(in) :: vector1
        type(vector3d), intent(in) :: vector2

        type(vector3d) :: res
        res%x = vector1%x + vector2%x
        res%y = vector1%y + vector2%y
        res%z = vector1%z + vector2%z
    end function sumvv

    ! vector - point
    function subvp(vector, point) result(res)
        implicit none
        type(vector3d), intent(in) :: vector
        type(point3d),  intent(in) :: point

        type(vector3d) :: res
        res%x = vector%x - point%x
        res%y = vector%y - point%y
        res%z = vector%z - point%z
    end function subvp

    ! point - vector
    function subpv(point, vector) result(res)
        implicit none
        type(point3d),  intent(in) :: point
        type(vector3d), intent(in) :: vector

        type(point3d) :: res
        res%x = point%x - vector%x
        res%y = point%y - vector%y
        res%z = point%z - vector%z
    end function subpv

    ! vector - vector = vector
    function subvv(vector1, vector2) result(res)
        implicit none
        type(vector3d), intent(in) :: vector1
        type(vector3d), intent(in) :: vector2

        type(vector3d) :: res
        res%x = vector1%x - vector2%x
        res%y = vector1%y - vector2%y
        res%z = vector1%z - vector2%z
    end function subvv

    ! point - point = vector
    function subpp(point1, point2) result(res)
        implicit none
        type(point3d), intent(in) :: point1
        type(point3d), intent(in) :: point2

        type(vector3d) :: res
        res%x = point1%x - point2%x
        res%y = point1%y - point2%y
        res%z = point1%z - point2%z
    end function subpp

    ! real * vector
    function mulrv(num, vector) result(res)
        implicit none
        real(dp), intent(in) :: num
        type(vector3d),  intent(in) :: vector

        type(vector3d) :: res
        res%x = num * vector%x
        res%y = num * vector%y
        res%z = num * vector%z
    end function mulrv

    ! vector * real
    function mulvr(vector, num) result(res)
        implicit none
        real(dp), intent(in) :: num
        type(vector3d),  intent(in) :: vector

        type(vector3d) :: res
        res%x = vector%x * num
        res%y = vector%y * num
        res%z = vector%z * num
    end function mulvr

    ! vector / real
    function divvr(vector, num) result(res)
        implicit none
        real(dp), intent(in) :: num
        type(vector3d),  intent(in) :: vector

        type(vector3d) :: res
        res%x = vector%x / num
        res%y = vector%y / num
        res%z = vector%z / num
    end function divvr

    ! Calculates the distance between two points
    function distance(point1, point2) result(res)
        implicit none
        type(point3d), intent(in) :: point1
        type(point3d), intent(in) :: point2
        
        real(dp) :: res
        res = sqrt( (point2%x - point1%x)**2 + &
                    (point2%y - point1%y)**2 + &
                    (point2%z - point1%z)**2 )

    end function distance


    ! Angle between two vectors
    function angle(vector1, vector2) result(res)
        type(vector3d), intent(in) :: vector1, vector2
        real(dp) :: res
        real(dp) :: scalar_prod, mod1, mod2, cos_angle

        scalar_prod = vector1%x * vector2%x + &
                      vector1%y * vector2%y + &
                      vector1%z * vector2%z

        mod1 = sqrt(vector1%x**2 + vector1%y**2 + vector1%z**2)
        mod2 = sqrt(vector2%x**2 + vector2%y**2 + vector2%z**2)

        if (mod1 == 0.0_dp .or. mod2 == 0.0_dp) then
            res = 0.0_dp
        else
            cos_angle = scalar_prod / (mod1 * mod2)
            res = acos(cos_angle)
        end if

    end function angle

    ! Returns the unitary vector in the same direction
    function normalize(vector) result(res)
        type(vector3d), intent(in) :: vector
        type(vector3d) :: res
        real(dp) :: mod

        mod = sqrt(vector%x**2 + vector%y**2 + vector%z**2)
        res = vector / mod
    
    end function normalize

    ! Cross product
    function cross_product(vector1, vector2) result(res)
        type(vector3d), intent(in) :: vector1, vector2
        type(vector3d) :: res

        res%x = vector1%y * vector2%z - vector1%z * vector2%y
        res%y = -(vector1%x * vector2%z - vector1%z * vector2%x)
        res%z = vector1%x * vector2%y - vector1%y * vector2%x
    
    end function cross_product

end module geometry

