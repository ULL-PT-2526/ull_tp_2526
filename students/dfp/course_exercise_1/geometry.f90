module geometry
    implicit none
    private
    public :: vector3d, point3d, dp
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: distance, angle, normalize, crossv
    integer, parameter :: dp = selected_real_kind(15, 307)

    type :: vector3d
        real(kind=dp) :: x, y, z
    end type vector3d

    type :: point3d
        real(kind=dp) :: x, y, z
    end type point3d

    interface operator(+)
        module procedure sumvv, sumpp, sumpv, sumvp
    end interface
    interface operator(-)
        module procedure subvv, subpp, subpv, subvp
    end interface
    interface operator(*)
        module procedure mulrv, mulvr
    end interface
    interface operator(/)
        module procedure divvr
    end interface

contains

    ! vector + vector
    function sumvv(a, b) result(r)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: r
        r%x = a%x + b%x
        r%y = a%y + b%y
        r%z = a%z + b%z
    end function

    ! point + point
    function sumpp(a, b) result(r)
        type(point3d), intent(in) :: a, b
        type(point3d) :: r
        r%x = a%x + b%x
        r%y = a%y + b%y
        r%z = a%z + b%z
    end function

    ! point + vector
    function sumpv(a, b) result(r)
        type(point3d), intent(in) :: a
        type(vector3d), intent(in) :: b
        type(point3d) :: r
        r%x = a%x + b%x
        r%y = a%y + b%y
        r%z = a%z + b%z
    end function

    ! vector + point
    function sumvp(a, b) result(r)
        type(vector3d), intent(in) :: a
        type(point3d), intent(in) :: b
        type(point3d) :: r
        r%x = a%x + b%x
        r%y = a%y + b%y
        r%z = a%z + b%z
    end function

    ! vector - vector
    function subvv(a, b) result(r)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: r
        r%x = a%x - b%x
        r%y = a%y - b%y
        r%z = a%z - b%z
    end function

    ! point - point
    function subpp(a, b) result(r)
        type(point3d), intent(in) :: a, b
        type(vector3d) :: r
        r%x = a%x - b%x
        r%y = a%y - b%y
        r%z = a%z - b%z
    end function

    ! point - vector
    function subpv(a, b) result(r)
        type(point3d), intent(in) :: a
        type(vector3d), intent(in) :: b
        type(point3d) :: r
        r%x = a%x - b%x
        r%y = a%y - b%y
        r%z = a%z - b%z
    end function

    ! vector - point
    function subvp(a, b) result(r)
        type(vector3d), intent(in) :: a
        type(point3d), intent(in) :: b
        type(vector3d) :: r
        r%x = a%x - b%x
        r%y = a%y - b%y
        r%z = a%z - b%z
    end function

    ! real * vector
    function mulrv(a, b) result(r)
        real(kind=dp), intent(in) :: a
        type(vector3d), intent(in) :: b
        type(vector3d) :: r
        r%x = a * b%x
        r%y = a * b%y
        r%z = a * b%z
    end function

    ! vector * real
    function mulvr(a, b) result(r)
        type(vector3d), intent(in) :: a
        real(kind=dp), intent(in) :: b
        type(vector3d) :: r
        r%x = a%x * b
        r%y = a%y * b
        r%z = a%z * b
    end function

    ! vector / real
    function divvr(a, b) result(r)
        type(vector3d), intent(in) :: a
        real(kind=dp), intent(in) :: b
        type(vector3d) :: r
        r%x = a%x / b
        r%y = a%y / b
        r%z = a%z / b
    end function

    ! distance between two points
    function distance(a, b) result(d)
        type(point3d), intent(in) :: a, b
        real(kind=dp) :: d
        real(kind=dp) :: dx, dy, dz
        dx = a%x - b%x
        dy = a%y - b%y
        dz = a%z - b%z
        d = sqrt(dx*dx + dy*dy + dz*dz)
    end function

    ! dot product of two vectors (internal helper)
    pure function dotv(a, b) result(d)
        type(vector3d), intent(in) :: a, b
        real(kind=dp) :: d
        d = a%x*b%x + a%y*b%y + a%z*b%z
    end function

    ! length (norm) of a vector (internal helper)
    pure function normv(a) result(n)
        type(vector3d), intent(in) :: a
        real(kind=dp) :: n
        n = sqrt(dotv(a,a))
    end function

    ! angle between two vectors in radians
    function angle(a, b) result(theta)
        type(vector3d), intent(in) :: a, b
        real(kind=dp) :: theta
    real(kind=dp) :: na, nb, quot
        na = normv(a)
        nb = normv(b)
        if (na <= 0.0_dp .or. nb <= 0.0_dp) then
            theta = 0.0_dp
        else
            quot = dotv(a,b) / (na*nb)
            ! Clamp to [-1,1] to avoid domain errors
            if (quot > 1.0_dp) quot = 1.0_dp
            if (quot < -1.0_dp) quot = -1.0_dp
            theta = acos(quot)
        end if
    end function

    ! normalize: return vector divided by its length (zero vector if length==0)
    function normalize(a) result(r)
        type(vector3d), intent(in) :: a
        type(vector3d) :: r
        real(kind=dp) :: n
        n = normv(a)
        if (n <= 0.0_dp) then
            r%x = 0.0_dp; r%y = 0.0_dp; r%z = 0.0_dp
        else
            r = a / n
        end if
    end function

    ! cross product of two vectors
    function crossv(a, b) result(r)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: r
        r%x = a%y*b%z - a%z*b%y
        r%y = a%z*b%x - a%x*b%z
        r%z = a%x*b%y - a%y*b%x
    end function

end module geometry