module geometry
    implicit none

    ! Making private internal operators and functions # Lecture 4
    private

    ! Needed types are public # Lecture 4
    public :: dp, vector3d, point3d
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: angle, normalize, cross_product, distance

    
        !--- Parameters and types ---!

        ! Defining the precision of the values to 64 bits # Lecture 3
        integer, parameter :: dp = selected_real_kind(p=15)

        ! Defining vector3d and point3d types # Lecture 5
        type vector3d
            real (kind=dp) :: x,y,z
        end type vector3d

        type point3d
            real (kind=dp) :: x,y,z
        end type point3d


        !--- Custom operators ---! 

        ! Defining some custom operators # Lecture 5

        ! Adding operator
        interface operator(+)
            module procedure sumvp, sumpv, sumvv
        end interface

        ! Substracting operator
        interface operator(-)
            module procedure subvp, subpv, subvv, subpp
        end interface

        ! Product operator
        interface operator(*)
            module procedure mulrv, mulvr
        end interface

        ! Dividing operator
        interface operator(/)
            module procedure divvr
        end interface

    contains

        !--- Functions for arimethic operators ---!

        ! Function sumvp for adding vectors and points
        function sumvp(v,p) result(res)

            ! Defining v and p as special type
            ! with three components each
            type(vector3d), intent(in) :: v
            type(point3d), intent(in) :: p
            type(point3d) :: res

            ! Accesing p and v components by using % # Lecture 5
            res = point3d(v%x + p%x, v%y + p%y, v%z + p%z)

        end function sumvp


        ! Function sumvp for adding points and vectors
        function sumpv(p,v) result(res)

            type(point3d), intent(in) :: p
            type(vector3d), intent(in) :: v
            type(point3d) :: res

            res = point3d(p%x + v%x, p%y + v%y, p%z + v%z)

        end function sumpv

        ! Function to add two vectors
        function sumvv(v1, v2) result(res)

            type(vector3d), intent(in) :: v1, v2
            type(vector3d) :: res

            res = vector3d(v1%x + v2%x, v1%y + v2%y, v1%z + v2%z)

        end function sumvv


        ! Function subvp for substracting vectors and points
        function subvp(v,p) result(res)

            type(vector3d), intent(in) :: v
            type(point3d), intent(in) :: p
            type(point3d) :: res

            res = point3d(v%x - p%x, v%y - p%y, v%z - p%z)

        end function subvp


        ! Function subpv for subtstracting points and vectors
        function subpv(p,v) result(res)   
            
            type(point3d), intent(in) :: p
            type(vector3d), intent(in) :: v
            type(point3d) :: res

            res = point3d(p%x - v%x, p%y - v%y, p%z - v%z)

        end function subpv

        ! Function to subtract two vectors
        function subvv(v1, v2) result(res)

            type(vector3d), intent(in) :: v1, v2
            type(vector3d) :: res

            res = vector3d(v1%x - v2%x, v1%y - v2%y, v1%z - v2%z)
            
        end function subvv

        ! Function to subtract two points
        function subpp(p1, p2) result(res)

            type(point3d), intent(in) :: p1, p2
            type(vector3d) :: res

            res = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
            
        end function subpp

        ! Function mulrv for multiplying reals and vectors
        function mulrv(r,v) result(res)

            real(kind=dp),  intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: res

            res = vector3d(r * v%x, r * v%y, r * v%z)

        end function mulrv


        ! Function mulvr for multiplying vectors and reals
        function mulvr(v,r) result(res)   
            
            type(vector3d), intent(in) :: v
            real(kind=dp),  intent(in) :: r
            type(vector3d) :: res

            res = vector3d(v%x * r, v%y * r, v%z * r)

        end function mulvr


        ! Function divvr for dividing vectors with reals
        function divvr(v,r) result(res)
            
            type(vector3d), intent(in) :: v
            real(kind=dp),  intent(in) :: r
            type(vector3d) :: res

            res = vector3d(v%x / r, v%y / r, v%z / r)

        end function divvr


        !--- Functions for geometric operators ---!

        ! Function distance between two points
        ! Computing by vectors module
        function distance(p1,p2) result(d)   

            type(point3d), intent(in) :: p1, p2
            real(kind=dp) :: d

            d = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z)**2)

        end function distance


        ! Function angle between two vectors
        ! Computing by dot product
        function angle(v1,v2) result(ang)   
            
            type(vector3d), intent(in) :: v1, v2
            real(kind=dp) :: ang
            real(kind=dp) :: dot_prod, mod_v1, mod_v2
            
            dot_prod = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
            mod_v1 = sqrt(v1%x**2 + v1%y**2 + v1%z**2)
            mod_v2 = sqrt(v2%x**2 + v2%y**2 + v2%z**2)
            
            ang = acos(dot_prod / (mod_v1 * mod_v2))

        end function angle  


        ! Function normalize a vector
        function normalize(v) result(res)
            
            type(vector3d), intent(in) :: v
            type(vector3d) :: res
            real(kind=dp) :: mod
        
            mod = sqrt(v%x**2 + v%y**2 + v%z**2)
            res = vector3d(v%x / mod, v%y / mod, v%z / mod)   

        end function normalize 

        ! Function cross product between two vectors
        function cross_product(v1,v2) result(res) 
            
            type(vector3d), intent(in) :: v1, v2
            type(vector3d) :: res
        
            res%x = v1%y * v2%z - v1%z * v2%y
            res%y = v1%z * v2%x - v1%x * v2%z
            res%z = v1%x * v2%y - v1%y * v2%x      

        end function cross_product 

end module geometry