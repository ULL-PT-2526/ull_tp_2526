module geometry
    implicit none
    
    integer, parameter :: bit64 = SELECTED_REAL_KIND(p=15, r=307) ! Defining our desired precision (64 bits) allowing full portability, instead of just using "kind=8"
      
    type point3d
     real(kind=bit64) :: x,y,z ! Using our defined precision
    end type point3d
    type vector3d
     real(kind=bit64) :: x,y,z 
    end type vector3d   

    ! Definition of operators to simplify the code
    interface operator(+)
        module procedure sumvp, sumpv, sumvv
    end interface
    interface operator(-)
        module procedure subvp, subpv, subvv
    end interface
    interface operator(*)
        module procedure mulrv, mulvr
    end interface
    interface operator(/)
        module procedure divvr
    end interface

contains
  
  ! In all of the functions below we use the types "point3d" and "vector3d" as points and vectors with 3 components, making the code easier to read (and write).
  
  ! Sum Vector + Point
  pure function sumvp(v, p) result(res)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in)  :: p
    type(point3d) :: res

    res%x = p%x + v%x
    res%y = p%y + v%y
    res%z = p%z + v%z
  end function sumvp

  ! Sum Point + Vector (same result, different input order)
  pure function sumpv(p, v) result(res)
    type(point3d), intent(in)  :: p
    type(vector3d), intent(in) :: v
    type(point3d) :: res

    res%x = p%x + v%x
    res%y = p%y + v%y
    res%z = p%z + v%z
  end function sumpv

  ! Subtract Vector - Point
  pure function subvp(v, p) result(res)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in)  :: p
    type(vector3d) :: res

    res%x = v%x - p%x
    res%y = v%y - p%y
    res%z = v%z - p%z
  end function subvp

  ! Substract Point - Vector
  pure function subpv(p, v) result(res)
    type(point3d), intent(in)  :: p
    type(vector3d), intent(in) :: v
    type(point3d) :: res

    res%x = p%x - v%x
    res%y = p%y - v%y
    res%z = p%z - v%z
  end function subpv

  ! Product of real number and vector
  pure function mulrv(r, v) result(res)
    real(kind=bit64), intent(in)  :: r
    type(vector3d), intent(in) :: v
    type(vector3d) :: res
    res%x = r*v%x
    res%y = r*v%y
    res%z = r*v%z
  end function mulrv
  
! Product of vector and real number (same result, different input order)
  pure function mulvr(v, r) result(res)
    real(kind=bit64), intent(in)  :: r
    type(vector3d), intent(in) :: v
    type(vector3d) :: res
    res%x = r*v%x
    res%y = r*v%y
    res%z = r*v%z
  end function mulvr

! Division of vector and real
  pure function divvr(v, r) result(res)
    real(kind=bit64), intent(in)  :: r
    type(vector3d), intent(in) :: v
    type(vector3d) :: res
    res%x = v%x / r
    res%y = v%y / r
    res%z = v%z / r
  end function divvr
  
  ! Vector + Vector (useful for later)
  pure function sumvv(a, b) result(res)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: res
    res%x = a%x + b%x
    res%y = a%y + b%y
    res%z = a%z + b%z
  end function sumvv

  ! Vector - Vector (useful for later)
  pure function subvv(a, b) result(res)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: res
    res%x = a%x - b%x
    res%y = a%y - b%y
    res%z = a%z - b%z
  end function subvv

  
  ! Distance between two points
  pure function distance(p1, p2) result(d) ! distance in same units as given
    type(point3d), intent(in) :: p1, p2
    real(kind=bit64) :: d 

    d = sqrt((p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2)
  end function distance
  
  !Two more useful functions are mulvv and norm
  pure function mulvv(a,b) result(res) ! scalar product of two vectors
    type(vector3d), intent(in) :: a,b
    real(kind=bit64) :: res
    res = a%x*b%x + a%y*b%y + a%z*b%z
  end function mulvv

  pure function norm(v) result(res)
    type(vector3d), intent(in) :: v
    real(kind=bit64) :: res
    res = sqrt(mulvv(v,v))
  end function norm
  
  ! Angle (radians) between two vectors
  pure function angle(a,b) result(theta) ! in radians
    type(vector3d), intent(in) :: a,b
    real(kind=bit64) :: theta
    if (norm(a)*norm(b) /= 0.0d0) then
        theta = acos( mulvv(a,b) / (norm(a)*norm(b)) ) !using the formula of the angle between two vectors
    else
        theta = 0.0d0
    end if
  end function angle
  
  ! Normalization of a vector
  pure function normalize(v) result(v_norm)
    type(vector3d), intent(in) :: v
    type(vector3d) :: v_norm
    if (norm(v) /= 0.0d0) then
        v_norm = divvr(v, norm(v)) ! using the divvr function with the "/" operator
    else
        v_norm = v  ! or could be zero vector
    end if
  end function normalize
  
  ! Cross product between two vectors
  pure function cross_product(a,b) result(c)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: c
    c%x = a%y*b%z - a%z*b%y
    c%y = a%z*b%x - a%x*b%z
    c%z = a%x*b%y - a%y*b%x
  end function cross_product

end module geometry
