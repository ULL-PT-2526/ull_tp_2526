module geometry


  use iso_fortran_env, only: real64
  implicit none


  !=======================
  ! Type definitions
  !=======================
  type :: vector3d
    real(real64) :: x, y, z
  end type vector3d


  type :: point3d
    real(real64) :: x, y, z
  end type point3d


  !=======================
  ! Operator interfaces 
  !=======================
  interface operator(+)
    module procedure :: sumvp, sumpv, sumvv, sumpp
  end interface


  interface operator(-)
    module procedure :: subvp, subpv, subvv
  end interface


  interface operator(*)
    module procedure :: mulrv, mulvr, mulrp, mulpr
  end interface


  interface operator(/)
    module procedure :: divvr, divpr
  end interface


contains


  !=======================
  ! Function implementations 
  !=======================

  ! === Addition: vector + point ===
  function sumvp(v, p) result(r)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in)  :: p
    type(point3d) :: r
    r = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  ! === Addition: point + vector ===
  function sumpv(p, v) result(r)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    type(point3d) :: r
    r = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

  ! === Addition: point + point ===
  function sumpp(p1, p2) result(r)
    type(point3d), intent(in) :: p1, p2
    type(point3d) :: r
    r = point3d(p1%x + p2%x, p1%y + p2%y, p1%z + p2%z)
  end function sumpp
  
  ! === Subtraction: point - point ===
  function subvp(p1, p2) result(r)
    type(point3d), intent(in) :: p1, p2
    type(vector3d) :: r
    r = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
  end function subvp

  ! === Subtraction: point - vector ===
  function subpv(p, v) result(r)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    type(point3d) :: r
    r = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
  end function subpv

  ! === Addition: vector + vector ===
  function sumvv(a, b) result(r)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: r
    r = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumvv

  ! === Subtraction: vector - vector ===
  function subvv(a, b) result(r)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: r
    r = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
  end function subvv

  ! === Multiplication: scalar * vector ===
  function mulrv(r_in, v) result(res)
    real(real64), intent(in) :: r_in
    type(vector3d), intent(in) :: v
    type(vector3d) :: res
    res = vector3d(r_in * v%x, r_in * v%y, r_in * v%z)
  end function mulrv

  ! === Multiplication: vector * scalar ===
  function mulvr(v, r_in) result(res)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r_in
    type(vector3d) :: res
    res = vector3d(v%x * r_in, v%y * r_in, v%z * r_in)
  end function mulvr

  ! === Multiplication: scalar * point ===
  function mulrp(r_in, p) result(res)
    real(real64), intent(in) :: r_in
    type(point3d), intent(in) :: p
    type(point3d) :: res
    res = point3d(r_in * p%x, r_in * p%y, r_in * p%z)
  end function mulrp

  ! === Multiplication: point * scalar ===
  function mulpr(p, r_in) result(res)
    type(point3d), intent(in) :: p
    real(real64), intent(in) :: r_in
    type(point3d) :: res
    res = point3d(p%x * r_in, p%y * r_in, p%z * r_in)
  end function mulpr

  ! === Division: vector / scalar ===
  function divvr(v, r_in) result(res)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r_in
    type(vector3d) :: res
    res = vector3d(v%x / r_in, v%y / r_in, v%z / r_in)
  end function divvr

  ! === Division: point / scalar ===
  function divpr(p, r_in) result(res)
    type(point3d), intent(in) :: p
    real(real64), intent(in) :: r_in
    type(point3d) :: res
    res = point3d(p%x / r_in, p%y / r_in, p%z / r_in)
  end function divpr

  ! === Distance between two points ===
  function distance(p1, p2) result(d)
    type(point3d), intent(in) :: p1, p2
    real(real64) :: d
    d = sqrt((p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2)
  end function distance

  ! === Squared distance between two points ===
  function distance2(p1, p2) result(d2)
    type(point3d), intent(in) :: p1, p2
    real(real64) :: d2
    d2 = (p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2
  end function distance2

  ! === Norm (magnitude) of a vector ===
  function norm(a) result(n)
    type(vector3d), intent(in) :: a
    real(real64) :: n
    n = sqrt(a%x**2 + a%y**2 + a%z**2)
  end function norm

  ! === Angle between two vectors ===
  function angle(a, b) result(theta)
    type(vector3d), intent(in) :: a, b
    real(real64) :: theta
    real(real64) :: dotprod, na, nb
    dotprod = a%x * b%x + a%y * b%y + a%z * b%z
    na = norm(a)
    nb = norm(b)
    theta = acos(dotprod / (na * nb))
  end function angle

  ! === Normalize a vector ===
  function normalize(a) result(res)
    type(vector3d), intent(in) :: a
    type(vector3d) :: res
    real(real64) :: n
    n = norm(a)
    if (n > 0.0) then
      res = a / n
    else
      res = vector3d(0.0, 0.0, 0.0)
    end if
  end function normalize

  ! === Cross product of two vectors ===
  function cross(a, b) result(c)
    type(vector3d), intent(in) :: a, b
    type(vector3d) :: c
    c = vector3d(a%y*b%z - a%z*b%y, a%z*b%x - a%x*b%z, a%x*b%y - a%y*b%x)
  end function cross


end module geometry
