module mymath
    implicit none
    real, parameter :: pi =3.14159
contains
    real function circle_area(r)
        real, intent(in) :: r
        circle_area = pi * r**2
    end function circle_area
end module mymath
