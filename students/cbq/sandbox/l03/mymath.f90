module mymath
    implicit none
    real , parameter :: pi = 3.14
    
    contains
    real function circle_area ( a )
    real , intent ( in ) :: a
    circle_area = pi*  a **2 
    end function circle_area
end module mymath