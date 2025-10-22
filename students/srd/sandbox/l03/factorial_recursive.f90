program factorial_recursive
    implicit none
    integer :: n

    print *, "Introduce un número entero no negativo:"
    read *, n

    if (n < 0) then
        print *, "Error: el factorial no está definido para negativos."
    else if (n > 20) then
        ! 20! = 2.43e18 → upper limit that can be solved in 64 bits
        print *, "Error: el número es demasiado grande para calcular el factorial en 64 bits."
    else
        print *, "Factorial de", n, "es", factorial(n)
    end if

contains

    recursive function factorial(n) result(f)
        integer, intent(in) :: n
        integer(kind=8) :: f
        if (n <= 1) then
            f = 1
        else
            f = n * factorial(n - 1)
        end if
    end function factorial

end program factorial_recursive

