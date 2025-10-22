program factor
    implicit none
    integer :: n, i
    integer :: fact

    print *, "Insert the number of imprints of n:"
    read(*,*) n

    do i = 0, n
    if (i == 0) then
      fact = 1
    else
      fact = fact * i
    end if
        print * , i, " ! = " , fact 
    end do
end program factor

