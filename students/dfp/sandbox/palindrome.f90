program palindrome
    implicit none
    integer :: n, n_dig, i, n_signed
    character(len=:), allocatable :: n_str ! dynamic string
    logical :: is_palindrome

    print *, "Enter a value for n"
    read *, n_signed

    if (n_signed < 0) then 
        n = abs(n_signed)
    end if

    if (n == 0) then
        n_dig = 1
    else
        n_dig = floor(log10(real(n))) + 1
    end if

    ! string with the exact len
    allocate(character(len=n_dig) :: n_str)

    ! convert number into a string
    write(n_str, '(I0)') n 

    ! verify if it is a palindrome
    is_palindrome = .true.
    do i = 1, n_dig/2
        if (n_str(i:i) /= n_str(n_dig-i+1:n_dig-i+1)) then
            print *, n_signed, 'is not a palindrome :('
            is_palindrome = .false.
            exit
        end if
    end do

    if (is_palindrome) then
        print *, n_signed, 'is a palindrome :)'
    end if

end program palindrome
