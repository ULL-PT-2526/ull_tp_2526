program ex3
    implicit none
    integer , parameter :: n=3
    integer :: a(n,n)
    a = reshape ((/1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9/) , (/n,n/))
    print *, 'Element =', a(3 ,2)
    print *, 'Column 1 =', a(: ,1)
    print *, 'Subarray =', a (:2 ,:2)
    print *, 'Whole =', a
    print *, 'Transposed =', transpose (a)
end program ex3