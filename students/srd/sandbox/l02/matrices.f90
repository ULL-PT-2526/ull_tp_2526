program matrices
    implicit none
    ! integer, dimension(3,4) :: A  ! 3 rows, 4 columns ! Better to not have them fixed size
    ! integer, dimension(3,3) :: B  ! 3 rows, 3 columns
    integer, allocatable :: A(:,:)
    integer, allocatable :: B(:,:)
    
    integer, allocatable :: C(:,:)  ! Fixed size matrix being C the multiplication of B and A
    
    integer :: i, j
    integer :: nrowsA, ncolsA, nrowsB, ncolsB ! Useful to check later the possible multiplications

    allocate(A(3,4)) ! Done like this the compiler doesn't check if all the matmuls are possible 
    allocate(B(3,3))
    
    A = reshape([3,2,4,1, &
                 2,4,2,2, &
                 1,2,3,7], shape(A))

    B = reshape([3,2,4, &
                 3,1,2, &
                 3,0,2], shape(B))

    
    ! Good way to print them is using a do loop over their elements
    
    print *, "Matrix A:"
    do i = 1, size(A,1)
        print *, (A(i,j), j=1,size(A,2))
    end do

    print *, "Matrix B:"
    do i = 1, size(B,1)
        print *, (B(i,j), j=1,size(B,2))
    end do

    nrowsA = size(A,1)
    ncolsA = size(A,2)
    nrowsB = size(B,1)
    ncolsB = size(B,2)
    
    ! A*B Possible? Number of columns in A must be equal to number of rows in B
    if (ncolsA == nrowsB) then
        print *, "A*B is possible, resulting size: ", nrowsA, "x", ncolsB
        allocate(C(nrowsA,ncolsB)) ! For any reason even if the condition is not satisfied the program execute these lines and gives error since matmul is not possible
        C = matmul(A,B)
        deallocate(C)
    else
        print *, "A*B is NOT possible with current dimensions."
    end if

    ! Same for B*A 
    
    if (ncolsB == nrowsA) then
        print *, "B*A is possible, resulting size: ", nrowsB, "x", ncolsA
        allocate(C(nrowsB,ncolsA))
        C = matmul(B,A)
        print *, "The resulting matrix C is:"
        do i = 1, 3
            print *, (C(i,j), j=1,4)
        end do
        deallocate(C)
    else
        print *, "B*A is NOT possible with current dimensions."
    end if
end program matrices

