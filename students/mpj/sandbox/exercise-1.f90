! Write a program that computes and prints the matrix multiplication of two real arrays
!
program matrix_multiplication
    integer :: A(3, 4)
    integer :: B(3, 3)
    integer :: subA(1, 3)     ! We cannot multiply A and B, due to the dimensions of the matrices
    integer :: R(size(subA, 1), size(B, 2))
    integer :: i, j, k, suma
    A = reshape ((/3, 2, 1, 2, 4, 2, 4, 2, 3, 1, 2, 7/), (/3, 4/))
    B = reshape ((/3, 2, 3, 2, 1, 0, 4, 2, 2/), (/3, 3/))
    subA = reshape(A(1,1:3), (/1,3/))
    do i=1, size(subA, 1)
        do k=1, size(B, 2)
            suma = 0
            do j=1, size(subA, 2)
                suma = suma + subA(i, j) * B(j, k)
            end do
            R(i, k) = suma
        end do
    end do
    print *, "Nueva matriz A:", R
end program matrix_multiplication