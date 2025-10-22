program matrix_mult_parallel

    use omp_lib
    implicit none

    ! Matrix dimensions
    integer :: na, ma, nb, mb
    real :: A(3,4), B(3,3), C(3,4)
    integer :: i, j, k
    real :: tmp_sum ! Temporary private variable for summation

    ! Define matrix sizes
    na = 3
    ma = 4
    nb = 3
    mb = 3

    ! Define matrices A and B
    A = reshape((/3.0,2.0,1.0, &
                  2.0,4.0,2.0, &
                  4.0,2.0,3.0, &
                  1.0,2.0,7.0/), (/na,ma/))

    B = reshape((/3.0,2.0,3.0, &
                  2.0,1.0,0.0, &
                  4.0,2.0,2.0/), (/nb,mb/))

    ! Check if B*A is possible (columns of B must equal rows of A)
    if (mb /= na) then
        print *, "ERROR: Cannot multiply B*A. Incompatible dimensions."
        print *, "B is ", nb, "x", mb, " and A is ", na, "x", ma
        stop
    end if

    ! Parallelize the two outer loops that iterate over the elements of C.
    ! 'collapse(2)' treats these two loops as a single large iteration space.
    ! 'i', 'j', 'k', and 'tmp_sum' must be private for each thread.
    ! Loop variables 'i' and 'j' are private by default.
    ! Matrices A, B, C and the dimension variables are shared.

    !$omp parallel do collapse(2) private(k, tmp_sum) shared(A, B, C)
    do i = 1, nb
        do j = 1, ma
            tmp_sum = 0.0
            do k = 1, mb
                tmp_sum = tmp_sum + B(i,k) * A(k,j)
            end do
            C(i,j) = tmp_sum
        end do
    end do
    !$omp end parallel do

    ! Print the resulting matrix
    print *, "Matrix B*A ="
    do i = 1, nb
        print '(4F6.1)', C(i,1:ma)
    end do

end program matrix_mult_parallel