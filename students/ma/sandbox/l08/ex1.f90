program mult_omp
    !$ use omp_lib 
    ! To avoid compatibility issues
    implicit none

    ! We take the exercise of matrix in lecture l03
    ! We declare the matrices
    real, dimension(3,4) :: A
    real, dimension(3,3) :: B
    real, dimension(3,4) :: C ! The matrix result
    integer :: i,j,k, nt = 1, tid = 0 ! the counters for the loops

    ! We can initialize the matrix using the reshape function
    ! The order is to specify that fills the matrix by rows

    A = reshape([3.0, 2.0, 4.0, 1.0, &
             2.0, 4.0, 2.0, 2.0, &
             1.0, 2.0, 3.0, 7.0], shape=[3,4], order=[2,1])

    B =reshape([3.0, 2.0, 4.0, &
             2.0, 1.0, 2.0, &
             3.0, 0.0, 2.0], shape=[3,3], order=[2,1])
    
    
    ! We can't multiply A and B, but B and A yes

    ! We can multiply matrices using matmul or loops

    C = 0.0 ! We initialize the matrix to 0
    
    
    !$omp parallel private(nt,tid,i,j,k) ! with collapse we parallize better and private is to avoid conflicts with variables
    !$ nt = omp_get_num_threads()
    !$ tid = omp_get_thread_num()
    !$omp do collapse(3)
    ! Collapse unravel the nested loops
    do i=1,3
        do j=1,4
            do k=1,3
                C(i,j) = C(i,j) + B(i,k)*A(k,j)
                !print '(xA,i3,xA,i2,xA,i3)', "Do loop, i =", i, & ! If we want to look the threads in the running 
                !& "thread ", tid, " of ", nt
            end do
        end do
    end do

    !$omp end do
    !$omp end parallel
    print *, "The result of the matrix multiplication is", C


end program mult_omp