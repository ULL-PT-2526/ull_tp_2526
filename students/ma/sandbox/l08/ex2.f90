program mult_seq

    use omp_lib 
    implicit none

    ! We take the exercise of matrix in lecture l03 and do the version of 200 elements
    ! We declare the matrices
    real, dimension(200,200) :: A, B, C
    integer :: i,j,k ! the counters for the loops
    real :: t1, t2, t3, t4

    ! We can initialize the matrix using the reshape function
    ! The order is to specify that fills the matrix by rows


  ! We initialize the random generator
    call random_seed()
    call random_number(A)
    call random_number(B)

    C = 0.0

    t1 = omp_get_wtime()
    ! Sequential multiplication
    do i = 1, 200
        do j = 1, 200
            do k = 1, 200
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do
    t2 = omp_get_wtime() ! In order to take the time of the procesess

    print *, 'Sequential time (s): ', t2 - t1

    C = 0.0

    t3 = omp_get_wtime()

    !$omp parallel
    !$omp do collapse(3) private(i,j,k)! with collapse we parallize better and private is to avoid conflicts with variables
    ! Collapse unravel the nested loops
     do i = 1, 200
        do j = 1, 200
            do k = 1, 200
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do

    !$omp end do
    !$omp end parallel
    t4 = omp_get_wtime()
    print *, 'Parallel time (s): ', t4 - t3

end program mult_seq