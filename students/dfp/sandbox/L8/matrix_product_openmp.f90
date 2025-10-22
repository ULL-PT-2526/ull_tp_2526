program matrix_product
    use omp_lib
    implicit none
    
    ! Matrix dimensions
    integer, parameter :: n = 3, m = 3, p = 3
    
    ! Matrix declarations
    real :: A(n, m), B(m, p), C(n, p)
    
    ! Loop indices
    integer :: i, j, k
    
    ! OpenMP variables
    integer :: num_threads, thread_id

    ! Initialize matrix A with default values
    A(1, 1) = 1.0; A(1, 2) = 2.0; A(1, 3) = 3.0
    A(2, 1) = 4.0; A(2, 2) = 5.0; A(2, 3) = 6.0
    A(3, 1) = 7.0; A(3, 2) = 8.0; A(3, 3) = 9.0
    
    ! Initialize matrix B with default values
    B(1, 1) = 2.0; B(1, 2) = 0.0; B(1, 3) = 1.0
    B(2, 1) = 1.0; B(2, 2) = 3.0; B(2, 3) = 0.0
    B(3, 1) = 0.0; B(3, 2) = 1.0; B(3, 3) = 2.0
    
    ! Initialize result matrix C to zero
    C = 0.0
    
    ! Print number of available threads
    write(*,'(A,I0)') 'Number of OpenMP threads: ', omp_get_max_threads()
    
    ! Print matrix A
    write(*,*) 'Matrix A:'
    do i = 1, n
        write(*,'(3F8.2)') (A(i, j), j = 1, m)
    end do
    
    ! Print matrix B
    write(*,*) 'Matrix B:'
    do i = 1, m
        write(*,'(3F8.2)') (B(i, j), j = 1, p)
    end do
    
    ! Matrix multiplication using nested do loops with OpenMP parallelization
    ! C(i,j) = sum of A(i,k) * B(k,j) for k = 1 to m
    ! Each element C(i,j) is calculated by a separate thread
    write(*,*) 'Performing parallel matrix multiplication...'
    
    !$omp parallel do collapse(2) private(k) shared(A,B,C)
    do i = 1, n
        do j = 1, p
            do k = 1, m
                C(i, j) = C(i, j) + A(i, k) * B(k, j)
            end do
        end do
    end do
    !$omp end parallel do
    
    ! Print result matrix C
    write(*,*) 'Result matrix C = A * B:'
    do i = 1, n
        write(*,'(3F8.2)') (C(i, j), j = 1, p)
    end do
    
end program matrix_product