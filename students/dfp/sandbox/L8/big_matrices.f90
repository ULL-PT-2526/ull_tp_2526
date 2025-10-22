program big_matrices
    use omp_lib
    implicit none
    
    ! Matrix dimensions
    integer, parameter :: n = 3000, m = 3000, p = 3000
    
    ! Matrix declarations (allocatable for large matrices)
    real, allocatable :: A(:,:), B(:,:), C(:,:)
    
    ! Loop indices
    integer :: i, j, k
    
    ! Timing variables
    real(8) :: start_time, end_time, execution_time
    
    write(*,*) 'Matrix Multiplication with OpenMP'
    write(*,'(A,I0,A,I0)') 'Matrix size: ', n, ' x ', p
    write(*,'(A,I0)') 'Number of OpenMP threads: ', omp_get_max_threads()
    
    ! Allocate matrices
    write(*,*) 'Allocating matrices...'
    allocate(A(n, m))
    allocate(B(m, p))
    allocate(C(n, p))
    
    ! Initialize matrices A and B with random values
    write(*,*) 'Initializing matrices...'
    
    call random_number(A)
    call random_number(B)
    C = 0.0
    
    ! Matrix multiplication with OpenMP
    write(*,*) 'Starting matrix multiplication...'
    
    start_time = omp_get_wtime()
    
    !$omp parallel do collapse(2) private(k) shared(A,B,C)
    do i = 1, n
        do j = 1, p
            do k = 1, m
                C(i, j) = C(i, j) + A(i, k) * B(k, j)
            end do
        end do
    end do
    !$omp end parallel do
    
    end_time = omp_get_wtime()
    execution_time = end_time - start_time
    
    write(*,'(A,F8.4,A)') 'Execution time: ', execution_time, ' seconds'
    
    ! Deallocate matrices
    deallocate(A, B, C)
    
end program big_matrices
