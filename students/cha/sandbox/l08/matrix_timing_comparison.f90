program matrix_timing_comparison

    use omp_lib
    use iso_fortran_env, only: real64 

    implicit none

    integer, parameter :: dp = real64

    ! Program parameters and variables
    integer, parameter :: matrix_size = 1000  
    real(kind=dp), allocatable :: A(:,:), B(:,:), C(:,:) 
    integer :: i, j, k
    
    ! Timing variables, omp_get_wtime returns a 64-bit real
    real(kind=dp) :: start_time, end_time
    real(kind=dp) :: t_seq, t_par, speedup

    ! Allocate memory for matrices (on the heap)
    allocate(A(matrix_size, matrix_size))
    allocate(B(matrix_size, matrix_size))
    allocate(C(matrix_size, matrix_size))

    ! Fill matrices with random values
    call random_seed()
    call random_number(A)
    call random_number(B)

    ! --- SEQUENTIAL COMPUTATION ---
    print *, "--- Starting SEQUENTIAL computation ---"
    C = 0.0_dp ! Initialize using the correct type
    start_time = omp_get_wtime()

    do i = 1, matrix_size
        do j = 1, matrix_size
            do k = 1, matrix_size
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do

    end_time = omp_get_wtime()
    t_seq = end_time - start_time
    print '("Sequential time: ", F10.6, " seconds")', t_seq
    print *

    ! --- PARALLEL COMPUTATION ---
    print *, "--- Starting PARALLEL computation ---"
    C = 0.0_dp ! Reset the result matrix
    start_time = omp_get_wtime()

    !$omp parallel do collapse(2) private(k)
    do i = 1, matrix_size
        do j = 1, matrix_size
            do k = 1, matrix_size
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do
    !$omp end parallel do

    end_time = omp_get_wtime()
    t_par = end_time - start_time
    print '("Parallel time:   ", F10.6, " seconds")', t_par

    ! --- SPEEDUP CALCULATION ---
    speedup = t_seq / t_par
    print *
    print '("By parallelizing this program, it runs ", F6.2, " times faster.")', speedup

    ! Free allocated memory
    deallocate(A, B, C)

end program matrix_timing_comparison

