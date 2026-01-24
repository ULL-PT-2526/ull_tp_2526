! module for matrix operations
module lecture_08_ex_matrix_ops
    implicit none
    
    ! using a portable kind for double precision
    integer, parameter :: dp = selected_real_kind(p=15, r=307)

    contains

    ! sequential matrix-matrix multiplication
    ! calculates c = a * b
    subroutine matmul_seq(n, a, b, c)
        integer, intent(in) :: n
        real(kind=dp), dimension(n,n), intent(in) :: a, b
        real(kind=dp), dimension(n,n), intent(out) :: c
        
        integer :: i, j, k
        real(kind=dp) :: sum_val

        ! iterate over rows of a
        do i = 1, n
            ! iterate over columns of b
            do j = 1, n
                sum_val = 0.0_dp
                ! perform dot product
                do k = 1, n
                    sum_val = sum_val + a(i,k) * b(k,j)
                end do
                c(i,j) = sum_val
            end do
        end do

    end subroutine matmul_seq


    ! parallel matrix-matrix multiplication using openmp
    ! the outer loops are parallelised so each element (or block of elements)
    ! is calculated by a separate thread.
    subroutine matmul_par(n, a, b, c)
        integer, intent(in) :: n
        real(kind=dp), dimension(n,n), intent(in) :: a, b
        real(kind=dp), dimension(n,n), intent(out) :: c
        
        integer :: i, j, k
        real(kind=dp) :: sum_val

        ! start a parallel region. 
        ! variables a, b, c, n are shared by default.
        ! loop indices i, j, k and sum_val must be private to each thread.
        !$omp parallel do private(i, j, k, sum_val) schedule(static) collapse(2)
        do i = 1, n
            do j = 1, n
                sum_val = 0.0_dp
                do k = 1, n
                    sum_val = sum_val + a(i,k) * b(k,j)
                end do
                c(i,j) = sum_val
            end do
        end do
        !$omp end parallel do

    end subroutine matmul_par

end module lecture_08_ex_matrix_ops