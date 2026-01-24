program lecture_08_ex
    use lecture_08_ex_matrix_ops
    use omp_lib
    implicit none

    ! define matrix size
    integer, parameter :: n = 200

    ! arrays for matrices
    real(kind=dp), dimension(n,n) :: a, b, c_seq, c_par
    
    ! variables for timing
    integer :: t_start, t_end, t_rate
    real(kind=dp) :: time_seq, time_par

    ! Threads variables
    integer :: requested_threads, actual_threads

    ! Set the number of threads directly in the code
    requested_threads = 8
    call omp_set_num_threads(requested_threads)

    ! initialize random seed
    call random_seed()

    ! fill matrices a and b with random numbers [0, 1]
    call random_number(a)
    call random_number(b)

    print *, "matrix size: ", n, "x", n


    ! 1. run sequential version
    call system_clock(t_start, t_rate)
    
    call matmul_seq(n, a, b, c_seq)
    
    call system_clock(t_end)
    
    time_seq = real(t_end - t_start, dp) / real(t_rate, dp)
    print '(a, f10.6, a)', "sequential time: ", time_seq, " s"


    ! 2. run parallel version

    ! open a small parallel region 
    ! to obtain the number of threads
    !$omp parallel
    !$omp master
    actual_threads = omp_get_num_threads()
    print *, "running parallel version with", actual_threads, "threads"
    !$omp end master
    !$omp end parallel

    call system_clock(t_start, t_rate)
    
    call matmul_par(n, a, b, c_par)
    
    call system_clock(t_end)
    
    time_par = real(t_end - t_start, dp) / real(t_rate, dp)
    print '(a, f10.6, a)', "parallel time:   ", time_par, " s"


    ! 3. verify results
    if (maxval(abs(c_seq - c_par)) < 1.0e-10_dp) then
        print *, "check: results match!"
    else
        print *, "check: error! results do not match."
    end if

    ! print speedup
    if (time_par > 0.0_dp) then
        print '(a, f6.2, a)', "speedup: ", time_seq / time_par, "x"
    end if

end program lecture_08_ex