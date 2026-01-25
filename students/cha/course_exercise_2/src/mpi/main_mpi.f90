program main_mpi

    use mpi
    use iso_fortran_env, only: real64
    use geometry
    use particle
    use barnes_hut_mpi_module
    use io_module

    ! conditional use only for the OMP library
#ifdef USE_OMP
    use omp_lib
#endif
    implicit none

    ! --- MPI VARIABLES ---
    integer :: ierr, my_rank, n_procs
    integer, allocatable :: sendcounts(:), displs(:)
    integer :: my_n, my_start, my_end, remainder, i_p

    integer :: i, n
    real(real64) :: dt, t_end, t, dt_out, t_out, disk_radius
    real(real64) :: mass_temp
    type(point3d) :: pos_temp
    type(vector3d) :: vel_temp

    real(real64), dimension(:), allocatable :: a_x, a_y, a_z
    type(vector3d) :: acc_vector

    ! TEMPORARY BUFFER ONLY FOR INITIAL BCAST 
    real(real64), allocatable :: init_buf(:)

    ! --- DIAGNOSTICS AND TIMING VARIABLES ---
    integer :: n_threads_omp
    integer(8) :: count_start, count_end, count_rate
    real(real64) :: time_elapsed


    ! ==================
    ! MPI INITIALIZATION
    ! ==================
    call MPI_INIT(ierr)                                     ! initialize MPI environment
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)       ! get current process rank
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr)       ! get total number of processes

    ! --- HEADER OUTPUT (rank 0 only) ---
    n_threads_omp = 1
#ifdef USE_OMP
    n_threads_omp = omp_get_max_threads()    ! detect OMP threads if enabled
#endif

    if (my_rank == 0) then
        print*, ""
        print*, "====================================="
        if (n_threads_omp > 1) then
            print*, "barnes-hut simulation (HYBRID MPI+OpenMP)"
        else
            print*, "barnes-hut simulation (MPI Pure)"
        end if
        print*, "====================================="
        print*, ""
        print*, "MPI Processes:    ", n_procs
        if (n_threads_omp > 1) then
            print*, "OpenMP Threads/Proc:", n_threads_omp
            print*, "Total Cores Used: ", n_procs * n_threads_omp
        end if
        print*, ""
    end if


    ! ======================================
    ! READING INPUT PARAMETERS (Rank 0 only)
    ! ======================================
    if (my_rank == 0) then
        print*, "reading input parameters..."
        read*, dt
        read*, dt_out
        read*, t_end
        read*, disk_radius
        read*, n

        print*, "  dt          = ", dt
        print*, "  dt_out      = ", dt_out
        print*, "  t_end       = ", t_end
        print*, "  disk_R      = ", disk_radius
        print*, "  n_particles = ", n
        print*, ""
    end if

    ! broadcast simulation parameters to all ranks
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(disk_radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


    ! =====================
    ! WORKLOAD DISTRIBUTION
    ! =====================
    allocate(sendcounts(n_procs), displs(n_procs))
    remainder = mod(n, n_procs)
    do i_p = 0, n_procs - 1
        sendcounts(i_p + 1) = n / n_procs
        if (i_p < remainder) sendcounts(i_p + 1) = sendcounts(i_p + 1) + 1
    end do
    displs(1) = 0
    do i_p = 2, n_procs
        displs(i_p) = displs(i_p-1) + sendcounts(i_p-1)
    end do
    
    my_n = sendcounts(my_rank + 1)                          ! particles handled by this rank
    my_start = displs(my_rank + 1) + 1                      ! global start index
    my_end = my_start + my_n - 1                            ! global end index


    ! ===============================
    ! ALLOCATION & INITIAL CONDITIONS
    ! ===============================
    n_particles = n                                         ! set module global variable
    allocate(particles(n))                                  ! all ranks store all particles
    allocate(a_x(n), a_y(n), a_z(n))
    allocate(init_buf(n*7))                                 ! temporary buffer for initial IO

    if (my_rank == 0) then
        print*, "reading particle initial conditions..."
        do i = 1, n
            read*, mass_temp, pos_temp%x, pos_temp%y, pos_temp%z, &
                   vel_temp%x, vel_temp%y, vel_temp%z
            particles(i)%m = mass_temp
            particles(i)%p = pos_temp
            particles(i)%v = vel_temp
            
            ! pack into flat buffer for initial BCAST
            init_buf((i-1)*7 + 1) = particles(i)%p%x
            init_buf((i-1)*7 + 2) = particles(i)%p%y
            init_buf((i-1)*7 + 3) = particles(i)%p%z
            init_buf((i-1)*7 + 4) = particles(i)%v%x
            init_buf((i-1)*7 + 5) = particles(i)%v%y
            init_buf((i-1)*7 + 6) = particles(i)%v%z
            init_buf((i-1)*7 + 7) = particles(i)%m
        end do
        print*, "read ", n, " particles"
        print*, ""
    end if

    ! distribute initial state to all processes
    call MPI_BCAST(init_buf, n*7, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    if (my_rank /= 0) then    ! unpack buffer on other ranks
        do i = 1, n
            particles(i)%p%x = init_buf((i-1)*7 + 1)
            particles(i)%p%y = init_buf((i-1)*7 + 2)
            particles(i)%p%z = init_buf((i-1)*7 + 3)
            particles(i)%v%x = init_buf((i-1)*7 + 4)
            particles(i)%v%y = init_buf((i-1)*7 + 5)
            particles(i)%v%z = init_buf((i-1)*7 + 6)
            particles(i)%m   = init_buf((i-1)*7 + 7)
        end do
    end if
    deallocate(init_buf)  ! free temporary buffer


    ! ===================
    ! SETUP OUTPUT & TREE
    ! ===================
    if (my_rank == 0) then 
        print*, "opening output file..."
        call io_open_output(n, dt, dt_out, t_end, disk_radius)
        print*, ""
    end if
    
    if (my_rank == 0) print*, "initializing barnes-hut tree..."
    call barnes_hut_initialize()    ! build initial tree locally
    if (my_rank == 0) print*, ""

    if (my_rank == 0) print*, "calculating initial accelerations..."
    a_x = 0.0_real64; a_y = 0.0_real64; a_z = 0.0_real64
    call barnes_hut_calculate_forces(a_x, a_y, a_z, my_start, my_end) ! compute local forces
    if (my_rank == 0) print*, ""


    ! ====================
    ! MAIN SIMULATION LOOP
    ! ====================
    if (my_rank == 0) then
        print*, "starting main simulation loop..."
        print*, ""
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)      ! sync before timing
    call system_clock(count_start, count_rate)

    t = 0.0_real64
    t_out = 0.0_real64

    do while (t < t_end)
        
        ! Leapfrog Step 1 (velocity kick)
        do i = my_start, my_end
            acc_vector = vector3d(a_x(i), a_y(i), a_z(i))
            particles(i)%v = particles(i)%v + acc_vector * (dt * 0.5_real64)
        end do

        ! update positions (drift)
        do i = my_start, my_end
            particles(i)%p = particles(i)%p + particles(i)%v * dt
        end do

        ! sync positions (version 1: Encapsulated inside module)
        ! This gathers updated positions from all procs to all procs
        call barnes_hut_sync_particles(my_start, my_n, sendcounts, displs)

        ! rebuild tree (replicated on every process)
        call barnes_hut_rebuild()

        ! compute forces (local subset)
        a_x = 0.0_real64; a_y = 0.0_real64; a_z = 0.0_real64
        call barnes_hut_calculate_forces(a_x, a_y, a_z, my_start, my_end)

        ! leapfrog step 2 (velocity kick)
        do i = my_start, my_end
            acc_vector = vector3d(a_x(i), a_y(i), a_z(i))
            particles(i)%v = particles(i)%v + acc_vector * (dt * 0.5_real64)
        end do

        ! output (master rank only)
        t_out = t_out + dt
        if (t_out >= dt_out) then
            if (my_rank == 0) call io_write_simulation_state(t, particles, n)
            t_out = 0.0_real64
        end if
        
        t = t + dt
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! sync after loop
    call system_clock(count_end)
    time_elapsed = real(count_end - count_start, real64) / real(count_rate, real64)

    if (my_rank == 0) then
        print*, ""
        print*, "simulation finished"
        print*, "-------------------------------------"
        write(*, '(A, F10.2, A)') "Total execution time: ", time_elapsed, " seconds"
        print*, "-------------------------------------"
        print*, ""
        call io_close_output()
    end if

    ! --- Cleanup ---
    deallocate(particles, a_x, a_y, a_z)

    if (my_rank == 0) then
        print*, "all resources released"
        print*, "====================================="
        print*, ""
    end if

    call MPI_FINALIZE(ierr)  ! shut down MPI

end program main_mpi