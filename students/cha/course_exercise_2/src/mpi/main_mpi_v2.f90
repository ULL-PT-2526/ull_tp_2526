program main_mpi_v2

    use mpi
    use iso_fortran_env, only: real64
    use geometry
    use particle
    use barnes_hut_mpi_module_v2
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
    
    ! --- COMMUNICATION BUFFERS ---
    real(real64), allocatable :: sendbuf(:), recvbuf(:)

    integer :: i, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    real(real64) :: disk_radius
    real(real64) :: mass_temp
    type(point3d) :: pos_temp
    type(vector3d) :: vel_temp

    real(real64), dimension(:), allocatable :: a_x, a_y, a_z
    type(vector3d) :: acc_vector

    ! --- TIMING VARIABLES ---
    integer :: n_threads_omp
    integer(8) :: count_start, count_end, count_rate
    real(real64) :: time_elapsed
    
    

    ! ==================
    ! MPI INITIALIZATION
    ! ==================
    call MPI_INIT(ierr)                                     ! initialize MPI environment
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)       ! get current process rank
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr)       ! get total number of processes

    ! detect hybrid mode vs pure mode
    n_threads_omp = 1                                       ! default to 1 (serial within process)
#ifdef USE_OMP
    n_threads_omp = omp_get_max_threads()                   ! detect OMP threads if enabled
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

    ! broadcast global parameters
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(disk_radius, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


    ! =====================
    ! WORKLOAD DISTRIBUTION
    ! =====================
    allocate(sendcounts(n_procs))
    allocate(displs(n_procs))

    remainder = mod(n, n_procs)
    do i_p = 0, n_procs - 1
        sendcounts(i_p + 1) = n / n_procs
        if (i_p < remainder) then
            sendcounts(i_p + 1) = sendcounts(i_p + 1) + 1
        end if
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
    n_particles = n 
    allocate(particles(n))
    allocate(a_x(n), a_y(n), a_z(n))
    
    allocate(sendbuf(my_n * 7))                             ! buffer for local particles
    allocate(recvbuf(n * 7))                                ! buffer for all particles

    if (my_rank == 0) then
        print*, "reading particle initial conditions..."
        do i = 1, n
            read*, mass_temp, &
                   pos_temp%x, pos_temp%y, pos_temp%z, &
                   vel_temp%x, vel_temp%y, vel_temp%z
            particles(i)%m = mass_temp
            particles(i)%p = pos_temp
            particles(i)%v = vel_temp
        end do
        print*, "read ", n, " particles"
        print*, ""
        
        call pack_subset_particles(particles, 1, n, recvbuf) ! pack initial state
    end if

    call MPI_BCAST(recvbuf, n*7, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (my_rank /= 0) call unpack_all_particles(recvbuf, n, particles)


    ! ===================
    ! SETUP OUTPUT & TREE
    ! ===================
    if (my_rank == 0) then
        print*, "opening output file..."
        call io_open_output(n, dt, dt_out, t_end, disk_radius)
        print*, ""
    end if

    if (my_rank == 0) print*, "initializing barnes-hut tree..."
    call barnes_hut_initialize()  ! build initial tree locally
    if (my_rank == 0) print*, ""

    if (my_rank == 0) print*, "calculating initial accelerations..."
    a_x = 0.0_real64
    a_y = 0.0_real64
    a_z = 0.0_real64
    
    call barnes_hut_calculate_forces(a_x, a_y, a_z, my_start, my_end) ! compute local forces
    
    if (my_rank == 0) print*, ""


    ! ====================
    ! MAIN SIMULATION LOOP
    ! ====================
    if (my_rank == 0) then
        print*, "starting main simulation loop..."
        print*, ""
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! sync before timing
    call system_clock(count_start, count_rate)
    
    t_out = 0.0_real64
    t = 0.0_real64

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

        ! sync positions (version 2: explicit communication)
        call pack_subset_particles(particles, my_start, my_n, sendbuf)
        call MPI_ALLGATHERV(sendbuf, my_n*7, MPI_DOUBLE_PRECISION, &
                            recvbuf, sendcounts*7, displs*7, MPI_DOUBLE_PRECISION, &
                            MPI_COMM_WORLD, ierr)
        call unpack_all_particles(recvbuf, n, particles)

        ! rebuild tree (replicated on every process)
        call barnes_hut_rebuild()

        ! compute forces (local subset)
        a_x = 0.0_real64
        a_y = 0.0_real64
        a_z = 0.0_real64
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

    deallocate(particles, a_x, a_y, a_z, sendbuf, recvbuf)
    if (associated(head)) deallocate(head)

    if (my_rank == 0) then
        print*, "all resources released"
        print*, "====================================="
        print*, ""
    end if

    call MPI_FINALIZE(ierr)    ! shut down MPI

contains

    subroutine pack_subset_particles(p_array, start_idx, count, buf)
        type(particle3d), intent(in) :: p_array(:)
        integer, intent(in) :: start_idx, count
        real(real64), intent(out) :: buf(:)
        integer :: j, p_idx
        do j = 1, count
            p_idx = start_idx + j - 1
            buf((j-1)*7 + 1) = p_array(p_idx)%p%x
            buf((j-1)*7 + 2) = p_array(p_idx)%p%y
            buf((j-1)*7 + 3) = p_array(p_idx)%p%z
            buf((j-1)*7 + 4) = p_array(p_idx)%v%x
            buf((j-1)*7 + 5) = p_array(p_idx)%v%y
            buf((j-1)*7 + 6) = p_array(p_idx)%v%z
            buf((j-1)*7 + 7) = p_array(p_idx)%m
        end do
    end subroutine pack_subset_particles

    subroutine unpack_all_particles(buf, count, p_array)
        real(real64), intent(in) :: buf(:)
        integer, intent(in) :: count
        type(particle3d), intent(out) :: p_array(:)
        integer :: j
        do j = 1, count
            p_array(j)%p%x = buf((j-1)*7 + 1)
            p_array(j)%p%y = buf((j-1)*7 + 2)
            p_array(j)%p%z = buf((j-1)*7 + 3)
            p_array(j)%v%x = buf((j-1)*7 + 4)
            p_array(j)%v%y = buf((j-1)*7 + 5)
            p_array(j)%v%z = buf((j-1)*7 + 6)
            p_array(j)%m   = buf((j-1)*7 + 7)
        end do
    end subroutine unpack_all_particles

end program main_mpi_v2