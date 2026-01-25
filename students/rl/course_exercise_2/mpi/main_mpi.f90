program main_mpi
    use mpi
    use geometry
    use particle
    use barnes_hut
    implicit none

    ! MPI variables
    integer :: my_rank, p, ierr
    integer :: my_n, my_start, my_end

    ! Simulation parameters
    real(dp) :: dt, dt_out, t_end, t, t_out
    real(dp) :: t_i, t_f 
    integer :: i
    integer :: output_unit

    !---------------------------------------------------------
    ! MPI initialization
    !---------------------------------------------------------
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_i = MPI_Wtime()

    !---------------------------------------------------------
    ! Open output file
    !---------------------------------------------------------

    if (my_rank == 0) then
        output_unit = 10
        open(unit=output_unit, file="output.dat", status="replace", action="write")
    end if

    !---------------------------------------------------------
    ! Read input on master
    !---------------------------------------------------------
    if (my_rank == 0) then
        read(*,*) dt
        read(*,*) dt_out
        read(*,*) t_end
        read(*,*) n
    end if

    !---------------------------------------------------------
    ! Broadcast parameters to the others
    !---------------------------------------------------------
    call MPI_Bcast(n,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dt,     1, MPI_REAL8,   0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dt_out, 1, MPI_REAL8,   0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(t_end,  1, MPI_REAL8,   0, MPI_COMM_WORLD, ierr)

    !---------------------------------------------------------
    ! Allocate arrays on all processes
    !---------------------------------------------------------
    allocate(parts(n))
    allocate(a(n))

    !---------------------------------------------------------
    ! Read particles on master
    !---------------------------------------------------------
    if (my_rank == 0) then
        do i = 1, n
            read(*,*) parts(i)%m, &
                       parts(i)%p%x, parts(i)%p%y, parts(i)%p%z, &
                       parts(i)%v%x, parts(i)%v%y, parts(i)%v%z
        end do
    end if

    !---------------------------------------------------------
    ! Broadcast particle data to the others
    !---------------------------------------------------------
    call MPI_Bcast(parts(:)%m,   n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%p%x, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%p%y, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%p%z, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%v%x, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%v%y, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(parts(:)%v%z, n, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    !---------------------------------------------------------
    ! Work partitioning; from here on there is no master/slave distinction,
    ! all processors collaborate equally
    ! (except when printing and performing the time check)
    !---------------------------------------------------------
    my_n     = n/p
    my_start = my_rank*my_n+1
    my_end   = my_start+my_n-1

    !---------------------------------------------------------
    ! Build initial Barnes–Hut tree
    !---------------------------------------------------------
    allocate(head)
    call calculate_ranges(head)
    head%type = 0
    call nullify_pointers(head)

    do i = 1, n
        call find_cell(head, temp_cell, parts(i)%p)
        call place_cell(temp_cell, parts(i)%p, i)
    end do

    call borrar_empty_leaves(head)
    call calculate_masses(head)

    !---------------------------------------------------------
    ! Initial accelerations (local block only)
    !---------------------------------------------------------
    a(my_start:my_end) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
    call calculate_forces_mpi(head, my_start, my_end)

    !---------------------------------------------------------
    ! Time integration loop
    !---------------------------------------------------------
    t_out = 0.0_dp
    t = 0.0_dp

    do while (t < t_end)

        ! Velocity half step
        do i = my_start, my_end
            parts(i)%v = parts(i)%v + a(i)*(dt/2.0_dp)
        end do

        ! Position update
        do i = my_start, my_end
            parts(i)%p = parts(i)%p + parts(i)%v*dt
        end do

        ! Synchronize positions
        call MPI_Allgather(parts(my_start)%p%x, my_n, MPI_REAL8, &
                           parts(1)%p%x,       my_n, MPI_REAL8, MPI_COMM_WORLD, ierr)
        call MPI_Allgather(parts(my_start)%p%y, my_n, MPI_REAL8, &
                           parts(1)%p%y,       my_n, MPI_REAL8, MPI_COMM_WORLD, ierr)
        call MPI_Allgather(parts(my_start)%p%z, my_n, MPI_REAL8, &
                           parts(1)%p%z,       my_n, MPI_REAL8, MPI_COMM_WORLD, ierr)

        ! Rebuild tree
        call borrar_tree(head)
        call calculate_ranges(head)
        head%type = 0
        call nullify_pointers(head)

        do i = 1, n
            call find_cell(head, temp_cell, parts(i)%p)
            call place_cell(temp_cell, parts(i)%p, i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head)

        ! New accelerations (local)
        a(my_start:my_end) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        call calculate_forces_mpi(head, my_start, my_end)

        ! Velocity half step
        do i = my_start, my_end
            parts(i)%v = parts(i)%v + a(i) * (dt / 2.0_dp)
        end do

        ! Output (master only)
        t_out = t_out + dt
        if (my_rank == 0) then
            if (t_out >= dt_out) then

                write(output_unit, '(ES20.10)', advance='no') t
                do i = 1, n
                    write(output_unit, '(3ES20.10)', advance='no') &
                        parts(i)%p%x, parts(i)%p%y, parts(i)%p%z
                end do
                write(output_unit, *)   ! newline

                t_out = 0.0_dp
            end if
        end if

        t = t + dt
    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_f = MPI_Wtime()

    ! Print execution time 
    if (my_rank == 0) then
        close(output_unit)
        print*, "Simulation complete. Output written to output.dat"
        print *, "Total execution time:", t_f - t_i
    end if
    
    call MPI_Finalize(ierr)

end program main_mpi
