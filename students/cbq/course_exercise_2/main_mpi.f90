! Barnes-Hut algorithm with MPI parallelization

program nbody_simulation_mpi
    use iso_fortran_env, only: real64
    use mpi
    use geometry
    use particle
    use barnes_hut
    implicit none
    
    integer :: i, n, output_unit, ierr, rank, nprocs
    integer :: local_start, local_end, local_n
    real(real64) :: dt, t_end, t, dt_out, t_out
    type(particle3d), dimension(:), allocatable :: particles
    type(vector3d), dimension(:), allocatable :: accelerations
    type(cell), pointer :: head, temp_cell
    real(real64) :: mass_temp, start_time, end_time
    type(point3d) :: pos_temp
    type(vector3d) :: vel_temp
    
    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    
    if (rank == 0) start_time = MPI_Wtime()
    
    ! Process reads input
    if (rank == 0) then
        read*, dt
        read*, dt_out
        read*, t_end
        read*, n
    end if
    
    call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    ! Allocate arrays
    allocate(particles(n))
    allocate(accelerations(n))
    
    ! Divide work among processes
    local_n = n / nprocs
    local_start = rank * local_n + 1
    if (rank == nprocs - 1) then
        local_end = n  
    else
        local_end = local_start + local_n - 1
    end if
    
    ! Root reads particle data
    if (rank == 0) then
        do i = 1, n
            read*, mass_temp, pos_temp%x, pos_temp%y, pos_temp%z, &
                  vel_temp%x, vel_temp%y, vel_temp%z
            particles(i)%m = mass_temp
            particles(i)%p = pos_temp
            particles(i)%v = vel_temp
        end do
        
        print*, "Barnes-Hut N-body (MPI)"
        print*, "Particles:", n
        print*, "MPI processes:", nprocs

        open(newunit=output_unit, file='output_mpi.dat', status='replace', action='write')
    end if
    
    ! Broadcast all particle data
    do i = 1, n
        call MPI_Bcast(particles(i)%m, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%p%x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%p%y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%p%z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%v%x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%v%y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(particles(i)%v%z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
    
    ! Build initial tree (all processes)
    allocate(head)
    call calculate_ranges(head, particles, n)
    head%type = 0
    call nullify_pointers(head)
    
    do i = 1, n
        call find_cell(head, temp_cell, particles(i)%p)
        call place_cell(temp_cell, particles(i)%p, i)
    end do
    call borrar_empty_leaves(head)
    call calculate_masses(head, particles)
    
    ! Calculate initial forces (distributed)
    accelerations = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
    do i = local_start, local_end
        call calculate_forces_aux(i, head, particles, accelerations)
    end do
    
    ! Sum accelerations from all processes
    call MPI_Allreduce(MPI_IN_PLACE, accelerations, n*3, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Main integration loop
    t_out = 0.0_real64
    t = 0.0_real64
    
    do while (t <= t_end)
        ! Leapfrog step 1
        do i = 1, n
            particles(i)%v = particles(i)%v + (dt / 2.0_real64) * accelerations(i)
        end do
        
        ! Update positions
        do i = 1, n
            particles(i)%p = particles(i)%p + dt * particles(i)%v
        end do
        
        ! Rebuild tree
        call borrar_tree(head)
        call calculate_ranges(head, particles, n)
        head%type = 0
        call nullify_pointers(head)
        
        do i = 1, n
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do
        call borrar_empty_leaves(head)
        call calculate_masses(head, particles)
        
        ! Calculate forces (distributed)
        accelerations = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
        do i = local_start, local_end
            call calculate_forces_aux(i, head, particles, accelerations)
        end do
        
        ! Sum accelerations
        call MPI_Allreduce(MPI_IN_PLACE, accelerations, n*3, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        
        ! Leapfrog step 2
        do i = 1, n
            particles(i)%v = particles(i)%v + (dt / 2.0_real64) * accelerations(i)
        end do
        
        t_out = t_out + dt
        
        ! Output (only root process)
        if (t_out >= dt_out .and. rank == 0) then
            write(output_unit, '(ES15.7)', advance='no') t
            do i = 1, n
                write(output_unit, '(3ES15.7)', advance='no') &
                    particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
            end do
            write(output_unit, *)
                   
            t_out = 0.0_real64
        end if
        
        t = t + dt
    end do
    
    ! Cleanup
    call borrar_tree(head)
    deallocate(head)
    deallocate(particles)
    deallocate(accelerations)
    
    if (rank == 0) then
        close(output_unit)
        end_time = MPI_Wtime()
        print*, ""
        print*, "Simulation completed!"
        print*, "Output: output_mpi.dat"
        print*, "Wall time:", end_time - start_time, "seconds"
    end if
    
    call MPI_Finalize(ierr)
    
end program nbody_simulation_mpi