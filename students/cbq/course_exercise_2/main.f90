! Main Program - N-body Simulation

! Barnes-Hut N-body gravitational simulation
! Supports serial and OpenMP parallel execution

program nbody_simulation
    use iso_fortran_env, only: real64
    use geometry
    use particle
    use barnes_hut
!$ use omp_lib
    implicit none
    
    integer :: i, n, output_unit
    real(real64) :: dt, t_end, t, dt_out, t_out
    type(particle3d), dimension(:), allocatable :: particles
    type(vector3d), dimension(:), allocatable :: accelerations
    type(cell), pointer :: head, temp_cell
    real(real64) :: mass_temp
    type(point3d) :: pos_temp
    type(vector3d) :: vel_temp
!$ real(real64) :: start_time, end_time, force_time
!$ integer :: num_threads
    
!$  start_time = omp_get_wtime()
    
    ! Read input parameters
    read*, dt
    read*, dt_out
    read*, t_end
    read*, n
    
    ! Allocate arrays
    allocate(particles(n))
    allocate(accelerations(n))
    
    ! Read particle data
    do i = 1, n
        read*, mass_temp, pos_temp%x, pos_temp%y, pos_temp%z, &
              vel_temp%x, vel_temp%y, vel_temp%z
        particles(i)%m = mass_temp
        particles(i)%p = pos_temp
        particles(i)%v = vel_temp
    end do
    
    ! Print simulation info
    print*, "Barnes-Hut N-body Simulation"
    print*, "Particles:", n
    print*, "Time step:", dt
    print*, "End time:", t_end
!$ num_threads = 1
!$ !$omp parallel
!$ num_threads = omp_get_num_threads()
!$ !$omp end parallel
!$ print*, "OpenMP threads:", num_threads
    print*, ""
    
    ! Open output file
    open(newunit=output_unit, file='output.dat', status='replace', action='write')
    
    ! Initialize tree
    allocate(head)
    call calculate_ranges(head, particles, n)
    head%type = 0
    call nullify_pointers(head)
    
    ! Build initial tree
    do i = 1, n
        call find_cell(head, temp_cell, particles(i)%p)
        call place_cell(temp_cell, particles(i)%p, i)
    end do
    call borrar_empty_leaves(head)
    call calculate_masses(head, particles)
    
    ! Calculate initial forces
    accelerations = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
    call calculate_forces(head, particles, accelerations, n)
    
    ! Main time integration loop
    t_out = 0.0_real64
    t = 0.0_real64
  !$  force_time = 0.0_real64
    
    do while (t <= t_end)
        ! Leapfrog step 1: advance velocities by half step
        do i = 1, n
            particles(i)%v = particles(i)%v + (dt / 2.0_real64) * accelerations(i)
        end do
        
        ! Update positions
        do i = 1, n
            particles(i)%p = particles(i)%p + dt * particles(i)%v
        end do
        
        ! Rebuild tree with new positions
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
        
        ! Calculate new forces (timed section)
!$      start_time = omp_get_wtime()
        accelerations = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
        call calculate_forces(head, particles, accelerations, n)
!$      force_time = force_time + (omp_get_wtime() - start_time)
        
        ! Leapfrog step 2: advance velocities by half step
        do i = 1, n
            particles(i)%v = particles(i)%v + (dt / 2.0_real64) * accelerations(i)
        end do
        
        ! Output data at intervals
        t_out = t_out + dt
        if (t_out >= dt_out) then
            ! Write to file: time x1 y1 z1 x2 y2 z2 ...
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
    close(output_unit)
    
!$  end_time = omp_get_wtime()
    
    ! Print timing results
    print*, "Simulation completed!"
    print*, "Output written to: output.dat"
!$  print*, "Total wall time:", end_time - start_time, "seconds"
!$  print*, "Force calculation time:", force_time, "seconds"
    
end program nbody_simulation