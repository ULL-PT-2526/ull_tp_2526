program main

    use iso_fortran_env, only: real64
    use geometry
    use particle
    use barnes_hut_module

    ! conditional use only for the OMP library
#ifdef USE_OMP
    use omp_lib
#endif

    use io_module
    implicit none

    integer :: i, n                                  ! loop index and number of particles
    real(real64) :: dt, t_end, t, dt_out, t_out      ! time variables
    real(real64) :: disk_radius                      ! disk radius
    real(real64) :: mass_temp                        ! temporary mass
    type(point3d) :: pos_temp                        ! temporary position
    type(vector3d) :: vel_temp                       ! temporary velocity

    ! acceleration arrays (required by barnes_hut_calculate_forces)
    real(real64), dimension(:), allocatable :: a_x, a_y, a_z

    ! temporary vector used for vectorized operations
    type(vector3d) :: acc_vector

    ! wall-clock timing variables
    integer(8) :: count_start, count_end, count_rate
    real(real64) :: time_elapsed


    print*, ""
    print*, "====================================="
    print*, "barnes-hut n-body simulation"
    print*, "====================================="
    print*, ""

#ifdef USE_OMP
    print*, "Parallel execution: YES"
    print*, "Number of threads: ", omp_get_max_threads()
#else
    print*, "Parallel execution: NO (Serial)"
#endif
    print*, ""

    ! we read input parameters
    print*, "reading input parameters..."
    read*, dt                                          ! integration timestep
    read*, dt_out                                      ! output timestep
    read*, t_end                                       ! end time
    read*, disk_radius                                 ! disk radius
    read*, n                                           ! number of particles

    print*, "  dt          = ", dt
    print*, "  dt_out      = ", dt_out
    print*, "  t_end       = ", t_end
    print*, "  disk_R      = ", disk_radius
    print*, "  n_particles = ", n
    print*, ""

    n_particles = n                                   ! set global particle count
    allocate(particles(n))                            ! allocate particle array
    allocate(a_x(n), a_y(n), a_z(n))                  ! allocate acceleration arrays

    ! we read initial conditions
    print*, "reading particle initial conditions..."
    do i = 1, n
        read*, mass_temp,                             & ! particle mass
               pos_temp%x, pos_temp%y, pos_temp%z,    & ! initial position
               vel_temp%x, vel_temp%y, vel_temp%z       ! initial velocity
        particles(i)%m = mass_temp                    ! store mass
        particles(i)%p = pos_temp                     ! store position
        particles(i)%v = vel_temp                     ! store velocity
    end do
    print*, "read ", n, " particles"
    print*, ""

    ! we open output file
    print*, "opening output file..."
    call io_open_output(n, dt, dt_out, t_end, disk_radius)
    print*, ""

    ! we initialize Barnes–Hut tree
    print*, "initializing barnes-hut tree..."
    call barnes_hut_initialize()                      ! build initial octree
    print*, ""

    ! we compute initial accelerations
    print*, "calculating initial accelerations..."
    a_x = 0.0_real64                                  ! reset accelerations
    a_y = 0.0_real64
    a_z = 0.0_real64
    call barnes_hut_calculate_forces(a_x, a_y, a_z)   ! compute gravitational forces
    print*, ""

    ! MAIN SIMULATION LOOP
    print*, "starting main simulation loop..."
    print*, ""

    ! --- START TIMER (wall clock) ---
    call system_clock(count_start, count_rate)
    ! --------------------------------

    t_out = 0.0_real64       ! output time accumulator
    t = 0.0_real64           ! simulation time

    do while (t < t_end)

        ! leapfrog: first half-step velocity update
        ! v(t+dt/2) = v(t) + a(t)*dt/2
        do i = 1, n
            acc_vector = vector3d(a_x(i), a_y(i), a_z(i))
            particles(i)%v = particles(i)%v + acc_vector * (dt * 0.5_real64)
        end do

        ! update positions
        ! r(t+dt) = r(t) + v(t+dt/2)*dt
        do i = 1, n
            particles(i)%p = particles(i)%p + particles(i)%v * dt
        end do

        ! we rebuild Barnes–Hut tree after particle motion
        call barnes_hut_rebuild()

        ! we compute new accelerations
        a_x = 0.0_real64
        a_y = 0.0_real64
        a_z = 0.0_real64
        call barnes_hut_calculate_forces(a_x, a_y, a_z)

        ! leapfrog: second half-step velocity update
        ! v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2
        do i = 1, n
            acc_vector = vector3d(a_x(i), a_y(i), a_z(i))
            particles(i)%v = particles(i)%v + acc_vector * (dt * 0.5_real64)
        end do

        ! output state 
        t_out = t_out + dt
        if (t_out >= dt_out) then
            call io_write_simulation_state(t, particles, n)
            t_out = 0.0_real64
        end if

        t = t + dt                                   ! we advance simulation time

    end do

    ! --- STOP TIMER ---
    call system_clock(count_end)
    time_elapsed = real(count_end - count_start, real64) / &
                   real(count_rate, real64)


    print*, ""
    print*, "simulation finished"
    print*, "-------------------------------------"
    write(*, '(A, F10.2, A)') "Total execution time: ", time_elapsed, " seconds"
    print*, "-------------------------------------"
    print*, ""

    ! we close output file
    call io_close_output()

    ! and finaly we cleanup memory
    deallocate(particles, a_x, a_y, a_z)
    if (associated(head)) deallocate(head)

    print*, "all resources released"
    print*, "====================================="
    print*, ""


end program main