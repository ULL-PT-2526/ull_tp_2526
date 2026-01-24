! main program for exercise 2 - MPI Version
program ex2_mpi
    use geometry
    use particle
    use barnes_hut_module
    implicit none

    include 'mpif.h'    

    ! MPI Variables
    integer :: rank, nprocs, ierr

    ! variables for simulation setup
    integer :: n_particles, i, step
    real(kind=dp) :: dt, t_end, t, dt_out, t_out
    
    ! arrays
    type(particle3d), dimension(:), allocatable :: particles
    type(vector3d), dimension(:), allocatable :: acc_local
    type(vector3d), dimension(:), allocatable :: acc_global
    
    ! work distribution variables
    integer :: n_local, i_start, i_end, remainder
    
    ! octree root
    type(cell_t), pointer :: root => null()

    ! i/o variablesariables
    integer, parameter :: u_in = 10, u_out = 20
    character(len=256) :: filename_in = 'input.dat'
    character(len=256) :: filename_out = 'output.dat'
    character(len=256) :: arg_string
    integer :: num_args, arg_idx
    logical :: input_file_set = .false.
    real(kind=dp) :: m, rx, ry, rz, vx, vy, vz
    
    ! simulation time variables
    double precision :: t_start, t_final
    real(kind=dp) :: total_time

    t_start = 0.0d0
    t_final = 0.0d0
    
    ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    
    ! Rank 0 reads input and broadcasts parameters
    if (rank == 0) then

        filename_in = 'input.dat'
        filename_out = 'output.dat'
        num_args = command_argument_count()
        arg_idx = 1

        do while (arg_idx <= num_args)
            call get_command_argument(arg_idx, arg_string)
            if (trim(arg_string) == '-i') then
                arg_idx = arg_idx + 1
                call get_command_argument(arg_idx, filename_in)
            else if (trim(arg_string) == '-o') then
                arg_idx = arg_idx + 1
                call get_command_argument(arg_idx, filename_out)
            end if
            arg_idx = arg_idx + 1
        end do

        print *, "MPI Master (Rank 0) reading from ", trim(filename_in)
        
        inquire(file=trim(filename_in), exist=input_file_set)
        if (.not. input_file_set) then
            print *, "Error: input file not found: ", trim(filename_in)
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if

        open(unit=u_in, file=trim(filename_in), status='old', action='read')
        read(u_in, *) dt
        read(u_in, *) dt_out
        read(u_in, *) t_end
        read(u_in, *) n_particles
        
        print *, "Simulation setup: N=", n_particles, " Procs=", nprocs
        
        allocate(particles(n_particles))
        
        do i = 1, n_particles
            read(u_in, *) m, rx, ry, rz, vx, vy, vz
            particles(i)%m = m
            particles(i)%p = point3d(rx, ry, rz)
            particles(i)%v = vector3d(vx, vy, vz)
        end do
        close(u_in)
        
        open(unit=u_out, file=trim(filename_out), status='replace', action='write')
    end if

    ! broadcast simulation parameters to all processes
    call MPI_BCAST(n_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory on worker processes
    if (rank /= 0) allocate(particles(n_particles))
    allocate(acc_local(n_particles))  ! Only stores forces calculated by this rank
    allocate(acc_global(n_particles)) ! Stores summed forces from all ranks

    
    ! broadcast initial particle data
    call MPI_BCAST(particles%m, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    ! temporary arrays to send point3d/vector3d components 
    call MPI_BCAST(particles%p%x, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(particles%p%y, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(particles%p%z, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(particles%v%x, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(particles%v%y, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(particles%v%z, n_particles, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    
    ! define workload distribution 

    ! divide particles evenly among ranks
    n_local = n_particles / nprocs
    remainder = mod(n_particles, nprocs)
    
    
    if (rank < remainder) then
        n_local = n_local + 1
        i_start = rank * n_local + 1
    else
        i_start = rank * n_local + remainder + 1
    end if
    i_end = i_start + n_local - 1

    ! start timer
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank == 0) t_start = MPI_WTIME()

    
    ! simulation Loop
    t = 0.0_dp
    t_out = 0.0_dp
    step = 0

    ! initial output (rank 0 only)
    if (rank == 0) then
        call write_output(u_out, t, particles, n_particles)
        t_out = t_out + dt_out
    end if

    do while (t < t_end)
        
        !--- build tree for every process ---!
        call build_tree(root, particles)

        !--- calculate the distributed forces ---!
        acc_local = vector3d(0.0_dp, 0.0_dp, 0.0_dp) 
        
        !--- each rank only calculates forces for its subset ---!
        do i = i_start, i_end
            call calculate_force_recursive(root, particles(i), acc_local(i))
        end do
        
        !--- reduce forces ---!
        call MPI_ALLREDUCE(acc_local%x, acc_global%x, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc_local%y, acc_global%y, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc_local%z, acc_global%z, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

        !--- velocity verlet ---!
        do i = 1, n_particles
            particles(i)%v = particles(i)%v + acc_global(i) * (0.5_dp * dt)
            particles(i)%p = particles(i)%p + particles(i)%v * dt
        end do
        
        !--- re-built tree ---!
        call delete_tree(root)
        call build_tree(root, particles)
        
        !--- calculate forces again ---!
        acc_local = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        do i = i_start, i_end
            call calculate_force_recursive(root, particles(i), acc_local(i))
        end do
        
        !--- reduce forces again ---!
        call MPI_ALLREDUCE(acc_local%x, acc_global%x, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc_local%y, acc_global%y, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(acc_local%z, acc_global%z, n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

        !--- verlet velocity again ---!
        do i = 1, n_particles
            particles(i)%v = particles(i)%v + acc_global(i) * (0.5_dp * dt)
        end do

        !--- output ---!
        t = t + dt
        step = step + 1

        if (rank == 0) then
            if (t >= t_out - 1.e-8) then
                call write_output(u_out, t, particles, n_particles)
                t_out = t_out + dt_out
            end if
            if (mod(step, 100) == 0) print *, "Step:", step, " Time:", t
        end if

        call delete_tree(root)
    end do

    
    ! simulation total time
    if (rank == 0) then
        t_final = MPI_WTIME()
        total_time = t_final - t_start
        print '(a, f15.6, a)', " Computation Time: ", total_time, " seconds"
        close(u_out)
    end if

    deallocate(particles)
    deallocate(acc_local)
    deallocate(acc_global)

    call MPI_FINALIZE(ierr)

contains

    ! write output subrutine
    subroutine write_output(unit_num, current_time, p_array, num_p)
        integer, intent(in) :: unit_num, num_p
        real(kind=dp), intent(in) :: current_time
        type(particle3d), dimension(:), intent(in) :: p_array
        integer :: k
        
        write(unit_num, '(f10.5)', advance='no') current_time
        do k = 1, num_p
            write(unit_num, '(3(1x, es12.5))', advance='no') p_array(k)%p%x, p_array(k)%p%y, p_array(k)%p%z
        end do
        write(unit_num, *)
    end subroutine write_output

end program ex2_mpi