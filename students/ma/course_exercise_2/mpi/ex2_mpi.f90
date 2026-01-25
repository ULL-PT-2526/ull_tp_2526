
program ex2_mpi

    use geometry
    use particles
    use tree_mpi

    implicit none

    integer :: ierr ! this is to check and deallocate
    character(len=100) :: infile   ! name of the input file
    character(len=200) :: in ! argument in the terminal for the input file
    real(dp) :: t1, t2

    !! mpi initialization
    !!!!!!!!!!!!!!!!!!!!!!!!
    call mpi_init(error)
    call mpi_comm_size(mpi_comm_world, p, error)
    call mpi_comm_rank(mpi_comm_world, my_rank, error)

    call get_command_argument(1,in)

    if (len_trim(in)==0) then ! if input is not given 
        print *, "Input example: ./ex2 <inputfile>"
        stop
    endif

    infile = trim(in) ! remove empty space

    if (my_rank == 0) then 
        print *, "Reading input file: ", trim(infile)
    end if

    ! we use the subroutine to read the file, and initialize mpi
    call read_input(infile,in, n,dt,dt_out,t_end,part,a)
    

    !! compute the block i have to compute.
    !! from here on there is no distinction master / worker except for printing.
    !! all processors collaborate equally.
    !! it is assumed that n is divisible by p, and we compute the start of my
    !! block (my_start) and the end of my block (my_end)

    my_n = n / p
    my_start = (my_n * my_rank) + 1
    my_end = my_start + my_n - 1


    !! initialize head node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(head)
  call calculate_ranges(head)
  head%type = 0
  call nullify_pointers(head)

    
 
   !! create initial tree
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1, n
    call find_cell(head, temp_cell, part(i)%p)
    call place_cell(temp_cell, part(i)%p, i)
  end do
  call borrar_empty_leaves(head)
  call calculate_masses(head)

     !! compute initial accelerations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! here we start the time calculation

    t1 = mpi_wtime()

    a(my_start:my_end) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)


    call calculate_forces(head, my_start, my_end)
    call main_loop(part, a, n, dt, dt_out, t_end)

    if (my_rank == 0) then
        ! in order to print only once in the terminal
        print *, "Simulation complete :) Particles positions save in output.dat"
    end if

        
    t2 = mpi_wtime()

    if (my_rank == 0) then
        print *, "Total execution time (s): ", t2 - t1
    end if

    
    ! finally we check the status of the arrays and deallocate
    if (allocated(part)) deallocate(part, stat=ierr)
    if (allocated(a)) deallocate(a, stat=ierr)

    contains

    subroutine read_input(infile, in,  n, dt, dt_out, t_end, part, a)
        character(len=100) :: infile   
        character(len=200) :: in 
        integer, intent(out) :: n
        real(dp), intent(out) :: dt, dt_out, t_end
        type(particle3d), allocatable, intent(out) :: part(:)
        type(vector3d), allocatable, intent(out) :: a(:)
        integer :: i, ierr 

    
        open(unit=10, file=infile, status='old', action='read', iostat=ierr)
        if (my_rank == 0) then
            if (ierr /= 0) then 
                print *, "Error: the file doesn't exist or it cannot be open: ", trim(infile)
            stop
            endif
        endif

    
        !! matrix initialization
        !! the master reads from file and broadcasts
        !! all variables to all slaves
        !!
        if ( my_rank == 0 ) then

            read(10,*) dt
            read(10,*) dt_out
            read(10,*) t_end
            read(10,*) n
        
            call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, error)
            call mpi_bcast(dt, 1, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(dt_out, 1, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(t_end, 1, mpi_double_precision, 0, mpi_comm_world, error)
        
            allocate(part(n))
            allocate(a(n)) 

            do i = 1,n
                read(10,*) part(i)%m, & 
                part(i)%p%x, part(i)%p%y, part(i)%p%z, &
                part(i)%v%x, part(i)%v%y, part(i)%v%z
            end do

            close(10)


            call mpi_bcast(part%m, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%x, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%y, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%z, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%x, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%y, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%z, n, mpi_double_precision, 0, mpi_comm_world, error)
            
        else
            call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, error)
            call mpi_bcast(dt, 1, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(dt_out, 1, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(t_end, 1, mpi_double_precision, 0, mpi_comm_world, error)

            allocate(part(n))
            allocate(a(n))
        
            call mpi_bcast(part%m, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%x, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%y, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%p%z, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%x, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%y, n, mpi_double_precision, 0, mpi_comm_world, error)
            call mpi_bcast(part%v%z, n, mpi_double_precision, 0, mpi_comm_world, error)
        end if

   end subroutine read_input



end program ex2_mpi