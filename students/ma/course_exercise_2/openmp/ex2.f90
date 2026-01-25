program ex2

use geometry
    use particles
    use tree
    !$ use omp_lib

    implicit none

    integer :: ierr ! This is to check and deallocate
    character(len=100) :: infile   ! Name of the input file
    character(len=200) :: in ! argument in the terminal for the input file
    real(dp) :: t1, t2
    ! Variables for system_clock (alternative to omp_get_wtime for serial running)
    integer(kind=8) :: count_start, count_end, count_rate

    call get_command_argument(1,in)

    if (len_trim(in)==0) then ! If input is not given 
        print *, "Input example: ./ex2 <inputfile>"
        stop
    endif

    infile = trim(in) ! remove empty space
    print *, "Reading input file: ", trim(infile)


    ! We use the subroutine to read the file
    call read_input(infile,in, n,dt,dt_out,t_end,part,a)
    
    !! Head node initialization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(head)
    call calculate_ranges(head)
    head%type = 0
    call nullify_pointers(head)
    
    !! Initial tree creation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,n
        call find_cell(head,temp_cell,part(i)%p)
        call place_cell(temp_cell,part(i)%p,i)
    end do

    call borrar_empty_leaves(head)
    call calculate_masses(head)
    
    !! Compute initial accelerations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Here we start the time calculation

    ! We initialize with -1.0 to know if OpenMP is active
    t1 = -1.0_dp 
    t2 = -1.0_dp

    !$ t1 = omp_get_wtime() ! Only if we compile with -fopenmp
    
    ! If we do not use OpenMP:
    if (t1 == -1.0_dp) then 
        call system_clock(count_start, count_rate)
        print *, "Running in serial mode (using system_clock)..."
    else
        print *, "Running in parallel mode (using OpenMP wtime)..."
    end if


    do i = 1, n
      a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
    end do
    
    call calculate_forces(head)

    call main_loop(part, a, n, dt, dt_out, t_end)

    !$ t2 = omp_get_wtime()

    print *, "Simulation complete :) Particles positions save in output.dat"

    ! If t1 changed, we used OpenMP 
    if (t1 /= -1.0_dp) then
        print *, "Total execution time (OPMP): ", t2 - t1, " s"
    else
        call system_clock(count_end)
        print *, "Total execution time (serial): ", &
             real(count_end - count_start, dp) / real(count_rate, dp), " s"
    end if

    ! Finally we check the status of the arrays and deallocate
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
        if (ierr /= 0) then 
            print *, "Error: The file doesn't exist or it cannot be open: ", trim(infile)
        stop
        endif

  
        read(10,*) dt
        read(10,*) dt_out
        read(10,*) t_end
        read(10,*) n

    ! We are going to associate the components of the particle type for simplicity, and now we read the last lines, n times to allocate n particles

        allocate(part(n))
        allocate(a(n)) 

        do i = 1,n
        read(10,*) part(i)%m, & 
                part(i)%p%x, part(i)%p%y, part(i)%p%z, &
                part(i)%v%x, part(i)%v%y, part(i)%v%z
        end do

        close(10)

        do i = 1, n
            a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        end do

   end subroutine read_input



end program ex2