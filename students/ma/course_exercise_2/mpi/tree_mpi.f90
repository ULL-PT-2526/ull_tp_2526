
! we choose the first option of the notes
! since i tried with the second option and the code was slower because of the communnications
!between the processes

module tree_mpi
    use MPI
    use geometry
    use particles

    implicit none

    !! mpi variables (rank, number of processors, etc.)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: my_rank, p, error
    integer, dimension(mpi_status_size) :: status
    integer :: my_n, my_start, my_end
    
    !! n-body problem variables
    !!
 
    integer :: i,j,k,n
    real(dp) :: dt, t_end, t, dt_out, t_out
    real(dp), parameter :: theta = 1
    type(particle3d), dimension(:), allocatable :: part
    type(vector3d), dimension(:), allocatable :: a, total_a

    type range ! limits of a cell
        real(dp), dimension(3) :: min,max
    end type range
    
    type cptr
        type(cell), pointer :: ptr
    end type cptr
    
    type cell
        type (range) :: range
        type(particle3d) :: part
        integer :: pos
        integer :: type !! 0 = no particule; 1 = particle; 2 = group of particles
        type(point3d) :: c_o_m !centre of mass
        type (cptr), dimension(2,2,2) :: subcell
    end type cell
    
    type (cell), pointer :: head, temp_cell

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate_ranges
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculate_ranges(goal)
        type(cell),pointer :: goal
        real(dp), dimension(3) :: mins,maxs,center
        real(dp) :: span
        real(dp), allocatable :: xs(:), ys(:), zs(:)

        xs = [ (part(i)%p%x, i=1,n) ]
        ys = [ (part(i)%p%y, i=1,n) ]
        zs = [ (part(i)%p%z, i=1,n) ]

        mins = [ minval(xs), minval(ys), minval(zs) ]
        maxs = [ maxval(xs), maxval(ys), maxval(zs) ]

        ! when computing span add 10% so particles don't land exactly on the boundary
        span = maxval(maxs - mins) * 1.1
        center = (maxs + mins) / 2.0
        goal%range%min = center - span/2.0
        goal%range%max = center + span/2.0

    end subroutine calculate_ranges
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! find_cell
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine find_cell(root, goal, partp)
        type(point3d) :: partp
        type(cell),pointer :: root,goal,temp
        integer :: i,j,k
        select case (root%type)
        case (2)
            out: do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (belongs(partp,root%subcell(i,j,k)%ptr)) then
                            call find_cell(root%subcell(i,j,k)%ptr,temp,partp)
                            goal => temp
                            exit out
                        end if
                    end do
                end do
            end do out
        case default
            goal => root
        end select
    end subroutine find_cell
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! place_cell
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine place_cell(goal, partp, n)
        type(cell), pointer :: goal, temp
        type(point3d) :: partp
        integer :: n
        select case (goal%type)
        case (0)
            goal%type = 1
            goal%part%p = partp
            goal%pos = n
        case (1)
            call crear_subcells(goal)
            call find_cell(goal,temp,partp)
            call place_cell(temp,partp,n)
        case default
            print*,"Should not be here. Error!"
        end select
    end subroutine place_cell
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! subcells
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine crear_subcells(goal)
        type(cell), pointer :: goal
        type(particle3d) :: part
        integer :: i, j, k, n
        integer, dimension(3) :: octant
        
        part = goal%part
        goal%type = 2
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = calcular_range(0, goal, octant)
                    goal%subcell(i,j,k)%ptr%range%max = calcular_range(1, goal, octant)
                    
                    if (belongs(part%p, goal%subcell(i,j,k)%ptr)) then
                        goal%subcell(i,j,k)%ptr%part%p = part%p
                        goal%subcell(i,j,k)%ptr%type = 1
                        goal%subcell(i,j,k)%ptr%pos = goal%pos
                    else
                        goal%subcell(i,j,k)%ptr%type = 0
                    end if
                    call nullify_pointers(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine crear_subcells
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! nullify_pointers
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nullify_pointers(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    nullify(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine nullify_pointers
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! belongs
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function belongs (partp, goal)
        type(point3d) :: partp
        type(cell), pointer :: goal
        logical :: belongs
        
        if (partp%x >= goal%range%min(1) .and. &
            partp%x <= goal%range%max(1) .and. &
            partp%y >= goal%range%min(2) .and. &
            partp%y <= goal%range%max(2) .and. &
            partp%z >= goal%range%min(3) .and. &
            partp%z <= goal%range%max(3)) then
            belongs = .true.
        else
            belongs = .false.
        end if
    end function belongs
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calcular_range
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function calcular_range (what, goal, octant)
        integer :: what, n
        type(cell), pointer :: goal
        integer, dimension(3) :: octant
        real(dp), dimension(3) :: calcular_range, valor_medio
        
        valor_medio = (goal%range%min + goal%range%max) / 2.0
        
        select case (what)
        case (0)
            where (octant == 1)
                calcular_range = goal%range%min
            elsewhere
                calcular_range = valor_medio
            endwhere
        case (1)
            where (octant == 1)
                calcular_range = valor_medio
            elsewhere
                calcular_range = goal%range%max
            endwhere
        end select
    end function calcular_range
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! borrar_empty_leaves
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_empty_leaves(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        
        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate (goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine borrar_empty_leaves
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! borrar_tree
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_tree(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate (goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine borrar_tree
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate_masses
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_masses(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        real(dp), dimension(3) :: c_o_m
        real(dp) :: m_temp
        
        goal%part%m = 0.0_dp
        goal%c_o_m = point3d(0.0_dp, 0.0_dp, 0.0_dp)

        select case (goal%type)
        case (1)
            goal%part%m = part(goal%pos)%m
            goal%c_o_m = part(goal%pos)%p
        case (2)
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (associated(goal%subcell(i,j,k)%ptr)) then
                            call calculate_masses(goal%subcell(i,j,k)%ptr)
                            m_temp = goal%part%m
                            goal%part%m = m_temp + goal%subcell(i,j,k)%ptr%part%m
                            goal%c_o_m = (m_temp * goal%c_o_m + &
                                goal%subcell(i,j,k)%ptr%part%m * goal%subcell(i,j,k)%ptr%c_o_m) / goal%part%m
                        end if
                    end do
                end do
            end do
        end select
    end subroutine calculate_masses
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate_forces
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_forces(head, start, end)
        type(cell), pointer :: head
        integer :: i, j, k, start, end
        do i = start, end
            call calculate_forces_aux(i, head)
        end do
    end subroutine calculate_forces
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! calculate_forces_aux
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_forces_aux(number, tree)
        type(cell), pointer :: tree
        integer :: i, j, k, number
        real(dp) :: l, rs, r3
        type(vector3d) :: rji
        
        select case (tree%type)
        case (1)
            if (number .ne. tree%pos) then
               rji = tree%c_o_m - part(number)%p
                rs = distance(tree%c_o_m, part(number)%p)
                r3 = rs**3
                a(number) = a(number) + part(tree%pos)%m * rji / r3
            end if
        case (2)
            l = tree%range%max(1) - tree%range%min(1) !! range has the same span in 3 dimensions
            rji = tree%c_o_m - part(number)%p
            rs = distance(tree%c_o_m, part(number)%p)
            
            if (l/rs < theta) then
                !! if cluster, we need to check if l/d < theta
                r3 = rs**3
                a(number) = a(number) + tree%part%m * rji / r3
            else
                do i = 1, 2
                    do j = 1, 2
                        do k = 1, 2
                            if (associated(tree%subcell(i,j,k)%ptr)) then
                                call calculate_forces_aux(number, tree%subcell(i,j,k)%ptr)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine calculate_forces_aux

    subroutine main_loop(part, a, n, dt, dt_out, t_end)
        type(particle3d), intent(inout) :: part(:)
        type(vector3d), intent(inout) :: a(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: dt, dt_out, t_end
        real(dp) :: t, t_out
        integer :: i
        
        ! we create the output file
        open(unit=12, file='output.dat', status='replace', action='write')

        t_out = 0.0

        ! if we change the number of particles the code could be slower, therefore it is convenient to write messages in order to clarify that all is 
        ! running smoothly
        if (my_rank == 0) then
            print *, "Starting the simulation..."
            print *, "Number of particles:", n
            print *, "Total time:", t_end, "Timestep:", dt
            print *, "---------------------------------------------"
        end if

        do t = 0.0, t_end, dt
            
            if (my_rank == 0) then
            ! with the following we print the time every 1.0 time units
                if (mod(t, 1.0_dp) < dt) then
                    print *, "t = ", t
                end if
            end if
    
            part(my_start:my_end)%v = part(my_start:my_end)%v + a(my_start:my_end) * dt/2.0_dp
            part(my_start:my_end)%p = part(my_start:my_end)%p + part(my_start:my_end)%v * dt

            ! We split the lines since Fortran has a limit of 132 characters per line
            call mpi_allgather( & 
                                part(my_start)%p%x, my_n, mpi_double_precision, &
                                part(1)%p%x, my_n, mpi_double_precision, &
                                mpi_comm_world, error )

            call mpi_allgather( &
                                part(my_start)%p%y, my_n, mpi_double_precision, &
                                part(1)%p%y, my_n, mpi_double_precision, &
                                mpi_comm_world, error )

            call mpi_allgather( &
                                part(my_start)%p%z, my_n, mpi_double_precision, &
                                part(1)%p%z, my_n, mpi_double_precision, &
                                mpi_comm_world, error )
        
            !! positions have changed, so we need to delete
            !! and reinitialize the tree
            call borrar_tree(head)
            call calculate_ranges(head)
            head%type = 0
            call nullify_pointers(head)
        
            do i = 1, n
                call find_cell(head, temp_cell, part(i)%p)
                call place_cell(temp_cell, part(i)%p, i)
            end do
        
            call borrar_empty_leaves(head)
            call calculate_masses(head)

            a(my_start:my_end) = vector3d(0.0_dp, 0.0_dp, 0.0_dp) 

            call calculate_forces(head, my_start, my_end)
            part(my_start:my_end)%v = part(my_start:my_end)%v + a(my_start:my_end) * dt/2.0_dp
        
            !! we only print if we are the master
            
            if (my_rank == 0) then
                t_out = t_out + dt
                if (t_out >= dt_out) then
                    write(12,'(es20.10)', advance ='no') t ! with, "advance=no" we avoid to make a split of line
                    do i = 1,n
                        write(12,'(3(1x,es20.10))', advance='no') part(i)%p
                    end do
                    write(12,*) 
                    t_out = 0.0
                end if
            end if
        end do
    
        call mpi_finalize (error)

    end subroutine main_loop

end module tree_mpi