module tree

    use geometry
    use particles
    !$ use omp_lib 
    ! To avoid compatibility issues

    implicit none

    integer :: i,j,k,n
    real(dp) :: dt, t_end, t, dt_out, t_out
    real(dp), parameter :: theta = 1
    type(particle3d), dimension(:), allocatable :: part
    type(vector3d), dimension(:), allocatable :: a
    
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
    !! Calculate_Ranges !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calculates the ranges of the particles in the
    !! array r in the 3 dimensions and stores them in
    !! the variable pointed to by goal
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine calculate_ranges(goal)
        type(cell),pointer :: goal
        real(dp), dimension(3) :: mins,maxs,medios
        real(dp) :: span
        real(dp), allocatable :: xs(:), ys(:), zs(:)

        xs = [ (part(i)%p%x, i=1,n) ]
        ys = [ (part(i)%p%y, i=1,n) ]
        zs = [ (part(i)%p%z, i=1,n) ]

        mins = [ minval(xs), minval(ys), minval(zs) ]
        maxs = [ maxval(xs), maxval(ys), maxval(zs) ]

        ! When computing span add 10% so particles don't land exactly on the boundary
        span = maxval(maxs - mins) * 1.1
        medios = (maxs + mins) / 2.0
        goal%range%min = medios - span/2.0
        goal%range%max = medios + span/2.0

    end subroutine calculate_ranges

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Find_Cell !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Finds the cell where we will place the particle.
    !! If the cell we are considering has no particle or has a single particle,
    !! that is the cell where we will place the particle.
    !! If the cell we are considering is a "group", we use the function BELONGS
    !! to determine which of the 8 subcells the particle belongs to and then
    !! call Find_Cell recursively on that subcell.
    !!
    !! NOTE: When a "group" cell is created its 8 subcells are allocated,
    !! so we can assume all 8 exist. Empty cells are deleted at the end,
    !! after the whole tree has been built.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine find_cell(root,goal,partp) ! In the input we put part%p, the part position that is a point3d
        type(point3d) :: partp
        type(cell),pointer :: root,goal,temp
        integer :: i,j,k
        select case (root%type)
        case (2)
            ! Here we dont use opnmp since we have exit in the do loop and to avoid cell assignation problems
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
    !! Place_Cell !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Executed after Find_Cell, in the cell returned by that function,
    !! so it is always a cell of type 0 (empty) or type 1 (contains a particle).
    !! If the cell is type 1 we must subdivide the cell and place both particles
    !! (the existing one and the new one).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine place_cell(goal,partp,n) ! Input is part%p, the position of the particle
        type(cell),pointer :: goal,temp
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
            print*,"SHOULD NOT BE HERE. ERROR!"
        end select
    end subroutine place_cell

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Crear_Subcells !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! This routine is called from Place_Cell and is only called when there
    !! is already a particle in the cell, so we need to subdivide it. It creates
    !! 8 subcells hanging from goal and places the particle that was in goal
    !! into the appropriate one of the 8 new subcells.
    !!
    !! To create the subcells use the routines CALCULAR_RANGE, BELONGS and NULLIFY_POINTERS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine crear_subcells(goal)
        type(cell), pointer :: goal
        type(particle3d) :: part
        integer :: i,j,k,n
        integer, dimension(3) :: octant
        part = goal%part
        goal%type=2
        ! We do not use openmp to construct the tree since it is tricky with the pointers
        do i = 1,2
            do j = 1,2
                do k = 1,2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = calcular_range (0,goal,octant)
                    goal%subcell(i,j,k)%ptr%range%max = calcular_range (1,goal,octant)
                    if (belongs(part%p,goal%subcell(i,j,k)%ptr)) then
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
    !! Nullify_Pointers !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Simply nullifies the pointers of the 8 subcells of cell "goal"
    !!
    !! Used in the main loop and by CREAR_SUBCELLS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nullify_pointers(goal)
        type(cell), pointer :: goal
        integer :: i,j,k
        do i = 1,2
            do j = 1,2
                do k = 1,2
                    nullify(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine nullify_pointers

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Belongs !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Returns TRUE if particle "part" is inside the range of cell "goal"
    !!
    !! Used by FIND_CELL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function belongs (partp,goal) ! The input is part%p, the position of the particle
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
    !! Calcular_Range !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Given an octant "octant" (1,1,1, 1,1,2 ... 2,2,2),
    !! computes its ranges based on the ranges of "goal".
    !! If "what" = 0 compute the minimums. If what=1 compute the maximums.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function calcular_range (what,goal,octant)
        integer :: what,n
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
    !! Borrar_empty_leaves !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Called once the tree is complete to delete (DEALLOCATE)
    !! the empty cells (i.e. without a particle).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_empty_leaves(goal)
        type(cell),pointer :: goal
        integer :: i,j,k
        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1,2
                do j = 1,2
                    do k = 1,2
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
    !! Borrar_tree !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Deletes the entire tree, except the "head".
    !!
    !! The tree must be regenerated continuously, so we need to delete the old one
    !! to avoid memory leaks.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_tree(goal)
        type(cell),pointer :: goal
        integer :: i,j,k
        do i = 1,2
            do j = 1,2
                do k = 1,2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate (goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine borrar_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_masses !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Computes for all cells hanging from "goal" their mass and center-of-mass.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_masses(goal)
        type(cell),pointer :: goal
        integer :: i,j,k
        real(dp), dimension(3) :: c_o_m
        real(dp) :: m_temp

        goal%part%m = 0.0_dp
        goal%c_o_m = point3d(0.0_dp, 0.0_dp, 0.0_dp)

        select case (goal%type)
        case (1)
            goal%part%m = part(goal%pos)%m
            goal%c_o_m = part(goal%pos)%p
        case (2)
            do i = 1,2
                do j = 1,2
                    do k = 1,2
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
    !! Calculate_forces !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calculates the forces on all particles from "head".
    !! Uses the function Calculate_forces_aux which performs the actual
    !! calculations for each particle.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine calculate_forces(head)
        type(cell),pointer :: head
        integer :: i,j,k,start,end
        !schedule(dynamic, 10) assigns particles in small chunks to threads as they become available, keeping 
        !all cores busy and preventing to wait for others threads.
        !$omp parallel do private(i) schedule(dynamic, 10) 
        do i = 1,n
            call calculate_forces_aux(i,head)
        end do
        !$omp end parallel do
    end subroutine calculate_forces

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces_aux !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Given a particle "goal" computes the forces on it from the cell "tree".
    !! If "tree" is a cell that contains a single particle the case is simple,
    !! since it is just two particles interacting.
    !!
    !! If "tree" is a group cell, first check whether l/D < theta. That is,
    !! if the side of the cell (l) divided by the distance from particle goal
    !! to the center_of_mass of cell tree (D) is less than theta.
    !! If so, treat the cell as a single particle. If not, consider all the
    !! subcells of tree and call Calculate_forces_aux recursively for each.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    recursive subroutine calculate_forces_aux(goal,tree)
        type(cell),pointer :: tree
        integer :: i,j,k,goal
        real(dp) :: l, rs, r3
        type(vector3d) :: rji
        select case (tree%type)
        case (1)
            if (goal .ne. tree%pos) then
                rji = tree%c_o_m - part(goal)%p
                rs = distance(tree%c_o_m, part(goal)%p)
                r3 = rs**3
                a(goal) = a(goal) + part(tree%pos)%m * rji / r3
            end if
        case (2)
            !! The range has the same span in the 3 dimensions
            !! so we can consider any dimension to compute the cell side
            !! (here dimension 1)
            l = tree%range%max(1) - tree%range%min(1)
            rji = tree%c_o_m - part(goal)%p
            rs = distance(tree%c_o_m, part(goal)%p)
            if (l/rs < theta) then
                !! If a group cell, check if l/D < theta
                r3 = rs**3
                a(goal) = a(goal) + tree%part%m * rji / r3
            else
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (associated(tree%subcell(i,j,k)%ptr)) then
                                call calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr)
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

        ! We create the output file
        open(unit=12, file='output.dat', status='replace', action='write')

        t_out = 0.0

        ! If we change the number of particles the code could be slower, therefore it is convenient to write messages in order to clarify that all is 
        ! running smoothly

        print *, "Starting the simulation..."
        print *, "Number of particles:", n
        print *, "Total time:", t_end, "Timestep:", dt
        print *, "---------------------------------------------"


        t_out = 0.0
        do t = 0.0, t_end, dt
            
           ! With the following we print the time every 1.0 time units
            if (mod(t, 1.0_dp) < dt) then
                 print *, "t = ", t
            end if

            !$omp parallel do private(i)
            do i= 1, n
                part(i)%v = part(i)%v + a(i) * dt/2.0_dp
                part(i)%p = part(i)%p + part(i)%v * dt
            end do
            !$omp end parallel do
            !! Positions have changed, so we must delete and reinitialize the tree
            call borrar_tree(head)
            call calculate_ranges(head)
            head%type = 0
            call nullify_pointers(head)
            do i = 1,n
                call find_cell(head,temp_cell,part(i)%p)
                call place_cell(temp_cell,part(i)%p,i)
            end do
            call borrar_empty_leaves(head)
            call calculate_masses(head)
            !$omp parallel do private(i)
            do i = 1, n
                a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
            end do
            !$omp end parallel do
            call calculate_forces(head)
            !$omp parallel do private(i)
            do i= 1, n
                part(i)%v = part(i)%v + a(i) * dt/2.0_dp
            end do
            !$omp end parallel do
                   ! Now we write in the output data the position of the particles each time
            t_out = t_out + dt
            if (t_out >= dt_out) then
                write(12,'(ES20.10)', advance ='no') t ! With, "advance=no" we avoid to make a split of line
                do i = 1,n
                    write(12,'(3(1X,ES20.10))', advance='no') part(i)%p
                end do
                write(12,*) 
                t_out = 0.0
            end if
        end do

        close(12)
    end subroutine main_loop


end module tree