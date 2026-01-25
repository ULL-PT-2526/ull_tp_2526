module barnes_hut_mpi_module

    use mpi
    use iso_fortran_env, only: real64
    use geometry
    use particle
    !$ use omp_lib
    implicit none
    private

    ! ------ PUBLIC SUBROUTINES (NOT NECESSARY, BUT A GOOD HABIT) ------

    public :: barnes_hut_initialize, barnes_hut_rebuild, barnes_hut_calculate_forces
    public :: barnes_hut_sync_particles
    public :: get_global_boundaries, find_cell, place_cell, calculate_masses
    public :: borrar_empty_leaves, borrar_tree

    ! ------ BARNES-HUT'S MAIN PARAMETER (OPENING CRITERION, THETA) ------
  
    real(real64), parameter, public :: theta = 1.0_real64 
    
    ! ------ DERIVED TYPES ------

    type, public :: range_type 
        type(point3d) :: min_corner
        type(point3d) :: max_corner
    end type range_type

    type, public :: cptr 
        type(cell), pointer :: ptr => null()
    end type cptr

    type, public :: cell 
        type(range_type) :: range 
        type(point3d) :: part 
        integer :: pos 
        integer :: type  
        real(real64) :: mass 
        type(point3d) :: c_o_m 
        type(cptr), dimension(2,2,2) :: subcell 
    end type cell

    ! ------ MODULE-LEVEL VARIABLES (GLOBAL TO THE MODULE) ------

    type(cell), pointer, public :: head, temp_cell  
    type(particle3d), dimension(:), allocatable, public :: particles
    integer, public :: n_particles

    ! --- COMMUNICATION BUFFERS ---
    real(real64), allocatable, private :: comm_sendbuf(:), comm_recvbuf(:)



contains



    ! =========================================================
    ! MAIN ROUTINES 
    ! =========================================================

    subroutine barnes_hut_initialize()
        integer :: i
        allocate(head)
        call get_global_boundaries(head)
        head%type = 0
        call nullify_pointers(head)

        do i = 1, n_particles
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head)
    end subroutine barnes_hut_initialize

    subroutine barnes_hut_rebuild()
        integer :: i
        call borrar_tree(head)
        
        call get_global_boundaries(head)
        head%type = 0
        call nullify_pointers(head)

        do i = 1, n_particles
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head)
    end subroutine barnes_hut_rebuild

    subroutine barnes_hut_calculate_forces(ax, ay, az, istart, iend)
        real(real64), dimension(:), intent(out) :: ax, ay, az
        integer, intent(in) :: istart, iend
        integer :: i
        type(vector3d) :: f_total

        ! Parallel region for force calculation (Hybrid support)
        !$omp parallel do default(none) &
        !$omp shared(ax, ay, az, particles, head, istart, iend) &
        !$omp private(i, f_total) &
        !$omp schedule(dynamic)
        do i = istart, iend
            f_total = vector3d(0.0_real64, 0.0_real64, 0.0_real64)
            call calculate_forces_aux(i, head, f_total)
            
            ax(i) = f_total%x / particles(i)%m
            ay(i) = f_total%y / particles(i)%m
            az(i) = f_total%z / particles(i)%m
        end do
        !$omp end parallel do
    end subroutine barnes_hut_calculate_forces

    subroutine barnes_hut_sync_particles(my_start, my_n, sendcounts, displs)
        integer, intent(in) :: my_start, my_n
        integer, dimension(:), intent(in) :: sendcounts, displs
        integer :: ierr

        ! Allocate buffers only once if not already allocated
        if (.not. allocated(comm_sendbuf)) then
             allocate(comm_sendbuf(my_n * 7))       
             allocate(comm_recvbuf(n_particles * 7)) 
        end if

        call pack_subset(my_start, my_n, comm_sendbuf)

        ! Encapsulated MPI communication
        call MPI_ALLGATHERV(comm_sendbuf, my_n*7, MPI_DOUBLE_PRECISION, &
                            comm_recvbuf, sendcounts*7, displs*7, MPI_DOUBLE_PRECISION, &
                            MPI_COMM_WORLD, ierr)

        call unpack_all(comm_recvbuf)
    end subroutine barnes_hut_sync_particles




    ! =========================================================
    ! ALGORITHM SUBROUTINES
    ! =========================================================

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Get Global Boundaries !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_global_boundaries(goal)
        type(cell), pointer :: goal
        type(point3d) :: mins, maxs, medios
        real(real64) :: span, half_span
        integer :: i
        
        mins = particles(1)%p
        maxs = particles(1)%p

        do i = 2, n_particles
            if (particles(i)%p%x < mins%x) mins%x = particles(i)%p%x
            if (particles(i)%p%x > maxs%x) maxs%x = particles(i)%p%x
            if (particles(i)%p%y < mins%y) mins%y = particles(i)%p%y
            if (particles(i)%p%y > maxs%y) maxs%y = particles(i)%p%y
            if (particles(i)%p%z < mins%z) mins%z = particles(i)%p%z
            if (particles(i)%p%z > maxs%z) maxs%z = particles(i)%p%z
        end do

        ! Add 10% margin so particles don't fall exactly on the edge
        span = max(maxs%x - mins%x, maxs%y - mins%y, maxs%z - mins%z) * 1.1_real64
        medios = (maxs + mins) * 0.5_real64
        
        ! Explicit constructor to avoid mixed-mode arithmetic errors
        half_span = span * 0.5_real64
        
        goal%range%min_corner = point3d(medios%x - half_span, &
                                        medios%y - half_span, &
                                        medios%z - half_span)
                                        
        goal%range%max_corner = point3d(medios%x + half_span, &
                                        medios%y + half_span, &
                                        medios%z + half_span)
    end subroutine get_global_boundaries

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Calculate Forces Aux !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_forces_aux(goal_idx, tree, force_accum)
        type(cell), pointer :: tree
        integer, intent(in) :: goal_idx
        type(vector3d), intent(inout) :: force_accum
        
        real(real64) :: d, r2, r3, l
        type(vector3d) :: rji, force_temp
        integer :: i, j, k
        real(real64), parameter :: epsilon = 0.01_real64

        select case (tree%type)
        case (1) ! Particle (Leaf Node)
            if (goal_idx /= tree%pos) then
                rji = particles(tree%pos)%p - particles(goal_idx)%p
                r2 = (rji%x**2 + rji%y**2 + rji%z**2)
                r3 = (r2 + epsilon**2) * sqrt(r2 + epsilon**2)
                
                force_temp = (particles(tree%pos)%m * particles(goal_idx)%m * rji) / r3
                force_accum = force_accum + force_temp
            end if
        case (2) ! Conglomerate (Internal Node)
            l = tree%range%max_corner%x - tree%range%min_corner%x
            rji = tree%c_o_m - particles(goal_idx)%p
            r2 = (rji%x**2 + rji%y**2 + rji%z**2)
            d = sqrt(r2)

            ! If conglomerate is far enough (l/d < theta), use approximation
            if (l/d < theta) then
                r3 = r2 * d
                force_temp = (tree%mass * particles(goal_idx)%m * rji) / r3
                force_accum = force_accum + force_temp
            else
                ! Otherwise, recurse into children
                do i = 1, 2
                    do j = 1, 2
                        do k = 1, 2
                            if (associated(tree%subcell(i, j, k)%ptr)) then
                                call calculate_forces_aux(goal_idx, tree%subcell(i, j, k)%ptr, force_accum)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine calculate_forces_aux

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Find Cell !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine find_cell(root, goal, part)
        type(point3d), intent(in) :: part
        type(cell), pointer, intent(in) :: root
        type(cell), pointer, intent(out) :: goal
        type(cell), pointer :: temp
        integer :: i, j, k

        select case (root%type)
        case (2)
            out: do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (belongs(part, root%subcell(i, j, k)%ptr)) then
                            call find_cell(root%subcell(i, j, k)%ptr, temp, part)
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

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Place Cell !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine place_cell(goal, part, n)
        type(cell), pointer :: goal, temp
        type(point3d), intent(in) :: part
        integer, intent(in) :: n

        select case (goal%type)
        case (0)
            goal%type = 1
            goal%part = part
            goal%pos = n
        case (1)
            call create_subcells(goal)
            call find_cell(goal, temp, part)
            call place_cell(temp, part, n)
        case (2)
             call find_cell(goal, temp, part)
             call place_cell(temp, part, n)
        end select
    end subroutine place_cell

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Create Subcells !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine create_subcells(goal)
        type(cell), pointer :: goal 
        type(point3d) :: part 
        integer :: i, j, k
        integer, dimension(3) :: octant 

        part = goal%part 
        goal%type = 2 

        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    octant = (/i, j, k/)
                    allocate(goal%subcell(i, j, k)%ptr)
                    
                    goal%subcell(i, j, k)%ptr%range%min_corner = get_subcell_limits(0, goal, octant)
                    goal%subcell(i, j, k)%ptr%range%max_corner = get_subcell_limits(1, goal, octant)

                    if (belongs(part, goal%subcell(i, j, k)%ptr)) then
                        goal%subcell(i, j, k)%ptr%type = 1 
                        goal%subcell(i, j, k)%ptr%part = part 
                        goal%subcell(i, j, k)%ptr%pos = goal%pos 
                    else 
                        goal%subcell(i, j, k)%ptr%type = 0 
                    end if
                    call nullify_pointers(goal%subcell(i, j, k)%ptr)
                end do
            end do
        end do
    end subroutine create_subcells

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Nullify Pointers !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine nullify_pointers(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    nullify(goal%subcell(i, j, k)%ptr)
                end do
            end do
        end do
    end subroutine nullify_pointers

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Belongs !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function belongs(part, goal)
        type(point3d), intent(in) :: part
        type(cell), pointer, intent(in) :: goal
        logical :: belongs

        if (part%x >= goal%range%min_corner%x .and. &
            part%x <= goal%range%max_corner%x .and. &
            part%y >= goal%range%min_corner%y .and. &
            part%y <= goal%range%max_corner%y .and. &
            part%z >= goal%range%min_corner%z .and. &
            part%z <= goal%range%max_corner%z) then
            belongs = .true.
        else
            belongs = .false.
        end if
    end function belongs

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Get Subcell Limits !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_subcell_limits(what, goal, octant)
        integer, intent(in) :: what
        type(cell), pointer, intent(in) :: goal
        integer, dimension(3), intent(in) :: octant 
        type(point3d) :: get_subcell_limits, valor_medio

        valor_medio = (goal%range%min_corner + goal%range%max_corner) * 0.5_real64

        select case (what)
        case (0) ! Min corner
            if (octant(1) == 1) then
                get_subcell_limits%x = goal%range%min_corner%x
            else
                get_subcell_limits%x = valor_medio%x
            end if
            
            if (octant(2) == 1) then
                get_subcell_limits%y = goal%range%min_corner%y
            else
                get_subcell_limits%y = valor_medio%y
            end if

            if (octant(3) == 1) then
                get_subcell_limits%z = goal%range%min_corner%z
            else
                get_subcell_limits%z = valor_medio%z
            end if

        case (1) ! Max corner
            if (octant(1) == 1) then
                get_subcell_limits%x = valor_medio%x
            else
                get_subcell_limits%x = goal%range%max_corner%x
            end if

            if (octant(2) == 1) then
                get_subcell_limits%y = valor_medio%y
            else
                get_subcell_limits%y = goal%range%max_corner%y
            end if

            if (octant(3) == 1) then
                get_subcell_limits%z = valor_medio%z
            else
                get_subcell_limits%z = goal%range%max_corner%z
            end if
        end select
    end function get_subcell_limits

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Calculate Masses !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine calculate_masses(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        real(real64) :: mass, total_mass

        goal%mass = 0.0_real64                                  
        goal%c_o_m = point3d(0.0_real64, 0.0_real64, 0.0_real64)

        select case (goal%type)
        case (1) 
            goal%mass = particles(goal%pos)%m
            goal%c_o_m = particles(goal%pos)%p
        case (2) 
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (associated(goal%subcell(i, j, k)%ptr)) then 
                            call calculate_masses(goal%subcell(i, j, k)%ptr) 
                            mass = goal%mass 
                            total_mass = mass + goal%subcell(i, j, k)%ptr%mass 
                            
                            goal%c_o_m = (goal%c_o_m * mass + &
                                          goal%subcell(i, j, k)%ptr%c_o_m * goal%subcell(i, j, k)%ptr%mass) / total_mass 
                            goal%mass = total_mass 
                        end if
                    end do
                end do
            end do
        end select
    end subroutine calculate_masses

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Erase Empty Leaves !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_empty_leaves(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        if (associated(goal%subcell(1, 1, 1)%ptr)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call borrar_empty_leaves(goal%subcell(i, j, k)%ptr)
                        if (goal%subcell(i, j, k)%ptr%type == 0) then
                            deallocate(goal%subcell(i, j, k)%ptr)
                            nullify(goal%subcell(i, j, k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine borrar_empty_leaves

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !! Erase Tree !!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive subroutine borrar_tree(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    if (associated(goal%subcell(i, j, k)%ptr)) then
                        call borrar_tree(goal%subcell(i, j, k)%ptr)
                        deallocate(goal%subcell(i, j, k)%ptr)
                        nullify(goal%subcell(i, j, k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine borrar_tree

    ! --- Private helpers for communication ---
    subroutine pack_subset(start_idx, count, buf)
        integer, intent(in) :: start_idx, count
        real(real64), intent(out) :: buf(:)
        integer :: j, p_idx
        do j = 1, count
            p_idx = start_idx + j - 1
            buf((j-1)*7 + 1) = particles(p_idx)%p%x
            buf((j-1)*7 + 2) = particles(p_idx)%p%y
            buf((j-1)*7 + 3) = particles(p_idx)%p%z
            buf((j-1)*7 + 4) = particles(p_idx)%v%x
            buf((j-1)*7 + 5) = particles(p_idx)%v%y
            buf((j-1)*7 + 6) = particles(p_idx)%v%z
            buf((j-1)*7 + 7) = particles(p_idx)%m
        end do
    end subroutine pack_subset

    subroutine unpack_all(buf)
        real(real64), intent(in) :: buf(:)
        integer :: j
        do j = 1, n_particles
            particles(j)%p%x = buf((j-1)*7 + 1)
            particles(j)%p%y = buf((j-1)*7 + 2)
            particles(j)%p%z = buf((j-1)*7 + 3)
            particles(j)%v%x = buf((j-1)*7 + 4)
            particles(j)%v%y = buf((j-1)*7 + 5)
            particles(j)%v%z = buf((j-1)*7 + 6)
            particles(j)%m   = buf((j-1)*7 + 7)
        end do
    end subroutine unpack_all

end module barnes_hut_mpi_module