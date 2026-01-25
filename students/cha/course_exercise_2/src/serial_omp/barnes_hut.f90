module barnes_hut_module

    use iso_fortran_env, only: real64
    use geometry
    use particle
    !$ use omp_lib
    implicit none
    private

    ! ------ PUBLIC SUBROUTINES (NOT NECESSARY, BUT A GOOD HABIT) ------

    public :: barnes_hut_initialize, barnes_hut_rebuild, barnes_hut_calculate_forces
    public :: get_global_boundaries, find_cell, place_cell, calculate_masses
    public :: borrar_empty_leaves, borrar_tree


    ! ------ BARNES-HUT'S MAIN PARAMETER (OPENING CRITERION, THETA) ------
    ! if a cell is sufficiently far from a particle (cell-size/distance < theta),
    ! one can treat it as a single point mass

    real(real64), parameter, public :: theta = 1.0_real64   ! theta = 1 is a standard 'safe' value
                                                            ! theta = 0 recovers the N-body case
                                                            ! smaller theta = more accurate but slower
    ! ------ DERIVED TYPES ------

    type, public :: range_type   ! defines a 3D cubic region using two opposite corners
        type(point3d) :: min_corner
        type(point3d) :: max_corner
    end type range_type

    type, public :: cptr    ! pointer to cell
        type(cell), pointer :: ptr => null()
    end type cptr

    type, public :: cell    ! octree node
        type(range_type) :: range   ! spatial region covered by this cell
        type(point3d) :: part       ! particle position (if type = 1)
        integer :: pos              ! particle index within the global array
        integer :: type             ! 0 = no particle, 1 = particle, 2 = multiple particles subdivided into 8 subcells
        real(real64) :: mass        ! total mass of the particles inside a cell
        type(point3d) :: c_o_m      ! center of mass of the cell
        type(cptr), dimension(2,2,2) :: subcell    ! 8 "children" (octree <-> 2x2x2)
    end type cell

    ! ------ MODULE-LEVEL VARIABLES (GLOBAL TO THE MODULE) ------
    ! head: pointer to the tree root (entire space)
    ! temp_cell: temporary cell used in recursive searches
    ! particles: dynamic array containing all particles in the system
    ! n_particles: total number of particles

    type(cell), pointer, public :: head, temp_cell  
    type(particle3d), dimension(:), allocatable, public :: particles
    integer, public :: n_particles



! ##########################################################################################




contains

    !==================================
    ! MAIN BARNES-HUT DRIVERS
    !==================================

    ! ------ barnes_hut_initialize ------
    ! creates the Barnes-Hut tree from scratch. Called once at the beginning of the simulation

    subroutine barnes_hut_initialize()
        integer :: i

        ! initialize head node
        allocate(head)
        call get_global_boundaries(head)   ! computes the size of the box that contains all particles
        head%type = 0                      ! starts empty

        ! initial tree construction
        do i = 1, n_particles
            call find_cell(head, temp_cell, particles(i)%p) ! finds where this particle should go in the tree
            call place_cell(temp_cell, particles(i)%p, i)   ! inserts the particle at that position
                                                            ! if the cell was already occupied, it splits into 8 children
        end do

        call borrar_empty_leaves(head)  ! removes empty leaves
        call calculate_masses(head)     ! traverses the tree bottom-up
                                        ! computing total mass and COM of each cell

    end subroutine barnes_hut_initialize



    ! ------ barnes_hut_rebuild ------
    ! called at each timestep. Since particles move, the old tree is no longer valid and must be rebuilt

    subroutine barnes_hut_rebuild()
        integer :: i

        call borrar_tree(head)               ! frees memory of the entire old tree but keeps the head structure
        call get_global_boundaries(head)     ! recomputes root box boundaries

        ! rest is identical to 'barnes_hut_initialize'
        head%type = 0

        do i = 1, n_particles
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head)

    end subroutine barnes_hut_rebuild



    ! ------ barnes_hut_calculate_forces ------
    ! computes the acceleration experienced by each particle due to the gravity of all others

    subroutine barnes_hut_calculate_forces(ax, ay, az)
        real(real64), dimension(:), intent(out) :: ax, ay, az  ! acceleration arrays
        integer :: i
        type(vector3d) :: f_total   ! temporary vector accumulating all forces acting on i

        !====== PARALELLISATION ======
        !$omp parallel do default(none) &
        !$omp shared(ax, ay, az, particles, n_particles, head) &
        !$omp private(i, f_total) &
        !$omp schedule(dynamic)
        do i = 1, n_particles
            f_total = vector3d(0.0_real64, 0.0_real64, 0.0_real64)  ! before starting for a new particle
                                                                    ! clear the accumulator
            ! compute total force on particle i
            call calculate_forces_aux(i, head, f_total)    ! function defined later
            
            ! convert force to acceleration: a = F / m
            ax(i) = f_total%x / particles(i)%m
            ay(i) = f_total%y / particles(i)%m
            az(i) = f_total%z / particles(i)%m
        end do
        !$omp end parallel do

    end subroutine barnes_hut_calculate_forces




    !=============================================================
    ! TREE MANAGEMENT (CREATE/MODIFY/DESTROY CELLS)
    !=============================================================

        ! ------ place_cell ------
    ! decides what to do when attempting to insert a particle into a specific 'goal' cell

    recursive subroutine place_cell(goal, part, n)
        type(cell), pointer :: goal, temp
        type(point3d), intent(in) :: part
        integer, intent(in) :: n

        select case (goal%type)  ! placement depends on the current case

        case (0)  ! "leaf" case with nothing inside
            goal%type = 1      ! change cell type to "occupied"
            goal%part = part   ! store particle coordinates
            goal%pos = n       ! and its index

        case (1)  ! "leaf" case that already had a particle inside
            
            ! split the cell into 8 children.
            ! IMPORTANT: 'create_subcells' automatically moves the existing (old) 
            ! particle into the correct child. We don't need to do it manually
            call create_subcells(goal) 
            
            ! 2. Now we just insert the NEW particle into the correct child.
            call find_cell(goal, temp, part)          ! find destination for new particle
            call place_cell(temp, part, n)            ! insert new particle
        
        case (2) ! "branch" case (already split)
             ! if we land here directly (rare but possible via recursion), just pass it down
             call find_cell(goal, temp, part)
             call place_cell(temp, part, n)

        end select

    end subroutine place_cell



    ! ------ create_subcells ------
    ! this is our "cellular mitosis": it takes a cell with a single particle (parent) and divides it into 8 children,
    ! relocating the parent particle into the corresponding child

    subroutine create_subcells(goal)
        type(cell), pointer :: goal        ! node to be split
        type(point3d) :: part              ! temporary variable to store the parent's particle
        integer :: i, j, k
        integer, dimension(3) :: octant    ! to identify the child
        logical :: placed                  ! flag to check if the particle has already been placed

        part = goal%part   ! safety copy
        goal%type = 2      ! change goal from type 1 to type 2
        placed = .false.   ! initialize flag as false (not placed yet)

        ! creation loop
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2

                    octant = (/i, j, k/)
                    allocate(goal%subcell(i, j, k)%ptr)

                    ! compute octant boundaries using 'get_subcell_limits'
                    goal%subcell(i, j, k)%ptr%range%min_corner = get_subcell_limits(0, goal, octant)
                    goal%subcell(i, j, k)%ptr%range%max_corner = get_subcell_limits(1, goal, octant)

                    ! only place if it belongs and hasn't been placed yet
                    if (.not. placed .and. belongs(part, goal%subcell(i, j, k)%ptr)) then
                        goal%subcell(i, j, k)%ptr%type = 1        ! child is born occupied
                        goal%subcell(i, j, k)%ptr%part = part     ! insert particle
                        goal%subcell(i, j, k)%ptr%pos = goal%pos  ! keep original properties
                        placed = .true.                           ! mark as placed

                    else  ! if it is not its "home"
                        goal%subcell(i, j, k)%ptr%type = 0        ! child is born empty
                    end if
                                                            
                end do
            end do
        end do

    end subroutine create_subcells



    ! ------ find_cell ------
    ! finds the cell in which a particle is located
    ! navigates from the root downwards until reaching the deepest node containing the particle position

    recursive subroutine find_cell(root, goal, part)
        type(point3d), intent(in) :: part
        type(cell), pointer, intent(in) :: root
        type(cell), pointer, intent(out) :: goal
        type(cell), pointer :: temp
        integer :: i, j, k

        select case (root%type)  ! behavior depends on the current cell type

        case (2) ! if root%type is 2, the cell is already subdivided into 8 children
                 ! the particle cannot stay here; it must descend into one child
            out: do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        ! the 'belongs' function is defined later
                        if (belongs(part, root%subcell(i, j, k)%ptr)) then
                            call find_cell(root%subcell(i, j, k)%ptr, temp, part)
                            goal => temp
                            exit out  ! once found, stop searching
                        end if
                    end do
                end do
            end do out

        case default      ! if we reach a node that is not subdivided (1 particle or empty cell)
            goal => root  ! this is the desired node

        end select

    end subroutine find_cell



    ! ------ borrar_empty_leaves ------
    ! removes nodes that were created but ended up empty (type = 0)
    ! this is similar to 'borrar_tree', but the latter resets the simulation at each timestep,
    ! while this one optimizes the current tree (removes unnecessary nodes before force computation)
    
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



    ! ------ borrar_tree ------
    ! same as above but for the entire tree

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



    !=======================================================
    ! SYSTEM PHYSICS (MASS AND FORCE CALCULATIONS)
    !=======================================================

    ! ------ calculate_masses ------
    ! after placing all particles geometrically, the node must contain
    ! information about its total enclosed mass and its center of mass

    recursive subroutine calculate_masses(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        real(real64) :: mass, total_mass

        ! initialize mass and COM to zero
        goal%mass = 0.0_real64                                    
        goal%c_o_m = point3d(0.0_real64, 0.0_real64, 0.0_real64)

        select case (goal%type)

        case (1)  ! leaf case (1 particle)
            goal%mass = particles(goal%pos)%m
            goal%c_o_m = particles(goal%pos)%p

        case (2)  ! branch case (multiple particles)
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2

                        if (associated(goal%subcell(i, j, k)%ptr)) then        ! if child exists...
                            call calculate_masses(goal%subcell(i, j, k)%ptr)    ! first we calculate the child's mass

                            mass = goal%mass  ! the accumulated mass of the children found so far (will be the parent’s mass)
                            total_mass = mass + goal%subcell(i, j, k)%ptr%mass  ! new total mass = parent + new child
                            
                            ! calculate weighted center of mass
                            ! [GEOMETRY MODULE]: we use point*scalar + point*scalar / scalar
                            ! formula: (COM_actual * Mass_actual + COM_child * Mass_child) / Mass_total
                            
                            goal%c_o_m = (goal%c_o_m * mass + &
                                              goal%subcell(i, j, k)%ptr%c_o_m * goal%subcell(i, j, k)%ptr%mass) / total_mass 
                            
                            goal%mass = total_mass ! and we update the parent’s mass, which will be the total mass
                        end if

                    end do
                end do
            end do

        end select

    end subroutine calculate_masses




    ! ------ calculate_forces_aux ------
    ! depending on distance, decides whether to apply Barnes-Hut approximation
    ! or treat particles individually to compute the force

    recursive subroutine calculate_forces_aux(goal, tree, force_accum)
        type(cell), pointer :: tree
        integer, intent(in) :: goal
        type(vector3d), intent(inout) :: force_accum
        
        real(real64) :: d, r2, r3, size_box
        type(vector3d) :: rji, force_temp
        integer :: i, j, k
        real(real64), parameter :: epsilon = 0.01_real64  ! softening length

        select case (tree%type)

        case (1)  ! leaf case (1 particle): particle-particle interaction

            if (goal /= tree%pos) then

                ! [GEOMETRY]: we use 'distance2' to calculate r^2 between points
                r2 = distance2(particles(tree%pos)%p, particles(goal)%p)
                    
                rji = particles(tree%pos)%p - particles(goal)%p
                
                r3 = (r2 + epsilon**2) * sqrt(r2 + epsilon**2)   ! small softening to avoid numerical blow-ups

                !we calculte the force F = (m1*m2 * r_vec) / r^3_softened
                force_temp = (particles(tree%pos)%m * particles(goal)%m * rji) / r3
                force_accum = force_accum + force_temp   ! and we sum up every force contribution to the particle
                
            end if

        case (2)  ! branch case (multiple particles)

            r2 = distance2(tree%c_o_m, particles(goal)%p)
            
            ! [GEOMETRY]: we use 'distance' for the opening criterion
            d = sqrt(r2) 
            
            size_box = tree%range%max_corner%x - tree%range%min_corner%x

            ! criterion of opening: size / distance < theta
            if (size_box / d < theta) then
                
                rji = tree%c_o_m - particles(goal)%p
                
                r3 = r2 * d
                
                force_temp = (tree%mass * particles(goal)%m * rji) / r3
                force_accum = force_accum + force_temp

            else
                do i = 1, 2
                    do j = 1, 2
                        do k = 1, 2
                            if (associated(tree%subcell(i, j, k)%ptr)) then
                                call calculate_forces_aux(goal, tree%subcell(i, j, k)%ptr, force_accum)
                            end if
                        end do
                    end do
                end do
            end if

        end select

    end subroutine calculate_forces_aux



    !=======================================================
    ! GEOMETRY AND NUMERICAL UTILITIES
    !=======================================================

    ! ------ get_global_boundaries ------
    ! this subroutine is used by 'barnes_hut_initialize' to define the root cell of the tree

    subroutine get_global_boundaries(goal)
        type(cell), pointer :: goal
        type(point3d) :: mins, maxs, medios
        real(real64) :: span
        integer :: i
        real(real64) :: minx, miny, minz, maxx, maxy, maxz

        ! we want to find the extremes
        ! first, take the position of the first particle as initial reference
        minx = particles(1)%p%x
        maxx = particles(1)%p%x
        miny = particles(1)%p%y
        maxy = particles(1)%p%y
        minz = particles(1)%p%z
        maxz = particles(1)%p%z

        ! update limits if we find any particle farther away
        do i = 2, n_particles
            if (particles(i)%p%x < minx) minx = particles(i)%p%x
            if (particles(i)%p%x > maxx) maxx = particles(i)%p%x
            if (particles(i)%p%y < miny) miny = particles(i)%p%y
            if (particles(i)%p%y > maxy) maxy = particles(i)%p%y
            if (particles(i)%p%z < minz) minz = particles(i)%p%z
            if (particles(i)%p%z > maxz) maxz = particles(i)%p%z
        end do

        ! at the end of this loop, define the minimum bounding box containing all particles
        mins = point3d(minx, miny, minz)
        maxs = point3d(maxx, maxy, maxz)

        ! force cubic shape
        span = max(maxx - minx, maxy - miny, maxz - minz) * 1.1_real64

        ! compute geometric center of the particle cloud
        medios = point3d((maxx + minx) / 2.0_real64, &
                         (maxy + miny) / 2.0_real64, &
                         (maxz + minz) / 2.0_real64)

        ! define root cell box centered at medios with total width span
        goal%range%min_corner = point3d(medios%x - span/2.0_real64, &
                                        medios%y - span/2.0_real64, &
                                        medios%z - span/2.0_real64)
        goal%range%max_corner = point3d(medios%x + span/2.0_real64, &
                                        medios%y + span/2.0_real64, &
                                        medios%z + span/2.0_real64)

    end subroutine get_global_boundaries



    ! ------ get_subcell_limits ------
    ! determines the exact coordinates of a child's corners
    ! based on its parent and which octant it belongs to

    function get_subcell_limits(what, goal, octant)
        integer, intent(in) :: what
        type(cell), pointer, intent(in) :: goal
        integer, dimension(3), intent(in) :: octant   
        type(point3d) :: get_subcell_limits, valor_medio

        ! first, split the parent cube in half along all three axes
        valor_medio = point3d((goal%range%min_corner%x + goal%range%max_corner%x) / 2.0_real64, &
                              (goal%range%min_corner%y + goal%range%max_corner%y) / 2.0_real64, &
                              (goal%range%min_corner%z + goal%range%max_corner%z) / 2.0_real64)

        select case (what)

        case (0)  ! define minimum corner of child cell

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


        case (1)  ! define maximum corner of child cell

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



    ! ------ belongs ------
    ! checks whether a 3D point belongs to a specific cell

    function belongs(part, goal)
        type(point3d), intent(in) :: part
        type(cell), pointer, intent(in) :: goal
        logical :: belongs

        ! to belong, all 6 conditions must be satisfied
        if (part%x >= goal%range%min_corner%x .and. &
            part%x <= goal%range%max_corner%x .and. &
            part%y >= goal%range%min_corner%y .and. &
            part%y <= goal%range%max_corner%y .and. &
            part%z >= goal%range%min_corner%z .and. &
            part%z <= goal%range%max_corner%z) then
            belongs = .true.
        ! Note: we use inclusive comparisons (<=, >=), meaning cell boundaries are included
        ! this is why in 'create_subcells' we added the 'placed' flag, to prevent a particle
        ! that lies exactly on a boundary from being placed twice
        else
            belongs = .false.
        end if

    end function belongs

    
end module barnes_hut_module

