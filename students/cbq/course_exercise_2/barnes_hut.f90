!!!! Barnes-Hut Tree Algorithm Module

module barnes_hut
    use iso_fortran_env, only: real64
    use geometry
    use particle
!$ use omp_lib
    implicit none
    
    ! Opening angle parameter for Barnes-Hut approximation
    real(real64), parameter :: theta = 1.0_real64
    
    ! Range type: stores min and max coordinates of a cell
    type range
        type(point3d) :: min, max
    end type range
    
    ! Cell pointer type (needed for recursive tree structure)
    type cptr
        type(cell), pointer :: ptr
    end type cptr
    
    ! Cell type: represents a node in the octree
    type cell
        type(range) :: range              ! Spatial boundaries
        type(point3d) :: part             ! Particle position (if leaf)
        integer :: pos                    ! Particle index (if leaf)
        integer :: type                   ! 0=empty, 1=leaf, 2=internal node
        real(real64) :: mass              ! Total mass in subtree
        type(point3d) :: c_o_m            ! Center of mass
        type(cptr), dimension(2,2,2) :: subcell  ! 8 children (octants)
    end type cell
    
contains


    ! Calculate spatial range for root cell

    subroutine calculate_ranges(goal, particles, n)
        type(cell), pointer :: goal
        type(particle3d), dimension(:), intent(in) :: particles
        integer, intent(in) :: n
        type(point3d) :: mins, maxs, medios
        real(real64) :: span
        integer :: i
        
        mins = particles(1)%p
        maxs = particles(1)%p
        
        ! Find bounding box
        do i = 2, n
            mins%x = min(mins%x, particles(i)%p%x)
            mins%y = min(mins%y, particles(i)%p%y)
            mins%z = min(mins%z, particles(i)%p%z)
            maxs%x = max(maxs%x, particles(i)%p%x)
            maxs%y = max(maxs%y, particles(i)%p%y)
            maxs%z = max(maxs%z, particles(i)%p%z)
        end do
        
        ! Create cubic domain with 10% margin
        span = max(maxs%x - mins%x, maxs%y - mins%y, maxs%z - mins%z) * 1.1_real64
        medios = point3d((maxs%x + mins%x) / 2.0_real64, &
                        (maxs%y + mins%y) / 2.0_real64, &
                        (maxs%z + mins%z) / 2.0_real64)
        
        goal%range%min = point3d(medios%x - span/2.0_real64, &
                                 medios%y - span/2.0_real64, &
                                 medios%z - span/2.0_real64)
        goal%range%max = point3d(medios%x + span/2.0_real64, &
                                 medios%y + span/2.0_real64, &
                                 medios%z + span/2.0_real64)
    end subroutine calculate_ranges


    ! Find appropriate cell for a particle

    recursive subroutine find_cell(root, goal, part)
        type(point3d), intent(in) :: part
        type(cell), pointer :: root, goal, temp
        integer :: i, j, k
        
        select case (root%type)
        case (2)  ! Internal node: search children
            out: do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (belongs(part, root%subcell(i,j,k)%ptr)) then
                            call find_cell(root%subcell(i,j,k)%ptr, temp, part)
                            goal => temp
                            exit out
                        end if
                    end do
                end do
            end do out
        case default  ! Empty or leaf: found target cell
            goal => root
        end select
    end subroutine find_cell


    ! Place particle in tree

    recursive subroutine place_cell(goal, part, n)
        type(cell), pointer :: goal, temp
        type(point3d), intent(in) :: part
        integer, intent(in) :: n
        
        select case (goal%type)
        case (0)  ! Empty cell: place particle here
            goal%type = 1
            goal%part = part
            goal%pos = n
        case (1)  ! Leaf cell: subdivide and place both particles
            call crear_subcells(goal)
            call find_cell(goal, temp, part)
            call place_cell(temp, part, n)
        end select
    end subroutine place_cell


    ! Create 8 subcells (octants)

    subroutine crear_subcells(goal)
        type(cell), pointer :: goal
        type(point3d) :: part
        integer :: i, j, k
        integer, dimension(3) :: octant
        
        part = goal%part
        goal%type = 2  ! Mark as internal node
        
        ! Create all 8 children
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    octant = (/i, j, k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = calcular_range(0, goal, octant)
                    goal%subcell(i,j,k)%ptr%range%max = calcular_range(1, goal, octant)
                    
                    ! Place existing particle in appropriate child
                    if (belongs(part, goal%subcell(i,j,k)%ptr)) then
                        goal%subcell(i,j,k)%ptr%part = part
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


    ! Initialize cell pointers to null

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


    ! Check if particle belongs to cell

    function belongs(part, goal)
        type(point3d), intent(in) :: part
        type(cell), pointer :: goal
        logical :: belongs
        
        if (part%x >= goal%range%min%x .and. part%x <= goal%range%max%x .and. &
            part%y >= goal%range%min%y .and. part%y <= goal%range%max%y .and. &
            part%z >= goal%range%min%z .and. part%z <= goal%range%max%z) then
            belongs = .true.
        else
            belongs = .false.
        end if
    end function belongs


    ! Calculate octant boundaries

    function calcular_range(what, goal, octant)
        integer, intent(in) :: what
        type(cell), pointer :: goal
        integer, dimension(3), intent(in) :: octant
        type(point3d) :: calcular_range, valor_medio
        
        valor_medio = point3d((goal%range%min%x + goal%range%max%x) / 2.0_real64, &
                             (goal%range%min%y + goal%range%max%y) / 2.0_real64, &
                             (goal%range%min%z + goal%range%max%z) / 2.0_real64)
        
        select case (what)
        case (0)  ! Minimum bounds
            calcular_range%x = merge(goal%range%min%x, valor_medio%x, octant(1) == 1)
            calcular_range%y = merge(goal%range%min%y, valor_medio%y, octant(2) == 1)
            calcular_range%z = merge(goal%range%min%z, valor_medio%z, octant(3) == 1)
        case (1)  ! Maximum bounds
            calcular_range%x = merge(valor_medio%x, goal%range%max%x, octant(1) == 1)
            calcular_range%y = merge(valor_medio%y, goal%range%max%y, octant(2) == 1)
            calcular_range%z = merge(valor_medio%z, goal%range%max%z, octant(3) == 1)
        end select
    end function calcular_range


    ! Remove empty leaf cells

    recursive subroutine borrar_empty_leaves(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        
        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate(goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine borrar_empty_leaves


    ! Delete entire tree (except root)

    recursive subroutine borrar_tree(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate(goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine borrar_tree


    ! Calculate masses and center of mass

    recursive subroutine calculate_masses(goal, particles)
        type(cell), pointer :: goal
        type(particle3d), dimension(:), intent(in) :: particles
        integer :: i, j, k
        real(real64) :: mass
        
        goal%mass = 0.0_real64
        goal%c_o_m = point3d(0.0_real64, 0.0_real64, 0.0_real64)
        
        select case (goal%type)
        case (1)  ! Leaf: copy particle data
            goal%mass = particles(goal%pos)%m
            goal%c_o_m = particles(goal%pos)%p
        case (2)  ! Internal node: sum children
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (associated(goal%subcell(i,j,k)%ptr)) then
                            call calculate_masses(goal%subcell(i,j,k)%ptr, particles)
                            mass = goal%mass
                            goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                            goal%c_o_m%x = (mass * goal%c_o_m%x + &
                                goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%mass
                            goal%c_o_m%y = (mass * goal%c_o_m%y + &
                                goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%mass
                            goal%c_o_m%z = (mass * goal%c_o_m%z + &
                                goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%mass
                        end if
                    end do
                end do
            end do
        end select
    end subroutine calculate_masses


    ! Calculate gravitational forces
    ! Parallelized with OpenMP if compiled with -fopenmp

    subroutine calculate_forces(head, particles, accelerations, n)
        type(cell), pointer :: head
        type(particle3d), dimension(:), intent(in) :: particles
        type(vector3d), dimension(:), intent(inout) :: accelerations
        integer, intent(in) :: n
        integer :: i
        
        ! OpenMP parallel loop: each thread calculates forces for different particles
        !$omp parallel do default(shared) private(i) schedule(dynamic)
        do i = 1, n
            call calculate_forces_aux(i, head, particles, accelerations)
        end do
        !$omp end parallel do
    end subroutine calculate_forces


    ! Auxiliary recursive force calculation

    recursive subroutine calculate_forces_aux(goal, tree, particles, accelerations)
        type(cell), pointer :: tree
        type(particle3d), dimension(:), intent(in) :: particles
        type(vector3d), dimension(:), intent(inout) :: accelerations
        integer, intent(in) :: goal
        integer :: i, j, k
        real(real64) :: l, D, r2, r3
        type(vector3d) :: rji
        
        select case (tree%type)
        case (1)  ! Leaf: direct particle-particle interaction
            if (goal /= tree%pos) then
                rji = tree%c_o_m - particles(goal)%p
                r2 = rji%x**2 + rji%y**2 + rji%z**2
                r3 = r2 * sqrt(r2)
                accelerations(goal) = accelerations(goal) + &
                    (particles(tree%pos)%m / r3) * rji
            end if
        case (2)  ! Internal node: check Barnes-Hut criterion
            l = tree%range%max%x - tree%range%min%x
            rji = tree%c_o_m - particles(goal)%p
            r2 = rji%x**2 + rji%y**2 + rji%z**2
            D = sqrt(r2)
            
            if (l/D < theta) then
                ! Cell far enough: treat as single particle
                r3 = r2 * D
                accelerations(goal) = accelerations(goal) + &
                    (tree%mass / r3) * rji
            else
                ! Cell too close: recurse into children
                do i = 1, 2
                    do j = 1, 2
                        do k = 1, 2
                            if (associated(tree%subcell(i,j,k)%ptr)) then
                                call calculate_forces_aux(goal, tree%subcell(i,j,k)%ptr, &
                                    particles, accelerations)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine calculate_forces_aux

end module barnes_hut