module barnes_hut
    use iso_fortran_env, only : real64
    use geometry
    use particle
    implicit none
    private

    real(real64), parameter :: theta = 1.0_real64
    ! softening parameter to avoid problems if r2=0
    real(real64), parameter :: eps2 = 1.0e-6_real64

    !==================================================
    ! Public types and subroutines
    !==================================================
    public :: cell
    public :: bh_init, bh_build_tree, bh_compute_forces, bh_destroy, borrar_tree

    !==================================================
    ! Data types
    !==================================================
    type range
        type(point3d) :: min,max
    end type range

    type cptr
        type(cell), pointer :: ptr
    end type cptr
    
    type cell
        type (range) :: range
        type(point3d) :: part
        integer :: pos
        integer :: type !! 0 = empty; 1 = particle; 2 = node
        real(real64) :: mass
        type(point3d) :: c_o_m
        type (cptr), dimension(2,2,2) :: subcell
    end type cell

contains

    !====================================================
    ! Subroutine: bh_init (public)
    !====================================================
    subroutine bh_init(head, p)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)

        allocate(head)

        call calculate_ranges(head, p)
        head%type = 0
        call nullify_pointers(head)

    end subroutine bh_init

    !====================================================
    ! Subroutine: bh_build_tree (public)
    !====================================================
    subroutine bh_build_tree(head,p)
        type (cell), pointer :: head, temp_cell
        type(particle3d), intent(in) :: p(:)

        integer :: i,n
        n = size(p)

        do i = 1,n
            call find_cell(head,temp_cell,p(i)%p)
            call place_cell(temp_cell,p(i)%p,i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head,p)
    
    end subroutine bh_build_tree

    !====================================================
    ! Subroutine: bh_compute_forces (public)
    !====================================================
    subroutine bh_compute_forces(head,p,a)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)

        a = vector3d(0.0, 0.0, 0.0)

        call calculate_forces(head,p,a)
    
    end subroutine bh_compute_forces

    !====================================================
    ! Subroutine: bh_destroy (public)
    !====================================================
    ! Suroutine to delete the tree (so we can create a new
    ! one with the updated positions)
    subroutine bh_destroy(head,p)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)

        call borrar_tree(head)
        deallocate(head)
        allocate(head)

        call calculate_ranges(head,p)
        head%type = 0
        call nullify_pointers(head)
    
    end subroutine bh_destroy

    !====================================================
    ! Subroutine: calculate_ranges
    !====================================================
    subroutine calculate_ranges(goal, p)
        type(cell),pointer :: goal
        type(particle3d), intent(in) :: p(:)

        type(point3d) :: mins,maxs,medios
        real(real64) :: span
        integer :: i

        ! Inicialization of the 1st particle
        mins = p(1)%p
        maxs = p(1)%p

        ! Loop over the particles
        do i = 2, size(p)
            mins%x = min(mins%x, p(i)%p%x)
            mins%y = min(mins%y, p(i)%p%y)
            mins%z = min(mins%z, p(i)%p%z)

            maxs%x = max(maxs%x, p(i)%p%x)
            maxs%y = max(maxs%y, p(i)%p%y)
            maxs%z = max(maxs%z, p(i)%p%z)
        end do

        ! Maximum span (the biggest dimension)
        ! When span is calculated, we add a 10% to ensure
        ! that particles do not fall exactly on the edge
        span = max( maxs%x - mins%x, &
                    maxs%y - mins%y, &
                    maxs%z - mins%z ) * 1.1_real64

        medios = 0.5_real64 * (maxs + mins)

        goal%range%min = medios - (span * 0.5_real64)
        goal%range%max = medios + (span * 0.5_real64)
    
    end subroutine calculate_ranges

    !====================================================
    ! Subroutine: find_cell
    !====================================================
    recursive subroutine find_cell(root,goal,part)
        type(cell),pointer :: root,goal,temp
        type(point3d) :: part
        integer :: i,j,k

        select case (root%type)
            case (2)
                out: do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (belongs(part,root%subcell(i,j,k)%ptr)) then
                                call find_cell(root%subcell(i,j,k)%ptr,temp,part)
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

    !====================================================
    ! Subroutine: place_cell
    !====================================================
    recursive subroutine place_cell(goal,part,n)
        type(cell),pointer :: goal,temp
        type(point3d) :: part
        integer :: n

        select case (goal%type)
            case (0)
                goal%type = 1
                goal%part = part
                goal%pos = n
            case (1)
                call crear_subcells(goal)
                call find_cell(goal,temp,part)
                call place_cell(temp,part,n)
            case default
                print*,"SHOULD NOT BE HERE. ERROR!"
        end select
    end subroutine place_cell

    !====================================================
    ! Subroutine: crear_subcells
    !====================================================
    subroutine crear_subcells(goal)
        type(cell), pointer :: goal
        type(point3d) :: part
        integer :: i,j,k
        integer, dimension(3) :: octant
        logical :: placed ! <--- NUEVA VARIABLE BANDERA

        part = goal%part
        goal%type=2
        
        placed = .false. ! Inicializamos a falso

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = calcular_range(0,goal,octant)
                    goal%subcell(i,j,k)%ptr%range%max = calcular_range(1,goal,octant)
                    
                    ! Inicializamos por defecto a vacío (0)
                    goal%subcell(i,j,k)%ptr%type = 0 

                    ! Solo intentamos colocarla si no ha sido colocada ya
                    if (.not. placed) then
                        if (belongs(part,goal%subcell(i,j,k)%ptr)) then
                            goal%subcell(i,j,k)%ptr%part = part
                            goal%subcell(i,j,k)%ptr%type = 1
                            goal%subcell(i,j,k)%ptr%pos = goal%pos
                            placed = .true. ! <--- MARCAMOS COMO COLOCADA
                        end if
                    end if
                    
                    call nullify_pointers(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
        
        ! Verificación de seguridad (opcional):
        if (.not. placed) then
             print *, "Error critico: Particula perdida en subdivision"
             ! Esto podria pasar solo por errores graves de punto flotante
        end if

    end subroutine crear_subcells

    !====================================================
    ! Subroutine: nullify_pointers
    !====================================================
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

    !====================================================
    ! Subroutine: belongs
    !====================================================
    function belongs (part,goal)
        type(point3d) :: part
        type(cell), pointer :: goal
        logical :: belongs

        belongs = (part%x >= goal%range%min%x .and. part%x <= goal%range%max%x) .and. &
          (part%y >= goal%range%min%y .and. part%y <= goal%range%max%y) .and. &
          (part%z >= goal%range%min%z .and. part%z <= goal%range%max%z)


    end function belongs

    !====================================================
    ! Subroutine: calcular_range
    !====================================================
    function calcular_range (what,goal,octant)
        integer :: what
        type(cell), pointer :: goal
        integer, dimension(3) :: octant
        type(point3d) :: calcular_range, valor_medio

        valor_medio = (goal%range%min + goal%range%max) * 0.5_real64
        
        select case (what)
        case (0)
            calcular_range%x = merge(goal%range%min%x, valor_medio%x, octant(1)==1)
            calcular_range%y = merge(goal%range%min%y, valor_medio%y, octant(2)==1)
            calcular_range%z = merge(goal%range%min%z, valor_medio%z, octant(3)==1)
        case (1)
            calcular_range%x = merge(valor_medio%x, goal%range%max%x, octant(1)==1)
            calcular_range%y = merge(valor_medio%y, goal%range%max%y, octant(2)==1)
            calcular_range%z = merge(valor_medio%z, goal%range%max%z, octant(3)==1)
        end select

    end function calcular_range

    !====================================================
    ! Subroutine: borrar_empty_leaves
    !====================================================
    recursive subroutine borrar_empty_leaves(goal)
        type(cell),pointer :: goal
        integer :: i,j,k

        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        call borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate(goal%subcell(i,j,k)%ptr)
                            nullify(goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine borrar_empty_leaves

    !====================================================
    ! Subroutine: borrar_tree (public)
    !====================================================
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

    !====================================================
    ! Subroutine: calculate_masses
    !====================================================
    recursive subroutine calculate_masses(goal,p)
        type(cell),pointer :: goal
        type(particle3d), intent(in) :: p(:)
        integer :: i,j,k
        real(real64) :: mass

        goal%mass = 0
        goal%c_o_m = point3d(0.0, 0.0, 0.0)

        select case (goal%type)
            case (1)
                goal%mass = p(goal%pos)%m
                goal%c_o_m = p(goal%pos)%p
            case (2)
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (associated(goal%subcell(i,j,k)%ptr)) then
                                call calculate_masses(goal%subcell(i,j,k)%ptr,p)
                                mass = goal%mass
                                goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                                goal%c_o_m = (mass * goal%c_o_m + &
                                    goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
                            end if
                        end do
                    end do
                end do
        end select
    end subroutine calculate_masses

    !====================================================
    ! Subroutine: calculate_forces
    !====================================================
    subroutine calculate_forces(head,p,a)
        type(cell),pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(inout) :: a(:)

        integer :: i,n
        n = size(p)
        
        !$omp parallel do schedule(dynamic) private(i) shared(head,p,a)
        do i = 1,n
            call calculate_forces_aux(i,head,p,a)
        end do
        !$omp end parallel do

    end subroutine calculate_forces

    !====================================================
    ! Subroutine: calculate_forces_aux
    !====================================================
    recursive subroutine calculate_forces_aux(goal,tree,p,a)
        integer, intent(in) :: goal
        type(cell), pointer :: tree
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(inout) :: a(:)
        type(vector3d) :: rji
        real(real64) :: r2, r3, l, D
        integer :: i,j,k

        select case (tree%type)
            case (1)
                if (goal .ne. tree%pos) then
                    rji = tree%c_o_m - p(goal)%p
                    r2 = rji%x*rji%x + rji%y*rji%y + rji%z*rji%z + eps2
                    r3 = r2 * sqrt(r2)
                    a(goal) = a(goal) + p(tree%pos)%m * rji / r3
                end if
            case (2)
                !! El rango tiene el mismo span en las 3 dimensiones
                !! por lo que podemos considerar una dimension cualquiera
                !! para calcular el lado de la celda (en este caso la
                !! dimension x)
                l = tree%range%max%x - tree%range%min%x
                rji = tree%c_o_m - p(goal)%p
                r2 = rji%x*rji%x + rji%y*rji%y + rji%z*rji%z + eps2
                D  = sqrt(r2)
                if (l/D < theta) then
                    !! Si conglomerado, tenemos que ver si se cumple l/D < @
                    r3 = r2 * D
                    a(goal) = a(goal) + tree%mass * rji / r3
                else
                    do i = 1,2
                        do j = 1,2
                            do k = 1,2
                                if (associated(tree%subcell(i,j,k)%ptr)) then
                                    call calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr,p,a)
                                end if
                            end do
                        end do
                    end do
                end if
        end select
    end subroutine calculate_forces_aux

end module barnes_hut