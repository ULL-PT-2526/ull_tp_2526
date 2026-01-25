module barnes_hut
    use geometry
    use particle
    !$ use omp_lib
    implicit none

    ! Barnes–Hut opening angle
    real(dp), parameter :: theta = 1.0_dp

    ! Number of particles
    integer :: n

    ! Particle array 
    type(particle3d), allocatable :: parts(:)

    ! Accelerations
    type(vector3d), allocatable :: a(:)

    !=========================================================
    ! Spatial range
    !=========================================================
    type :: range
        type(point3d) :: min, max
    end type range

    type :: cptr
        type(cell), pointer :: ptr
    end type cptr

    !=========================================================
    ! Barnes–Hut cell 
    !=========================================================
    type :: cell
        type(range) :: range
        type(point3d) :: part
        integer :: pos
        integer :: type        ! 0 = empty, 1 = particle, 2 = conglomerate
        real(dp) :: mass
        type(point3d) :: c_o_m
        type(cptr), dimension(2,2,2) :: subcell
    end type cell

    type(cell), pointer :: head => null(), temp_cell => null()

contains

!===========================================================
! Calculate_Ranges
!===========================================================
    subroutine calculate_ranges(goal)
        type(cell), pointer :: goal
        real(dp) :: span
        type(point3d) :: mins, maxs, medios
        integer :: i

        mins = parts(1)%p
        maxs = parts(1)%p

        do i = 2, n
            mins%x = min(mins%x, parts(i)%p%x)
            mins%y = min(mins%y, parts(i)%p%y)
            mins%z = min(mins%z, parts(i)%p%z)

            maxs%x = max(maxs%x, parts(i)%p%x)
            maxs%y = max(maxs%y, parts(i)%p%y)
            maxs%z = max(maxs%z, parts(i)%p%z)
        end do

        span = max( maxs%x-mins%x, max(maxs%y-mins%y, maxs%z-mins%z) ) * 1.1_dp
        medios = point3d( (maxs%x+mins%x)/2, (maxs%y+mins%y)/2, (maxs%z+mins%z)/2 )

        goal%range%min = point3d( medios%x-span/2, medios%y-span/2, medios%z-span/2 )
        goal%range%max = point3d( medios%x+span/2, medios%y+span/2, medios%z+span/2 )
    end subroutine calculate_ranges

!===========================================================
! Find_Cell
!===========================================================
    recursive subroutine find_cell(root, goal, part)
        type(cell), pointer :: root, goal, temp
        type(point3d), intent(in) :: part
        integer :: i, j, k

        select case (root%type)
        case (2)
            out: do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (belongs(part, root%subcell(i,j,k)%ptr)) then
                            call find_cell(root%subcell(i,j,k)%ptr, temp, part)
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

!===========================================================
! Place_Cell
!===========================================================
    recursive subroutine place_cell(goal, part, idx)
        type(cell), pointer :: goal, temp
        type(point3d), intent(in) :: part
        integer, intent(in) :: idx

        select case (goal%type)
        case (0)
            goal%type = 1
            goal%part = part
            goal%pos  = idx
        case (1)
            call crear_subcells(goal)
            call find_cell(goal, temp, part)
            call place_cell(temp, part, idx)
        case default
            print *, "ERROR in place_cell"
        end select
    end subroutine place_cell

!===========================================================
! Create_Subcells
!===========================================================
    subroutine crear_subcells(goal)
        type(cell), pointer :: goal
        type(point3d) :: part
        integer :: i, j, k
        integer, dimension(3) :: octant

        part = goal%part
        goal%type = 2

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    octant = (/i,j,k/)
                    allocate(goal%subcell(i,j,k)%ptr)
                    goal%subcell(i,j,k)%ptr%range%min = calcular_range(0, goal, octant)
                    goal%subcell(i,j,k)%ptr%range%max = calcular_range(1, goal, octant)

                    if (belongs(part, goal%subcell(i,j,k)%ptr)) then
                        goal%subcell(i,j,k)%ptr%part = part
                        goal%subcell(i,j,k)%ptr%pos  = goal%pos
                        goal%subcell(i,j,k)%ptr%type = 1
                    else
                        goal%subcell(i,j,k)%ptr%type = 0
                    end if

                    call nullify_pointers(goal%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine crear_subcells

!===========================================================
! Nullify_Pointers
!===========================================================
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

!===========================================================
! Belongs
!===========================================================
    function belongs(part, goal)
        type(point3d) :: part
        type(cell), pointer :: goal
        logical :: belongs
        
        if (part%x >= goal%range%min%x .and. &
            part%x <= goal%range%max%x .and. &
            part%y >= goal%range%min%y .and. &
            part%y <= goal%range%max%y .and. &
            part%z >= goal%range%min%z .and. &
            part%z <= goal%range%max%z) then
            belongs = .true.
        else
            belongs = .false.
        end if
        
    end function belongs

!===========================================================
! Calculate_Range
!===========================================================
    function calcular_range(what, goal, octant) result(r)
        integer, intent(in) :: what
        type(cell), pointer :: goal
        integer, dimension(3), intent(in) :: octant
        type(point3d) :: r, mid

        mid = point3d( (goal%range%min%x+goal%range%max%x)/2, &
                       (goal%range%min%y+goal%range%max%y)/2, &
                       (goal%range%min%z+goal%range%max%z)/2 )

        select case (what)
        case (0)
            r%x = merge(goal%range%min%x, mid%x, octant(1)==1)
            r%y = merge(goal%range%min%y, mid%y, octant(2)==1)
            r%z = merge(goal%range%min%z, mid%z, octant(3)==1)
        case (1)
            r%x = merge(mid%x, goal%range%max%x, octant(1)==1)
            r%y = merge(mid%y, goal%range%max%y, octant(2)==1)
            r%z = merge(mid%z, goal%range%max%z, octant(3)==1)
        end select
    end function calcular_range

!===========================================================
! Delete_empty_leaves
!===========================================================
    recursive subroutine borrar_empty_leaves(goal)
        type(cell), pointer :: goal
        integer :: i,j,k

        if (associated(goal%subcell(1,1,1)%ptr)) then
            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        call borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        if (goal%subcell(i,j,k)%ptr%type == 0) then
                            deallocate(goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end if
    end subroutine borrar_empty_leaves

!===========================================================
! Delete_tree
!===========================================================
    recursive subroutine borrar_tree(goal)
        type(cell), pointer :: goal
        integer :: i,j,k

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    if (associated(goal%subcell(i,j,k)%ptr)) then
                        call borrar_tree(goal%subcell(i,j,k)%ptr)
                        deallocate(goal%subcell(i,j,k)%ptr)
                    end if
                end do
            end do
        end do
    end subroutine borrar_tree

!===========================================================
! Calculate_masses
!===========================================================
    recursive subroutine calculate_masses(goal)
        type(cell), pointer :: goal
        integer :: i, j, k
        real(dp) :: mass
        !type(point3d) :: c_o_m
        
        goal%mass = 0.0_dp
        goal%c_o_m = point3d(0.0_dp, 0.0_dp, 0.0_dp)
        
        select case (goal%type)
        case (1)
            goal%mass = parts(goal%pos)%m
            goal%c_o_m = parts(goal%pos)%p
        case (2)
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        if (associated(goal%subcell(i,j,k)%ptr)) then
                            call calculate_masses(goal%subcell(i,j,k)%ptr)
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

!===========================================================
! Calculate_forces
!===========================================================
    subroutine calculate_forces(head)
        type(cell), pointer :: head
        integer :: i

        a = vector3d(0.0_dp,0.0_dp,0.0_dp)

        !$omp parallel default(shared) private(i)
        !$omp do
        do i = 1,n
            call calculate_forces_aux(i,head)
        end do
        !$omp end do
        !$omp end parallel 
    end subroutine calculate_forces

!===========================================================
! Calculate_forces_mpi
!===========================================================
    subroutine calculate_forces_mpi(head, istart, iend)
        type(cell), pointer :: head
        integer, intent(in) :: istart, iend
        integer :: i

        do i = istart, iend
            call calculate_forces_aux(i, head)
        end do
    end subroutine calculate_forces_mpi

!===========================================================
! Calculate_forces_aux
!===========================================================
    recursive subroutine calculate_forces_aux(idx, tree)
        integer, intent(in) :: idx
        type(cell), pointer :: tree
        integer :: i,j,k
        real(dp) :: l, D, r2, r3
        type(vector3d) :: rji

        select case (tree%type)
        case (1)
            if (idx /= tree%pos) then
                rji = tree%c_o_m - parts(idx)%p
                r2 = dot_prod(rji, rji)
                r3 = r2*sqrt(r2)
                a(idx) = a(idx) + (parts(tree%pos)%m / r3)*rji
            end if

        case (2)
            l = tree%range%max%x - tree%range%min%x
            rji = tree%c_o_m - parts(idx)%p
            r2 = dot_prod(rji, rji)
            D  = sqrt(r2)

            if (l/D < theta) then
                r3 = r2 * D
                a(idx) = a(idx) + (tree%mass / r3) * rji
            else
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (associated(tree%subcell(i,j,k)%ptr)) then
                                call calculate_forces_aux(idx, tree%subcell(i,j,k)%ptr)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine calculate_forces_aux

end module barnes_hut
