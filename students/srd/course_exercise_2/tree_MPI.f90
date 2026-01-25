PROGRAM tree
use mpi
use geometry
use particle
IMPLICIT NONE

! MPI control variables
integer :: rank                   ! Process ID
integer :: nprocs                 ! Total number of processes
integer :: ierr                   ! MPI error flag
integer :: status(MPI_STATUS_SIZE) ! MPI status array
integer, allocatable :: counts(:), offsets(:)
character(len=20) :: mode_name = "MPI (Dist.)"
! Important for the distributed version:
integer :: n_total                ! Global number of particles 
integer :: n_local                ! Number of particles in a given process
integer :: i_start, i_end         ! Global indices range for a given process
real(kind=bit64) :: tmp_pos(3), tmp_vel(3), tmp_mass ! Temporary buffers for distribution
real(kind=bit64), allocatable :: tmp_pos_list(:), full_pos_list(:)
real(kind=bit64), allocatable :: tmp_mass_list(:), full_mass_list(:)
type(particle3d), allocatable :: p_global(:)
real(kind=bit64), allocatable :: global_positions(:)


INTEGER :: i,j,k
INTEGER, PARAMETER :: out_unit = 10 ! Logic number for the output file
real(kind=bit64) :: dt, t_end, t, dt_out, t_out, rs, r2, r3 ! Maintaining the precision type between modules
real(kind=bit64), parameter :: theta = 1.0_bit64
real(kind=bit64), parameter :: epsilon = 0.01_bit64 ! softening parameter 
type(particle3d), allocatable :: p(:) ! Using the mass, position and velocities from the particle module
type(vector3d) :: rji 
type(vector3d), allocatable :: a(:)
! Time control variables:
integer(kind=8)  :: c_ini, c_fin, c_rate
real(kind=bit64) :: t_tree, t_forces, t_update, t_total

TYPE RANGE
    REAL(kind=bit64), DIMENSION(3) :: min,max
END TYPE RANGE

TYPE CPtr
    TYPE(CELL), POINTER :: ptr
END TYPE CPtr

TYPE CELL
    TYPE(RANGE) :: range
    TYPE(particle3d) :: part
    INTEGER :: pos
    INTEGER :: type        !! 0 = no particle; 1 = particle; 2 = conglomerado
    REAL(kind=bit64) :: mass
    TYPE(vector3d) :: c_o_m          ! Changed to be a vector
    TYPE(CPtr), DIMENSION(2,2,2) :: subcell
END TYPE CELL

TYPE(CELL), POINTER :: head, temp_cell

! MPI initialization
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
! Starting the global timer
if (rank == 0) call system_clock(c_ini, c_rate)
! Only the first rank reads the input file
if (rank == 0) then
    ! Same order as the input.dat
    read(*, *) dt
    read(*, *) dt_out
    read(*, *) t_end
    read(*, *) n_total
end if

! Variables that we want to share across all processes
call MPI_BCAST(n_total, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! Split the n_total particles among the nprocs
n_local = n_total / nprocs
i_start = rank * n_local + 1
i_end   = (rank + 1) * n_local
    
! Fix the last process in case of remainder
if (rank == nprocs - 1) i_end = n_total
! Recalculate n_local for the last process
n_local = i_end - i_start + 1

! Allocate the local memory
allocate(p(n_local))
allocate(a(n_local))

! Since we use "< input.dat", only Rank 0 can read the data

if (rank == 0) then
        ! Root reads its own local particles
    do i = 1, n_local
        read(*, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
                     p(i)%v%x, p(i)%v%y, p(i)%v%z
    end do
        
    ! Root reads and sends data to each worker process
     do j = 1, nprocs - 1
        ! Calculate how many particles process 'j' expects
        ! (Assuming same distribution logic as in Step 5)
        k = n_total / nprocs
        if (j == nprocs - 1) k = n_total - (j * k)
            
        do i = 1, k
	    ! Read from stdin into temporary buffers to avoid overwriting simulation parameters
            read(*, *) tmp_mass, tmp_pos(1), tmp_pos(2), tmp_pos(3), &
                         tmp_vel(1), tmp_vel(2), tmp_vel(3)
	    call MPI_SEND(tmp_mass, 1, MPI_DOUBLE_PRECISION, j, 10, MPI_COMM_WORLD, ierr) 
            call MPI_SEND(tmp_pos,  3, MPI_DOUBLE_PRECISION, j, 11, MPI_COMM_WORLD, ierr) 
            call MPI_SEND(tmp_vel,  3, MPI_DOUBLE_PRECISION, j, 12, MPI_COMM_WORLD, ierr) 
        end do
    end do
else
    ! The other processes receive their assigned particles from Rank 0
    do i = 1, n_local
        call MPI_RECV(p(i)%m, 1, MPI_DOUBLE_PRECISION, 0, 10, MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(tmp_pos, 3, MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, status, ierr)
        p(i)%p%x = tmp_pos(1)
        p(i)%p%y = tmp_pos(2)
        p(i)%p%z = tmp_pos(3)
        call MPI_RECV(tmp_vel, 3, MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, status, ierr)
        p(i)%v%x = tmp_vel(1)
        p(i)%v%y = tmp_vel(2)
        p(i)%v%z = tmp_vel(3)
    end do
end if

! Wait for the distribution to finish
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! We need to store all positions to build the tree locally
allocate(global_positions(n_total))
allocate(counts(nprocs), offsets(nprocs))

! Share the n_local of each process so everyone knows the distribution
call MPI_ALLGATHER(n_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

! Calculate displacements
offsets(1) = 0
do i = 2, nprocs
    offsets(i) = offsets(i-1) + counts(i-1)
end do

! We need all masses to build the tree, but they don't change over time -> one list
allocate(tmp_mass_list(n_local))
allocate(full_mass_list(n_total))

do i = 1, n_local
    tmp_mass_list(i) = p(i)%m
end do

! We use 'counts' and 'offsets' because ranks might have different number of particles
call MPI_ALLGATHERV(tmp_mass_list, n_local, MPI_DOUBLE_PRECISION, &
                    full_mass_list, counts, offsets, MPI_DOUBLE_PRECISION, &
                    MPI_COMM_WORLD, ierr)

! Store them in our p_global array (which we use for the tree)
allocate(p_global(n_total))
do i = 1, n_total
    p_global(i)%m = full_mass_list(i)
end do

! Clean up temporary mass buffers
deallocate(tmp_mass_list, full_mass_list)

! Final sync before the big loop
call MPI_BARRIER(MPI_COMM_WORLD, ierr)




if (rank == 0) then
    OPEN(unit = out_unit, file = "output.dat", status = "replace", action = "write") ! Only root opens the file
end if

! Initialization of the time variables
t_tree   = 0.0_bit64
t_forces = 0.0_bit64
t_update = 0.0_bit64
! Obtaining the tick frequency
call system_clock(count_rate=c_rate)

!! Bucle principal
t = 0.0_bit64
t_out = 0.0_bit64

nullify(head) 
allocate(head)
call Nullify_Pointers(head)

! First, we need initial accelerations
allocate(tmp_pos_list(3*n_local))
allocate(full_pos_list(3*n_total))

! Initial force calculation before loop
call update_global_tree_and_times(t_tree)
call calculate_local_forces_and_times(t_forces)

DO WHILE (t <= t_end)

    call system_clock(count=c_ini) ! Clock initialization
    DO i = 1,n_local
        p(i)%v = p(i)%v + a(i) * (dt/2.0_bit64)
        p(i)%p = p(i)%p + p(i)%v * dt
    END DO
    call system_clock(count=c_fin) ! Clock stop
    t_update = t_update + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

    ! Update global tree
    call update_global_tree_and_times(t_tree)

    ! Calculate forces
    call calculate_local_forces_and_times(t_forces) ! Packaging the subrutines

    call system_clock(count=c_ini)
    DO i = 1, n_local
        p(i)%v = p(i)%v + a(i) * (dt/2.0_bit64)
    END DO
    
    t_out = t_out + dt
    IF (t_out >= dt_out) THEN
        do i = 1, n_local
            tmp_pos_list(3*i-2:3*i) = (/ p(i)%p%x, p(i)%p%y, p(i)%p%z /)
        end do
        ! Gather the rank 0
        call MPI_GATHERV(tmp_pos_list, 3*n_local, MPI_DOUBLE_PRECISION, &
                         full_pos_list, 3*counts, 3*offsets, MPI_DOUBLE_PRECISION, &
                         0, MPI_COMM_WORLD, ierr)
        
        IF (rank == 0) THEN
            WRITE(out_unit, *) t, ( full_pos_list(i), i = 1, 3*n_total )
        END IF
        t_out = 0.0_bit64
    END IF
    
    call system_clock(count=c_fin)
    t_update = t_update + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

    t = t + dt
END DO

CLOSE(out_unit)

t_total = t_tree + t_forces + t_update

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if (rank == 0) then
    PRINT*, ""
    PRINT*, "==============================================="
    PRINT*, "      PERFORMANCE REPORT (MPI)        "
    PRINT*, "==============================================="
    PRINT '(A, I4)',       " Number of Processes:      ", nprocs
    PRINT '(A, F12.4, A)', " Total Calculation Time:   ", t_total, " s"
    PRINT*, "-----------------------------------------------"
    PRINT '(A, F12.4, A, F6.2, A)', " 1. Tree Management:      ", t_tree,   " s | ", (t_tree/t_total)*100.0, "%"
    PRINT '(A, F12.4, A, F6.2, A)', " 2. Force Calculation:    ", t_forces, " s | ", (t_forces/t_total)*100.0, "%"
    PRINT '(A, F12.4, A, F6.2, A)', " 3. Integration & I/O:    ", t_update, " s | ", (t_update/t_total)*100.0, "%"
    PRINT*, "==============================================="
end if

if (rank == 0) then
    CLOSE(out_unit)
    PRINT*, "==============================================="
    PRINT*, "Simulation complete. Data saved."
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_FINALIZE(ierr)

CALL EXIT(0)

CONTAINS

SUBROUTINE Calculate_Ranges(goal, p, n)
    TYPE(CELL), POINTER :: goal
    type(particle3d), intent(in) :: p(:)
    integer, intent(in) :: n

    real(kind=bit64), dimension(3) :: mins, maxs, medios
    real(kind=bit64) :: span
    integer :: i

    ! Inicializamos con la primera partícula
    mins = [ p(1)%p%x, p(1)%p%y, p(1)%p%z ]
    maxs = mins

    ! Recorremos el resto de partículas
    DO i = 2, n
        mins = MIN( mins, [ p(i)%p%x, p(i)%p%y, p(i)%p%z ] )
        maxs = MAX( maxs, [ p(i)%p%x, p(i)%p%y, p(i)%p%z ] )
    END DO

    span   = MAXVAL(maxs - mins) * 1.1_bit64
    medios = (maxs + mins) / 2.0_bit64

    goal%range%min = medios - span/2.0_bit64
    goal%range%max = medios + span/2.0_bit64
END SUBROUTINE Calculate_Ranges


RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
    TYPE(particle3d), INTENT(IN) :: part
    TYPE(CELL),POINTER :: root,goal,temp
    INTEGER :: i,j,k

    SELECT CASE (root%type)
    CASE (2)
out:    DO i = 1,2
            DO j = 1,2
                DO k = 1,2
                    IF (Belongs(part,root%subcell(i,j,k)%ptr)) THEN
                        CALL Find_Cell(root%subcell(i,j,k)%ptr,temp,part)
                        goal => temp
                        EXIT out
                    END IF
                END DO
            END DO
        END DO out
    CASE DEFAULT
        goal => root
    END SELECT
END SUBROUTINE Find_Cell

RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
    TYPE(CELL),POINTER :: goal,temp
    TYPE(particle3d), INTENT(IN) :: part
    INTEGER :: n

    SELECT CASE (goal%type)
    CASE (0)
        goal%type = 1
        goal%part = part
        goal%pos = n
    CASE (1)
	IF (ABS(goal%part%p%x - part%p%x) < 1.e-8_bit64 .AND. &
            ABS(goal%part%p%y - part%p%y) < 1.e-8_bit64) THEN
            RETURN
	END IF

        CALL Crear_Subcells(goal)
        CALL Find_Cell(goal,temp,part)
        CALL Place_Cell(temp,part,n)
    CASE DEFAULT
        PRINT*,"SHOULD NOT BE HERE. ERROR!"
    END SELECT
END SUBROUTINE Place_Cell

SUBROUTINE Crear_Subcells(goal)
    TYPE(CELL), POINTER :: goal
    TYPE(particle3d) :: part
    INTEGER :: i,j,k
    INTEGER, DIMENSION(3) :: octant

    part = goal%part
    goal%type = 2

    DO i = 1,2
        DO j = 1,2
            DO k = 1,2
                octant = (/i,j,k/)
                ALLOCATE(goal%subcell(i,j,k)%ptr)
                goal%subcell(i,j,k)%ptr%range%min = Calcular_Range(0,goal,octant)
                goal%subcell(i,j,k)%ptr%range%max = Calcular_Range(1,goal,octant)

                IF (Belongs(part,goal%subcell(i,j,k)%ptr)) THEN
                    goal%subcell(i,j,k)%ptr%part = part
                    goal%subcell(i,j,k)%ptr%type = 1
                    goal%subcell(i,j,k)%ptr%pos = goal%pos
                ELSE
                    goal%subcell(i,j,k)%ptr%type = 0
                END IF

                CALL Nullify_Pointers(goal%subcell(i,j,k)%ptr)
            END DO
        END DO
    END DO
END SUBROUTINE Crear_Subcells

SUBROUTINE Nullify_Pointers(goal)
    TYPE(CELL), POINTER :: goal
    INTEGER :: i,j,k

    DO i = 1,2
        DO j = 1,2
            DO k = 1,2
                NULLIFY(goal%subcell(i,j,k)%ptr)
            END DO
        END DO
    END DO
END SUBROUTINE Nullify_Pointers

FUNCTION Belongs (part,goal)
    TYPE(particle3d), INTENT(IN) :: part
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs

    IF (part%p%x >= goal%range%min(1) .AND. part%p%x <= goal%range%max(1) .AND. &
        part%p%y >= goal%range%min(2) .AND. part%p%y <= goal%range%max(2) .AND. &
        part%p%z >= goal%range%min(3) .AND. part%p%z <= goal%range%max(3)) THEN
        Belongs = .TRUE.
    ELSE
        Belongs = .FALSE.
    END IF
END FUNCTION Belongs

FUNCTION Calcular_Range (what,goal,octant)
    INTEGER :: what
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    REAL(kind=bit64), DIMENSION(3) :: Calcular_Range, valor_medio

    valor_medio = (goal%range%min + goal%range%max) / 2.0_bit64
    SELECT CASE (what)
    CASE (0)
        WHERE (octant == 1)
            Calcular_Range = goal%range%min
        ELSEWHERE
            Calcular_Range = valor_medio
        END WHERE
    CASE (1)
        WHERE (octant == 1)
            Calcular_Range = valor_medio
        ELSEWHERE
            Calcular_Range = goal%range%max
        END WHERE
    END SELECT
END FUNCTION Calcular_Range

RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k

    IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
        DO i = 1,2
            DO j = 1,2
                DO k = 1,2
                    CALL Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                    IF (goal%subcell(i,j,k)%ptr%type == 0) THEN
                        DEALLOCATE(goal%subcell(i,j,k)%ptr)
                    END IF
                END DO
            END DO
        END DO
    END IF
END SUBROUTINE Borrar_empty_leaves

RECURSIVE SUBROUTINE Borrar_tree(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k

    DO i = 1,2
        DO j = 1,2
            DO k = 1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                    CALL Borrar_tree(goal%subcell(i,j,k)%ptr)
                    DEALLOCATE(goal%subcell(i,j,k)%ptr)
                END IF
            END DO
        END DO
    END DO
END SUBROUTINE Borrar_tree

RECURSIVE SUBROUTINE Calculate_masses(goal, p)
    TYPE(CELL), POINTER :: goal
    type(particle3d), intent(in) :: p(:)
    INTEGER :: i, j, k

    ! Inicialización
    goal%mass   = 0.0_bit64
    goal%c_o_m  = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64)

    SELECT CASE (goal%type)
    CASE (1)   ! Nodo hoja con partícula
        goal%mass   = p(goal%pos)%m
        goal%c_o_m = vector3d(p(goal%pos)%p%x, p(goal%pos)%p%y, p(goal%pos)%p%z)

    CASE (2)   ! Nodo interno
        DO i = 1,2
            DO j = 1,2
                DO k = 1,2
                    IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                        CALL Calculate_masses(goal%subcell(i,j,k)%ptr, p)

                        ! Centro de masas ponderado directamente
                        IF (goal%mass + goal%subcell(i,j,k)%ptr%mass /= 0.0_bit64) THEN
                            goal%c_o_m = (goal%c_o_m * goal%mass + &
                                          goal%subcell(i,j,k)%ptr%c_o_m * &
                                          goal%subcell(i,j,k)%ptr%mass) / &
                                          (goal%mass + goal%subcell(i,j,k)%ptr%mass)
                        END IF

                        ! Actualizamos la masa total
                        goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                    END IF
                END DO
            END DO
        END DO
    END SELECT
END SUBROUTINE Calculate_masses


!SUBROUTINE Calculate_forces(head, p, a)
!    TYPE(CELL), POINTER :: head
!    type(particle3d), intent(in) :: p(:)
!    type(vector3d), intent(inout) :: a(:)
!    INTEGER :: i
!
!    DO i = 1, SIZE(p)
!        CALL Calculate_forces_aux(i, head, p, a)
!    END DO
!END SUBROUTINE Calculate_forces


RECURSIVE SUBROUTINE Calculate_forces_aux(idx_global,idx_local, tree, p, a)
    TYPE(CELL), POINTER :: tree
    INTEGER, INTENT(IN) :: idx_global  
    INTEGER, INTENT(IN) :: idx_local   
    type(particle3d), intent(in) :: p(:)
    type(vector3d), intent(inout) :: a(:)
    INTEGER :: i, j, k
    real(kind=bit64) :: l, D, r2, r3
    type(vector3d) :: rji

    SELECT CASE (tree%type)
    CASE (1)
        IF (idx_global /= tree%pos) THEN
            rji = tree%c_o_m - p(idx_global)%p
            r2 = mulvv(rji, rji) + epsilon**2 ! including the softening parameter
            r3 = r2 * sqrt(r2)
            a(idx_local) = a(idx_local) + p(tree%pos)%m * rji / r3
        END IF
    CASE (2)
        l = tree%range%max(1) - tree%range%min(1)
        rji = tree%c_o_m - p(idx_global)%p
        r2 = mulvv(rji, rji) + epsilon**2
        D = sqrt(r2)

        IF (l/D < theta) THEN
            r3 = r2 * D
            a(idx_local) = a(idx_local) + tree%mass * rji / r3
        ELSE
            DO i = 1,2
                DO j = 1,2
                    DO k = 1,2
                        IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                            CALL Calculate_forces_aux(idx_global, idx_local, tree%subcell(i,j,k)%ptr, p, a)
                        END IF
                    END DO
                END DO
            END DO
        END IF
    END SELECT
END SUBROUTINE Calculate_forces_aux

SUBROUTINE update_global_tree_and_times(timer_accum)
    real(kind=bit64), intent(inout) :: timer_accum
    integer(kind=8) :: t1, t2
    
    call system_clock(count=t1)
    
    do i = 1, n_local
        tmp_pos_list(3*i-2:3*i) = (/ p(i)%p%x, p(i)%p%y, p(i)%p%z /)
    end do
    
    call MPI_ALLGATHERV(tmp_pos_list, 3*n_local, MPI_DOUBLE_PRECISION, &
                        full_pos_list, 3*counts, 3*offsets, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, ierr)
    
    do i = 1, n_total
        p_global(i)%p%x = full_pos_list(3*i-2)
        p_global(i)%p%y = full_pos_list(3*i-1)
        p_global(i)%p%z = full_pos_list(3*i)
    end do

    CALL Borrar_tree(head)
    CALL Calculate_Ranges(head, p_global, n_total)
    head%type = 0
    CALL Nullify_Pointers(head)
    DO i = 1, n_total
        CALL Find_Cell(head, temp_cell, p_global(i))
        CALL Place_Cell(temp_cell, p_global(i), i)
    END DO
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head, p_global)
    
    call system_clock(count=t2)
    timer_accum = timer_accum + real(t2 - t1, bit64) / real(c_rate, bit64)
END SUBROUTINE update_global_tree_and_times

SUBROUTINE calculate_local_forces_and_times(timer_accum)
    real(kind=bit64), intent(inout) :: timer_accum
    integer(kind=8) :: t1, t2
    integer :: global_idx
    
    call system_clock(count=t1)
    DO i = 1, n_local
        a(i) = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64)
        global_idx = i + offsets(rank+1)
        if (global_idx > n_total) exit
        CALL Calculate_forces_aux(global_idx, i, head, p_global, a)
    END DO
    call system_clock(count=t2)
    timer_accum = timer_accum + real(t2 - t1, bit64) / real(c_rate, bit64)
END SUBROUTINE calculate_local_forces_and_times

END PROGRAM tree

