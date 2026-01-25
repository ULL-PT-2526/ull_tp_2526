PROGRAM tree
use geometry
use particle
!$ use omp_lib ! OpenMP library. Only activates with the -fopenmp when compiling
IMPLICIT NONE

! Assume serial. Only applies the OpenMP changes with the flag at compilation
character(len=20) :: mode_name = "SERIAL" ! Relevant when printing the performance
INTEGER :: nt = 1

INTEGER :: i,j,k,n
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

!! Lectura de datos
READ*, dt
READ*, dt_out
READ*, t_end
READ*, n

OPEN(unit = out_unit, file = "output.dat", status = "replace", action = "write") ! Opening the file


! Changes in the variables allocated
ALLOCATE(p(n))
ALLOCATE(a(n))

!Reading all the "p" properties from the module
DO i = 1, n
    READ*, p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
                 p(i)%v%x, p(i)%v%y, p(i)%v%z
END DO


!! Inicialización head node
ALLOCATE(head)
CALL Calculate_Ranges(head,p,n)
head%type = 0
CALL Nullify_Pointers(head)

!! Creación del árbol inicial
DO i = 1,n
    CALL Find_Cell(head, temp_cell, p(i))
    CALL Place_Cell(temp_cell, p(i), i)
END DO

CALL Borrar_empty_leaves(head)
CALL Calculate_masses(head,p)

!! Calcular aceleraciones iniciales
DO i = 1,n
    a(i) = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64) ! Now using the vector3d acceleration
END DO

CALL Calculate_forces(head,p,a)

! Initialization of the time variables
t_tree   = 0.0_bit64
t_forces = 0.0_bit64
t_update = 0.0_bit64
! Obtaining the tick frequency
call system_clock(count_rate=c_rate)

!! Bucle principal
t = 0.0_bit64
t_out = 0.0_bit64
DO WHILE (t <= t_end)

    call system_clock(count=c_ini) ! Clock initialization
    !$omp parallel do private(i)
    DO i = 1,n
        p(i)%v = p(i)%v + a(i) * (dt/2.0_bit64)
        p(i)%p = p(i)%p + p(i)%v * dt
    END DO
    !$omp end parallel do
    call system_clock(count=c_fin) ! Clock stop
    t_update = t_update + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

    call system_clock(count=c_ini) ! Clock initialization
    ! We DON'T paralellize the tree construction using OpenMP
    ! Multiple cores would try to construct the same cell at the same time -> program breaks

    !! Las posiciones han cambiado, reinicializamos el árbol
    CALL Borrar_tree(head)
    CALL Calculate_Ranges(head,p,n)
    head%type = 0
    CALL Nullify_Pointers(head)

    DO i = 1,n
        CALL Find_Cell(head, temp_cell, p(i)) ! Using particle
        CALL Place_Cell(temp_cell, p(i), i)
    END DO

    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head,p)
    call system_clock(count=c_fin) ! Clock stop
    t_tree = t_tree + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

    call system_clock(count=c_ini) ! Clock initialization
    !$omp parallel do private(i)
    DO i = 1,n
        a(i) = vector3d(0.0_bit64, 0.0_bit64, 0.0_bit64)
    END DO
    !$omp end parallel do

    CALL Calculate_forces(head,p,a)
    call system_clock(count=c_fin) ! Clock stop
    t_forces = t_forces + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

    call system_clock(count=c_ini) ! Clock initialization
    !$omp parallel do private(i)
    DO i = 1,n
        p(i)%v = p(i)%v + a(i) * (dt/2.0_bit64)
    END DO
    !$omp end parallel do

    t_out = t_out + dt
    IF (t_out >= dt_out) THEN
        WRITE(out_unit, *) t, ( p(i)%p%x, p(i)%p%y, p(i)%p%z, i = 1, n )
        t_out = 0.0
    END IF

    t = t + dt
    call system_clock(count=c_fin) ! Clock stop
    t_update = t_update + real(c_fin - c_ini, bit64) / real(c_rate, bit64)

END DO

CLOSE(out_unit)

t_total = t_tree + t_forces + t_update

! Print time info
PRINT*, ""
PRINT*, "==============================================="
!$ mode_name = "OpenMP"
!$ nt = omp_get_max_threads()
PRINT*, "          PERFORMANCE REPORT (", trim(mode_name), ") "
!$ PRINT*, "          THREADS USED: ", nt
    
    PRINT*, "==============================================="
PRINT '(A, F12.4, A)', " Total Calculation Time:  ", t_total, " s"
PRINT*, "-----------------------------------------------"
PRINT '(A, F12.4, A, F6.2, A)', " 1. Tree Management:      ", t_tree,   " s | ", (t_tree/t_total)*100.0, "%"
PRINT '(A, F12.4, A, F6.2, A)', " 2. Force Calculation:    ", t_forces, " s | ", (t_forces/t_total)*100.0, "%"
PRINT '(A, F12.4, A, F6.2, A)', " 3. Integration & I/O:    ", t_update, " s | ", (t_update/t_total)*100.0, "%"
PRINT*, "===============================================" 

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


SUBROUTINE Calculate_forces(head, p, a)
    TYPE(CELL), POINTER :: head
    type(particle3d), intent(in) :: p(:)
    type(vector3d), intent(inout) :: a(:)
    INTEGER :: i
    
    ! Dynamic parallelization: particle assignment meanwhile the cores are free
    !$omp parallel do private(i) schedule(dynamic)
    DO i = 1, SIZE(p)
        CALL Calculate_forces_aux(i, head, p, a)
    END DO
    !$omp end parallel do
END SUBROUTINE Calculate_forces


RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, p, a)
    TYPE(CELL), POINTER :: tree
    INTEGER :: goal
    type(particle3d), intent(in) :: p(:)
    type(vector3d), intent(inout) :: a(:)
    INTEGER :: i, j, k
    real(kind=bit64) :: l, D, r2, r3
    type(vector3d) :: rji

    SELECT CASE (tree%type)
    CASE (1)
        IF (goal /= tree%pos) THEN
            rji = tree%c_o_m - p(goal)%p
            r2 = mulvv(rji, rji) + epsilon**2 ! including the softening parameter
            r3 = r2 * sqrt(r2)
            a(goal) = a(goal) + p(tree%pos)%m * rji / r3
        END IF
    CASE (2)
        l = tree%range%max(1) - tree%range%min(1)
        rji = tree%c_o_m - p(goal)%p
        r2 = mulvv(rji, rji) + epsilon**2
        D = sqrt(r2)

        IF (l/D < theta) THEN
            r3 = r2 * D
            a(goal) = a(goal) + tree%mass * rji / r3
        ELSE
            DO i = 1,2
                DO j = 1,2
                    DO k = 1,2
                        IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                            CALL Calculate_forces_aux(goal, tree%subcell(i,j,k)%ptr, p, a)
                        END IF
                    END DO
                END DO
            END DO
        END IF
    END SELECT
END SUBROUTINE Calculate_forces_aux


END PROGRAM tree

