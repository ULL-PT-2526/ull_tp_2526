! Barnes-Hut algorithm in Fortran90 parallelized with MPI
PROGRAM tree_mpi
    use geometry
    use particle
    use mpi
    IMPLICIT NONE

    ! MPI variables
    INTEGER :: ierr, my_rank, n_procs, i_start, i_end
    
    ! Simulation variables
    INTEGER :: i, j, k, n
    REAL(dp) :: dt, t_end, t, dt_out, t_out
    REAL(dp), PARAMETER :: theta = 1.0_dp
    REAL(dp), PARAMETER :: epsilon = 0.2_dp

    TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: p
    TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a
    
    ! Barnes-Hut Tree Structures
    TYPE RANGE
        TYPE(point3d) :: min, max
    END TYPE RANGE

    TYPE CPtr
        TYPE(CELL), POINTER :: ptr
    END TYPE CPtr

    TYPE CELL
        TYPE (RANGE) :: range
        TYPE(point3d) :: part
        INTEGER :: pos
        INTEGER :: type
        REAL(dp) :: mass
        TYPE(vector3d) :: c_o_m
        TYPE (CPtr), DIMENSION(2,2,2) :: subcell
    END TYPE CELL

    TYPE (CELL), POINTER :: head, temp_cell

    ! Variables for controlling time
    integer(kind=8) :: c_ini, c_fin, c_rate
    real(dp) :: t_tree, t_forces, t_update, t_total
    INTEGER :: ios

    ! MPI initialization
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr)

    ! reading data
    IF (my_rank == 0) THEN
        READ*, dt
        READ*, dt_out
        READ*, t_end
        READ*, n
    END IF

    ! Broadcast basic parameters to all processes
    CALL MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ALLOCATE(p(n), a(n))

    IF (my_rank == 0) THEN
        DO i = 1, n
            READ*, p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
                   p(i)%v%x, p(i)%v%y, p(i)%v%z
        END DO
    END IF

    ! Initial syncronization
    CALL Sincronizar_Particulas_MPI(p, n)

    ! Load balancing in the parallelization
    i_start = (my_rank * n / n_procs) + 1
    i_end   = ((my_rank + 1) * n / n_procs)
    if (my_rank == n_procs - 1) i_end = n

    ALLOCATE(head)

    t_tree = 0.0
    t_forces = 0.0
    t_update = 0.0
    call system_clock(count_rate = c_rate)

    IF (my_rank == 0) THEN
        OPEN(UNIT=10, FILE='output_mpi.dat', STATUS='REPLACE', ACTION='WRITE')
    END IF

    !! Main loop
    t_out = 0.0_dp
    DO t = 0.0_dp, t_end, dt
        call system_clock(count = c_ini)

        ! Uploading velocities and position
        IF (my_rank == 0) THEN
            DO i = 2, n    ! We skip the 1 (central body)
                p(i)%v = p(i)%v + a(i) * (dt/2.0_dp)
                p(i)%p = p(i)%p + p(i)%v * dt
            END DO
        END IF

        CALL Sincronizar_Particulas_MPI(p, n)

        ! Tree reconstruction
        CALL Borrar_tree(head)
        CALL Calculate_ranges(head)
        head%type = 0
        CALL Nullify_Pointers(head)

        DO i = 1, n
            CALL Find_Cell(head, temp_cell, p(i)%p)
            CALL Place_Cell(temp_cell, p(i)%p, i)
        END DO

        CALL Borrar_empty_leaves(head)
        CALL Calculate_masses(head)

        call system_clock(count=c_fin)
        t_tree = t_tree + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini)

        ! Computation of accelerations in parallel
        a = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        DO i = i_start, i_end
            CALL Calculate_forces_aux(p(i), a(i), head)
        END DO

        ! Full acceleration vector across all processes
        CALL Sincronizar_Aceleraciones_MPI(a, n)

        call system_clock(count=c_fin)
        t_forces = t_forces + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini)

        ! Final velocity update
        IF (my_rank == 0) THEN
            DO i = 1, n
                p(i)%v = p(i)%v + a(i) * (dt/2.0_dp)
            END DO
            
            t_out = t_out + dt

            ! Writing data to output file
            IF (t_out >= dt_out) THEN
                WRITE(10,'(E15.7,1X)', ADVANCE='NO') t
                DO i = 1, n
                    WRITE(10,'(3E15.7,1X)', ADVANCE='NO') p(i)%p%x, p(i)%p%y, p(i)%p%z
                END DO
                WRITE(10,*)
                t_out = 0.0_dp
            END IF
        END IF

        call system_clock(count=c_fin)
        t_update = t_update + real(c_fin - c_ini, dp) / real(c_rate, dp)
    END DO

    CLOSE(10)
    CALL MPI_FINALIZE(ierr)

    IF (my_rank == 0) THEN
        t_total = t_tree + t_forces + t_update

        ! Print performance results for MPI
        print*
        print '(A)', "  ===================================================="
        
        ! Title
        print '(A)', "            ANALYSIS OF EXECUTION - Mode: MPI"
        print '(A, I2)', "             Active MPI processes:     ", n_procs

        print '(A)', "  ===================================================="
        print '(A, F10.3, A)', "   Total execution time:    ", t_total, " s"
        print '(A)', "  ----------------------------------------------------"
        
        ! Breakdown of the phases
        print '(A, F10.3, A, F5.1, A)', "   - Octree reconstruction:  ", t_tree,   " s  (", (t_tree/t_total)*100.0, "%)"
        print '(A, F10.3, A, F5.1, A)', "   - Gravity calculations:   ", t_forces, " s  (", (t_forces/t_total)*100.0, "%)"
        print '(A, F10.3, A, F5.1, A)', "   - Integration and output: ", t_update, " s  (", (t_update/t_total)*100.0, "%)"
        print '(A)', "  ____________________________________________________"
        
        if (t_total > 0) then
            print '(A, F12.2, A)', "   Throughput: ", (real(n, dp) / t_total), " particles/s"
        end if
        print '(A)', "  ===================================================="
        print*
    END IF

CONTAINS

    ! --- INCLUDING SUBROUTINES FOR MPI COMMUNICATION ---

    SUBROUTINE Sincronizar_Particulas_MPI(p_vec, n_particles)
        TYPE(particle3d), DIMENSION(:) :: p_vec
        INTEGER, INTENT(IN) :: n_particles
        INTEGER :: ierr
        CALL MPI_BCAST(p_vec, n_particles * 56, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
    END SUBROUTINE

    SUBROUTINE Sincronizar_Aceleraciones_MPI(a_vec, n_particles)
        TYPE(vector3d), DIMENSION(:) :: a_vec
        INTEGER, INTENT(IN) :: n_particles
        REAL(dp), DIMENSION(3*n_particles) :: buf_send, buf_recv
        INTEGER :: i, ierr
        DO i = 1, n_particles
            buf_send(3*i-2) = a_vec(i)%x
            buf_send(3*i-1) = a_vec(i)%y
            buf_send(3*i)   = a_vec(i)%z
        END DO
        CALL MPI_ALLREDUCE(buf_send, buf_recv, 3*n_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        DO i = 1, n_particles
            a_vec(i)%x = buf_recv(3*i-2)
            a_vec(i)%y = buf_recv(3*i-1)
            a_vec(i)%z = buf_recv(3*i)
        END DO
    END SUBROUTINE

    SUBROUTINE Calculate_Ranges(goal)
        TYPE(CELL), POINTER :: goal
        TYPE(point3d) :: mins, maxs, medios
        REAL(dp) :: span
        mins%x = MINVAL(p%p%x); mins%y = MINVAL(p%p%y); mins%z = MINVAL(p%p%z)
        maxs%x = MAXVAL(p%p%x); maxs%y = MAXVAL(p%p%y); maxs%z = MAXVAL(p%p%z)
        span = MAX(maxs%x-mins%x, MAX(maxs%y-mins%y, maxs%z-mins%z)) * 1.1_dp
        medios%x = (mins%x+maxs%x)/2.0_dp; medios%y = (mins%y+maxs%y)/2.0_dp; medios%z = (mins%z+maxs%z)/2.0_dp
        goal%range%min = medios - vector3d(span/2, span/2, span/2)
        goal%range%max = medios + vector3d(span/2, span/2, span/2)
    END SUBROUTINE

    RECURSIVE SUBROUTINE Find_Cell(root, goal, part)
        TYPE(CELL), POINTER :: root, goal, temp
        TYPE(point3d), INTENT(IN) :: part
        INTEGER :: i, j, k
        IF (root%type == 2) THEN
            DO i=1,2; DO j=1,2; DO k=1,2
                IF (Belongs(part, root%subcell(i,j,k)%ptr)) THEN
                    CALL Find_Cell(root%subcell(i,j,k)%ptr, temp, part)
                    goal => temp; RETURN
                END IF
            END DO; END DO; END DO
        ELSE
            goal => root
        END IF
    END SUBROUTINE

    RECURSIVE SUBROUTINE Place_Cell(goal, part, n_idx)
        TYPE(CELL), POINTER :: goal, temp
        TYPE(point3d), INTENT(IN) :: part
        INTEGER, INTENT(IN) :: n_idx
        IF (goal%type == 0) THEN
            goal%type = 1; goal%part = part; goal%pos = n_idx
        ELSE IF (goal%type == 1) THEN
            IF (MAXVAL(ABS([goal%part%x-part%x, goal%part%y-part%y, goal%part%z-part%z])) < epsilon) RETURN
            CALL Crear_Subcells(goal)
            CALL Find_Cell(goal, temp, part)
            CALL Place_Cell(temp, part, n_idx)
        END IF
    END SUBROUTINE

    SUBROUTINE Crear_Subcells(goal)
        TYPE(CELL), POINTER :: goal
        TYPE(point3d) :: part_old
        INTEGER :: i, j, k, p_old
        part_old = goal%part; p_old = goal%pos; goal%type = 2
        DO i=1,2; DO j=1,2; DO k=1,2
            ALLOCATE(goal%subcell(i,j,k)%ptr)
            goal%subcell(i,j,k)%ptr%range%min = Calcular_Range(0, goal, [i,j,k])
            goal%subcell(i,j,k)%ptr%range%max = Calcular_Range(1, goal, [i,j,k])
            IF (Belongs(part_old, goal%subcell(i,j,k)%ptr)) THEN
                goal%subcell(i,j,k)%ptr%type = 1
                goal%subcell(i,j,k)%ptr%part = part_old
                goal%subcell(i,j,k)%ptr%pos = p_old
            ELSE
                goal%subcell(i,j,k)%ptr%type = 0
            END IF
            CALL Nullify_Pointers(goal%subcell(i,j,k)%ptr)
        END DO; END DO; END DO
    END SUBROUTINE

    FUNCTION Belongs(part, goal)
        TYPE(point3d), INTENT(IN) :: part
        TYPE(CELL), POINTER :: goal
        LOGICAL :: Belongs
        Belongs = (part%x >= goal%range%min%x .AND. part%x <= goal%range%max%x .AND. &
                   part%y >= goal%range%min%y .AND. part%y <= goal%range%max%y .AND. &
                   part%z >= goal%range%min%z .AND. part%z <= goal%range%max%z)
    END FUNCTION

    FUNCTION Calcular_Range(what, goal, octant)
        INTEGER, INTENT(IN) :: what
        TYPE(CELL), POINTER :: goal
        INTEGER, DIMENSION(3), INTENT(IN) :: octant
        TYPE(point3d) :: Calcular_Range, mid
        mid%x = (goal%range%min%x + goal%range%max%x)/2.0_dp
        mid%y = (goal%range%min%y + goal%range%max%y)/2.0_dp
        mid%z = (goal%range%min%z + goal%range%max%z)/2.0_dp
        IF (what == 0) THEN
            Calcular_Range%x = MERGE(goal%range%min%x, mid%x, octant(1)==1)
            Calcular_Range%y = MERGE(goal%range%min%y, mid%y, octant(2)==1)
            Calcular_Range%z = MERGE(goal%range%min%z, mid%z, octant(3)==1)
        ELSE
            Calcular_Range%x = MERGE(mid%x, goal%range%max%x, octant(1)==1)
            Calcular_Range%y = MERGE(mid%y, goal%range%max%y, octant(2)==1)
            Calcular_Range%z = MERGE(mid%z, goal%range%max%z, octant(3)==1)
        END IF
    END FUNCTION

    RECURSIVE SUBROUTINE Calculate_masses(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i, j, k
        REAL(dp) :: m_sub
        goal%mass = 0.0; goal%c_o_m = vector3d(0.0, 0.0, 0.0)
        IF (goal%type == 1) THEN
            goal%mass = p(goal%pos)%m
            goal%c_o_m = vector3d(p(goal%pos)%p%x, p(goal%pos)%p%y, p(goal%pos)%p%z)
        ELSE IF (goal%type == 2) THEN
            DO i=1,2; DO j=1,2; DO k=1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                    CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
                    m_sub = goal%subcell(i,j,k)%ptr%mass
                    goal%c_o_m = (goal%mass * goal%c_o_m + m_sub * goal%subcell(i,j,k)%ptr%c_o_m) / (goal%mass + m_sub)
                    goal%mass = goal%mass + m_sub
                END IF
            END DO; END DO; END DO
        END IF
    END SUBROUTINE

    RECURSIVE SUBROUTINE Calculate_forces_aux(goal_p, goal_a, tree)
        TYPE(particle3d), INTENT(IN) :: goal_p
        TYPE(vector3d), INTENT(INOUT) :: goal_a
        TYPE(CELL), POINTER :: tree
        TYPE(vector3d) :: rji
        REAL(dp) :: r2, D, l
        INTEGER :: i, j, k

        IF (tree%type == 1) THEN
            rji = tree%c_o_m - goal_p%p
            r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2
            IF (r2 > epsilon**2 + 1e-10) THEN
                goal_a = goal_a + (tree%mass * rji) / (r2 * SQRT(r2))
            END IF
        ELSE IF (tree%type == 2) THEN
            l = tree%range%max%x - tree%range%min%x
            rji = tree%c_o_m - goal_p%p
            r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2
            D = SQRT(r2)
            IF (l/D < theta) THEN
                goal_a = goal_a + (tree%mass * rji) / (r2 * D)
            ELSE
                DO i=1,2; DO j=1,2; DO k=1,2
                    IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) CALL Calculate_forces_aux(goal_p, goal_a, tree%subcell(i,j,k)%ptr)
                END DO; END DO; END DO
            END IF
        END IF
    END SUBROUTINE

    SUBROUTINE Nullify_Pointers(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i, j, k
        DO i=1,2; DO j=1,2; DO k=1,2; NULLIFY(goal%subcell(i,j,k)%ptr); END DO; END DO; END DO
    END SUBROUTINE

    RECURSIVE SUBROUTINE Borrar_tree(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i, j, k
        DO i=1,2; DO j=1,2; DO k=1,2
            IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                CALL Borrar_tree(goal%subcell(i,j,k)%ptr)
                DEALLOCATE(goal%subcell(i,j,k)%ptr)
            END IF
        END DO; END DO; END DO
    END SUBROUTINE

    RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i, j, k
        IF (goal%type == 2) THEN
            DO i=1,2; DO j=1,2; DO k=1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                    CALL Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                    IF (goal%subcell(i,j,k)%ptr%type == 0) DEALLOCATE(goal%subcell(i,j,k)%ptr)
                END IF
            END DO; END DO; END DO
        END IF
    END SUBROUTINE

END PROGRAM tree_mpi