! Barnes-Hut algorithm in Fortran90 
! Serial code for N-body problem with Barnes-Hut algorithm
PROGRAM tree
    use geometry
    use particle
    IMPLICIT NONE

    INTEGER :: i,j,k,n
    REAL(dp) :: dt, t_end, t, dt_out, t_out, r2, r3, D, l
    REAL(dp), PARAMETER :: theta = 1.0_dp
    real(dp), parameter :: epsilon = 0.2_dp

    TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: p
    TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a
    TYPE(vector3d) :: rji

    character(len=20) :: mode_name = "SERIAL"
    integer(kind=8)  :: c_ini, c_fin, c_rate
    real(dp) :: t_tree, t_forces, t_update, t_total

    TYPE RANGE
        TYPE(point3d) :: min,max
    END TYPE RANGE

    TYPE CPtr
        TYPE(CELL), POINTER :: ptr
    END TYPE CPtr

    TYPE CELL
        TYPE (RANGE) :: range
        TYPE(point3d) :: part
        INTEGER :: pos
        INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerate
        REAL(dp) :: mass
        TYPE(vector3d) :: c_o_m
        TYPE (CPtr), DIMENSION(2,2,2) :: subcell
    END TYPE CELL

    TYPE (CELL), POINTER :: head, temp_cell

    ! New variables for output
    INTEGER :: ios
    INTEGER :: unit_out

    !! Data reading
    READ*, dt
    READ*, dt_out
    READ*, t_end
    READ*, n

    ALLOCATE(p(n))
    ALLOCATE(a(n))

    DO i = 1, n
        READ*, p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
               p(i)%v%x, p(i)%v%y, p(i)%v%z
    END DO

    !! Head node initialization
    ALLOCATE(head)
    CALL Calculate_ranges(head)
    head%type = 0
    CALL Nullify_Pointers(head)

    !! Initial tree creation
    DO i = 1,n
        CALL Find_Cell(head,temp_cell,p(i)%p)
        CALL Place_Cell(temp_cell,p(i)%p,i)
    END DO

    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head)

    !! Initial accelerations
    a = vector3d(0.0_dp,0.0_dp,0.0_dp)
    CALL Calculate_forces(head)

    ! Initialization of the time variables
    t_tree   = 0.0_dp
    t_forces = 0.0_dp
    t_update = 0.0_dp

    call system_clock(count_rate = c_rate)

    ! Opening the output file
    OPEN(UNIT=10, FILE='output.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT*, 'Error abriendo output.dat'
        STOP
    END IF

    !! Main loop
    t_out = 0.0_dp

    DO t = 0.0_dp, t_end, dt
        call system_clock(count = c_ini) ! INITIAL MEASUREMENT OF TREE + POSITIONS

        ! Velocities and positions update
        DO i=2,n  ! We start at 2 to avoid moving the central object
            p(i)%v = p(i)%v + a(i) * (dt/2.0_dp)
            p(i)%p = p(i)%p + p(i)%v * dt
        END DO

        ! Tree reconstruction
        CALL Borrar_tree(head)
        CALL Calculate_ranges(head)
        head%type = 0
        CALL Nullify_Pointers(head)

        DO i = 1,n
            CALL Find_Cell(head,temp_cell,p(i)%p)
            CALL Place_Cell(temp_cell,p(i)%p,i)
        END DO

        CALL Borrar_empty_leaves(head)
        CALL Calculate_masses(head)

        call system_clock(count=c_fin) 
        t_tree = t_tree + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini) ! Start forces measurement

        ! Computation of accelerations
        a = vector3d(0.0_dp,0.0_dp,0.0_dp)
        CALL Calculate_forces(head)

        call system_clock(count=c_fin)
        t_forces = t_forces + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini) ! Start integration and output measurement

        ! Actualización final de velocidades
        DO i=1,n
            p(i)%v = p(i)%v + a(i) * (dt/2)
        END DO

        t_out = t_out + dt

        ! Writing data to file
        IF (t_out >= dt_out) THEN
            WRITE(10,'(E15.7,1X)', ADVANCE='NO') t

            ! Write positions of all particles 
            ! This format maintains the order x1, y1, z1, x2, y2, z2...
            DO i=1,n
                WRITE(10,'(3E15.7,1X)', ADVANCE='NO') p(i)%p%x, p(i)%p%y, p(i)%p%z
            END DO
            WRITE(10,*)
            t_out=0.0_dp
        END IF

        call system_clock(count=c_fin)
        t_update = t_update + real(c_fin - c_ini, dp) / real(c_rate, dp)
    END DO

    CLOSE(10)

    t_total = t_tree + t_forces + t_update

    print*
    print '(A)', "  ===================================================="
    
    ! Title
    print '(A)', "            ANALYSIS OF EXECUTION - Mode: Serial"
    print '(A)', "             Resource usage: 1 Core (Single)"

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



    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_Ranges !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calculates the ranges of the particles in the 
    !! matrix p in the 3 dimensions and stores it in the 
    !! variable pointed to by goal
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE Calculate_Ranges(goal)
            TYPE(CELL), POINTER :: goal
            TYPE(point3d) :: mins, maxs, medios
            REAL(dp) :: span

            mins%x = MINVAL([p(:)%p%x])
            mins%y = MINVAL([p(:)%p%y])
            mins%z = MINVAL([p(:)%p%z])

            maxs%x = MAXVAL([p(:)%p%x])
            maxs%y = MAXVAL([p(:)%p%y])
            maxs%z = MAXVAL([p(:)%p%z])
            ! When calculating span, I add 10% so that the 
            ! particles do not fall exactly on the boundary
            span = MAX(maxs%x - mins%x, MAX(maxs%y - mins%y, maxs%z - mins%z)) * 1.1_dp


            medios%x = (mins%x + maxs%x)/2.0_dp
            medios%y = (mins%y + maxs%y)/2.0_dp
            medios%z = (mins%z + maxs%z)/2.0_dp

            goal%range%min = medios - vector3d(span/2, span/2, span/2)
            goal%range%max = medios + vector3d(span/2, span/2, span/2)
        END SUBROUTINE Calculate_Ranges
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Find_Cell !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Finds the cell where we will place the particle. 
    !! If the cell we are considering has no 
    !! particle or has one particle, it is this cell where 
    !! we will place the particle. 
    !! If the cell we are considering is a "conglomerate", 
    !! we use the BELONGS function to find which of the 8 
    !! possible subcells it belongs to and call Find_Cell 
    !! again with this subcell. 
    !! 
    !! NOTE: When a "conglomerate" cell is created, all 
    !! 8 subcells are created, so we can assume they always 
    !! exist. Empty cells are deleted at the very end, once 
    !! the entire tree has already been created.
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
        TYPE(CELL), POINTER :: root, goal, temp
        TYPE(point3d), INTENT(IN) :: part
        INTEGER :: i,j,k

        SELECT CASE (root%type)
        CASE (2)
            out: DO i = 1,2
                DO j = 1,2
                    DO k = 1,2
                        IF (Belongs(part, root%subcell(i,j,k)%ptr)) THEN
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Place_Cell !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Executed after Find_Cell, in the cell that 
    !! function returns, so it is always a cell of 
    !! type 0 (no particle) or type 1 (with one particle). 
    !! If it is a type 1 cell, the cell must be subdivided 
    !! and the two particles (the original and the new one) 
    !! placed in its stead.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
        TYPE(CELL), POINTER :: goal, temp
        TYPE(point3d), INTENT(IN) :: part
        INTEGER, INTENT(IN) :: n

        SELECT CASE (goal%type)
        CASE (0)
            goal%type = 1
            goal%part = part
            goal%pos = n
        CASE (1)
            ! If particles are virtually in the same spot, do not subdivide further
            IF (ABS(goal%part%x - part%x) < epsilon .AND. &
                ABS(goal%part%y - part%y) < epsilon .AND. &
                ABS(goal%part%z - part%z) < epsilon) THEN
                RETURN
            END IF

            CALL Crear_Subcells(goal)
            CALL Find_Cell(goal,temp,part)
            CALL Place_Cell(temp,part,n)
        CASE DEFAULT
            print*,"SHOULD NOT BE HERE. ERROR!"
        END SELECT
    END SUBROUTINE Place_Cell

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Crear_Subcells !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! This function is called from Place_Cell and 
    !! is only called when there is already a particle 
    !! in the cell, thus requiring subdivision. It creates 
    !! 8 subcells that "hang" from goal and places the 
    !! particle that was in goal into the corresponding 
    !! subcell among the 8 new ones created. 
    !! 
    !! To create the subcells, use the functions 
    !! CALCULAR_RANGE, BELONGS, and NULLIFY_POINTERS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Crear_Subcells(goal)
        TYPE(CELL), POINTER :: goal
        TYPE(point3d) :: part
        INTEGER :: i,j,k
        INTEGER, DIMENSION(3) :: octant

        part = goal%part
        goal%type = 2

        DO i=1,2
            DO j=1,2
                DO k=1,2
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Nullify_Pointers !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Simply nullifies the pointers of the 8 
    !! subcells of the "goal" cell. 
    !! 
    !! Used in the main loop and by CREAR_SUBCELLS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Nullify_Pointers(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i,j,k
        DO i=1,2
            DO j=1,2
                DO k=1,2
                    NULLIFY(goal%subcell(i,j,k)%ptr)
                END DO
            END DO
        END DO
    END SUBROUTINE Nullify_Pointers

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Belongs !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Returns TRUE if the particle "part" is 
    !! within the range of the cell "goal" 
    !! 
    !! Used by FIND_CELL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Belongs(part,goal)
        TYPE(point3d), INTENT(IN) :: part
        TYPE(CELL), POINTER :: goal
        LOGICAL :: Belongs

        IF (part%x >= goal%range%min%x .AND. part%x <= goal%range%max%x .AND. &
            part%y >= goal%range%min%y .AND. part%y <= goal%range%max%y .AND. &
            part%z >= goal%range%min%z .AND. part%z <= goal%range%max%z) THEN
            Belongs = .TRUE.
        ELSE
            Belongs = .FALSE.
        END IF
    END FUNCTION Belongs

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calcular_Range !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Given an octant "octant" (1,1,1, 1,1,2 ... 2,2,2), 
    !! calculates its ranges based on the ranges of "goal". 
    !! If "what" = 0, calculates minimums. If what = 1, 
    !! calculates maximums.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Calcular_Range(what,goal,octant)
        INTEGER, INTENT(IN) :: what
        TYPE(CELL), POINTER :: goal
        INTEGER, DIMENSION(3), INTENT(IN) :: octant
        TYPE(point3d) :: Calcular_Range, mid

        mid%x = (goal%range%min%x + goal%range%max%x)/2.0_dp
        mid%y = (goal%range%min%y + goal%range%max%y)/2.0_dp
        mid%z = (goal%range%min%z + goal%range%max%z)/2.0_dp

        SELECT CASE (what)
        CASE (0)   ! min
            Calcular_Range%x = MERGE(goal%range%min%x, mid%x, octant(1)==1)
            Calcular_Range%y = MERGE(goal%range%min%y, mid%y, octant(2)==1)
            Calcular_Range%z = MERGE(goal%range%min%z, mid%z, octant(3)==1)

        CASE (1)   ! max
            Calcular_Range%x = MERGE(mid%x, goal%range%max%x, octant(1)==1)
            Calcular_Range%y = MERGE(mid%y, goal%range%max%y, octant(2)==1)
            Calcular_Range%z = MERGE(mid%z, goal%range%max%z, octant(3)==1)
        END SELECT
    END FUNCTION Calcular_Range


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Borrar_empty_leaves !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Called once the tree is complete to delete 
    !! (DEALLOCATE) empty cells (i.e., those without 
    !! a particle).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
        TYPE(CELL),POINTER :: goal
        INTEGER :: i,j,k

        IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
            DO i = 1,2
                DO j = 1,2
                    DO k = 1,2
                        CALL Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                        IF (goal%subcell(i,j,k)%ptr%type == 0) THEN
                            DEALLOCATE (goal%subcell(i,j,k)%ptr)
                        END IF
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE Borrar_empty_leaves

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Borrar_tree !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Deletes the entire tree except for the "head". 
    !! 
    !! The tree must be continuously regenerated, 
    !! so we must delete the old one to avoid memory leaks.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Borrar_tree(goal)
        TYPE(CELL),POINTER :: goal
        INTEGER :: i,j,k
        
        DO i = 1,2
            DO j = 1,2
                DO k = 1,2
                    IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                        CALL Borrar_tree(goal%subcell(i,j,k)%ptr)
                        DEALLOCATE (goal%subcell(i,j,k)%ptr)
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE Borrar_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_masses !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calculates the mass and center-of-mass for all 
    !! cells hanging from "goal".
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Calculate_masses(goal)
        TYPE(CELL), POINTER :: goal
        INTEGER :: i,j,k
        TYPE(particle3d) :: part
        TYPE(vector3d) :: c_o_m
        REAL(dp) :: total_mass

        goal%mass = 0.0_dp
        goal%c_o_m = vector3d(0.0_dp,0.0_dp,0.0_dp)

        SELECT CASE(goal%type)
        CASE(1)
            goal%mass = p(goal%pos)%m
            goal%c_o_m = vector3d(p(goal%pos)%p%x, p(goal%pos)%p%y, p(goal%pos)%p%z)
        CASE(2)
            DO i=1,2
                DO j=1,2
                    DO k=1,2
                        IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                            CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
                            total_mass = goal%mass
                            goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                            goal%c_o_m = (total_mass * goal%c_o_m + &
                                          goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
                        END IF
                    END DO
                END DO
            END DO
        END SELECT
    END SUBROUTINE Calculate_masses

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calculates the forces of all particles against "head". 
    !! Uses the function Calculate_forces_aux which 
    !! actually performs the calculations for each particle.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Calculate_forces(head)
    TYPE(CELL), POINTER :: head
    INTEGER :: i

    DO i=1, n
        CALL Calculate_forces_aux(p(i), a(i), head)
    END DO
    END SUBROUTINE Calculate_forces

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces_aux !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Given a particle "goal_particle", calculates the forces 
    !! exerted on it by the cell "tree". If "tree" is a 
    !! cell containing a single particle, the case is simple, 
    !! as they are treated as two particles. 
    !! 
    !! If "tree" is a conglomerate cell, we first check 
    !! if l/D < theta. That is, if the cell side (l) 
    !! divided by the distance from the goal particle to 
    !! the center_of_mass of the tree cell (D) is less than theta. 
    !! If so, we treat the cell as a single particle. 
    !! If it is not less than theta, then we must consider 
    !! all subcells of tree and recursively call 
    !! Calculate_forces_aux for each of them.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RECURSIVE SUBROUTINE Calculate_forces_aux(goal_particle, goal_accel, tree)
        TYPE(particle3d), INTENT(IN) :: goal_particle
        TYPE(vector3d), INTENT(INOUT) :: goal_accel
        TYPE(CELL), POINTER :: tree
        INTEGER :: i,j,k
        TYPE(vector3d) :: rji
        REAL(dp) :: r2, r3, D, l

        SELECT CASE(tree%type)
        CASE(1)
            ! Verify that it is not the same particle by position
            IF (ABS(goal_particle%p%x - tree%c_o_m%x) > 1.e-10_dp .OR. &
                ABS(goal_particle%p%y - tree%c_o_m%y) > 1.e-10_dp .OR. &
                ABS(goal_particle%p%z - tree%c_o_m%z) > 1.e-10_dp) THEN

                rji = tree%c_o_m - goal_particle%p
                r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2
                r3 = r2 * SQRT(r2)
                ! NOW, WE SUM THE ACCELERATION
                goal_accel%x = goal_accel%x + (tree%mass * rji%x) / r3
                goal_accel%y = goal_accel%y + (tree%mass * rji%y) / r3
                goal_accel%z = goal_accel%z + (tree%mass * rji%z) / r3
            END IF

        CASE(2)
            l = tree%range%max%x - tree%range%min%x
            rji = tree%c_o_m - goal_particle%p
            r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2
            D = SQRT(r2)

            IF (l/D < theta) THEN
                r3 = r2 * D
                goal_accel%x = goal_accel%x + (tree%mass * rji%x) / r3
                goal_accel%y = goal_accel%y + (tree%mass * rji%y) / r3
                goal_accel%z = goal_accel%z + (tree%mass * rji%z) / r3
            ELSE
                DO i=1,2
                    DO j=1,2
                        DO k=1,2
                            IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                                CALL Calculate_forces_aux(goal_particle, goal_accel, tree%subcell(i,j,k)%ptr)
                            END IF
                        END DO
                    END DO
                END DO
            END IF
        END SELECT
    END SUBROUTINE Calculate_forces_aux

END PROGRAM tree