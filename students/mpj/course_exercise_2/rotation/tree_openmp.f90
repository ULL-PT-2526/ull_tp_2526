! Algoritmo Barnes-Hut en Fortran90
! 5.2 Codigo serie para problema N-body con algoritmo Barnes-Hut
PROGRAM tree_omp
    use geometry
    use particle
    use omp_lib
    IMPLICIT NONE

    character(len=20) :: mode_name = "SERIAL"

    INTEGER :: i,j,k,n
    REAL(dp) :: dt, t_end, t, dt_out, t_out, r2, r3, D, l
    REAL(dp), PARAMETER :: theta = 1.0_dp
    real(dp), parameter :: epsilon = 0.2_dp

    TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: p
    TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a
    TYPE(vector3d) :: rji

    ! Time control variables:
    integer(kind=8)  :: c_ini, c_fin, c_rate
    real(dp) :: t_tree, t_forces, t_update, t_total
    INTEGER :: nt = 1

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
        INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
        REAL(dp) :: mass
        TYPE(vector3d) :: c_o_m
        TYPE (CPtr), DIMENSION(2,2,2) :: subcell
    END TYPE CELL

    TYPE (CELL), POINTER :: head, temp_cell

    ! Nuevas variables para salida
    INTEGER :: ios
    INTEGER :: unit_out

    !! Lectura de datos
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

    !! Inicializacion head node
    ALLOCATE(head)
    CALL Calculate_ranges(head)
    head%type = 0
    CALL Nullify_Pointers(head)

    !! Creacion del arbol inicial
    DO i = 1,n
        CALL Find_Cell(head,temp_cell,p(i)%p)
        CALL Place_Cell(temp_cell,p(i)%p,i)
    END DO

    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head)

    !! Aceleraciones iniciales
    !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i)
    DO i=1,n
        a(i) = vector3d(0.0_dp,0.0_dp,0.0_dp)
    END DO
    !$OMP END PARALLEL DO

    CALL Calculate_forces(head)

    ! Initialization of the time variables
    t_tree   = 0.0_dp
    t_forces = 0.0_dp
    t_update = 0.0_dp

    call system_clock(count_rate = c_rate)
    call system_clock(count = c_ini)


    ! Apertura del archivo de salida
    OPEN(UNIT=10, FILE='output_openmp.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT*, 'Error abriendo output_openmp.dat'
        STOP
    END IF

    !! Bucle principal
    t_out = 0.0_dp

    DO t = 0.0_dp, t_end, dt
        call system_clock(count = c_ini)

        ! Actualización de velocidades y posiciones
        !$OMP PARALLEL DO PRIVATE(i)
        DO i=2,n  ! EMPEZAMOS EN 2 PARA NO MOVER EL AGUJERO NEGRO
            p(i)%v = p(i)%v + a(i) * (dt/2.0_dp)
            p(i)%p = p(i)%p + p(i)%v * dt
        END DO
        !$OMP END PARALLEL DO

        ! Reconstrucción del árbol
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

        call system_clock(count=c_fin) ! Clock stop
        t_tree = t_tree + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini) ! Clock initialization

        ! Cálculo de aceleraciones
        !$omp parallel do private(i)
        DO i = 1,n
            a(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        END DO
        !$omp end parallel do
        CALL Calculate_forces(head)

        call system_clock(count=c_fin) ! Clock stop
        t_forces = t_forces + real(c_fin - c_ini, dp) / real(c_rate, dp)

        call system_clock(count=c_ini) ! Clock initialization

        ! Actualización final de velocidades
        !$OMP PARALLEL DO PRIVATE(i)
        DO i=1,n
            p(i)%v = p(i)%v + a(i) * (dt/2.0_dp)
        END DO
        !$OMP END PARALLEL DO


        t_out = t_out + dt

        ! Escritura de los datos al archivo
        IF (t_out >= dt_out) THEN
            WRITE(10,'(E15.7,1X)', ADVANCE='NO') t

            ! Escribir posiciones de todas las partículas
            ! Esta forma mantiene el orden x1, y1, z1, x2, y2, z2...
            DO i=1,n
                WRITE(10,'(3E15.7,1X)', ADVANCE='NO') p(i)%p%x, p(i)%p%y, p(i)%p%z
            END DO
            WRITE(10,*)
            t_out=0.0_dp
        END IF

        call system_clock(count=c_fin) ! Clock stop
        t_update = t_update + real(c_fin - c_ini, dp) / real(c_rate, dp)

    END DO

    CLOSE(10)

    t_total = t_tree + t_forces + t_update

    ! Print time info
    print*
    print '(A)', "  ===================================================="
    
    ! Centramos el título del modo manualmente
    print '(A)', "          ANALYSIS OF EXECUTION - Mode: OpenMP"
    
    !$ nt = omp_get_max_threads()
    !$ print '(A, I2)', "           Processing threads active: ", nt

    print '(A)', "  ===================================================="
    print '(A, F10.3, A)', "   Total execution time: ", t_total, " s"
    print '(A)', "  ----------------------------------------------------"
    
    ! Desglose de las fases
    print '(A, F10.3, A, F5.1, A)', "   - Octree reconstruction:     ", t_tree,   " s  (", (t_tree/t_total)*100.0, "%)"
    print '(A, F10.3, A, F5.1, A)', "   - Gravity calculations:      ", t_forces, " s  (", (t_forces/t_total)*100.0, "%)"
    print '(A, F10.3, A, F5.1, A)', "   - Integration and output:    ", t_update, " s  (", (t_update/t_total)*100.0, "%)"
    print '(A)', "  ____________________________________________________"
    
    ! Mantenemos solo el Throughput como dato de interés científico
    if (t_total > 0) then
        print '(A, F12.2, A)', "   Throughput: ", (n / t_total), " particles/s"
    end if
    print '(A)', "  ===================================================="
    print*



    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_Ranges !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Calcula los rangos de las partıculas en la
    !! matriz r en las 3 dimensiones y lo pone en la
    !! variable apuntada por goal
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
            ! Al calcular span le sumo un 10% para que las
            ! particulas no caigan justo en el borde
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
    !! Encuentra la celda donde colocaremos la particula.
    !! Si la celda que estamos considerando no tiene
    !! particula o tiene una particula, es esta celda donde
    !! colocaremos la particula.
    !! Si la celda que estamos considerando es un "conglomerado",
    !! buscamos con la funcion BELONGS a que subcelda de las 8
    !! posibles pertenece y con esta subcelda llamamos de nuevo
    !! a Find_Cell
    !!
    !! NOTA: Cuando se crea una celda "conglomerado" se crean las
    !! 8 subceldas, por lo que podemos asumir que siempre existen
    !! las 8. Las celdas vacıas se borran al final del todo, cuando
    !! todo el arbol ha sido ya creado.
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
    !! Se ejecuta tras Find_Cell, en la celda que
    !! esa funcion nos devuelve, por lo que siempre
    !! es una celda de tipo 0 (sin particula) o de tipo 1
    !! (con una particula). En el caso de que es una celda
    !! de tipo 1 habra que subdividir la celda y poner en
    !! su lugar las dos particulas (la que originalmente
    !! estaba, y la nueva).
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
            ! Si las partículas están virtualmente en el mismo sitio, no subdividas más
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
    !! Esta funcion se llama desde Place_Cell y
    !! solo se llama cuando ya hay una particula
    !! en la celda, con lo que la tenemos que
    !! subdividir. Lo que hace es crear 8 subceldas
    !! que "cuelgan" de goal y la particula que
    !! estaba en goal la pone en la subcelda que
    !! corresponda de la 8 nuevas creadas.
    !!
    !! Para crear las subceldas utilizar las funciones
    !! CALCULAR_RANGE, BELONGS y NULLIFY_POINTERS
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
    !! Simplemente me NULLIFYca los punteros de
    !! las 8 subceldas de la celda "goal"
    !!
    !! Se utiliza en el bucle principal y por
    !! CREAR_SUBCELLS
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
    !! Devuelve TRUE si la particula "part" esta
    !! dentro del rango de la celda "goal"
    !!
    !! Utilizada por FIND_CELL
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
    !! Dado un octante "otctant" (1,1,1, 1,1,2 ... 2,2,2),
    !! calcula sus rangos en base a los rangos de
    !! "goal". Si "what" = 0 calcula los minimos. Si what=1
    !! calcula los maximos.
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
    !! Se llama una vez completado el arbol para
    !! borrar (DEALLOCATE) las celdas vacıas (i.e.
    !! sin partıcula).
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
    !! Borra el arbol completo, excepto la "head".
    !!
    !! El arbol se ha de regenerar continuamente,
    !! por lo que tenemos que borrar el antiguo
    !! para evitar "memory leaks".
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
    !! Nos calcula para todas las celdas que cuelgan
    !! de "goal" su masa y su center-of-mass.
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
    !! Calcula las fuerzas de todas las particulas contra "head".
    !! Se sirve de la funcion Calculate_forces_aux que es la
    !! que en realidad hace los calculos para cada particula
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Calculate_forces(head)
        TYPE(CELL), POINTER :: head
        INTEGER :: i

        !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i) SCHEDULE(dynamic)
        DO i=1, n
            CALL Calculate_forces_aux(p(i), a(i), head)
        END DO
        !$OMP END PARALLEL DO

    END SUBROUTINE Calculate_forces


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Calculate_forces_aux !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Dada una particula "goal" calcula las fuerzas
    !! sobre ella de la celda "tree". Si "tree" es una
    !! celda que contiene una sola particula el caso
    !! es sencillo, pues se tratan de dos particulas.
    !!
    !! Si "tree" es una celda conglomerado, hay que ver primero
    !! si l/D < theta. Es decir si el lado de la celda (l)
    !! dividido entre la distancia de la particula goal
    !! al center_of_mass de la celda tree (D) es menor que theta.
    !! En caso de que asi sea, tratamos a la celda como una
    !! sola particula. En caso de que no se menor que theta,
    !! entonces tenemos que considerar todas las subceldas
    !! de tree y para cada una de ellas llamar recursivamente
    !! a Calculate_forces_aux
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
            ! Verificamos que no sea la misma partícula por posición
            IF (ABS(goal_particle%p%x - tree%c_o_m%x) > 1.e-10_dp .OR. &
                ABS(goal_particle%p%y - tree%c_o_m%y) > 1.e-10_dp .OR. &
                ABS(goal_particle%p%z - tree%c_o_m%z) > 1.e-10_dp) THEN

                rji = tree%c_o_m - goal_particle%p
                r2 = rji%x**2 + rji%y**2 + rji%z**2 + epsilon**2
                r3 = r2 * SQRT(r2)
                ! AHORA SUMA A LA ACELERACIÓN
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

END PROGRAM tree_omp