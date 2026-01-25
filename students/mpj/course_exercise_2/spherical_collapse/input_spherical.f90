PROGRAM generador_colapso
    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(p=15, r=307)

    INTEGER :: i, n, k
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    INTEGER :: values(1:8)
    REAL(dp) :: mass, rx, ry, rz, radio_sq, r_min, r_max
    REAL(dp) :: dt, dt_out, t_end

    ! --- Parameters of simulation ---
    dt = 0.002_dp       ! Time step
    dt_out = 1.0_dp     ! Output frequency
    t_end = 5.0_dp      ! Total time
    n = 50000           ! Number of particles
    mass = 1.0_dp       ! Total mass of the system
    
    r_min = 0.5_dp      ! Internal radius
    r_max = 1.0_dp      ! External radius
    ! -----------------------------------

    CALL date_and_time(values=values)
    CALL random_seed(size=k)
    ALLOCATE(seed(k))
    seed(:) = values(8)
    CALL random_seed(put=seed)

    OPEN(unit=10, file="input_spherical.dat", status="replace", action="write")

    WRITE(10,'(F10.5)') dt
    WRITE(10,'(F10.5)') dt_out
    WRITE(10,'(F10.5)') t_end
    WRITE(10,'(I8)') n

    DO i = 1, n
        DO
            ! Generating coordinates between -r_max and r_max
            CALL random_number(rx); rx = (rx * 2.0_dp - 1.0_dp) * r_max
            CALL random_number(ry); ry = (ry * 2.0_dp - 1.0_dp) * r_max
            CALL random_number(rz); rz = (rz * 2.0_dp - 1.0_dp) * r_max
            
            radio_sq = rx**2 + ry**2 + rz**2
            
            ! The particle is only accepted if it is inside the spherical shell (r_min < r < r_max)
            IF (radio_sq <= r_max**2 .AND. radio_sq >= r_min**2) EXIT
        END DO

        ! Format: Mass, Pos(x,y,z), Vel(x,y,z)
        WRITE(10,'(F12.8, 3F12.8, 3F12.8)') mass/n, rx, ry, rz, 0.0_dp, 0.0_dp, 0.0_dp
    END DO

    CLOSE(10)
    PRINT*, "File input.txt generated succesfully for the spherical collapse."

END PROGRAM generador_colapso