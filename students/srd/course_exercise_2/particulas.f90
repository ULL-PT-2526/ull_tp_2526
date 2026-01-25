PROGRAM particulas
  IMPLICIT NONE

  INTEGER :: i, n, k
  INTEGER :: values(1:8)
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  REAL :: mass, rx, ry, rz
  REAL :: dt, dt_out, t_end
  REAL :: R
  INTEGER :: outunit

  ! Inicialización de semilla aleatoria
  CALL date_and_time(values = values)
  CALL random_seed(size = k)
  ALLOCATE(seed(1:k))
  seed(:) = values(8)
  CALL random_seed(put = seed)

  ! Número de partículas
  PRINT*, "Number of bodies?"
  READ*, n

  R = 8.0
  mass = 1.0 / n

  ! Parámetros de tiempo
  dt = 0.01
  dt_out = 0.1
  t_end = 30.0

  ! Abrir archivo para escribir
  OPEN(NEWUNIT=outunit, FILE="input.dat", STATUS="REPLACE", ACTION="WRITE")

  ! Escribir parámetros iniciales
  WRITE(outunit,*) dt
  WRITE(outunit,*) dt_out
  WRITE(outunit,*) t_end
  WRITE(outunit,*) n

  ! Escribir partículas
  DO i = 1, n
     DO
        ! Generamos las 3 coordenadas de golpe en el cubo [-R, R]
        CALL random_number(rx)
        CALL random_number(ry)
        CALL random_number(rz)
        
        rx = (rx * 2.0 - 1.0) * R
        ry = (ry * 2.0 - 1.0) * R
        rz = (rz * 2.0 - 1.0) * R
        
        ! Si el punto cae dentro de la esfera, lo aceptamos
        IF (rx**2 + ry**2 + rz**2 <= R**2) EXIT
     END DO

     WRITE(outunit,'(F10.6,3F12.6,3I2)') mass, rx, ry, rz, 0, 0, 0
  END DO

  CLOSE(outunit)
  PRINT*, "Archivo input.dat generado con exito (Radio =", R, ")"
END PROGRAM particulas


