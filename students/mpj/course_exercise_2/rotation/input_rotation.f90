PROGRAM galaxia_bh
  IMPLICIT NONE
  ! Definimos el tipo de precisión igual que en tu código
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER :: i, n, outunit
  REAL(dp) :: mass_disk, mass_center, r, theta, dist, v_mag
  REAL(dp) :: rx, ry, rz, vx, vy, vz
  REAL(dp) :: dt, dt_out, t_end

  ! 1. Configuración idéntica a tu READ
  n = 10000
  dt = 0.01_dp
  dt_out = 0.1_dp
  t_end = 30.0_dp
  
  mass_center = 100.0_dp  ! Agujero negro central
  mass_disk = 0.001_dp   ! Estrellas

  OPEN(NEWUNIT=outunit, FILE="input_rotation.dat", STATUS="REPLACE", ACTION="WRITE")
  
  ! 2. ESCRITURA EN EL ORDEN QUE PIDE TU PROGRAMA:
  WRITE(outunit,*) dt        ! Primero dt
  WRITE(outunit,*) dt_out    ! Segundo dt_out
  WRITE(outunit,*) t_end     ! Tercero t_end
  WRITE(outunit,*) n         ! Y AL FINAL el entero n

  ! --- 1. EL AGUJERO NEGRO ---
  ! Masa, px, py, pz, vx, vy, vz
  WRITE(outunit,'(7F14.8)') mass_center, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp

  ! --- 2. EL DISCO ---
  DO i = 2, n
     CALL RANDOM_NUMBER(r)
     ! Distancia entre 10.0 y 20.0
     dist = 10.0_dp + (r * 10.0_dp) 
     
     CALL RANDOM_NUMBER(theta)
     theta = theta * 2.0_dp * 3.1415926535_dp
     
     rx = dist * COS(theta)
     ry = dist * SIN(theta)
     rz = (r - 0.5_dp) * 0.1_dp ! Un disco muy fino

     ! VELOCIDAD CIRCULAR CORREGIDA (considerando el epsilon del árbol)
     ! Si usas epsilon = 0.2 en el árbol, úsalo aquí también.
     v_mag = SQRT(mass_center * dist**2 / (dist**2 + 0.2_dp**2)**1.5_dp)
     
     vx = -v_mag * SIN(theta)
     vy =  v_mag * COS(theta)
     vz = 0.0_dp

     WRITE(outunit,'(7F14.8)') mass_disk, rx, ry, rz, vx, vy, vz
  END DO

  CLOSE(outunit)
  PRINT*, "input_rotation.dat generado respetando tu orden de lectura."
END PROGRAM galaxia_bh