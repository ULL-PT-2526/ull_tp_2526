PROGRAM galaxia_bh
  IMPLICIT NONE

  ! Precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

  INTEGER :: i, n, outunit
  REAL(dp) :: mass_disk, mass_center, r, theta, dist, v_mag
  REAL(dp) :: rx, ry, rz, vx, vy, vz
  REAL(dp) :: dt, dt_out, t_end

  ! First lines
  dt = 0.01_dp
  dt_out = 0.1_dp
  t_end = 30.0_dp
  n = 10000
  
  mass_center = 100.0_dp  ! Central object
  mass_disk = 0.001_dp    ! Stars

  OPEN(NEWUNIT=outunit, FILE="input_rotation.dat", STATUS="REPLACE", ACTION="WRITE")
  
  ! Writing the first lines in the input file
  WRITE(outunit,*) dt
  WRITE(outunit,*) dt_out
  WRITE(outunit,*) t_end
  WRITE(outunit,*) n         

  ! --- Central object ---
  ! Mass, px, py, pz, vx, vy, vz
  WRITE(outunit,'(7F14.8)') mass_center, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp

  ! --- DISK ---
  DO i = 2, n
     CALL RANDOM_NUMBER(r)
     ! Distance between 10.0 and 20.0
     dist = 10.0_dp + (r * 10.0_dp) 
     
     CALL RANDOM_NUMBER(theta)
     theta = theta * 2.0_dp * 3.1415926535_dp
     
     rx = dist * COS(theta)
     ry = dist * SIN(theta)
     rz = (r - 0.5_dp) * 0.1_dp ! Very thin disk

     ! CIRCULAR VELOCITY
     v_mag = SQRT(mass_center * dist**2 / (dist**2 + 0.2_dp**2)**1.5_dp)
     
     vx = -v_mag * SIN(theta)
     vy =  v_mag * COS(theta)
     vz = 0.0_dp

     WRITE(outunit,'(7F14.8)') mass_disk, rx, ry, rz, vx, vy, vz
  END DO

  CLOSE(outunit)
  PRINT*, "input_rotation.dat generated."
END PROGRAM galaxia_bh