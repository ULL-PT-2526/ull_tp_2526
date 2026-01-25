program N_body_generator
! This program is from the notes
  use geometry ! here we define the dp
  implicit none
  integer :: I, N
  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed
  real(dp) :: mass, rx, ry, rz
  real(dp) :: dt = 0.001
  real(dp) :: dt_im = 0.1
  real(dp) :: t = 10.0 

  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  print*, "Number of bodies:"
  read*, N
  mass = 1.0 / N


! We create the input file with the particle initial positions

  open(unit=12, file='input.dat', status='replace', action='write')

! We put the input parameters in the file like in the input.dat in ex1

  write(12,'(F6.3)')dt
  write(12,'(F6.3)')dt_im
  write(12,'(F10.3)')t
  write(12,'(I6)')N

  do I = 1, N
    call random_number(rx)
    do
      call random_number(ry)
      if ((rx**2 + ry**2) .le. 1) exit
    end do

    do
      call random_number(rz)
      if ((rx**2 + ry**2 + rz**2) .le. 1) exit
    end do

    ! Finally we write the positions in the input.dat

    write(12,'(F12.6, 6F12.6)') mass, rx, ry, rz, 0.0_dp, 0.0_dp, 0.0_dp! Velocities are 0

  end do

end program N_body_generator