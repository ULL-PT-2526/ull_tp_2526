program leapfrog
    use iso_fortran_env, only: real64
    use geometry
    use particle
    implicit none

    integer :: i, ierr, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    character(len=300) :: setup_file_name
    type(particle3d), allocatable :: p(:)
    type(vector3d), allocatable :: a(:)

    ! Read setup file name
    print *, 'Insert name of the simulation setup file:'
    read *, setup_file_name

    ! Read input parameters and initialize particles
    call read_setup(setup_file_name, dt, dt_out, t_end, n, p)

    ! Allocate acceleration array
    allocate(a(n), stat=ierr)
    if (ierr /= 0) stop 'Error: allocation of acceleration array failed'

    ! Initial acceleration computation
    call compute_acceleration(p, a, n)

    ! Create and write file for the output of the simulation
    open(unit=20, file='output.dat', status='replace', action='write')

    !======================
    ! Main integration loop
    !======================
    t = 0.0
    t_out = 0.0

    do while (t <= t_end)
        call update_positions(p, a, dt, n)
        call compute_acceleration(p, a, n)
        call update_velocities(p, a, dt, n)

        t_out = t_out + dt
        if (t_out >= dt_out) then
            write(20, *) t, (p(i)%p, i=1,n)
            t_out = 0.0
        end if

        t = t + dt
    end do

    ! Close output file
    close(20)

    ! Deallocate allocated variables
    if (allocated(p)) deallocate(p)
    if (allocated(a)) deallocate(a)

    ! Print final comment to check that everything went well
    print *, 'Simulation completed successfully. Results saved in "output.dat".'

contains
    !====================================================
    ! Subroutine: read_setup
    !====================================================
    subroutine read_setup(filename, dt, dt_out, t_end, n, p)
        character(len=*), intent(in) :: filename
        real(real64), intent(out) :: dt, dt_out, t_end
        integer, intent(out) :: n
        type(particle3d), allocatable, intent(out) :: p(:)
        integer :: i, ierr

        ! Open and read file with the input data
        open(unit=10, file=filename, status='old', action='read')

        read(10, *) dt
        read(10, *) dt_out
        read(10, *) t_end
        read(10, *) n

        allocate(p(n), stat=ierr)
        if (ierr /= 0) stop 'Error: allocation of particle array failed'
        
        do i = 1, n
            read(10, *) p(i)%m, p(i)%p, p(i)%v
        end do
        
        ! Close input file
        close(10)
    end subroutine read_setup

    !====================================================
    ! Subroutine: compute_acceleration
    !====================================================
    subroutine compute_acceleration(p, a, n)
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)
        integer, intent(in) :: n
        type(vector3d) :: rji
        real(real64) :: dist, r3
        integer :: i, j

        a = vector3d(0.0, 0.0, 0.0)  ! Initialize all elements of 'a' to zero

        do i = 1, n
            do j = i+1, n
                rji = p(j)%p - p(i)%p
                dist = distance(p(i)%p, p(j)%p)
                r3 = dist * dist**2
                a(i) = a(i) + p(j)%m * rji / r3  ! Acceleration on particle i due to particle j
                a(j) = a(j) - p(i)%m * rji / r3  ! Acceleration on particle j due to particle i
            end do
        end do
    end subroutine compute_acceleration

    !====================================================
    ! Subroutine: update_positions (half velocity + full position step)
    !====================================================
    subroutine update_positions(p, a, dt, n)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            p(i)%v = p(i)%v + a(i) * (dt / 2.0)
            p(i)%p = p(i)%p + p(i)%v * dt
        end do
    end subroutine update_positions

    !====================================================
    ! Subroutine: update_velocities (second half-step)
    !====================================================
    subroutine update_velocities(p, a, dt, n)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            p(i)%v = p(i)%v + a(i) * (dt / 2.0)
        end do
    end subroutine update_velocities

end program leapfrog
