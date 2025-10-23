program leapfrog
    use geometry
    use particle
    implicit none

    integer :: i, n
    real(dp) :: dt, t_end, t, dt_out, t_out
    character(len=100) :: infile
    type(particle3d), allocatable :: p(:)
    type(vector3d), allocatable :: a(:)
    integer :: input_unit, output_unit, ios

    ! Read input file with initial conditions
    call get_command_argument(1, infile)
    if (len_trim(infile) == 0) then
        print *, "Usage: ./ex1 <input_file>"
        stop
    end if

    call read_input(trim(infile), dt, dt_out, t_end, n, p)
    allocate(a(n))

    output_unit = 20
    open(unit=output_unit, file="output.dat", status="replace", action="write")

    ! --- Initialization ---
    call compute_accelerations(p, a, n) ! Compute acceleration
    t_out = 0.0 
    t = 0.0

    ! --- Time integration loop ---
    do while (t <= t_end)
        call update_velocities(p, a, dt, n, 0.5_dp) ! Half step velocity update
        call update_positions(p, dt, n) ! Full step position update
        call compute_accelerations(p, a, n) ! Recompute acceleration
        call update_velocities(p, a, dt, n, 0.5_dp) ! Finish velocity update


        ! Output: write out all the particle positions at intervals in a .dat file
        ! for each time, there will be a line in the file with: time p1x p1y p1z p2x p2y p2z ... pnx pny pnz
        t_out = t_out + dt
        if (t_out >= dt_out) then
            call write_output(output_unit, t, p, n)
            t_out = 0.0
        end if

        t = t + dt
    end do

    close(output_unit)
    print *, "Simulation complete. Results written to output.dat"

    ! Deallocate before exiting
    deallocate(p, a)

contains

    !---------------------------------------------------
    subroutine read_input(filename, dt, dt_out, t_end, n, p)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: dt, dt_out, t_end
        integer, intent(out) :: n
        type(particle3d), allocatable, intent(out) :: p(:)
        integer :: input_unit, ios, i

        input_unit = 10
        open(unit=input_unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening input file ", filename
            stop
        end if

        read (input_unit, *) dt ! time step 
        read (input_unit, *) dt_out ! output interval
        read (input_unit, *) t_end ! end time
        read (input_unit, *) n ! number of particles
        
        allocate(p(n))
        do i = 1, n
            read(input_unit, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, p(i)%v%x, p(i)%v%y, p(i)%v%z
        end do
        close(input_unit)
    end subroutine read_input

    !---------------------------------------------------
    subroutine compute_accelerations(p, a, n)
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)
        integer, intent(in) :: n
        integer :: i, j
        real(dp) :: rs, r3
        type(vector3d) :: rji

        a = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
        do i = 1, n
            do j = i+1, n
                ! for each unique pair of particles (i,j) (no double-counting):
                rji = p(j)%p - p(i)%p ! vector from i to j
                rs = distance(p(i)%p, p(j)%p)
                r3 = rs**3
                a(i) = a(i) + (p(j)%m / r3) * rji
                a(j) = a(j) - (p(i)%m / r3) * rji
            end do
        end do
    end subroutine compute_accelerations

    !---------------------------------------------------
    subroutine update_velocities(p, a, dt, n, factor)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(dp), intent(in) :: dt, factor
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            p(i)%v = p(i)%v + factor * dt * a(i)
        end do
    end subroutine update_velocities

    !---------------------------------------------------
    subroutine update_positions(p, dt, n)
        type(particle3d), intent(inout) :: p(:)
        real(dp), intent(in) :: dt
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            p(i)%p = p(i)%p + dt * p(i)%v
        end do
    end subroutine update_positions

    !---------------------------------------------------
    subroutine write_output(output_unit, t, p, n)
        integer, intent(in) :: output_unit, n
        real(dp), intent(in) :: t
        type(particle3d), intent(in) :: p(:)
        integer :: i

        write(output_unit, '(F10.4)', advance='no') t
        do i = 1, n
            write(output_unit, '(3(1X,F12.6))', advance='no') p(i)%p%x, p(i)%p%y, p(i)%p%z
        end do
        write(output_unit, *)
    end subroutine write_output

end program leapfrog
