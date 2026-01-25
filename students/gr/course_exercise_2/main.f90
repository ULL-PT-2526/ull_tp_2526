program main
    use iso_fortran_env, only: real64

#ifdef _OPENMP
    use omp_lib
#endif

    use barnes_hut
    use geometry
    use particle
    implicit none

    integer :: ierr, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    character(len=300) :: setup_file_name

    type(particle3d), allocatable :: p(:)
    type(vector3d), allocatable :: a(:)
    type(cell), pointer :: head

    !---- TIMING ----
    real(real64) :: time_start, time_end, elapsed
    character(len=10) :: mode
    integer :: nthreads

    ! Read setup file name
    print *, 'Insert name of the simulation setup file:'
    read *, setup_file_name

    ! Read input parameters and initialize particles
    call read_setup(setup_file_name, dt, dt_out, t_end, n, p)

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
    if (nthreads > 1) then
        mode = "openmp"
    else
        mode = "serial"
    end if
#else
    nthreads = 1
    mode = "serial"
#endif

    ! Allocate acceleration array
    allocate(a(n), stat=ierr)
    if (ierr /= 0) stop 'Error: allocation of acceleration array failed'

    ! Initialize Barnes–Hut and compute initial accelerations
    call initialize_bh(head, p, a)

    ! Create and write file for the output of the simulation
    ! with a dynamic name
    block
        character(len=400) :: out_name
        ! Write variable out_name: folder/output_mode_nVALUE.dat
        write(out_name, '(A, A, A, I0, A)') &
            "outputs/output_", trim(mode), "_n", n, ".dat"
        
        call execute_command_line("mkdir -p outputs")
        open(unit=20, file=trim(out_name), status='replace', action='write')
        print *, 'Output file created: ', trim(out_name)
    end block

    !======================
    ! TIMING START
    !======================
    call cpu_time(time_start)

    !======================
    ! Main integration loop
    !======================
    t = 0.0_real64
    t_out = 0.0_real64

    do while (t <= t_end)

        call update_positions(p, a, dt)
        call rebuild_tree_and_forces(head, p, a)
        call update_velocities(p, a, dt)

        t_out = t_out + dt
        if (t_out >= dt_out) then
            call write_output(20, t, p)
            t_out = 0.0_real64
        end if

        t = t + dt
    end do

    !======================
    ! TIMING END
    !======================
    call cpu_time(time_end)
    elapsed = time_end - time_start

    ! Close output file
    close(20)

    !---- SAVE THE TIME ----
    open(unit=99, file='timing.dat', position='append', action='write')

    write(99,'(A,1X,I8,1X,I4,1X,F12.6)') &
     trim(mode), n, nthreads, elapsed

    close(99)

    ! Deallocate allocated variables
    if (allocated(p)) deallocate(p)
    if (allocated(a)) deallocate(a)

    ! Delete and deallocate tree
    call borrar_tree(head)

    print *, 'Simulation time:', elapsed, 'seconds'
    print *, 'Results saved in "outputs" folder'
    print *, 'Timing saved in "timing.dat"'

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
        
        close(10)
    end subroutine read_setup

    !====================================================
    ! Subroutine: initialize_bh
    !====================================================
    subroutine initialize_bh(head, p, a)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)

        call bh_init(head, p)
        call bh_build_tree(head, p)
        call bh_compute_forces(head, p, a)
    end subroutine initialize_bh

    !====================================================
    ! Subroutine: update_positions
    !====================================================
    subroutine update_positions(p, a, dt)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt

        p%v = p%v + a * (dt * 0.5_real64)
        p%p = p%p + p%v * dt
    end subroutine update_positions

    !====================================================
    ! Subroutine: rebuild_tree_and_forces
    !====================================================
    subroutine rebuild_tree_and_forces(head, p, a)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)

        call bh_destroy(head, p)
        call bh_build_tree(head, p)
        call bh_compute_forces(head, p, a)
    end subroutine rebuild_tree_and_forces

    !====================================================
    ! Subroutine: update_velocities
    !====================================================
    subroutine update_velocities(p, a, dt)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt

        p%v = p%v + a * (dt * 0.5_real64)
    end subroutine update_velocities

    !====================================================
    ! Subroutine: write_output
    !====================================================
    subroutine write_output(unit, t, p)
        integer, intent(in) :: unit
        real(real64), intent(in) :: t
        type(particle3d), intent(in) :: p(:)
        integer :: i, n

        n = size(p)
        write(unit, *) t, (p(i)%p, i = 1, n)
    end subroutine write_output

end program main
