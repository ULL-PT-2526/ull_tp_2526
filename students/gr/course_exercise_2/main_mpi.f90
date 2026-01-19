program main
    use MPI
    use iso_fortran_env, only: real64
    use barnes_hut_mpi
    use geometry
    use particle
    implicit none

    integer :: my_rank, n_procs, error
    integer :: my_n, my_start, my_end, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    character(len=300) :: setup_file_name

    type(particle3d), allocatable :: parts(:)
    type(vector3d),  allocatable :: a(:)
    type(cell), pointer :: head

    integer :: out_unit = 20

    !---- TIMING ----
    real(real64) :: time_start, time_end, elapsed

    !=======================
    ! MPI INIT
    !=======================
    call MPI_INIT(error)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, error)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, error)

    !=======================
    ! Read setup (only rank 0)
    !=======================
    if (my_rank == 0) then
        print *, 'Insert name of the simulation setup file:'
        read *, setup_file_name
    end if

    call read_setup_mpi(setup_file_name, dt, dt_out, t_end, n, parts, my_rank)

    allocate(a(n))

    my_n     = n / n_procs
    my_start = (my_n * my_rank) + 1
    my_end   = my_start + my_n - 1

    !=======================
    ! BH INIT
    !=======================
    allocate(head)
    call initialize_bh_mpi(head, parts, a, my_start, my_end)

    !=======================
    ! Output file
    !=======================
    if (my_rank == 0) then
        block
            character(len=400) :: out_name
            ! Generate name: outputs/output_mpi_nVALUE.dat
            write(out_name, '(A, I0, A)') &
                "outputs/output_mpi_n", n, ".dat"
            
            call execute_command_line("mkdir -p outputs")
            open(unit=out_unit, file=trim(out_name), status='replace', action='write')
            print *, 'MPI Output file created: ', trim(out_name)
        end block
    end if

    !======================
    ! TIMING START
    !======================
    time_start = MPI_WTIME()

    !=======================
    ! Main integration loop
    !=======================
    t = 0.0_real64
    t_out = 0.0_real64

    do while (t <= t_end)

        call update_positions_mpi(parts, a, dt, my_start, my_end)

        call sync_positions(parts, my_start, my_end, n)

        call rebuild_tree_and_forces_mpi(head, parts, a, my_start, my_end)

        call update_velocities_mpi(parts, a, dt, my_start, my_end)

        t_out = t_out + dt

        if (t_out >= dt_out) then
            if (my_rank == 0) then
                call write_output(out_unit, t, parts)
            end if
            t_out = 0.0_real64
        end if

        t = t + dt
    end do

    !======================
    ! TIMING END
    !======================
    time_end = MPI_WTIME()
    elapsed = time_end - time_start

    if (my_rank == 0) then
        close(out_unit)

        open(unit=99, file='timing.dat', position='append', action='write')

        write(99,'(A,1X,I8,1X,I4,1X,F12.6)') &
             "mpi", n, n_procs, elapsed

        close(99)

        print *, 'MPI Simulation time:', elapsed, 'seconds'
        print *, 'Results saved in "outputs" folder'
        print *, 'Timing saved in "timing.dat"'
    end if

    call MPI_Finalize(error)

contains

    !====================================================
    ! Subroutine: read_setup_mpi
    !====================================================
    subroutine read_setup_mpi(filename, dt, dt_out, t_end, n, parts, rank)
        character(len=*), intent(in) :: filename
        real(real64), intent(out) :: dt, dt_out, t_end
        integer, intent(out) :: n
        type(particle3d), allocatable, intent(out) :: parts(:)
        integer, intent(in) :: rank
        
        real(real64), allocatable :: buffer_all(:,:)
        integer :: i, error

        if (rank == 0) then
            open(10,file=filename,status='old',action='read')
            read(10,*) dt
            read(10,*) dt_out
            read(10,*) t_end
            read(10,*) n
        end if

        call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
        call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
        call MPI_BCAST(dt_out,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
        call MPI_BCAST(t_end,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)

        allocate(parts(n))
        ! 7 for each particle: m, x, y, z, vx, vy, vz
        allocate(buffer_all(7, n)) 

        if (rank == 0) then
            do i=1,n
                read(10,*) parts(i)%m, parts(i)%p, parts(i)%v
                buffer_all(1,i) = parts(i)%m
                buffer_all(2,i) = parts(i)%p%x
                buffer_all(3,i) = parts(i)%p%y
                buffer_all(4,i) = parts(i)%p%z
                buffer_all(5,i) = parts(i)%v%x
                buffer_all(6,i) = parts(i)%v%y
                buffer_all(7,i) = parts(i)%v%z
            end do
            close(10)
        end if

        call MPI_BCAST(buffer_all, 7*n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)

        do i=1,n
            parts(i)%m   = buffer_all(1,i)
            parts(i)%p%x = buffer_all(2,i)
            parts(i)%p%y = buffer_all(3,i)
            parts(i)%p%z = buffer_all(4,i)
            parts(i)%v%x = buffer_all(5,i)
            parts(i)%v%y = buffer_all(6,i)
            parts(i)%v%z = buffer_all(7,i)
        end do
        
        deallocate(buffer_all)

    end subroutine read_setup_mpi

    !====================================================
    ! Subroutine: initialize_bh_mpi
    !====================================================
    subroutine initialize_bh_mpi(head, p, a, s, e)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)
        integer, intent(in) :: s,e

        call bh_init(head, p)
        call bh_build_tree(head, p)
        call bh_compute_forces(head, s, e, p, a)
    end subroutine

    !====================================================
    ! Subroutine: update_positions_mpi
    !====================================================
    subroutine update_positions_mpi(p, a, dt, s, e)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt
        integer, intent(in) :: s,e

        p(s:e)%v = p(s:e)%v + a(s:e)*(dt*0.5_real64)
        p(s:e)%p = p(s:e)%p + p(s:e)%v*dt
    end subroutine

    !====================================================
    ! Subroutine: sync_positions
    !====================================================
    subroutine sync_positions(parts, s, e, n)
        type(particle3d), intent(inout) :: parts(:)
        integer, intent(in) :: s, e, n
        integer :: error
        
        real(real64), allocatable :: send_buf(:,:), recv_buf(:,:)
        integer :: my_n, i, k

        my_n = e - s + 1
        
        allocate(send_buf(3, my_n))
        allocate(recv_buf(3, n))
        
        do k = 1, my_n
            i = s + k - 1
            send_buf(1, k) = parts(i)%p%x
            send_buf(2, k) = parts(i)%p%y
            send_buf(3, k) = parts(i)%p%z
        end do
        
        call MPI_ALLGATHER( &
            send_buf, 3*my_n, MPI_DOUBLE_PRECISION, &
            recv_buf, 3*my_n, MPI_DOUBLE_PRECISION, &
            MPI_COMM_WORLD, error)

        do i = 1, n
            parts(i)%p%x = recv_buf(1, i)
            parts(i)%p%y = recv_buf(2, i)
            parts(i)%p%z = recv_buf(3, i)
        end do

        deallocate(send_buf)
        deallocate(recv_buf)

    end subroutine sync_positions

    !====================================================
    ! Subroutine: rebuild_tree_and_forces_mpi
    !====================================================
    subroutine rebuild_tree_and_forces_mpi(head,p,a,s,e)
        type(cell), pointer :: head
        type(particle3d), intent(in) :: p(:)
        type(vector3d), intent(out) :: a(:)
        integer, intent(in) :: s,e

        call bh_destroy(head, p)
        call bh_build_tree(head, p)
        call bh_compute_forces(head, s, e, p, a)
    end subroutine

    !====================================================
    ! Subroutine: update_velocities_mpi
    !====================================================
    subroutine update_velocities_mpi(p,a,dt,s,e)
        type(particle3d), intent(inout) :: p(:)
        type(vector3d), intent(in) :: a(:)
        real(real64), intent(in) :: dt
        integer, intent(in) :: s,e

        p(s:e)%v = p(s:e)%v + a(s:e)*(dt*0.5_real64)
    end subroutine

    !====================================================
    ! Subroutine: write_output
    !====================================================
    subroutine write_output(unit, t, p)
        integer, intent(in) :: unit
        real(real64), intent(in) :: t
        type(particle3d), intent(in) :: p(:)
        integer :: i,n

        n = size(p)
        write(unit,*) t, (p(i)%p, i=1,n)
    end subroutine

end program main




