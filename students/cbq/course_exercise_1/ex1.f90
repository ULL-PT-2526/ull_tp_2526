program leapfrog
    use geometry
    use particle

    implicit none

  ! We define the integers and vectors
    integer :: i, j, n, k, n_steps
    real(8) :: dt, t_end, t, dt_out, t_out
    real(8) :: r2, r3
    type(particle3d), dimension(:), allocatable :: particles
    type(vector3d) :: rji
    type(vector3d), dimension(:), allocatable :: accelerations
    integer :: u, uout, iarg, stat  !for the input and output files
    character(len=256) :: filename

    iarg = command_argument_count()
    if (iarg /= 1) then
        print *, 'Usage: ./ex1 <input_file>'
    stop
    end if
    
    call get_command_argument(1, filename)
    
    ! Open input file
    open(newunit=u, file=trim(filename), status='old', iostat=stat)

    ! Read input parameters
    read (u, *) dt
    read (u, *) dt_out  
    read (u, *) t_end
    read (u, *) n
    
    allocate(particles(n))
    allocate(accelerations(n))

    ! First read particle data
    do i = 1, n
        read(u, *) particles(i)%m, &
                  particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                  particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    end do
    
    close(u)

    accelerations = vector3d(0.0, 0.0, 0.0)


    do i= 1, n
        do j= i+1, n
        rji= particles(j)%p - particles(i)%p
        r2= rji%x**2 + rji%y**2 + rji%z **2
        r3= r2*sqrt(r2)
        accelerations(i) = accelerations(i) + (rji / r3) * particles(j)%m
        accelerations(j) = accelerations(j) - (rji / r3) * particles(i)%m
        end do
    end do

    !Output file
    open(newunit=uout, file='output.dat', status='replace', action='write', iostat=stat)
    
    
    n_steps = int(t_end / dt)
    t_out = 0.0_8
    t = 0.0_8

    do k = 0, n_steps
        t = k * dt
        do i = 1, n
            particles(i)%v = particles(i)%v + accelerations(i) * (dt/2.0_8)
        end do
        
        do i = 1, n
            particles(i)%p = particles(i)%p + particles(i)%v * dt
        end do
        
        do i = 1, n
            accelerations(i) = vector3d(0.0_8, 0.0_8, 0.0_8)
        end do

     do i = 1,n
        do j = i+1,n
            rji= particles(j)%p - particles(i)%p
            r2= rji%x**2 + rji%y**2 + rji%z **2
            r3= r2*sqrt(r2)
            accelerations(i) = accelerations(i) + (rji / r3) * particles(j)%m
            accelerations(j) = accelerations(j) - (rji / r3) * particles(i)%m
        end do
     end do

     
     do i = 1, n
            particles(i)%v = particles(i)%v + accelerations(i) * (dt/2.0)
        end do
     
     t_out = t_out + dt
        if (t_out >= dt_out) then
            write(uout, '(F12.6)', advance='no') t
            do i = 1, n
                write(uout, '(3F15.8)', advance='no') particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
            end do
            write(uout, *) 
            t_out = 0.0
        end if
        end do

    close(uout)

deallocate(particles)
deallocate(accelerations)

end program leapfrog