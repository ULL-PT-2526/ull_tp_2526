program main
    use geometry
    use particle
    use barnes_hut
    !$ use omp_lib
    implicit none

    !---------------------------------------------------------
    ! Simulation parameters
    !---------------------------------------------------------
    real(dp) :: dt, dt_out, t_end, t, t_out
    real(dp) :: t_i, t_f 
    integer  :: i
    integer, parameter :: output_unit = 10

    t_i = omp_get_wtime() 

    !---------------------------------------------------------
    ! Read input (same order as original code)
    !---------------------------------------------------------
    read(*,*) dt
    read(*,*) dt_out
    read(*,*) t_end
    read(*,*) n

    allocate(parts(n))
    allocate(a(n))

    do i = 1, n
        read(*,*) parts(i)%m, &
                   parts(i)%p%x, parts(i)%p%y, parts(i)%p%z, &
                   parts(i)%v%x, parts(i)%v%y, parts(i)%v%z
    end do

    ! Open output file
    open(unit=output_unit, file="output.dat", status="replace", action="write")
    
    !---------------------------------------------------------
    ! Initialize tree
    !---------------------------------------------------------
    allocate(head)
    call calculate_ranges(head)
    head%type = 0
    call nullify_pointers(head)

    do i = 1, n
        call find_cell(head, temp_cell, parts(i)%p)
        call place_cell(temp_cell, parts(i)%p, i)
    end do

    call borrar_empty_leaves(head)
    call calculate_masses(head)

    !---------------------------------------------------------
    ! Initial accelerations
    !---------------------------------------------------------
    a = vector3d(0.0_dp, 0.0_dp, 0.0_dp) !check calculate_forces, this might be redundant
    call calculate_forces(head)

    !---------------------------------------------------------
    ! Main time integration loop (leapfrog)
    !---------------------------------------------------------
    t_out = 0.0_dp
    t     = 0.0_dp

    do while (t < t_end)

        do i = 1, n
            parts(i)%v = parts(i)%v + a(i) * (dt/2.0_dp)
        end do

        do i = 1, n
            parts(i)%p = parts(i)%p + parts(i)%v * dt
        end do

        ! Rebuild tree
        call borrar_tree(head)
        call calculate_ranges(head)
        head%type = 0
        call nullify_pointers(head)

        do i = 1, n
            call find_cell(head, temp_cell, parts(i)%p)
            call place_cell(temp_cell, parts(i)%p, i)
        end do

        call borrar_empty_leaves(head)
        call calculate_masses(head)

        ! Compute new aelerations
        a = vector3d(0.0_dp, 0.0_dp, 0.0_dp) !NOTE: check calculate_forces, this might be redundant
        call calculate_forces(head)

        do i = 1, n
            parts(i)%v = parts(i)%v + a(i) * (dt/2.0_dp)
        end do

        !-----------------------------------------------------
        ! Output 
        !-----------------------------------------------------
        t_out = t_out + dt
        if (t_out >= dt_out) then

            write(output_unit, '(ES20.10)', advance='no') t
            do i = 1, n
                write(output_unit, '(3ES20.10)', advance='no') &
                    parts(i)%p%x, parts(i)%p%y, parts(i)%p%z
            end do
            write(output_unit, *)   ! newline

            t_out = 0.0_dp
        end if

        t = t + dt
    end do

    !---------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------
    call borrar_tree(head)
    deallocate(head)
    deallocate(parts)
    deallocate(a)

    print*, "Simulation complete. Output written to output.dat"

    t_f = omp_get_wtime()

    print *, "Total execution time:", t_f - t_i

end program main


