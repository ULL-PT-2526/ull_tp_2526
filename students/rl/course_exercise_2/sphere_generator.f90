program generate_input
    implicit none

    integer :: i, n
    integer :: values(8), k
    integer, allocatable :: seed(:)
    real(8) :: mass, rx, ry, rz
    real(8) :: dt, dt_out, t_end

    ! -----------------------------
    ! Simulation parameters
    ! -----------------------------
    dt     = 1.0d-3
    dt_out = 1.0d-2
    t_end  = 1.0d1

    ! -----------------------------
    ! Number of particles
    ! -----------------------------
    read(*,*) n

    mass = 1.0d0 / n

    ! -----------------------------
    ! Random seed initialization
    ! -----------------------------
    call date_and_time(values=values)
    call random_seed(size=k)
    allocate(seed(k))
    seed = values(8)
    call random_seed(put=seed)

    ! -----------------------------
    ! Write header (input.dat format)
    ! -----------------------------
    write(*,'(ES20.10)') dt
    write(*,'(ES20.10)') dt_out
    write(*,'(ES20.10)') t_end
    write(*,'(I10)') n

    ! -----------------------------
    ! Generate particles in sphere
    ! -----------------------------
    do i = 1, n

        call random_number(rx)
        call random_number(ry)
        call random_number(rz)

        ! Rejection sampling to stay inside unit sphere
        do while (rx*rx + ry*ry + rz*rz > 1.0d0)
            call random_number(rx)
            call random_number(ry)
            call random_number(rz)
        end do

        write(*,'(ES20.10,3ES20.10,3ES20.10)') &
            mass, rx, ry, rz, 0.0d0, 0.0d0, 0.0d0
    end do

end program generate_input