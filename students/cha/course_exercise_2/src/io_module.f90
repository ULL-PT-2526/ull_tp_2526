module io_module

    use iso_fortran_env, only: real64
    use particle
    use geometry
    implicit none
    private

    ! ------ PUBLIC SUBROUTINES ------
    public :: io_open_output, io_close_output, io_write_simulation_state

    ! ------ MODULE-LEVEL VARIABLES ------
    integer, public, save :: output_unit = -1
    character(len=256), public, save :: output_filename = ''
    logical, public, save :: file_opened = .false.

contains

    ! ------ io_open_output ------
    ! Opens the output file and writes the simulation parameters at the top
    subroutine io_open_output(n_particles, dt, dt_out, t_end, disk_radius)
        integer, intent(in) :: n_particles
        real(real64), intent(in) :: dt, dt_out, t_end, disk_radius
        
        integer :: io_stat
        character(len=32) :: n_str

        ! Convert number of particles to string (for filename)
        write(n_str, '(I0)') n_particles - 1

        ! Construct the filename: "output_disk_" + N + ".dat"
        output_filename = 'output_disk_' // trim(n_str) // '.dat'

        ! Open file (always replace)
        open(newunit=output_unit, file=trim(output_filename), &
             status='replace', action='write', iostat=io_stat)

        ! Minimal error check
        if (io_stat /= 0) then
            print*, 'Error opening file: ', trim(output_filename)
            stop
        end if

        file_opened = .true.
        print*, 'Opened output file: ', trim(output_filename)

        ! --- WRITE HEADER ---
        write(output_unit, '(ES12.3)') dt
        write(output_unit, '(ES12.3)') dt_out
        write(output_unit, '(F10.3)')  t_end
        write(output_unit, '(F10.3)')  disk_radius
        write(output_unit, '(I0)')     n_particles

    end subroutine io_open_output

    ! ------ io_close_output ------
    subroutine io_close_output()
        if (file_opened) then
            close(output_unit)
            print*, 'Closed output file: ', trim(output_filename)
            file_opened = .false.
        end if
    end subroutine io_close_output

    ! ------ io_write_simulation_state ------
    ! Writes [Time] [x1] [y1] [z1] ... [xN] [yN] [zN]
    subroutine io_write_simulation_state(time, particles_array, n_particles)
        real(real64), intent(in) :: time
        type(particle3d), dimension(:), intent(in) :: particles_array
        integer, intent(in) :: n_particles
        integer :: i
        character(len=32) :: fmt_string

        if (.not. file_opened) return

        ! Build dynamic format string
        write(fmt_string, '(a,i0,a)') '(es15.8,', 3*n_particles, '(1x,es14.6))'

        write(output_unit, fmt=trim(fmt_string)) time, &
            (particles_array(i)%p%x, particles_array(i)%p%y, particles_array(i)%p%z, i=1, n_particles)

    end subroutine io_write_simulation_state

end module io_module






