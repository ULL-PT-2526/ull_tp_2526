program main_barnes_hut
	use geometry
	use particle
	use barnes_hut
	implicit none

		integer :: k, n, nstep
		real(kind=dp) :: dt, t_end, dt_out, t_out, t
		type(particle3d), allocatable :: particles(:)
		type(vector3d), allocatable :: a(:)
		integer :: outunit
	
		! open output file
		open(newunit=outunit, file='output.dat', status='replace', action='write')

		! read input
		call read_input(dt, dt_out, t_end, n, particles)

		! allocate accelerations
		allocate(a(n))

		! compute initial accelerations and half-step velocities
		call compute_accelerations(particles, n, a)
		call update_velocities(particles, n, a, dt/2.0_dp)

		! init sim variables
		nstep = int(t_end/dt)
		t_out = 0.0_dp
		t = 0.0_dp

		! write initial state
		call write_output(outunit, t, particles, n)

		! main simulation loop
		do k = 1, nstep
			! leapfrog step
			call update_positions(particles, n, dt)
			call compute_accelerations(particles, n, a)
			call update_velocities(particles, n, a, dt)

			! update time and periodic write
			t = t + dt
			t_out = t_out + dt
			if (t_out >= dt_out .or. k == nstep) then
				call write_output(outunit, t, particles, n)
				t_out = 0.0_dp
			end if
		end do

		close(outunit)

		! cleanup
		deallocate(particles)
		deallocate(a)

	contains

		subroutine read_input(dt, dt_out, t_end, n, particles)
			implicit none
			real(kind=dp), intent(out) :: dt, dt_out, t_end
			integer, intent(out) :: n
			type(particle3d), allocatable, intent(out) :: particles(:)
			integer :: i
			read*, dt
			read*, dt_out
			read*, t_end
			read*, n
			allocate(particles(n))
			do i = 1, n
				read*, particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
						particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
			end do
		end subroutine read_input

		subroutine compute_accelerations(particles, n, accelerations)
			implicit none
			type(particle3d), intent(in) :: particles(:)
			integer, intent(in) :: n
			type(vector3d), intent(out) :: accelerations(:)
			real(kind=dp), parameter :: theta = 1.0_dp
			integer :: i
			type(cell), pointer :: head, temp
			! build tree and compute accelerations via Barnes–Hut
			allocate(head)
			call bh_calculate_ranges(head, particles)
			head%type = 0
			call bh_nullify_pointers(head)
			! insert all particles
			do i = 1, n
				call bh_find_cell(head, temp, particles(i)%p)
				call bh_place_cell(temp, particles(i)%p, i)
			end do
			call bh_borrar_empty_leaves(head)
			call bh_calculate_masses(head, particles)
			! zero accelerations
			!$omp parallel do default(shared) private(i) schedule(static)
			do i = 1, n
				accelerations(i)%x = 0.0_dp
				accelerations(i)%y = 0.0_dp
				accelerations(i)%z = 0.0_dp
			end do
			!$omp end parallel do
			call bh_calculate_forces(head, particles, accelerations, theta)
			! teardown tree
			call bh_borrar_tree(head)
			deallocate(head)
		end subroutine compute_accelerations

		subroutine update_positions(particles, n, dt)
			implicit none
			type(particle3d), intent(inout) :: particles(:)
			integer, intent(in) :: n
			real(kind=dp), intent(in) :: dt
			integer :: i
			!$omp parallel do default(shared) private(i) schedule(static)
			do i = 1, n
				particles(i)%p = particles(i)%p + particles(i)%v * dt
			end do
			!$omp end parallel do
		end subroutine update_positions

		subroutine update_velocities(particles, n, accelerations, dt)
			implicit none
			type(particle3d), intent(inout) :: particles(:)
			integer, intent(in) :: n
			type(vector3d), intent(in) :: accelerations(:)
			real(kind=dp), intent(in) :: dt
			integer :: i
			!$omp parallel do default(shared) private(i) schedule(static)
			do i = 1, n
				particles(i)%v = particles(i)%v + accelerations(i) * dt
			end do
			!$omp end parallel do
		end subroutine update_velocities

		subroutine write_output(outunit, t, particles, n)
			implicit none
			integer, intent(in) :: outunit
			real(kind=dp), intent(in) :: t
			type(particle3d), intent(in) :: particles(:)
			integer, intent(in) :: n
			integer :: i
			write(outunit, '(F10.5)', advance='no') t
			do i = 1, n
				write(outunit, '(3F15.8)', advance='no') particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
			end do
			write(outunit, *)
		end subroutine write_output

	end program main_barnes_hut
