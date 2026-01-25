program main_barnes_hut_mpi
	use mpi
	use geometry, only: dp
	use bh_mpi_module
	implicit none

	integer :: my_rank, p, ierr
	integer :: n, my_n, my_start, my_end
	real(kind=dp) :: dt, t_end, dt_out, t_out, t
	real(kind=dp), allocatable :: m(:)
	real(kind=dp), allocatable :: r(:,:), v(:,:), a(:,:)
	integer :: outunit
	integer, allocatable :: counts(:), displs(:), displs0(:)
	integer :: i
	real(kind=dp), parameter :: theta = 1.0_dp
	type(cell), pointer :: head, temp

	! mpi init
	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world, p, ierr)
	call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

	! rank 0 reads inputs and broadcasts
	! rank 0 reads and broadcasts input
	if (my_rank == 0) then
		read*, dt
		read*, dt_out
		read*, t_end
		read*, n
	end if
	call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
	call mpi_bcast(dt, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
	call mpi_bcast(dt_out, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
	call mpi_bcast(t_end, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
	allocate(m(n), r(n,3), v(n,3), a(n,3))
	if (my_rank == 0) then
		do i = 1, n
			read*, m(i), r(i,1), r(i,2), r(i,3), v(i,1), v(i,2), v(i,3)
		end do
	end if
	call mpi_bcast(m, n, mpi_double_precision, 0, mpi_comm_world, ierr)
	call mpi_bcast(r, n*3, mpi_double_precision, 0, mpi_comm_world, ierr)
	call mpi_bcast(v, n*3, mpi_double_precision, 0, mpi_comm_world, ierr)
	call mpi_bcast(a, n*3, mpi_double_precision, 0, mpi_comm_world, ierr)

	! compute block decomposition
	allocate(counts(p), displs(p), displs0(p))
	do i = 1, p
		counts(i) = n / p
	end do
	do i = 1, mod(n,p)
		counts(i) = counts(i) + 1
	end do
	displs(1) = 1
	do i = 2, p
		displs(i) = displs(i-1) + counts(i-1)
	end do
	do i = 1, p
		displs0(i) = displs(i) - 1
	end do
	my_n = counts(my_rank+1)
	my_start = displs(my_rank+1)
	my_end = my_start + my_n - 1

	! Initialize Barnes–Hut node pool once (conservative size ~8*N)
	call bhm_pool_init(8*n + 8)

	! initial accelerations for my block and half-step velocities
	! build tree and initial accelerations on my block
	call bhm_pool_reset(head)
	call bhm_compute_ranges(head, r)
	head%type = 0
	call bhm_nullify_pointers(head)
	do i = 1, n
		call bhm_find_cell(head, temp, r(i,:))
		call bhm_place_cell(temp, r(i,:), i)
	end do
	call bhm_remove_empty_leaves(head)
	call bhm_compute_masses(head, m, r)
	a(my_start:my_end,:) = 0.0_dp
	call bhm_compute_forces_range(head, m, r, a, theta, my_start, my_end)
	! No dynamic allocation for head
	v(my_start:my_end,:) = v(my_start:my_end,:) + a(my_start:my_end,:) * (dt/2.0_dp)

	! initialize simulation variables
	t_out = 0.0_dp
	t = 0.0_dp

	! open output on master and write initial state
	if (my_rank == 0) then
		open(newunit=outunit, file='output.dat', status='replace', action='write')
		write(outunit, '(F10.5)', advance='no') t
		do i = 1, n
			write(outunit, '(3F15.8)', advance='no') r(i,1), r(i,2), r(i,3)
		end do
		write(outunit, *)
	end if

	! main loop
	do while (t < t_end)
		! update positions for my block
		r(my_start:my_end,:) = r(my_start:my_end,:) + v(my_start:my_end,:) * dt
		! gather positions from all ranks so everyone has full positions
	   call mpi_allgatherv(r(my_start,1), my_n, mpi_double_precision, r(1,1), counts, displs0, &
		   mpi_double_precision, mpi_comm_world, ierr)
	   call mpi_allgatherv(r(my_start,2), my_n, mpi_double_precision, r(1,2), counts, displs0, &
		   mpi_double_precision, mpi_comm_world, ierr)
	   call mpi_allgatherv(r(my_start,3), my_n, mpi_double_precision, r(1,3), counts, displs0, &
		   mpi_double_precision, mpi_comm_world, ierr)
		! compute accelerations for my block and full-step velocities
		call bhm_pool_reset(head)
		call bhm_compute_ranges(head, r)
		head%type = 0
		call bhm_nullify_pointers(head)
		do i = 1, n
			call bhm_find_cell(head, temp, r(i,:))
			call bhm_place_cell(temp, r(i,:), i)
		end do
		call bhm_remove_empty_leaves(head)
		call bhm_compute_masses(head, m, r)
		a(my_start:my_end,:) = 0.0_dp
		call bhm_compute_forces_range(head, m, r, a, theta, my_start, my_end)
		! Pool-based head
		v(my_start:my_end,:) = v(my_start:my_end,:) + a(my_start:my_end,:) * dt

		! update time and output
		t_out = t_out + dt
		if (my_rank == 0) then
			if (t_out >= dt_out .or. t + dt >= t_end) then
				write(outunit, '(F10.5)', advance='no') t
				do i = 1, n
					write(outunit, '(3F15.8)', advance='no') r(i,1), r(i,2), r(i,3)
				end do
				write(outunit, *)
				t_out = 0.0_dp
			end if
		end if
		t = t + dt
	end do

	if (my_rank == 0) close(outunit)
	call mpi_finalize(ierr)

end program main_barnes_hut_mpi
