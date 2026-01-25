module barnes_hut
	use geometry
	use particle
	implicit none
	private
	
	! public types and procedures
	public :: cell, cptr, range
	public :: bh_compute_ranges, bh_find_cell, bh_place_cell, bh_create_subcells
	public :: bh_nullify_pointers, bh_belongs, bh_compute_range
	public :: bh_remove_empty_leaves, bh_remove_tree, bh_compute_masses
	public :: bh_compute_forces
	! Legacy aliases expected by drivers
	public :: bh_calculate_ranges, bh_calculate_masses, bh_calculate_forces
	public :: bh_borrar_empty_leaves, bh_borrar_tree

	type :: range
		type(point3d) :: min, max
	end type range

	type :: cptr
		type(cell), pointer :: ptr
	end type cptr

	type :: cell
		type(range) :: range
		type(point3d) :: part
		integer :: pos
		integer :: type !! 0 = no particle; 1 = particle; 2 = aggregate
		real(kind=dp) :: mass
		type(point3d) :: c_o_m
		type(cptr), dimension(2,2,2) :: subcell
	end type cell

contains

	! =============================
	! bh_compute_ranges
	! =============================
	subroutine bh_compute_ranges(goal, particles)
		type(cell), pointer :: goal
		type(particle3d), intent(in) :: particles(:)
		type(point3d) :: mins, maxs, means
		real(kind=dp) :: span, dx, dy, dz
		integer :: i
		! initialize mins/maxs with first particle position
		mins = particles(1)%p
		maxs = particles(1)%p
		do i = 2, size(particles)
			if (particles(i)%p%x < mins%x) mins%x = particles(i)%p%x
			if (particles(i)%p%y < mins%y) mins%y = particles(i)%p%y
			if (particles(i)%p%z < mins%z) mins%z = particles(i)%p%z
			if (particles(i)%p%x > maxs%x) maxs%x = particles(i)%p%x
			if (particles(i)%p%y > maxs%y) maxs%y = particles(i)%p%y
			if (particles(i)%p%z > maxs%z) maxs%z = particles(i)%p%z
		end do
		dx = maxs%x - mins%x
		dy = maxs%y - mins%y
		dz = maxs%z - mins%z
		span = max(dx, max(dy, dz)) * 1.1_dp  ! add 10% padding
		means%x = (maxs%x + mins%x) / 2.0_dp
		means%y = (maxs%y + mins%y) / 2.0_dp
		means%z = (maxs%z + mins%z) / 2.0_dp
		goal%range%min%x = means%x - span/2.0_dp
		goal%range%min%y = means%y - span/2.0_dp
		goal%range%min%z = means%z - span/2.0_dp
		goal%range%max%x = means%x + span/2.0_dp
		goal%range%max%y = means%y + span/2.0_dp
		goal%range%max%z = means%z + span/2.0_dp
	end subroutine bh_compute_ranges

	! =============================
	! bh_find_cell
	! =============================
	recursive subroutine bh_find_cell(root, goal, part)
		type(cell), pointer :: root, goal, temp
		type(point3d), intent(in) :: part
		integer :: i, j, k
		select case (root%type)
		case (2)
			out: do i = 1, 2
				do j = 1, 2
					do k = 1, 2
						! only test belonging if the subcell pointer exists
						if (associated(root%subcell(i,j,k)%ptr)) then
							if (bh_belongs(part, root%subcell(i,j,k)%ptr)) then
								call bh_find_cell(root%subcell(i,j,k)%ptr, temp, part)
								goal => temp
								exit out
							end if
						end if
					end do
					end do
				end do out
		case default
			goal => root
		end select
	end subroutine bh_find_cell

	! =============================
	! bh_place_cell
	! =============================
	recursive subroutine bh_place_cell(goal, part, n)
		type(cell), pointer :: goal, temp
		type(point3d), intent(in) :: part
		integer, intent(in) :: n
		select case (goal%type)
		case (0)
			goal%type = 1
			goal%part = part
			goal%pos = n
		case (1)
			call bh_create_subcells(goal)
			call bh_find_cell(goal, temp, part)
			call bh_place_cell(temp, part, n)
		case default
			print*, 'SHOULD NOT BE HERE. ERROR!'
		end select
	end subroutine bh_place_cell

	! =============================
	! bh_create_subcells
	! =============================
	subroutine bh_create_subcells(goal)
		type(cell), pointer :: goal
		type(point3d) :: part
		integer :: i, j, k
		logical :: placed
		integer, dimension(3) :: octant
		part = goal%part
		goal%type = 2
		placed = .false.
		do i = 1, 2
			do j = 1, 2
				do k = 1, 2
					octant = (/ i, j, k /)
					allocate(goal%subcell(i,j,k)%ptr)
					goal%subcell(i,j,k)%ptr%range%min = bh_compute_range(0, goal, octant)
					goal%subcell(i,j,k)%ptr%range%max = bh_compute_range(1, goal, octant)
					! Assign the existing particle to exactly one child to avoid duplication
					if (.not. placed .and. bh_belongs(part, goal%subcell(i,j,k)%ptr)) then
						goal%subcell(i,j,k)%ptr%part = part
						goal%subcell(i,j,k)%ptr%type = 1
						goal%subcell(i,j,k)%ptr%pos = goal%pos
						placed = .true.
					else
						goal%subcell(i,j,k)%ptr%type = 0
					end if
					call bh_nullify_pointers(goal%subcell(i,j,k)%ptr)
				end do
			end do
		end do
	end subroutine bh_create_subcells

	! =============================
	! bh_nullify_pointers
	! =============================
	subroutine bh_nullify_pointers(goal)
		type(cell), pointer :: goal
		integer :: i, j, k
		do i = 1, 2
			do j = 1, 2
				do k = 1, 2
					nullify(goal%subcell(i,j,k)%ptr)
				end do
			end do
		end do
	end subroutine bh_nullify_pointers

	! =============================
	! bh_belongs
	! =============================
	logical function bh_belongs(part, goal)
		type(point3d), intent(in) :: part
		type(cell), pointer :: goal
		bh_belongs = &
			(part%x >= goal%range%min%x .and. part%x <= goal%range%max%x .and. &
			 part%y >= goal%range%min%y .and. part%y <= goal%range%max%y .and. &
			 part%z >= goal%range%min%z .and. part%z <= goal%range%max%z)
	end function bh_belongs

	! =============================
	! bh_compute_range
	! =============================
	function bh_compute_range(what, goal, octant) result(r)
		integer, intent(in) :: what
		type(cell), pointer :: goal
		integer, dimension(3), intent(in) :: octant
		type(point3d) :: r, mean_value
		mean_value%x = (goal%range%min%x + goal%range%max%x) / 2.0_dp
		mean_value%y = (goal%range%min%y + goal%range%max%y) / 2.0_dp
		mean_value%z = (goal%range%min%z + goal%range%max%z) / 2.0_dp
		select case (what)
		case (0)
			! mins
			r%x = merge(goal%range%min%x, mean_value%x, octant(1) == 1)
			r%y = merge(goal%range%min%y, mean_value%y, octant(2) == 1)
			r%z = merge(goal%range%min%z, mean_value%z, octant(3) == 1)
		case (1)
			! maxs
			r%x = merge(mean_value%x, goal%range%max%x, octant(1) == 1)
			r%y = merge(mean_value%y, goal%range%max%y, octant(2) == 1)
			r%z = merge(mean_value%z, goal%range%max%z, octant(3) == 1)
		end select
	end function bh_compute_range

		! Legacy wrappers (aliases) to keep driver code working
		subroutine bh_calculate_ranges(goal, particles)
			type(cell), pointer :: goal
			type(particle3d), intent(in) :: particles(:)
			call bh_compute_ranges(goal, particles)
		end subroutine bh_calculate_ranges

		subroutine bh_borrar_empty_leaves(goal)
			type(cell), pointer :: goal
			call bh_remove_empty_leaves(goal)
		end subroutine bh_borrar_empty_leaves

		subroutine bh_borrar_tree(goal)
			type(cell), pointer :: goal
			call bh_remove_tree(goal)
		end subroutine bh_borrar_tree

		subroutine bh_calculate_masses(goal, particles)
			type(cell), pointer :: goal
			type(particle3d), intent(in) :: particles(:)
			call bh_compute_masses(goal, particles)
		end subroutine bh_calculate_masses

		subroutine bh_calculate_forces(head, particles, a, theta)
			type(cell), pointer :: head
			type(particle3d), intent(in) :: particles(:)
			type(vector3d), intent(inout) :: a(:)
			real(kind=dp), intent(in) :: theta
			call bh_compute_forces(head, particles, a, theta)
		end subroutine bh_calculate_forces

		

	! =============================
	! bh_remove_empty_leaves
	! =============================
	recursive subroutine bh_remove_empty_leaves(goal)
		type(cell), pointer :: goal
		integer :: i, j, k
		if (associated(goal%subcell(1,1,1)%ptr)) then
			do i = 1, 2
				do j = 1, 2
					do k = 1, 2
						call bh_remove_empty_leaves(goal%subcell(i,j,k)%ptr)
						if (goal%subcell(i,j,k)%ptr%type == 0) then
							deallocate(goal%subcell(i,j,k)%ptr)
						end if
					end do
				end do
			end do
		end if
	end subroutine bh_remove_empty_leaves

	! =============================
	! bh_remove_tree
	! =============================
	recursive subroutine bh_remove_tree(goal)
		type(cell), pointer :: goal
		integer :: i, j, k
		do i = 1, 2
			do j = 1, 2
				do k = 1, 2
					if (associated(goal%subcell(i,j,k)%ptr)) then
						call bh_remove_tree(goal%subcell(i,j,k)%ptr)
						deallocate(goal%subcell(i,j,k)%ptr)
					end if
				end do
			end do
		end do
	end subroutine bh_remove_tree

	! =============================
	! bh_compute_masses
	! =============================
	recursive subroutine bh_compute_masses(goal, particles)
		type(cell), pointer :: goal
		type(particle3d), intent(in) :: particles(:)
		integer :: i, j, k
		real(kind=dp) :: mass
		goal%mass = 0.0_dp
		goal%c_o_m%x = 0.0_dp
		goal%c_o_m%y = 0.0_dp
		goal%c_o_m%z = 0.0_dp
		select case (goal%type)
		case (1)
			goal%mass = particles(goal%pos)%m
			goal%c_o_m = particles(goal%pos)%p
		case (2)
			do i = 1, 2
				do j = 1, 2
					do k = 1, 2
						if (associated(goal%subcell(i,j,k)%ptr)) then
							call bh_compute_masses(goal%subcell(i,j,k)%ptr, particles)
							mass = goal%mass
							goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
							if (goal%mass > 0.0_dp) then
								goal%c_o_m%x = (mass * goal%c_o_m%x + &
									goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%mass
								goal%c_o_m%y = (mass * goal%c_o_m%y + &
									goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%mass
								goal%c_o_m%z = (mass * goal%c_o_m%z + &
									goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%mass
							end if
						end if
					end do
				end do
			end do
		end select
	end subroutine bh_compute_masses

	! =============================
	! bh_compute_forces
	! =============================
	subroutine bh_compute_forces(head, particles, a, theta)
		type(cell), pointer :: head
		type(particle3d), intent(in) :: particles(:)
		type(vector3d), intent(inout) :: a(:)
		real(kind=dp), intent(in) :: theta
		integer :: i
		!$omp parallel do default(shared) private(i) schedule(static)
		do i = 1, size(particles)
			call bh_compute_forces_aux(i, head, particles, a, theta)
		end do
		!$omp end parallel do
	end subroutine bh_compute_forces



	! =============================
	! bh_compute_forces_aux
	! =============================
	recursive subroutine bh_compute_forces_aux(goal, tree, particles, a, theta)
		integer, intent(in) :: goal
		type(cell), pointer :: tree
		type(particle3d), intent(in) :: particles(:)
		type(vector3d), intent(inout) :: a(:)
		real(kind=dp), intent(in) :: theta
		integer :: i, j, k
		real(kind=dp) :: r2, r3, D, l
		type(vector3d) :: rji
		select case (tree%type)
		case (1)
			if (goal /= tree%pos) then
				rji = tree%c_o_m - particles(goal)%p
				r2 = rji%x*rji%x + rji%y*rji%y + rji%z*rji%z
				r3 = r2 * sqrt(r2)
				a(goal) = a(goal) + (tree%mass * (rji / r3))
			end if
		case (2)
			! cell side length (cube) from any dimension (use x)
			l = tree%range%max%x - tree%range%min%x
			rji = tree%c_o_m - particles(goal)%p
			r2 = rji%x*rji%x + rji%y*rji%y + rji%z*rji%z
			D = sqrt(r2)
			if (l / D < theta) then
				r3 = r2 * D
				a(goal) = a(goal) + (tree%mass * (rji / r3))
			else
				do i = 1, 2
					do j = 1, 2
						do k = 1, 2
							if (associated(tree%subcell(i,j,k)%ptr)) then
								call bh_compute_forces_aux(goal, tree%subcell(i,j,k)%ptr, particles, a, theta)
							end if
						end do
					end do
				end do
			end if
		end select
	end subroutine bh_compute_forces_aux

end module barnes_hut
