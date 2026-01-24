! module containing the barnes-hut algorithm types and procedures
module barnes_hut_module
    use geometry
    use particle
    implicit none

        private
        public :: cell_t, barnes_hut_node
        public :: build_tree, delete_tree
        public :: calculate_forces, calculate_force_recursive

        ! barnes-hut accuracy parameter
        real(kind=dp), parameter :: theta = 0.5_dp

        ! type to define the spatial boundaries of a cell
        type :: range_t
            type(point3d) :: min_p, max_p
        end type range_t

        ! wrapper for the recursive pointer structure
        type :: barnes_hut_node
            type(cell_t), pointer :: ptr => null()
        end type barnes_hut_node

        ! the main octree cell structure
        type :: cell_t
            type(range_t) :: range       ! spatial bounds
            type(point3d) :: com         ! center of mass
            real(kind=dp) :: total_mass  ! total mass in the cell
            integer :: particle_idx = 0  ! index of the particle 
            integer :: num_particles = 0 ! 0: empty, 1: leaf, >1: branch
            ! 8 children subcells (2x2x2)
            type(barnes_hut_node), dimension(2,2,2) :: subcells 
        end type cell_t

    contains


        ! subroutine: build_tree
        subroutine build_tree(root, particles)
            type(cell_t), pointer :: root
            type(particle3d), dimension(:), intent(in) :: particles
            
            type(point3d) :: min_bound, max_bound
            integer :: i

            ! 1. determine the bounding box of the system
            call get_boundaries(particles, min_bound, max_bound)

            ! 2. allocate and initialize root
            allocate(root)
            root%range%min_p = min_bound
            root%range%max_p = max_bound
            root%total_mass = 0.0_dp
            root%com = point3d(0.0_dp, 0.0_dp, 0.0_dp)
            root%num_particles = 0
            root%particle_idx = 0

            ! 3. insert all particles into the tree
            do i = 1, size(particles)
                call insert_particle(root, particles, i)
            end do

        end subroutine build_tree


        ! subroutine: insert_particle
        ! inserts a particle index into the tree
        recursive subroutine insert_particle(node, particles, p_idx)
            type(cell_t), pointer :: node
            type(particle3d), dimension(:), intent(in) :: particles
            integer, intent(in) :: p_idx
            
            integer :: old_idx
            type(point3d) :: p_pos

            p_pos = particles(p_idx)%p

            ! update node mass and center of mass
            ! mass * com = old_mass * old_com + new_mass * new_pos
            ! new_com = (mass * com) / new_total_mass
            if (node%total_mass > 0.0_dp) then
                node%com%x = (node%com%x * node%total_mass + p_pos%x * particles(p_idx)%m) / (node%total_mass + particles(p_idx)%m)
                node%com%y = (node%com%y * node%total_mass + p_pos%y * particles(p_idx)%m) / (node%total_mass + particles(p_idx)%m)
                node%com%z = (node%com%z * node%total_mass + p_pos%z * particles(p_idx)%m) / (node%total_mass + particles(p_idx)%m)
            else
                node%com = p_pos
            end if
            node%total_mass = node%total_mass + particles(p_idx)%m
            node%num_particles = node%num_particles + 1

            ! if node was empty, it becomes a leaf
            if (node%num_particles == 1) then
                node%particle_idx = p_idx
                return
            end if

            ! if node was a leaf, we must subdivide it
            if (node%num_particles == 2) then
                old_idx = node%particle_idx
                node%particle_idx = 0 
                
                ! push the old particle to a child
                call push_to_child(node, particles, old_idx)
            end if

            ! insert the new particle into the correct child
            call push_to_child(node, particles, p_idx)

        end subroutine insert_particle


        ! subroutine: push_to_child
        recursive subroutine push_to_child(node, particles, idx)
            type(cell_t), pointer :: node
            type(particle3d), dimension(:), intent(in) :: particles
            integer, intent(in) :: idx
            
            integer :: ix, iy, iz
            type(point3d) :: mid
            
            ! calculate midpoint
            mid%x = 0.5_dp * (node%range%min_p%x + node%range%max_p%x)
            mid%y = 0.5_dp * (node%range%min_p%y + node%range%max_p%y)
            mid%z = 0.5_dp * (node%range%min_p%z + node%range%max_p%z)

            ! determine octant
            ix = 1; iy = 1; iz = 1
            if (particles(idx)%p%x > mid%x) ix = 2
            if (particles(idx)%p%y > mid%y) iy = 2
            if (particles(idx)%p%z > mid%z) iz = 2

            ! create child if it doesn't exist
            if (.not. associated(node%subcells(ix,iy,iz)%ptr)) then
                allocate(node%subcells(ix,iy,iz)%ptr)
                call init_child_range(node%subcells(ix,iy,iz)%ptr, node%range, ix, iy, iz)
            end if

            call insert_particle(node%subcells(ix,iy,iz)%ptr, particles, idx)

        end subroutine push_to_child


        ! subroutine: init_child_range
        ! calculates the bounding box for a child cell
        subroutine init_child_range(child, parent_range, ix, iy, iz)
            type(cell_t), pointer :: child
            type(range_t), intent(in) :: parent_range
            integer, intent(in) :: ix, iy, iz
            
            type(point3d) :: mid

            mid%x = 0.5_dp * (parent_range%min_p%x + parent_range%max_p%x)
            mid%y = 0.5_dp * (parent_range%min_p%y + parent_range%max_p%y)
            mid%z = 0.5_dp * (parent_range%min_p%z + parent_range%max_p%z)

            if (ix == 1) then
                child%range%min_p%x = parent_range%min_p%x
                child%range%max_p%x = mid%x
            else
                child%range%min_p%x = mid%x
                child%range%max_p%x = parent_range%max_p%x
            end if
            
            if (iy == 1) then
                child%range%min_p%y = parent_range%min_p%y
                child%range%max_p%y = mid%y
            else
                child%range%min_p%y = mid%y
                child%range%max_p%y = parent_range%max_p%y
            end if
            
            if (iz == 1) then
                child%range%min_p%z = parent_range%min_p%z
                child%range%max_p%z = mid%z
            else
                child%range%min_p%z = mid%z
                child%range%max_p%z = parent_range%max_p%z
            end if
            
            ! init other properties
            child%total_mass = 0.0_dp
            child%num_particles = 0
            child%particle_idx = 0
            child%com = point3d(0.0_dp, 0.0_dp, 0.0_dp)

        end subroutine init_child_range


        ! subroutine: calculate_forces
        ! computes acceleration for all particles
        ! OMP parallelization
        subroutine calculate_forces(root, particles, accelerations)
        type(cell_t), pointer :: root
        type(particle3d), dimension(:), intent(in) :: particles
        type(vector3d), dimension(:), intent(out) :: accelerations
        
        integer :: i

        !$omp parallel do private(i) shared(root, particles, accelerations)
        do i = 1, size(particles)
            accelerations(i) = vector3d(0.0_dp, 0.0_dp, 0.0_dp)
            call calculate_force_recursive(root, particles(i), accelerations(i))
        end do
        !$omp end parallel do

        end subroutine calculate_forces

        ! subroutine: calculate_force_recursive
        ! barnes-hut criterion logic
        recursive subroutine calculate_force_recursive(node, p, acc)
            type(cell_t), pointer :: node
            type(particle3d), intent(in) :: p
            type(vector3d), intent(inout) :: acc
            
            real(kind=dp) :: d, s, ratio
            type(vector3d) :: r_vec
            integer :: i, j, k

            if (.not. associated(node)) return
            if (node%num_particles == 0) return

            ! vector from particle to node COM
            r_vec = node%com - p%p
            d = sqrt(r_vec%x**2 + r_vec%y**2 + r_vec%z**2)

            ! width of the cell (assuming cube-ish, take max dim)
            s = max(node%range%max_p%x - node%range%min_p%x, &
                    node%range%max_p%y - node%range%min_p%y)

            ! if node is a leaf (and not the particle itself)
            if (node%num_particles == 1) then
                if (d > epsilon(0.0_dp)) then ! avoid self-interaction
                    acc = acc + (node%total_mass / d**3) * r_vec
                end if
            else 
                ! checks if node is far enough
                if (d > epsilon(0.0_dp)) then
                    ratio = s / d
                    if (ratio < theta) then
                        ! use approximation
                        acc = acc + (node%total_mass / d**3) * r_vec
                    else
                        ! recurse into children
                        do i = 1, 2
                            do j = 1, 2
                                do k = 1, 2
                                    call calculate_force_recursive(node%subcells(i,j,k)%ptr, p, acc)
                                end do
                            end do
                        end do
                    end if
                end if
            end if

        end subroutine calculate_force_recursive


        ! subroutine: get_boundaries
        subroutine get_boundaries(particles, min_p, max_p)
            type(particle3d), dimension(:), intent(in) :: particles
            type(point3d), intent(out) :: min_p, max_p
            
            integer :: i
            
            ! init with first particle
            min_p = particles(1)%p
            max_p = particles(1)%p
            
            do i = 2, size(particles)
                if (particles(i)%p%x < min_p%x) min_p%x = particles(i)%p%x
                if (particles(i)%p%y < min_p%y) min_p%y = particles(i)%p%y
                if (particles(i)%p%z < min_p%z) min_p%z = particles(i)%p%z
                
                if (particles(i)%p%x > max_p%x) max_p%x = particles(i)%p%x
                if (particles(i)%p%y > max_p%y) max_p%y = particles(i)%p%y
                if (particles(i)%p%z > max_p%z) max_p%z = particles(i)%p%z
            end do
            
            ! add small buffer to avoid boundary issues
            min_p%x = min_p%x - 0.001_dp
            max_p%x = max_p%x + 0.001_dp

        end subroutine get_boundaries

        ! subroutine: delete_tree
        recursive subroutine delete_tree(node)
            type(cell_t), pointer :: node
            integer :: i, j, k
            
            if (.not. associated(node)) return
            
            do i = 1, 2
                do j = 1, 2
                    do k = 1, 2
                        call delete_tree(node%subcells(i,j,k)%ptr)
                    end do
                end do
            end do
            
            deallocate(node)
            nullify(node)

        end subroutine delete_tree

end module barnes_hut_module