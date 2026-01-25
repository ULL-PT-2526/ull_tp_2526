module bh_mpi_module
  use geometry, only: dp
  implicit none
  private

  public :: range, cptr, cell
  public :: bhm_compute_ranges, bhm_find_cell, bhm_place_cell, bhm_create_subcells
  public :: bhm_nullify_pointers, bhm_belongs, bhm_compute_range
  public :: bhm_remove_empty_leaves, bhm_remove_tree, bhm_compute_masses
  public :: bhm_compute_forces_range
  public :: bhm_pool_init, bhm_pool_reset

  type :: range
    real(kind=dp), dimension(3) :: min, max
  end type range

  type :: cptr
    type(cell), pointer :: ptr
  end type cptr

  type :: cell
    type(range) :: range
    real(kind=dp), dimension(3) :: part
    integer :: pos
    integer :: type !! 0 = empty; 1 = particle; 2 = aggregate
    real(kind=dp) :: mass
    real(kind=dp), dimension(3) :: c_o_m
    type(cptr), dimension(2,2,2) :: subcell
  end type cell

  ! Preallocated pool for cells to avoid frequent allocate/deallocate
  type(cell), target, allocatable, save :: bhm_node_pool(:)
  integer, save :: bhm_pool_size = 0
  integer, save :: bhm_node_count = 0

contains

  subroutine bhm_pool_init(max_nodes)
    integer, intent(in) :: max_nodes
  integer :: i
  type(cell), pointer :: p
    if (allocated(bhm_node_pool)) deallocate(bhm_node_pool)
    allocate(bhm_node_pool(max_nodes))
    bhm_pool_size = max_nodes
    bhm_node_count = 0
    ! Ensure all child pointers are nullified initially
    do i = 1, max_nodes
      p => bhm_node_pool(i)
      call bhm_nullify_pointers(p)
      bhm_node_pool(i)%type = 0
      bhm_node_pool(i)%pos = 0
      bhm_node_pool(i)%mass = 0.0_dp
      bhm_node_pool(i)%part = 0.0_dp
      bhm_node_pool(i)%c_o_m = 0.0_dp
      bhm_node_pool(i)%range%min = 0.0_dp
      bhm_node_pool(i)%range%max = 0.0_dp
    end do
  end subroutine bhm_pool_init

  subroutine bhm_init_cell(node)
    type(cell), pointer :: node
    node%type = 0
    node%pos = 0
    node%mass = 0.0_dp
    node%part = 0.0_dp
    node%c_o_m = 0.0_dp
    node%range%min = 0.0_dp
    node%range%max = 0.0_dp
    call bhm_nullify_pointers(node)
  end subroutine bhm_init_cell

  subroutine bhm_pool_reset(root)
    type(cell), pointer :: root
    ! Start from first node in the pool as root
    if (.not. allocated(bhm_node_pool)) then
      stop 'bhm_pool_reset called before bhm_pool_init'
    end if
    bhm_node_count = 1
    root => bhm_node_pool(1)
    call bhm_init_cell(root)
  end subroutine bhm_pool_reset

  subroutine bhm_new_node(node)
    type(cell), pointer :: node
    if (bhm_node_count >= bhm_pool_size) then
      stop 'Barnes-Hut node pool exhausted. Increase pool size (bhm_pool_init)'
    end if
    bhm_node_count = bhm_node_count + 1
    node => bhm_node_pool(bhm_node_count)
    call bhm_init_cell(node)
  end subroutine bhm_new_node

  subroutine bhm_compute_ranges(goal, r)
    type(cell), pointer :: goal
    real(kind=dp), intent(in) :: r(:,:)  ! (n,3)
    real(kind=dp), dimension(3) :: mins, maxs, medios
    real(kind=dp) :: span
    mins(1) = minval(r(:,1)); mins(2) = minval(r(:,2)); mins(3) = minval(r(:,3))
    maxs(1) = maxval(r(:,1)); maxs(2) = maxval(r(:,2)); maxs(3) = maxval(r(:,3))
    span = maxval(maxs - mins) * 1.1_dp
    medios = (maxs + mins) / 2.0_dp
    goal%range%min = medios - span/2.0_dp
    goal%range%max = medios + span/2.0_dp
  end subroutine bhm_compute_ranges

  recursive subroutine bhm_find_cell(root, goal, part)
    type(cell), pointer :: root, goal, temp
    real(kind=dp), intent(in) :: part(3)
    integer :: i,j,k
    select case (root%type)
    case (2)
      out: do i = 1,2
        do j = 1,2
          do k = 1,2
            if (associated(root%subcell(i,j,k)%ptr)) then
              if (bhm_belongs(part, root%subcell(i,j,k)%ptr)) then
                call bhm_find_cell(root%subcell(i,j,k)%ptr, temp, part)
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
  end subroutine bhm_find_cell

  recursive subroutine bhm_place_cell(goal, part, n)
    type(cell), pointer :: goal, temp
    real(kind=dp), intent(in) :: part(3)
    integer, intent(in) :: n
    select case (goal%type)
    case (0)
      goal%type = 1
      goal%part = part
      goal%pos = n
    case (1)
      call bhm_create_subcells(goal)
      call bhm_find_cell(goal, temp, part)
      call bhm_place_cell(temp, part, n)
    case default
      print*, 'ERROR: invalid cell state'
    end select
  end subroutine bhm_place_cell

  subroutine bhm_create_subcells(goal)
    type(cell), pointer :: goal
    real(kind=dp), dimension(3) :: part
    integer :: i,j,k
    logical :: placed
    integer, dimension(3) :: octant
    part = goal%part
    goal%type = 2
    placed = .false.
    do i = 1,2
      do j = 1,2
        do k = 1,2
          octant = (/ i, j, k /)
          call bhm_new_node(goal%subcell(i,j,k)%ptr)
          goal%subcell(i,j,k)%ptr%range%min = bhm_compute_range(0, goal, octant)
          goal%subcell(i,j,k)%ptr%range%max = bhm_compute_range(1, goal, octant)
          ! Assign the existing particle to exactly one child to avoid duplication
          if ((.not. placed) .and. bhm_belongs(part, goal%subcell(i,j,k)%ptr)) then
            goal%subcell(i,j,k)%ptr%part = part
            goal%subcell(i,j,k)%ptr%type = 1
            goal%subcell(i,j,k)%ptr%pos = goal%pos
            placed = .true.
          else
            goal%subcell(i,j,k)%ptr%type = 0
          end if
          call bhm_nullify_pointers(goal%subcell(i,j,k)%ptr)
        end do
      end do
    end do
  end subroutine bhm_create_subcells

  subroutine bhm_nullify_pointers(goal)
    type(cell), pointer :: goal
    integer :: i,j,k
    do i = 1,2
      do j = 1,2
        do k = 1,2
          nullify(goal%subcell(i,j,k)%ptr)
        end do
      end do
    end do
  end subroutine bhm_nullify_pointers

  logical function bhm_belongs(part, goal)
    real(kind=dp), intent(in) :: part(3)
    type(cell), pointer :: goal
    bhm_belongs = (part(1) >= goal%range%min(1) .and. part(1) <= goal%range%max(1) .and. &
                   part(2) >= goal%range%min(2) .and. part(2) <= goal%range%max(2) .and. &
                   part(3) >= goal%range%min(3) .and. part(3) <= goal%range%max(3))
  end function bhm_belongs

  function bhm_compute_range(what, goal, octant) result(r)
    integer, intent(in) :: what
    type(cell), pointer :: goal
    integer, dimension(3), intent(in) :: octant
    real(kind=dp), dimension(3) :: r, valor_medio
    valor_medio = (goal%range%min + goal%range%max) / 2.0_dp
    select case (what)
    case (0)
      where (octant == 1)
        r = goal%range%min
      elsewhere
        r = valor_medio
      endwhere
    case (1)
      where (octant == 1)
        r = valor_medio
      elsewhere
        r = goal%range%max
      endwhere
    end select
  end function bhm_compute_range

  recursive subroutine bhm_remove_empty_leaves(goal)
    type(cell), pointer :: goal
    integer :: i,j,k
    if (associated(goal%subcell(1,1,1)%ptr)) then
      do i = 1,2
        do j = 1,2
          do k = 1,2
            call bhm_remove_empty_leaves(goal%subcell(i,j,k)%ptr)
            if (goal%subcell(i,j,k)%ptr%type == 0) then
              nullify(goal%subcell(i,j,k)%ptr)
            end if
          end do
        end do
      end do
    end if
  end subroutine bhm_remove_empty_leaves

  recursive subroutine bhm_remove_tree(goal)
    type(cell), pointer :: goal
    integer :: i,j,k
    do i = 1,2
      do j = 1,2
        do k = 1,2
          if (associated(goal%subcell(i,j,k)%ptr)) then
            call bhm_remove_tree(goal%subcell(i,j,k)%ptr)
            nullify(goal%subcell(i,j,k)%ptr)
          end if
        end do
      end do
    end do
  end subroutine bhm_remove_tree

  recursive subroutine bhm_compute_masses(goal, m, r)
    type(cell), pointer :: goal
    real(kind=dp), intent(in) :: m(:)
    real(kind=dp), intent(in) :: r(:,:)  ! (n,3)
    integer :: i,j,k
    real(kind=dp) :: mass
    goal%mass = 0.0_dp
    goal%c_o_m = 0.0_dp
    select case (goal%type)
    case (1)
      goal%mass = m(goal%pos)
      goal%c_o_m(1) = r(goal%pos,1)
      goal%c_o_m(2) = r(goal%pos,2)
      goal%c_o_m(3) = r(goal%pos,3)
    case (2)
      do i = 1,2
        do j = 1,2
          do k = 1,2
            if (associated(goal%subcell(i,j,k)%ptr)) then
              call bhm_compute_masses(goal%subcell(i,j,k)%ptr, m, r)
              mass = goal%mass
              goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
              if (goal%mass > 0.0_dp) then
                goal%c_o_m = (mass * goal%c_o_m + goal%subcell(i,j,k)%ptr%mass * &
                               goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
              end if
            end if
          end do
        end do
      end do
    end select
  end subroutine bhm_compute_masses

  subroutine bhm_compute_forces_range(head, m, r, a, theta, i_start, i_end)
    type(cell), pointer :: head
    real(kind=dp), intent(in) :: m(:)
    real(kind=dp), intent(in) :: r(:,:)  ! (n,3)
    real(kind=dp), intent(inout) :: a(:,:) ! (n,3)
    real(kind=dp), intent(in) :: theta
    integer, intent(in) :: i_start, i_end
    integer :: i
    do i = i_start, i_end
      call bhm_compute_forces_aux(i, head, m, r, a, theta)
    end do
  end subroutine bhm_compute_forces_range

  recursive subroutine bhm_compute_forces_aux(goal, tree, m, r, a, theta)
    integer, intent(in) :: goal
    type(cell), pointer :: tree
    real(kind=dp), intent(in) :: m(:)
    real(kind=dp), intent(in) :: r(:,:)  ! (n,3)
    real(kind=dp), intent(inout) :: a(:,:) ! (n,3)
    real(kind=dp), intent(in) :: theta
    integer :: i,j,k
    real(kind=dp) :: r2, r3, D, l
    real(kind=dp), dimension(3) :: rji
    select case (tree%type)
    case (1)
      if (goal /= tree%pos) then
        rji = tree%c_o_m - r(goal,:)
        r2 = rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3)
        r3 = r2 * sqrt(r2)
        a(goal,:) = a(goal,:) + m(tree%pos) * rji / r3
      end if
    case (2)
      l = tree%range%max(1) - tree%range%min(1)
      rji = tree%c_o_m - r(goal,:)
      r2 = rji(1)*rji(1) + rji(2)*rji(2) + rji(3)*rji(3)
      D = sqrt(r2)
      if (l / D < theta) then
        r3 = r2 * D
        a(goal,:) = a(goal,:) + tree%mass * rji / r3
      else
        do i = 1,2
          do j = 1,2
            do k = 1,2
              if (associated(tree%subcell(i,j,k)%ptr)) then
                call bhm_compute_forces_aux(goal, tree%subcell(i,j,k)%ptr, m, r, a, theta)
              end if
            end do
          end do
        end do
      end if
    end select
  end subroutine bhm_compute_forces_aux

end module bh_mpi_module
