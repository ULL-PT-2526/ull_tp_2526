module slist

  implicit none

  !=============================================================
  ! TYPE DEFINITIONS
  !=============================================================
  !--------- CELL ---------
  ! A single node (cell) in the list
  ! Each node stores a string and a pointer to the next node
  type :: Cell
    character(len=:), allocatable :: data  ! String stored in the cell
    type(Cell), pointer :: next => null()  ! Pointer to the next cell
  end type Cell

  !--------- SORTEDLIST ---------
  ! A sorted singly linked list.
  ! Contains a pointer to the first cell and the current length.
  type :: SortedList
    integer :: length = 0                  ! Number of elements
    type(Cell), pointer :: head => null()  ! Pointer to first cell
  end type SortedList

contains

  !=============================================================
  ! INIT_CELL:: creates and initialize a new Cell
  !=============================================================
  subroutine init_cell(new_cell, str)
    !! Allocates and initializes a new cell with a given string
    type(Cell), pointer, intent(out) :: new_cell
    character(len=*), intent(in) :: str
    integer :: istat

    allocate(new_cell, stat=istat)
    if (istat /= 0) stop 'Error: unable to allocate memory for Cell'

    if (len_trim(str) > 0) then
      allocate(character(len=len_trim(str)) :: new_cell%data, stat=istat)
      if (istat /= 0) stop 'Error: unable to allocate memory for cell data'
      new_cell%data = trim(str)
    else
      allocate(character(len=1) :: new_cell%data)
      new_cell%data = ' '
    end if

    new_cell%next => null()
  end subroutine init_cell

  !=============================================================
  ! EXTEND: adds a new string to the sorted list
  !=============================================================
  subroutine extend(list, str)
    ! Inserts a string into the list while maintaining ascending order
    type(SortedList), intent(inout) :: list
    character(*), intent(in) :: str
    type(Cell), pointer :: new_cell, current

    call init_cell(new_cell, str)

    ! Case 1: the list is empty
    if (.not. associated(list%head)) then
        list%head => new_cell
    else
        ! Case 2: the list is not empty. Check if new element should be the head
        if (new_cell%data < list%head%data) then
            new_cell%next => list%head
            list%head => new_cell
        ! Case 3: find the correct insertion spot within the list
        else
            current => list%head
            do while (associated(current%next))
                if (new_cell%data > current%next%data) then
                    current => current%next
                else
                    exit ! Found the spot
                end if
            end do
            
            ! Insert the cell
            new_cell%next => current%next
            current%next => new_cell
        end if
    end if

    list%length = list%length + 1
  end subroutine extend


  !=============================================================
  ! POP: removes and return the last element (tail) of the list
  !=============================================================
  function pop(list) result(popped_cell)
    ! Removes the last element in the list and returns a pointer to it
    ! If the list is empty, returns a null pointer
    type(SortedList), intent(inout) :: list
    type(Cell), pointer :: popped_cell, current, previous

    popped_cell => null()
    if (.not. associated(list%head)) return  ! Empty list

    previous => null()
    current => list%head

    ! Traverse to the last node
    do while (associated(current%next))
      previous => current
      current => current%next
    end do

    ! Detach the last node
    popped_cell => current
    if (associated(previous)) then
      previous%next => null()
    else
      list%head => null()  ! List had only one element
    end if

    list%length = list%length - 1
  end function pop

  !=============================================================
  ! PRINT_FORWARD: prints the list in forward order (head → tail)
  !=============================================================
  recursive subroutine print_forward(node)
    ! Recursively prints each node from head to tail
    type(Cell), pointer, intent(in) :: node
    if (.not. associated(node)) return
    print '("  ", A)', trim(node%data)
    if (associated(node%next)) call print_forward(node%next)
  end subroutine print_forward

  !=============================================================
  ! PRINT_REVERSE: prints the list in reverse order (tail → head)
  !=============================================================
  recursive subroutine print_reverse(node)
    ! Recursively prints each node from tail to head
    type(Cell), pointer, intent(in) :: node
    if (.not. associated(node)) return
    if (associated(node%next)) call print_reverse(node%next)
    print '("  ", A)', trim(node%data)
  end subroutine print_reverse

  !=============================================================
  ! EMPTY_LIST: delete all elements and release memory
  ! This is a personal implementation, I think is important in order
  ! to realease memory when the process is finished
  !=============================================================
  subroutine empty_list(list)
    ! Frees all nodes in the list and resets it to empty state
    type(SortedList), intent(inout) :: list
    type(Cell), pointer :: current, next_cell

    current => list%head
    do while (associated(current))
      next_cell => current%next
      if (allocated(current%data)) deallocate(current%data)
      deallocate(current)
      current => next_cell
    end do

    list%head => null()
    list%length = 0
  end subroutine empty_list

end module slist

