! module containing the linked list definition and procedures
! step 1 of the exercise
module lecture_10_ex_linkedlist_slist
    implicit none

    ! definition of the node type
    type :: cell
        character(len=:), allocatable :: data
        type(cell), pointer :: next => null()
    end type cell

    ! definition of the list container
    type :: sortedlist
        integer :: length = 0
        type(cell), pointer :: head => null()
    end type sortedlist

contains

    ! subroutine: extend
    ! adds a new string to the list, sorted by the first letter.
    subroutine extend(list, str)
        type(sortedlist), intent(inout) :: list
        character(len=*), intent(in) :: str
        
        type(cell), pointer :: new_node, current, previous

        ! allocate memory for the new node
        allocate(new_node)
        new_node%data = str
        new_node%next => null()

        ! list is empty or new node goes before head
        if (.not. associated(list%head)) then
            list%head => new_node
        else if (str(1:1) < list%head%data(1:1)) then
            new_node%next => list%head
            list%head => new_node
        else
            ! traverse to find insertion point
            current => list%head
            previous => null()
            
            do while (associated(current))
                if (str(1:1) < current%data(1:1)) exit
                previous => current
                current => current%next
            end do
            
            ! insert new node between previous and current
            new_node%next => current
            if (associated(previous)) then
                previous%next => new_node
            end if
        end if

        ! increment list length
        list%length = list%length + 1

    end subroutine extend


    ! subroutine: empty_list
    ! deallocates all memory used by the list
    subroutine empty_list(list)
        type(sortedlist), intent(inout) :: list
        type(cell), pointer :: current, temp

        current => list%head
        
        do while (associated(current))
            temp => current
            current => current%next
            deallocate(temp)
        end do
        
        nullify(list%head)
        list%length = 0

    end subroutine empty_list


    ! function: pop
    ! returns the last cell of the list and removes it from the structure
    function pop(list) result(res)
        type(sortedlist), intent(inout) :: list
        type(cell) :: res ! copy of the cell content
        
        type(cell), pointer :: current, previous

        ! initialize result with empty/default
        res%data = ""
        res%next => null()

        if (.not. associated(list%head)) then
            print *, "warning: popping from empty list"
            return
        end if

        current => list%head
        previous => null()

        ! traverse to the end of the list
        do while (associated(current%next))
            previous => current
            current => current%next
        end do

        ! copy data to result
        res%data = current%data
        
        ! remove the node from list logic
        if (associated(previous)) then
            ! if there was a previous node, cut the link
            previous%next => null()
        else
            ! if there was no previous node, we removed the head
            list%head => null()
        end if

        ! free the memory of the removed node
        deallocate(current)
        list%length = list%length - 1

    end function pop


    ! recursive printing
    
    ! prints from last to first
    subroutine print_forward(list)
        type(sortedlist), intent(in) :: list
        call print_recursive_last_to_first(list%head)
        print * ! new line for cleanliness
    end subroutine print_forward

    ! prints from first to last
    subroutine print_reverse(list)
        type(sortedlist), intent(in) :: list
        call print_recursive_first_to_last(list%head)
        print * ! new line for cleanliness
    end subroutine print_reverse


    ! internal private recursive implementations
    
    recursive subroutine print_recursive_last_to_first(current)
        type(cell), pointer, intent(in) :: current
        
        if (associated(current)) then
            ! recursion first (go to the end)
            call print_recursive_last_to_first(current%next)
            ! then print (on the way back)
            write(*, '(a, 1x)', advance='no') current%data
        end if
    end subroutine print_recursive_last_to_first

    recursive subroutine print_recursive_first_to_last(current)
        type(cell), pointer, intent(in) :: current
        
        if (associated(current)) then
            ! print first
            write(*, '(a, 1x)', advance='no') current%data
            ! then recurse
            call print_recursive_first_to_last(current%next)
        end if
    end subroutine print_recursive_first_to_last

end module lecture_10_ex_linkedlist_slist