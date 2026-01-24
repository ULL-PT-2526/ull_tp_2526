program lecture_10_ex_main
    ! use the module defined in the previous step
    use lecture_10_ex_linkedlist_slist
    implicit none

    ! define the list variable
    type(sortedlist) :: my_list
    
    ! array with the test words provided
    character(len=11), dimension(11) :: words = [character(len=11) :: &
        'river', 'spark', 'melody', 'whisper', 'canyon', &
        'drift', 'lantern', 'echo', 'quartz', 'breeze', 'fortran']
        
    integer :: i

    ! step 1: adding words to the list
    print *, "step 1: adding words to the list..."
    
    ! loop through the words array and add them one by one
    do i = 1, 11
        print '(a, a)', " inserting: ", trim(words(i))
        call extend(my_list, trim(words(i)))
    end do
    
    print *
    print '(a, i3)', "final list length: ", my_list%length
    print *

    ! step 2: print in ascending order (a-z)
    print *, "step 2: ascending order (a -> z)"
    call print_reverse(my_list)
    print *


    ! step 3: print in descending order (z-a)
    print *, "step 3: descending order (z -> a)"
    call print_forward(my_list)
    print *

    ! step 4: clean up memory
    print *, "cleaning up memory..."
    call empty_list(my_list)
    
    print '(a, i3)', "list length after cleanup: ", my_list%length

end program lecture_10_ex_main