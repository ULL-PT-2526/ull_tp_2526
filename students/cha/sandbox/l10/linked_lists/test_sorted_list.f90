program test_sorted_list

  use slist
  implicit none

  type(SortedList) :: my_list
  character(len=20), dimension(10) :: words
  integer :: i

  ! Sample words
  words = [character(len=20) :: &
       "river", "spark", "melody", "whisper", "canyon", &
       "drift", "lantern", "echo", "quartz", "breeze"]

  ! Initialize the list structure
  my_list%length = 0
  my_list%head => null()

  ! Insert each word into the sorted list
  do i = 1, size(words)
    call extend(my_list, trim(words(i)))
  end do

  ! Print list in ascending order (A → Z)
  print *, "List in ascending order (A-Z):"
  call print_forward(my_list%head)
  print *, ""

  ! Print list in descending order (Z → A)
  print *, "List in descending order (Z-A):"
  call print_reverse(my_list%head)
  print *, ""

  ! Show final length of the list
  print *, "Final list length:", my_list%length

  ! Empty and free all memory
  call empty_list(my_list)
  print *, "List emptied successfully."

end program test_sorted_list
