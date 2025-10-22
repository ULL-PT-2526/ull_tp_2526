program palindrome
  implicit none
  integer :: num, original, reversed, digit
  
  print *, "Enter a number:"
  read(*,*) num
  
  original = num
  reversed = 0
  
  do while (num > 0)
     digit = mod(num, 10)
     reversed = reversed * 10 + digit
     num = num / 10
  end do
  
  print *, "Original number:", original
  print *, "Reversed number:", reversed
  
  if (original == reversed) then
     print *, "It's a palindrome"
  else
     print *, "It's not a palindrome"
  end if
  
end program palindrome