program array
  implicit none
  integer :: A(3,4)  
  integer :: B(3,3)
  integer, allocatable :: C(:,:)

  
  A = reshape( (/ 3,2,1,  2,4,2,  4,2,3,  1,2,7 /), shape(A))

  B = reshape( (/ 3,2,3,  2,1,0,  4,2,2 /), shape(B))

  print *, "Matrix A:"
  call print_matrix(A)

  print *, "Matrix B:"
  call print_matrix(B)


  if (size(B,2) == size(A,1)) then
     allocate(C(size(B,1), size(A,2)))
     C = matmul(B,A)
     print *, "Resulting B x A = C:"
     call print_matrix(C)
     deallocate(C)
  else
     print *, "Cannot perform B x A: dimensions not compatible"
     print *, "B columns:", size(B,2), "A rows:", size(A,1)
  end if

  contains

  subroutine print_matrix(mat)
    integer, intent(in) :: mat(:,:)
    integer :: i, j
    do i = 1, size(mat,1)
       do j = 1, size(mat,2)
          write(*,'(I4)',advance="no") mat(i,j)
       end do
       print *
    end do
  end subroutine print_matrix

end program array


