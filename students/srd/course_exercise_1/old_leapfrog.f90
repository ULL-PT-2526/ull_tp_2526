program leapfrog
  implicit none
  integer :: i, j, k
  integer :: n
  real :: dt, t_end, t, dt_out, t_out
  real :: rs, r2, r3

  real, dimension(:), allocatable :: m
  real, dimension(:,:), allocatable :: r,v,a
  real, dimension(3) :: rji

  read *, dt
  read *, dt_out
  read *, t_end
  read *, n

  allocate(m(n))
  allocate(r(n,3))
  allocate(v(n,3))
  allocate(a(n,3))

  do i = 1, n
     read *, m(i), r(i,:),v(i,:)
  end do

  a = 0.0
  do i = 1,n
     do j = i+1,n
        rji = r(j,:) - r(i,:)
        r2 = sum(rji**2)
        r3 = r2 * sqrt(r2)
        a(i,:) = a(i,:) + m(j) * rji / r3
        a(j,:) = a(j,:) - m(i) * rji / r3
     end do
  end do

  t_out = 0.0
  do t = 0.0, t_end, dt
     v = v + a * dt/2
     r = r + v * dt
     a = 0.0
     do i = 1,n
        do j = i+1,n
           rji = r(j,:) - r(i,:)
           r2 = sum(rji**2)
           r3 = r2 * sqrt(r2)
           a(i,:) = a(i,:) + m(j) * rji / r3
           a(j,:) = a(j,:) - m(i) * rji / r3
        end do
     end do
     
     v = v + a * dt/2
     
     t_out = t_out + dt
     if (t_out >= dt_out) then
        do i = 1,n
           print*, r(i,:)
        end do
        t_out = 0.0
     end if
  end do
  
end program leapfrog
