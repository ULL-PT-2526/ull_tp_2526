program test
    use iso_fortran_env, only: real64
    use geometry
    implicit none
    type(vector3d) :: a, b, cross, mul1, mul2, mulvv_res, div1, normal, sumvv_res, subvv_res
    type(point3d) :: d, e, sum1, sum2, sub1, sub2
    type(vector3d) :: subpp_res
    real(real64) :: s = 2.3
    real(real64) :: dist, ang

    !=============================
    ! Initialization
    !=============================
    a = vector3d(7.0, 6.0, 8.0)
    b = vector3d(4.0, 3.0, 6.0)

    d = point3d(4.0, 6.0, 9.0)
    e = point3d(15.0, 8.0, -3.0)

    !=============================
    ! Operations to test
    !=============================
    sum1 = a + d              ! vector + point
    sum2 = e + b              ! point + vector
    sumvv_res = a + b         ! vector + vector

    sub1 = a - d              ! vector - point
    sub2 = e - b              ! point - vector
    subvv_res = a - b         ! vector - vector
    subpp_res = e - d         ! point - point → vector

    mul1 = a * s              ! vector * scalar
    mul2 = s * a              ! scalar * vector
    mulvv_res = a * b         ! vector * vector (component-wise)
    div1 = a / s              ! vector / scalar

    dist = distance(d, e)     ! distance between two points
    ang = angle(a, b)         ! angle between vectors
    normal = normalize(a)     ! normalized vector
    cross = cross_product(a, b) ! cross product

    !=============================
    ! Print results
    !=============================
    print *, '--- SUMAS ---'
    print *, 'a + d = ', sum1
    print *, 'e + b = ', sum2
    print *, 'a + b = ', sumvv_res

    print *, '--- RESTAS ---'
    print *, 'a - d = ', sub1
    print *, 'e - b = ', sub2
    print *, 'a - b = ', subvv_res
    print *, 'e - d = ', subpp_res

    print *, '--- MULTIPLICACIONES ---'
    print *, 'a * s = ', mul1
    print *, 's * a = ', mul2
    print *, 'a * b (component-wise) = ', mulvv_res

    print *, '--- DIVISION ---'
    print *, 'a / s = ', div1

    print *, '--- FUNCIONES GEOMÉTRICAS ---'
    print *, 'distance(d, e) = ', dist
    print *, 'angle(a, b) [radians] = ', ang
    print *, 'normalize(a) = ', normal
    print *, 'cross_product(a, b) = ', cross

end program test
