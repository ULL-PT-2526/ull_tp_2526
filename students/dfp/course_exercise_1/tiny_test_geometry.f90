program tiny_test_geometry
    use geometry
    implicit none

    type(vector3d) :: v, w, c, n
    type(vector3d) :: sumvw, diffvw, scaled, halfv
    type(point3d) :: p, q
    real(kind=dp) :: d, a

    v = vector3d(1.0_dp, 1.0_dp, 1.0_dp)
    w = vector3d(0.0_dp, 1.0_dp, 2.0_dp)
    p = point3d(0.0_dp, 0.0_dp, 0.0_dp)
    q = point3d(1.0_dp, 2.0_dp, 1.0_dp)

    print *, 'v =', v%x, v%y, v%z
    print *, 'w =', w%x, w%y, w%z

    sumvw = v + w
    diffvw = v - w
    scaled = 2.0_dp * v
    halfv = v / 2.0_dp
    print *, 'v + w =', sumvw%x, sumvw%y, sumvw%z
    print *, 'v - w =', diffvw%x, diffvw%y, diffvw%z
    print *, '2*v =', scaled%x, scaled%y, scaled%z
    print *, 'v/2 =', halfv%x, halfv%y, halfv%z

    d = distance(p,q)
    print *, 'distance(p,q)=', d

    a = angle(v,w)
    print *, 'angle(v,w)=', a

    n = normalize(vector3d(3.0_dp,4.0_dp,0.0_dp))
    print *, 'normalize(3,4,0)=', n%x, n%y, n%z

    c = crossv(v,w)
    print *, 'cross(v,w)=', c%x, c%y, c%z

end program tiny_test_geometry
