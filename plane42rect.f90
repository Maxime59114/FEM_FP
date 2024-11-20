module plane42rect

    !! This module contains subroutines specific to the plane42 element
    !!
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc
    !! `E1` = element face 1, `E2` = element face 2, etc

    use types
    implicit none
    save

    private
    public :: plane42rect_ke, plane42rect_re, plane42rect_ss

contains

    subroutine plane42rect_ke(xe, young, nu, thk, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        real(wp) :: cmat(3,3), fact, aa, bb
        real(wp) :: t1, t2, t3, t4, t5, t6, t9, t11, t13, t15, &
                        t16, t17, t18, t19, t22, t25, t29, t33, &
                        t34, t38, t42, t43, t47, t51, t55, t56, &
                        t57, t58, t59, t62, t66, t70, t74, t78, &
                        t82, t86, t90, t94, t98, t102, t106, t110
        real(wp) :: d11, d12, d13, d22, d23, d33

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        d11 = cmat(1,1)
        d12 = cmat(1,2)
        d13 = cmat(1,3)
        d22 = cmat(2,2)
        d23 = cmat(2,3)
        d33 = cmat(3,3)

        t1 = bb**2
        t2 = t1*d11
        t3 = 2*t2
        t4 = aa**2
        t5 = t4*d33
        t6 = 2*t5
        t9 = 3*d13*aa*bb
        t11 = 1/aa
        t13 = 1/bb
        t15 = (t3+t6+t9)*t11*t13/6
        t16 = d13*t1
        t17 = 4*t16
        t18 = t4*d23
        t19 = 4*t18
        t22 = 3*d12*aa*bb
        t25 = 3*aa*d33*bb
        t29 = (t17+t19+t22+t25)*t11*t13/12
        t33 = (-t3+t5)*t11*t13/6
        t34 = 2*t18
        t38 = (-t17+t34+t22-t25)*t11*t13/12
        t42 = (t2+t5+t9)*t11*t13/6
        t43 = 2*t16
        t47 = (t43+t34+t22+t25)*t11*t13/12
        t51 = (-t2+t6)*t11*t13/6
        t55 = (-t43+t19+t22-t25)*t11*t13/12
        t56 = d33*t1
        t57 = 2*t56
        t58 = t4*d22
        t59 = 2*t58
        t62 = 3*aa*d23*bb
        t66 = (t57+t59+t62)*t11*t13/6
        t70 = (-t17+t34-t22+t25)*t11*t13/12
        t74 = (-t57+t58)*t11*t13/6
        t78 = (t56+t58+t62)*t11*t13/6
        t82 = (-t43+t19-t22+t25)*t11*t13/12
        t86 = (-t56+t59)*t11*t13/6
        t90 = (t3+t6-t9)*t11*t13/6
        t94 = (t17+t19-t22-t25)*t11*t13/12
        t98 = (t2+t5-t9)*t11*t13/6
        t102 = (t43+t34-t22-t25)*t11*t13/12
        t106 = (t57+t59-t62)*t11*t13/6
        t110 = (t56+t58-t62)*t11*t13/6
        ke(1,1) = t15
        ke(1,2) = t29
        ke(1,3) = t33
        ke(1,4) = t38
        ke(1,5) = -t42
        ke(1,6) = -t47
        ke(1,7) = -t51
        ke(1,8) = -t55
        ke(2,1) = t29
        ke(2,2) = t66
        ke(2,3) = t70
        ke(2,4) = t74
        ke(2,5) = -t47
        ke(2,6) = -t78
        ke(2,7) = -t82
        ke(2,8) = -t86
        ke(3,1) = t33
        ke(3,2) = t70
        ke(3,3) = t90
        ke(3,4) = t94
        ke(3,5) = -t51
        ke(3,6) = -t82
        ke(3,7) = -t98
        ke(3,8) = -t102
        ke(4,1) = t38
        ke(4,2) = t74
        ke(4,3) = t94
        ke(4,4) = t106
        ke(4,5) = -t55
        ke(4,6) = -t86
        ke(4,7) = -t102
        ke(4,8) = -t110
        ke(5,1) = -t42
        ke(5,2) = -t47
        ke(5,3) = -t51
        ke(5,4) = -t55
        ke(5,5) = t15
        ke(5,6) = t29
        ke(5,7) = t33
        ke(5,8) = t38
        ke(6,1) = -t47
        ke(6,2) = -t78
        ke(6,3) = -t82
        ke(6,4) = -t86
        ke(6,5) = t29
        ke(6,6) = t66
        ke(6,7) = t70
        ke(6,8) = t74
        ke(7,1) = -t51
        ke(7,2) = -t82
        ke(7,3) = -t98
        ke(7,4) = -t102
        ke(7,5) = t33
        ke(7,6) = t70
        ke(7,7) = t90
        ke(7,8) = t94
        ke(8,1) = -t55
        ke(8,2) = -t86
        ke(8,3) = -t102
        ke(8,4) = -t110
        ke(8,5) = t38
        ke(8,6) = t74
        ke(8,7) = t94
        ke(8,8) = t106

        ke = ke*thk
    end subroutine plane42rect_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42rect_re(xe, eface, fe, thk, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: aa, bb, nface(2,8), f(2)

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        nface = 0
        f = 0
        if (eface == 1) then
            nface(1,1) = aa
            nface(1,3) = aa
            nface(2,2) = aa
            nface(2,4) = aa
            f(2) = -fe
        elseif (eface == 2) then
            nface(1,3) = bb
            nface(1,5) = bb
            nface(2,4) = bb
            nface(2,6) = bb
            f(1) = fe
        elseif (eface == 3) then
            nface(1,5) = aa
            nface(1,7) = aa
            nface(2,6) = aa
            nface(2,8) = aa
            f(2) = fe
        elseif (eface == 4) then
            nface(1,1) = bb
            nface(1,7) = bb
            nface(2,2) = bb
            nface(2,8) = bb
            f(1) = -fe
        endif
        re = matmul(transpose(nface), f) * thk
        !print *, 'ERROR in plane42rect/plane42rect_re'
        !print *, 'subroutine incomplete -- you need to add some code in this subroutine'
        !stop
    end subroutine plane42rect_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42rect_ss(xe, de, young, nu, estress, estrain, eprincipals)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)
        real(wp) :: bmat(3, 8), cmat(3, 3)
        real(wp) :: principalstresses(2), principaldirections(2)
        real(wp), dimension(:), intent(out) :: eprincipals

        real(wp) :: x, y, a, b !the coordinates in element for which the strain is evaluated and the size of the element
        real(wp) :: n1x, n1y, n2x, n2y, n3x, n3y, n4x, n4y

        x = 0
        y = 0

        a = (xe(3)-xe(1))/2
        b = (xe(8)-xe(2))/2

        n1x = -1*(b-y)
        n1y = -1*(a-x)
        n2x = (b-y)
        n2y = -1*(a+x)
        n3x = (b+y)
        n3y = (a+x)
        n4x = -1*(b+y)
        n4y = (a-x)


        ! Build strain-displacement matrix
        bmat(1,:) = [ real(wp):: n1x, 0.0, n2x, 0.0, n3x, 0.0, n4x, 0.0]
        bmat(2,:) = [ real(wp):: 0.0, n1y, 0.0, n2y, 0.0, n3y, 0.0, n4y]
        bmat(3,:) = [ real(wp):: n1y, n1x, n2y, n2x, n3y, n3x, n4y, n4x]
        bmat = bmat/(4*a*b)
        !print *, 'this is bmat'
        !print *, bmat(1,1:8)

        ! Compute element strain
        estrain = matmul(bmat, de)
        !print *, 'Element displacements'
        !print *, de
        !print *, estrain

        ! Build constitutive matrix (plane stress)

        cmat = young/(1.0-nu**2.0)*transpose(reshape([real(wp):: 1.0, nu, 0.0, nu, 1.0, &
         0.0, 0.0, 0.0, (1.0-nu)**2.0 ], shape(cmat)))
        !print *, 'this is cmat'
        !print *, cmat

        ! Compute element stress
        estress = matmul(cmat, estrain) !in a from sigma11, sigma22, sigma12
        !print *, estress

        ! Compute principal stress and direction
        principalstresses = [0.5*(estress(1)+estress(2))+sqrt(((estress(1)-estress(2))/2.0)**2.0+estress(3)**2.0), &
        0.5*(estress(1)+estress(2))-sqrt(((estress(1)-estress(2))/2.0)**2.0+estress(3)**2.0)]

        principaldirections = [(estress(1)-estress(2))/(principalstresses(1)-principalstresses(2)), &
        (-2.0*estress(3))/(principalstresses(1)-principalstresses(2))]

        eprincipals(1) = principalstresses(1)
        eprincipals(2) = principalstresses(2)
        eprincipals(3) = atan2(principaldirections(2),principaldirections(1))/2
        print *, 'epri'
        print *, eprincipals


        !print *, 'WARNING in plane42rect/plane42rect_ss: subroutine incomplete -- you need to' &
         !   // 'add some code in this subroutine'
    end subroutine plane42rect_ss

end module plane42rect
