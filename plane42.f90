module plane42

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
    public :: plane42_ke, plane42_re, plane42_ss, plane42_ke_plastic

contains

    subroutine plane42_ke(xe, young, nu, thk, ke)

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

        integer :: gpn


        real(wp), allocatable :: gauss_points(:), weights(:)



        real(wp) :: cmat(3,3), fact, aa, bb, volume
        real(wp) :: n1ze, n2ze, n3ze, n4ze, n1et, n2et, n3et, n4et, zeta, eta, wi, wj, dxdzeta, dydzeta, dxdeta, dydeta, det_J
        real(wp) :: J_(2,2), n_tylde(4, 8), G_tylde(4, 4), L(3, 4), Bmat(3, 8), test(4,8)
        integer :: i, j

        gpn = 2
        allocate(gauss_points(gpn))
        allocate(weights(gpn))


        ke = 0
        volume = 0

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


        if (gpn == 1) then
            gauss_points = [real(wp):: 0.0]
            weights = [real(wp):: 2.0]
        elseif (gpn == 2) then
            gauss_points = [real(wp):: 1.0/sqrt(3.0), -1.0/sqrt(3.0)]
            weights = [real(wp):: 1.0, 1.0]
        elseif (gpn == 3) then
            gauss_points = [real(wp):: 0, sqrt(0.6), -1.0*sqrt(0.6)]
            weights = [real(wp):: 8.0/9.0, 5.0/9.0, 5.0/9.0]
        else
            print *, 'Invalid number of gauss points'
            stop
        end if

        !print *, 'no gauss points'
        !print *, gpn

        do i=1, gpn
            zeta = gauss_points(i)
            wi = weights(i)
            do j=1, gpn
                eta = gauss_points(j)
                wj = weights(j)

                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)

                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)

                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                det_J = dxdzeta*dydeta-dydzeta*dxdeta

                !print *, det_J
                !print *, 'and matrix'
                !print *, J_

                n_tylde(1, 1:8) = [real(wp) :: n1ze, 0.0, n2ze, 0.0, n3ze, 0.0, n4ze, 0.0]
                n_tylde(2, 1:8) = [real(wp) :: n1et, 0, n2et, 0, n3et, 0, n4et, 0.0]
                n_tylde(3, 1:8) = [real(wp) :: 0.0, n1ze, 0.0, n2ze, 0.0, n3ze, 0.0, n4ze]
                n_tylde(4, 1:8) = [real(wp) :: 0.0, n1et, 0.0, n2et, 0.0, n3et, 0.0, n4et]


                G_tylde = 0
                G_tylde(1, 1:4) = [real(wp) :: J_(2,2), -J_(1,2), 0.0,0.0]
                G_tylde(2, 1:4) = [real(wp) :: -J_(2,1), J_(1,1), 0.0,0.0]
                G_tylde(3, 1:4) = [real(wp) :: 0.0, 0.0, J_(2,2), -J_(1,2)]
                G_tylde(4, 1:4) = [real(wp) :: 0.0, 0.0, -J_(2,1), J_(1,1)]

                G_tylde = (1.0/det_J)*G_tylde

                L(1, 1:4) = [real(wp) :: 1.0, 0.0, 0.0, 0.0]
                L(2, 1:4) = [real(wp) :: 0.0, 0.0, 0.0, 1.0]
                L(3, 1:4) = [real(wp) :: 0.0, 1.0, 1.0, 0.0]

                test = matmul(G_tylde, n_tylde)

                !print *, test(2, 1:8)


                Bmat = matmul(L, test)

                volume = volume + thk*det_J*wi*wj

                ke = ke + thk*matmul(transpose(Bmat), matmul(cmat, Bmat))*det_J*wi*wj
                !print*, ke(1:8,1)
            end do
        end do

    end subroutine plane42_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_re(xe, eface, fe, thk, re)

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
        integer :: gpn, i
        real(wp), allocatable :: gauss_points(:), weights(:)
        real(wp) :: n1ze, n2ze, n3ze, n4ze, n1et, n2et, n3et, n4et, zeta, eta, wi, wj, dxdzeta, dydzeta, dxdeta, dydeta, det_J, p, w
        real(wp) :: J_(2,2), n_tylde(4, 8), G_tylde(4, 4), L(3, 4), Bmat(3, 8), j(2,1)
        real(wp) :: test(8,1), cast(8)

        gpn = 2
        allocate(gauss_points(gpn))

        if (gpn == 1) then
            gauss_points = [real(wp):: 0.0]
            weights = [real(wp):: 2.0]
        elseif (gpn == 2) then
            gauss_points = [real(wp):: 1.0/sqrt(3.0), -1.0/sqrt(3.0)]
            weights = [real(wp):: 1.0, 1.0]
        elseif (gpn == 3) then
            gauss_points = [real(wp):: 0, sqrt(0.6), -1.0*sqrt(0.6)]
            weights = [real(wp):: 8.0/9.0, 5.0/9.0, 5.0/9.0]
        else
            print *, 'Invalid number of gauss points'
            stop
        end if

        nface = 0
        f = 0
        if (eface == 1) then
            eta = -1
            p = -fe
            do i=1, gpn
                nface = 0
                zeta = gauss_points(i)
                w = weights(i)
                nface(1,1) = 0.5*(1-zeta)
                nface(1,3) = 0.5*(1+zeta)
                nface(2,2) = 0.5*(1-zeta)
                nface(2,4) = 0.5*(1+zeta)
                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)
                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)
                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                j(1,1) = -J_(1,2)
                j(2,1) =  J_(1,1)
                test = matmul(transpose(nface), j)
                cast = test(1:8,1)
                re = re + p*thk*w*cast
            end do

        elseif (eface == 2) then
            zeta = 1
            p = fe
            do i=1, gpn
                nface = 0
                eta = gauss_points(i)
                w = weights(i)
                nface(1,3) = 0.5*(1-eta)
                nface(1,5) = 0.5*(1+eta)
                nface(2,4) = 0.5*(1-eta)
                nface(2,6) = 0.5*(1+eta)
                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)
                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)
                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                j(1,1) = -J_(2,2)
                j(2,1) =  J_(2,1)
                test = matmul(transpose(nface), j)
                cast = test(1:8,1)
                re = re + p*thk*w*cast
            end do
        elseif (eface == 3) then
            eta = 1
            p = fe
            do i=1, gpn
                nface = 0
                zeta = gauss_points(i)
                w = weights(i)
                nface(1,5) = 0.5*(1.0-eta)
                nface(1,7) = 0.5*(1.0+eta)
                nface(2,6) = 0.5*(1.0-eta)
                nface(2,8) = 0.5*(1.0+eta)
                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)
                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)
                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                j(1,1) = J_(1,2)
                j(2,1) =  -J_(1,1)
                test = matmul(transpose(nface), j)
                cast = test(1:8,1)
                re = re + p*thk*w*cast
            end do

        elseif (eface == 4) then
            zeta = -1
            p = -fe
            do i=1, gpn
                nface = 0
                eta = gauss_points(i)
                w = weights(i)
                nface(1,1) = 0.5*(1-eta)
                nface(1,7) = 0.5*(1+eta)
                nface(2,2) = 0.5*(1-eta)
                nface(2,8) = 0.5*(1+eta)
                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)
                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)
                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                j(1,1) = J_(2,2)
                j(2,1) =  -J_(2,1)
                test = matmul(transpose(nface), j)
                cast = test(1:8,1)
                re = re + p*thk*w*cast
            end do
        endif
        !print *, 'ERROR in plane42rect/plane42rect_re'
        !print *, 'subroutine incomplete -- you need to add some code in this subroutine'
        !stop
    end subroutine plane42_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_ss(xe, de, young, nu, estress, estrain, eprincipals)

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

        real(wp) :: n1ze, n2ze, n3ze, n4ze, n1et, n2et, n3et, n4et, zeta, eta, wi, wj, dxdzeta, dydzeta, dxdeta, dydeta, det_J
        real(wp) :: J(2,2), n_tylde(4, 8), G_tylde(4, 4), L(3, 4)

        !Gaussian point must be chose - chosen in then middle

        zeta = 0
        eta = 0
        wi = 2
        wj = 2

        n1ze = -0.25*(1.0-eta)
        n1et = -0.25*(1.0-zeta)
        n2ze = 0.25*(1.0-eta)
        n2et = -0.25*(1.0+zeta)
        n3ze = 0.25*(1.0+eta)
        n3et = 0.25*(1.0+zeta)
        n4ze = -0.25*(1.0+eta)
        n4et = 0.25*(1.0-zeta)

        dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
        dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
        dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
        dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)

        J(1, 1:2) = [dxdzeta, dydzeta]
        J(2, 1:2) = [dxdeta, dydeta]
        det_J = dxdzeta*dydeta-dxdeta*dydzeta

        n_tylde(1, 1:8) = [real(wp) :: n1ze, 0, n2ze, 0, n3ze, 0, n4ze, 0]
        n_tylde(2, 1:8) = [real(wp) :: n1et, 0, n2et, 0, n3et, 0, n4et, 0]
        n_tylde(3, 1:8) = [real(wp) :: 0, n1ze, 0, n2ze, 0, n3ze, 0, n4ze]
        n_tylde(4, 1:8) = [real(wp) :: 0, n1et, 0, n2et, 0, n3et, 0, n4et]

        G_tylde = 0
        G_tylde(1, 1:4) = [real(wp) :: J(2,2), -J(1,2), 0.0,0.0]
        G_tylde(2, 1:4) = [real(wp) :: -J(2,1), J(1,1), 0.0,0.0]
        G_tylde(3, 1:4) = [real(wp) :: 0.0, 0.0, J(2,2), -J(1,2)]
        G_tylde(4, 1:4) = [real(wp) :: 0.0, 0.0, -J(2,1), J(1,1)]
        G_tylde = (1/det_J)*G_tylde

        L(1, 1:4) = [real(wp) :: 1.0, 0.0, 0.0, 0.0]
        L(2, 1:4) = [real(wp) :: 0.0, 0.0, 0.0, 1.0]
        L(3, 1:4) = [real(wp) :: 0.0, 1.0, 1.0, 0.0]

        bmat = matmul(L, matmul(G_tylde, N_tylde))

        ! Compute element strain
        estrain = matmul(bmat, de)
        !print *, 'Element displacements'
        !print *, de
        !print *, estrain

        ! Build constitutive matrix (plane stress)

        cmat(1, 1:3) = [real(wp):: 1.0, nu, 0]
        cmat(2, 1:3) = [real(wp):: nu, 1.0, 0]
        cmat(3, 1:3) = [real(wp):: 0.0, 0.0, (1-nu)/2]

        cmat = young/(1.0-nu**2.0)*cmat
        !print *, 'this is cmat'
        !print *, cmat

        ! Compute element stress
        estress = matmul(cmat, estrain) !in a from sigma11, sigma22, sigma12
        if((xe(1) - 0.14 < 10e-8 .or. xe(3) - 0.14 < 10e-8 .or. xe(5) - 0.14 < 10e-8 .or. xe(7) - 0.14 < 10e-8 ) &
         .and. (xe(2) <10e-8 .or. xe(4) <10e-8 .or. xe(6) <10e-8 .or. xe(8) <10e-8 )) then
            print *, 'Stress A'
            print *, xe(5)
            print *, xe(6)
            print *, sqrt(estress(1)**2.0 + estress(2)**2.0 - estress(1)*estress(2) + 3.0*estress(3)**2.0)
        end if


        !print *, estress

        ! Compute principal stress and direction
        principalstresses = [0.5*(estress(1)+estress(2))+sqrt(((estress(1)-estress(2))/2.0)**2.0+estress(3)**2.0), &
        0.5*(estress(1)+estress(2))-sqrt(((estress(1)-estress(2))/2.0)**2.0+estress(3)**2.0)]

        principaldirections = [(estress(1)-estress(2))/(principalstresses(1)-principalstresses(2)), &
        (-2.0*estress(3))/(principalstresses(1)-principalstresses(2))]

        eprincipals(1) = principalstresses(1)
        eprincipals(2) = principalstresses(2)
        eprincipals(3) = atan2(principaldirections(2),principaldirections(1))/2
        !print *, 'epri'
        !print *, eprincipals


        !print *, 'WARNING in plane42rect/plane42rect_ss: subroutine incomplete -- you need to' &
         !   // 'add some code in this subroutine'
    end subroutine plane42_ss
!
!--------------------------------------------------------------------------------------------------
!


    subroutine plane42_ke_plastic(xe, young, youngy, sigma_Y_p, nu, thk, ke, stress_e_p)

        !! This subroutine constructs the plastic stiffness matrix for
        !! a rectangular 4-noded quad element.

        real(wp), intent(in) :: young, youngy, sigma_Y_p
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
        real(wp), dimension(3), intent(in) :: stress_e_p

        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix

        integer :: gpn


        real(wp), allocatable :: gauss_points(:), weights(:)



        real(wp) ::  fact, aa, bb, volume, F


        real(wp) :: n1ze, n2ze, n3ze, n4ze, n1et, n2et, n3et, n4et
        real(wp) :: zeta, eta, wi, wj, dxdzeta, dydzeta, dxdeta, dydeta, det_J, sigma_e, h
        real(wp) :: J_(2,2), n_tylde(4, 8), G_tylde(4, 4), L(3, 4), Bmat(3, 8), test(4,8), cmat_ep(3,3), c_inter(1,3), cmat_epd(1,1)
        integer :: i, j

        real(wp), dimension(3,3) :: cmat
        real(wp), dimension(1,3) :: dFdsigma

        gpn = 2
        allocate(gauss_points(gpn))
        allocate(weights(gpn))


        ke = 0
        volume = 0

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

        sigma_e = sqrt(stress_e_p(1)**2 + stress_e_p(2)**2 - stress_e_p(1)*stress_e_p(2) + 3*stress_e_p(3)**2)
        h = (youngy*young)/(young - youngy)

        F = sigma_e - sigma_Y_p
        dFdsigma(1,1) = (2*stress_e_p(1) - stress_e_p(2))/sigma_e
        dFdsigma(1,2) = (2*stress_e_p(2) - stress_e_p(1))/sigma_e
        dFdsigma(1,3) = 6*stress_e_p(3)/sigma_e

        print*,'dFdsigma', dFdsigma
        print*,'cmat',cmat
        print*,'transposedF',transpose(dFdsigma)

        cmat_epd = matmul(matmul(dFdsigma,cmat),transpose(dFdsigma))
        print*,'cmat_epd',cmat_epd
        stop
        cmat_ep = cmat - (matmul(matmul(cmat,transpose(dFdsigma)),matmul(dFdsigma,cmat)))/(cmat_epd(1,1))


        if (gpn == 1) then
            gauss_points = [real(wp):: 0.0]
            weights = [real(wp):: 2.0]
        elseif (gpn == 2) then
            gauss_points = [real(wp):: 1.0/sqrt(3.0), -1.0/sqrt(3.0)]
            weights = [real(wp):: 1.0, 1.0]
        elseif (gpn == 3) then
            gauss_points = [real(wp):: 0, sqrt(0.6), -1.0*sqrt(0.6)]
            weights = [real(wp):: 8.0/9.0, 5.0/9.0, 5.0/9.0]
        else
            print *, 'Invalid number of gauss points'
            stop
        end if


        !print *, 'no gauss points'
        !print *, gpn

        do i=1, gpn
            zeta = gauss_points(i)
            wi = weights(i)
            do j=1, gpn
                eta = gauss_points(j)
                wj = weights(j)

                n1ze = -0.25*(1.0-eta)
                n1et = -0.25*(1.0-zeta)
                n2ze = 0.25*(1.0-eta)
                n2et = -0.25*(1.0+zeta)
                n3ze = 0.25*(1.0+eta)
                n3et = 0.25*(1.0+zeta)
                n4ze = -0.25*(1.0+eta)
                n4et = 0.25*(1.0-zeta)

                dxdzeta = n1ze*xe(1) + n2ze*xe(3) + n3ze*xe(5) + n4ze*xe(7)
                dydzeta = n1ze*xe(2) + n2ze*xe(4) + n3ze*xe(6) + n4ze*xe(8)
                dxdeta = n1et*xe(1) + n2et*xe(3) + n3et*xe(5) + n4et*xe(7)
                dydeta = n1et*xe(2) + n2et*xe(4) + n3et*xe(6) + n4et*xe(8)

                J_(1, 1:2) = [dxdzeta, dydzeta]
                J_(2, 1:2) = [dxdeta, dydeta]
                det_J = dxdzeta*dydeta-dydzeta*dxdeta

                !print *, det_J
                !print *, 'and matrix'
                !print *, J_

                n_tylde(1, 1:8) = [real(wp) :: n1ze, 0.0, n2ze, 0.0, n3ze, 0.0, n4ze, 0.0]
                n_tylde(2, 1:8) = [real(wp) :: n1et, 0, n2et, 0, n3et, 0, n4et, 0.0]
                n_tylde(3, 1:8) = [real(wp) :: 0.0, n1ze, 0.0, n2ze, 0.0, n3ze, 0.0, n4ze]
                n_tylde(4, 1:8) = [real(wp) :: 0.0, n1et, 0.0, n2et, 0.0, n3et, 0.0, n4et]


                G_tylde = 0
                G_tylde(1, 1:4) = [real(wp) :: J_(2,2), -J_(1,2), 0.0,0.0]
                G_tylde(2, 1:4) = [real(wp) :: -J_(2,1), J_(1,1), 0.0,0.0]
                G_tylde(3, 1:4) = [real(wp) :: 0.0, 0.0, J_(2,2), -J_(1,2)]
                G_tylde(4, 1:4) = [real(wp) :: 0.0, 0.0, -J_(2,1), J_(1,1)]

                G_tylde = (1.0/det_J)*G_tylde

                L(1, 1:4) = [real(wp) :: 1.0, 0.0, 0.0, 0.0]
                L(2, 1:4) = [real(wp) :: 0.0, 0.0, 0.0, 1.0]
                L(3, 1:4) = [real(wp) :: 0.0, 1.0, 1.0, 0.0]

                test = matmul(G_tylde, n_tylde)

                !print *, test(2, 1:8)


                Bmat = matmul(L, test)

                volume = volume + thk*det_J*wi*wj

                if (F < 0) then
                    ke = ke + thk*matmul(transpose(Bmat), matmul(cmat, Bmat))*det_J*wi*wj
                else
                    ke = ke + thk*matmul(transpose(Bmat), matmul(cmat_ep, Bmat))*det_J*wi*wj

                end if


            end do
        end do

    end subroutine plane42_ke_plastic
!
!--------------------------------------------------------------------------------------------------
!


end module plane42
