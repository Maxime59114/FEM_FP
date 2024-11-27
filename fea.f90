module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, plastic_iterator_1, plastic_enforce

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial

        !! This subroutine is mainly used to allocate vectors and matrices

        use fedata
        use link1
        use plane42rect
        use plane42

        integer:: e

! Hint for continuum elements:
        !integer, parameter :: mdim = 8
!        integer, dimension(mdim) :: edof

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            bw = 0
            do e = 1, ne
                if (bw < (maxval(element(e)%ix) - minval(element(e)%ix) + 1)*2) then
                    bw = (maxval(element(e)%ix) - minval(element(e)%ix) + 1)*2
                end if
            end do
            allocate(kmat(bw,neqn))
            print *, 'allocation implemented for band'
        end if
        allocate (p(neqn), d(neqn), del_p(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
        allocate (principals(ne, 3))
        allocate (sigma_n(ne, 3))
        allocate (sigma_yield(ne))

        ! Initial stress and strain
        strain = 0
        stress = 0
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e
        real(wp), dimension(:), allocatable :: plotval, psv1, psv2, psvang
        real(wp), allocatable :: kmat_old(:,:)

        if (plasticity) then
            call plastic_iterator_1
            stop
        end if

        allocate(kmat_old(bw,neqn))

        ! Build load-vector
        call buildload
        print *, 'Buildload done'

        ! Build stiffness matrix
        call buildstiff
        print *, 'Buildstiff done'
        kmat_old = kmat

        ! Remove rigid body modes
        call enforce
        print *, 'Enforce done'

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
            call bfactor(kmat)
            call bsolve(kmat, p)
            !print *, 'ERROR in fea/displ'
            !print *, 'Band form not implemented -- you need to add your own code here'
            !stop
        end if
        print *, 'Solving done'

        ! Transfer results
        d(1:neqn) = p(1:neqn)

        call compliance(kmat_old)

        ! Recover stress
        call recover
        print *, 'Recovery done'

        ! Output results
        call output

        ! Plot deformed structure
        call plotmatlabdef('Deformed')

        ! Plot element values
        allocate (plotval(ne))
        allocate (psv1(ne))
        allocate (psv2(ne))
        allocate (psvang(ne))
        psv1 = 0
        psv2 = 0
        psvang = 0
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                !print *, 'Finding vm stresses for continuum elements'
                !print*, stress(e, 1:3)
                plotval(e) = sqrt(stress(e,1)**2.0 + stress(e,2)**2.0 - stress(e,1)*stress(e,2) + 3.0*stress(e,3)**2.0) !THIS HAS TO BE ADDED
                psv1(e) = principals(e,1)
                psv2(e) = principals(e,2)
                psvang(e) = principals(e,3)
                !print *, 'WARNING in fea/displ: Plot value not set -- you need to add your own code here'
            end if
        end do
        call plotmatlabeval('Stresses',plotval)
        !Here we need to plot the thing for principals
        call plotmatlabevec('Principal stresses and direction',psv1,psv2,psvang)

    end subroutine displ

    subroutine plastic_iterator_1

        use fedata
        use numeth
        use processor

        integer :: n_increments, i
        real(wp) :: p_n(size(p)), delta_p(size(p)), deld(size(p))

        n_increments = 10

        call buildload !Makes sure vector p is correct

        p_n = 0.0
        del_p = 0.0
        delta_p = p/REAL(n_increments)
        deld = 0.0
        d = 0.0

        sigma_n = 0.0
        sigma_yield = 0.1 !This has to be initialized in agreement with data

        DO i = 1, n_increments
            PRINT *, "Iteration number: ", i

            ! Update P_n
            p_n = p_n + delta_p
            del_p = delta_p

            ! Compute stiffness matrix
            call buildstiff

            ! Enforce boundary conditions
            call plastic_enforce

            ! Solve for delta_D
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, del_p)
            ! Transfer results
            deld(1:neqn) = del_p(1:neqn)
            ! Update D_n
            d = d + deld
            stop
            ! Recover stress and strain
            call recover
        END DO
        !!Process the results

   end subroutine plastic_iterator_1
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42rect
        use plane42

        integer :: i
        integer :: j
        integer :: e
        integer :: eface
        integer :: nen
        integer :: idx
        integer :: k
        real(wp) :: fe
        real(wp) :: thk
        integer :: face
! Hint for continuum elements:
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re

        ! Build load vector
        re = 0.0
        p(1:neqn) = 0

        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	idx = (2*loads(i, 2) + loads(i, 3) - 2)
            	p(idx) = p(idx) + loads(i, 4)
            case( 2 )
                print *, 'surface traction'
                e = loads(i, 2)
                nen = element(e)%numnode
                do j = 1, nen
                    xe(2*j-1) = x(element(e)%ix(j),1)
                    xe(2*j  ) = x(element(e)%ix(j),2)
                    edof(2*j-1) = 2 * element(e)%ix(j) - 1
                    edof(2*j)   = 2 * element(e)%ix(j)
                end do
                fe = loads(i, 4)
                thk = mprop(element(e)%mat)%thk
                eface = loads(i, 3)
                call plane42_re(xe, eface, fe, thk, re)
                do k = 1, nen*2
                    p(edof(k)) = p(edof(k)) + re(k)
                end do
            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rect
        use plane42

        integer :: e, i, j
        integer :: nen
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        kmat = 0

        do e = 1, ne
            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 if (plasticity) then
                    !call link1_ke_plastic()
                    print *, "No plastic implementation"
                    stop
                 else
                    call link1_ke(xe, young, area, ke)
                 end if
            case( 2 )
                young = mprop(element(e)%mat)%young
                nu = mprop(element(e)%mat)%nu
                thk = mprop(element(e)%mat)%thk
                if (plasticity) then
                    call plane42_ke_plastic(xe, young, sigma_yield(e), nu, thk, ke, sigma_n(e, 1:3))
                 else
                    call plane42_ke(xe, young, nu, thk, ke)
                 end if
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do
! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
            else
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        !print *, edof(i)
                        !print *, edof(j)
                        if (edof(i) - edof(j) >= 0) then
                            if (edof(i) - edof(j) < bw) then
                                kmat(edof(i) - edof(j) + 1, edof(j)) = kmat(edof(i) - edof(j) + 1, edof(j)) + ke(i, j)
                            end if
                        end if
                    end do
                end do
            end if
        end do
    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof, j, k, spc, spc2
        real(wp) :: penal

        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))

                    spc = idof+bw-1

                    if (spc >bw) then
                        spc = bw
                    end if

                    spc2 = spc - idof + 1
                    p(idof:spc) = p(idof:spc) - kmat(1:spc2, idof) * bound(i, 3)

                    p(idof) = bound(i, 3)
                    kmat(1:bw,idof) = 0

                    do j = 1, idof
                        kmat(idof-j+1, j) = 0
                    end do
                    kmat(1, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(1, idof) = kmat(1, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
            print *, 'Implemented the boundary condition enforcement - check for correctness'
        end if
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain,
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42rect
        use plane42

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area
! Hint for continuum elements:
        real(wp):: nu, dens, thk
        real(wp), dimension(3) :: estrain, estress, eprincipals

        ! Reset force vector
        p = 0

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                if (plasticity) then
                    !call link1_ke_plastic()
                    print *, "NO implementation for truss"
                    stop
                    p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                    !call link1_ss_plastic()
                else
                    call link1_ke(xe, young, area, ke)
                    p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                    call link1_ss(xe, de, young, estress, estrain)
                end if
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
            case( 2 )
                young = mprop(element(e)%mat)%young
                nu = mprop(element(e)%mat)%nu
                if (plasticity) then
                    youngt = mprop(element(e)%mat)%youngt
                    yieldstress = mprop(element(e)%mat)%yieldstress

                    do i = 1, nen
                        delta_de_n(2*i-1,1) = del_d(edof(2*i-1))
                        delta_de_n(2*i,1)   = del_d(edof(2*i))

                        estress_p(1,1) = stress_p(element(e)%ix(i),1)
                        estress_p(2,1) = stress_p(element(e)%ix(i),2)
                        estress_p(3,1) = stress_p(element(e)%ix(i),3)

                        estrain_p(1,1) = strain_p(element(e)%ix(i),1)
                        estrain_p(2,1) = strain_p(element(e)%ix(i),2)
                        estrain_p(3,1) = strain_p(element(e)%ix(i),3)

                    end do

                    esigma_Y_p = sigma_yield(e)

                    call plane42_ss_plastic(xe, delta_de_n, young, youngt, nu,estress_p,estress_n, estrain_p, &
                                            esigma_Y_p, esigma_Y_n)
                    stress(e,1) = stress(e,1) + estress_n(1,1)
                    stress(e,2) = stress(e,2) + estress_n(2,1)
                    stress(e,3) = stress(e,3) + estress_n(3,1)

                    stress_p = stress
                    sigma_yield(e) = esigma_Y_n
                 else
                    call plane42_ss(xe, de, young, nu, estress, estrain, eprincipals)
                 end if
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
                principals(e, 1:3) = eprincipals
            end select
        end do
    end subroutine recover

    subroutine compliance(kmat_old)
        use fedata
        real(wp), dimension(:,:), intent(in) :: kmat_old
        real (wp) :: C
        C = dot_product(d, matmul(d, kmat_old))
        !k matrix not factorized and no boundary conditions
        print *, 'compliance'
        print *, C

    end subroutine


    subroutine plastic_enforce

        !! This subroutine enforces the support boundary conditions - for plastic case

        use fedata

        integer :: i, idof, j, k, spc, spc2
        real(wp) :: penal

        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    del_p(1:neqn) = del_p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    del_p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    del_p(idof) = penal * bound(i, 3)
                end do
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))

                    spc = idof+bw-1

                    if (spc >bw) then
                        spc = bw
                    end if

                    spc2 = spc - idof + 1
                    del_p(idof:spc) = del_p(idof:spc) - kmat(1:spc2, idof) * bound(i, 3)

                    del_p(idof) = bound(i, 3)
                    kmat(1:bw,idof) = 0

                    do j = 1, idof
                        kmat(idof-j+1, j) = 0
                    end do
                    kmat(1, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(1, idof) = kmat(1, idof) + penal
                    del_p(idof) = penal * bound(i, 3)
                end do
            end if
        end if
    end subroutine plastic_enforce

end module fea
