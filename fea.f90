module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

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
        allocate (p(neqn), d(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
        allocate (principals(ne, 3))

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

    subroutine plastic_iterator_1(P, n_increments, X, IX, ne, mprop, bound, &
                            D_n, stress_n, strain_n, elem_stress, elem_strain)
        ! Input variables
        REAL, INTENT(IN) :: P(:), X(:,:), mprop(:), bound(:)
        INTEGER, INTENT(IN) :: n_increments, ne, IX(:,:)

        ! Output variables
        REAL, INTENT(OUT) :: D_n(:), stress_n(:), strain_n(:), elem_stress(:), elem_strain(:)

        ! Local variables
        INTEGER :: i, elem_interest
        REAL :: delta_P(size(P)), P_n(size(P)), stress_previous(3,ne), strain_previous(3,ne)
        LOGICAL :: elastic(ne)
        REAL :: delta_D(size(P))
        REAL, ALLOCATABLE :: Ktmatr(:,:), K_num(:,:)

        REAL, dimension(3,ne) ::  stress_np, strain_np

        ! Initialize variables
        elem_interest = 2
        delta_P = P / REAL(n_increments)
        P_n = 0.0
        D_n = 0.0
        stress_n = 0.0
        strain_n = 0.0
        elastic = .TRUE.
        elem_strain = 0.0
        elem_stress = 0.0

        stress_previous = 0.0
        strain_previous = 0.0


        DO i = 1, n_increments
            PRINT *, "Iteration number: ", i
            ! Update P_n

            P_n = P_n + delta_P

            ! Compute stiffness matrix
            CALL plasticStiffness(X, ne, IX, mprop , elastic, K_num,stress_previous, strain_previous,i)
            
            ! Enforce boundary conditions
            CALL enforce(K_num, delta_P, bound, Ktmatr, delta_P)

            ! Solve for delta_D
            CALL solveSystem(Ktmatr, delta_P, delta_D)

            ! Update D_n
            D_n = D_n + delta_D

            ! Recover stress and strain
            CALL plasticRecover(IX, ne, X, mprop, delta_D, strain_n, stress_n, elastic)

            ! Store strain and stress for the element of interest
            elem_strain(i) = strain_n(elem_interest)
            elem_stress(i) = stress_n(elem_interest)
        END DO

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
        !print *, np

        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
            	!print *, 'the loads come here'
            	!print *, loads(i, :)
            	idx = (2*loads(i, 2) + loads(i, 3) - 2)
            	!print *, idx
            	p(idx) = p(idx) + loads(i, 4)
                !p(6) = -0.1_wp
                !print *, 'WARNING in fea/buildload: You need to replace hardcoded nodal load with your code'
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
                !print *, 're'
                !print *, re
                !print *, 'p'
                !print *, p

                do k = 1, nen*2
                    p(edof(k)) = p(edof(k)) + re(k)
                end do

            	! Build uniformly distributed surface (pressure) load contribution
                !print *, 'ERROR in fea/buildload'
                !print *, 'Distributed loads not defined -- you need to add your own code here'
                !stop !restore this stop if needed - idk what it was for
            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
        !print *, 'This is p'
        !print *, p
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
! Hint for system matrix in band form:
!        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
! Hint for modal analysis:
!        real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, area
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0
        else
            kmat = 0
            print *, 'Stiffnes band matrix reset'
            !print*,'ERROR in fea/buildstiff'
            !print*,'Band form not implemented -- you need to add your own code here'
            !stop
        end if

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
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                young = mprop(element(e)%mat)%young
                nu = mprop(element(e)%mat)%nu
                thk = mprop(element(e)%mat)%thk
                call plane42_ke(xe, young, nu, thk, ke)
                 !print *, 'ERROR in fea/buildstiff:'
                 !print *, 'Stiffness matrix for plane42rect elements not implemented -- you need to add your own code here'
                 !stop
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

                !stop
                print *, 'The global matrix assembly was implemented - check correctness'
                !print *, 'ERROR in fea/buildstiff'
                !print *, 'Band form not implemented -- you need to add our own code here'
                !stop
            end if
        end do
        !print *, kmat(8,1:12)
        !stop
        !print *, 'end buildstiff'
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
                    !print *, 'spec2'
                    !print *, kmat(1:spc2, idof)
                    !print *, 'ps'
                    !print *, p(idof:spc)
                    p(idof:spc) = p(idof:spc) - kmat(1:spc2, idof) * bound(i, 3)

                    p(idof) = bound(i, 3)
                    kmat(1:bw,idof) = 0

                    do j = 1, idof
                        kmat(idof-j+1, j) = 0
                    end do
                    kmat(1, idof) = 1
                end do
                !print *, 'The changed matrix'
                !print *, kmat(8, 1:12)
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(1, idof) = kmat(1, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
            print *, 'Implemented the boundary condition enforcement - check for correctness'
            !print *, 'ERROR in fea/enforce'
            !print *, 'Band form not implemented -- you need to add your own code here'
            !stop
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
                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
            case( 2 )
                young = mprop(element(e)%mat)%young
                nu = mprop(element(e)%mat)%nu
                call plane42_ss(xe, de, young, nu, estress, estrain, eprincipals)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
                principals(e, 1:3) = eprincipals
                !print *, 'WARNING in fea/recover: Stress and strain not calculated for continuum' &
                    !// 'elements -- you need to add your own code here'
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
!
!--------------------------------------------------------------------------------------------------
!

    subroutine plasticStiffness(elastic, K_num,stress_previous, strain_previous,i)

        use fedata
        use link1
        use plane42rect
        use plane42

        integer :: e, j
        integer :: nen

        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke

        integer, intent(in) :: i



        real(wp), dimension(:,:), intent(inout) :: K_num

        real(wp) :: young, youngt, nu, dens, thk, sigma_Y




        real(wp) :: stress_e_p(3)
        K_num = 0


        do e = 1,ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode

            do j = 1, nen
                 xe(2*j-1) = x(element(e)%ix(j),1)
                 xe(2*j  ) = x(element(e)%ix(j),2)
                 edof(2*j-1) = 2 * element(e)%ix(j) - 1
                 edof(2*j)   = 2 * element(e)%ix(j)
            end do

            young = mprop(element(e)%mat)%young
            youngt = mprop(element(e)%mat)%youngt
            nu = mprop(element(e)%mat)%nu
            thk = mprop(element(e)%mat)%thk
            sigma_Y = mprop(element(e)%mat)%sigma_Y


            stress_e_p(1) = stress_previous(1,e)
            stress_e_p(2) = stress_previous(2,e)
            stress_e_p(3) = stress_previous(3,e)


            call plane42_ke_plastic(xe, young, youngt, sigma_Y, nu, thk, ke, stress_e_p)

        end do

end module fea




