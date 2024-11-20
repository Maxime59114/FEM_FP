module link1

    !! This module contains subroutines specific to the link1 element
    !!
    !! The link1 element has 2 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the x- and y-coordinate directions.
    !!
    !!             o N2
    !!            /
    !!           / E1
    !!          /
    !!         o N1
    !!
    !! `N1` = element node 1
    !! `N2` = element node 2
    !! `E1` = element edge 1
    use types
    implicit none
    save

    private
    public :: link1_ke, link1_ss

contains

    subroutine link1_ke(xe, young, area, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a truss element

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: area
            !! Cross-sectional area of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinate of this element in undeformed configuration
            !!
            !! * `xe(1)` = \(x\)-coordinate of element node 1
            !! * `xe(2)` = \(y\)-coordinate of element node 1
            !! * `xe(3)` = \(x\)-coordinate of element node 2
            !! * `xe(4)` = \(y\)-coordinate of element node 2
            !!
            !! See also [[link1]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Element stiffness matrix

        real(wp) :: delx, dely, l0, delx2, dely2, delxy, coef

        delx = xe(3) - xe(1)
        dely = xe(4) - xe(2)
        delx2 = delx * delx
        dely2 = dely * dely
        delxy = delx * dely
        l0 = sqrt(delx2 + dely2)

        coef = young*area/l0**3

        ke(1,1) = coef * delx2
        ke(1,2) = coef * delxy
        ke(1,3) = -coef * delx2
        ke(1,4) = -coef * delxy
        ke(2,1) = coef * delxy
        ke(2,2) = coef * dely2
        ke(2,3) = -coef * delxy
        ke(2,4) = -coef * dely2
        ke(3,1) = -coef * delx2
        ke(3,2) = -coef * delxy
        ke(3,3) = coef * delx2
        ke(3,4) = coef * delxy
        ke(4,1) = -coef * delxy
        ke(4,2) = -coef * dely2
        ke(4,3) = coef * delxy
        ke(4,4) = coef * dely2
    end subroutine link1_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine link1_ss(xe, de, young, estress, estrain)

        !! This subrotuine constructs the element stress and strain

        real(wp), intent(in) :: young
            !! Young's Modulus of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinate of this element in undeformed configuration (see alse [[link1_ke]])
        real(wp), dimension(:), intent(in) :: de
            !! Displacement field
            !!
            !! * de(1) = displacement along x-axis of element node 1
            !! * de(2) = displacement along y-axis of element node 1
            !! * de(3) = displacement along x-axis of element node 2
            !! * de(4) = displacement along y-axis of element node 2
        real(wp), dimension(:), intent(out) :: estress
            !! Element stress
            !!
            !! * estress(1) = stress on section perpendicular to element axis
            !! * estress(2:) = unused
        real(wp), dimension(:), intent(out) :: estrain
            !! Element strain
            !!
            !! * estrain(1) = strain measured along element axis
            !! * estrain(2:) = unused
        real(wp) :: delx, dely, l0, delx2, dely2, delu, delv

        delx = xe(3) - xe(1)
        dely = xe(4) - xe(2)
        delx2 = delx * delx
        dely2 = dely * dely
        l0 = sqrt(delx2 + dely2)

        delu = de(3) - de(1)
        delv = de(4) - de(2)

        estrain(1) = (delu*delx + delv*dely)/l0**2
        estrain(2:3) = 0
        estress(1) = young * estrain(1)
        estress(2:3) = 0
    end subroutine link1_ss

end module link1
