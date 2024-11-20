module types

    !! In this module we define so-called "kind"-parameters, which can be used
    !! to control the precision of the floating point reals in the program. 
    !!
    !! Also we define the value of some constants
    !!
    !! ** Note:  It is not necessary to alter in this module **
    
    ! Get "kind" parameters for single, double, and quad precision floating-point reals
    ! (ref: http://fortranwiki.org/fortran/show/Real+precision)
    integer, parameter ::                              &
        sp = kind(1.0),                                &
        dp = selected_real_kind(2*precision(1.0_sp)),  &
        qp = selected_real_kind(2*precision(1.0_dp))
    !     
    ! Select working precision
    integer, parameter :: wp = dp
    
    ! Constants
!    real(wp), parameter ::  pi = 4*atan(1.0_wp) ! Fortran 2003+
    real(wp), parameter ::  pi = 3.14159265358979323848_wp ! Fortran 90/95
    
end module types
