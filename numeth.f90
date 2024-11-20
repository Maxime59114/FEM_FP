module numeth

    !! This module contains subroutines for numerical solution
    !! of a system of linear equations

    use types
    implicit none
    save

    private
    public :: factor, solve, bfactor, bsolve

contains

    subroutine factor(a)

        !! This subroutine computes the LU-factorization of a matrix. 
        !!
        !! The following assumptions are made
        !!
        !! * the matrix is square
        !! * the matrix is non-singular
        !!
        !! If the assumptions are not meet the behavior is undefined -- probably the
        !! program will crash.
        !!
        !! The LU-factorization is stored in-place, thereby replacing the original
        !! content of the given matrix.

        real(wp), dimension(:, :), intent(inout) :: a
            !! Matrix to be LU-factored
            
        integer :: i, j, neqn
        neqn = size(a, 1)

        do j = 2, neqn
            do i = 1, j-1
                a(j, i) = a(j, i) - dot_product(a(j, 1:i-1), a(1:i-1, i))
                a(i, j) = a(i, j) - dot_product(a(i, 1:i-1), a(1:i-1, j))
                a(j, i) = a(j, i)/a(i, i)
            end do
            a(j, j) = a(j, j) - dot_product(a(j, 1:j-1), a(1:j-1, j))
        end do
    end subroutine factor
!
!--------------------------------------------------------------------------------------------------
!
    subroutine solve(a, b)

        !! This subroutine solves the linear system of equations [M]{x} = {b}, where
        !! [a] is the LU-factoriztion of [M], as obtained from the subroutine [[factor]].  
        !!   
        !! The subroutine returns the solution {x} by overwriting {b}.
        
        real(wp), dimension(:, :), intent(in) :: a
            !! LU-factorization of system matrix
        real(wp), dimension(:), intent(inout) :: b
            !! Right-hand side to be overwritten by the computed solution

        integer :: i, neqn
        neqn = size(a, 1)

        ! Forward substitution
        do i = 2, neqn
            b(i) = b(i) - dot_product(a(i, 1:i-1), b(1:i-1))
        end do

        ! Backward substitution
        b(neqn) = b(neqn)/a(neqn, neqn)
        do i = neqn-1, 1, -1
            b(i) = (b(i) - dot_product(a(i, i+1:neqn), b(i+1:neqn)))/a(i, i)
        end do
    end subroutine solve
 !
!--------------------------------------------------------------------------------------------------
!
    subroutine bfactor(a)

        !! This subroutine computes the LU-factorization of a matrix.
        !!
        !! The following assumptions are made
        !!
        !! * the matrix is stored in banded form (rectangular matrix)
        !! * the matrix is non-singular
        !!
        !! If the assumptions are not meet the behavior is undefined -- probably the
        !! program will crash.
        !!
        !! The LU-factorization is stored in-place, thereby replacing the original
        !! content of the given matrix.

        real(wp), dimension(:, :), intent(inout) :: a        
            !! Matrix in banded form to be LU-factored
        real(wp) :: c
        integer :: i, j, k, m, n, bw, neqn

        bw = size(a, 1)
        neqn = size(a, 2)

        do n = 1, neqn
            do m = 2, bw
                if (a(m, n) == 0.0_wp) cycle
                i = n+m-1
                c = a(m, n)/a(1, n)
                j = 0
                do k = m, bw
                    j = j+1
                    a(j, i) = a(j, i) - c*a(k, n)
                end do
                a(m, n) = c
            end do
        end do 
    end subroutine bfactor
!
!--------------------------------------------------------------------------------------------------
!
    subroutine bsolve(a, b)
       
        !! This subroutine solves the linear system of equations [M]{x} = {b}, where
        !! [a] is the LU-factoriztion of [M], as obtained from the subroutine [[bfactor]].  
        !!   
        !! The subroutine returns the solution {x} by overwriting {b}.

        real(wp), dimension(:, :), intent(in) :: a
            !! LU-factorization of system matrix
        real(wp), dimension(:), intent(inout) :: b
            !! Right-hand side to be overwritten by the computed solution

        integer :: i, j, k, m, n, bw, neqn
        bw = size(a, 1)
        neqn = size(a, 2)

        do n = 1,neqn
            do j = 2, bw
                if (a(j, n) == 0.0_wp) cycle 
                i = n+j-1
                b(i) = b(i)-a(j, n)*b(n)
            end do
            b(n) = b(n)/a(1, n)
        end do

        do m = 2, neqn
            n = neqn+1-m
            do j = 2, bw
                if (a(j, n) == 0.0_wp) cycle 
                k = n+j-1
                b(n) = b(n)-a(j, n)*b(k)
            end do
        end do
    end subroutine bsolve

end module numeth

