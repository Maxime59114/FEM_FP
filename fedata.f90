module fedata
    use types
    implicit none
    save

    !! This module is used to define global values and matrices

    integer, parameter :: mdim = 8
        !! Maximum number of element degrees-of-freedom

    integer :: ne
        !! Total number of elements
    integer :: nn
        !! Total number of nodes
    integer :: nb
        !! Number of boundary conditions
    integer :: np
        !! Number of loads
    integer :: nm
        !! Number of materials
    integer :: neqn
        !! Number of equations
    integer :: bw
        !! Bandwidth of system matrix

    real(wp), dimension(:,:), allocatable :: x
        !! Nodal coordinates:
        !! * _i_-th row: coordinates for node _i_
        !! * column 1: \(x\)-coordinate
        !! * column 2: \(y\)-coordinate

    type element_def
        !! The "element" data type, which is used to
        !! describe a single element of a structure
        integer, dimension(4) :: ix
            !! Node list
        integer :: id
            !! Element id:
            !! * 1: truss element
            !! * 2: continuum element
        integer :: mat
            !! Identification number for material used by this element
        integer :: numnode
            !! Number of nodes on this element
    end type element_def
    type(element_def), dimension(:), allocatable :: element
        !! The elements of a structure

    type matprop
        !! The "material" data type, which is used to
        !! describe a material. Note that a given material
        !! might not use all of the properties defined below
        real(wp) :: young
            !! Young's Modulus
        real(wp) :: nu
            !! Poisson's Ratio
        real(wp) :: thk
            !! The thickness of the structure where this material is used
        real(wp) :: area
            !! The cross-sectional area of the structure where this material is used
        real(wp) :: dens
            !! Density
        real(wp) :: youngy
            !! Young's Modulus in the transverse direction (used for orthotropic materials)
        real(wp) :: shear
            !! Shear Modulus (used for orthotropic materials)
    end type matprop
    type(matprop), dimension(:), allocatable :: mprop
        !! The materials of a structure

    ! Boundary conditions
    real(wp), dimension(:,:), allocatable :: bound
        !! Boundary conditions
        !!
        !! * _i_-th row: boundary condition number _i_
        !! * column 1: index of node where boundary condition is applied
        !! * column 2: degree of freedom _i_ affected by boundary condition
        !!     * _i_ = 1: \(x\)-component
        !!     * _i_ = 2: \(y\)-component
        !! * column 3: prescribed value for the displacement
    real(wp), dimension(:,:), allocatable :: loads
        !! External loading
        !!
        !! * _i_-th row: external load number _i_
        !! * column 1: type _j_ of external loading
        !!     * _j_ = 1: point load
        !!     * _j_ = 2: pressure load
        !!
        !! The meaning of columns 2 - 4 depend on the type of external loading:
        !!
        !! * For point loads:
        !!     * column 2: index of node where load is applied
        !!     * column 3: degree of freedom _k_ affected by load
        !!         - _k_ = 1: \(x\)-component
        !!         - _k_ = 2: \(y\)-component
        !!      * column 4: value of applied force
        !!
        !! * For surface loads:
        !!     * column 2: index of element where pressure load is applied
        !!     * column 3: index of element face where pressure is applied
        !!     * column 4: value of applied pressure (always normal to surface)
    real(wp) :: accel(2)
        !! Acceleration forces
        !!
        !! * index 1: acceleration along \(x\)-component
        !! * index 2: acceleration along \(y\)-component

    ! Working arrays:
    real(wp), dimension(:,:), allocatable :: kmat
        !! Stiffness matrix
    real(wp), dimension(:,:), allocatable :: strain
        !! Strains at different places in the structure
        !!
        !! * _i_-th row: strain in element _i_
        !! * column 1: \(\epsilon_{11}\)
        !! * column 2: \(\epsilon_{22}\)
        !! * column 3: \(\epsilon_{12}\)
    real(wp), dimension(:,:), allocatable :: stress
        !! Stresses at different places in the structure
        !!
        !! * _i_-th row: stress in element _i_
        !! * column 1: \(\sigma_{11}\)
        !! * column 2: \(\sigma_{22}\)
        !! * column 3: \(\sigma_{12}\)
    real(wp), dimension(:),   allocatable :: p
        !! Force vector
    real(wp), dimension(:),   allocatable :: del_p
    real(wp), dimension(:),   allocatable :: d
        !! Displacement vector

    real(wp), dimension(:,:), allocatable :: principals
        !! Stresses at different places in the structure
        !!
        !! * _i_-th row: stress in element _i_
        !! * column 1: principal stress 1
        !! * column 2: principal stress 2
        !! * column 3: angle between principal stresses

    ! Input/Output
    character(len=50) :: filename
        !! Name of input file

    ! Analysis type
    character(len=20) :: antype
        !! The type of analysis to be carried out
        !!
        !! STATIC
        !! : Static linear analysis (default)
        !!
        !! STATIC_NL
        !! : Static geometrically nonlinear analysis (HINT: not implemented yet)
        !!
        !! MODAL
        !! : Analysis of eigenfrequencies and mode shapes (HINT: not implemented yet)
        !!
        !! ANGLE
        !! : Optimization of fiber orientation of anisotropic material (HINT: not implemented yet)
        !!
        !! TRANS
        !! : Transient analysis (HINT: not implemented yet)

    ! Parameters
    real(wp), parameter :: scale_def = 1.0_wp
        !! Scale deformations when plotting
    real(wp), parameter :: scale_vec = 1.0_wp
        !! Scale length of vectors for plotting
    real(wp), parameter :: scale_thk = 1.0_wp
        !! Scale thickness of lines
    logical, parameter :: banded = .false.
    !logical, parameter :: banded = .true.
        !! Indicate whether the system matrix is in banded form or not (full matrix)
    logical, parameter :: plasticity = .true.
        !! Indicate if the analysis is also plastic
    logical, parameter :: penalty = .false.
        !! Indicate whether boundary conditions are imposed by the penalty method or not (zero-one method)
    real(wp), parameter :: penalty_fac = 1.0e10_wp
        !! Scaling factor for boundary conditions imposed by the penalty method

end module fedata
