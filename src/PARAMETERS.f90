MODULE parameters

use constants

!**********************
!!INITIAL DEFINITIONS *
!**********************

!*Define chain topology*!:
!* Check one *!
logical, parameter :: linear = .true.
logical, parameter :: star = .false.
logical, parameter :: lattice = .false.
logical, parameter :: squared = .false.

!*Define couplings scheme*!:
!* Check one *!
logical, parameter :: uniform = .true.
logical, parameter :: pst = .false.
logical, parameter :: ssh_a = .false.
logical, parameter :: ssh_b = .false.
logical, parameter :: abc = .false.
logical, parameter :: kitaev = .false.

!*You want to calculate dynamical figures?*!:
logical, parameter :: dynamics = .true. !calculation of dynamics
logical, parameter :: single = .false.  !single point calculation

!*You want to read the initial state from file?*!:
logical, parameter :: read_state = .false. !read from state, time needs to be higher than previous

!*What method do you wish to use to solve the Schrodinger eq.?*!:
logical, parameter :: integration = .false.
logical, parameter :: diagonalisation = .true.

!*Define presence of disorder*!:
logical, parameter :: random_J = .false. !Off-diagonal disorder
logical, parameter :: random_D = .false. !Diagonal disorder

!*Manage files*!:
logical, parameter :: output = .true.
logical, parameter :: files = .true.
logical, parameter :: graphical = .true.

!**************************************
!!Basic characteristics of the system *
!**************************************

integer, parameter :: N = 3           !Total number of sites
integer, parameter :: exno = 1         !Total number of excitations

integer, parameter :: numI = 1         !Total number of initial injected states
integer, dimension(numI) :: initialVec
!GO TO END OF FILE TO SPECIFY INITIAL VECTORS

integer, parameter :: branches = 1     !Number of branches, if linear set to 1.

!**********************
!!Coupling parameters *
!**********************

real(kind=dbl), parameter :: J_max = 1.0    !Maximum coupling in the middle
                                            !for PST. Used also when uniform
real(kind=dbl), parameter :: J_strong = 1.0 !Strong versus weak coupling for
real(kind=dbl), parameter :: J_weak = 0.1   !SSH-like schemes.

!************************************
!!Disorder and tolerance parameters *
!************************************

!*For average purposes*!:
integer, parameter        :: num_realisations = 1 !How many noise realisations (1 by default)

real(kind=dbl), parameter :: E_J = 1.0_dbl  !scale of the disorder on the couplings, in units of J_max
real(kind=dbl), parameter :: E_D = 0.10_dbl !scale of the disorder on the sites, in units of J_max

real(kind=dbl), parameter :: error=0.0001_dbl !allowed error for integration method and normalisation checks

!**********************
!!Dynamics parameters *
!**********************

integer, parameter :: steps = 50000
real(kind=dbl), parameter :: totalTime = 8
real(kind=dbl), parameter :: t_A = 500   !time for single point calculation (set single option)


!****************
!!Entanglements *
!****************

!choose your measure
logical, parameter :: eof = .true.
logical, parameter :: max_eof = .true. !calculates the maximum eof over
                                      !a time window = totalTime

!qubits to trace for EOF
integer, parameter :: Q1 = 1
integer, parameter :: Q2 = N


!****************
!!INITIAL STATE *
!****************

contains

subroutine initialState(initialVec)

integer, dimension(numI), intent(inout) :: initialVec

!Comment vectors not needed


initialVec(1) = 2 !Initial state 1
!initialVec(2) = 0  !Initial state 2
!initialVec(3) = 0  !Initial state 3
!initialVec(4) = 0  !Initial state 4

!!Keep adding vectors in this fashion and with this same numenclature
!!Normalisation is done automatically

!!Overall initial State = norm*[initialVec(1)+...+intialVec(numI)]

end subroutine initialState

END MODULE

