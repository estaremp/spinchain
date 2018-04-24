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
logical, parameter :: uniform = .false.
logical, parameter :: pst = .false.
logical, parameter :: ssh_a = .false.
logical, parameter :: ssh_b = .false.
logical, parameter :: abc = .true.
logical, parameter :: kitaev = .false.

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

integer, parameter :: N = 7           !Total number of sites
integer, parameter :: exno = 2         !Total number of excitations

integer, parameter :: numI = 1         !Total number of initial injected states
integer, dimension(numI) :: initialVec
!GO TO END OF FILE TO SPECIFY INITIAL VECTORS

integer, parameter :: branches = 3     !Number of branches

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

integer, parameter :: num_realisations = 1

real(kind=dbl), parameter :: E_J = 0.0_dbl !scale of the disorder on the couplings
real(kind=dbl), parameter :: E_D = 1.0_dbl !scale of the disorder on the sites

real(kind=dbl), parameter :: error=0.000001_dbl
!**********************
!!Dynamics parameters *
!**********************

integer, parameter :: steps = 5000
real(kind=dbl), parameter :: totalTime = 500

!****************
!!Entanglements *
!****************

!choose your measure
logical, parameter :: eof = .true.

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
!!Keep adding vectors in this fasion and with this same numenclature

end subroutine initialState

END MODULE

