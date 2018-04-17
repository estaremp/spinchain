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
logical, parameter :: pst = .true.
logical, parameter :: ssh_a = .false.
logical, parameter :: ssh_b = .false.
logical, parameter :: abc = .false.
logical, parameter :: kitaev = .false.

!*Define presence of disorder*!:
logical, parameter :: random_J = .false. !Off-diagonal disorder
logical, parameter :: random_D = .true. !Diagonal disorder

!*Manage files*!:
logical, parameter :: output = .true.
logical, parameter :: files = .true.
logical, parameter :: graphical = .true.

!**************************************
!!Basic characteristics of the system *
!**************************************

integer, parameter :: N = 10           !Total number of sites
integer, parameter :: exno = 1         !Total number of excitations

integer, parameter :: numI = 1         !Total number of initial injected states
integer, parameter :: initialVec1 = 2  !Initial state 1
integer, parameter :: initialVec2 = 0  !Initial state 2
integer, parameter :: initialVec3 = 0  !Initial state 3
integer, parameter :: initialVec4 = 0  !Initial state 4

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
real(kind=dbl), parameter :: totalTime = 20

!****************
!!Entanglements *
!****************

!qubits to trace for EOF
integer, parameter :: Q1 = 1
integer, parameter :: Q2 = N


END MODULE

