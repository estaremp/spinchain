MODULE parameters

use constants

!*********************
!!Initial parameters *
!*********************
integer, parameter :: N = 4            !Total number of sites
integer, parameter :: exno = 1         !Total number of excitations

integer, parameter :: numI = 1         !Total number of initial injected states
integer, parameter :: initialVec1 = 3  !Initial state 1
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
real(kind=dbl), parameter :: r_noise = 0.0_dbl
real(kind=dbl), parameter :: error = 0.00000001_dbl

!**********************
!!Dynamics parameters *
!**********************
integer, parameter :: steps = 50000
real(kind=dbl), parameter :: totalTime = 100

!***********************
!!Other specifications *
!***********************

!*Define chain topology*!:
logical, parameter :: linear = .false.
logical, parameter :: star = .true.
logical, parameter :: lattice = .false.
!logical, parameter :: crossed_three, crossed_four, crossed_five, crossed_six
logical, parameter :: squared = .false.

!*Define couplings scheme*!:
logical, parameter :: uniform = .true.
logical, parameter :: pst = .false.
logical, parameter :: ssh_a = .false.
logical, parameter :: ssh_b = .false.
logical, parameter :: abc = .false.
logical, parameter :: kitaev = .false.

!*Define presence of disorder*!:
logical, parameter :: random = .false.

!*Manage files*!:
logical, parameter :: output = .true.
logical, parameter :: graphical = .true.
logical, parameter :: files = .true.
logical, parameter :: clean = .true.


END MODULE

