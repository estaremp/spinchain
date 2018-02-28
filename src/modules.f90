!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE constants
implicit none

!!Precision
integer, parameter :: sbl = 4
integer, parameter :: dbl = 8

!!Pi
real(kind=dbl), parameter :: pi = 3.141592654_dbl

!!Imaginary i
complex(kind=dbl), parameter :: im = cmplx(0._dbl,1._dbl,kind=dbl)

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  PARAMETERS - modify THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE parameters

use constants

!*********************
!!Initial parameters *
!*********************
integer, parameter :: N = 4
integer, parameter :: branches = 3
integer, parameter :: exno = 1

!**********************
!!Coupling parameters *
!**********************
real(kind=dbl), parameter :: J_max = 1.0_dbl
real(kind=dbl), parameter :: J_strong = 1.0_dbl
real(kind=dbl), parameter :: J_weak = 0.1_dbl

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
logical, parameter :: linear = .true.
logical, parameter :: branched = .false.
logical, parameter :: crossed = .false.
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


END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  SUBROUTINES DEPENDENCIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE dependencies
implicit none

CONTAINS

    include 'BaseGenerator.f90'
    include 'Hamiltonian_linear.f90'
    include 'Hamiltonian_crosses.f90'
    include 'Dynamics.f90'

END MODULE
