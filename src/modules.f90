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

MODULE parameters

use constants

!!Initial parameters
integer, parameter :: N = 4
integer, parameter :: branches = 3
integer, parameter :: exno = 1

!!Coupling parameters
real(kind=dbl), parameter :: J_max = 1.0_dbl
real(kind=dbl), parameter :: J_strong = 1.0_dbl
real(kind=dbl), parameter :: J_weak = 0.1_dbl

!!Disorder and tolerance parameters
real(kind=dbl), parameter :: r_noise = 0.0_dbl
real(kind=dbl), parameter :: error = 0.00000001_dbl

!!Dynamics parameters
integer, parameter :: steps = 50000
real(kind=dbl), parameter :: totalTime = 100


END MODULE

MODULE dependencies
implicit none

CONTAINS

    include 'BaseGenerator.f90'
    include 'Hamiltonian_linear.f90'
    include 'Hamiltonian_crosses.f90'
    include 'Dynamics.f90'

END MODULE
