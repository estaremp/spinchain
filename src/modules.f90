MODULE constants
implicit none

integer, parameter :: dbl = 8

real(kind=dbl), parameter :: pi     = 3.141592654_dbl

complex(kind=dbl), parameter :: im = cmplx(0._dbl,1._dbl,kind=dbl)

END MODULE

MODULE dependencies
implicit none

CONTAINS

    include 'BaseGenerator.f90'
    include 'Hamiltonian_linear.f90'
    include 'Hamiltonian_crosses.f90'
    include 'Dynamics.f90'

END MODULE
