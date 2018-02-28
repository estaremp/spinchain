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

