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

complex(kind=dbl), dimension (4,4) :: pauli = cmplx(0.0_dbl,0.0_dbl,kind=dbl)


contains
    subroutine calc_pauli(pauli)

        complex(kind=dbl), dimension (4,4), intent(inout) :: pauli
        !!Pauli mattrices (for EOF)
        pauli(1,4)=(-1._dbl,0._dbl)
        pauli(2,3)=(1._dbl,0._dbl)
        pauli(3,2)=(1._dbl,0._dbl)
        pauli(4,1)=(-1._dbl,0._dbl)
    end subroutine calc_pauli

END MODULE

