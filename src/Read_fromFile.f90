!=========================================================================!
! Subroutine that reads the initial state from                            !
! a state.data file and normalizes it                                     !
! REMEMBER: time                                                          !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   TODO                                                                  !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 10/05/2018                           !
!=========================================================================!

subroutine readFromFile(vectorstotal,initial_state,time)

!modules
use constants
use parameters

!inputs
integer, intent(in) :: vectorstotal

!outputs
real(kind=dbl), intent(out) :: time
complex(kind=dbl), dimension(vectorstotal), intent(out) :: initial_state

complex(kind=dbl), dimension(vectorstotal) :: c_i

!local variables
integer :: i,j,k,p

real(kind=dbl) :: norm

!open files
open (unit=60,file='final_state.data',status='old')
read(60,*)
read(60,*) time
do i=1,vectorstotal
    read (60,*) c_i(i)
enddo
close(60)

!********************
!!RENORMALIZATION ***
!********************

norm = 0._dbl
do i=1,vectorstotal
    norm = norm + dconjg(c_i(i))*(c_i(i))
enddo

do i=1,vectorstotal
    initial_state(i)=(1._dbl/(sqrt(norm)))*c_i(i)
enddo


end subroutine
