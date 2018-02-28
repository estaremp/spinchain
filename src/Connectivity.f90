!=========================================================================!
! Subroutine that defines the connectivity of the system.                 !
!                                                                         !
! Available options:                                                      !
!   + Linear                                                              !
!   + Cross                                                               !
!        - Three branches                                                 !
!        - Four branches                                                  !
!        - Five branches                                                  !
!        - Six branches                                                   !
!   + More to be included... circular, lattices,...                       !
!                                                                         !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!    => Parameters that determine connectivity                            !
!           + hub                                                         !
!           + len branch
!-------------------------------------------------------------------------!
! Written by  Marta Estarellas, v0.1, 14/02/2016                          !
!=========================================================================!
subroutine couplings(Js)

use constants
use parameters

real(kind=dbl), dimension(N-1), intent(inout) :: Js
real(kind=dbl) :: J_0              !characteristic coupling (energy scale) for PST
integer :: i,j,k

SELECT CASE ()


end subroutine
