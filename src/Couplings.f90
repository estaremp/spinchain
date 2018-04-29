!=========================================================================!
! Sbroutine that returns an array with all the couplings distributed      !
! according the chosen connectivity (linear, cross, branched...)          !
!                                                                         !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   Js => array with all the couplings                                    !
!-------------------------------------------------------------------------!
! Written by  Marta Estarellas, v0.1, 14/02/2016                          !
!=========================================================================!
subroutine couplings(Js,len_branch,hub,limits)

use constants
use parameters


real(kind=dbl), dimension(N-1), intent(inout) :: Js
integer, intent(in) :: hub
integer, intent(in) :: len_branch
integer, dimension(branches-1), intent(in) :: limits

real(kind=dbl) :: J_0              !characteristic coupling (energy scale) for PST
integer :: i,j,k

!Uniform:

if (uniform) then
    do i=1,N-1
        Js(i) = J_max
    enddo
endif

!PST

if (pst) then
    J_0 = (2._dbl*J_max)/(N*sqrt(1._dbl-(1._dbl/(N**2))))
    do i=1,N-1
        Js(i) = J_0*sqrt(float(i*(N-i)))
    enddo
endif

!SSH


!type (a)
if (ssh_a) then
    do i=1,N-1
        if (i<=(N/2)) then !Induce alpha-configuration 
            if (MOD(i,2)==0) then
                Js(i)=J_weak
            else
                Js(i)=J_strong
            endif
        endif

        if (i>(N/2)) then !Induce beta-configuration
            if (MOD(i,2)==0) then
                Js(i)=J_strong
            else
                Js(i)=J_weak
            endif
        endif
    enddo
endif

!type (b)
if (ssh_b) then
    do i=1,N-1
        if (i<=(N/2)) then !Induce alpha-configuration
            if (MOD(i,2)==0) then
                Js(i)=J_strong
            else
                Js(i)=J_weak
            endif
        endif

        if (i>(N/2)) then !Induce beta-configuration
            if (MOD(i,2)==0) then
                Js(i)=J_weak
            else
                Js(i)=J_strong
            endif
        endif
    enddo
endif

!abc
if (abc) then
    if (linear) then
    do i=1,N-1
        if (i<=(N/2)) then
            if (MOD(i,2)==0) then
                Js(i)=J_strong
            else
                Js(i)=J_weak
            endif
        else
            if (MOD(i,2)==0) then
                Js(i)=J_weak
            else
                Js(i)=J_strong
            endif
        endif
    enddo
    endif

if (star) then
Js=0
do i=1,N-1
    if (i<len_branch) then
        if (i<(hub)) then
            if (MOD(i,2)==0) then
                Js(i)=J_strong
            else
                Js(i)=J_weak
            endif
        else
            if (MOD(i,2)==0) then
                Js(i)=J_weak
            else
                Js(i)=J_strong
            endif
        endif
    endif
    if (i>=len_branch) then
        if (i<(limits(1)+((len_branch-1)/2))) then
            if (MOD(i,2)==0) then
                Js(i)=J_strong
            else
                Js(i)=J_weak
            endif
        else
            if (MOD(i,2)==0) then
                Js(i)=J_weak
            else
                Js(i)=J_strong
            endif
        endif
    endif
enddo
endif

endif

!Kitaev
if (kitaev) then
    do i=1,N-1
        if (MOD(i,2)==0) then
            Js(i)=J_strong
        else
            Js(i)=J_weak
        endif
    enddo
endif

end subroutine
