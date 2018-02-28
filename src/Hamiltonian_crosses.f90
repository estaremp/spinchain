!=========================================================================!
! Subroutine that builds up the Hamiltonian depending on the initial      !
! conditions for the chain:                                               !
!  - CROSSES                                                              !
!-------------------------------------------------------------------------!
! Input values:                                                           !
!   HT : matrix with total hilbert space used                             !
!   Js : array with the coupling pattern                                  !
!   N  : number of sites                                                  !
!   vectorstotal : total number of basis vectors                          !
!   hami : hamiltonian matrix to be returned                              !
!   branches : number of brances
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   hami => matrix with the build hamiltonian                             !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 04/04/2017                           !
!=========================================================================!
subroutine build_hamiltonian_crosses(HT,Js,N,vectorstotal,hami,branches)

use constants

!input variables
integer, intent(in) :: vectorstotal
integer, intent(in) :: N
integer, intent(in) :: branches

real(kind=dbl), dimension(N), intent(in) :: Js
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT

real(kind=dbl), dimension(vectorstotal,vectorstotal), intent(inout) :: hami

integer :: i,j,k
integer :: l,ll,m
integer :: hub, len_branch
integer, dimension(N) :: test
integer, dimension(branches-1) :: limits

!initialize
limits = 0

!*************************************!
!******+ DEFINE CONNECTIVITY *********!
!*************************************!

len_branch = ((2*(N - 1)/branches)+1)
hub = ((len_branch-1)/2) + 1

limits(1) = len_branch
limits(2) = len_branch + 1
do i=3,branches-1
    limits(i) = limits(i-1) + (len_branch-2)
enddo

hami=0.0_dbl
do i=1,vectorstotal
    do j=1,vectorstotal
        l=0
        ll=0
        m=0
        test=0
        do k=1,N
            if (HT(i,k).eq.1) then
                l=l+1 !number of 1s
            endif
            if (HT(j,k).eq.1) then
                ll=ll+1 !number of 1s
            endif
        enddo
        if (l.ne.ll) cycle !if we are not in the same hilbert subspace, cycle
        do k=1,N
            if (HT(i,k).eq.1) then
                if (HT(j,k).ne.1) then
                    m=m+1 !number of 1's moved
                endif
            endif
        enddo
        !cutre salchichero method (aka IFS party) - add two vectors to compare that only moved 1ex, if the new vector
        !has two consecutive 1s (i.e. 10100 + 10010 = 20110) the hopping term betwen both 1s goes to the Hamiltonian
        !yes, I know
        if (m==1) then
            do k=1,N
                test(k)=HT(i,k)+HT(j,k)
            enddo
            !building first linear branch (horitzontal branch)
            do k=1,limits(1)-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                endif
                !check for connectivity between hub and first sites of upper/lower branches
                if (k==hub) then
                if (test(hub).eq.1) then
                    if (test(limits(1)+(len_branch-1)/2).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                    if ((branches==4).or.(branches==5).or.(branches==6)) then
                    if (test(limits(2)+(len_branch-1)/2).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                    if ((branches==5).or.(branches==6)) then
                    if (test(limits(3)+(len_branch-1)/2).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                    endif
                    if (branches==6) then
                    if (test(limits(4)+(len_branch-1)/2).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                    endif
                    endif
                endif
                endif
            enddo
            !building second branch
            do k=limits(2),(limits(1)+(len_branch-1)/2)-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k-1)
                    endif
                endif
            enddo
            !building third branch
            if ((branches==4).or.(branches==5).or.(branches==6)) then
            do k=(limits(2)+(len_branch-1)/2),limits(3)-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                endif
            enddo
            endif
            !building forth branch
            if ((branches==5).or.(branches==6)) then
            do k=(limits(3)+(len_branch-1)/2),limits(4)-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                endif
            enddo
            endif
            !building fifth branch
            if (branches==6) then
                do k=(limits(4)+(len_branch-1)/2),limits(5)-1
                    if ((test(k)).eq.1) then
                        if (test(k+1).eq.1) then
                            hami(j,i)=Js(k)
                        endif
                    endif
                enddo
            endif
        endif
    enddo
enddo


end subroutine
