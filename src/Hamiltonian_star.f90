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
subroutine build_hamiltonian_star(HT,Js,N,vectorstotal,hami,branches,limits,hub,len_branch)

use constants

!input variables
integer, intent(in) :: vectorstotal
integer, intent(in) :: N
integer, intent(in) :: branches

real(kind=dbl), dimension(N), intent(in) :: Js
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT

real(kind=dbl), dimension(vectorstotal,vectorstotal), intent(inout) :: hami

integer :: i,j,k
integer :: l,ll,m,b
integer, dimension(N) :: test

integer, intent(in) :: hub
integer, intent(in) :: len_branch
integer, dimension(branches-1), intent(in) :: limits

!*************************************!
!******+ BUILD HAMILTONIAN ***********!
!*************************************!

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
                    do b=3,branches
                    if (test(limits(b-2)+(len_branch-1)/2).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                    enddo
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
            !building consecutive branches
            do b=4,branches
            do k=(limits(b-2)+(len_branch-1)/2),limits(b-1)-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                endif
            enddo
            enddo
        endif
    enddo
enddo

end subroutine
