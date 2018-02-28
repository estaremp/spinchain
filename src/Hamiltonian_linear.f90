!=========================================================================!
! Subroutine that builds up the Hamiltonian depending on the initial      !
! conditions for the chain:                                               !
!  - LINEAR                                                               !
!-------------------------------------------------------------------------!
! Input values:                                                           !
!   HT : matrix with total hilbert space used                             !
!   Js : array with the coupling pattern                                  !
!   N  : number of sites                                                  !
!   vectorstotal : total number of basis vectors                          !
!   hami : hamiltonian matrix to be returned                              !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   hami => matrix with the build hamiltonian                             !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 05/09/2017                           !
!=========================================================================!

subroutine build_hamiltonian_linear(HT,Js,N,vectorstotal,hami)

use constants

integer, intent(in) :: vectorstotal
integer, intent(in) :: N
real(kind=dbl), dimension(N), intent(in) :: Js
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT
real(kind=dbl), dimension(vectorstotal,vectorstotal), intent(inout) :: hami

integer :: i,j,k
integer :: l,ll,m
integer, dimension(N) :: test

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
            !building the hamiltonian
            do k=1,N-1
                if ((test(k)).eq.1) then
                    if (test(k+1).eq.1) then
                        hami(j,i)=Js(k)
                    endif
                endif
            enddo
        endif
    enddo
enddo

end subroutine
