!=========================================================================!
! Subroutine that calculates the reduced density matrix by tracing out    !
! qubit Q1 and Q2.                                                        !
!                                                                         !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!  red_rho                                                                !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 07/04/2017                           !
!=========================================================================!

subroutine reduced_density_matrix(HT,vectorstotal,red_rho,c_i)

use constants
use parameters

integer, parameter :: Nrho = N-2

integer :: nit,Ninit,ex
integer, dimension(N) :: vec
integer :: vectors1ex = Nrho
integer :: vectors2ex = Nrho
integer :: vectors3ex = Nrho
integer :: vectorstotalrho

integer, allocatable, dimension(:,:)  :: Hrho1,Hrho2,Hrho3,HrhoT
integer, allocatable, dimension (:,:) :: HTrhoP
integer, allocatable, dimension (:,:) :: rho

integer, intent(in) :: vectorstotal
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT
complex(kind=dbl), dimension(vectorstotal), intent(in) :: c_i
complex(kind=dbl), dimension (4,4), intent(inout) :: red_rho

integer :: i,j,k,v,w,jj,vv
logical :: ok = .false.


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!TODO - GENERALIZE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!WRITE BASIS SET FOR THE TRACE
!Allocate matrices that will contain all the vectors:
if (exno==1) then
    vectorstotalrho = vectors1ex+1
else if (exno==2) then
    vectors2ex = (Nrho*(Nrho-1)/2)
    vectorstotalrho = vectors1ex+vectors2ex+1
else if (exno==3) then
    vectors2ex = (Nrho*(Nrho-1)/2)
    vectors3ex = (Nrho*(Nrho-1)*(N-2)/6)
    vectorstotalrho = vectors1ex+vectors2ex+vectors3ex+1
end if


allocate(Hrho1(vectors1ex,Nrho))
allocate(Hrho2(vectors2ex,Nrho))
allocate(Hrho3(vectors3ex,Nrho))
allocate(HrhoT(vectorstotalrho,Nrho))

Hrho1 = 0  !1ex subspace matrix
Hrho2 = 0  !2ex subspace matrix
Hrho3 = 0  !3ex subspace matrix
!... keep adding matrices
HrhoT = 0  !total vectors

!Create the subsectors matrices through a recursive call to Permutations
!First subsector (including ground state - all spins down):

do i=1,Nrho
    do j=1,Nrho
        if (i.eq.j) then
            Hrho1(i,j)=1
        endif
    enddo
enddo

HrhoT(2:,:) = Hrho1

!Second subsector (two excitations):

if (exno>1) then
nit=1
Ninit=1
vec=0
k=1
ex=2
call permutations(ex,nit,vec,Nrho,Ninit,Hrho2,vectors2ex,k)

HrhoT(vectors1ex+2:,:) = Hrho2
endif

!Third subsector (three excitations):

if (exno>2) then
nit=1
Ninit=1
vec=0
k=1
ex=3
call permutations(ex,nit,vec,Nrho,Ninit,Hrho3,vectors3ex,k)

HrhoT(vectors1ex+vectors2ex+2:,:) = Hrho3
endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allocate (HTrhoP(vectorstotal,Nrho))
allocate (rho(vectorstotalrho,4))

HTrhoP=0
rho=0
!reconstruction of the hami matrix taking Q1 and Q2 out
do i =1, vectorstotal
    w=1
    do j=1,Q1-1
        HTrhoP(i,w)=HT(i,j)
        w=w+1
    enddo
    do j=Q1+1,Q2-1
        HTrhoP(i,w)=HT(i,j)
        w=w+1
    enddo
    do j=Q2+1,N
        HTrhoP(i,w)=HT(i,j)
    enddo
enddo

!build rho
    do i=1,vectorstotalrho
        w=0
        do j=1,vectorstotal
            ok=.true.
            do v=1,Nrho
                if (HTrhoP(j,v).ne.HrhoT(i,v)) then
                ok=.false.
                exit
                endif
            enddo
        if (ok) then
            w=w+1
            rho(i,w)=j
            if (w.gt.4) print*, 'There are supposedly more than', 4,'vectors with the same middle...'
        endif
        enddo
    enddo

    do i=1,vectorstotalrho
        do vv=1,4
            if (rho(i,vv).ne.0) then
                do jj=vv,4
                if (rho(i,jj).ne.0) then
                    red_rho(vv,jj)=red_rho(vv,jj)+c_i(rho(i,vv))*conjg(c_i(rho(i,jj)))
                endif
                enddo
            endif
        enddo
    enddo

    ! Fill in other triangle of rho_mat2
    do k=1,4
        do j=k+1,4
        red_rho(j,k)=conjg(red_rho(k,j))
        end do
    end do

deallocate(HTrhoP)
deallocate(rho)
deallocate(Hrho1)
deallocate(Hrho2)
deallocate(Hrho3)
deallocate(HrhoT)

end subroutine
