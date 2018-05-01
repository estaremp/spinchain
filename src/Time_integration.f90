!=========================================================================!
! Subroutine that computes the dynamics of the system given a initial     !
! vector through a time integration method                                !
!                                                                         !
! REMEMBER: each column of the returned matrix corresponds to the         !
! fidelity of the initial state against the vector with the same number   !
! than that column. First column is time                                  !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   Fidelity => matrix with the dynamics of the system                    !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 31/03/2017                           !
!=========================================================================!

subroutine time_integration(HT,hami,vectorstotal,c_i)

!modules
use constants
use parameters

!inputs
integer, intent(in) :: vectorstotal
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT

real(kind=dbl), dimension(vectorstotal,vectorstotal), intent(in) :: hami
complex(kind=dbl), dimension(vectorstotal), intent(inout) :: c_i

!local variables
integer :: i,j,k,p

real(kind=dbl) :: norm
real(kind=dbl) :: sum
real(kind=dbl) :: Dt
real(kind=dbl) :: time

!arrays
real(kind=dbl), dimension(vectorstotal) :: fidelity
real(kind=dbl), dimension(N) :: siteProb

complex(kind=dbl), dimension(vectorstotal) :: Ac

!initialize fidelity
fidelity=0._dbl

!open files
open (unit=44,file='dynamics.data',status='unknown')
open (unit=45,file='exmap.data',status='unknown')

!writte comments on the oytputs
write(44,*) '#DYNAMICS. TIME (1st COL) AND FIDELITIES AGAINST ALL THE BASIS VECTORS ORDERED BY INDEX NUMBER'
write(45,*) '#DYNAMICS. TIME (1st COL) AND SITE OCCUPATION PROBABILITES ORDERED BY SITE'

Dt = real(totalTime/steps) !stepsize

!normalisation factor
norm=(1._dbl/sqrt(float(numI)))

!get initial vector from PARAMETERS file inputs
call initialState(initialVec)

!initial stat
do i=1,numI
    c_i(initialVec(i))=norm
enddo

!start iterations
time = 0._dbl

do while (time<=totalTime)

    !sum coef
    Ac = (0.0_dbl,0.0_dbl)
    do j=1,vectorstotal
        do k=1,vectorstotal
            Ac(j)=Ac(j)+c_i(k)*hami(j,k)
        enddo
    enddo

    !evolve coef
    do p=1,vectorstotal
        c_i(p) = c_i(p) - im*Dt*Ac(p)
    enddo

    !evaluate error
    sum=0
    do i=1,vectorstotal
        sum=sum+(abs(c_i(i))**2)
    enddo
    if ((sum-1.0_dbl).gt.error) then
        print*, "Error too big at time ",time, " and step_size ", Dt
        STOP
    endif

    !Fidelities
    do i=1,vectorstotal
        fidelity(i) = (abs(c_i(i)))**2
    enddo

    siteProb=0
    do i=1,vectorstotal
        do k=1,N
            if (HT(i,k)==1) then
                siteProb(k) = siteProb(k) + fidelity(i)
            endif
        enddo
    enddo


    if (pst) then
        write(44,*) time*J_max, fidelity
        write(45,*) time*J_max, siteProb
    else
        write(44,*) time*J_max, fidelity
        write(45,*) time*J_max, siteProb
    endif

    time = time + Dt
enddo

!close files
close(44)
close(45)

end subroutine
