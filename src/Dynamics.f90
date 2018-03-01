!=========================================================================!
! Subroutine that computes the dynamics of the system given a initial     !
! vector.                                                                 !
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

subroutine injection_dynamics(HT,hami,eigvals,vectorstotal,initialVec,norm,c_i)

use constants
use parameters

!inputs
integer, intent(in) :: vectorstotal
integer, intent(in) :: initialVec
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT

real(kind=dbl), intent(in) :: norm

real(kind=dbl), dimension(vectorstotal), intent(in) :: eigvals
complex(kind=dbl), dimension(vectorstotal,vectorstotal), intent(in) :: hami
complex(kind=dbl), dimension(vectorstotal), intent(inout) :: c_i

integer :: i,j,k

real(kind=dbl) :: step_size, time
real(kind=dbl), allocatable, dimension(:) :: fidelity
real(kind=dbl), dimension(N) :: siteProb
complex(kind=dbl), allocatable, dimension(:) :: a_m
complex(kind=dbl), allocatable, dimension(:,:) :: Y_o, Y_t
complex(kind=dbl) :: sum_vec

allocate(a_m(vectorstotal))
allocate(fidelity(vectorstotal))
allocate(Y_o(vectorstotal,vectorstotal))
allocate(Y_t(vectorstotal,vectorstotal))

fidelity=0._dbl

!open files
open (unit=44,file='dynamics.data',status='unknown')
open (unit=45,file='exmap.data',status='unknown')

!TODO: to be converted in function
!Define |¥(0)> = \sum{a_m|m>} and a_m = <¥(0)|m>
!being <¥(0)|=norm*(<initialVec1|+<initialVec2|+...)
do i=1,vectorstotal
    a_m(i) = norm*dconjg(hami(initialVec,i))
enddo

step_size = totalTime/steps
time = 0._dbl

do while (time<=totalTime)

    do i=1,vectorstotal
        do j=1,vectorstotal
            Y_o(i,j) = hami(i,j)*a_m(j)
            Y_t(i,j) = Y_o(i,j)*exp(-1._dbl*im*time*eigvals(j))
        enddo
    enddo

    !|¥(t)> = \sum{c_i|i>}
    c_i=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
    do i=1,vectorstotal
        sum_vec = cmplx(0._dbl,0._dbl,kind=dbl)
        do j=1,vectorstotal
            sum_vec = sum_vec + Y_t(i,j)
        enddo
        c_i(i) = sum_vec
    enddo

    siteProb=0
    do i=1,vectorstotal
        do k=1,N
            if (HT(i,k)==1) then
                siteProb(k) = siteProb(k) + fidelity(i)
            endif
        enddo
    enddo

    !Fidelity
    do i=1,vectorstotal
        fidelity(i) = (abs(c_i(i)))**2
    enddo

    if (pst) then
        write(44,*) time*J_max, fidelity
        write(45,*) time*J_max, siteProb
    else
        write(44,*) time*J_max, fidelity
        write(45,*) time*J_max, siteProb
    endif

    time = time + step_size

enddo

close(44)
close(45)

deallocate(a_m)
deallocate(Y_o)
deallocate(Y_t)

end subroutine
