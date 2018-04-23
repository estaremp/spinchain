!=========================================================================!
! Subroutine that computes the EOF of the system given an initial pair    !
! of sites and the reduced density matrix.                                !
!                                                                         !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   EOF => eof values at every time.                                      !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 31/03/2018                           !
!=========================================================================!

subroutine entanglement_of_formation(red_rho,eof_rho)

use constants
use parameters

complex(kind=dbl), dimension(4,4), intent(in) :: red_rho
real(kind=dbl), intent(inout) :: eof_rho


integer :: info
integer :: i, j, ii

real(kind=dbl) :: tau_rho
real(kind=dbl) :: x_rho
real(kind=dbl) :: log_rho

integer, dimension (1) :: eigvalsarr_rho

real(kind=dbl), dimension (4)      :: eigvals_rhob

!complex(kind=dbl), dimension (4,4) :: pauli
complex(kind=dbl), dimension (4,4) :: m_right
complex(kind=dbl), dimension (4,4) :: m_left
complex(kind=dbl), dimension (4,4) :: rho_conjg
complex(kind=dbl), dimension (4,4) :: rho_rho
complex(kind=dbl), dimension (4)   :: eigvals_rho
complex(kind=dbl), dimension (4)   :: eigvals_rho2
complex(kind=dbl), dimension (8)   :: work_rho
complex(kind=dbl), dimension (8)   :: rwork_rho

call calc_pauli(pauli)

rho_conjg = cmplx(0.0_dbl,0.0_dbl,kind=dbl)
m_right=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
m_left=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
rho_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)

!!Calculate CONJG(rho)
do i=1,4
    do j=1,4
        rho_conjg(i,j)=conjg(red_rho(i,j))
    enddo
enddo

!!Calculate m_left=rho*sigma
do i=1,4
    do j=1,4
        do ii=1,4
            m_left(i,j)=m_left(i,j)+(red_rho(i,ii)*pauli(ii,j))
        enddo
    enddo
enddo

!!Calculate m_right=m_left*CONJG(rho)
do i=1,4
    do j=1,4
        do ii=1,4
            m_right(i,j)=m_right(i,j)+(m_left(i,ii)*rho_conjg(ii,j))
        enddo
    enddo
enddo

!!Calculate rho_rho=m_right*sigma
do i=1,4
    do j=1,4
        do ii=1,4
            rho_rho(i,j)=rho_rho(i,j)+(m_right(i,ii)*pauli(ii,j))
        enddo
    enddo
enddo

eigvals_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
m_right=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
m_left=cmplx(0.0_dbl,0.0_dbl,kind=dbl)

!*********************************************************************************************
! LAPACK SUBROUTINE ©                                                                       !*
! ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the                           !*
!*  eigenvalues and, optionally, the left and/or right eigenvectors.                        !*
!*                                                                                          !*
call zgeev('N','N',4,rho_rho,size(rho_rho,1),eigvals_rho,m_left,4,m_right,&
&4,work_rho,size(work_rho,1),rwork_rho,info)
if(info/=0) stop 'Error in zgeev'                                                           !*
!*                                                                                          !*
!*********************************************************************************************


eigvals_rhob=dble(eigvals_rho) !Take the double precision real part

!Sort eigenvalues in decreasing order and take square root
do i=1,4
    eigvals_rho2(i)=sqrt(maxval(eigvals_rhob))
    eigvalsarr_rho=maxloc(eigvals_rhob)
    j=eigvalsarr_rho(1)
    eigvals_rhob(j)=0._dbl
enddo

tau_rho=eigvals_rho2(1)
do i=2,4
    tau_rho=tau_rho-eigvals_rho2(i)
enddo

if (tau_rho.le.0) then
    tau_rho=0._dbl
else
    tau_rho=tau_rho**2
endif



x_rho=(1._dbl+sqrt(1._dbl-tau_rho))*0.5_dbl


log_rho=1._dbl/log(2._dbl)
if ((1._dbl-x_rho).le.(error)) then
    eof_rho=0._dbl
else
    eof_rho=-x_rho*(log(x_rho)*log_rho)-(1._dbl-x_rho)*(log(1._dbl-x_rho)*log_rho)
endif

end subroutine
