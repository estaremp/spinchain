subroutine eof(rho_mat2,eof_rho)

real(kind=dbl), parameter :: error  = 0.01_dbl !This is the allowed error

integer :: n=2 !size of the reduced density matrix
integer :: i,j,ii,v,w
integer :: status !for eigenvalues
integer, dimension (1)                                :: eigvalsarr_rho

real(kind=dbl) :: tau_rho
real(kind=dbl) :: x_rho
real(kind=dbl) :: eof_rho
real(kind=dbl) :: log_rho
real(kind=dbl) :: eof_rho

real(kind=dbl), dimension (qtrace2)                   :: eigvals_rhob
real(kind=dbl), dimension (qtrace2)                   :: eigvals_rhob


complex(kind=dbl), dimension (qtrace2,qtrace2)        :: rho_mat3
complex(kind=dbl), dimension (qtrace2,qtrace2)        :: vl_rho
complex(kind=dbl), dimension (qtrace2,qtrace2)        :: vr_rho
complex(kind=dbl), dimension (qtrace2,qtrace2)        :: rho_mat2b
complex(kind=dbl), dimension (qtrace2,qtrace2)        :: pauli
complex(kind=dbl), dimension (qtrace2)                :: eigvals_rho
complex(kind=dbl), dimension (qtrace2)                :: eigvals_rho2
complex(kind=dbl), dimension (2*qtrace2)              :: work_rho !This is for eigenvector computation
complex(kind=dbl), dimension (2*qtrace2)              :: rwork_rho



    pauli=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    pauli(1,4)=(-1._dbl,0._dbl)
    pauli(2,3)=(1._dbl,0._dbl)
    pauli(3,2)=(1._dbl,0._dbl)
    pauli(4,1)=(-1._dbl,0._dbl)

    rho_mat3=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    vl_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    vr_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)

    rho_mat2b = cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    do i=1,n
        do j=1,n
            rho_mat2b(i,j)=conjg(rho_mat2(i,j))
        enddo
    enddo

    do i=1,n
        do j=1,n
            do ii=1,n
                vl_rho(i,j)=vl_rho(i,j)+(rho_mat2(i,ii)*pauli(ii,j))
            enddo
        enddo
    enddo

    do i=1,n
        do j=1,n
            do ii=1,n
                vr_rho(i,j)=vr_rho(i,j)+(vl_rho(i,ii)*rho_mat2b(ii,j))
            enddo
        enddo
    enddo

    do i=1,n
        do j=1,n
            do ii=1,n
                rho_mat3(i,j)=rho_mat3(i,j)+(vr_rho(i,ii)*pauli(ii,j))
            enddo
        enddo
    enddo

    vl_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    vr_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    eigvals_rho=cmplx(0.0_dbl,0.0_dbl,kind=dbl)
    v=n
    w=n
    status=0

    !rho diagonalization
    call zgeev('N','N',qtrace2,rho_mat3,size(rho_mat3,1),eigvals_rho,vl_rho,v,vr_rho,&
    &w,work_rho,size(work_rho,1),rwork_rho,status)
    if(status/=0) stop 'Error in zgeev'

    eigvals_rhob=dble(eigvals_rho) !Take the double precision real part

    !Sort eigenvalues in decreasing order and take square root
    do i=1,n
        eigvals_rho2(i)=sqrt(maxval(eigvals_rhob))
        eigvalsarr_rho=maxloc(eigvals_rhob)
        j=eigvalsarr_rho(1)
        eigvals_rhob(j)=0._dbl
    enddo


    tau_rho=eigvals_rho2(1)
    do i=2,n
        tau_rho=tau_rho-eigvals_rho2(i)
    enddo

    if (tau_rho.le.0) then
        tau_rho=0._dbl
    else
        tau_rho=tau_rho**2
    endif

    x_rho=(1._dbl+sqrt(1._dbl-tau_rho))*0.5_dbl

    log_rho=1._dbl/log(2._dbl)
    if ((1._dbl-x_rho).le.(error*0.00000000000001_dbl)) then
        eof_rho=0._dbl
    else
        eof_rho=-x_rho*(log(x_rho)*log_rho)-(1._dbl-x_rho)*(log(1._dbl-x_rho)*log_rho)
    endif


end subroutine
