!=========================================================================!
! Subroutine that computes the Entanglement of formation between two      !
! specific qubits passed as argument Q1 and Q2                            !
!                                                                         !
! REMEMBER: each column of the returned matrix corresponds to the         !
! fidelity of the initial state against the vector with the same number   !
! than that column. First column is time                                  !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!  EOF => matrix with the dynamics of the system                          !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 07/04/2017                           !
!=========================================================================!

subroutine entanglement_formation(HT,N,vectorstotal,Q1,Q2)



integer, intent(in) :: Q1, Q2

integer :: i,j,k,v,w


!reconstruction of the hami matrix taking Q1 and Q2 out
do i =1, vectorstotal
    w=1
    do j=1,Q1-1
        HT_rho(i,w)=HT(i,j)
        w=w+1
    enddo
    do j=Q1+1,Q2-1
        HT_rho(i,w)=HT(i,j)
        w=w+1
    enddo
    do j=Q2+1,N
        HT_rho(i,w)=HT(i,j)
    enddo
enddo

!build rho
do i=1,vectors_rho
    w=0
    do j=1,vectors






end subroutine
