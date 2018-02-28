!=========================================================================!
! Recursive subroutine that returns a matrix with all the possible        !
! vectors of the system depending on N for the required Hilbert space     !
! sector                                                                  !
!                                                                         !
! REMEMBER: this still needs to be programmed on its generalised way      !
! so it is able to output the entire matrix with all the subsectors       !
! together                                                                !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   HX => matrix with all the X sector vectors                            !
!-------------------------------------------------------------------------!
! Written by Pablo H. Sampedro and Marta Estarellas, v0.1, 14/02/2016     !
!=========================================================================!
recursive subroutine permutations(ex,nit,vec,N,Ninit,H,vectors,k)


integer, intent(in) :: ex, nit, N, Ninit, vectors
integer, dimension(N), intent(inout) :: vec
integer, dimension(vectors,N), intent(inout) :: H
integer :: i,j,k

if (nit<ex) then
    do i=Ninit, N-ex+nit
        vec(i)=1
        call permutations(ex,nit+1,vec,N,i+1,H,vectors,k)
        vec(i)=0
    enddo
else
    do i= Ninit, N
        vec(i)=1
        H(k,:)=vec
        k=k+1
        vec(i)=0
    enddo
endif

end subroutine
