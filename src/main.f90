module noisemod

integer, parameter :: dubp=kind(1.0d0)

! For random numbers
integer, parameter :: bit_32=kind(4)
integer(bit_32) :: ix,iy

real(dubp) :: Noise_SDev=0.5_dubp ! For Normal distribution

end module

program spinchain

!!load subroutines
use dependencies
!!load constants
use constants
!!load initial parameters
use parameters

use noisemod

implicit none

!!This program solves an XY spin chain problem. It diagonalizes the hamiltonian obtaining
!!for chains or networks of different geometries and distributions, and obtains its
!!eigenvectors and eigenvalues as well as evolves it dynamically with a defined set of
!!initial conditions. Properties such Fidelity, Entropy and Entanglement of Formation
!!can be computed by using the relevant subroutines
!a change

!NOTES:


!***************************************************!
!******************  VARIABLES *********************!
!***************************************************!

!integers
integer :: i,j,k,v,w   !loop dummies
integer :: nit,Ninit,ex      !subroutine Permutations variables
integer :: seed

integer :: vectors1ex = N       !Initially set to N, reallocate later if needed (E.I.)
integer :: vectors2ex = N       !Initially set to N, reallocate later if needed (E.I.)
integer :: vectors3ex = N       !Initially set to N, reallocate later if needed (E.I.)
integer :: vectorstotal         !Sum of all the vectors
integer :: info, liwork  !Info in lapack subroutines

integer,dimension(8) :: values !array with date
integer, dimension(N) :: vec
integer, allocatable, dimension (:) :: iwork

integer, allocatable, dimension(:,:) :: H1,H2,H3,HT !Hilbert subspaces matrices

!floats
real(kind=dbl) :: norm,normal,orto !normaliztion constant
real(kind=dbl) :: r !random number

real(kind=dbl), dimension(N-1) :: Js = 0.0_dbl

real(kind=dbl), allocatable, dimension(:) :: eigvals
real(kind=dbl), allocatable, dimension(:) :: rwork

real(kind=dbl), allocatable, dimension(:) :: Noise

!complex
complex(kind=dbl), allocatable, dimension(:) :: work
complex(kind=dbl), allocatable, dimension(:) :: c_i

complex(kind=dbl), allocatable, dimension(:,:) :: red_rho
complex(kind=dbl), allocatable, dimension(:,:) :: hamiD !diagonalized Hamiltonian

!matrices
real(kind=dbl), allocatable, dimension(:,:) :: hami !Hamiltonian

!strings
character :: a
character(len=32) :: tmp
character(len=500) :: fmt1,fmt2 !format descriptors

!random number generator
real(dubp), external :: algor_uniform_random


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! START PROGRAM AND WRITE OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (output) then
!retrieve date
call date_and_time(VALUES=values)
101 FORMAT (1X,59("*"))
102 FORMAT (1X,16("*")," SPIN CHAIN PROGRAM OUTPUT ",16("*"))
103 FORMAT (20X,I2,"/",I2.1,"/",I4,2X,I2,":",I2)
104 FORMAT (1X,59("-"))
open (unit=40,file='spinchain.out',status='replace')
write(40,101)
write(40,102)
write(40,101)
write(40,*) '           © Marta P. Estarellas, 27/07/2016              '
write(40,*) '                   University of York                     '
write(40,103) values(3),values(2),values(1),values(5),values(6)
write(40,104)
write(40,*)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINING THE DESIRED TYPE OF CHAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**********************************************
!this is done in the module called PARAMETERS
!you should ONLY modify that module to set the
!conditions and structure of the chain.
!**********************************************

write(*,*) '>> Defining System'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIAL CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (linear) then
if (ssh_a.or.ssh_b) then
if (MOD(N-1,4)/=0) then
STOP 'ERROR: for type (a) ssh chain N needs to be odd and N-1 needs to be divisible by 4.'
endif
endif

if (abc) then
if (MOD(N-3,4)/=0) then
STOP 'ERROR: for type ABC chain N needs to be odd and N-3 needs to be divisible by 4.'
endif
endif

if (kitaev) then
if (MOD(N,2)/=0) then
STOP 'ERROR: for a kitaev chain N needs to be even.'
endif
endif
endif

if (star) then
if (branches==3) then
if (MOD((N-1),3)/=0) then
STOP 'ERROR: Triple branched networks need to have EVEN number of sites and (N-1) needs to be divisible by 3.'
endif
endif

if (branches==4) then
if (MOD((N-1),4)/=0) then
STOP 'ERROR: Four branched networks need to have ODD number of sites and (N-1) needs to be divisible by 4.'
endif
endif

if (branches==5) then
if (MOD((N-1),5)/=0) then
STOP 'ERROR: Five branched networks need to have EVEN number of sites and (N-1) needs to be divisible by 5.'
endif
endif

if (branches==6) then
if (MOD((N-1),6)/=0) then
STOP 'ERROR: Six branched networks need to have ODD number of sites and (N-1) needs to be divisible by 6.'
endif
endif
endif

write(*,*) '>> Initial checks'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINING BASIS VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!TODO - GENERALIZE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Calculate number of vectors for each excitation N!/exno!(N-exno)! subspace and the total number
!this is done progressively, sector by sector for sake of efficiency:

if (exno==1) then
    vectorstotal = vectors1ex+1
else if (exno==2) then
    vectors2ex = (N*(N-1)/2)
    vectorstotal = vectors1ex+vectors2ex+1
else if (exno==3) then
    vectors2ex = (N*(N-1)/2)
    vectors3ex = (N*(N-1)*(N-2)/6)
    vectorstotal = vectors1ex+vectors2ex+vectors3ex+1
end if

!Allocate matrices that will contain all the vectors:

allocate(H1(N,N))
allocate(H2(vectors2ex,N))
allocate(H3(vectors3ex,N))
allocate(HT(vectorstotal,N))

H1 = 0  !1ex subspace matrix
H2 = 0  !2ex subspace matrix
H3 = 0  !3ex subspace matrix
!... keep adding matrices
HT = 0  !total vectors

!GENERALIZE THIS
!Create the subsectors matrices through a recursive call to Permutations
!First subsector (including ground state - all spins down):

do i=1,N
    do j=1,N
        if (i.eq.j) then
            H1(i,j)=1
        endif
    enddo
enddo

HT(2:,:) = H1

!Second subsector (two excitations):

if (exno>1) then
    nit=1
    Ninit=1
    vec=0
    k=1
    ex=2
    call permutations(ex,nit,vec,N,Ninit,H2,vectors2ex,k)
    HT(vectors1ex+2:,:) = H2
endif

!Third subsector (three excitations):

if (exno>2) then
    nit=1
    Ninit=1
    vec=0
    k=1
    ex=3
    call permutations(ex,nit,vec,N,Ninit,H3,vectors3ex,k)
    HT(vectors1ex+vectors2ex+2:,:) = H3
endif

!**(NOTE: Add extra subsectors in the same fashion if needed)**

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Stdout vectors martix
if (output) then
    201 FORMAT (/A)
    write(40,FMT=201) 'BASIS VECTORS:'
    do i=1,vectorstotal
        write(40,*) i,'-->',(HT(i,j),j=1,N)
    enddo
endif

write(*,*) '>> Basis vectors defined'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIAL STATE NORMALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) '>> Defining initial injection'

!normalization factor dependenig
!on the number of initial injections

norm=(1._dbl/sqrt(float(numI)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINE CONNECTIVITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!call connectivity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINE COUPLING PATTERN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call couplings(Js)

!Stdout coupling pattern
301 FORMAT ("(spin",I3,")-(spin",I3,") -->",F6.2)
if (output) then
    write(40,FMT=201) 'COUPLING PATTERN:'
    do i=1,N-1
        write(40,FMT=301) i, i+1, Js(i)
    enddo
endif

write(*,*) '>> Coupling pattern defined'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ADD PERTURBATION FACTORS TO THE COUPLINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

seed = 0
call noise_sub1 (seed)

if (random_J) then
    do i=1,N-1
        r=algor_uniform_random()
        Js(i)=Js(i)+(r*E_J*J_max)
    enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! BUILD HAMILTONIAN IN THE SPIN BASIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(hami(vectorstotal,vectorstotal))

hami=0.0_dbl
if (linear) then
    call build_hamiltonian_linear(HT,Js,N,vectorstotal,hami)
else if (star) then
    call build_hamiltonian_star(HT,Js,N,vectorstotal,hami,branches)
endif

!Stdout Hamiltonian
if (output) then
    write(40,FMT=201) 'HAMILTONIAN MATRIX:'
    do i=1,vectorstotal
        write(40,*) (hami(i,j),j=1,vectorstotal)
    enddo
endif

write(*,*) '>> Hamiltonian Build'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ADD PERTURBATION FACTORS TO THE DIAGONAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Noise(vectorstotal))

Noise=0
if (random_D) then
    do i=2,vectorstotal
    w = 0
        do j=1,N
            if (HT(i,j).eq.1) then
                w=w+1 !w counts the number of excitations in a vector
            endif
        enddo

    if (w.eq.1) then
        r=algor_uniform_random()
        Noise(i)=(E_D*r*J_max)
    else
        do j=2,N+1
            do k=1,N
                if ((HT(i,k).eq.1).and.(HT(i,k).eq.HT(j,k))) then
                    Noise(i)=Noise(i)+Noise(j)
                endif
            enddo
        enddo
    endif
    hami(i,i) = hami(i,i) + Noise(i)
    enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TRANSLATE THE HAMILTONIAN IN THE MJ BASIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DIAGONALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(hamiD(vectorstotal,vectorstotal))
allocate(eigvals(vectorstotal))
allocate(rwork((2*(vectorstotal**2))+5*vectorstotal+1))
allocate(work((vectorstotal**2)+2*vectorstotal))

liwork=5*vectorstotal+3

allocate(iwork(liwork))

hamiD=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
do i=1,vectorstotal
    do j=1,vectorstotal
        hamiD(i,j)=(hami(i,j))
    enddo
enddo

!save Hamiltonian matrix
if (files) then
    open(unit=89,file='hami.data',status='unknown')
    do i=1,vectorstotal
        write(89,*) (hami(i,j),j=1,vectorstotal)
    enddo
    close(89)
endif


!***********************************************************************************************************************
! LAPACK SUBROUTINE ©                                                                                                 !*
! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix                          !*
!call zheev('V','U',vectorstotal,hamiD,size(hamiD,1),eigvals,work,size(work,1),rwork,info)                            !*
                                                                                                                      !*
call zheevd('V','U',vectorstotal,hamiD,size(hamiD,1),eigvals,work,size(work,1),rwork,size(rwork,1),iwork,liwork,info) !*
if(info/=0) stop 'ERROR in ZHEEV diagonalization'                                                                     !*
                                                                                                                      !*
!***********************************************************************************************************************

!check normalisation eigenvectors
do i=1,vectorstotal
    normal=0.
    do j=1,vectorstotal
        normal=normal+abs(hamiD(i,j))**2
    enddo
    if (abs(1.-normal)>=error) then
        print*, 'ERROR: your eigenvectors are not well normalized'
        STOP
    endif
enddo

!check eigenvectors orthogonality
do v=1,vectorstotal
    do i=1,vectorstotal
        orto=0.
        do j=1,vectorstotal
            orto=orto+real(hamiD(i,j)*dconjg(hamiD(v,j)))
        enddo
        if ((orto>error).and.(v/=i)) then
            print*, 'ERROR: your eigenvectors are not orthogonal'
            STOP
        endif
    enddo
enddo

!!Stdout Eigenvalues
if (output) then

    !set formats
    write(tmp,'(i3.1)') vectorstotal
    fmt1='(1X,i3.1,1X,'//tmp//'("(",f7.3,f7.3,")"))'

    !print the headins
    !if N is small,
    !otherwise the
    !file gets too
    !messy
    if (N.le.20) then
    fmt2='(6X,'
        do i=1,vectorstotal
            write(tmp,'(i3.1)') i
            fmt2=trim(fmt2)//'"Eigenvector'//trim(adjustl(tmp))//':",3X,'
        enddo

        fmt2=trim(fmt2)//")"
    endif

        !Eigenvalues
        write(40,FMT=201) 'EIGENVALUES:'
        do i=1,vectorstotal
            write(40,*) eigvals(i)
        enddo

        !Eigenvectors
        write(40,FMT=201) 'EIGENVECTORS'

        if (N.le.20) then
        write(40,fmt2)
        endif
        do i=1,vectorstotal
            write(40,fmt1) i ,(hamiD(i,:))
        enddo

endif

!!Save data in files
if (files) then

    open (unit=41,file='coefficients.data',status='unknown')
    open (unit=42,file='probabilities.data',status='unknown')
    open (unit=43,file='eigenvalues.data',status='unknown')


    do i=1,vectorstotal
        write(41,*) real(hamiD(i,:))
        write(42,*) (abs(dconjg(hamiD(i,:))*(hamiD(i,:))))
        if (eigvals(i)==0._dbl) cycle
        write(43,*) eigvals(i)
    enddo

    close(41)
    close(42)
    close(43)

endif

write(*,*) '>> Hamiltonian Diagonalization'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DYNAMICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(c_i(vectorstotal))

c_i=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)

call injection_dynamics(HT,hamiD,eigvals,vectorstotal,initialVec1,norm,c_i)

write(*,*) '>> Dynamics'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ENTANGLEMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(red_rho(4,4))

red_rho=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)

call reduced_density_matrix(HT,vectorstotal,red_rho,c_i)

!!Save reduced density matrix
if (files) then
    open(unit=89,file='reduced_rho.data',status='unknown')
    do i=1,4
        write(89,*) (red_rho(i,j),j=1,4)
    enddo
    close(89)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ENTROPY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PLOTTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !Writes in a file data needed for plots
    open(unit=46,file='info.data',status='unknown')

    401 FORMAT ("GRAPHICAL=",L)
    write(46,401) graphical

    501 FORMAT ("REALISATIONS=",A)
    write(tmp,'(i5.4)') num_realisations
    write(46,501) adjustl(trim(tmp))

    601 FORMAT ("N=",A)
    write(tmp,'(i5.4)') N
    write(46,601) adjustl(trim(tmp))

    701 FORMAT ("VECTORS=",A)
    write(tmp,'(i5.4)') vectorstotal
    write(46,701) adjustl(trim(tmp))

    801 FORMAT ("TOTALTIME=",A)
    write(tmp,'(f6.2)') totaltime
    write(46,801) adjustl(trim(tmp))

    901 FORMAT ("INITIALVEC=",A)
    write(tmp,'(i5.2)') initialVec1
    write(46,901) adjustl(trim(tmp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FREE SPACE AND CLOSE FILES AND CLEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (output) then
    close(unit=40)
endif

if (graphical) then
    close(unit=46)
endif

!!deallocation **VERY IMPORTANT**
deallocate(H1)
deallocate(H2)
deallocate(H3)
deallocate(HT)
deallocate(hami)
deallocate(hamiD)
deallocate(eigvals)
deallocate(rwork)
deallocate(work)
deallocate(red_rho)
deallocate(c_i)
deallocate(Noise)

end program

!*********************TO BE REPROGRAMED**********************!

subroutine noise_sub1 (seed)

use noisemod
implicit none
integer :: seed             !may or may not be 32 bit...
!32 bit integer version of seed
integer(kind=bit_32) :: iseed              !random number seed

!f90 intrinsic time call
character(len=10) :: system_time           !length is crucial ...
real (kind=dubp)    :: rtime

if (bit_size(iseed)/=32) call io_abort('Error in algor_set_random_seed - 32-bit integer kind type problem')

if (seed == 0 ) then
call date_and_time(time=system_time)  !character string hhmmss.xxx
read (system_time,*) rtime            !convert to real
rtime = rtime * 1000.0_dubp             !0<rtime<235959999.0 which fits within huge(1)
else
rtime = real(abs(seed),kind=dubp)          !convert seed to real
end if

!and then convert to bit_32 size integer
iseed = int(rtime,kind=bit_32)              !must fit within huge(1)

ix=ieor(777755555_bit_32,iseed)                   !Marsaglia generator
iy=ior(ieor(888889999_bit_32,iseed),1_bit_32)     !Parks-Miller generator

seed = int(rtime)                                 !return the seed that was used in default integer


return
end subroutine noise_sub1

function algor_uniform_random()
!=========================================================================!
! Return a single random deviate ~ uniform [0,1] or [-1,1] or [-0.5,0.5]  !
! Based on Park-Miller "minimal standard" generator with Schrage's method !
!  to do 32-bit multiplication without requiring higher precision, plus   !
!  additional Marsaglia shift to suppress any weaknesses & correlations.  !
! Using two independent methods greatly increases the period of the       !
!  generator, s.t. resulting period ~2*10^18                              !
! NB Routine is only set to work with 32 bit integers!                    !
!-------------------------------------------------------------------------!
! References:                                                             !
!   S.K. Park and K.W. Miller, Commun. ACM, 31, p1192-1201 (1988)         !
!   L. Schrage, ACM Trans. Math. Soft., 5, p132-138 (1979)                !
!   G. Marsaglia, Linear Algebra and its Applications, 67, p147-156 (1985)!
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   algor_uniform_random => required random deviate                       !
!-------------------------------------------------------------------------!
! Parent module variables used:                                           !
!   ix as next 32-bit integer in Marsaglia generator (updated)            !
!   iy as next 32-bit integer in Park-Miller generator (updated)          !
!-------------------------------------------------------------------------!
! Written by Matt Probert, v0.1, 01/07/2000                               !
!=========================================================================!
use noisemod
implicit none
real(kind=dubp)                 :: algor_uniform_random

!NB We use 3 logical XOR operations to do Marsaglia shift
!=> hard-wire 3 shift values (not all triplets any good)
!=> entire routine preset to only work with 32 bit integers.

!local variables ...
integer(kind=bit_32)            :: iy_tmp       !working value to force integer division
integer(kind=bit_32), parameter :: iy_max=2147483647 !2^31-1
real(kind=dubp), parameter        :: inv_iy_max=1.0_dubp/2147483647.0_dubp
real(kind=dubp)                   :: notright

!do Marsaglia shift sequence, period 2^32-1, to get new ix
ix=ieor(ix,ishft(ix, 13_bit_32))
ix=ieor(ix,ishft(ix,-17_bit_32))
ix=ieor(ix,ishft(ix,  5_bit_32))

!Park-Miller sequence, period iy_max-1, to get new iy
iy_tmp=iy/127773_bit_32                         !NB integer division
iy=16807_bit_32*(iy-iy_tmp*127773_bit_32)-iy_tmp*2836_bit_32  !next value of iy
if (iy < 0_bit_32) iy=iy+iy_max                 !integer modulo iy_max

!Combine ix and iy to get new random number, rescale onto range [0,1]
notright=inv_iy_max*ior(iand(iy_max,ieor(ix,iy)),1_bit_32)       !with masking to ensure non-zero

!Between -1 to 1
!algor_uniform_random=(notright*2.0_dubp)-1_dubp when {0-1}: !notright

!Between -0.5 and 0.5
algor_uniform_random=(notright*1.0_dubp)-0.5_dubp!notright

!Between 0 and 1
!algor_uniform_random=notright

return
end function algor_uniform_random

!**********************************************************************************!