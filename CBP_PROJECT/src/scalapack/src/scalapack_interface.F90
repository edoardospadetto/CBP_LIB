#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX

#define pxgemm pcgemm
#define pxgeqrf pcgeqrf
#define pxelset pcelset
#define pxelget pcelget

#define pxxxmqr pcunmqr
#define pxgesvd pcgesvd
#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define CTYPE0 0.d0
#define CTYPE1 1.d0
#define MPI_CTYPE MPI_DOUBLE_PRECISION

#define pxgemm pdgemm
#define pxgeqrf pdgeqrf
#define pxelset pdelset
#define pxelget pdelget

#define pxxxmqr pdormqr

#define pxgesvd pdgesvd

#endif

#ifdef CBPZMPLX 
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX

#define pxgemm pzgemm
#define pxgeqrf pzgeqrf
#define pxelset pzelset
#define pxelget pzelget

#define pxxxmqr pzunmqr

#define pxgesvd pzgesvd
#endif

#ifdef CBPREAL 
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL

#define pxgemm psgemm
#define pxgeqrf psgeqrf
#define pxelset pselset
#define pxelget pselget

#define pxxxmqr psormqr
#define pxgesvd psgesvd
#endif


#define BF_ 30

!##################################################
!Module to handle matrices in a distribute fashion
!#################################################
module scalapack_interface


use debug_module

implicit none

contains
INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )


INTEGER              EXTRABLKS, MYDIST, NBLOCKS
INTEGER 		N, NB, IPROC, ISRCPROC, NPROCS
INTRINSIC            MOD

MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )

NBLOCKS = N / NB

NUMROC = (NBLOCKS/NPROCS) * NB

EXTRABLKS = MOD( NBLOCKS, NPROCS )

IF( MYDIST.LT.EXTRABLKS ) THEN
NUMROC = NUMROC + NB

ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
NUMROC = NUMROC + MOD( N, NB )
END IF

RETURN

END


!Computes the sum for distribute matrices, safely.
subroutine dsum(A, descA, B, descB, C, descC)

integer, dimension(:) :: descA, descB, descC
double complex, dimension(:,:) :: A, B, C
double complex :: alpha, beta
integer ii, jj

call breakifn("Wrong dimensions", all(descA(3:4) .eq. descB(3:4)), .true.)
call breakifn("Wrong dimensions", all(descA(3:4) .eq. descC(3:4)), .true.)

if (all((descA .eq. descB) .and. (descB .eq. descC))) then 
C = A + B	
!else
!	do ii = 1 , descA(3)
!		do jj = 1, descA(4)
!		call pzelget('A','I',A,ii, jj ,descA, alpha)
!		call pzelget('A','I',B,ii, jj ,descB, beta)
!		 
!		
!		call pzelset(C, ii ,jj , descC, alpha + beta)
!		end do 
!	end do 
end if

end subroutine

! Distributed matrix builder
subroutine set_matrix(M, iam, A, desca)

	use mpi
	use matrix_interface

	implicit none

	! input
	CTYPE , dimension(:,:), intent(INOUT) :: M, A
	integer, intent(INOUT) ::  iam
	integer, dimension(:), intent(INOUT) ::  desca
	! variables
	integer , dimension(2) :: sizeG 
	integer :: ierr
	integer :: ii, jj 

	! Check if M is square
	sizeG = shape(M)
	! Use M from the first process
	if(iam .eq. 0) then 
		call MPI_BCAST(m,product(sizeg), MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr)
	else 
		call MPI_BCAST(m,product(sizeg), MPI_CTYPE, 0, MPI_COMM_WORLD,ierr )
	end if
	

	if (ierr .ne. 0) then 
		print*, "Error building matrix!"
		STOP
	end if

	! Distribute M 
	DO ii = 1, sizeg(1)
		DO jj = 1, sizeg(2)
		CALL pxelset(A, ii, jj, DESCA, M(ii,jj))	
		END DO
	END DO

end subroutine 

subroutine get_matrix(M, A, desca)
integer :: ii, jj 
integer , dimension(2):: sizeG 
CTYPE , dimension(:,:) :: A, M 
integer, dimension(9) :: descA


sizeG = shape(M)


DO ii = 1, sizeg(1)
DO jj = 1, sizeg(2)

CALL pxelget('A', 'D' ,M(ii,jj) ,A,ii, jj, DESCA)	
END DO
END DO

end subroutine

subroutine sca_identity(A, desca)

	use mpi
	use matrix_interface

	implicit none

	! input
	CTYPE , dimension(:,:), intent(INOUT) :: A
	integer, dimension(9), intent(INOUT) ::  desca
	! variables
	
	integer :: ii
	A = CTYPE0
	if (descA(3) .ne. descA(4)) then 
		print*, "SCA_IDENTITY : Error not square Matrix"
	end if 
	DO ii = 1, descA(3)
		CALL pxelset(A, ii, ii, DESCA, CTYPE1 )	
	END DO

end subroutine 

! Simple interface to print a distributed matrix
subroutine printmat(A, descA)

	integer, dimension(:) :: desca 
	double complex, dimensioN(:,:) :: A
	double complex, dimension(30) :: work

	CALL PZLAPRNT(descA(3), descA(4), A, 1, 1, DESCA, 0, 0, 'A', 6, WORK)

end subroutine 

! Distributed matrix diagonalization 
subroutine ddzm(A, desca, Z, descz, W)

	use matrix_interface

	implicit none

	double complex, dimension(:,:), intent(INOUT) :: A, Z
	integer, dimension(:), intent(INOUT) :: desca, descz
	double precision, dimension(:), intent(INOUT) :: W
	double complex, dimension(:), allocatable :: rwork
	integer :: lwork, lrwork
	integer :: info
	double complex, dimension(:), allocatable :: work

	allocate(work(1))
	allocate(rwork(1))

	CALL pzheev('V', 'U', size(w), A, 1, 1, DESCA, W, Z, 1, 1, DESCZ, WORK, -1, RWORK, -1, INFO)

	lwork = work(1)
	lrwork = int(rwork(1))

	deallocate(work)
	deallocate(rwork)

	allocate(work(lwork))
	allocate(rwork(lrwork))	

	CALL pzheev('V', 'U', size(w), A, 1, 1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, INFO)

end subroutine

!diagonalize matrix
subroutine ddzm2(A, desca, Z, descz, W)

	use matrix_interface

	implicit none

	double complex, dimension(:,:), intent(INOUT) :: A, Z
	integer, dimension(:), intent(INOUT) :: desca, descz
	double precision, dimension(:), intent(INOUT) :: W
	double precision, dimension(:), allocatable :: rwork
	integer :: lwork, lrwork, liwork, m, nz
	integer :: info
	double complex, dimension(:), allocatable :: work
	integer, dimension(:), allocatable :: iwork


	allocate(work(1))
	allocate(iwork(1))
	allocate(rwork(1))

	!CALL pzheevr('V', 'A', 'U', size(w), A, 1, 1, DESCA, -100.d0, 100.d0, 0, 100, m, &
	!			 & nz, W, Z, 1, 1, DESCZ, WORK, -1, RWORK, -1, IWORK, -1, INFO)

	CALL pzheevr('V', 'I', 'U', size(w), A, 1, 1, DESCA, -100.d0, 100.d0, 1, 3, m, &
	& nz, W, Z, 1, 1, DESCZ, WORK, -1, RWORK, -1, IWORK, -1, INFO)

	lwork = work(1)
	lrwork = rwork(1)
	liwork = iwork(1)

	deallocate(work)
	deallocate(rwork)
	deallocate(iwork)

	allocate(work(lwork))
	allocate(rwork(lrwork))	
	allocate(iwork(liwork))

	!CALL pzheevr('V', 'A', 'U', size(w), A, 1, 1, DESCA, -100.d0, 100.d0, 0, 100, m, &
	!			 & nz, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

	CALL pzheevr('V', 'I', 'U', size(w), A, 1, 1, DESCA, -100.d0, 100.d0, 1, 3, m, &
	& nz, W, Z, 1, 1, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)

	deallocate(work)
	deallocate(rwork)

end subroutine

function locr(K, rowproc) result (val)
	integer :: K, rowproc
	integer :: val , ii 
	val = 0
	do ii= 0, rowproc-1
		val = max(  NUMROC( K,  BF_ , ii , 0, rowproc ), val)
	end do  

end function 

function locc(K, colproc) result (val)
	integer :: K , colproc 
	integer:: val , ii 
	val = 0
	do ii= 0, colproc-1
		val = max(  NUMROC( K,  BF_ , ii , 0, colproc ), val)
	end do  

end function

function locrc(K, procs) result (val)

	integer, dimension(2) :: val,  procs , K 
	integer :: ii
	val = 0 
	do ii = 0, procs(1)-1 

		val(1) = max(  NUMROC( K(1),  BF_ , ii , 0, procs(1) ), val(1))

	end do
	do ii= 0, procs(2)-1
		val(2) = max(  NUMROC( K(2),  BF_ , ii , 0, procs(2) ), val(2))
	end do  

	end function


function g2l(global, descA, dimA) result(local) 

	integer , dimension(9) :: descA 
	integer, dimension(2) :: npxsca, coord_sca  , global, blockidx, idx, dimA
	integer, dimension(4) :: local 

	call BLACS_GRIDINFO(descA(2), npxsca(1) , npxsca(2) , coord_sca(1) , coord_sca(2) )

	if (any(global .gt. dimA))  then 
	local  = -1
	return
	else 

	global = blockidx*npxsca*descA(6:7) + coord_sca*descA(6:7) + idx 

	blockidx = global/(npxsca*descA(6:7))
	coord_sca = mod(global,npxsca*descA(6:7))/descA(6:7) 
	idx = mod( mod(global,npxsca*descA(6:7)),descA(6:7) )

	local(1:2) = blockidx*descA(6:7)+idx
	local(3:4) = coord_sca
	return 
	end if 


end function 

!correct if i know in advance the process
function l2g(local, descA, dimA) result(global)

	integer , dimension(9) :: descA 
	integer, dimension(2) :: npxsca, coord_sca, global, blockidx, idx, dimA
	integer, dimension(4) :: local 


	call BLACS_GRIDINFO(descA(2), npxsca(1) , npxsca(2) , coord_sca(1) , coord_sca(2) )
	coord_sca= local(3:4)

	blockidx= local(1:2) / descA(6:7) 
	idx = mod(local(1:2),descA(6:7) )

	global = blockidx*npxsca*descA(6:7) + coord_sca*descA(6:7) + idx 

	if (any(global .gt. dimA)) global = -1

	return

end function

!subroutine sca_get_row(A, descA, v, descv, globalstart, globalend, dimA ) 
!	integer, dimension(9) :: descA, descv 
!	integer :: globalstart, globalend , ii, jj , gridblock
!	CTYPE, dimension (:) :: A, v 
!	integer, dimension(2) :: npxsca, coord_sca , dimA
!	
!	call BLACS_GRIDINFO(descA(2), npxsca(1) , npxsca(2) , coord_sca(1) , coord_sca(2) )
!	
!	gridblock = dimA(1)/ ( npsca(1)*descA(6))
!	
!	if ( mod(dimA(1), npsca(1)*descA(6) ) .gt. 0 ) gridblock = gridblock +  1 
!	
!	!find first local row globally ge globalstart 
!	!g2l(globalstart,globalend)
!	!who owns start 
!	
!	do ii = 1, gridblock    !number of "grid block"
!		do jj = descA(6)*coord_sca(1)*1, descA(6)
!			!put in the same index 
!			v(:) =
!			 
!	end do 
!	 
!	
!end subroutine 

subroutine sca_qr(A,Q,descA,descQ)
	use mpi
	implicit none 


	integer :: IA, JA, LWORK, INFO
	CTYPE, dimension(:), allocatable ::  TAU
	CTYPE, dimension(:), allocatable :: WORK
	integer, dimension(2) :: npxsca, coord_sca, dimA
	integer, dimension(9) :: descA, descQ
	integer :: ii,  t
	CTYPE, dimension(:,:)::  A, Q
	CTYPE, dimension(:,:), allocatable :: Qtmp, Atemp,A2, Qt, Q2
	CTYPE :: v, va
	!Compute the svd 

	call BLACS_GRIDINFO(descA(2), npxsca(1) , npxsca(2) , coord_sca(1) , coord_sca(2) ) 
	
	if (all(locrc([descA(4), descA(4)],npxsca) .lt. shape(Q)) ) print*, "sca_QR: Q matrix seems too small"	
	t = maxval(locrc([descA(4), descA(4)], npxsca)) 
	
	call DESCINIT(descQ, descA(4), descA(4), BF_, BF_, 0, 0, descA(2), & 
					T, info)
	
	
	
	
	allocate(work(1))
	allocate (TAU(locc(minval(descA(3:4)), npxsca(2))) ) 
	
	lwork = -1
	call pxgeqrf( descA(3), descA(4), A, 1, 1, DESCA, TAU, WORK, LWORK, INFO )

	lwork = int(work(1))
	deallocate(work)
	allocate(work(lwork))
	
	call pxgeqrf( descA(3), descA(4), A, 1, 1, DESCA, TAU, WORK, LWORK, INFO )
	
	call sca_identity(Q, descQ)

	lwork = -1 
	call pxxxmqr( 'L', 'N', descA(3), descA(3), descA(4)-1, A, 1, 1, DESCA, TAU, &
                         Q, 1, 1, DESCQ, WORK, LWORK, INFO )
                         
    lwork = int(work(1))
    deallocate(work)
	allocate(work(lwork))
    call pxxxmqr( 'L', 'N', descA(3), descA(3), descA(4)-1, A, 1, 1, DESCA, TAU, &
                         Q, 1, 1, DESCQ, WORK, LWORK, INFO )
	deallocate(TAU) 
	deallocate(work)
	!START FROM Q  
	!number of times blocks
	!number of rowproc



end subroutine 

! Distributed matrix multiplication
subroutine sca_matmul(A, descA, B, descB, C, descC)

	use mpi 

	implicit none

	CTYPE, dimension(:,:), intent(INOUT) :: A, B, C
	integer, dimension(:), intent(INOUT) :: descA, descB, descC 

	
		if (descA(4) .ne. descB(3)) then 
		write(*,*) "Invalid matrix dimensions!"
		STOP
		end if
		if (descC(4) .ne. descB(4))  then 
		write(*,*) "Invalid matrix dimensions!"
		STOP
		end if
		if (descC(3) .ne. descA(3)) then 
		write(*,*) "Invalid matrix dimensions!"
		STOP
		end if

	C = dcmplx(0.d0,0.d0)

	call pxgemm('N', 'N', descA(3), descB(4), descA(4), CTYPE1, & 
				& A, 1, 1, descA, B, 1, 1, descB, CTYPE0, C, 1, 1, descC)

end subroutine

subroutine sca_svd(A,U,VT,descA,descU,descVT, JOBU,JOBVT)
	use mpi 
	implicit none
	CTYPE, dimension(:,:), intent(INOUT) :: A, U, VT
	integer, dimension(9) :: descA, descU, descVT
	character*1 :: jobu, jobvt
	integer :: lwork , info
	CTYPE, dimension(:) , allocatable :: WORK
#if defined(CBPCMPLX) | defined(CBPREAL)	
	real, dimension(:) , allocatable :: S
#endif
#if defined(CBPZMPLX) | defined(CBPDOUBLE)	
	double precision, dimension(:) , allocatable :: S
#endif
#ifdef CBPCMPLX
	real , dimension(:), allocatable ::RWORK
	allocate(rwork(1+4*maxval(descA(3:4))))
#endif
#ifdef CBPZMPLX
	double precision , dimension(:), allocatable ::RWORK
	allocate(rwork(1+4*maxval(descA(3:4))))
#endif 

	lwork = -1	
	allocate(S(minval(descA(3:4))))
	allocate(work(1))

#if defined(CBPZMPLX) | defined(CBPCMPLX)
	call pxgesvd(JOBU,JOBVT,DESCA(3),DESCA(4),A,1,1,DESCA,S,U,1,1,DESCU, & 
                      VT,1,1,DESCVT,WORK,LWORK,RWORK,info)
    

#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)  
  	call pxgesvd(JOBU,JOBVT,DESCA(3),DESCA(4),A,1,1,DESCA,S,U,1,1,DESCU, & 
                      VT,1,1,DESCVT,WORK,LWORK,info)  
#endif

    print*, info
    lwork = ceiling(abs(work(1)))
    deallocate(work) 
    allocate(work(lwork))
    
#if defined(CBPZMPLX) | defined(CBPCMPLX)
	call pxgesvd(JOBU,JOBVT,DESCA(3),DESCA(4),A,1,1,DESCA,S,U,1,1,DESCU, & 
                      VT,1,1,DESCVT,WORK,LWORK,RWORK,info)
    deallocate(rwork)
#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)  
  	call pxgesvd(JOBU,JOBVT,DESCA(3),DESCA(4),A,1,1,DESCA,S,U,1,1,DESCU, & 
                      VT,1,1,DESCVT,WORK,LWORK,info) 
    print*, info  
#endif
    

	deallocate(S)
	deallocate(work)
end subroutine




! -----------------------------------------------------------
! SUBROUTINES/FUNCTIONS TO BE TESTED OR NOT WORKING
! -----------------------------------------------------------

!Must be tested
subroutine distributedvsfull(A, descA, B, iam )
use mpi
implicit none

CTYPE, dimension(:,:) :: A
CTYPE, dimension (:,:), allocatable :: Bprime
CTYPE , dimension(:,:) :: B
integer , dimension(:) :: descA
integer, dimension(9) :: descB
double complex :: alpha, tot
integer :: ii, jj ,iam,ierr

descB = descA

allocate(Bprime(size(A, dim = 1), size(A, dim= 2)))



Bprime = CTYPE0

if(iam .eq. 0) then 

call MPI_BCAST(B,product(shape(B)), MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr)
else 
call MPI_BCAST(B,product(shape(B)), MPI_CTYPE, 0, MPI_COMM_WORLD,ierr )
end if



do ii = 1 , descA(3)
do jj = 1, descA(4)
call pxelset(Bprime, ii ,jj , descB, B(ii,jj))
end do 
end do

print*,"dvsfull", sum(A-Bprime),sum(A), sum(Bprime), descA(3:4)
deallocate(Bprime)

end subroutine

! Get column
function getcol(A, descA, col) result(res)

double complex, dimension(:,:) :: A
integer:: dimr, dimc, descA(9), col, ii
double complex :: alpha
!parameter(dimr = descA(3), dimc = descA(4))
double complex, dimension(descA(3)) :: res

do ii = 1, descA(3)
call pzelget('A', ' ', res(ii), A, ii, col, descA)
end do

return

end function


end module scalapack_interface
