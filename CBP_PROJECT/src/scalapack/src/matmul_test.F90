#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX

#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define CTYPE0 0.0
#define CTYPE1 1.0
#define MPI_CTYPE MPI_DOUBLE_PRECISION

#endif

#ifdef CBPZMPLX 
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX

#endif

#ifdef CBPREAL 
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL

#endif

#define BF_ 30

program matmul_test 
	use misc_mod
	use scalapack_interface
	use mpi
	implicit none 
	CTYPE, dimension(:,:) , allocatable :: A, B, C
	integer , dimension(9) :: descA, descB, descC
	integer, dimension(2) :: dimA ,t 
	integer :: MYROW, MYCOL, IAM, CONTEXT , NPROW, NPCOL,NPROCS, stat, info 
	character*8 :: args(2)
	integer :: ii, jj , ierr
	real :: start, finish
    
	  ! ---- INITIALIZE BLACS -------------------------------------------------------
	    call get_command_argument(1,args(1))
	    read(args(1),*,iostat=stat)  NPROW
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if
	    
	    call get_command_argument(2,args(2))
	    read(args(2),*,iostat=stat)  NPCOL
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if
	    
	    !NPROW = 3

	    CALL BLACS_PINFO(IAM, NPROCS)
	
	    
	    IF (NPROCS.LT.1) THEN
		CALL BLACS_SETUP(IAM, NPROW*NPCOL)
	    END IF
	    
	    CALL BLACS_GET(-1, 0, CONTEXT)
	    CALL BLACS_GRIDINIT(CONTEXT, 'R', NPROW, NPCOL)
	    CALL BLACS_GRIDINFO(CONTEXT, NPROW, NPCOL, MYROW, MYCOL)
		
	    IF (MYROW .eq. -1) THEN
		WRITE(*,*) "ERROR, blacs context not valid."
		CALL BLACS_EXIT(0) 
		STOP
	    END IF 
    ! -----------------------------------------------------------------------------
    	!if(iam.eq.0)  open(unit=127, file='scalapack_matmul.dat')
    	
    	do ii = 1,100,10
			dimA = [500*ii,500*ii]

		    allocate(A(locr(dimA(1),NPROW),locr(dimA(2),NPCOL)))
		    allocate(B(locr(dimA(1),NPROW),locr(dimA(2),NPCOL)))
			allocate(C(locr(dimA(1),NPROW),locr(dimA(2),NPCOL)))
		   
			t = maxval( locrc(DIMa,[NPROW,NPCOL]))
			
			call DESCINIT(descA, dimA(1), dimA(2), BF_, BF_, 0, 0, CONTEXT, & 
					T, info)
			call DESCINIT(descB, dimA(1), dimA(2), BF_, BF_, 0, 0, CONTEXT, & 
					T, info)
			call DESCINIT(descC, dimA(1), dimA(2), BF_, BF_, 0, 0, CONTEXT, & 
					T, info)
	
			call c_random(A)
			call c_random(B)
			
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call cpu_time(start)
			call sca_matmul(A, descA, B, descB, C, descC)
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call cpu_time(finish)
			deallocate(A)
			deallocate(B)
			deallocate(C)
			
			if(iam.eq.0)  print*, dimA(1) , finish-start
		end do 
	
    ! ---- EXIT BLACS -------------------------------------------------------------
    CALL BLACS_GRIDEXIT(CONTEXT)
    CALL BLACS_EXIT(0)
    ! -----------------------------------------------------------------------------
end program 

