#ifdef CBPCMPLX
#define CTYPE complex
#define MPI_CTYPE MPI_COMPLEX
#define CZERO cmplx(0.0,0.0)
#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#define CZERO 0.0
#endif

#ifdef CBPZMPLX 
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#define CZERO dcmplx(0.0)
#endif

#ifdef CBPREAL 
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#define CZERO 0.0
#endif


program passm
    use cbp_DistributeMOD
    use mpi
   

    implicit none
    ! #### MPI Variables ####
     CTYPE , dimension(:,:) , allocatable :: A,A0,R,L
    
    integer :: rank ! the process ID
    integer :: nprocs ! number of processes
    integer , dimension(:,:) , allocatable :: process_grid 
    integer :: rowproc, colproc
    integer :: ierr , idxcol, idxrow, ii, jj,kk
    !Global
    
   
    integer :: dimA(2), dimL(2), dimR(2)
    
    integer , dimension(4):: descriptorA,descriptorL, descriptorR!, descriptorR, descriptorRtemp
    integer , dimension(4):: descriptorQtemp, descriptorQitemp, descriptorQinv, descriptorQres
    
    integer :: localdimA(2), localdimL(2), localdimR(2)
    

    double precision :: start, end, flops 
    
    integer :: zeros, zerosall, expected
    
    type(c_ptr) :: context
    
#ifdef CBPCMPLX
        real , dimension (:,:), allocatable  ::  rai, rar
#endif
#ifdef CBPZMPLX
        double  precision , dimension(:,:) , allocatable :: rar, rai
#endif

    
    CTYPE , allocatable, dimension(:,:) :: localA, localL, localQ, localR, localQinv, localQres
    
    rowproc = 2
    colproc = 1
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
do kk = 1,1
    dimA=(/8,8/)
    dimR=[dimA(2), dimA(2)]
    dimL=(/dimA(1),dimA(1)/)

   
    !dimA = kk*2000 + dimA

 
    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    !allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call cbp_init_local_matrix_descriptor(descriptorA, [dimA(1), dimA(2)], localdimA, rowproc, colproc, rank)  
    call cbp_init_local_matrix_descriptor(descriptorL, [dimA(1), dimA(1)], localdimL, rowproc, colproc, rank)
    call cbp_init_local_matrix_descriptor(descriptorR, [dimA(2), dimA(2)], localdimR, rowproc, colproc, rank)
 
    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if

    !------------------- Local Matrices
        
   
    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localL(localdimL(1),localdimL(2)))
    allocate(localR(localdimR(1),localdimR(2)))
    
   ! allocate(localA0(localdimA(1),localdimA(2)))
   
    !allocate(localQ(localdimQ(1),localdimQ(2)))
    !allocate(localQinv(localdimQ(1),localdimQ(2)))
    !allocate(localQres(localdimQ(1),localdimQ(2)))
    
    localA= CZERO
    localL = CZERO
    localR = CZERO
    
     
    !-------------------! Global Matrices

    if(rank .eq. 0 ) then 
       
        allocate(A(dimA(1), dimA(2))) ! good one
        allocate(A0(dimA(1), dimA(2))) 
        allocate(L(dimL(1) , dimL(1)))
        allocate(R(dimR(1), dimR(2)))
    
    end if


#if defined(CBPCMPLX) |defined (CBPZMPLX) 
    allocate(rar(localdimA(1), localdimA(2)))
    allocate(rai(localdimA(1), localdimA(2)))
   do ii = 1, rank+1
   call random_number(rar)
   call random_number(rai)
   end do 
#ifdef CBPCMPLX       
   localA = cmplx(rar, rai)
#endif
#ifdef CBPZMPLX
   localA = dcmplx(rar, rai)
#endif
   deallocate(rar)
   deallocate(rai) 
#endif


#if defined(CBPREAL) | defined(CBPDOUBLE)
    !Fill A, B with random numbers
     do ii = 1, rank+1
    call random_number(localA)
	end do 
#endif

    ! R = A
    locala = (locala- 0.5)*14.0
    
    !localR=localA
    !descriptorR = descriptorA
    !do ii = 1, dimA(1) 
    !	print*, localA(ii,:) 
    !end do 
    !print*, " "
    call init_cublasXT(1,context)
    if(rank .eq. 0) then 
            start = MPI_Wtime();
    endif
    !Needed for get back matrix old implementation sorry :( 
    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
    call GetBackMatrix(localA, descriptorA, A0, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
         	  	
  
         do jj = 0, nprocs-1
	if (rank.eq.jj)then 
	print*,"RANK" ,rank
	print*, " " 
	do ii = 1,descriptorA(1)

	print "(*('('sf6.2xspf6.2x'i)':x))", localA(ii,:descriptorA(2))
	!write(*,*), M(ii,:)
	end do 

	print*, " " 
	end if
	CALL SLEEP(1)
	call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
	print*, " "
	end do 
    call cbp_Bidiag(localA, descriptorA, localR, descriptorR, localL, & 
    			descriptorL, dimA ,rowproc,colproc, rank, nprocs, context)
    			
 
    
    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
         	  	
    call GetBackMatrix(localL, descriptorL, L, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
    call GetBackMatrix(localR, descriptorR, R, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
       
       
	
    if (rank.eq.0) then 
    
    print*, " " 
    print*, " " 
    
    do ii = 1, dimA(1) 
    	print  "(*('('sf6.2xspf6.2x'i)':x))", A(ii,:)  
    end do 
    A(1,1) = CZERO
    do ii = 2, minval(dimA) 
    		A(ii,ii ) =  CZERO
    		A(ii-1,ii) =  CZERO 
    	 
    end do 
      
    !end do 
   	 print*,"this ->", maxval(abs(A)) 
    deallocate(A0) 
    deallocate(R) 
    deallocate(L) 
    end if 
    end do 
    call MPI_FINALIZE( ierr )
  
end program 




