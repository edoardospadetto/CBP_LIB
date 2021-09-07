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
    
    integer :: rank ! the process ID
    integer :: nprocs ! number of processes
    integer , dimension(:,:) , allocatable :: process_grid 
    integer :: rowproc, colproc
    integer :: ierr , idxcol, idxrow, ii, jj,kk
    !Global
    
    CTYPE , dimension(:,:) , allocatable :: A,A0,R,Q, Qres
    
    integer :: dimA(2), dimQ(2), dimR(2)
    
    integer , dimension(4):: descriptorA,descriptorA0, descriptorQ, descriptorR, descriptorRtemp
    integer , dimension(4):: descriptorQtemp, descriptorQitemp, descriptorQinv, descriptorQres
    
    integer :: localdimA(2), localdimQ(2), localdimR(2)
    

    double precision :: start, end, flops 
    
    integer :: zeros, zerosall, expected
    
    type(c_ptr) :: context
    
#ifdef CBPCMPLX
        real , dimension (:,:), allocatable  ::  rai, rar
#endif
#ifdef CBPZMPLX
        double  precision , dimension(:,:) , allocatable :: rar, rai
#endif

    
    CTYPE , allocatable, dimension(:,:) :: localA, localA0, localQ, localR, localQinv, localQres
    
    rowproc = 1
    colproc = 1
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
do kk = 1,20
    dimA=(/30,30/)
    dimR=dimA
    dimQ=(/dimA(1),dimA(1)/)

   
    dimA = kk*100 + dimA

 
    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    !allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call cbp_init_local_matrix_descriptor(descriptorA, [dimA(1), dimA(2)], localdimA, rowproc, colproc, rank)  
    call cbp_init_local_matrix_descriptor(descriptorQ, [dimA(1), dimA(1)], localdimQ, rowproc, colproc, rank)
   
 
    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if

    !------------------- Local Matrices
        
   
    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localA0(localdimA(1),localdimA(2)))
    allocate(localR(localdimA(1),localdimA(2)))
    allocate(localQ(localdimQ(1),localdimQ(2)))
    allocate(localQinv(localdimQ(1),localdimQ(2)))
    allocate(localQres(localdimQ(1),localdimQ(2)))
    
    localA= CZERO
    localQ = CZERO
    localR = CZERO
    
     
    !-------------------! Global Matrices

    if(rank .eq. 0 ) then 
       
        allocate(A(dimA(1), dimA(2))) ! good one
        allocate(A0(dimA(1), dimA(2))) 
        allocate(Q(dimA(1) , dimA(1)))
        allocate(Qres(dimA(1) , dimA(1)))
        allocate(R(dimA(1), dimA(2)))
    
    end if


#if defined(CBPCMPLX) |defined (CBPZMPLX) 
    allocate(rar(localdimA(1), localdimA(2)))
    allocate(rai(localdimA(1), localdimA(2)))
    
   call random_number(rar)
   call random_number(rai)

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

    call random_number(localA)

#endif

   
    locala = (locala-0.5)*14.0
    localR=localA
    descriptorR = descriptorA

    call init_cublasXT(1,context)
    if(rank .eq. 0) then 
            start = MPI_Wtime();
    endif


    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)


    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
    
    if (rank.eq.0) then 
    	print*, "start QR decomposition"
    end if 
    
    
    call cbp_QR_decompose_Givens(localR, descriptorR, localQ, descriptorQ, dimA, & 
    		[dimA(1), dimA(1)],rowproc, colproc,rank,nprocs, context )
    		
  
    zeros = 0 
    
    !-----------CHECK R ---------------
    
   ! do ii = 1, descriptorR(1)
   ! 	!
   ! 	do jj = 1, descriptorR(2) 
   ! 		if (abs(localR(ii,jj)) .lt. 1e-3) then 
   ! 			zeros = zeros +1
   ! 		end if 	
   ! 	end do 
   ! end do 
   call MPI_Reduce(zeros, zerosall, 1, MPI_INTEGER,&
               MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

   
   expected = 0
  ! do ii = 2, dimA(1)
  ! 	if ((ii-1) .lt. dimA(2)) then 
  ! 		expected = (ii-1) + expected
  ! 	else 
  ! 		expected = expected + dimA(2)
  ! 	end if 
  ! end do 
  ! if (rank.eq.0) then 
  ! 	print *, "Obtained and expected zero values in R matrix",  zerosall , expected 
  ! end if
   
   
   !------------------CHECK QR = A ---------------------
   
   descriptorQtemp= descriptorQ
   descriptorRtemp = descriptorR
   
    allocate(process_grid(max(rowproc, colproc), max(rowproc,colproc)))
   
    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
    
    call cbp_matmul_descriptors(descriptorQ, descriptorR, descriptorA0, dimQ, dimR, dimA, & 
         localdimQ, localdimR, localdimA , rowproc, colproc, process_grid,rank) 
    
    call cbp_redispatch(localQ, descriptorQtemp,descriptorQ,rank,nprocs)
    call cbp_redispatch(localR, descriptorRtemp,descriptorR,rank,nprocs)
    
    call cbp_matmul(localQ, localR, localA0, descriptorQ, descriptorR ,& 
                 descriptorA0, rowproc, colproc, rank, process_grid,  context)
    
   
    call GetBackMatrix(localA0, descriptorA0, A0, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
    
   
   
   if (rank .eq. 0) then 
   ! print *, "if low means QR = A" ,  maxval (  abs( A0 - A ) ), sum(A) , sum(A0) 
   end if 
    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
    
    localQinv = localQ 
    descriptorQinv = descriptorQ
 
    call cbp_adjoint(localQinv, descriptorQinv) 
    
    call cbp_matmul_descriptors(descriptorQtemp, descriptorQitemp, descriptorQres, dimQ, dimQ, dimQ, & 
         localdimQ, localdimQ, localdimQ , rowproc, colproc, process_grid,rank) 
    
    call cbp_redispatch(localQinv, descriptorQinv, descriptorQitemp, rank, nprocs) 
    call cbp_redispatch(localQ, descriptorQ, descriptorQtemp, rank, nprocs) 
    
    descriptorQinv = descriptorQitemp
    descriptorQ = descriptorQtemp
    
     localQres = CZERO  
      
    call cbp_matmul(localQ, localQinv, localQres, descriptorQ, descriptorQinv ,& 
                 descriptorQres, rowproc, colproc, rank, process_grid,  context)
    
   
    
    call GetBackMatrix(localQres, descriptorQres, Qres, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
    call GetBackMatrix(localQ, descriptorQ, Q, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
    
   if ( rank .eq. 0 ) then 

#if defined(CBPREAL) | defined(CBPDOUBLE)
     print*, dimA, sum ( matmul(Q, (transpose(Q))  )) -  dble(dimQ(1)),& 
                                                  maxval (  abs( A0 - A ) ), sum(A) , sum(A0) , sum(qres)
endif
#if defined(CBPCMPLX) | defined(CBPZMPLX)
     print*, dimA, sum ( matmul(Q, conjg(transpose(Q)) ))-  dble(dimQ(1)),& 
                                                  maxval (  abs( A0 - A ) ), sum(A) , sum(A0) , sum(qres) 
#endif
   end if 
    
    call destroy_cublasXt(context)

    deallocate (localA)
    deallocate (localA0) 
    deallocate (localQ) 
    deallocate (localR)
    deallocate (localQinv)
    deallocate (localQres)
  
    
    if(rank .eq. 0 ) then 
	 deallocate(Q) 
	 deallocate(R) 
	 deallocate(A) 
	 deallocate(A0) 
    end if
   
    deallocate (process_grid)
 
 end do 
    call MPI_FINALIZE(ierr)
end program 
