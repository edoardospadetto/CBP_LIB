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
    use misc_mod
   

    implicit none
    ! #### MPI Variables ####
    
    integer :: rank ! the process ID
    integer :: nprocs ! number of processes
    integer :: rowproc, colproc
    integer :: ierr , idxcol, idxrow, ii, jj,kk
    !Global

    
    integer :: dimA(2), dimQ(2), dimR(2)
    
    integer , dimension(4):: descriptorQ, descriptorR
    integer :: localdimQ(2), localdimR(2), stat, ngpu
   
    double precision :: start, finish, flops 
	character*8 :: args(3)
    
    type(c_ptr) :: context
    
    CTYPE , allocatable, dimension(:,:) :: localQ, localR
    
        
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
    !-----------------------------------------------	
    	call get_command_argument(1,args(1))
	    read(args(1),*,iostat=stat)  rowproc
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if
	    
	    call get_command_argument(2,args(2))
	    read(args(2),*,iostat=stat)  colproc
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if

	    call get_command_argument(3,args(3))
	    read(args(3),*,iostat=stat)  ngpu
	    if (stat .ne. 0 ) then 
	    	print *, "error"
	    	STOP 
	    end if
	!----------------------------------------------------------    
    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)

    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if
    
  
 
do kk = 1,100
    dimA=(/50*kk,50*kk/)
    dimR=dimA
    dimQ=[dimA(1),dimA(1)]


    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    !allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call cbp_init_local_matrix_descriptor(descriptorR, [dimA(1), dimA(2)], localdimR, rowproc, colproc, rank)  
    call cbp_init_local_matrix_descriptor(descriptorQ, [dimA(1), dimA(1)], localdimQ, rowproc, colproc, rank)
   
  

    !------------------- Local Matrices
    allocate(localR(localdimR(1),localdimR(2)))
    allocate(localQ(localdimQ(1),localdimQ(2)))

    localQ = CZERO
    localR = CZERO
    call c_random(localR)
    localr = (localr-0.5)*14.0
    !------------------- 

	!----------------------------
    call init_cublasXT(4,context)
    !----------------------------

 
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(rank .eq. 0) start = MPI_Wtime();
	 
    call cbp_QR_decompose_Householder(localR, descriptorR, localQ, descriptorQ, dimA, & 
    								  rowproc, colproc,rank,nprocs, context )
   
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if(rank .eq. 0) finish = MPI_Wtime();
    
 	
    call destroy_cublasXt(context)


    deallocate (localQ) 
    deallocate (localR)
 	
 	if(rank .eq. 0) print*, dimA(1) , finish -start
 
 end do 
    call MPI_FINALIZE(ierr)
end program 
