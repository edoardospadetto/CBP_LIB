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
    !real , dimension(:,:),  allocatable :: tmpa , tmpb, tmpc
    CTYPE , dimension(:,:) , allocatable :: A, B, C,D
    CTYPE , dimension(:) , allocatable :: lineA
    integer :: dimA(2), dimB(2), dimC(2)
    
    integer , dimension(4):: descriptorA(4), descriptorB(4)
    integer :: descriptorC(4)
    integer :: localdimA(2), localdimB(2), localdimC(2)
    integer :: testr(2), testc(2) , testf(2)

    double precision :: start, end, flops 
    
    type(c_ptr) :: context
    
#ifdef CBPCMPLX
        real , dimension (:,:), allocatable  ::  rar, rbr, rai, rbi
#endif
#ifdef CBPZMPLX
        double  precision , dimension(:,:) , allocatable :: rar, rbr, rai, rbi 
#endif

    
    CTYPE , allocatable, dimension(:,:) :: localA, localB
    CTYPE , allocatable, dimension(:,:) :: localC
    
    rowproc = 2
    colproc = 2
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
    dimA = [5,4]
    
    call cbp_init_local_matrix_descriptor & 
        	(descriptorA,dimA,localdimA,&  
        	rowproc, colproc, rank)  
        	
    allocate(localA(localdimA(1), localdimA(2)))
    allocate(A(dimA(1), dimA(2))) 
    allocate(lineA(dimA(2)))
    print*, descriptorA,'|', localdimA
    
    localA = CZERO
    lineA = rank
    !print*, lineA
    A=CZERO
    
    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
    
    call cbp_setgetprcline(localA, descriptorA, lineA, [dimA(2),1,1,1], 1, rank, nprocs, 'set')
    
    call cbp_setgetprcval(localA, descriptorA, cmplx(7.0,0.0), [1,1], 1, rank, nprocs, 'set')
    
    !print*, "here"
    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,& 
                    	 idxrow, idxcol,rank)
    !print*, "here1"
    if(rank.eq.0) then
    do ii = 1, dimA(1) 
    	print*, A(ii,:)
    end do 
    end if 
    
    deallocate(A)
    deallocate(localA) 
    deallocate(lineA)
    
    print*, "done"
    
    call MPI_FINALIZE(ierr)
end program 
