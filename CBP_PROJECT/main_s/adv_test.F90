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
    integer :: dimA(2), dimB(2), dimC(2)
    
    integer , dimension(4,2):: descriptorA(4,2), descriptorB(4,2)
    integer :: descriptorC(4)
    integer :: localdimA(2), localdimB(2), localdimC(2)
    integer :: testr(2), testc(2) , testf(2)

    double precision :: start, end, flops, t1, t2 
    
    type(c_ptr) :: context
    
#ifdef CBPCMPLX
        real , dimension (:,:), allocatable  ::  rar, rbr, rai, rbi
#endif
#ifdef CBPZMPLX
        double  precision , dimension(:,:) , allocatable :: rar, rbr, rai, rbi 
#endif

    
    CTYPE , allocatable, dimension(:,:,:) :: localA, localB
    CTYPE , allocatable, dimension(:,:) :: localC
    
    rowproc = 4
    colproc = 4
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
do kk = 1,25
    dimA=(/320,240/)
    dimB=(/240,434/)
    dimC=(/dimA(1), dimB(2)/)
   
    dimA = kk*20000 + dimA
    dimB = kk*20000 + dimB
    dimC = kk*20000 + dimC



 
    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
    call cbp_matmul_descriptors(descriptorA, descriptorB, descriptorC, dimA, dimB, dimC, & 
         localdimA, localdimB, localdimC , rowproc, colproc, process_grid,rank) 
    
    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if

    !-------------------!

    allocate(localA(localdimA(1),localdimA(2),2))
    allocate(localB(localdimB(1),localdimB(2),2))
    allocate(localC(localdimC(1),localdimC(2)))
    
    localC = CZERO
    localA = CZERO
    localB = CZERO
    
     
    !-------------------!

    if(rank .eq. 0 ) then 
        print*, "start"
    
        !allocate(A(dimA(1), dimA(2)))
        !allocate(B(dimB(1), dimB(2)))
        !allocate(C(dimC(1), dimC(2)))
        !allocate(D(dimC(1), dimC(2)))

    end if

#ifdef CBPCMPLX  
    allocate(rar(localdimA(1), localdimA(2)))
    allocate(rai(localdimA(1), localdimA(2)))
    allocate(rbr(localdimB(1), localdimB(2)))
    allocate(rbi(localdimB(1), localdimB(2)))

   call random_number(rar)
   call random_number(rai)
   call random_number(rbr)
   call random_number(rbi)
        
   localA(:,:,1) = cmplx(rar, rai)
   localB(:,:,1) = cmplx(rbr, rbi)

   deallocate(rbr)
   deallocate(rar)
   deallocate(rai) 
   deallocate(rbi)
#endif
#ifdef CBPZMPLX
    allocate(rar(localdimA(1), localdimA(2)))
    allocate(rai(localdimA(1), localdimA(2)))
    allocate(rbr(localdimB(1), localdimB(2)))
    allocate(rbi(localdimB(1), localdimB(2)))
    
   call random_number(rar)
   call random_number(rai)
   call random_number(rbr)
   call random_number(rbi)
        
   localA(:,:,1) = dcmplx(rar, rai)
   localB(:,:,1) = dcmplx(rbr, rbi)

   deallocate(rbr)
   deallocate(rar)
   deallocate(rai) 
   deallocate(rbi)
#endif
#ifdef CBPREAL
    !Fill A, B with random numbers

    call random_number(localA(:,:,1))
    call random_number(localB(:,:,1))

#endif
#ifdef CBPDOUBLE
    !Fill A, B with random numbers

    call random_number(localA(:,:,1))
    call random_number(localB(:,:,1))

#endif   
    
    !GET BACK INPUT MATRICES
    !call GetBackMatrix(localA(:,:,1), descriptorA, A, rowproc, colproc,& 
    !                 idxrow, idxcol,rank,process_grid) 
    !call GetBackMatrix(localB(:,:,1), descriptorB, B, rowproc, colproc,& 
    !                 idxrow, idxcol,rank,process_grid)  
    
    if (rank .eq. 0) then
        !D = matmul(A,B) 
   
    print*, "done gen"
    end if
    
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
  
    do ii = 1, 4
        do jj = 0,0
        
    if (rank .eq. 0) then 
       start = MPI_Wtime();
    endif
    call init_cublasXT(ii,context)
    if (rank .eq. 0) then 
       t1 = MPI_Wtime();
    endif

     call ComputeMatmul(localA, localB, localC, descriptorA, descriptorB ,& 
                 descriptorC, rowproc, colproc, rank, process_grid, ii, jj*0.001, context)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) then 
       t2 = MPI_Wtime();
    endif
 

     call destroy_cublasXt(context)
   
    if(rank .eq. 0) then 
            end = MPI_Wtime();
    endif

    !FILL THE GLOBAL MATRIX
    !call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,& 
     !                idxrow, idxcol,rank,process_grid) 
    
    !CHECK RESULTS
    if(rank .eq. 0 ) then
            !print*, maxval(abs(C-D))*dble(dimA(1))*dble(dimA(2))/sum(C), sum(C) , sum(D) , & 
            print*,        t1-start, t2-t1, end-t2, 2*dble(dimA(2))*dble(dimA(1))*dble(dimB(2)) / ((end-start)*(10**9)), & 
            ii, jj*0.01,kk,  kk*dble(20000)
        
            !C=0.0
    end if    
    
    localC = 0.0  
    
    end do
    
    end do
   
     
     
    !ENDING..
    deallocate (localA) 
    deallocate (localB) 
    deallocate (localC)
    
    
    if(rank .eq. 0 ) then 
    !deallocate (A) 
    !deallocate (B) 
    !deallocate (C)
    !deallocate (D)
    end if
    deallocate (process_grid)
 end do 
    call MPI_FINALIZE(ierr)
end program 
