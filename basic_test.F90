#ifdef CBPCMPLX
#define CTYPE complex
#define MPI_CTYPE MPI_COMPLEX
#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#endif

#ifdef CBPZMPLX 
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#endif

#ifdef CBPREAL 
#define CTYPE real
#define MPI_CTYPE MPI_REAL
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
    

    
    CTYPE , allocatable, dimension(:,:,:) :: localA, localB
    CTYPE , allocatable, dimension(:,:) :: localC
    
    rowproc = 2
    colproc = 2
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
    

    dimA=(/320,240/)
    dimB=(/240,434/)
    dimC=(/dimA(1), dimB(2)/)
    
    if(rank .eq. 0 ) then 
    
    allocate(A(dimA(1), dimA(2)))
    allocate(B(dimB(1), dimB(2)))
    allocate(C(dimC(1), dimC(2)))
    allocate(D(dimC(1), dimC(2)))

    
    !Fill A, B with random numbers

    call random_number(A)
    call random_number(B)

    D = matmul(A,B)
    
    print*, " "
    !print*, ceiling(real(dimA(1))/real(rowproc)), ceiling(real(dimB(2))
    !do ii = 1, rowproc
    !    do jj = 1,colproc
    !        do kk = 1, min(rowproc,colproc)
    !            testr = (/(ii-1)*ceiling(real(dimA(1))/real(rowproc))+1,& 
    !                  ii*ceiling(real(dimA(1))/real(rowproc))/)
    !            testc = (/(jj-1)*ceiling(real(dimB(2))/real(colproc))+1,& 
    !                  jj*ceiling(real(dimB(2))/real(colproc))/)
    !            testf = (/(kk-1)*ceiling(real(dimB(1))/real(min(colproc,rowproc)))+1,& 
    !                  kk*ceiling(real(dimB(1))/real(min(colproc,rowproc)))/)
    !            if(testc(2) .gt. dimB(2)) then
    !                testc(2) = dimB(2)
    !            end if
    !            if(testr(2) .gt. dimA(1)) then 
    !                testr(2) = dimA(1)
    !            end if
    !            if(testf(2) .gt. dimB(1)) then 
    !                testf(2) = dimB(1)
    !            end if
    !    
    !        end do
    !    end do 
    !end do
    
    end if
    
    
    
    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
    call init_local_matrices_descriptors(descriptorA, descriptorB, descriptorC, dimA, dimB, dimC, & 
         localdimA, localdimB, localdimC , rowproc, colproc, process_grid,rank) 
    
    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if

        !-------------------!
    
    

    !print*, "dd" , idxrow, idxcol,   "| ", descriptorA(:,1)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !###############DEBUG##################
    !if(rank .eq. 0 ) then
    !
    !    do ii = 1, rowproc
    !        do jj = 1, colproc
    !         B( (abs(process_grid(ii,jj))-1)*localdimB(1)+1,(jj-1)*localdimB(2)+1 )  = cmplx( descriptorB(1,1), descriptorB(2,1) ) 
            ! if (process_grid(ii,jj) .gt. 0 ) then 
            
            
    !         A( (ii-1)*localdimA(1)+1,(abs(process_grid(ii,jj))-1)*localdimA(2)+1 )  = cmplx( descriptorA(1,1), descriptorA(2,1) ) 
             !!print *, "init",ii, jj,  (ii-1)*localdimA(1)+1, process_grid(ii,jj)*localdimA(2)+1
             !end if
    !        end do
    !    end do 
    
    !end if
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                
    !#####################################
    !And allocate space for the local matrices
    !Doubled the space, in the 3dim, one for computation and recv, the other queued to send.
    
    allocate(localA(localdimA(1),localdimA(2),2))
    allocate(localB(localdimB(1),localdimB(2),2))
    allocate(localC(localdimC(1),localdimC(2)))
    
    localC = cmplx(0.0,0.0)
    localA = cmplx(0.0,0.0)
    localB = cmplx(0.0,0.0)
    
    
    
    !Dispatch global matrix to processes ---------------------------
    
    call dispatch_global_matrix(A, localA(:,:,1), descriptorA(:,1), & 
        rowproc, colproc, idxrow, idxcol, rank, process_grid)
    call dispatch_global_matrix(B, localB(:,:,1), descriptorB(:,1), & 
         rowproc, colproc, idxrow, idxcol, rank, process_grid)
         
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !HERE CALL MATMUL 
    call ComputeMatmul(localA, localB, localC, descriptorA, descriptorB ,& 
                 descriptorC, rowproc, colproc, rank, process_grid)
    
    
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !FILL THE GLOBAL MATRIX
    call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,& 
                     idxrow, idxcol,rank,process_grid)
    
    
    
    !CHECK RESULTS
    
    if(rank .eq. 0 ) then
    print*, sum(C-D), sum(C), sum(D)
    end if
    
    
                
     
    !ENDING..
    deallocate (localA) 
    deallocate (localB) 
    deallocate (localC)
    
    
    if(rank .eq. 0 ) then 
    deallocate (A) 
    deallocate (B) 
    deallocate (C)
    deallocate (D)
    end if
    deallocate (process_grid)
    call MPI_FINALIZE(ierr)
end program 
