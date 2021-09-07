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
    CTYPE , dimension(:,:) , allocatable :: A,A0, B, C,D
    integer :: dimA(2), dimB(2), dimC(2)
    
    integer , dimension(4):: descriptorA, descriptorA0, descriptorB
    integer :: descriptorC(4)
    integer :: localdimA(2), localdimB(2), localdimC(2)
    integer :: testr(2), testc(2) , testf(2), tempdim

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
    
    rowproc = 3
    colproc = 2
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
do kk = 0,0
    dimA=(/321,247/)
    dimB=(/247,434/)
    dimC=(/dimA(1), dimB(2)/)
   
    dimA = kk*2000 + dimA
    dimB = kk*2000 + dimB
    dimC = kk*2000 + dimC

 
    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    
    call cbp_init_local_matrix_descriptor(descriptorA0, [dimA(2), dimA(1)], localdimA, rowproc, colproc, rank)  
    call cbp_init_local_matrix_descriptor(descriptorB, dimB, localdimB, rowproc, colproc, rank) 
    
  
    
    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
    call cbp_matmul_descriptors(descriptorA, descriptorB, descriptorC, dimA, dimB, dimC, & 
         localdimA, localdimB, localdimC , rowproc, colproc, process_grid,rank) 
    
   
    
    if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs 
        STOP 
    end if

    !-------------------!
        
    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localB(localdimB(1),localdimB(2)))
    allocate(localC(localdimC(1),localdimC(2)))
    
   
    
    localC = CZERO
    localA = CZERO
    localB = CZERO
    
     
    !-------------------!

    if(rank .eq. 0 ) then 
        print*, "start"
        allocate(A(dimA(1), dimA(2))) ! good one
        allocate(A0(dimA(2), dimA(1)))  !need to be tranposed
        allocate(B(dimB(1), dimB(2)))
        allocate(C(dimC(1), dimC(2)))
        allocate(D(dimC(1), dimC(2)))

    end if


#if defined(CBPCMPLX) |defined (CBPCMPLX) 
    allocate(rar(localdimA(1), localdimA(2)))
    allocate(rai(localdimA(1), localdimA(2)))
    allocate(rbr(localdimB(1), localdimB(2)))
    allocate(rbi(localdimB(1), localdimB(2)))

   call random_number(rar)
   call random_number(rai)
   call random_number(rbr)
   call random_number(rbi)
#ifdef CBPCMPLX       
   localA = cmplx(rar, rai)
   localB = cmplx(rbr, rbi)
#endif
#ifdef CBPZMPLX
   localA = dcmplx(rar, rai)
   localB = dcmplx(rbr, rbi)
#endif
   deallocate(rbr)
   deallocate(rar)
   deallocate(rai) 
   deallocate(rbi)
#endif


#if defined(CBPREAL) | defined(CBPDOUBLE)
    !Fill A, B with random numbers

    call random_number(localA)
    call random_number(localB)
    !do ii = 1, localdimA(2)
    !    localA(:,ii,1) = descriptorA0(4,1)+ii-1
    !end do 

#endif
  
    !GET BACK INPUT MATRICES
    A0 = CZERO
    
    call GetBackMatrix(localA, descriptorA0, A0, rowproc, colproc,& 
                     idxrow, idxcol,rank) 
    call GetBackMatrix(localB, descriptorB, B, rowproc, colproc,& 
                     idxrow, idxcol,rank)  
   
    ! Transpose A
    if(rank .eq. 0) then 
            start = MPI_Wtime();
    endif

    call cbp_transpose(localA, descriptorA0) !Mat, olddescriptor, newdescriptor
    
    call cbp_redispatch(localA, descriptorA0,descriptorA,rank,nprocs) 
    
    if(rank .eq. 0) then 
            end = MPI_Wtime();
    endif




    !GET BACK INPUT MATRICES
    
    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,& 
                     idxrow, idxcol,rank)

    
    if (rank .eq. 0) then  
                print*,  'time = ', end- start , 'zero is good = ' , sum(A-transpose(A0)), sum(A), sum(A0),& 
                 sum(transpose(A0))  , shape(A0) , shape(A)
    end if 
   
   

    if (rank .eq. 0) then
        D = matmul(transpose(A0),B) 
     
    print*, "done gen"
    end if
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
  
    do ii = 1, 1
     
         do jj = 0,0
      
    if (rank .eq. 0) then 
       
       start = MPI_Wtime();
    endif
    call init_cublasXT(ii,context)
     

     call cbp_matmul(localA, localB, localC, descriptorA, descriptorB ,& 
                 descriptorC, rowproc, colproc, rank, process_grid,  context)


     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call destroy_cublasXt(context)
   
    if(rank .eq. 0) then 
            end = MPI_Wtime();
    endif

    !FILL THE GLOBAL MATRIX
    call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,& 
                     idxrow, idxcol,rank)
    
    !CHECK RESULTS
    if(rank .eq. 0 ) then
            print*, maxval(abs(C-D))*dble(dimA(1))*dble(dimA(2))/sum(C), sum(C) , sum(D) , & 
                   end-start, 2*dble(dimA(2))*dble(dimA(1))*dble(dimB(2)) / (end-start), & 
            ii, jj*0.01, kk*dble(2000)
        
            C=0.0
    end if    
    
    localC = 0.0  
    
    end do
    
    end do
   
     
     
    !ENDING..
    deallocate (localA) 
    deallocate (localB) 
    deallocate (localC)
    
    
    if(rank .eq. 0 ) then 
    deallocate(A0)
    deallocate (A) 
    deallocate (B) 
    deallocate (C)
    deallocate (D)
    end if
    deallocate (process_grid)
 end do 
    call MPI_FINALIZE(ierr)
end program 
