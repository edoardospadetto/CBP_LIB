#define N_TEMPORARY_BYTES 200000000000
#define N_PERMANENT_BYTES 2000
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
    integer :: ierr , idxcol, idxrow, ii, jj,kk, mm,ll
    !Global
    !real , dimension(:,:),  allocatable :: tmpa , tmpb, tmpc
    CTYPE , dimension(:,:) , allocatable :: A, B, C,D
    integer :: dimA(2), dimB(2), dimC(2)
    
    integer , dimension(4):: descriptorA, descriptorB
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
    colproc = 3
    
    
    call MPI_INIT(ierr)
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
 
do kk = 2,2
    dimA=(/321,247/)
    dimB=(/247,434/)
    dimC=(/dimA(1), dimB(2)/)
   
    dimA = kk*2000 + dimA
    dimB = kk*2000 + dimB
    dimC = kk*2000 + dimC



 
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

    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localB(localdimB(1),localdimB(2)))
    allocate(localC(localdimC(1),localdimC(2)))
    
    localC = CZERO
    localA = CZERO
    localB = CZERO
    
     
    !-------------------!

    if(rank .eq. 0 ) then 
        print*, "start"
    
        allocate(A(dimA(1), dimA(2)))
        allocate(B(dimB(1), dimB(2)))
        allocate(C(dimC(1), dimC(2)))
        allocate(D(dimA(1), dimA(2)))

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
        
   localA = cmplx(rar, rai)
   localB = cmplx(rbr, rbi)

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
        
   localA = dcmplx(rar, rai)
   localB = dcmplx(rbr, rbi)

   deallocate(rbr)
   deallocate(rar)
   deallocate(rai) 
   deallocate(rbi)
#endif
#ifdef CBPREAL
    !Fill A, B with random numbers

    call random_number(localA)
    call random_number(localB)
    localB = 1.0
    localA = 2.0
#endif
#ifdef CBPDOUBLE
    !Fill A, B with random numbers

    call random_number(localA)
    call random_number(localB)

#endif   
    
    !GET BACK INPUT MATRICES
    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,& 
                     idxrow, idxcol,rank) 
    call GetBackMatrix(localB, descriptorB, B, rowproc, colproc,& 
                     idxrow, idxcol,rank)  
        
    print*, descriptorA
    if (rank .eq. 0) then
        !D = matmul(A,B) 
   
    print*, "done gen"
    end if
    
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
  
    do ii = 1, 4
       if(rank .eq. 0) then
                print*, "ii"    
        end if
         do jj = 0,0
      
    if (rank .eq. 0) then 
       
       start = MPI_Wtime();
    endif
    call init_cublasXT(ii,context)
     
    call ComputeMatmul(localA, localB, localC, descriptorA, descriptorB ,& 
                 descriptorC, rowproc, colproc, rank, process_grid, ii, 0.00, context)

     
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     call destroy_cublasXt(context)
   
    if(rank .eq. 0) then 
            end = MPI_Wtime();
    endif

    !FILL THE GLOBAL MATRIX
    call GetBackMatrix(localA, descriptorA, D, rowproc, colproc,& 
                     idxrow, idxcol,rank) 
   
    !print*, rank, descriptorA         
    !CHECK RESULTS

    if (rank.eq.0) then
    !do ll = 1, dimA(2)
    !    do mm = 1, dimA(1) 
    !            if( abs(A(mm,ll) - D(mm,ll)) .gt. 1e-4) then   
    !            print*, mm, ll, A(mm,ll), D(mm,ll) 
    !            end if
    !    end do 
    !end do 
    end if
    if(rank .eq. 0 ) then
            print*, maxval(abs(A-D))*dble(dimB(1))*dble(dimB(2)), sum(B) , sum(D) , & 
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
    deallocate (A) 
    deallocate (B) 
    deallocate (C)
    deallocate (D)
    end if
    deallocate (process_grid)
 end do 
    call MPI_FINALIZE(ierr)
end program 
