#ifdef CBPCMPLX
#define CUSTOM_ABS cabs
#define CTYPE complex
#define typestring  "cmplx"
#define MPI_CTYPE MPI_COMPLEX
#endif

#ifdef CBPDOUBLE 
#define CUSTOM_ABS dabs
#define typestring  "double"
#define CTYPE double precision
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#endif

#ifdef CBPZMPLX 
#define CUSTOM_ABS zabs
#define typestring  "zmplx"
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#endif

#ifdef CBPREAL 
#define typestring  "real"
#define CUSTOM_ABS abs
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#endif





program passm
    use cbp_DistributeMOD
    use mpi
    implicit none
    !### Perform TEST Vars###
    integer :: dims(6), read_stat, iargc
    character*128 :: testname, arg, command, & 
               sparsity_deg, files

    
    ! #### MPI Variables ####
     
    integer :: rank ! the process ID
    integer :: nprocs ! number of processes
    integer , dimension(:,:) , allocatable :: process_grid 
    integer :: rowproc, colproc
    integer :: ierr , idxcol, idxrow, ii, jj,kk
    !Global
    !real , dimension(:,:),  allocatable :: tmpa , tmpb, tmpc
    CTYPE , dimension(:,:) , allocatable :: A, B, C, result
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
    
       
     
    
    print*, " "
    
    if (iargc() .eq. 6) then 
        DO ii = 1, 6
            CALL getarg(ii, arg)
            read(arg, *) dims(ii)
        END DO
              
        write(testname ,"(I10.10, I10.10, I10.10)") dims(1:3)
        write(sparsity_deg, "(I10.10)") dims(4)

        testname=trim(testname)//typestring//sparsity_deg//".bin"
        
    else if(iargc() .eq. 3) then 
        call getarg(1, testname)
        call getarg(2, arg)
        read(arg,*) dims(5) !precision
        call getarg(3, arg)
        read(arg,*)dims(6) !number of gpus
        read(testname(1:10) ,  '(i10)')dims(1)
        read(testname(11:20),  '(i10)')dims(2)
        read(testname(21:30),  '(i10)')dims(3)
        if(typestring .ne. testname(31:len(trim(testname))-14)) then
                print*, "Error type"
                print*, testname(31:len(trim(testname))-14)
                STOP
        end if
        read(testname( len(trim(testname))-13:len(trim(testname))-4 ) , "(i10)")dims(4)
        

    else
        print*, "Error args "
        STOP
     end if


            
    
    dimA=(/dims(1),dims(3)/)
    dimB=(/dims(3),dims(2)/)
    dimC=(/dimA(1),dimB(2)/)


    if(rank .eq. 0) then
    allocate(  A(dimA(1), dimA(2)) )
    allocate(  B(dimB(1), dimB(2)) )
    allocate(  C(dimC(1), dimC(2)) )
    allocate(result(dimC(1),dimC(2)))
        
    command='ls -a tests > temp.txt'
    CALL system(command)
    OPEN(unit=75,file='temp.txt')

    DO WHILE(1.EQ.1)
       READ(75,FMT=*,iostat = read_stat) files
          if(read_stat .lt. 0) then
             print*, trim(testname) &
             , " ERROR, FILE NOT FOUND"
             STOP
          else
             if (trim(files) .eq. trim(testname)) then
               EXIT
             end if
          end if

     END DO
     CLOSE(75)
     CALL system("rm temp.txt")


     OPEN(unit=76,file="./tests/"//trim(files),form='unformatted')
     read(76) A
     read(76) B

    
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
    
    if(rank .eq. 0 ) then 
        deallocate (A) 
        deallocate (B) 
    end if

     
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !HERE CALL MATMUL 
    call ComputeMatmul(localA, localB, localC, descriptorA, descriptorB ,& 
                 descriptorC, rowproc, colproc, rank, process_grid)
    
    
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !FILL THE GLOBAL MATRIX
    call GetBackMatrix(localC, descriptorC, result, rowproc, colproc,& 
                     idxrow, idxcol,rank,process_grid)
    
    
    
    !CHECK RESULTS
    
   

     if (rank .eq. 0 ) then
        read(76) C
        close(76)
        print* , dims , trim(typestring),(maxval(real(CUSTOM_ABS(C-result))))/(maxval(real(CUSTOM_ABS(C)))) ,&
         ((maxval(CUSTOM_ABS(C-result))/(maxval(CUSTOM_ABS(C)))) .le. 10.0**(-1.0*dims(5))), &
        10.0**(-1.0*dims(5)), sum(C), sum(result)

        deallocate(C)
        deallocate(result)
     end if
    
    
                
     
    !ENDING..
    deallocate (localA) 
    deallocate (localB) 
    deallocate (localC)
    
    
    deallocate (process_grid)
    call MPI_FINALIZE(ierr)

end program 
