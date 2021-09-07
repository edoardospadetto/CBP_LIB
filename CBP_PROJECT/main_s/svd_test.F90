#ifdef CBPCMPLX
#define CTYPE complex
#define MPI_CTYPE MPI_COMPLEX
#define CZERO cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#define CZERO 0.0
#define CTYPE1  1.0
#endif

#ifdef CBPZMPLX
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#define CZERO dcmplx(0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#endif

#ifdef CBPREAL
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#define CZERO 0.0
#define CTYPE1 1.0
#endif


program passm
    use cbp_DistributeMOD
    use mpi
    use lapack_interface


    implicit none
    ! #### MPI Variables ####
     CTYPE , dimension(:,:) , allocatable :: A,A0,R,L, U,VT

    integer :: rank, idxrow, idxcol ! the process ID
    integer :: nprocs ! number of processes
    integer , dimension(:,:) , allocatable :: process_grid
    integer :: rowproc, colproc
    integer :: ierr ,ii, jj,kk
    !Global
    character*16 :: num1char, num2char
    integer,dimension(2) :: ld1, ld2, ld3
    integer :: dimA(2), dimL(2), dimR(2), descriptorSV(2)

    integer , dimension(4):: descriptorA,descriptorL, descriptorR, descriptorA0
    integer , dimension(4):: temp1, temp2, temp3

    integer :: localdimA(2), localdimL(2), localdimR(2)


    double precision :: start, end, flops

    integer :: zeros, zerosall, expected

    type(c_ptr) :: context
   double precision, dimension(:), allocatable :: D
    CTYPE, dimension(:), allocatable :: singular_values
    CTYPE :: trace, trace2
#ifdef CBPCMPLX
     real , dimension (:,:), allocatable  ::  rai, rar
#endif
#ifdef CBPZMPLX
     double  precision , dimension(:,:) , allocatable :: rar, rai
#endif


    CTYPE , allocatable, dimension(:,:) :: localA, localL, localQ, localR, localQinv, localQres, localA0

  	CALL getarg(1,num1char)   !first, read in the two values
	CALL getarg(2,num2char)

	READ(num1char,*)rowproc                   !then, convert them to REALs
	READ(num2char,*)colproc

   ! rowproc = 1
   ! colproc = 2


    call MPI_INIT(ierr)


    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    OPEN(unit=1111, file="svd_precision_ZMPLX.txt")
do kk = 1,50
    dimA=(/50,50/)
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

    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
    allocate(process_grid(rowproc,colproc))

    !------------------- Local Matrices


    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localL(localdimL(1),localdimL(2)))
    allocate(localR(localdimR(1),localdimR(2)))

    allocate(localA0(localdimA(1),localdimA(2)))

    !allocate(localQ(localdimQ(1),localdimQ(2)))
    !allocate(localQinv(localdimQ(1),localdimQ(2)))
    !allocate(localQres(localdimQ(1),localdimQ(2)))

    localA= CZERO
    localL = CZERO
    localR = CZERO


    !-------------------! Global Matrices



        allocate(A(dimA(1), dimA(2))) ! good one
        allocate(A0(dimA(1), dimA(2)))
        allocate(L(dimL(1) , dimL(1)))
        allocate(R(dimR(1), dimR(2)))




#if defined(CBPCMPLX) |defined (CBPZMPLX)
    allocate(rar(dimA(1), dimA(2)))
    allocate(rai(dimA(1), dimA(2)))


   call random_number(rar)
   call random_number(rai)


#ifdef CBPCMPLX
   A = cmplx(rar, rai)
#endif
#ifdef CBPZMPLX
   A = dcmplx(rar, rai)
#endif
   call dispatch_global_matrix(A, localA, descriptorA, rowproc, colproc,&
        		 idxrow, idxcol,rank)
   deallocate(rar)
   deallocate(rai)
#endif


#if defined(CBPREAL) | defined(CBPDOUBLE)
    !Fill A, B with random numbers

    call random_number(localA)

#endif

    ! R = A
    if (rank.eq.0) then
    print*, "start"
    end if
    locala = (locala- 0.5)*14.0


    call init_cublasXT(1,context)
    if(rank .eq. 0) then
            start = MPI_Wtime();
    endif



    !Needed for get back matrix old implementation sorry :(
    call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
    call GetBackMatrix(localA, descriptorA, A0, rowproc, colproc,&
         	  	idxrow, idxcol,rank)

    call cbp_Bidiag(localA, descriptorA, localR, descriptorR, localL, &
    			descriptorL, dimA ,rowproc,colproc, rank, nprocs, context)

     if (rank.eq.0) then
     print*, "done bidiag"
     end if

     call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,&
         	  	idxrow, idxcol,rank)

    allocate(VT(minval(dimA), dimA(2)) )
    allocate(U(dimA(1), minval(dimA))  )
    allocate(D(minval(dimA))  )

    call SVD_decompose_lapack(A0,U,D,VT)

    deallocate(U)
    deallocate(VT)
    allocate (singular_values(cbp_singular_value_array_size(dimA, descriptorA, rank, nprocs)))


    localL = CZERO
    call cbp_SVD_GulobReinsch( localA, descriptorA, dimA, singular_values, descriptorSV, rowproc, colproc,rank, nprocs)
    trace2 = cbp_sum(singular_values**2,1,descriptorSV(1),rank,nprocs)



    if (rank.eq.0) then
      A = matmul(A,conjg(transpose(A)))
      trace = CZERO
      do ii = 1, dimA(1)

    		trace = A(ii,ii) + trace !print "(*('('sf8.2xspf8.2x'i)':x))", L(ii, :)

     end do
    write(1111,*,advance="no") , abs(trace-trace2), abs((sum(D)-trace)) , abs((D-singular_values))
    end if

    deallocate(process_grid)
    deallocate(singular_values)
    deallocate(d)
    deallocate(localA)
    deallocate(localR)
    deallocate(localL)
    deallocate(localA0)
    deallocate(A)
    deallocate(A0)
    deallocate(R)
    deallocate(L)
    !end if
    end do
    close(1111)
    call MPI_FINALIZE( ierr )

end program
