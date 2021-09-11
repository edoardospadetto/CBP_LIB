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
    integer , dimension(:,:) , allocatable :: process_grid
    integer :: rowproc, colproc
    integer :: ierr , idxcol, idxrow, ii, jj,kk
    !Global
    !real , dimension(:,:),  allocatable :: tmpa , tmpb, tmpc
    CTYPE , dimension(:,:) , allocatable :: A, B, C,D
    integer :: dimA(2), dimB(2), dimC(2)
    character*8 , dimension(3) :: args
    integer , dimension(4):: descriptorA(4), descriptorB(4)
    integer :: descriptorC(4), ngpu, stat
    integer :: localdimA(2), localdimB(2), localdimC(2)
    integer :: testr(2), testc(2) , testf(2)

    double precision :: start, end, flops

    type(c_ptr) :: context


    CTYPE , allocatable, dimension(:,:) :: localA, localB
    CTYPE , allocatable, dimension(:,:) :: localC



    call MPI_INIT(ierr)


    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)


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

      print*, "rowproc", args(1), "colproc", args(2), "#gpu", args(3)

do kk = 1,20
    dimA=(/0,0/)
    dimB=(/0,0/)
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


        allocate(A(dimA(1), dimA(2)))
        allocate(B(dimB(1), dimB(2)))
        allocate(C(dimC(1), dimC(2)))
        allocate(D(dimC(1), dimC(2)))

    end if

    do ii = 1, descriptorA(1)
        do jj = 1, descriptorA(2)
          !  print*, "in",0.5**(descriptorA(3)+jj-2)
          localA(ii,jj) = 1.01010101 !1e5*0.5**(descriptorA(4)+jj-2)
        end do
    end do
    do ii = 1, descriptorB(1)
        do jj = 1, descriptorB(2)
          localB(ii,jj) = 1.01010101 !  1e5*0.5**(descriptorB(3)+ii-2)
        end do
    end do

    !print*, "gen", sum(localA), sum(localB), dimA(1)*(1-(0.5)**dimA(2))/(1-(0.5))
    !GET BACK INPUT MATRICES
    call GetBackMatrix(localA, descriptorA, A, rowproc, colproc,&
                     idxrow, idxcol,rank)
    call GetBackMatrix(localB, descriptorB, B, rowproc, colproc,&
                     idxrow, idxcol,rank)

    if (rank .eq. 0) then
        D = matmul(A,B)

    end if


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)



    call init_cublasXT(1,context)

    call cbp_matmul(localA, localB, localC, descriptorA, descriptorB ,&
                 descriptorC, rowproc, colproc, rank, process_grid,context)


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
            print*, nprocs,  abs(maxval(abs(C))), abs(maxval(abs(D)))
    end if

    !GET BACK INPUT MATRICES
    A = CZERO
    B = CZERO

    localC = 0.0




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
