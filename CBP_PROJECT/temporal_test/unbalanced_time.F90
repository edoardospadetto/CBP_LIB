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
    integer :: ierr , idxcol, idxrow, ii, jj,kk,ngpu
    !Global
    character*8 :: args(4)
    integer :: dimA(2), dimB(2), dimC(2), stat,op

    integer , dimension(4):: descriptorA,descriptorA0, descriptorB, descriptorC

    integer :: localdimA(2), localdimB(2), localdimC(2)


    double precision ::  flops, t1, t2

    type(c_ptr) :: context


    CTYPE , allocatable, dimension(:,:) :: localA, B
    CTYPE , allocatable, dimension(:,:) :: localC

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
	    read(args(3),*,iostat=stat) op
	    if (stat .ne. 0 ) then
	    	print *, "error"
	    	STOP
	    end if

	    call get_command_argument(4,args(4))
	    read(args(4),*,iostat=stat)  ngpu
	    if (stat .ne. 0 ) then
	    	print *, "error"
	    	STOP
	    end if

	!----------------------------------------------------------

 	 if(nprocs .ne. rowproc*colproc) then
        print *, "invalid process grid, fatal error ", nprocs
        STOP
    end if

do kk = 1,200,10
    dimA=(/500*kk,10/)
    dimB=(/10,10/)
    dimC=(/500*kk, 10/)




    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))

    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)




    !-------------------!

    call cbp_unbalanced_matmul_descriptors(descriptorA, dimA(1), dimA(2), rank, nprocs)


    allocate(B(localdimB(1),localdimB(2)))
    allocate(localC(localdimC(1),localdimC(2)))
    allocate(localA(localdimA(1),localdimA(2)))

    localC = CZERO
    localA = CZERO
    localB = CZERO


	call c_random(localA)
	call c_random(B)

  if ( rank.eq.0) then
			call MPI_BCAST(B,product(dimB), MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr )
	else
			call MPI_BCAST(B,product(dimB), MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr )
	end if



    call MPI_BARRIER(MPI_COMM_WORLD, ierr)


    do ii = 1, ngpu




    call init_cublasXT(ii,context)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (rank .eq. 0)  t1 = MPI_Wtime()
         if(op .eq. 1) then

        call cbp_init_local_matrix_descriptor(descriptorA0,dimA,&
                localdimA,rowproc, colproc, rank)
        call cbp_transpose(localA, descriptorA0)
        call cbp_redispatch(localA,descriptorA0, descriptorA, rank, nprocs)
        end if

       call cbp_unbalanced_matmul(localA, descriptorA, B, localC, descriptorC, context)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) t2 = MPI_Wtime()




     call destroy_cublasXt(context)



    !FILL THE GLOBAL MATRIX
    !call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,&
     !                idxrow, idxcol,rank,process_grid)

    !CHECK RESULTS
    if(rank .eq. 0 ) print*,ii, dimA(1), t2-t1
    localC = 0.0

    end do





    !ENDING..
    deallocate (localA)
    deallocate (localB)
    deallocate (localC)


    deallocate (process_grid)
 end do

    call MPI_FINALIZE(ierr)
end program
