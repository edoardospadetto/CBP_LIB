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
    integer :: ierr , idxcol, idxrow, ii, jj,kk,ngpu
    !Global
    character*8 :: args(4)
    integer :: dimA(2), dimB(2), dimC(2), stat,op

    integer , dimension(2):: descriptorv

    integer :: localdimA(2), localdimB(2), localdimC(2)


    double precision ::  flops, t1, t2

    type(c_ptr) :: context


    CTYPE , allocatable, dimension(:) :: v
    double precision :: vnorm

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

    dimv = kk*100
    call cbp_init_1D_descriptor(descriptorv, dimv, rank, nprocs)
    allocate(v(ceiling(dimv/nprocs)))

    v = CTYPE0

	  call c_random(V)

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)




    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (rank .eq. 0)  t1 = MPI_Wtime()

    vnorm = cbp_l2norm(v, 1, descriptorv1 , rank, nprocs)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) t2 = MPI_Wtime()

    deallocate(v)



    !FILL THE GLOBAL MATRIX
    !call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,&
     !                idxrow, idxcol,rank,process_grid)

    !CHECK RESULTS
    if(rank .eq. 0 ) print*,ii, dimv(1), t2-t1







    !ENDING..
    deallocate (localA)
    deallocate (localB)
    deallocate (localC)


    deallocate (process_grid)
 end do

    call MPI_FINALIZE(ierr)
end program
