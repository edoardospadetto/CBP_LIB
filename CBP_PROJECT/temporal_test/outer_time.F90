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
    integer , dimension(2):: descriptorv, descriptoru
    integer :: localdimA(2), localdimB(2), localdimC(2)


    double precision ::  flops, t1, t2,t3

    type(c_ptr) :: context


    CTYPE , allocatable, dimension(:,:) :: localA, localB, A
    CTYPE , allocatable, dimension(:,:) :: localC
    CTYPE , allocatable, dimension(:) :: localv, localu , res1 , res2, v

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



do kk = 1,20,1
  dimC=[500*kk,500*kk]
  dimA=[dimC(1),1]
  dimB=[1,dimC(2)]





    !-----------------------------------------------------------------
    !DIVIDE MATRIX IN THE PROCESS GRID
    allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))

    call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)


  call cbp_init_1D_descriptor(descriptoru, dimA(1), rank, nprocs)
  call cbp_init_1D_descriptor(descriptorv, dimA(2), rank, nprocs)
    !-------------------!

	 call cbp_matmul_descriptors(descriptorA, descriptorB, descriptorC, dimA, dimB, dimC, &
     localdimA, localdimB, localdimC , rowproc, colproc, process_grid,rank)


    allocate(localB(localdimB(1),localdimB(2)))
    allocate(localC(localdimC(1),localdimC(2)))
    allocate(localA(localdimA(1),localdimA(2)))
    allocate(localv(localdimA(2)))
    allocate(localu(localdimA(1)))

    localC = CZERO
    localA = CZERO
    localB = CZERO



    allocate(A(dimA(1), dimA(2)))
    allocate(v(dimA(2)))
    allocate(res1(dimA(1)))
    allocate(res2(dimA(1)))


	call c_random(localA)
	call c_random(localB)
  call c_random(localv)
  call c_random(localu)




    call MPI_BARRIER(MPI_COMM_WORLD, ierr)


    do ii = 1, 1




    call init_cublasXT(ngpu,context)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (rank .eq. 0)  t1 = MPI_Wtime()
         if(op .eq. 1) then

        call cbp_init_local_matrix_descriptor(descriptorA0,dimA,&
                localdimA,rowproc, colproc, rank)
        call cbp_transpose(localA, descriptorA0)

         end if
         call cbp_redispatch(localA,descriptorA0, descriptorA, rank, nprocs)
        call cbp_matmul(localA, localB, localC, descriptorA, descriptorB ,&
                 descriptorC, rowproc, colproc, rank, process_grid,context)

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if (rank .eq. 0) t2 = MPI_Wtime()

       call cbp_outer(localv, descriptorv, localu, descriptoru, localC, descriptorC, dimC, 0, rowproc, colproc, rank, nprocs, context  )


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) t3 = MPI_Wtime()




    !res2 = localu(:descriptorres(1))
    !call cbp_redispatch_1d(res2,descriptorres,[dimA(1),1],rank,nprocs)


     !res1 = matmul(A,v)


     call destroy_cublasXt(context)



    !FILL THE GLOBAL MATRIX
    !call GetBackMatrix(localC, descriptorC, C, rowproc, colproc,&
     !                idxrow, idxcol,rank,process_grid)

    !CHECK RESULTS
    if(rank .eq. 0 ) print*, dimA(1), t2-t1, t3-t2!, sum(abs(res1-res2))/dimA(1), sum(abs(res1)), sum(abs(res2)), sum(v)
    localC = 0.0

    end do




    !ENDING..
    deallocate(A)
    deallocate(res1)
    deallocate(res2)
    deallocate(v)
    deallocate(localu)
    deallocate(localv)
    deallocate (localA)
    deallocate (localB)
    deallocate (localC)


    deallocate (process_grid)
 end do

    call MPI_FINALIZE(ierr)
end program
