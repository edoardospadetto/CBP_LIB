#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX
#define CUSPARSETRIDIAG cusparse_cgtsv
#define CUSPARSEBUFFER  cusparseCgtsv2_bufferSizeExt
#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define CTYPE0 0.0
#define CTYPE1 1.0
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#define CUSPARSETRIDIAG cusparse_dgtsv
#define CUSPARSEBUFFER  cusparseDgtsv2_bufferSizeExt
#endif

#ifdef CBPZMPLX
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#define CUSPARSETRIDIAG cusparse_zgtsv
#define CUSPARSEBUFFER  cusparseZgtsv2_bufferSizeExt
#endif

#ifdef CBPREAL
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#define CUSPARSETRIDIAG cusparse_sgtsv
#define CUSPARSEBUFFER  cusparseSgtsv2_bufferSizeExt
#endif


module cbp_DistributeMOD
    use mpi
    use cbpfor_parallel

    implicit none
    contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------BASIC------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!Find first value equal to the requested in an array, if not found resturns -1
    function findval(array, val) result(ii)
    integer , dimension(:) :: array
    integer :: val,ii
    logical :: good

    good = .false.
    do ii = 1, size(array)
        if(array(ii) .eq. val) then
        good = .true.
        return
        exit
        end if
    end do

    ii = -1

    return


    end function


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------DEBUG------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !Just print a local Matrix on terminal
        SUBROUTINE printmatrix(mat)
          IMPLICIT NONE
          INTEGER :: ii,jj
          INTEGER, DIMENSION(:,:), INTENT(IN) :: mat
          INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
          sizem=SHAPE(mat)

          PRINT*, ""
          DO ii=1,sizem(1)
      		print*, mat(ii,:)

          END DO


        END SUBROUTINE printmatrix

  !Check error of mpi
    subroutine mpic_checkerror(ierr)
    implicit none
    character(512) ::what
    integer :: ierr, ierr0,ll
    if(ierr .ne. 0 ) then
    call MPI_ERROR_STRING(ierr, what, ll, ierr0)
    print*, trim(what)
    end if
    end subroutine

    !=========================================================================================

!Process #0 contains the global Matrix, request the descriptor from all the others, and according to it
!send back the correct chunk of the global matrix



subroutine dispatch_global_matrix(A, localA, descriptor, rowproc, colproc,&
         idxrow, idxcol,rank)
        use mpi
        implicit none
        !Use the descriptor and the process grid slepy dumb idiot.
        CTYPE, dimension(:,:) :: A
        CTYPE, dimension(:,:) :: localA
        integer, dimension(4) :: descriptor, extdescriptor
        integer :: rowproc, colproc, idxrow, idxcol, rank
        integer, dimension(MPI_STATUS_SIZE) :: status_

        integer :: ierr, ii,jj


        if(rank .eq. 0) then

        do ii = 1, rowproc
            do jj = 1, colproc
                if(ii*jj .ne. 1) then

                    call MPI_RECV(extdescriptor, 4, MPI_INTEGER, &
                              grid2rank(ii,jj,rowproc,colproc) &
                             ,3*(ii*(colproc+1)+jj)+1 , MPI_COMM_WORLD, status_, ierr)
                    call mpic_checkerror(ierr)

                    if(sum(extdescriptor) .ne. 0 ) then

                    call MPI_SSEND(A(extdescriptor(3):extdescriptor(3)+extdescriptor(1)-1 ,&
                            extdescriptor(4):extdescriptor(4)+extdescriptor(2)-1) , extdescriptor(1)*extdescriptor(2),&
                         MPI_CTYPE, grid2rank(ii,jj,rowproc,colproc) ,&
                         ii*(colproc+1)+jj , MPI_COMM_WORLD, ierr)

                    call mpic_checkerror(ierr)

                    end if
                else
                    localA = A(descriptor(3): descriptor(3)+descriptor(1)-1,&
                           descriptor(4):descriptor(4)+descriptor(2)-1)
                end if

            end do
        end do
        else




        call MPI_SSEND( descriptor , 4,&
                     MPI_INTEGER, 0 ,&
                     3*(idxrow*(colproc+1)+idxcol)+1 , MPI_COMM_WORLD, ierr)

        call mpic_checkerror(ierr)

        if(sum(descriptor) .ne. 0 ) then

        call MPI_RECV(localA(:descriptor(1),:descriptor(2)), descriptor(1)*descriptor(2), MPI_CTYPE, 0 &
            ,idxrow*(colproc+1)+idxcol , MPI_COMM_WORLD, status_, ierr)

        call mpic_checkerror(ierr)



        end if


        end if






end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!--------------------- DUMB FUNCTIONS FOR THE IMPLEMENTED ALGEBRA-----------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    !Computes the rank of the process block with coordinate merow+down, mecol+side
    !If the process grid pass the limits it restart from the other side
    !Pacman Effect :3
    function next2me(merow, mecol, down,side, totrows, totcols) result (idx)
        implicit none
        integer :: merow, mecol, side, down, totrows, totcols, itcol, itrow
        integer :: idx

        itcol = (mecol+side)
        itrow = (merow+down)

        if (itrow .gt. totrows)  then
            itrow = mod(itrow-1,totrows)+1
        else if( itrow .le. 0  ) then
            itrow = totrows- mod(itrow, totrows)
        end if

        if (itcol .gt. totcols)then
            itcol = mod(itcol-1,totcols)+1
        else if( itcol .le. 0  ) then
            itcol = totcols -mod(itcol, totcols)
        end if

        !get to the rank
        !print*, merow,mecol, "||", down, side,"||" ,itrow, itcol
        idx = (itcol-1)*totrows + itrow-1

        return

    end function

!Next step of the process grid to perform matmul

    subroutine update_process_grid(process_grid, rowproc, colproc, rank)
            use mpi
        implicit none
        integer , dimension(:,:) :: process_grid
        integer :: ld, fromc, jj, colproc, rowproc, rank
        integer, dimension(:), allocatable :: temp_
        integer ::ierr


        ld = size(process_grid, 1)

        if(rank .eq. 0 ) then
        allocate(temp_(ld));


        if (colproc .lt. rowproc) then

            temp_ = process_grid(:,1)

            do jj = 1, ld-1
                fromc = mod(jj,ld) +1
                process_grid(:,jj) = process_grid(:,fromc)
            end do
                process_grid(:,ld) = temp_



            do jj = 1, ld
                if (process_grid(jj,colproc) .ne. 0 ) then
                    if( findval(process_grid(jj,1:colproc), -abs(process_grid(jj,colproc)) ) .gt. 0)    then
                    process_grid(jj ,findval(process_grid(jj,1:colproc), -abs(process_grid(jj,colproc)) ) ) = 0
                    end if
                end if
                if (process_grid(jj,ld) .ne. 0 ) then
                    if( findval(process_grid(jj,1:colproc), 0)  .gt. 0)    then
                    process_grid(jj ,findval(process_grid(jj,1:colproc), 0) ) = - abs(process_grid(jj,ld))
                    end if
                    !process_grid(jj,ld) =  0
                end if

            end do

        else if (colproc .gt. rowproc) then

            temp_ = process_grid(1,:)

            do jj = 1, ld-1
                fromc = mod(jj,ld) +1
                process_grid(jj,:) = process_grid(fromc,:)
            end do
                process_grid(ld,:) = temp_


            do jj = 1, ld
                if (process_grid(rowproc,jj) .ne. 0 ) then
                    if( findval(process_grid(1:rowproc, jj), -abs(process_grid(rowproc, jj))) .gt. 0) then
                    process_grid(findval(process_grid(1:rowproc, jj), -abs(process_grid(rowproc, jj))), jj ) = 0
                    end if
                end if
                if (process_grid(ld,jj) .ne. 0 ) then
                       if( findval(process_grid(1:rowproc,jj), 0)  .gt. 0)    then
                    process_grid(findval(process_grid(1:rowproc,jj), 0), jj ) = - abs(process_grid(ld,jj))
                    end if
                end if

            end do
        else
            temp_ = process_grid(1,:)

            do jj = 1, ld-1
                fromc = mod(jj,ld) +1
                process_grid(jj,:) = process_grid(fromc,:)
            end do
                process_grid(ld,:) = temp_
        end if
            deallocate(temp_);

            call MPI_BCAST(process_grid, size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )

        else
            call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )

        end if


    end subroutine


    !INIT CUBLASXT CONTEXT
    subroutine init_cublasXT(ngpu, h)
            use iso_c_binding
            integer:: ngpu
            type(c_ptr) :: h

            call init_xt_context(ngpu, h)

    end subroutine

   subroutine destroy_cublasXT(h)
        use iso_c_binding
        type(c_ptr) :: h

        call destroy_xt_context(h)

   end subroutine

    !Inititalizes the process grid
    !Block Cij is located alwais in the ij process
    !Aik is in the i row
    !Bkj in the j column.
    !It selects the k to assing at each ij grid block

    subroutine init_process_grid(process_grid, rowproc, colproc,rank,idxrow, idxcol)
        implicit none
        integer , dimension(:,:), allocatable :: process_grid
        integer , dimension(:), allocatable :: setofpieces
        integer :: ii, jj, kk, rowproc, colproc, ierr, rank
        integer :: idxrow, idxcol

        idxcol= rank/rowproc+1
        idxrow= mod(rank,rowproc)+1

        !if  ( .not. allocated(process_grid)) then
        !	print*, "Error, Process grid not allocated "
        !end if

        allocate(setofpieces(min(rowproc, colproc)))

        process_grid = 0
        do ii = 1, max(rowproc,colproc)
            do jj = 1, max(rowproc,colproc) !columns
                process_grid(ii,jj) = mod( (ii+jj-2) , max(rowproc,colproc) )+1

                if( process_grid(ii,jj) .gt. min(colproc,rowproc) ) then
                    process_grid(ii,jj) = 0 !Invalid Pieces
                end if
            end do

        end do

        do kk = 1, min(rowproc, colproc)
            setofpieces(kk) = kk
        end do

        do ii = 1, min(rowproc, colproc)
        do jj = 1, max(rowproc, colproc)

        if(rowproc .gt. colproc) then
            if(process_grid(jj,ii) .eq. 0) then
            do kk = 1, min(rowproc, colproc)
                if(.not. any(abs(process_grid(jj,1:colproc)) == setofpieces(kk)) ) then
                    process_grid(jj,ii) = - setofpieces(kk)
                end if
            end do
            end if

        else if(rowproc .lt. colproc) then
            if(process_grid(ii,jj) .eq. 0) then
            do kk = 1, min(rowproc, colproc)
                if(.not. any(abs(process_grid(1:rowproc, jj)) == setofpieces(kk)) ) then
                    process_grid(ii,jj) = - setofpieces(kk)
                end if
            end do
            end if
        end if


        end do
        end do


        deallocate(setofpieces)

        !Dispatch Process Grid to all processes
        if(rank .eq. 0 ) then
        !call printmatrix(process_grid(:rowproc, :colproc))
        call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )
        else
        call MPI_BCAST(process_grid,size(process_grid,dim=1)**2, MPI_INTEGER, 0 ,MPI_COMM_WORLD,ierr )
        end if

    end subroutine



    !given row and column in the process grid computes the rank

    function grid2rank(idxrow, idxcol , rowproc, colproc) result(idx)
        implicit none
        integer :: idxrow, idxcol, idx, rowproc, colproc
    	if (idxcol .le. 0 ) then
    		print*, "grid2rank error, ivalid idx col or row "
    	endif
    	if (idxrow .le. 0 ) then
    		print*, "grid2rank error, ivalid idx col or row "
    	endif
    	if (idxcol .gt. colproc ) then
    		print*, "grid2rank error, ivalid idx col or row "
    	endif
    	if (idxrow .gt. rowproc ) then
    		print*, "grid2rank error, ivalid idx col or row "
    	endif
        idx = (idxcol-1)*rowproc + idxrow-1

        return

    end function


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!---------------------TO MOVE AROUND DATA IN THE MPI HELL ------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!Set/ get to/from the  process "who" the right line (row or column ) of a distributed matrix
!Linedescriptor is equal to a matrix descriptor in the form if line is row.. if it is column..

    subroutine cbp_setgetprcline(Mat, matdescriptor, line, linedescriptor, who, rank, nprocs, what)
       integer :: rank, nprocs, who
       CTYPE , dimension(:,:) :: Mat
       CTYPE, dimension(:) :: line
       character*3 :: what

       integer , dimension(4) :: linedescriptor
       integer , dimension(4) :: matdescriptor, other
       integer , dimension(4) :: row_lims, col_lims
       integer, dimension(2):: slice  , localrow, localcol

       integer :: ii , ierr, Rs, RR

       if ( any(matdescriptor .lt. 0) ) then
       	!print*, matdescriptor
       	print*, "error invalid descriptor"
       end if
       if( size(mat) .lt. product(matdescriptor(:2))) then
       	print*, "Matrix too small for descriptor, cbp_setgetprcline", size(mat), product(matdescriptor(:2))
       end if
         if(who .ge. nprocs) then
      		print*, "Invalid process cbp_setgetprcline"
       end if
       if( size(line) .lt. maxval(linedescriptor(:2)) .and.  maxval(linedescriptor(:2)) .ne. 1 ) then
       	print*, "Error, wrong line sizes", size(line, dim=1), maxval(linedescriptor(:2))
       end if
         if(size(shape(line)) .ne. 1) then
       	print*, "Error, wrong line sizes"
       end if


       do ii = 0, nprocs-1
       	!check me with me
       	if (who .eq. ii) then

       		other = matdescriptor
       		! mat start, mat end, chunk start , chunkend, all global
			row_lims = [ other(3), other(1)+other(3)-1, linedescriptor(3) , linedescriptor(3)+ linedescriptor(1) -1]
			col_lims = [ other(4), other(2)+other(4)-1, linedescriptor(4) , linedescriptor(4)+ linedescriptor(2) -1]


		!check me with the others
		else if (rank .eq. who) then
			!I send the line descriptor
			!Matr
			call MPI_RECV( other, 4 , MPI_INTEGER, ii, 10*ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
			!line
			call MPI_SEND( linedescriptor, 4, MPI_INTEGER, ii, 10*ii+1, MPI_COMM_WORLD, ierr )


			row_lims = [ other(3), other(1)+other(3)-1 ,linedescriptor(3) , linedescriptor(3)+ linedescriptor(1) -1 ]
			col_lims = [ other(4), other(2)+other(4)-1 ,linedescriptor(4) , linedescriptor(4)+ linedescriptor(2) -1 ]



		else if (rank .eq. ii) then
			!He send the full descriptor
			!Matr
			call MPI_SEND( matdescriptor, 4 , MPI_INTEGER, who, 10*ii, MPI_COMM_WORLD,  ierr )

			!line
			call MPI_RECV( other, 4, MPI_INTEGER, who, 10*ii+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

       		row_lims = [ matdescriptor(3), matdescriptor(3) + matdescriptor(1) -1, other(3), other(1)+other(3)-1 ]
			col_lims = [ matdescriptor(4), matdescriptor(4) + matdescriptor(2) -1, other(4), other(2)+other(4)-1 ]



       	end if
       	!Look for overlap

       	if ((rank .eq. ii) .or. (rank .eq. who)) then



		      if (( minval(row_lims(2:4:2)) - maxval(row_lims(1:3:2)) .ge. 0) .and. &
		           (minval(col_lims(2:4:2)) - maxval(col_lims(1:3:2)) .ge. 0 )) then

		      	       if      (row_lims(4)-row_lims(3) .eq. 0) then
		      	       	slice = [maxval(col_lims(1:3:2))-col_lims(3)+1 , minval(col_lims(2:4:2))-col_lims(3)+1]
		      	       else if (col_lims(4)-col_lims(3) .eq. 0) then
		      	       	slice = [maxval(row_lims(1:3:2))-row_lims(3)+1 , minval(row_lims(2:4:2))-row_lims(3)+1]
		      	       end if

       		       if (what .eq. 'set') then
					if (rank.eq.ii .and. rank.eq.who ) then
						localrow = [maxval(row_lims(1:3:2))-row_lims(1)+1 , minval(row_lims(2:4:2))-row_lims(1)+1]
					 	localcol = [maxval(col_lims(1:3:2))-col_lims(1)+1 , minval(col_lims(2:4:2))-col_lims(1)+1 ]
						!print*, localrow, '|', localcol
							if (localrow(2)-localrow(1) .eq. 0 ) then
								mat( localrow(1) , localcol(1) : localcol(2)) =  line(slice(1) : slice(2))
							else if (localcol(2)-localcol(1) .eq. 0 ) then
								mat( localrow(1):localrow(2) , localcol(1)) = line(slice(1) : slice(2))
							else
							 	mat( localrow(1) , localcol(1)) = line(slice(1))
							end if
					else if (rank .eq. who) then

						call MPI_SEND( line( slice(1) : slice(2)) , &
						       slice(2)-slice(1)+1, &
						       MPI_CTYPE, ii , 10*ii+3, MPI_COMM_WORLD,ierr )


					else if (rank .eq. ii) then
						localrow = [maxval(row_lims(1:3:2))-row_lims(1)+1 , minval(row_lims(2:4:2))-row_lims(1)+1]
						localcol = [maxval(col_lims(1:3:2))-col_lims(1)+1 , minval(col_lims(2:4:2))-col_lims(1)+1 ]

						call MPI_RECV(  mat( localrow(1) : localrow(2) , localcol(1) : localcol(2)), &
							         slice(2)-slice(1)+1, &
							         MPI_CTYPE, who, 10*ii+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

					end if

					!if (rank.eq.who.and.rank.ne.ii) then
					!	call  MPI_Wait(RS, MPI_STATUS_IGNORE, IERR)
					!
					!end if
					!if (rank.eq.ii.and.rank.ne.who) then
					!	call  MPI_Wait(RR, MPI_STATUS_IGNORE, IERR)

					!end if

	       		       if (rank.eq.who) then
	       		   		!print*,ii, who ,col_lims,'|' ,localcol, localcol(2) -localcol(1), '|',&
	       		   		!	 sum (mat( localrow(1) : localrow(2) , localcol(1) : localcol(2)))
	       		   	end if
	       		   	if (rank.eq.ii) then
	       		   		!print*,ii, who ,col_lims,'|',  slice , slice(2) -slice(1) , '|', sum (line(slice(1):slice(2)))
	       		   	end if


	       		else if (what .eq. 'get') then
					if (rank.eq.ii .and. rank.eq.who ) then
						localrow = [maxval(row_lims(1:3:2))-row_lims(1)+1 , minval(row_lims(2:4:2))-row_lims(1)+1]
						localcol = [maxval(col_lims(1:3:2))-col_lims(1)+1 , minval(col_lims(2:4:2))-col_lims(1)+1 ]

						if (localrow(2)-localrow(1) .eq. 0 ) then

							line(slice(1) : slice(2)) =  mat( localrow(1) , localcol(1) : localcol(2))


						else if (localcol(2)-localcol(1) .eq. 0 ) then

							line(slice(1) : slice(2)) =  mat( localrow(1):localrow(2) , localcol(1))
						else
							line(slice(1)) = mat( localrow(1) , localcol(1))
						end if

					else if (rank .eq. ii ) then

						localrow = [maxval(row_lims(1:3:2))-row_lims(1)+1 , minval(row_lims(2:4:2))-row_lims(1)+1]
						localcol = [maxval(col_lims(1:3:2))-col_lims(1)+1 , minval(col_lims(2:4:2))-col_lims(1)+1 ]

						call MPI_SEND( mat( localrow(1) : localrow(2) , localcol(1) : localcol(2)), &
									slice(2)-slice(1)+1, &
									MPI_CTYPE, who, 10*ii+5, MPI_COMM_WORLD ,ierr )

					else if (rank .eq. who) then


						call MPI_RECV( line( slice(1) : slice(2)), &
							slice(2)-slice(1)+1, &
							MPI_CTYPE, ii, 10*ii+5, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

					end if

					!if (rank.eq.ii.and.rank.ne.who) then
					!	call  MPI_Wait(RS, MPI_STATUS_IGNORE, IERR)
					!
					!end if
					!if (rank.eq.who.and.rank.ne.ii) then
					!	call  MPI_Wait(RR, MPI_STATUS_IGNORE, IERR)
					!end if

					if (rank.eq.ii) then
	       		   		!print*,ii, who ,col_lims,'|' ,localcol, localcol(2) -localcol(1), '|',&
	       		   		!	 sum (mat( localrow(1) : localrow(2) , localcol(1) : localcol(2)))
	       		   	end if
	       		   	if (rank.eq.who) then
	       		   		!print*,ii, who ,col_lims,'|',  slice , slice(2) -slice(1) , '|', sum (line(slice(1):slice(2)))
	       		   	end if

	       		else
	       			print*, "Error, mod not available"
       		        end if

			 end if !overlap

       	end if !processes

       end do

    end subroutine
!====================================================================================
!=====================================================================================

!	subroutine scalapack_to_cbp(A,descA, descriptorA, rank, nprocs, modo)
!	CTYPE, dimension(:,:) :: A
!	CTYPE, dimension(:,:), allocatable :: Atemp
!	integer , dimension(4) :: descriptorA, desc_chunk !cbp
!	integer, dimension(9) :: descA !scalapack
!	integer, dimension(2) :: np_sca, coord_sca , start, noslices, slicecoord, local_sca
!
!	integer :: ii, ll, kk , jj
!	character*3:: modo
!	!get block info
!	call BLACS_GRIDINFO(descA(2), np_sca(1) , np_sca(2) , coord_sca(1) , coord_sca(2) )
!	allocate(Atemp())
!	noslices = [dimA(1) / npsca(1)*descA(6), dimA(2) / npsca(2)*descA(7) ]
!
!	if ( mod(dimA(1),npsca(1)*descA(6) ) .gt. 0 ) noslices(1) = noslices(1) +1
!	if ( mod(dimA(2),npsca(2)*descA(7) ) .gt. 0 ) noslices(2) = noslices(2) +1
!
!	do kk = 1, noslices(1)
!	do ll = 1, noslices(2)
!		slicecoord = [kk*npsca(1)*descA(6)+1, ll*npsca(2)*descA(7)+1 ] !upper left of slice
!		do ii = 0, np_sca(1)-1
!			do jj = 0, np_sca(2)-1 !for every process of the scalapack  grid
!
!			start = slicecoord + [ii*descA(6), jj*descA(7) ] !coords inside slice no need to add 1
!
!			if ( all ( start .le. dimA) ) then
!				!global coordinates of the block
!				desc_chunk = [ min( descA(6), dimA(1) - start(1)+1) , min( descA(7), dimA(2) - start(2)+1), &
!					       start(1) , start(2)  ]
!				!where it is locally in A ?
!				local_sca = [kk*descA(6)+1, ll*descA(7)+1]
!
!				cbp_setgetprcblock(Atemp, descriptorA, &
!								   A(local_sca(1):local_sca(1)+desc_chunk(1), local_sca(2):local_sca(2)+desc_chunk(2) ),&
!						           desc_chunk, ii*np_sca(2)+ll, rank, nprocs, modo)
!
!			end if
!			end do
!		end do
!	end do
!	end do
!	end do

!	end subroutine
!Set the chunk from process who to the right distributed matrix
!or get at the process who

!Never tested surely it does not work
!But it is similar to the one above

    subroutine cbp_setgetprcblock(Mat, matdescriptor, chunk, chunkdescriptor, who, rank, nprocs, what)
       integer :: rank, nprocs, who
       CTYPE , dimension(:,:) :: Mat, chunk

       integer , dimension(4) :: chunkdescriptor
       integer , dimension(4) :: matdescriptor, other
       integer , dimension(4) :: row_lims, col_lims
       integer, dimension(2) :: overlap , localcoord
       character*3 :: what
       integer :: ii , ierr

       if(who .ge. nprocs) then
      		print*, "Invalid process cbp_setgetprcblock"
       end if


       do ii = 0, nprocs-1
       	!check me with me
       		if (who .eq. ii) then

       			other = matdescriptor
       			! mat start, mat end, chunk start , chunkend
			row_lims = [ other(3), other(1)+other(3)-1, chunkdescriptor(3) , chunkdescriptor(3)+ chunkdescriptor(1) -1]
			col_lims = [ other(4), other(2)+other(4)-1, chunkdescriptor(4) , chunkdescriptor(4)+ chunkdescriptor(2) -1]


		!check me with the others
		else if (rank .eq. who) then
			!I send the chunk descriptor
			!Matr
			call MPI_RECV( other, 4 , MPI_INTEGER, ii, 10*ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
			!Chunk
			call MPI_SEND( chunkdescriptor, 4, MPI_INTEGER, ii, 10*ii+1, MPI_COMM_WORLD, ierr )


			row_lims = [ other(3), other(1)+other(3)-1 ,chunkdescriptor(3) , chunkdescriptor(3)+ chunkdescriptor(1) -1 ]
			col_lims = [ other(4), other(2)+other(4)-1 ,chunkdescriptor(4) , chunkdescriptor(4)+ chunkdescriptor(2) -1 ]



		else if (rank .eq. ii) then
			!He send the full descriptor
			!Matr
			call MPI_SEND( matdescriptor, 4 , MPI_INTEGER, who, 10*ii, MPI_COMM_WORLD, ierr )
			!Chunk
			call MPI_RECV( other, 4, MPI_INTEGER, who, 10*ii+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

       		        row_lims = [ matdescriptor(3), matdescriptor(3) + matdescriptor(1) -1, other(3), other(1)+other(3)-1 ]
			col_lims = [ matdescriptor(4), matdescriptor(4) + matdescriptor(2) -1, other(4), other(2)+other(4)-1 ]



       		end if
       	!Look for overlap
       	if ((rank .eq. ii) .or. (rank .eq. who)) then

       		 overlap = [ minval(row_lims(2:4:2)) - maxval(row_lims(1:3:2)) , minval(col_lims(2:4:2)) - maxval(col_lims(1:3:2)) ]

       		 if (all(overlap.gt. 0)) then

       		      if (what .eq. 'set') then
       		      	       !------------------------------------------
       		      	      if (rank .eq. who .and. who .eq. ii ) then

       		      	      localcoord = [ maxval(row_lims(1:3:2))-row_lims(3)+1 , maxval(col_lims(1:3:2))-col_lims(3)+1 ]

       		      	      	mat( localcoord(1) : localcoord(1)+overlap(1), &
	       		      	     localcoord(2) : localcoord(1)+overlap(2) ) = &
							       		      	      chunk( localcoord(1) : localcoord(1)+overlap(1), &
							       		      	             localcoord(2) : localcoord(1)+overlap(2))

	       		      else if (rank .eq. who) then

	       		        localcoord = [ maxval(row_lims(1:3:2))-row_lims(3)+1 , maxval(col_lims(1:3:2))-col_lims(3)+1 ]

	       		      	call MPI_SEND( chunk( localcoord(1) : localcoord(1)+overlap(1), &
	       		      			      localcoord(2) : localcoord(1)+overlap(2)), &
	       		      	 		      (overlap(1)+1)*(overlap(2)+1), &
	       		      			      MPI_CTYPE, ii , 10*ii+2, MPI_COMM_WORLD, ierr )

	       		      else if (rank.eq.ii) then

				 localcoord = [ maxval(row_lims(1:3:2))-row_lims(1)+1 , maxval(col_lims(1:3:2))-col_lims(1)+1 ]


	       		      	 call MPI_RECV( mat( localcoord(1) : localcoord(1)+overlap(1) , &
	       		      			             localcoord(2) : localcoord(1)+overlap(2)), &
	       		      			        (overlap(1)+1)*(overlap(2)+1), &
	       		      	 		        MPI_CTYPE, who, 10*ii+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

	       		      end if
	       		        !------------------------------------------
	       		else if (what .eq. 'get') then
	       		    !------------------------------------------
					if (rank .eq. who .and. who .eq. ii ) then

						localcoord = [ maxval(row_lims(1:3:2))-row_lims(3)+1 , maxval(col_lims(1:3:2))-col_lims(3)+1 ]


						chunk( localcoord(1) : localcoord(1)+overlap(1), &
							   localcoord(2) : localcoord(1)+overlap(2)) = &
						   		      	      mat( localcoord(1) : localcoord(1)+overlap(1), &
												localcoord(2) : localcoord(1)+overlap(2) )

					else if (rank .eq. who) then

						localcoord = [ maxval(row_lims(1:3:2))-row_lims(3)+1 , maxval(col_lims(1:3:2))-col_lims(3)+1 ]

						call MPI_RECV( chunk( localcoord(1) : localcoord(1)+overlap(1), &
								  localcoord(2) : localcoord(1)+overlap(2)), &
								  (overlap(1)+1)*(overlap(2)+1), &
								  MPI_CTYPE, ii , 10*ii+2, MPI_COMM_WORLD,MPI_STATUS_IGNORE,  ierr )

					else if (rank.eq.ii) then

						localcoord = [ maxval(row_lims(1:3:2))-row_lims(1)+1 , maxval(col_lims(1:3:2))-col_lims(1)+1 ]


						call MPI_SEND( mat( localcoord(1) : localcoord(1)+overlap(1), &
										 localcoord(2) : localcoord(1)+overlap(2) ),  &
									(overlap(1)+1)*(overlap(2)+1), &
									MPI_CTYPE, who, 10*ii+2, MPI_COMM_WORLD, ierr )

					end if
	       		     !--------------------------------------------
	       		  end if ! set get
       	  end if !overlap

		end if !ranks



       end do



    end subroutine
!===========================================================================================
!===========================================================================================
function cbp_getval2d(Mat, matdescriptor, valpos1, rank, nprocs) result(val)
       integer :: rank, nprocs
       CTYPE , dimension(:,:) :: Mat
       CTYPE :: val
       integer , dimension(4) :: matdescriptor
       integer , dimension(2) :: valpos1
       integer :: ii

       do ii = 0 , nprocs-1
       	 call cbp_setgetprcval(Mat, matdescriptor, val, valpos1, ii, rank, nprocs, 'get')
       end do


end  function


!===========================================================================================
!===========================================================================================

!Set the chunk at process sender to the right distributed matrix
    subroutine cbp_setgetprcval(Mat, matdescriptor, val, valpos1, who, rank, nprocs, what)

       integer :: rank, nprocs, who
       CTYPE , dimension(:,:) :: Mat
       CTYPE :: val


       integer , dimension(4) :: matdescriptor, other
       integer , dimension(2) :: row_lims, col_lims, valpos, valpos1

       integer :: ii ,ierr
       character*3 :: what

       if(who .ge. nprocs) then
      		print*, "Invalid process cbp_setgetprcval"
       end if
       !print*, "hello"

       if(rank.ne.who) then
       	valpos =  0
       else
       	valpos = valpos1
       end if

       do ii = 0, nprocs-1
       	!check me with me
       	!print*, ii , who, rank
       	if (who .eq. ii) then

       		other = matdescriptor

			row_lims = [ other(3), other(1)+other(3)-1 ]
			col_lims = [ other(4), other(2)+other(4)-1 ]

		!check me with the others
		else if (rank .eq. who) then
			!I send the chunk descriptor
			!Matr
			!print*, "in who"
			call MPI_RECV( other, 4 , MPI_INTEGER, ii, 10*ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
			!Chunk
			!print*, "whorecv"
			call MPI_SEND( valpos, 2, MPI_INTEGER, ii, 10*ii+1, MPI_COMM_WORLD, ierr )


			row_lims = [ other(3), other(1)+other(3)-1 ]
			col_lims = [ other(4), other(2)+other(4)-1 ]


		else if (rank .eq. ii) then
			!He send the full descriptor
			!Matr
			!print*, "in ii"
			call MPI_SEND( matdescriptor, 4 , MPI_INTEGER, who, 10*ii, MPI_COMM_WORLD, ierr )
			!Chunk

			call MPI_RECV( valpos, 2, MPI_INTEGER, who, 10*ii+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
			!print*, "iirecv"
       		row_lims = [ matdescriptor(3), matdescriptor(3) + matdescriptor(1) -1 ]
			col_lims = [ matdescriptor(4), matdescriptor(4) + matdescriptor(2) -1 ]



       	end if
       	!Look for overlap
       	!print*,ii, who, "-->", rank ,'|', row_lims, '|', col_lims , '|', &
       	!	( ( (valpos(1)-row_lims(1))*(valpos(1)-row_lims(2)) .le. 0) .and. &
       	!	     ( (valpos(1)-col_lims(1))*(valpos(1)-col_lims(2)) .le. 0 ) ), &
       	!	     valpos
       	!print*,ii, who, "-->", rank'|', col_lims, '|', rank , ii, who, valpos,

       	if ((rank .eq. ii) .or. (rank .eq. who)) then

       		 if( ( (valpos(1) .ge. row_lims(1)) .and. (valpos(1) .le. row_lims(2)) ) .and. &
       		     ( (valpos(2) .ge. col_lims(1)) .and. (valpos(2) .le. col_lims(2))  ))then

       		      if (what.eq. 'set') then

	       		     if (rank .eq. who) then

	       		      	call MPI_SEND( val , 1 , MPI_CTYPE, ii , 10*ii+8, MPI_COMM_WORLD, ierr )
	       		      	!print*, "val", val
	       		     end if
	       		     if (rank .eq. ii) then



	       		     	call MPI_RECV( Mat( valpos(1) - row_lims(1)  +1 ,&
	       		     			     valpos(2) - col_lims(1)  +1 ) &
	       		     			   ,1 ,MPI_CTYPE,who,10*ii+8, MPI_COMM_WORLD,MPI_STATUS_IGNORE,  ierr )
	       		       ! print*, "old" , Mat( valpos(1) - row_lims(1)  +1 ,&
	       		     	!		    valpos(2) - col_lims(1)  +1 ) , &
	       		     	!		     valpos(1) - row_lims(1)  +1, &
	       		     	!		     valpos(2) - col_lims(1)  +1, row_lims(2), col_lims(2)
	       		     	!print*, "done", ii, who
	       		      end if

       		      else if (what .eq. 'get') then
       		         if (rank .eq. ii) then

	       		     	call MPI_SEND( Mat( valpos(1) - row_lims(1)  +1 ,&
	       		     			   valpos(2) - col_lims(1)  +1 ), &
	       		     			   1 ,MPI_CTYPE,who,10*ii+2,MPI_COMM_WORLD, ierr )
	       		      end if
       		        if (rank .eq. who) then

	       		      	call MPI_RECV( val , 1 , MPI_CTYPE, ii , 10*ii+2, MPI_COMM_WORLD,&
	       		      		MPI_STATUS_IGNORE, ierr )
	       		end if


       		      else
       		      	print*, "Error, mod not available"
       		      endif

			 end if


       	end if

       end do



    end subroutine

!===================================================================================
!===================================================================================


!Given a Distributed Matrix (i.e. many little matrices and their descriptors) move
!move data in order to make the matrix have the new descriptors

!At the moment it is unsafe since it does not check that localmat ( understood in
!terms of space to keep numbers ) are large enough to keep the data prescribed by the
!new descriptors

!Additionally it is unsafe because it does not implement batched MPI SEND and RECEIVE
!necessary in the case of matrices so big that the block size trepass huge(int*4)

 subroutine cbp_redispatch_2d(localMat, olddescriptor, newdescriptor, rank, nprocs)
 	 integer , dimension(4) :: olddescriptor , newdescriptor, ini, fin
         CTYPE , dimension(:,:) :: localMat
         CTYPE , dimension(:,:) , allocatable :: localtemp
         integer, dimension(2) :: localrow, localcol, starts , enbs
         integer :: ii, jj, nprocs, rank , ierr

         allocate(localtemp(size(localmat, dim=1), size(localmat, dim=2)))
         localtemp = CTYPE0

                do ii = 0, nprocs-1
                do jj = 0 , nprocs-1
        	if(rank.eq.ii .or. rank.eq.jj) then
        		!|1) EXCHANGE
        		if (rank .eq. ii .and. rank .eq. jj) then
        			ini = olddescriptor
        			fin = newdescriptor
        		else if ( rank .eq. ii) then
				ini = olddescriptor
				call MPI_SEND( olddescriptor, 4, MPI_INTEGER,jj,jj,MPI_COMM_WORLD, ierr )
				call MPI_RECV( fin, 4, MPI_INTEGER, jj, ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        		else if(rank .eq. jj) then
        			ini = olddescriptor
				call MPI_RECV( ini, 4, MPI_INTEGER, ii, jj, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				call MPI_SEND( fin, 4, MPI_INTEGER, ii,ii,MPI_COMM_WORLD, ierr )
        		end if
        		!|2) OVERLAP
        		!overlap rows .and. overlap cols
        		enbs = [ max(ini(3),fin(3)), max(ini(4),fin(4)) ]
        		starts = [ min( sum(ini(1:3:2)), sum(fin(1:3:2)) ) -1 ,  &
        			   min( sum(ini(2:4:2)), sum(fin(2:4:2)) ) -1 ]


        		if( (( enbs(1) - starts(1) )  .ge. 0 ) .and. (( enbs(2) - starts(2) ) .ge. 0)  ) then
        		 !|3) PASS DATA
        		 	if (rank.eq.ii.and.rank.eq.jj) then
        		 		localrow = [ starts(1)-olddescriptor(3)+1, enbs(1) -  olddescriptor(3) +1 ]
        		 		localcol = [ starts(2)-olddescriptor(4)+1, enbs(2) -  olddescriptor(4) +1 ]

        		 		localTEMP(localrow(1):localrow(2), localcol(1):localcol(2)) = &
        		 			localmat(localrow(1):localrow(2), localcol(1):localcol(2))


				else if (rank.eq.ii) then
					localrow = [ starts(1)-olddescriptor(3)+1, enbs(1) -  olddescriptor(3) +1 ]
        		 		localcol = [ starts(2)-olddescriptor(4)+1, enbs(2) -  olddescriptor(4) +1 ]

        		 	        call MPI_SEND( localmat(localrow(1):localrow(2), localcol(1):localcol(2)),&
        		 	       		 (localrow(2)-localrow(1)+1)*(localcol(2)-localcol(1)+1), &
        		 	       		  MPI_CTYPE, jj,ii,MPI_COMM_WORLD, ierr )


				else if(rank.eq.jj) then
					localrow = [ starts(1)-newdescriptor(3)+1, enbs(1) - newdescriptor(3) +1 ]
        		 		localcol = [ starts(2)-newdescriptor(4)+1, enbs(2) - newdescriptor(4) +1 ]

					  call MPI_RECV( localTEMP(localrow(1):localrow(2), localcol(1):localcol(2)),&
        		 	       		  	 (localrow(2)-localrow(1)+1)*(localcol(2)-localcol(1)+1), &
        		 	       		  	  MPI_CTYPE, ii,ii,MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

				end if
        		end if


        	end if
                end do
                end do
                localtemp = localmat
                deallocate(localtemp)

 end subroutine

 subroutine cbp_redispatch(localMat, olddescriptor, newdescriptor, rank, nprocs)
        integer , dimension(4) :: olddescriptor , newdescriptor, other
        CTYPE , dimension(:,:) :: localMat
        CTYPE , dimension(:,:), allocatable :: localtemp
        integer :: rank, nprocs, ii, jj, kk, ierr
        integer , dimension(4) :: orderr, orderc
        integer , dimension(4) :: row_lims, col_lims
        integer, dimension(MPI_STATUS_SIZE) :: stat
        logical :: colok, rowok
        integer :: tmp
        allocate(localtemp(size(localmat,dim = 1),size(localmat,dim = 2) ))
        !localMat(:,:,2) = localMat(:,:,1)
        !localMat(:,:,1) = 0.0

	localtemp = localmat
	localmat = CTYPE0

	if(newdescriptor(1) .gt. size(localMat, dim=1) ) then
		print*, "Error, matrix too small redispatch"
	end if

	if(newdescriptor(2) .gt. size(localMat, dim=2) ) then
			print*, "Error, matrix too small redispatch"
	end if

        do ii = 0, nprocs-1
                do jj = 0 , nprocs-1

                        if (rank .eq. ii) then

                                if(ii .ne. jj) then

                                call MPI_SEND( newdescriptor, 4, MPI_INTEGER,jj,ii,MPI_COMM_WORLD, ierr )
                                call MPI_RECV( other, 4, MPI_INTEGER, jj, jj, MPI_COMM_WORLD, stat, ierr)

                                else
                                        other = olddescriptor
                                end if
                                !compare overlap
                                row_lims = [ newdescriptor(3), newdescriptor(3)+newdescriptor(1)-1 , other(3), other(1)+other(3)-1 ]
                                col_lims = [ newdescriptor(4), newdescriptor(4)+newdescriptor(2)-1 , other(4), other(2)+other(4)-1 ]


                                !print*, ii, jj , newdescriptor, other, "r"


                                !print*, row_lims, col_lims


                                 if(( min(newdescriptor(3)+newdescriptor(1),other(3)+other(1))- &
                                      max(newdescriptor(3), other(3))  .gt. 0 ).and. &
                                    ( min(newdescriptor(4)+newdescriptor(2),other(4)+other(2))- &
                                      max(newdescriptor(4), other(4))  .gt. 0 ) ) then


                                   !    print*,ii, jj , "||"   ,max(newdescriptor(3),other(3)) , &
                                    !             min(newdescriptor(3)+newdescriptor(1),other(3)+other(1))-1,&
                                     !           max(newdescriptor(4),other(4)), &
                                      !           min(newdescriptor(4)+newdescriptor(2),other(4))+other(2)-1,"r"

                                        if (rank .ne. jj ) then
                                        call  MPI_RECV( &
                                        localMat(max(newdescriptor(3),other(3))-newdescriptor(3)+1 :&
                                                 min(newdescriptor(3)+newdescriptor(1),other(3)+other(1))-newdescriptor(3), &
                                                 max(newdescriptor(4),other(4))-newdescriptor(4)+1:&
                                                 min(newdescriptor(4)+newdescriptor(2),other(4)+other(2))-newdescriptor(4)) , &
                                                (min(newdescriptor(3)+newdescriptor(1), other(3)+other(1)) - &
                                                 max(newdescriptor(3),other(3))) * &
                                                (min(newdescriptor(4)+newdescriptor(2), other(4)+other(2)) - &
                                                max(newdescriptor(4),other(4))) ,&
                                                 MPI_CTYPE, jj, ii, MPI_COMM_WORLD, stat, ierr)
                                        else
                                        localMat(max(newdescriptor(3),olddescriptor(3))-newdescriptor(3)+1 :&
                                                 min(newdescriptor(3)+newdescriptor(1),&
                                                 olddescriptor(3)+olddescriptor(1))-newdescriptor(3), &
                                                 max(newdescriptor(4),olddescriptor(4))-newdescriptor(4)+1:&
                                                 min(newdescriptor(4)+newdescriptor(2),&
                                                 olddescriptor(4)+olddescriptor(2))-newdescriptor(4)) = &
                                         localtemp(max(olddescriptor(3),newdescriptor(3))-olddescriptor(3)+1:&
                                                 min(olddescriptor(3)+olddescriptor(1),&
                                                 newdescriptor(3)+newdescriptor(1))-olddescriptor(3), &
                                                 max(olddescriptor(4),newdescriptor(4))-olddescriptor(4)+1:&
                                                 min(olddescriptor(4)+olddescriptor(2),&
                                                 newdescriptor(4)+newdescriptor(2))-olddescriptor(4))

                                        end if

                                end if
                            end if





                        if (rank .eq. jj .and. jj .ne. ii ) then
                                !if different send otherwise do nothing, all handled above

                                !print*, ii, jj, "firstmet s"
                                call MPI_RECV( other, 4, MPI_INTEGER, ii, ii, MPI_COMM_WORLD,stat, ierr)
                                call MPI_SEND( olddescriptor, 4, MPI_INTEGER , ii , jj, MPI_COMM_WORLD, ierr )



                                row_lims = [ olddescriptor(3), olddescriptor(3)+newdescriptor(1)-1 , other(3), other(1)+other(3)-1 ]
                                col_lims = [ olddescriptor(4), olddescriptor(4)+newdescriptor(2)-1 , other(4), other(2)+other(4)-1 ]

                                !print *, ii, jj, other, olddescriptor, "s"
                                if( ( min(olddescriptor(3)+olddescriptor(1),other(3)+other(1))- &
                                      max(olddescriptor(3), other(3))  .gt. 0 ).and. &
                                   (  min(olddescriptor(4)+olddescriptor(2),other(4)+other(2))- &
                                      max(olddescriptor(4), other(4))  .gt. 0 )) then

                                       ! print*, ii, jj , "||",  max(olddescriptor(3),other(3)),&
                                        !         min(olddescriptor(3)+olddescriptor(1),other(3)+other(1))-1, &
                                         !        max(olddescriptor(4),other(4)),&
                                          !       min(olddescriptor(4)+olddescriptor(2),other(4)+other(2))-1,"s"
                                       ! print*, ii, jj , "overlap  send"
                                        call MPI_SEND( &
                                        localtemp(max(olddescriptor(3),other(3))-olddescriptor(3)+1:&
                                                 min(olddescriptor(3)+olddescriptor(1),other(3)+other(1))-olddescriptor(3), &
                                                 max(olddescriptor(4),other(4))-olddescriptor(4)+1:&
                                                 min(olddescriptor(4)+olddescriptor(2),other(4)+other(2))-olddescriptor(4)) , &
                                                (min(olddescriptor(3)+olddescriptor(1),other(3)+other(1)) - &
                                                 max(olddescriptor(3),other(3))) * &
                                                (min(olddescriptor(4)+olddescriptor(2), other(4)+other(2)) - &
                                                 max(olddescriptor(4),other(4))) ,&
                                                 MPI_CTYPE, ii, ii, MPI_COMM_WORLD, ierr)


                                end if


                        end if


                end do
        end do

    deallocate(localtemp)

    end subroutine


 !Not used at the moment since it is an old version, it is an MPI wrapper to pass matrices
!Blocks between two processes in the process grid , used in the MATMUL.
!Now to allow really big matrices I use the other version

!Local Matrix -> 2D array of geenric dimensions
!MatrixDescriptor ->
!Merow-mecol-> Process grid coordinates of the process who perform the action
!Otherrow-othercol->Process grid coordinates of eho receives (if i am sender), of send (if i am receiver)
!Mod -> 's' send 'r' receive
!Tag -> A number from 0 to 10 Don't worry its made unique for this step.


    subroutine MPI_comm_wrap_matmul2(LocalMatrix, MatrixDescriptor ,&
                      merow,mecol,otherrow,othercol,rowproc,colproc,&
                      mod,tag ,requestM, requestD)

        implicit none
        CTYPE, dimension(:,:) :: LocalMatrix
        integer, dimension(4) :: MatrixDescriptor

        integer :: dims(2)
        integer :: otherrow,othercol, ierr, tag, merow, mecol, rowproc, colproc
        integer :: requestM , requestD
        character*1 :: mod

        dims(1) = size(LocalMatrix, dim = 1)
        dims(2) = size(LocalMatrix, dim = 2)
        !print *,"dimension : ",  dims


        if (mod .eq. 's') then

            call MPI_ISEND( LocalMatrix, dims(1)*dims(2), MPI_CTYPE,&
                     grid2rank(otherrow,othercol,rowproc, colproc) , &
                     10*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag ,&  !tag is the  dest
                     MPI_COMM_WORLD, requestM, ierr)

            call MPI_ISEND( MatrixDescriptor, 4, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc, colproc) , &
                     10*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag+1 ,&  !tag is the  dest
                     MPI_COMM_WORLD, requestD, ierr)

            call mpic_checkerror(ierr)
        else if (mod .eq. 'r') then

            call MPI_IRECV( LocalMatrix, dims(1)*dims(2), MPI_CTYPE,&
                     grid2rank(otherrow,othercol,rowproc,colproc) , &
                     10*grid2rank(merow,mecol,rowproc,colproc)+2*tag ,&  !tag is the receiver, so me
                     MPI_COMM_WORLD, requestM, ierr)

            call MPI_IRECV( MatrixDescriptor, 4, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc,colproc) , &
                     10*grid2rank(merow,mecol,rowproc,colproc)+2*tag+1 ,&  !tag is the receiver, so me
                     MPI_COMM_WORLD, requestD, ierr)


            call mpic_checkerror(ierr)
        else
            print*, "Mod error"
        end if



    end subroutine

    !Send and Receive a data.
    !Old version, never used in practice and hence not guaranteed

    function MPI_comm_wrap(thedata,merow,mecol,otherrow,othercol,rowproc,colproc,mod,tag) result (request)
        implicit none
        integer, dimension(:) :: thedata
        integer :: dims , otherrow,othercol, request, ierr, tag, merow, mecol, rowproc, colproc
        character*1 :: mod

        dims = size(thedata, dim = 1)
        !print *,"dimension : ",  dims


        if (mod .eq. 's') then
            call MPI_ISEND( thedata, dims, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc, colproc) , &
                     4*grid2rank(otherrow,othercol,rowproc,colproc)+tag ,&  !tag is the  dest
                     MPI_COMM_WORLD, request, ierr)
            call mpic_checkerror(ierr)
        else if (mod .eq. 'r') then
            call MPI_IRECV( thedata, dims, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc,colproc) , &
                     4*grid2rank(merow,mecol,rowproc,colproc)+tag ,&  !tag is the receiver, so me
                     MPI_COMM_WORLD, request, ierr)
            call mpic_checkerror(ierr)
        else
            print*, "Mod error"
        end if


        return
    end function

!This function is friend to the one below, it is needed to
!compute outv the dimension of the MPI request array.
!In practice you call it
!then allocate the request array with dimension = outv
!then pass the request array and outv to the function below
!In general if the matrices are small outv is 1.

function mpi_norequest(LocalMatrix) result(outv)
integer*8 :: dimstmp(2)
integer :: cut , outv, ii , maxint
CTYPE, dimension(:,:) :: LocalMatrix

dimstmp(1) = size(LocalMatrix, dim = 1)
dimstmp(2) = size(LocalMatrix, dim = 2)
ii = 1
outv = 1

	 if(dimstmp(1) .gt. huge(maxint)) then
                        print*, "Matrix too big "
                end if
                cut = dimstmp(2)


                do while (cut*dimstmp(1) .gt. huge(maxint))
                        cut = ceiling(real(dimstmp(2))/real(ii))
                        outv = ii
                        ii = ii +1
                        !print*, ii , cut, out
                end do
end function

!Batched version of MPI_comm_wrap_matmul2

!merow, mecol -> coordinates of the process in the process_grid of
!otherrow, othercol -> coordinates of the process in the process_grid to
!			/from i will receive/send data.
!Tag is the tag of MPI, but look at the code before using...
!rowproc colproc are the sides of the process grid
!mod if you want to send or receive with merow/mecol

!Well, to understand it look at cbp_matmul where i use it..


subroutine MPI_comm_wrap_matmul(LocalMatrix, MatrixDescriptor ,&
                      merow,mecol,otherrow,othercol,rowproc,colproc,&
                      mod,tag ,requestM, requestD, outv)

        implicit none
        CTYPE, dimension(:,:) :: LocalMatrix
        integer, dimension(4) :: MatrixDescriptor

        integer*8 :: dimstmp(2)
        integer :: ii, cut , outv, tmp(2)
        integer :: dims(2)
        integer :: otherrow,othercol, ierr, tag, merow, mecol, rowproc, colproc
        integer, dimension(:):: requestM
        integer ::requestD
        character*1 :: mod
        integer :: maxint

        dimstmp(1) = size(LocalMatrix, dim = 1)
        dimstmp(2) = size(LocalMatrix, dim = 2)
        !print *,"dimension : ",  dims
        dims = dimstmp

        cut = ceiling(real(dimstmp(2))/real(outv))

               ! if(dimstmp(1) .gt. huge(maxint)) then
               !         print*, "Error at 415:cbp_distribute.F90 : Matrix too big "
               ! end if
               ! cut = dimstmp(2)


                !do while (cut*dimstmp(1) .gt. huge(maxint))
                 !       cut = ceiling(real(dimstmp(2))/real(ii))
                  !      outv = ii
                   !     ii = ii +1
                        !print*, ii , cut, out
                !end do
                !print*, "end" , cut*out , dimstmp(2), out, cut


        if (mod .eq. 's') then

            call MPI_ISEND( MatrixDescriptor, 4, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc, colproc) , &
                     10*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag+1 ,&  !tag is the  dest
                     MPI_COMM_WORLD, requestD, ierr)

            call mpic_checkerror(ierr)

            else if (mod .eq. 'r') then

            call MPI_IRECV( MatrixDescriptor, 4, MPI_INTEGER,&
                     grid2rank(otherrow,othercol,rowproc,colproc) , &
                     10*grid2rank(merow,mecol,rowproc,colproc)+2*tag+1 ,&  !tag is the receiver, so me
                     MPI_COMM_WORLD, requestD, ierr)


            call mpic_checkerror(ierr)
        else
            print*, "Mod error"
        end if


        tmp(1)=1

        if ( size(requestM) .lt. (outv)) then
         print*, "Error, RequestM to small", size(requestM), outv
        end if

        do  ii = 1, outv
            tmp(2) = min( dimstmp(2), cut*ii)

            if (mod .eq. 's') then

            call MPI_ISEND( LocalMatrix(:, tmp(1):tmp(2)), dims(1)*(tmp(2)-tmp(1)+1), MPI_CTYPE,&
                     grid2rank(otherrow,othercol,rowproc, colproc) , &
                     100*grid2rank(otherrow,othercol,rowproc,colproc)+2*tag+ii+11 ,&  !tag is the  dest
                     MPI_COMM_WORLD, requestM(ii), ierr)

            call mpic_checkerror(ierr)

            else if (mod .eq. 'r') then

            call MPI_IRECV( LocalMatrix(:,tmp(1):tmp(2)), dims(1)*(tmp(2)-tmp(1)+1), MPI_CTYPE,&
                     grid2rank(otherrow,othercol,rowproc,colproc) , &
                     100*grid2rank(merow,mecol,rowproc,colproc)+2*tag+ii+11 ,&  !tag is the receiver, so me
                     MPI_COMM_WORLD, requestM(ii), ierr)

            call mpic_checkerror(ierr)
        else
            print*, "Mod error"
        end if
        !print*, "good " ,  tmp(1), tmp(2) , "-", dimstmp(2),'||' , ii, "of" , out , '||'
        tmp(1)= tmp(2)+1

        end do


    end subroutine

!Given new process grid decide which processes will exchange pieces of the distributed matrix
!Various cases, sa send A, ra recv A etc..
!If rowproc > colproc you do not need copies of B. it is guaranteed that at least one process has one. hence negative id
!are not dispatched, the same for A in the case rowproc < colproc.
!The process ij looking at the new process grid understand who to send to  or receive from the piece of this step.
! ra is receive a , sa send a etc etc..
    subroutine Find_coop(old_proc_grid, new_proc_grid,idxrow,idxcol,otherrow, othercol, rowproc, colproc, mod )
    implicit none
    integer , dimension (:,:) :: old_proc_grid, new_proc_grid
    integer :: idxrow, idxcol, otherrow, othercol, rowproc, colproc
    character*2:: mod


    !A PART
    !A has to stay in the same column.
    if(mod .eq. "sa") then
        otherrow = idxrow
        if(rowproc .gt. colproc) then
            othercol = findval(abs(new_proc_grid(idxrow, 1:colproc)), abs(old_proc_grid(idxrow,idxcol)))

        else
            if (old_proc_grid(idxrow,idxcol) .gt. 0) then
                othercol = findval(new_proc_grid(idxrow, 1:colproc), old_proc_grid(idxrow,idxcol))
            else
                otherrow =0
                othercol =0
            end if
            !if(new_proc_grid(othercol, othercol) .lt. 0 ) then
            !    otherrow =0
            !    othercol =0
            !end if
            if (othercol .lt. 0) then
            print*, "error coop not found for sending A "
            end if
        end if





    else if(mod .eq. "ra") then
        otherrow = idxrow
        if(rowproc .gt. colproc) then
            othercol = findval(abs(old_proc_grid(idxrow, 1:colproc)), abs(new_proc_grid(idxrow,idxcol)))

        else
            if (new_proc_grid(idxrow,idxcol) .gt. 0) then
                othercol = findval(old_proc_grid(idxrow, 1:colproc), new_proc_grid(idxrow,idxcol))
            else
                otherrow = 0
                othercol = 0
            end if
            !if(new_proc_grid(othercol, othercol) .lt. 0 ) then
            !    otherrow =0
            !    othercol =0
            !end if
            if (othercol .lt. 0) then
            print*, "error coop not found for receiving A "
            end if
        end if




    !B part
    !B stays in the same row
    else if (mod .eq. "sb") then
        othercol = idxcol
        if(rowproc .lt. colproc) then
            otherrow = findval(abs(new_proc_grid(1:rowproc, idxcol)), abs(old_proc_grid(idxrow,idxcol)))
            if (otherrow .lt. 0) then
            print*, "error coop not found for receiving A "
            end if
        else
            if (old_proc_grid(idxrow,idxcol) .gt. 0) then
                otherrow = findval(new_proc_grid(1:rowproc, idxcol), old_proc_grid(idxrow,idxcol))
            else
                otherrow = 0
                othercol = 0
            end if

            if (otherrow .lt. 0) then
            print*, "error coop not found for receiving A "
            end if

        end if



    else if(mod .eq. "rb") then
        othercol = idxcol
        if(rowproc .lt. colproc) then
            otherrow = findval(abs(old_proc_grid(1:rowproc, idxcol)), abs(new_proc_grid(idxrow,idxcol)))

        else
            if (new_proc_grid(idxrow,idxcol) .gt. 0) then
                otherrow = findval(old_proc_grid(1:rowproc, idxcol), new_proc_grid(idxrow,idxcol))
            else
                otherrow =0
                othercol =0
            end if

            if (otherrow .lt. 0) then
            print*, "error coop not found for receiving B "
            end if

        end if


    else
        print*, "not valid mod to coop with a process"
    end if



    end subroutine



!it does this :
	! Matrix(row, col ) = val
!but matrix is distributed and val must be available to all processes
!if val is own by a specific process use cbp_setgetprcval

subroutine cbp_setval2D(localMat, descriptor, row,col, val)

        CTYPE, dimension(:,:) :: localMat
        CTYPE :: val
        integer, dimension(4) :: descriptor

        integer :: row, col

        if ( (descriptor(3) .le. row) .and. (descriptor(3)+descriptor(1) .gt. row)) then
                    if ( (descriptor(4) .le. col) .and. (descriptor(4)+descriptor(2).gt. col) ) then

                        localmat(row-descriptor(3)+1, col-descriptor(4)+1) = val
                    end if
        end if



end subroutine

!given an entry of the distributed matrix it makes it available to
!all processes in the variable val

subroutine cbp_getval1D(array,globalsize,descriptor,idx,val,rank,nprocs,root)
integer :: idx , ierr,globalsize,rank, nprocs
integer, optional :: root
integer , dimension(2) :: descriptor
CTYPE , dimension(:) :: array
CTYPE :: val
if(.not. present(root)) then
root = (idx-1)/ceiling(real(globalsize)/real(nprocs))
	if ( idx .ge. descriptor(2) .and.  idx .lt. sum(descriptor) ) then
		if (rank .ne. root) then
			print*, "Error getval1D, Unexpected process"
		end if

		call MPI_BCAST(array(idx - descriptor(2) +1 ),1, MPI_CTYPE, root ,MPI_COMM_WORLD,ierr )
		val = array(idx - descriptor(2) +1 )
	else
		call MPI_BCAST(val,1, MPI_CTYPE, root ,MPI_COMM_WORLD,ierr )
	end if

else
	if(rank.eq.root) then

		call MPI_BCAST(array(idx - descriptor(2) +1 ),1, MPI_CTYPE, root ,MPI_COMM_WORLD,ierr )
		val = array(idx - descriptor(2) +1 )
	else
		call MPI_BCAST(val,1, MPI_CTYPE, root ,MPI_COMM_WORLD,ierr )

	endif

end if

end subroutine



subroutine cbp_setval1D(array,descriptor,idx,val,rank,nprocs)
integer :: idx ,rank,nprocs
integer , dimension(2) :: descriptor
CTYPE , dimension(:) :: array
CTYPE :: val
if ( idx .ge. descriptor(2) .and.  idx .lt. sum(descriptor) ) then
	array(idx - descriptor(2) +1) = val
end if

end subroutine





!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!----------------------------------LINEAR ALGEBRA -------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!Transpose distributed matrix.

    subroutine cbp_transpose(localMat, descriptor)

        CTYPE , dimension(:,:) :: localMat
        !CTYPE,  dimension(:,:) allocatable :: tempmat
        integer , dimension(4) :: descriptor, descriptort
        integer :: row, col, ii


        if( size(localMat, dim = 1 ) .eq. size(localMat, dim= 2 ))  then
                !transpose matrix


                !Maybe in the future it is possible to parallelize this
                localMat= transpose(localmat)
                descriptort = descriptor

                !transpose descriptor dim
                descriptor(1) = descriptort(2)
                descriptor(2) = descriptort(1)

                !transpose descriptor starting indexes
                descriptor(3) = descriptort(4)
                descriptor(4) = descriptort(3)


        else
                print*, "Error, Transposition not safe"
        end if
    end subroutine




!Adjoint distributed matrix.

    subroutine cbp_adjoint(localMat, descriptor)
        implicit none
        CTYPE , dimension(:,:) :: localMat
        !CTYPE,  dimension(:,:) allocatable :: tempmat
        integer , dimension(4) :: descriptor, descriptortemp
        integer :: row, col, ii
       if( size(localMat, dim = 1 ) .eq. size(localMat, dim= 2 ))  then

#if defined(CBPCMPLX) | defined(CBPZMPLX)
                localMat = conjg(localMat)
#endif

                localMat = transpose(localMat)


                descriptortemp = descriptor
                !transpose descriptor dim

                descriptor(1) = descriptortemp(2)
                descriptor(2) = descriptortemp(1)

                !transpose descriptor starting indexes
                descriptor(3) = descriptortemp(4)
                descriptor(4) = descriptortemp(3)
		!print*, "desc" , descriptort,  descriptor

        else
                print*, "Warning, Adjoint not safe" , size(localMat, dim = 1 ),size(localMat, dim = 2 )
        end if
    end subroutine





!The matrix multipliction

!To use it be sure that the descriptors are suitable otherwise generate them with cbp_matmul_descriptors, and use cbp_redispatch
!to move the matrices according to the new descriptors
!then roproc, colproc are sides of the process grid, rank mpi rank, process grid is a matrix of integers

!Lets say that you have nprocs mpi processes, you choose rowproc and colproc s.t. rowproc*colproc = nprocs
!then generate the process_grid
!handle is the pointer to the cublas xt context

!the result C is a distributed matrix , described by descriptor C
!for the input is the same

subroutine cbp_matmul(localA, localB, localC, descriptorA, descriptorB , descriptorC, rowproc, colproc, rank, process_grid,&
                 handle)
    use iso_c_binding
    CTYPE, dimension(:,:):: localA , localB
    CTYPE, dimension(:,:), allocatable:: localAt , localBt
    CTYPE, dimension(:,:):: localC
    integer , dimension(4) :: descriptorA, descriptorB
    integer , dimension(4) :: descriptorAt, descriptorBt

    integer :: descriptorC(4)

    integer , dimension(:,:)  :: process_grid
    integer :: rowproc, colproc, rank
    integer ::nprocs
    type(c_ptr) :: handle




    integer , dimension(:,:) , allocatable :: old_process_grid

    integer :: idxrow, idxcol
    integer :: switch, withwhorow, withwhocol
    integer, dimension(:), allocatable :: requestsendA, requestrecvA ,requestsendB, requestrecvB
    integer :: requestsendDA, requestrecvDA ,requestsendDB, requestrecvDB
    logical :: sa, ra, sb , rb

    !integer, dimension(MPI_STATUS_SIZE) :: statussendA,statussendB, statusrecvA,statusrecvB
    !integer, dimension(MPI_STATUS_SIZE) :: statussendDA,statussendDB, statusrecvDA,statusrecvDB

    integer :: ngpu
    real :: cpuratio
    integer :: ii, pace
    integer :: ierr ! MPI error handler
    integer :: outv(2)

    integer , dimension(:,:), allocatable :: stat




    ierr = 0
    if ( size(process_grid, dim= 1) .ne. max(rowproc,colproc)) then
    	print*, "Error, invalid shape of process grid"
    end if
     if ( size(process_grid, dim= 2) .ne. max(rowproc,colproc)) then
    	print*, "Error, invalid shape of process grid"
    end if
    allocate (old_process_grid(max(rowproc,colproc), max(rowproc, colproc)))


    !end here
    nprocs = colproc*rowproc
    pace = max(colproc,rowproc)
    switch = 1

    allocate(stat(pace,8))
    allocate(localAt(size(localA, dim =1), size(localA, dim =2)))
    allocate(localBt(size(localB, dim =1), size(localB, dim =2)))

    !Each process understand its identity
    idxcol= rank/rowproc+1            !grid assigned by cols
    idxrow= mod(rank,rowproc)+1

    localC = CTYPE0
    if(process_grid(idxrow,idxcol) .gt. 0  ) then

    call cbpfor_auto_gemm(localA(:descriptorA(1),:descriptorA(2)),&
         		  localB(:descriptorB(1),:descriptorB(2)),&
          		  localC(:descriptorC(1), :descriptorC(2)), &
          		  1, 1 , 0.0 , handle)



    end if



     outv(1) = mpi_norequest(localA)!A
     outv(2) = mpi_norequest(localB)!B



    do ii = 2,pace

        !Step UPDATE-------------------------------
        old_process_grid = process_grid

        sa = .false.
        sb = .false.
        ra = .false.
        rb = .false.


        call update_process_grid(process_grid, rowproc, colproc,rank);


        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !FIND COOP and SEND
        !############################################################################################################!







        !send B -------------------------------------------------------------------------------------------------------
        call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "sb")
        if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
            if((withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol)) then
		    if( switch .eq. 1) then

		        LocalBt = LocalB
		        descriptorBt = descriptorB

		    else if ( switch .eq. 2 ) then

		    	LocalB = LocalBt
		        descriptorB = descriptorBt

		    end if

            else
            sb = .true.
            allocate(requestsendB(outv(2)))
            requestsendB = 0
		    if( switch .eq. 1) then
			    	call  MPI_comm_wrap_matmul(LocalB,descriptorB,idxrow,idxcol,&
				withwhorow, withwhocol,rowproc,colproc,'s',1, requestsendB, requestsendDB, outv(2))

		    else if ( switch .eq. 2 ) then
		    		call  MPI_comm_wrap_matmul(LocalBt,descriptorBt,idxrow,idxcol,&
				withwhorow, withwhocol,rowproc,colproc,'s',1, requestsendB, requestsendDB, outv(2))
		    end if

            end if
        !    print*, "sending"
        end if
        !############################################################################################################!


        !############################################################################################################!
        !send A--------------------------------------------------------------------------------------------------------
        call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "sa")
        !print*, idxrow, idxcol, "->", withwhorow, withwhocol ,"sa"
        if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then

            if((withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol)) then
                    if( switch .eq. 1) then

		        LocalAt = LocalA
		        descriptorAt = descriptorA

		    else if ( switch .eq. 2 ) then

		    	LocalA = LocalAt
		        descriptorA = descriptorAt

		    end if
                !print*, "Already Here send A"
            else
            sa = .true.
           allocate(requestsendA(outv(1)))
  	   requestsendA = 0

            	    if( switch .eq. 1) then

		                call  MPI_comm_wrap_matmul(LocalA,descriptorA,idxrow,idxcol,&
               		 withwhorow, withwhocol,rowproc,colproc,'s',0,requestsendA, requestsendDA, outv(1))

		    else if ( switch .eq. 2 ) then

		    		call  MPI_comm_wrap_matmul(LocalAt,descriptorAt,idxrow,idxcol,&
               		 withwhorow, withwhocol,rowproc,colproc,'s',0,requestsendA, requestsendDA, outv(1))


            	    end if
            end if
        end if
        !############################################################################################################!



        !@@@@@@@@@@@@@@@@@@@@@@
        switch = mod(switch,2)+1
        !@@@@@@@@@@@@@@@@@@@@@@

        !FIND COOP and RECV
        !############################################################################################################!
        !recv B--------------------------------------------------------------------------------------------------------
        call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "rb")
        if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
            if(.not. ( (withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol) ) ) then
                rb = .true.
                allocate(requestrecvB(outv(2)))
                requestrecvB = 0

               ! call MPI_comm_wrap_matmul(LocalB(:,:,switch),descriptorB(:,switch),idxrow,idxcol,&
                !    withwhorow, withwhocol,rowproc,colproc,'r',1, requestrecvB, requestrecvDB)

                 if( switch .eq. 1) then

		    	call  MPI_comm_wrap_matmul(LocalB,descriptorB,idxrow,idxcol,&
		        withwhorow, withwhocol,rowproc,colproc,'r',1, requestrecvB, requestrecvDB, outv(2))

                 else if ( switch .eq. 2 ) then

		    	 call  MPI_comm_wrap_matmul(LocalBt,descriptorBt,idxrow,idxcol,&
			 withwhorow, withwhocol,rowproc,colproc,'r',1, requestrecvB, requestrecvDB, outv(2))

                 end if



            end if
        end if
        !############################################################################################################!


        !############################################################################################################!
        !recv A--------------------------------------------------------------------------------------------------------
        call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "ra")

        if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then

            if(.not. ( (withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol) ) ) then
            ra = .true.
            allocate(requestrecvA(outv(1)))
            requestrecvA = 0


            	    if( switch .eq. 1) then

		                call  MPI_comm_wrap_matmul(LocalA,descriptorA,idxrow,idxcol,&
               		 withwhorow, withwhocol,rowproc,colproc,'r',0,requestrecvA, requestrecvDA, outv(1))

		    else if ( switch .eq. 2 ) then

		    		call  MPI_comm_wrap_matmul(LocalAt,descriptorAt,idxrow,idxcol,&
               		 withwhorow, withwhocol,rowproc,colproc,'r',0,requestrecvA, requestrecvDA, outv(1))


            	    end if
            end if
        end if
        !############################################################################################################!



        !if small matrices, but default for big, using the simplest communication wrapper
        !if(sb) then
        !    call MPI_WAIT(requestsendB, MPI_STATUS_IGNORE, ierr )
        !    call MPI_WAIT(requestsendDB, MPI_STATUS_IGNORE, ierr )
        !end if
        !if(sa) then
        !    call MPI_WAIT(requestsendA, MPI_STATUS_IGNORE, ierr )
        !    call MPI_WAIT(requestsendDA, MPI_STATUS_IGNORE, ierr )
        !end if
        !if (ra) then
        !    call MPI_WAIT(requestrecvA, MPI_STATUS_IGNORE, ierr )
        !    call MPI_WAIT(requestrecvDA, MPI_STATUS_IGNORE, ierr )
        !end if
        !if (rb) then
        !    call MPI_WAIT(requestrecvB, MPI_STATUS_IGNORE, ierr )
        !    call MPI_WAIT(requestrecvDB, MPI_STATUS_IGNORE, ierr )
        !end if
	!print*,'sb',rank,sb,  ii, requestsendB
	!print*,'sa',rank, sa,  ii, requestsendA
	!print*,'rb',rank, rb,  ii, requestrecvB
	!print*,'ra',rank,ra,  ii, requestrecvA
        !Using the chunked communication wrapper
        if(sb) then
            call MPI_WAITALL( size(requestsendB,dim=1),  requestsendB, MPI_STATUSES_IGNORE, ierr )
            call MPI_WAIT(requestsendDB, MPI_STATUS_IGNORE, ierr )
            deallocate(requestsendB)
        end if
        if(sa) then
            call MPI_WAITALL( size(requestsendA,dim=1) , requestsendA, MPI_STATUSES_IGNORE, ierr )
            call MPI_WAIT(requestsendDA, MPI_STATUS_IGNORE, ierr )
            deallocate(requestsendA)
        end if
        if (ra) then
            call MPI_WAITALL( size(requestrecvA,dim=1) , requestrecvA, MPI_STATUSES_IGNORE, ierr )
            call MPI_WAIT(requestrecvDA, MPI_STATUS_IGNORE, ierr )
            deallocate(requestrecvA)
        end if
        if (rb) then
            call MPI_WAITALL(size(requestrecvB,dim=1), requestrecvB, MPI_STATUSES_IGNORE, ierr )
            call MPI_WAIT(requestrecvDB, MPI_STATUS_IGNORE, ierr )
            deallocate(requestrecvB)
        end if





        !-------------------------End Step Update ----------------------------------
        !MATMUL
        if(process_grid(idxrow,idxcol) .gt. 0  ) then



            	    if( switch .eq. 1) then

		                 call cbpfor_auto_gemm(localA(:descriptorA(1),:descriptorA(2)),&
                                  localB(:descriptorB(1),:descriptorB(2)),&
                                  localC(:descriptorC(1), :descriptorC(2)), &
                                  1, ngpu, cpuratio, handle)

		    else if ( switch .eq. 2 ) then

		    		 call cbpfor_auto_gemm(localAt(:descriptorAt(1),:descriptorAt(2)),&
                                  localBt(:descriptorBt(1),:descriptorBt(2)),&
                                  localC(:descriptorC(1), :descriptorC(2)), &
                                  1, ngpu, cpuratio, handle)
            	    end if

        end if

        !--------------------------------------------------------------------------
    end do

    if(switch .eq. 2 .and. process_grid(idxrow,idxcol) .ne. 0 ) then

            descriptorA = descriptorAt
            descriptorB = descriptorBt
            localA = localAt
            localB = localBt

    end if

    deallocate(stat)
    deallocate(localAt)
    deallocate(localBt)

    deallocate(old_process_grid)




end subroutine
!=========================================================================================
!if you have a matrix B available in all processes and you want to compute
!distributed matrix A *B and store the result in a distributed matrix C

subroutine cbp_unbalanced_matmul(localA, descriptorA, B, localC, descriptorC, handle)
CTYPE, dimension(:,:) :: localA, B , localC
integer , dimension(4) :: descriptorA, descriptorC
type(c_ptr) :: handle
integer :: ierr, ii, jj


descriptorC = [descriptorA(1), size(B, dim = 2) , descriptorA(3), 1]
localC = 0

call cbpfor_auto_gemm(localA(:descriptorA(1),:descriptorA(2)),&
                      B,&
                      localC(:descriptorC(1), :descriptorC(2)), 1, 0,0.0, handle)


end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine GetBackMatrix(localC, descriptor, C, rowproc, colproc,&
         idxrow, idxcol,rank)

        use mpi
        implicit none


        CTYPE, dimension(:,:) :: C
        CTYPE, dimension(:,:) :: localC
        CTYPE, dimension(:,:),allocatable :: tempc

        integer, dimension(4) :: descriptor, extdescriptor
        integer :: rowproc, colproc, idxrow, idxcol, rank
        integer, dimension(MPI_STATUS_SIZE) :: status_

        integer :: ierr, ii,jj

        if(descriptor(1) .gt. size(localC, dim=1) .and. all(descriptor .ne.0) ) then
		print*, "Error, getbackmatrix, localmat too small,row", descriptor(1), size(localC, dim=1)
	end if

	if(descriptor(2) .gt. size(localC, dim=2) .and. all(descriptor .ne.0) ) then
			print*, "Error, getbackmatrix, localmat too small,col", descriptor(2), size(localC, dim=2)
	end if


        if(rank .eq. 0) then

        do ii = 1, rowproc
            do jj = 1, colproc
                if(ii*jj .ne. 1) then


                        call MPI_RECV(extdescriptor, 4, MPI_INTEGER, &
                                  grid2rank(ii,jj,rowproc,colproc), &
                                  3*(ii*(colproc+1)+jj) , MPI_COMM_WORLD, status_, ierr)

                        call mpic_checkerror(ierr)

                        allocate (tempc(extdescriptor(1),extdescriptor(2)))
                        tempc = cmplx(0.0,0.0)
                        call MPI_RECV(tempc,&
                                extdescriptor(1)*extdescriptor(2), MPI_CTYPE, &
                                grid2rank(ii,jj,rowproc,colproc), &
                                3*(ii*(colproc+1)+jj)+1 , MPI_COMM_WORLD,status_, ierr)

                        call mpic_checkerror(ierr)

                        C(extdescriptor(3):extdescriptor(3)+extdescriptor(1)-1 ,&
                          extdescriptor(4):extdescriptor(4)+extdescriptor(2)-1) =&
                                tempc

                        deallocate (tempc)


                else

                    if(descriptor(3)+descriptor(1)-1 .gt. size(C, dim=1) ) then
			print*, "Error, matrix too small redispatch", descriptor(3)+descriptor(1)-1, shape(C)
		    end if
		    if(descriptor(2)+descriptor(4)-1 .gt. size(C, dim=2) ) then
			print*, "Error, matrix too small redispatch", descriptor(2)+descriptor(4)-1 , shape(C)
		    end if


                    C(descriptor(3):descriptor(3)+descriptor(1)-1,&
                      descriptor(4):descriptor(4)+descriptor(2)-1) = localC(: descriptor(1), : descriptor(2))

                end if

            end do
        end do

        else



            call MPI_SSEND( descriptor , 4,&
                    MPI_INTEGER, 0 ,&
                    3*(idxrow*(colproc+1)+idxcol) , MPI_COMM_WORLD, ierr)

            call mpic_checkerror(ierr)


            call MPI_SSEND(localC(:descriptor(1), :descriptor(2)) , descriptor(1)*descriptor(2),&
                 MPI_CTYPE, 0 ,&
                 3*(idxrow*(colproc+1)+idxcol)+1 , MPI_COMM_WORLD, ierr)



            call mpic_checkerror(ierr)



        end if



end subroutine


!subroutine GetBackMatrix2(localC, descriptor, C, dimC rank)
!
!        use mpi
!        implicit none
!
!
!        CTYPE, dimension(:,:) :: C
!        CTYPE, dimension(:,:) :: localC
!        CTYPE, dimension(:,:),allocatable :: tempc
!
!        integer , dimension(2) :: dimC
!        integer , dimension(4) :: fakedesc
!        if rank.eq.0
!        fakedesc = [dimC(1), dimC(2), 1,1]
!        else
!        fakedesc = 0
!        end if
!        if(rank .eq. 0) then
!
!        end
!
!end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


subroutine cbp_init_local_matrix_descriptor &
        (descriptorA,dimA,localdimA,&
        rowproc, colproc, rank)

        implicit none

        integer , dimension(4)  :: descriptorA
        integer , dimension(2)  :: dimA
        integer , dimension(2)  :: localdimA
        integer :: rowproc, colproc
        integer :: idxcol, idxrow, rank

        idxcol= rank/rowproc+1
        idxrow= mod(rank,rowproc)+1

        localdimA(1) = ceiling(real(dimA(1))/real(rowproc))
        localdimA(2) = ceiling(real(dimA(2))/real(colproc))


        descriptorA =  (/maxval(localdimA), maxval(localdimA), (idxrow-1)*maxval(localdimA)+1, &
                (idxcol-1)*maxval(localdimA)+1 /)

        descriptorA =  (/localdimA(1), localdimA(2), (idxrow-1)*localdimA(1)+1, &
                (idxcol-1)*localdimA(2)+1 /)



        !---LIMIT CASES---

        !if(idxrow .eq. rowproc) then
        !        descriptorA(1) = dimA(1) - descriptorA(3) +1
        !end if
        !if(idxcol .eq. colproc) then
        !        descriptorA(2) = dimA(2) - descriptorA(4) +1
        !end if

        if(descriptorA(3)+descriptorA(1) -1 - dimA(1) .ge. 0)then
        	!print*, dimA , ([idxrow-1,idxcol-1])*(localdimA)
        	descriptorA(1) = dimA(1) - (idxrow-1)*localdimA(1)

        end if
         if(descriptorA(4)+descriptorA(2) -1 - dimA(2) .ge. 0)then
        	!print*, dimA , ([idxrow-1,idxcol-1])*(localdimA)
        	descriptorA(2) = dimA(2) - (idxcol-1)*(localdimA(2))

        end if

        if(any(descriptorA(3:) - dimA .gt. 0) )then
        	descriptorA(:2) = 0
        	descriptorA(3:) = dimA+1
        end if

        localdimA = ceiling(real(maxval(dimA)/real(min(rowproc,colproc))))



end subroutine


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!For matrices with low number of colums and hig number of rows
subroutine cbp_unbalanced_matmul_descriptors(descriptor2d, rows, cols, rank, nprocs)
integer :: rows, cols, rank, nprocs
integer, dimension(4) :: descriptor2d
descriptor2d=[ceiling(real(rows)/real(nprocs)), cols, rank*ceiling(real(rows)/real(nprocs))+1 ,1]

if (rank .eq. nprocs-1) then
	descriptor2d(1) = rows - rank*ceiling(real(rows)/real(nprocs))
end if


end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
	integer :: idxcol, idxrow, rowproc, colproc, rank
	idxcol= rank/rowproc+1
        idxrow= mod(rank,rowproc)+1
end subroutine
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine cbp_matmul_descriptors         &
        (descriptorA, descriptorB, descriptorC,&
        dimA, dimB, dimC,            &
        localdimA, localdimB, localdimC,       &
        rowproc, colproc, process_grid ,rank)

        implicit none

        integer, dimension(:,:) :: process_grid
        integer , dimension(4)  :: descriptorA, descriptorB
        integer , dimension(4)  ::descriptorC
        integer , dimension(2)  :: dimA, dimB, dimC
        integer , dimension(2)  :: localdimA, localdimB, localdimC
        integer :: rowproc, colproc
        integer :: idxcol, idxrow, rank

        idxcol= rank/rowproc+1
        idxrow= mod(rank,rowproc)+1

        if (dimA(2) .ne. dimB(1)) then
        	print*, "Error, wrong matmul dimensions"
        end if
        if (dimA(1) .ne. dimC(1)) then
        	print*, "Error, wrong matmul dimensions"
        end if
         if (dimB(2) .ne. dimC(2)) then
        	print*, "Error, wrong matmul dimensions"
        end if

        localdimA(1) = ceiling(real(dimA(1))/real(rowproc))
        localdimB(2) = ceiling(real(dimB(2))/real(colproc))
        localdimC(1) = localdimA(1)
        localdimC(2) = localdimB(2)
        localdimA(2) = ceiling(real(dimA(2))/real(min(colproc,rowproc)))
        localdimB(1) = localdimA(2)

        descriptorA =  (/localdimA(1), localdimA(2), (idxrow-1)*localdimA(1)+1, &
                (abs(process_grid(idxrow,idxcol))-1)*localdimA(2)+1 /)

        descriptorB =  (/localdimB(1), localdimB(2),&
                 (abs(process_grid(idxrow,idxcol))-1)*localdimB(1)+1, (idxcol-1)*localdimB(2)+1 /)

        descriptorC =  (/localdimC(1), localdimC(2), (idxrow-1)*localdimC(1)+1,&
                 (idxcol-1)*localdimC(2)+1/)





        !---LIMIT CASES---

        if(idxrow .eq. rowproc) then
                descriptorA(1) = dimA(1) - descriptorA(3) +1
                descriptorC(1) = dimC(1) - descriptorC(3) +1
        end if
        if(idxcol .eq. colproc) then
                descriptorB(2) = dimB(2) - descriptorB(4) +1
                descriptorC(2) = dimC(2) -descriptorC(4)+1
        end if
        if(abs(process_grid(idxrow,idxcol)) .eq. min(colproc,rowproc)) then
                descriptorA(2) = dimA(2) - descriptorA(4) +1
                descriptorB(1) = dimB(1) - descriptorB(3) +1
            if(descriptorA(2) - descriptorB(1) .ne. 0) then
                    print*, "Erroror Init_matmul_descriptors, &
                            Getting Wrong Sizes of local matrices blocks."
            end if
        end if

        if( process_grid(idxrow, idxcol) .le. 0 ) then
            if (rowproc .lt. colproc) then
                descriptorA(:2)  = 0
                descriptorA(3:)  = dimA+1
            else if (colproc .lt. rowproc) then
                descriptorB(:2)  = 0
                descriptorB(3:)  = dimA+1
            end if

        end if

        localdimA = ceiling(real(maxval(dimA)/real(min(rowproc,colproc))))
        localdimB = ceiling(real(maxval(dimB)/real(min(rowproc,colproc))))
        localdimc = ceiling(real(maxval(dimC)/real(min(rowproc,colproc))))

end subroutine




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!descriptorv & descriptor res descriptor for unidimensional array
!localv localres vectors, so matrices with 1column.

subroutine cbp_vecmatmul(localA, descriptorA, localv, descriptorv, localres, descriptorres, rank, rowproc,colproc, handle, bcast)

        CTYPE, dimension(:,:) ::localA
        CTYPE, dimension(:), target, allocatable :: localv , localres
        CTYPE, dimension(:,:) , pointer :: ptrv, ptrres
        CTYPE, dimension(:), allocatable :: tmp

        integer, dimension(2) :: descriptorv, descriptorres, descriptorvtmp !Vector descriptors
        integer :: idxrow, idxcol, rowproc, colproc, rank

        integer, dimension(4) ::  descriptorA !descriptor matrices
        integer :: rowcomm, ierr
        type(c_ptr) :: handle

        logical, optional :: bcast

        descriptorres =  [ descriptorA(1) , descriptorA(3) ]
        descriptorvtmp = [ descriptorA(2) , descriptorA(4) ]

        if (size(localv) .lt. size(localA,dim=2) ) then
          allocate(tmp(size(localA,dim=2)))
          tmp(:descriptorv(1)) = localv(descriptorv(1))
          call move_alloc(tmp,localv)
        endif
        if (size(localres) .lt. size(localA,dim=1) ) then
          allocate(tmp(size(localA,dim=1)))
          tmp(:descriptorres(1)) = localres(descriptorres(1))
          call move_alloc(tmp,localres)
        endif

        call cbp_redispatch_1d(localv, descriptorv, descriptorvtmp, rank, rowproc*colproc)

        !allocate(ptrv(descriptorv(1),1))
        !allocate(ptrres(descriptorres(1),1))

        ptrv(1:descriptorv(1),1:1) => localv(1:descriptorv(1))
        ptrres(1:descriptorres(1),1:1) => localres(1:descriptorres(1))

        descriptorv = descriptorvtmp
        localres = CTYPE0
      !  print*, "testptr", sum(abs(ptrres(:,1))), sum(abs(localres)), descriptorres, size(localres)

        call cbpfor_auto_gemm( localA(:descriptorA(1),:descriptorA(2)),&
                               ptrv, ptrres,1, 1, 0.0, handle, &
                               [descriptorA(1),descriptorA(2)], [descriptorv(1),1], [descriptorres(1),1])

        !allocate(tmp(descriptorA(1)))
        !tmp =  matmul(localA(:descriptorA(1),:descriptorA(2)), localv(1:descriptorv(1)))

        !print*, "testptr", sum(abs(ptrv(:,1))), sum(abs(localv))

        call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
        call MPI_Comm_split(MPI_COMM_WORLD, idxrow, idxcol, rowcomm, ierr)
        CALL mpi_allreduce(MPI_IN_PLACE, localres, descriptorres(1), MPI_CTYPE, MPI_SUM, rowcomm, ierr)


      !  print*, "testptr", sum(abs(ptrres(:,1))), sum(abs(localres))
        !print*, ierr
        !IF (idxcol .eq. 1) THEN
        !  CALL mpi_reduce(MPI_IN_PLACE, localres, descriptorres(1), MPI_CTYPE, MPI_SUM, 1, rowcomm, ierr)
        !  print*, "bad"
        !  if (present(bcast)) then
        !  if (bcast) call MPI_BCAST( localres, descriptorres(1),  MPI_CTYPE, 1 ,rowcomm,ierr )

        ! end if
        !ELSE
        ! CALL mpi_reduce(localres, localres, descriptorres(1), MPI_CTYPE, MPI_SUM, 1, rowcomm, ierr)
        ! print*, "good"
        !if (present(bcast)) then
        !if (bcast) call MPI_BCAST( localres, descriptorres(1),  MPI_CTYPE, 1 ,rowcomm,ierr )
        !end if
        !END IF

      !  deallocate(ptrv)
      !  deallocate(ptrres)
        return

end subroutine

!#############################################

!function to do

!dest( startdest : enddest : stridedest ) = src ( startsrc : endsrc )

!debugged, it works.
subroutine cbp_AssignRange(dest, startdest, enddest,stridedest, src, startsrc, endsrc, descriptordest, descriptorsrc,rank,nprocs)
	integer ::  startdest, startsrc, enddest, endsrc,  stridedest, rank, nprocs
	integer, dimension(2) :: otherdescriptor, correctlims
	CTYPE, dimension(:) :: src, dest
	integer :: ierr, localend , localstart
	integer :: idxdest, idxsrc, ii, jj
	logical :: strideok
	integer , dimension(2) :: descriptorsrc , descriptordest

	!descriptor contains dimension , starting_idx
	!@@@@@@@@@@@@@@@ CHECKS @@@@@@@@@@@@@@@@@@@@@@@@
	!if(size(dest) .lt. descriptordest(1)) then
	!	print*, "Error, dimension of local tensor too small. , " ,&
	!	 size(dest) , descriptordest(1)
	!end if
	!if(size(src) .lt. descriptorsrc(1)) then
	!	print*, "Error, dimension of local tensor too small.", &
	!	size(src) , descriptorsrc(1)
	!end if
	!if ( (enddest-startdest+1)/stridedest .ne. (endsrc-startsrc+1)) then
	!	print*, 'cbp assign_ranges , sure of the ranges? ' , (enddest-startdest+1)/stridedest &
	!	, (endsrc-startsrc+1)
	!end if


	!@@@@@@@@@@@@@ COPY @@@@@@@@@@@@@@@@@@@@@@@@@@@


	do ii = 0, nprocs-1
		do jj = 0, nprocs-1

			correctlims(1) = huge(correctlims(1))
			correctlims(2) = 0

			!#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			if ( rank .eq. ii ) then

				if ( jj .ne. ii ) then

					call MPI_SEND( descriptorsrc, 2, MPI_INTEGER, jj, nprocs, MPI_COMM_WORLD, ierr ) !send src of ii
					call MPI_RECV( otherdescriptor , 2, MPI_INTEGER, jj, nprocs+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) !Receive the dest of jj
				else
					otherdescriptor = descriptordest !me w me
				end if

				!print *, ii, jj, descriptorsrc, otherdescriptor
			!-----------------------------------------------
			!Check if src localmat contains something to pass
			!ie overlap between descriptor and source
			!print*, min( endsrc, sum(descriptorsrc)-1) - max( startsrc, descriptorsrc(2) )  &
			!, min( enddest, sum(otherdescriptor)-1) - max(startdest, otherdescriptor(2)),ii,jj,'u'

				if( min( endsrc, sum(descriptorsrc)-1) - max( startsrc, descriptorsrc(2) ) .ge. 0 ) then
				!check if dest_descriptor has overlap with global_destination
				!print*, ii, jj, 'u'
				if( min( enddest, sum(otherdescriptor)-1) - max(startdest, otherdescriptor(2)) .ge. 0 ) then
				!Good, check for overlap between the parts we contain and pass
				!local vars are rescaled to startsrc and startdest  = 1

				!localend =  min( min(endsrc, sum(descriptorsrc)-1 ) -startsrc +1 ,&
				!		 min(enddest, sum(otherdescriptor)-1 ) - startdest +1)
				!localstart = max( max(startsrc, descriptorsrc(2) )-startsrc +1 ,&
				!		 max(startdest, otherdescriptor(2) ) - startdest +1)

				localend =  min(endsrc, sum(descriptorsrc)-1)  -startsrc +1
				localstart =  max(startsrc, descriptorsrc(2)) -startsrc +1

				!print*, localstart, localend, "s", ii, jj
				!----
				!if ( localend - localstart .ge. 0) then
				!---
				if (stridedest .gt. 1) then
					idxsrc = 1
					strideok = .false.
					do idxdest = startdest, enddest, stridedest

						if ( idxdest .le. sum(otherdescriptor) -1 .and. &
						     idxdest .ge. otherdescriptor(2) .and. &
						     idxsrc .ge. localstart .and. &
						     idxsrc .le. localend ) then

						     correctlims(1) = min( correctlims(1), idxsrc)
						     correctlims(2) = max( correctlims(2), idxsrc)
						     strideok = .true.

					  	end if
					idxsrc = idxsrc +1
					end do
					!print*, correctlims, ii, jj, 'u'
				else if (stridedest .eq. 1) then
					correctlims(1) = localstart
					correctlims(2) = localend
					strideok = .true.
				else
					print*, "Error at 'assing_range' routine: Invalid Stride"
					strideok = .false.
				end if
				!---

				!---
				if ( jj .ne. ii .and. strideok ) then

					!to rescale to descriptor coordinates we need descriptor(2) = 1
					call MPI_SEND(src(correctlims(1) -(descriptorsrc(2) -startsrc) : correctlims(2) - (descriptorsrc(2) - startsrc)) ,  &
							correctlims(2)-correctlims(1)+1, &
							MPI_CTYPE, jj, ii, MPI_COMM_WORLD, ierr )

					!print*, correctlims(2)-correctlims(1)+1, "S",  ii, jj



				else if (strideok)  then
				!Different change of indexes
					dest(localstart-(otherdescriptor(2) -startdest) : localend-(otherdescriptor(2) -startdest) : stridedest) = &
						src ( localstart-(descriptorsrc(2) -startsrc) : localend-(descriptorsrc(2) -startsrc) )

				end if
				!---

				!end if
				!----
				end if


				end if
			end if

			!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			if( (rank .eq. jj) .and. (jj .ne. ii) ) then

				call MPI_RECV( otherdescriptor , 2, MPI_INTEGER, ii, nprocs, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
				call MPI_SEND( descriptordest, 2, MPI_INTEGER, ii, nprocs+1, MPI_COMM_WORLD, ierr ) !send src of ii

				!print*, min( endsrc, sum(otherdescriptor)-1) - max( startsrc, otherdescriptor(2) )  &
				!, min( enddest, sum(descriptordest)-1) - max(startdest, descriptordest(2)),ii,jj,'d'



			!Same checks as above
			if( min( endsrc, sum(otherdescriptor)-1) - max( startsrc, otherdescriptor(2) ) .ge. 0  ) then
			!check if dest_descriptor has overlap with global_destination

			if( min( enddest, sum(descriptordest)-1) - max(startdest, descriptordest(2) ) .ge. 0 ) then
			!-----------------------------
			!Check overlap of the right stuff to be sended
				!local vars are rescaled to startsrc and startdest  = 1
				!localend =  min( min(endsrc, sum(otherdescriptor)-1 ) - startsrc +1 ,&
				!		 min(enddest, sum(descriptordest)-1 ) - startdest +1)
				!localstart = max( max(startsrc, otherdescriptor(2) ) - startsrc +1 ,&
				!		 max(startdest, descriptordest(2) )  - startdest +1)

				localend =  min(endsrc, sum(otherdescriptor)-1) -startsrc +1
				localstart =  max(startsrc, otherdescriptor(2)) -startsrc +1

				!print*, endsrc - startsrc +1,   sum(otherdescriptor)-1 - startsrc +1  , &
				!        enddest - startdest +1,  sum(descriptordest)-1 - startdest +1
			        !print*, localstart, localend, "A", ii, jj

				!if ( localend - localstart .ge. 0) then
				!correct stride

					if (stridedest .gt. 1) then
					correctlims(1) = huge( correctlims(1))
					correctlims(2) = 0
					idxsrc = 1
					strideok = .false.
					do idxdest = startdest, enddest, stridedest

						if ( idxdest .le. sum(descriptordest) -1 .and. &
						     idxdest .ge. descriptordest(2) .and. &
						     idxsrc .ge. localstart .and. &
						     idxsrc .le. localend ) then

						     correctlims(1) = min( correctlims(1), idxdest)
						     correctlims(2) = max( correctlims(2), idxdest)
						     strideok = .true.

					  	end if
						idxsrc = idxsrc +1
					end do

						correctlims(1) = correctlims(1) - descriptordest(2) +1
						correctlims(2) = correctlims(2) - descriptordest(2) +1

					else if (stridedest .eq. 1) then
						correctlims(1) = localstart-(descriptordest(2) -startdest)
						correctlims(2) = localend-(descriptordest(2) -startdest)
						strideok = .true.
					else
						print*, "Error at 'assing_range' routine: Invalid Stride"
						strideok = .false.
					end if


					if (strideok) then
						!print*, ((correctlims(2)-correctlims(1))/stridedest)  +1 ,'A' , ii, jj

						call MPI_RECV( dest (  correctlims(1) : &
								  correctlims(2) : stridedest ) , &
								  ((correctlims(2)-correctlims(1))/stridedest)  +1, &
								  MPI_CTYPE, ii, ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )



					!end if

				end if
			!-------------------------
			end if
			end if
			end if
			!-------------------------

			CALL MPI_BARRIER (MPI_COMM_WORLD, ierr)


		end do
	end do
end subroutine


!===================================================================================
!===================================================================================

!QR decomposition 4 distributed matrices, A input, QR , G the givens rotation matrix
!A at output will be the R matrix

subroutine cbp_QR_decompose_Givens(A, descriptorA, Q, descriptorQ, dimA, dimQ,rowproc, colproc,rank,nprocs, context )

	CTYPE, dimension(:,:) :: A,Q
	CTYPE, dimension(:,:) ,allocatable :: Qt, G  ! Qtest ,Qtt, Gtest, Id, Qttest
	CTYPE, dimension(:,:), allocatable :: currentrows, newrows,currentrows2, newrows2!, newwrows
	CTYPE, dimension(:,:,:), allocatable :: mod_matrix
	CTYPE, dimension(:,:), allocatable :: mod_matrix2
	CTYPE :: valtest
	double precision :: rad


	integer, dimension(4) :: descriptorA, descriptorQ, descriptorG, descriptorQt
	integer, dimension(4) :: descriptorQtemp , descriptorGtemp
	integer, dimension(4) :: descriptorr , descriptori, descriptorres

	integer, dimension(2) :: dimA, dimQ, ld1, ld2, ld3
	integer, dimension(:,:), allocatable :: process_grid
	integer :: idxrow, idxcol, nprocs, colproc, rowproc, rank, saferow, safecol , ierr,kk, self



	integer :: ii , jj, maxdiag, howmanyparallel,parallelxproc, meparallel, ia, ja, yy

	type(c_ptr)  :: context

	if (any(descriptorA .lt. 0)) then
		print*, "error, invalid R/A descriptor:" , descriptorA
	end if
	if (any(descriptorQ .lt. 0)) then
		print*, "error, invalid Q descriptor :", descriptorQ
	end if

	!Checks
	if ( (dimA(1) .ne. dimQ(1) ).or. (dimA(1) .ne. dimQ(2))  ) then
		print*, "Wrong dimension of matrices in Qr "
	end if

	!-----------------------Init Params -------------------------------

	allocate (process_grid(max(rowproc,colproc), max(rowproc,colproc)))
    allocate(G(size(Q,dim=1), size(Q,dim=1))) ! Givens rotation Matrix
	descriptorG = descriptorQ
	allocate(Qt(size(Q,dim=1), size(Q,dim=2)))





	!allocate(Qtest(dimA(1), dimA(2)))
	!allocate(Gtest(dimQ(1), dimQ(2)))
	!allocate(Qttest(dimQ(1), dimQ(2)))


	!Id = CTYPE0
	!do ii = 1, dimQ(1)
	!Id(ii,ii) = CTYPE1
	!end do
	!allocate(Qtt(size(Q,dim=1), size(Q,dim=2))) !Temp space of updates to Q matrix
	descriptorG = descriptorQ


	!Q matrix at start as identity
	Q = CTYPE0
	G=CTYPE0

	!Qtt = CTYPE0
	!do ii = 1, dimQ(1)
	!Qtt(ii,ii) = CTYPE1
	!end do

	!see if is possible to change theese with cbp_identity
	do ii = 1,dimQ(1)
		call cbp_setval2D(Q, descriptorQ, ii , ii, CTYPE1 )
		call cbp_setval2D(G, descriptorG, ii , ii, CTYPE1 ) !G descriptor = Q descriptor
	end do





	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	do  ii = 3, minval(dimA) + dimA(1) !Big sequential cycle, ii represent a chunk of parallizable op.
					    !Basically each iteration is an antidiagonal of the matrix.


	    !-----------Handle parallelism----------------------------------------
	    !Need to understand each ii how many parallel operation contains
	    !and decide how many parallel op to each process.

	    howmanyparallel = (ii-1)/2
	    saferow = 0
	    safecol = 0

	    if ( ii-1 .gt. dimA(1)  ) then
	    	saferow = ( ii - dimA(1) -1 )
	    	howmanyparallel = howmanyparallel - saferow
	    end if
	    !maxcol = (ii -1)/2
	    if ( ii .gt. 2*dimA(2) ) then
	    	safecol = (ii - 2*dimA(2) -1 )/2
	    	howmanyparallel = howmanyparallel -safecol

	    end if


	    parallelxproc = ceiling(real(howmanyparallel)/real(nprocs))
	    meparallel = parallelxproc

	    if( parallelxproc*(rank+1) .gt. howmanyparallel ) then
	    	meparallel = howmanyparallel - parallelxproc*rank
	    end if

	    if ( meparallel .le. 0 ) then
	    	meparallel = 0
	    end if


	    allocate(currentrows(2*parallelxproc, dimA(2)))
	    allocate(newrows(2*parallelxproc, dimA(2)))
	   !  allocate(currentrows2(2*parallelxproc, dimA(2)))
	   ! allocate(newrows2(2*parallelxproc, dimA(2)))


		!print*, "hh", ii, rank ,howmanyparallel, meparallel
	    !----------------------------------------------------------------------
	    !Sequential cycle to get the elements for each process from the distributed matrix A
	    !exchange of lines.
	    !Once evry process has the line it needs and then we can compute the new lines
	    !----------------------------------------------------------------------
	    do jj = 1, howmanyparallel

	    	 ia = ii - jj -saferow
	    	 ja = jj +saferow



	    	 call cbp_setgetprcline(A, descriptorA, currentrows(2*mod((jj-1),parallelxproc)+1,: ), &
	    	 			[1, dimA(2), ja, 1], (jj-1)/parallelxproc, rank, nprocs,'get')

	    	 call cbp_setgetprcline(A, descriptorA, currentrows(2*mod((jj-1),parallelxproc)+2,: ), &
	    	 			[1, dimA(2), ia, 1], (jj-1)/parallelxproc, rank, nprocs,'get')



	    end do

	    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	    !do jj = 1, 2*meparallel
	    !	print*,rank, 's', currentrows(jj,:)
	    !end do
	    !parallel in the sense that this cicle is divided in the mpi processes
	    !Maybe split is in gpus

	    !allocate(mod_matrix(2,2,meparallel))
	    allocate(mod_matrix2(2*parallelxproc,2*parallelxproc))

	     newrows = CTYPE0
	    !------------------------------Comment-------------------------------------------------
	    !Parallel , each process perform a different for cycle and computes the QR coefficients,
	    !store them in mod matrix and compte new rows for the A/R matrix

	    !It can become interesting in the future to parallelize this cycle with openmp
	    !and use not cublasxt but with simple cublas this way each openmp thread use a single gpu
	    ! and compute a single cycle of this for loop
	    !----------------------------------------------------------------------------------------

	    mod_matrix2 = CTYPE0
	    do jj =  1, meparallel!(i/2)!range(int(i/2)+1) : #parallel

		ia = ii - jj - (parallelxproc*rank) - saferow !row
		ja = jj + (parallelxproc*rank) + saferow !column



		self = 2*(jj-1)
#if defined(CBPCMPLX) | defined(CBPZMPLX)

		rad = sqrt( currentrows(2*(jj-1)+2,ja)*conjg(currentrows(2*(jj-1)+2,ja))&
			+currentrows(2*(jj-1)+1,ja)*conjg(currentrows(2*(jj-1)+1,ja)) )


		!if (rad .ge. 1e-7) then
		!mod_matrix(1,1,jj) = conjg(currentrows(2*(jj-1)+1,ja))/rad !upperleft c
		!mod_matrix(2,2,jj) = currentrows(2*(jj-1)+1,ja)/rad !lowerright c
		!mod_matrix(1,2,jj) = conjg(currentrows(2*(jj-1)+2,ja))/rad !upperright -s
		!mod_matrix(2,1,jj) = -currentrows(2*(jj-1)+2,ja)/rad !lowerleft s
		!else
		!mod_matrix(:,:,jj) = CTYPE0
		!mod_matrix(1,1,jj) = CTYPE1
		!mod_matrix(2,2,jj) = CTYPE1
		!end if


		if (rad .ge. 1e-7) then
		mod_matrix2(1+self,1+self) = conjg(currentrows(2*(jj-1)+1,ja))/rad !upperleft c
		mod_matrix2(2+self,2+self) = currentrows(2*(jj-1)+1,ja)/rad !lowerright c
		mod_matrix2(1+self,2+self) = conjg(currentrows(2*(jj-1)+2,ja))/rad !upperright -s
		mod_matrix2(2+self,1+self) = -currentrows(2*(jj-1)+2,ja)/rad !lowerleft s
		else
		mod_matrix2(1+self:2+self,1+self:2+self) = CTYPE0
		mod_matrix2(1+self,1+self) = CTYPE1
		mod_matrix2(2+self,2+self) = CTYPE1
		end if
#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)
		rad = sqrt( currentrows(2*(jj-1)+2,ja)*(currentrows(2*(jj-1)+2,ja))&
			+currentrows(2*(jj-1)+1,ja)*(currentrows(2*(jj-1)+1,ja)) )

		if (rad .ge. 1e-7) then
		mod_matrix(1,1,jj) = (currentrows(2*(jj-1)+1,ja))/rad !upperleft c
		mod_matrix(2,2,jj) = currentrows(2*(jj-1)+1,ja)/rad !lowerright c
		mod_matrix(1,2,jj) = (currentrows(2*(jj-1)+2,ja))/rad !upperright -s
		mod_matrix(2,1,jj) = -currentrows(2*(jj-1)+2,ja)/rad !lowerleft s
		else
		mod_matrix(:,:,jj) = CTYPE0
		mod_matrix(1,1,jj) = CTYPE1
		mod_matrix(2,2,jj) = CTYPE1
		end if
#endif


	       !Set R

	       !newwrows = CTYPE0
	       !call cbpfor_auto_gemm(mod_matrix(:,:,jj),currentrows(2*(jj-1)+1:2*(jj-1)+2,:),&
          	!	  newrows(2*(jj-1)+1:2*(jj-1)+2,:), &
          	!	  0 , 0 , 0.0 , context)

	    end do





	     call cbpfor_auto_gemm(mod_matrix2,currentrows,&
          		  newrows, &
          		  0 , 0 , 0.0 , context)

        !do jj =  1, meparallel
	    !	self = 2*(jj-1)
	    !	mod_matrix(:,:,jj) = mod_matrix2(1+self:2+self,1+self:2+self)

	    !end do


	    !newrows = newrows2
	    !print*, "test", sum(newrows2-newrows)

	    !for cycle on all the parallel done to update the distributed G matrix and A/R matrix
	    !print*, howmanyparallel

	    do jj = 1, howmanyparallel

	        ia = ii - jj -saferow !row of G matrix (ia,ja) lower left sine
	    	ja = jj +saferow  !column


	       call cbp_setgetprcline(A, descriptorA, newrows(2*mod((jj-1),parallelxproc)+1,: ), &
	       			[1, dimA(2), ja, 1], (jj-1)/parallelxproc, rank, nprocs, 'set')

	       call cbp_setgetprcline(A, descriptorA, newrows(2*mod((jj-1),parallelxproc)+2,: ), &
	       			[1, dimA(2), ia, 1], (jj-1)/parallelxproc, rank, nprocs, 'set')

	        !print*, ia, ja
	        !print*, "mod mat -->" , sum(matmul(transpose(conjg(mod_matrix(:,:,mod((jj-1),parallelxproc)+1))) , &
	        !	mod_matrix(:,:,mod((jj-1),parallelxproc)+1))), mod((jj-1),parallelxproc)+1
	     	!G(ja,ja) = mod_matrix(1,1,mod((jj-1),parallelxproc)+1)
	     	!G(ia,ja) = mod_matrix(2,1,mod((jj-1),parallelxproc)+1)
	     	!G(ja,ia) = mod_matrix(1,2,mod((jj-1),parallelxproc)+1)
	     	!G(ia,ia) = mod_matrix(2,2,mod((jj-1),parallelxproc)+1)



		   self = 2*mod((jj-1),parallelxproc)
	       call  cbp_setgetprcval(G, descriptorG, mod_matrix2(1+self,1+self), &
	       				[ja,ja], (jj-1)/parallelxproc, rank, nprocs, 'set')
	       call  cbp_setgetprcval(G, descriptorG, mod_matrix2(2+self,1+self), &
	       				[ia,ja], (jj-1)/parallelxproc, rank, nprocs, 'set')
	       call  cbp_setgetprcval(G, descriptorG, mod_matrix2(1+self,2+self), &
	       				[ja,ia], (jj-1)/parallelxproc, rank, nprocs, 'set')
	       call  cbp_setgetprcval(G, descriptorG, mod_matrix2(2+self,2+self), &
	       				[ia,ia], (jj-1)/parallelxproc, rank, nprocs, 'set')


	       !call  cbp_setgetprcval(G, descriptorG, valtest, &
	       !				[ja,ja], (jj-1)/parallelxproc, rank, nprocs, 'get')

	       !if (rank .eq. (jj-1)/parallelxproc ) then
	       !print*,'1', valtest -mod_matrix(1,1,mod((jj-1),parallelxproc)+1)
	       !end if

	       !call  cbp_setgetprcval(G, descriptorG,valtest, &
	       !	[ia,ja], (jj-1)/parallelxproc, rank, nprocs, 'get')

	       !if (rank .eq. (jj-1)/parallelxproc ) then
	       !print*, '2', valtest -mod_matrix(2,1,mod((jj-1),parallelxproc)+1)
	       !end if

	       !call  cbp_setgetprcval(G, descriptorG, valtest, &
	       !		[ja,ia], (jj-1)/parallelxproc, rank, nprocs, 'get')

	       !if (rank .eq. (jj-1)/parallelxproc ) then
	       !print*, '3', valtest -mod_matrix(1,2,mod((jj-1),parallelxproc)+1)
	       !end if

	       !call  cbp_setgetprcval(G, descriptorG, valtest, &
	       !		[ia,ia], (jj-1)/parallelxproc, rank, nprocs, 'get')
	       !if (rank .eq. (jj-1)/parallelxproc ) then
	       !print*, '4', valtest -mod_matrix(2,2,mod((jj-1),parallelxproc)+1)
	       !end if

	    end do



	   ! deallocate(mod_matrix)
	    deallocate(currentrows)
  	    deallocate(newrows)
  	    deallocate(mod_matrix2)
	    !deallocate(currentrows2)
  	    !deallocate(newrows2)
  	   !deallocate(newwrows)

	 	!------------Update of Q--------------------------------------------------
		call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
		Qt = CTYPE0
		Qt = Q
		descriptorQt = descriptorQ
		descriptorQtemp = descriptorQ
        descriptorGtemp = descriptorG

                !Test G orthogonal --------------------------------------------------------------------
                 ! call cbp_adjoint(G, descriptorG)
                 ! Gtest = CTYPE0
                 ! call GetBackMatrix(G, descriptorG, Gtest, rowproc, colproc,&
         	 ! 	idxrow, idxcol,rank)
         	 ! if (rank.eq. 0 ) then
                 ! print*, "Test G", ( abs(sum(matmul(transpose(conjg(Gtest)), Gtest )-Id)+ &
                 ! sum(matmul(transpose(conjg(Gtest)), Gtest))&
                 ! -sum(Id)  + sum(matmul(Gtest, transpose(conjg(Gtest)) )-Id) ).lt. 1e-5)
         	 ! end if
         	 !  call cbp_adjoint(G, descriptorG)
		!--------------------------------------------------------------------------------------

		call cbp_matmul_descriptors(descriptorG, descriptorQt, descriptorQ, dimQ, dimQ, dimQ, &
		ld1, ld2, ld3 , rowproc, colproc, process_grid,rank)

		call cbp_redispatch(Qt,descriptorQtemp, descriptorQt, rank, nprocs)

		call cbp_redispatch(G,descriptorGtemp, descriptorG, rank, nprocs)

		Q = CTYPE0
		call cbp_matmul(G, Qt, Q, descriptorG, descriptorQt , descriptorQ, rowproc, colproc, rank, process_grid,&
                context)


               ! call GetBackMatrix(A, descriptorA, Qtest, rowproc, colproc,&
         	 ! 	idxrow, idxcol,rank)
         	 !if (rank.eq.0) then
         	 !do kk = 1, dimA(1)
         	 	!print*, Qtest(kk,:)
         	 !end do
         	 !print*, " "
         	!end if

         	 !call GetBackMatrix(Q, descriptorQ, Gtest, rowproc, colproc,&
         	 ! 	idxrow, idxcol,rank)
         	! Qtest = matmul(conjg(transpose(Gtest)), Qtest)
         	 !if (rank.eq.0) then
         	 !do kk = 1, dimQ(1)
         	 	!print*, Qtest(kk,:)
         	 !end do
         	 !print*, " "
         	!end if



                !Test Q orthogonal -------------------------------------------------------------------
                !  Qtest = CTYPE0
                !  Qttest = CTYPE0
                !  call GetBackMatrix(Q, descriptorQ, Qtest, rowproc, colproc,&
         	!	idxrow, idxcol,rank)
         	! call GetBackMatrix(Qt, descriptorQt, Qttest, rowproc, colproc,&
         	!	idxrow, idxcol,rank)
         	!  call GetBackMatrix(G, descriptorG, Gtest, rowproc, colproc,&
         	!  	idxrow, idxcol,rank)
         	!if(rank.eq.0 ) then
             	!  print*, "Test Q",  sum(matmul(transpose(conjg(Qttest)), Qttest )-Id), &
             	!  	sum(matmul(transpose(conjg(Gtest)), Gtest )-Id)
             	!
             	!  Qtest = matmul(Gtest, Qttest)
             	!  print*, "Test Q2",  sum(matmul(transpose(conjg(Qttest)), Qttest )-Id) !&    !+ sum(matmul(transpose(conjg(Qtest)), Qtest)) &
             	!  !- sum(Id)
             	!end if

		G = CTYPE0

		call cbp_identity(G, descriptorG)
		!do jj = 1,dimQ(1)
		!	call cbp_setval2D(G, descriptorG, jj , jj, CTYPE1 ) !G descriptor = Q descriptor
		!end do



	end do


		!deallocate(R)
#if defined(CBPCMPLX) | defined(CBPZMPLX)
		call cbp_adjoint(Q, descriptorQ)
#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)
		call cbp_transpose(Q, descriptorQ)
#endif

		deallocate(process_grid)
		deallocate(G)
		deallocate(Qt)
		!deallocate(Qtt)

end subroutine



!===================================================================================
!===================================================================================
subroutine cbp_QR_decompose_Householder(A, descriptorA, L, descriptorL, &
					dimA ,rowproc,colproc, rank, nprocs, context)
 	CTYPE , dimension(:,:) :: A, L
 	CTYPE, dimension(:,:) , allocatable :: H , Atemp, T, Atest
 	integer , dimension(4) :: descriptorA , descriptorR, descriptorL, &
 				  descriptorVt, descriptorH,&
 				  descriptorAtemp, old1, old2, helper, descriptorT
 	integer :: ierr
 	integer, dimension(2) :: ld1, ld2, ld3, dimL, dimR, dimA, descriptorV
 	integer :: ii ,jj, kk, dstr_dim, idxrow, idxcol, fromnorm, startglob,&
 		   startl, lenl, ngpu,rowa, cola, rowproc, colproc, rank, nprocs

 	CTYPE, dimension(:), allocatable :: v
 	CTYPE :: Aii
 	double precision :: lineanorm, vnorm

 	integer, dimension(:) , allocatable :: groupproc
 	integer, dimension(:,:), allocatable :: process_grid
 	type(c_ptr) :: context

 	rowA = dimA(1)
 	colA = dimA(2)
 	dimL = [dimA(1), dimA(1)]
 	dimR = [dimA(2), dimA(2)]

 	call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
 	allocate(process_grid(max(rowproc,colproc),max(rowproc,colproc)))
 	allocate (Atemp ( size(A, dim=1), size(A, dim = 2)))
 	Atemp = A
 	descriptorAtemp = descriptorA
 	call cbp_identity(L, descriptorL)


 	!---DEBUG---
 	allocate (Atest(dimA(1), dimA(2)))
 	!allocate (Htest(dimA(1), dimA(1)))
 	!Atest = CTYPE0
 	!---------


	do ii =  1,cola

	    if (  ii .le. rowA-1 ) then

		call cbp_init_local_matrix_descriptor(helper,dimL,ld1,&
							       rowproc, colproc, rank)
		!allocate(lineA(ld1(1)))
		allocate(v(ld1(1)))
		v = CTYPE0
		!lineA = CTYPE0


		do jj = 0, nprocs-1
			startglob = max(helper(3), ii)
			lenl =  helper(1)
			startl =  startglob - helper(3) +1

			call cbp_setgetprcline( Atemp , descriptorAtemp,V(startl:lenl), [lenl-startl+1,1,startglob,ii] , jj, rank, nprocs,'get')
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)

		end do



		allocate(groupproc(rowproc))

		do jj = 1, rowproc
			groupproc(jj) = grid2rank(jj,1,rowproc,colproc)
		end do


		lineAnorm = 0.0
		lineAnorm = cbp_l2norm(v,1, descriptorAtemp(1) ,rank, nprocs, groupproc )
		!print*, rank, lineAnorm
		if ( lineAnorm.gt. 1e-7 ) then


			!COMPUTE HOUSHOLDER VECTOR

			descriptorV(1) = helper(1) ! dimension used
			descriptorV(2) = helper(3) !start_point


			!v(:descriptorv(1)) = lineA(:descriptorv(1))

			Aii = CTYPE0
			Aii = cbp_getval2d(Atemp, descriptorAtemp, [ii,ii], rank, nprocs)



			if ( sum(descriptorv) -1 .ge. ii .and. descriptorv(2) .le. ii) then
				v(ii-descriptorv(2) +1 ) = v(ii-descriptorv(2) +1 ) - lineAnorm*Aii/abs(Aii)
			end if



			vnorm = cbp_l2norm(v(:descriptorv(1)),1, descriptorv(1),rank, nprocs, groupproc )
			v = v/vnorm

			!-------------COMPUTE HOSHOLDER MATRIX---------------------

			call cbp_init_local_matrix_descriptor(descriptorH,dimL,ld1,&
							       rowproc, colproc, rank)
			allocate(H(ld1(1),ld1(2)))
			H = CTYPE0

			call cbp_identity(H,descriptorH)
			H = -0.5*(H)


			call cbp_outer(v, descriptorV, v, descriptorV, H, descriptorH, dimL, 1, &
						rowproc, colproc, rank, nprocs, context  )

			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			H = -2*H


			deallocate(v)

			!-------------------------------UPDATE L --------------------------

			allocate(T(dimA(1),dimA(1)))
			T = CTYPE0

			old1 = descriptorL
			old2 = descriptorH
			process_grid = 0

			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)

			call cbp_matmul_descriptors(descriptorH, descriptorL, descriptorT,&
						    dimL, dimL, dimL, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)


			call cbp_redispatch(L, old1, descriptorL, rank, nprocs)
			call cbp_redispatch(H, old2, descriptorH, rank, nprocs)


			call cbp_matmul(H,L,T,descriptorH,descriptorL,descriptorT, &
					rowproc, colproc, rank, process_grid, context)

			L = T
			descriptorL = descriptorT


			deallocate(T)

			!-------------------------------UPDATE RESULT MATRIX A --------------------------

			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
			old1 = 0
			old2 = 0
			old1 = descriptorA
			old2 = descriptorL
			call cbp_matmul_descriptors(descriptorL, descriptorA, descriptorAtemp,&
						    dimL, dimA, dimA, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)


			call cbp_redispatch(A,old1,descriptorA, rank, nprocs)
			call cbp_redispatch(L,old2,descriptorL, rank, nprocs)


			call cbp_matmul(L,A,Atemp,descriptorL,descriptorA,descriptorAtemp, &
					rowproc, colproc, rank, process_grid, context)

			!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

			! ---DEBUG : column annihilated?---
			!Atest = CTYPE0
			!call GetBackMatrix(Atemp, descriptorAtemp, Atest, rowproc, colproc,&
         	!  		idxrow, idxcol,rank)
         	  	!----------------------------------

			!if (rank.eq.0) then
			!print*, "Annihilated ?  ", ii,  sum(Atest(ii+1:rowA,ii))
			!end if

			deallocate(H)


		else
			deallocate(v)
		end if

		!deallocate(lineA)
		deallocate(groupproc)
	    end if
	     !if (ii+1 < colA-1 and i < rowA)


	    end do

	A = Atemp
	descriptorA = descriptorAtemp
	deallocate(Atemp)
	deallocate(process_grid)

	deallocate(Atest)

end subroutine



!===================================================================================
!Has to be cleandeeed..
!===================================================================================
subroutine cbp_bidiag(A, descriptorA, R, descriptorR, L, descriptorL, &
					dimA ,rowproc,colproc, rank, nprocs, context)
 	CTYPE , dimension(:,:) :: A, R, L
 	CTYPE, dimension(:,:) , allocatable :: H , Atemp, T!, Htest, Atest, outtest
 	integer , dimension(4) :: descriptorA , descriptorR, descriptorL, &
 				              descriptorH, descriptorAtemp, descriptorT, &
 				              old1, old2, helper

 	integer, dimension(2) :: ld1, ld2, ld3, dimL, dimR, dimA, descriptorV
 	integer :: ii ,jj, idxrow, idxcol,  kk,ierr,ngpu, jju, iiu
 	integer :: startl, lenl, startglob, rowa, cola, rowproc, colproc, rank, nprocs
 	CTYPE, dimension(:), allocatable :: v!, vtest
 	CTYPE :: Aii

 	double precision :: lineanorm, vnorm
 	CTYPE, dimension(:) , allocatable :: lineA
 	integer, dimension(:) , allocatable :: groupproc
 	integer, dimension(:,:), allocatable :: process_grid
 	type(c_ptr) :: context

 	rowA = dimA(1)
 	colA = dimA(2)
 	dimL = [dimA(1), dimA(1)]
 	dimR = [dimA(2), dimA(2)]
 	allocate(process_grid(max(rowproc,colproc),max(rowproc,colproc)))
 	!---DEBUG---
 	!allocate (Atest(dimA(1), dimA(2)))
 	!allocate (Htest(dimA(2), dimA(2)))
 	!allocate (outtest(dimA(1), dimA(1)))
 	!allocate(vtest(dimA(2)))
 	!---------
 	!Atest = CTYPE0


 	call cbp_identity(L, descriptorL)
 	call cbp_identity(R, descriptorR)

 	call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)


 	allocate (Atemp ( size(A, dim=1), size(A, dim = 2)))

 	Atemp = A

 	descriptorAtemp = descriptorA
 	!do jj = 1, rowA
	!	print*, A(jj,:)
	!end do


	do ii =  1,max(rowa,cola)
	     !============================ SET COLS TO ZERO ==============================
	    if ( ii .le. colA .and. ii .le. rowA-1 ) then

		call cbp_init_local_matrix_descriptor(helper,dimL,ld1,&
							       rowproc, colproc, rank)
		allocate(lineA(ld1(1)))
		lineA = CTYPE0

		!Debug vtest
		!vtest = CTYPE0
		!do jj = 0, nprocs-1
		!call cbp_setgetprcline( Atemp , descriptorAtemp,vtest(ii:), [rowA-ii+1,1,ii,ii] , jj, rank, nprocs,'get')
		!end do

		startglob = max(helper(3), ii)
		lenl =  helper(1)
		startl =  startglob - helper(3) +1

		do jj = 0, nprocs-1
			call cbp_setgetprcline( Atemp , descriptorAtemp,lineA(startl:lenl), [lenl-startl+1,1,startglob,ii] , jj, rank, nprocs,'get')
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		end do


		allocate(groupproc(rowproc))
		do jj = 1, rowproc
			groupproc(jj) = grid2rank(jj,1,rowproc,colproc)
		end do
		lineAnorm = 0.0
		lineAnorm = cbp_l2norm(lineA,1, descriptorAtemp(1) ,rank, nprocs, groupproc )
		!print*, rank, lineAnorm
		if ( lineAnorm.gt. 1e-7 ) then


			!COMPUTE HOUSHOLDER VECTOR

			allocate(v(size(linea,dim=1))) !,2 ))
			descriptorV(1) = helper(1) ! dimension used
			descriptorV(2) = helper(3) !start_point

			v = CTYPE0
			v(:descriptorv(1)) = lineA(:descriptorv(1))

			Aii = CTYPE0
			Aii = cbp_getval2d(Atemp, descriptorAtemp, [ii,ii], rank, nprocs)



			if ( sum(descriptorv) -1 .ge. ii .and. descriptorv(2) .le. ii) then
				v(ii-descriptorv(2) +1 ) = v(ii-descriptorv(2) +1 ) - lineAnorm*Aii/abs(Aii)
			end if



			vnorm = cbp_l2norm(v(:descriptorv(1)),1, descriptorv(1),rank, nprocs, groupproc )
			v = v/vnorm


			!Debug vtest
			!vtest(ii)  = vtest(ii) - Aii*sqrt(dot_product(vtest,vtest))/abs(Aii)
			!vtest = vtest/sqrt(dot_product(vtest,vtest))
			!vnorm = 0.0
			!vnorm = cbp_l2norm(v,1, descriptorv(1),rank, nprocs, groupproc )

			!COMPUTE HOSHOLDER MATRIX

			call cbp_init_local_matrix_descriptor(descriptorH,dimL,ld1,&
							       rowproc, colproc, rank)
			allocate(H(ld1(1),ld1(2)))
			H = CTYPE0

			call cbp_identity(H,descriptorH)
			H = -0.5*(H)


			call cbp_outer(v, descriptorV, v, descriptorV, H, descriptorH, dimL, 1, &
						rowproc, colproc, rank, nprocs, context  )

			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			H = -2*H


			!DEBUG :

			!1)GET THE CORRECT HOUSHOLDER MATRIX FROM V TEST

			!call GetBackMatrix(Atemp, descriptorAtemp, Atest, rowproc, colproc,&
         	  	!	idxrow, idxcol,rank)
         	  	!call GetBackMatrix(H, descriptorH, outtest, rowproc, colproc,&
         	  	!	idxrow, idxcol,rank)

			!Htest = CTYPE0
			!do jj = 1, rowA
			!	do kk = 1, rowA
			!		Htest(jj,kk) =-2*  vtest(jj) * conjg(vtest(kk))
			!		if (jj .eq. kk) then
			!		Htest(jj,kk) = Htest(jj,kk)+1
			!		end if
			!	end do
			!end do

			!2) VERIFY IF ANNIHILATES A COLUMN OF A
			!Atest = matmul(Htest,Atest)

			!if (rank.eq.0) then
			!print*, "A ANNIHILATED FROM TEST ?  if 0.0 good : ",maxval(abs(Atest(ii+1:rowa,ii)))
			!end if




			deallocate(v)

			!-------------------------------UPDATE L --------------------------

			allocate(T(size(L,dim=1),size(L,dim=2)))
			T = CTYPE0

			old1 = descriptorL
			old2 = descriptorH
			process_grid = 0


			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)

			call cbp_matmul_descriptors(descriptorH, descriptorL, descriptorT,&
						    dimL, dimL, dimL, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)



			call cbp_redispatch(L, old1, descriptorL, rank, nprocs)
			call cbp_redispatch(H, old2, descriptorH, rank, nprocs)

			call cbp_matmul(H,L,T,descriptorH,descriptorL,descriptorT, &
					rowproc, colproc, rank, process_grid, context)


			!call printMatrix(process_grid)



			L = T
			descriptorL = descriptorT


			deallocate(T)

			!-------------------------------UPDATE RESULT MATRIX A --------------------------

			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
			old1 = 0
			old2 = 0

			A = Atemp
			descriptorA = descriptorAtemp

			old1 = descriptorA
			old2 = descriptorH
			call cbp_matmul_descriptors(descriptorH, descriptorA, descriptorAtemp,&
						    dimL, dimA, dimA, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)


			call cbp_redispatch(A,old1,descriptorA, rank, nprocs)
			call cbp_redispatch(H,old2,descriptorH, rank, nprocs)

			!Atest = CTYPE0

			call cbp_matmul(H,A,Atemp,descriptorH,descriptorA,descriptorAtemp, &
					rowproc, colproc, rank, process_grid, context)

			!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

			! DEBUG : column annihilated?
			!call GetBackMatrix(Atemp, descriptorAtemp, Atest, rowproc, colproc,&
         	!  		idxrow, idxcol,rank)
			!if (rank.eq.0) then
			!print*, "Annihilated ?  ", ii,  sum(Atest(ii+1:rowA,ii))
			!end if

			deallocate(H)
			!deallocate(Htest)
		!+++++++++++++++++++++++++++++++++++++++++
		end if !LINEANORM
		deallocate(lineA)
		deallocate(groupproc)
	    end if !condition annihilate columns
	    !+++++++++++++++++++++++++++++++++++++++++


	    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	     !============================ SET ROWS TO ZERO ==============================
	     if (ii+1 .le. colA-1 .and. ii .le. rowA) then

	     	call cbp_init_local_matrix_descriptor(helper,dimR,ld1,&
							  rowproc, colproc, rank)
		allocate(lineA(ld1(1)))
		lineA = CTYPE0

		!------------------------------- DEBUG
		!vtest = CTYPE0
		!do jj = 0, nprocs-1
		!call cbp_setgetprcline( Atemp , descriptorAtemp, vtest(ii+1:), [1,colA-ii,ii,ii+1] , jj, rank, nprocs,'get')
		!end do
	    !------------------------------
		!get ii-th rows starting from ii+1
		!and in the process grid put it in a single row, hence more rows more copies.

		startglob = max(helper(3), ii+1)
		lenl =  helper(1)
		startl =  startglob - helper(3) +1

		do jj = 0, nprocs-1
			call cbp_setgetprcline( Atemp , descriptorAtemp,lineA(startl:lenl), &
				                    [1,lenl-startl+1,ii,startglob] , jj, rank, nprocs,'get')
			!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		end do
			!if (rank.eq.jj) then
				!print*, startglob, lenl, startl,sum(helper(1:3:2))-1, helper(3), colA
				!print*, sum(lineA(startl:lenl) )
				!print*, sum(vtest(helper(3):sum(helper(1:3:2))-1))
				!print*, " "
				!print*, sum(Atemp(ii, ii+1:colA))
			!end if

		allocate(groupproc(rowproc))  !now a row owns the data

		do jj = 1, rowproc
			groupproc(jj) = grid2rank(jj,1,rowproc,colproc)
		end do

		lineAnorm = 0.0
		!print*, rank, sum(lineA), descriptorAtemp
		lineAnorm = cbp_l2norm(lineA,1, descriptorAtemp(1) ,rank, nprocs, groupproc )
		!print*, rank, lineAnorm, descriptorAtemp
		if ( lineAnorm.gt. 1e-7 ) then


			!COMPUTE HOUSHOLDER VECTOR

			allocate(v(size(linea,dim=1))) !,2 ))
			descriptorV(1) = helper(1) ! dimension used
			descriptorV(2) = helper(3) !start_point

			v = CTYPE0
			v(:descriptorv(1)) = lineA(:descriptorv(1))

			Aii = CTYPE0
			Aii = cbp_getval2d(Atemp, descriptorAtemp, [ii,ii+1], rank, nprocs)

			if ( sum(descriptorv) -1 .ge. ii+1 .and. descriptorv(2) .le. ii+1) then
				v(ii-descriptorv(2) +2 ) = v(ii-descriptorv(2) +2 ) - lineAnorm*Aii/abs(Aii)
			end if

			!vtest(ii+1) = vtest(ii+1) - vtest(ii+1)*sqrt(dot_product(vtest,vtest))/abs(vtest(ii+1))
			!vtest = vtest / sqrt(dot_product(vtest,vtest))

			vnorm = cbp_l2norm(v(:descriptorv(1)),1, descriptorv(1),rank, nprocs, groupproc )

			v = v/vnorm


			!print*, idxrow, idxcol, sum( v(:descriptorv(1)) ) , sum(vtest(descriptorv(2):sum(descriptorv)-1 ))
			v = conjg(v)
			!COMPUTE HOSHOLDER MATRIX

			call cbp_init_local_matrix_descriptor(descriptorH,dimR,ld1,&
							                      rowproc, colproc, rank)
			allocate(H(ld1(1),ld1(2)))

			H = CTYPE0
			call cbp_identity(H,descriptorH)
			H = -0.5*(H)
			call cbp_outer(v, descriptorV, v, descriptorV, H, descriptorH, dimR, 1, &
						rowproc, colproc, rank, nprocs, context  )

			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			H = -2*H

			deallocate(v)

			!Atest = CTYPE0

			! DEBUG : column annihilated?
			!call GetBackMatrix(Atemp, descriptorAtemp, Atest, rowproc, colproc,&
         	 ! 		idxrow, idxcol,rank)



         	  	!Htest = CTYPE0
			!do jj = 1, colA
			!	do kk = 1, colA
			!		Htest(kk,jj) =-2*  vtest(jj) * conjg(vtest(kk))
			!		if (jj .eq. kk) then
			!		Htest(jj,kk) = Htest(jj,kk)+1
			!		end if
			!	end do
			!end do

			!if (rank.eq.0) then
			!	Htest = matmul(Atest,Htest)
			!	print*, sum(Atest-Htest)
			!	print*," "
			!	do jj = 1, colA
			!		print*, Htest(jj,:)
			!	end do
			!	print*, " "
			!end if


			!-------------------------------UPDATE R --------------------------

			allocate(T(size(R,dim=1),size(R,dim=2)))
			T = CTYPE0

			old1 = descriptorR
			old2 = descriptorH
			process_grid = 0

			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)

			call cbp_matmul_descriptors(descriptorR, descriptorH, descriptorT,&
						    dimR, dimR, dimR, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)


			call cbp_redispatch(R, old1, descriptorR, rank, nprocs)
			call cbp_redispatch(H, old2, descriptorH, rank, nprocs)


			call cbp_matmul(R,H,T,descriptorR,descriptorH,descriptorT, &
					rowproc, colproc, rank, process_grid, context)

			R = T
			descriptorR = descriptorT

			deallocate(T)

			!-------------------------------UPDATE RESULT MATRIX A --------------------------

			call init_process_grid(process_grid, rowproc, colproc,rank, idxrow, idxcol)
			old1 = 0
			old2 = 0
			A = Atemp
			descriptorA = descriptorAtemp

			old1 = descriptorA
			old2 = descriptorH

			call cbp_matmul_descriptors(descriptorA, descriptorH, descriptorAtemp,&
						    dimA, dimR, dimA, &
						    ld1, ld2, ld3 , rowproc, &
						    colproc, process_grid,rank)


			call cbp_redispatch(A,old1,descriptorA, rank, nprocs)
			call cbp_redispatch(H,old2,descriptorH, rank, nprocs)


			call cbp_matmul(A,H,Atemp,descriptorA,descriptorH,descriptorAtemp, &
					rowproc, colproc, rank, process_grid, context)

			call MPI_BARRIER(MPI_COMM_WORLD, ierr)

			!Atest = CTYPE0

			! DEBUG : column annihilated?
			!call GetBackMatrix(Atemp, descriptorAtemp, Atest, rowproc, colproc,&
         	  	!	idxrow, idxcol,rank)

			!if (rank.eq.0) then
			!print*, "Annihilated col?  ", ii,  sum(Atest(ii,ii+2:colA))
			!end if

			deallocate(H)
		!+++++++++++++++++++++++++++++++++++++++++
		end if
		deallocate(lineA)
		deallocate(groupproc)
	    end if
	    !+++++++++++++++++++++++++++++++++++++++++



	    end do

	A = Atemp
	descriptorA = descriptorAtemp

	deallocate(Atemp)
	deallocate(process_grid)
	!deallocate(Atest)
	!deallocate(Htest)
end subroutine

!===================================================================================
!===================================================================================

subroutine check_descmat(matrix, desc, funct)
	CTYPE , dimension(:,:) :: Matrix
	integer, dimension(4) :: desc
	character(*) :: funct

	if(desc(1) .gt. size(Matrix,dim=1) )then
		print*, "Error, matrix too small with respect to descriptor", funct
	end if
	if(desc(2) .gt. size(Matrix,dim=2) )then
		print*, "Error, matrix too small with respect to descriptor", funct

	end if
end subroutine
!===================================================================================
!===================================================================================

subroutine cbp_identity(I, descriptorI, add)
	CTYPE , dimension(:,:) :: I
	integer , dimension(4) :: descriptorI
	integer :: ii, localrow, localcol, overlap
	integer , optional :: add

	if ( all( descriptorI .ne. 0 )) then
	if(.not. present(add)) then
		I = CTYPE0
	else if ( add .eq. 0)  then
		I = CTYPE0
	end if

	call check_descmat(I, descriptorI, "cbp_identity")

	overlap =  min(sum(descriptorI(1:3:2)), sum(descriptorI(2:4:2)))-1 - &
		   max( descriptorI(3) , descriptorI(4) )

	if ( overlap .ge. 0 ) then

	 do ii = 1,overlap+1

	   	localrow = max(descriptorI(3),descriptorI(4))-descriptorI(3) +ii
	   	localcol = max(descriptorI(3),descriptorI(4))-descriptorI(4) +ii

	   	if(localrow .gt. descriptorI(1)) then
	   		print*, "Error, bug inside cbp_identity"
	   	end if
	   	if(localcol.gt. descriptorI(2)) then
	   		print*, "Error, bug inside cbp_identity"
	   	end if


		I (localrow, localcol ) = I (localrow, localcol ) + CTYPE1

	  end do
	end if
	end if



end subroutine

!===================================================================================
!===================================================================================
function cbp_l2norm(v, from, to , rank, nprocs, procs) result(norm)

	integer, optional , dimension(:), intent(in) :: procs
	CTYPE, dimension(:) :: v
	CTYPE, dimension(2) :: descriptorv
	double precision :: halfnorm2 , temp, norm

	integer :: ierr, rank, nprocs, from ,to ,ii

	halfnorm2 = dot_product(v(from:to),v(from:to))

	norm = 0.0
	!print*, "halfnorm2", halfnorm2
	if (to .lt. from) then
		halfnorm2 = 0.0
	end if


	if (.not. present(procs)) then
		do ii = 0, nprocs-1
			if (rank .eq. ii .and. ii .eq. 0 ) then
				norm = norm +halfnorm2
			else if ( rank .eq. ii) then
				call MPI_SEND(halfnorm2,1, MPI_DOUBLE_PRECISION, 0 ,ii, MPI_COMM_WORLD,ierr )
				!norm = norm + halfnorm2
			else if ( rank .eq. 0) then
				call MPI_RECV(temp,1, MPI_DOUBLE_PRECISION, ii , ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

				norm = norm + temp
			end if
		end do
	else

		do ii = 1, size(procs)

			if ( rank .eq. procs(ii) .and. rank .eq. 0 ) then
				norm = norm +halfnorm2

			else if ( rank .eq. procs(ii)) then

				call MPI_SEND(halfnorm2,1, MPI_DOUBLE_PRECISION,0, procs(ii),MPI_COMM_WORLD,ierr )
				!norm = norm + halfnorm2
			else if ( rank.eq.0) then

				call MPI_RECV(temp,1, MPI_DOUBLE_PRECISION, procs(ii) ,procs(ii), MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				norm = norm + temp


			end if
		end do
	end if


	if ( rank.eq.0) then
			norm = sqrt(norm)

			call MPI_BCAST(norm,1, MPI_DOUBLE_PRECISION, 0 ,MPI_COMM_WORLD,ierr )

	else
			call MPI_BCAST(norm,1, MPI_DOUBLE_PRECISION, 0 ,MPI_COMM_WORLD,ierr )

	end if


end function


function cbp_sum(v, from, to , rank, nprocs, procs) result(norm)

	integer, optional , dimension(:), intent(in) :: procs
	CTYPE, dimension(:) :: v
	CTYPE, dimension(2) :: descriptorv
	CTYPE :: halfnorm2 , temp, norm

	integer :: ierr, rank, nprocs, from ,to ,ii

	halfnorm2 = sum(v(from:to))

	norm = 0.0
	!print*, "halfnorm2", halfnorm2
	if (to .lt. from) then
		halfnorm2 = 0.0
	end if


	if (.not. present(procs)) then
		do ii = 0, nprocs-1
			if (rank .eq. ii .and. ii .eq. 0 ) then
				norm = norm +halfnorm2
			else if ( rank .eq. ii) then
				call MPI_SEND(halfnorm2,1, MPI_CTYPE, 0 ,ii, MPI_COMM_WORLD,ierr )
				!norm = norm + halfnorm2
			else if ( rank .eq. 0) then
				call MPI_RECV(temp,1, MPI_CTYPE, ii , ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )

				norm = norm + temp
			end if
		end do
	else

		do ii = 1, size(procs)

			if ( rank .eq. procs(ii) .and. rank .eq. 0 ) then
				norm = norm +halfnorm2

			else if ( rank .eq. procs(ii)) then

				call MPI_SEND(halfnorm2,1, MPI_CTYPE,0, procs(ii),MPI_COMM_WORLD,ierr )
				!norm = norm + halfnorm2
			else if ( rank.eq.0) then

				call MPI_RECV(temp,1, MPI_CTYPE, procs(ii) ,procs(ii), MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
				norm = norm + temp


			end if
		end do
	end if


	if ( rank.eq.0) then

			call MPI_BCAST(norm,1, MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr )

	else
			call MPI_BCAST(norm,1, MPI_CTYPE, 0 ,MPI_COMM_WORLD,ierr )

	end if


end function

!FROM HERE MUST BE TESTED
!===================================================================================
!===================================================================================

subroutine cbp_outer(u, descriptorU, v, descriptorV, R, descriptorR, dimR, add, rowproc, colproc, rank, nprocs, context  )

	CTYPE, dimension(:) :: v, u
	CTYPE, dimension(:), allocatable :: urow, vcol
	CTYPE, dimension(:,:) :: R
	integer, dimension(4) :: descriptorR, temp
	integer, dimension(2) :: descriptorv,descriptoru, dimR, ldim

	integer :: idxrow, idxcol, rowproc, colproc, rank, nprocs, dimv, add , dimu
	type(c_ptr) :: context
	dimu = dimR(1)
	dimv = dimR(2)

	call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
	call cbp_init_local_matrix_descriptor(temp, [dimu, dimv], ldim, rowproc, colproc , rank )
	!print*,"outer", max( ldim(1) , descriptoru(1)), max( ldim(2) , descriptorv(1))
	allocate(urow(max( ldim(1) , descriptoru(1))))
	allocate(vcol(max( ldim(2) , descriptorv(1))))

	urow = CTYPE0
	vcol = CTYPE0



	if (descriptoru(1) .gt. 0) then
	urow(:descriptoru(1)) = u(:descriptoru(1))
	end if
	if (descriptorv(1) .gt. 0) then
#if defined(CBPCMPLX) | defined(CBPZMPLX)
	vcol(:descriptorv(1)) = conjg(v(:descriptorv(1)))
#endif
#if defined(CBPREAL) | defined(CBPDOUBLE)
	vcol(:descriptorv(1)) = v(:descriptorv(1))
#endif
	end if


	call cbp_redispatch_1d(vcol, descriptorv, temp(2:4:2), rank, nprocs)
	call cbp_redispatch_1d(urow, descriptoru, temp(1:3:2), rank, nprocs)


	if (add.eq.1) then
		call cbp_redispatch(R,descriptorR,&
				temp,rank,nprocs)
	end if
	descriptorR = temp

	call cbpfor_outer(urow,vcol,R,add,0,0.0,context,temp(1), temp(2), descriptorR(:2))

	deallocate(vcol)
	deallocate(urow)
end subroutine

!===================================================================================
!===================================================================================

subroutine cbp_redispatch_1d(v,olddesc,newdesc, rank, nprocs, ppp)
	integer, dimension(2) :: olddesc, newdesc, fin, ini
	integer, dimension(2)  ::  ends , starts !) 1 ini 2 fin
	CTYPE , dimension(:) :: v
	CTYPE , dimension(:), allocatable :: vtemp
	integer, optional :: ppp
	integer :: ppp1
	integer :: newdim, ii, jj,  localstart, localend, rank, nprocs, ierr, RS, RR

	if (.not. present(ppp) ) then
		ppp1 = 0
	else
		ppp1 = 1
	end if


	allocate(vtemp(size(v)))
	if ( any(olddesc .eq. 0)) then
		olddesc = 0
	end if
	if ( any(newdesc .eq. 0)) then
		newdesc = 0
	end if
	vtemp = CTYPE0
	do ii = 0, nprocs-1
		do jj = 0, nprocs-1

			if (rank.eq.ii .or. rank.eq.jj) then

				!1) check overlap
				if (rank.eq. ii .and. rank.eq.jj) then
					ini = olddesc
					fin = newdesc
				else if (rank.eq.ii) then
					ini = olddesc
					call MPI_SEND( olddesc, 2, MPI_INTEGER, jj, ii, MPI_COMM_WORLD, ierr )
					call MPI_RECV( fin, 2, MPI_INTEGER, jj, jj, MPI_COMM_WORLD,  MPI_STATUS_IGNORE, ierr )
				else if(rank.eq.jj) then
					fin = newdesc
					call MPI_RECV( ini, 2, MPI_INTEGER, ii, ii, MPI_COMM_WORLD,  MPI_STATUS_IGNORE, ierr )
					call MPI_SEND( newdesc, 2, MPI_INTEGER, ii, jj, MPI_COMM_WORLD, ierr )
				end if
				ends = [sum(ini)-1, sum(fin)-1]
				starts = [ini(2) , fin(2) ]
				!if (ppp1.eq.1) then
				!print*, ii, jj , minval(ends)-maxval(starts) , ini, fin
				!end if
				!2) pass data
				if ( minval(ends)-maxval(starts) .ge. 0 ) then

					if (rank.eq.ii .and. rank.eq.jj) then
						vtemp(maxval(starts) -  newdesc(2) +1:minval(ends) -  newdesc(2) +1) =  &
						    v(maxval(starts) -  olddesc(2) +1:minval(ends) -  olddesc(2) +1)
					else if (rank.eq.ii) then
						localstart = maxval(starts) - olddesc(2) +1
						localend = minval(ends) - olddesc(2) +1

						call MPI_SEND(v(localstart:localend), &
								minval(ends)-maxval(starts)+1, MPI_CTYPE, jj, ii,MPI_COMM_WORLD,ierr)

					else if (rank .eq. jj ) then
						localstart = maxval(starts) -  newdesc(2) +1
						localend = minval(ends) -  newdesc(2) +1
						call MPI_RECV(vtemp(localstart:localend) ,&
								minval(ends)-maxval(starts)+1, MPI_CTYPE, ii, ii,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )



					end if
					!if(rank.eq.ii) then
					!	call  MPI_Wait(RS, MPI_STATUS_IGNORE, IERR)
					!end if
					!if(rank.eq.jj) then
					!	call  MPI_Wait(RR, MPI_STATUS_IGNORE, IERR)
					!end if

				end if
			end if


		end do
	end do

	v = vtemp
	deallocate(vtemp)

end subroutine

!===================================================================================
!===================================================================================
!SVD Decomposition
!A must be dispathced in blocks " basic descriptor or matmul descriptors are ok"
!Also out from transpose or adjoint
#define CUDACOMP
#ifdef CUDACOMP
subroutine cbp_SVD_GulobReinsch(A, descriptorA, dimA, d,descriptorD ,rowproc, colproc, rank, nprocs)
	!use svd
	use kernels
  	use cudafor
  	use cusparse
  	use iso_c_binding

	CTYPE , dimension(:,:) :: A, U0,V0
	CTYPE , dimension(:) :: d
	CTYPE,  dimension(:) , allocatable ::  u, uu, dd, evec, testevect

	CTYPE, dimension(2) ::  ker_condition

	integer, dimension(4) :: descriptorA, descriptorU0, descriptorV0, svddesc, tmpdesc
	CTYPE :: tempu, tempd, tempd2, tempuu

	CTYPE,dimension(:),allocatable, device :: d_d, u_d, dd_d, uu_d
	CTYPE :: temp, trace
	integer :: ii  , it , maxdim, localpivot, localmaxit , localoddeven,jj,ll

	integer :: rank, nprocs, ierr, ngpu, cntr, ldiag, rowproc, colproc

	integer, dimension(2) ::  dimA, descriptord, descriptoru, descriptordd, descriptorevec
	integer :: side, oddeven, maxit,pivot,workonpivot, tag_factor

	double precision :: e, sens
	type(dim3) :: grid, tBlock
	logical :: PerformStep
	integer :: ownspivot, kk

	!Update of U and V Matrices
	CTYPE , dimension(:,:), allocatable :: M, Tstep, Ttemp, ortho
	integer , dimension(:,:) , allocatable :: process_grid
   	integer, dimension(2) :: ld1, ld2, ld3 , dd_end
   	integer, dimension(4) :: temp_descriptor, step_descriptor, descriptorTstep, descriptorTtemp, &
		 		 oldtemp, oldstep, old_in
	integer :: idxrow, idxcol
	integer :: buffersize

   	type(c_ptr) :: context
   	type(cusparseHandle) :: cusparse_context


  	tag_factor = 9
	e = 1e-10
	sens = 1e-50
	cntr = 0


	call cbp_getcolrowidx(idxrow, idxcol, rank, rowproc, colproc)
	descriptord = descriptorA(1:3:2)
	descriptoru = descriptorA(1:3:2)

	allocate(process_grid(rowproc, colproc) )

	! A is a upper bidiagonal matrix, redispatch it as blocks overlapping on the diagonal
	! with a 2x2 overlap on the diagonal
	!---------------------------------------------------------


	ldiag = ceiling(real(min(dimA(1), dimA(2)))/ real(nprocs))
	svddesc = [ ldiag+2 , ldiag+2 , ldiag*rank+1, ldiag*rank+1 ]

	if (rank .eq. nprocs-1)  then
		svddesc(1:2) = minval(dimA) - ldiag*rank
		if ( dimA(2) .gt. dimA(1) ) then
			svddesc(2) = svddesc(2) +1 ! last bidiagonal element
		end if
	end if

	maxdim = maxval(svddesc(1:2))

	tBlock = dim3( 2,1,1)
	grid = dim3(20,1,1)

	!allocate(d(max(maxval(descriptorA(1:2)),maxdim+1)))
	if ( size(d) .lt. max(maxval(descriptorA(1:2)),maxdim+1)) then
		print*, "size singularvalues array is unsafe, &
			suggested max(maxval(descriptorA(1:2)),ceiling(real(min(dimA(1), dimA(2)))/ real(nprocs)) +4)"
	end if
	allocate(u(max(maxval(descriptorA(1:2)),maxdim+1)))
	allocate(uu(max(maxval(descriptorA(1:2)),maxdim+1)))
	allocate(dd(max(maxval(descriptorA(1:2)),maxdim+1)))
	allocate(evec(max(maxval(descriptorA(1:2)),maxdim+1)))
	!allocate(testevect(dimA(2)))
	!allocate(M(dimA(1), dimA(2)))
	!allocate(Ortho(dimA(1), dimA(1)))

	!device pointers
	allocate(d_d(maxdim+1))
	allocate(dd_d(maxdim+1))
	allocate(uu_d(maxdim+1))
	allocate(u_d(maxdim+1))



	allocate(Tstep(size(A,dim=1),size(A,dim=2)))
	allocate(Ttemp(size(A,dim=1),size(A,dim=2)))

	call cbp_setget_diagonal(A,descriptorA, d, descriptord, 0, 'get')
	call cbp_setget_diagonal(A,[descriptorA(1),descriptorA(2),descriptorA(3),descriptorA(4)-1],&
				 u, descriptoru, 0, 'get')



	call cbp_redispatch_1d(d, descriptord, [svddesc(1), svddesc(3)],rank,nprocs)
	call cbp_redispatch_1d(u, descriptoru, [svddesc(1), svddesc(3)],rank,nprocs)



	u(maxdim:) = CTYPE0
	d_d(:maxdim) = d(:maxdim)
	u_d(:maxdim) = u(:maxdim)
	uu_d = CTYPE0
	dd_d = CTYPE0


	pivot = 1




	it = 0

	if (pivot .ge. svddesc(3) .and. pivot .le. svddesc(3)+maxdim -1 ) then
	  	  localpivot = pivot-svddesc(3)+1
	else
		localpivot =0
	end if
	PerformStep = .TRUE.
	ownspivot = 0


	!------------------------------------------------------
	do while(ownspivot .ne. nprocs)


		  maxit = max( (it+2)/2  , 0 )  ! - svddesc(3) works only if the upper left corner of the block is on the global diagonal
		  oddeven =  1-mod(it,4)/2 !first it is odd, second is odd, third even , 4th even , 5 odd.. Global ticking
		  side = mod(it,2) !global ticking
		  it = it+1

		  localoddeven =  (1-oddeven)*(1-mod(svddesc(3),2))+mod(svddesc(3),2)*oddeven  ! 0,0 -> 1 / 1,1->1 /  1,0 -> 0 / 0,1-> 1
		  localmaxit = max( maxit - svddesc(3) +1 , 0)

		  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

		  if (PerformStep) then

			  call GulobReinsch_step<<<tBlock,grid>>>(   d_d, u_d, &
		  						     dd_d, uu_d,    & !diagonal, abovediagonal,cosine, since
								     e, sens,      & !precision, sensibility to treat sine and cosine
								     localpivot,   &
								     localmaxit,   &
								     side,         & !multiply left or right, odd indeces in parallel or even indices inparallel
								     localoddeven, &
								     maxdim,&
								     max(0, -nprocs+2+rank))	   !maximum dimension
								     ! lastproc


			  !to check
			  if ( localpivot .gt. 0 ) then
			  	u(localpivot) = u_d(localpivot)
			  	uu(localpivot) = uu_d(localpivot)
			  	dd(localpivot+1) = dd_d(localpivot+1)
			  end if

			  !GPU to CPU

			u(maxdim-1:maxdim)    = u_d(maxdim-1:maxdim)
		  	d(maxdim-1:maxdim)    = d_d(maxdim-1:maxdim)
		 	dd(maxdim-1:maxdim+1) = dd_d(maxdim-1:maxdim+1)
		  	u(:2) =  u_d(:2)
		  	d(:2) = d_d(:2)
		  	dd(:2) = dd_d(:2)

			! print"(I3, F8.5,SP,F0.0,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i')",&
			! 	 rank ,u(1), d(1) , u(maxdim-1), d(maxdim-1)
			 !GPU to GPU
			! print*, " "
			end if
			 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

			  if(side .eq. 0 .and. rank.ge.ownspivot) then
			  	if (rank.lt.nprocs-1 .and. performStep) then

				  	call MPI_SEND( d(maxdim-1),1, MPI_CTYPE, rank+1, tag_factor*rank,   MPI_COMM_WORLD, ierr )
				  	call MPI_SEND( dd(maxdim) ,1, MPI_CTYPE, rank+1, tag_factor*rank+1, MPI_COMM_WORLD, ierr )

				  	if ( localoddeven .eq. mod(maxdim,2)  ) then !Half overlap in the two blocks some to send some to receive
				  	call MPI_RECV( u(maxdim-1) ,1, MPI_CTYPE, rank+1, tag_factor*rank+2, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
				  	call MPI_RECV( d(maxdim)   ,1, MPI_CTYPE, rank+1, tag_factor*rank+3, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
				  	else !all to send
				  	call MPI_SEND( u(maxdim-1),1, MPI_CTYPE, rank+1, tag_factor*rank+2, MPI_COMM_WORLD, ierr )
				  	call MPI_SEND( d(maxdim)  ,1, MPI_CTYPE, rank+1, tag_factor*rank+3,   MPI_COMM_WORLD, ierr )

				  	end if
			  	end if
			  	if (rank.gt.ownspivot ) then

				  	call MPI_RECV( d(1)   ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1),   MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
				  	call MPI_RECV( dd(2)  ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
				  	if (localoddeven .eq. 0) then
				  	call MPI_SEND( u(1)   ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+2,   MPI_COMM_WORLD,ierr )
				  	call MPI_SEND( d(2)   ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+3,   MPI_COMM_WORLD,ierr )
				  	else
				  	call MPI_RECV( u(1)   ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
				  	call MPI_RECV( d(2)   ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+3,   MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
				  	end if
			  	end if

			  else if(side .eq. 1) then
			  	if (rank.gt.ownspivot) then
			  		if (localoddeven .eq. 1) then !send all the shared square
				  	call MPI_SEND( d(1),  1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+4, MPI_COMM_WORLD, ierr )
				  	call MPI_SEND( u(1)  ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+5, MPI_COMM_WORLD, ierr )
				  	call MPI_SEND( uu(1) ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+6, MPI_COMM_WORLD, ierr )
				  	call MPI_SEND( dd(2) ,1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+7, MPI_COMM_WORLD, ierr )
				  	end if
				  	call MPI_SEND( d(2),1, MPI_CTYPE, rank-1, tag_factor*(rank-1)+8, MPI_COMM_WORLD, ierr )

			  	end if
			  	if (rank.lt.nprocs-1.and. performStep ) then

			  		if ( mod(maxdim-1,2) .eq. localoddeven) then !send all the shared square
				  	call MPI_RECV( d(maxdim-1)       ,1, MPI_CTYPE, rank+1, tag_factor*(rank)+4, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
				  	call MPI_RECV( u(maxdim-1)       ,1, MPI_CTYPE, rank+1, tag_factor*(rank)+5, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )

				  	call MPI_RECV( uu(maxdim-1)      ,1, MPI_CTYPE, rank+1, tag_factor*(rank)+6, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
				  	call MPI_RECV( dd(maxdim)        ,1, MPI_CTYPE, rank+1, tag_factor*(rank)+7, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
				  	end if
				  	call MPI_RECV( d(maxdim)         ,1, MPI_CTYPE, rank+1, tag_factor*(rank)+8, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )


				  	!print*, abs(tempd2 -  d(maxdim)), abs(tempd  -  d(maxdim-1)), abs(tempuu -  dd(maxdim) )
			  	end if
			  end if


		 !endif !PerformStep

		 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		! print"(I3, F8.5,SP,F0.0,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i',F8.5,SP,F8.5,'i')" ,&
		 !		rank , u(1), d(1) ,u(maxdim-1), d(maxdim-1)

		 !JUDGE PIVOT AFTER EXCHANGE
		 !Is my first element the pivot?

		 if (pivot .ge. svddesc(3) .and. pivot .le. svddesc(3)+maxdim-1 ) then
		  	  localpivot = pivot-svddesc(3)+1



			  if (abs(u(localpivot))+abs(uu(localpivot))+abs(dd(localpivot+1)) .lt. 3*e) then



				if (localpivot .eq. maxdim-1 + max(0,rank-nprocs+2))then
					ownspivot = rank+1
					localpivot = 0
					pivot = 0
					do ll =0, nprocs-1
			  		call MPI_SEND( rank+1, 1, MPI_INTEGER, ll, ll, MPI_COMM_WORLD, ierr )
			  		end do

				else
					do ll =0, nprocs-1
			  		call MPI_SEND( rank, 1, MPI_INTEGER, ll, ll, MPI_COMM_WORLD, ierr )
			  		end do
					pivot = pivot+1
					localpivot = localpivot+1

				end if


			  else
			  	do ll =0, nprocs-1
			  	call MPI_SEND( rank, 1, MPI_INTEGER, ll, ll, MPI_COMM_WORLD, ierr )
			  	end do

			  end if
		  else
		  	 call MPI_RECV(ownspivot, 1, MPI_INTEGER, ownspivot,rank ,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )
		  	 if(rank.eq.ownspivot) then
		  		localpivot = 2
		  		pivot = svddesc(3)+1
		  	 end if

		  end if
		  PerformStep = rank .ge. ownspivot


		  !CPU to GPU

		  u_d(maxdim-1:maxdim)    =  u(maxdim-1:maxdim)
		  d_d(maxdim-1:maxdim)    = d(maxdim-1:maxdim)
		  dd_d(maxdim-1:maxdim+1) = dd(maxdim-1:maxdim+1)
		  u_d(:2) =  u(:2)
		  d_d(:2) = d(:2)
		  dd_d(:2) = dd(:2)
		  !print*, "end H2D"

		  descriptorD = [min(ldiag,minval(dimA) - rank*ldiag), svddesc(3) ]


		  !-------------------------------------------------------------------------

  end do
  !DEBUG SINGULAR VALUES OK
	d(:maxdim) = d_d(:maxdim)
	deallocate(u)
	deallocate(uu)
	deallocate(dd)
	deallocate(evec)
	!deallocate(testevect(dimA(2)))
	!deallocate(M(dimA(1), dimA(2)))
	!deallocate(Ortho(dimA(1), dimA(1)))

	!device pointers
	deallocate(d_d)
	deallocate(dd_d)
	deallocate(uu_d)
	deallocate(u_d)





  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  deallocate (process_grid)
end subroutine
#endif

function judge_pivot_svd(localpivot, coldim, rank, nprocs) result(workonpivot)
	integer :: rank, nprocs, pivot, coldim, localpivot
	integer, dimension(4) :: svddesc
	integer :: workonpivot

	!I can judge the pivot if it is the first only if i am the process 1
	!otherwise it is judge by the one before

	if (rank .gt. 0 .and. rank .lt. nprocs) then
		if (localpivot .ge. 2 .and. localpivot .le. coldim-1) then
			workonpivot = 1
		else
			workonpivot = 0
		end if

	else if (rank.eq.0) then
		if (localpivot .ge. 1 .and. localpivot .le. coldim-1) then
			workonpivot = 1
		else
			workonpivot = 0
		end if
	end if
	return

end function

subroutine cbp_setget_diagonal(A,descriptorA,d,descriptord,add,op)
	CTYPE , dimension(:,:) :: A
	CTYPE , dimension(:) :: d
	integer , dimension(4) :: descriptorA
	integer , dimension(2) :: descriptord
	integer :: ii, localrow, localcol, overlap, rank, nprocs
	integer :: add
	character*3 :: op


	if (op .eq. 'set') then
	if ( add .eq. 0)  then
		A = CTYPE0
	end if
	end if
	if (op .eq. 'get') then
	if ( add .eq. 0)  then
		d = CTYPE0
	end if
	end if

	call check_descmat(A, descriptorA, "cbp_setget_diagonal")

	overlap =  min(sum(descriptorA(1:3:2)), sum(descriptorA(2:4:2)))-1 - &
		   max( descriptorA(3) , descriptorA(4) )
	if (op .eq. 'get' .and. add .eq. 0 ) then
		descriptord = [max(0,overlap+1), max( descriptorA(3) , descriptorA(4) ) ]
	else
		call cbp_redispatch_1d(d, descriptord, descriptorA(1:3:2),rank,nprocs)
		descriptord = [max(0,overlap+1), max( descriptorA(3) , descriptorA(4) ) ]
	end if



	if ( overlap .ge. 0 ) then

	 do ii = 1,overlap+1
	   	localrow = max(descriptorA(3),descriptorA(4))-descriptorA(3) +ii
	   	localcol = max(descriptorA(3),descriptorA(4))-descriptorA(4) +ii

	   	if(localrow .gt. descriptorA(1)) then
	   		print*, "Error, bug inside cbp_setget_diagonal"
	   	end if
	   	if(localcol.gt. descriptorA(2)) then
	   		print*, "Error, bug inside cbp_setget_diagonal"
	   	end if
	   	if (op .eq.'get') then
	   		!print*, localrow, localcol, size(d) , ii , descriptord(1)
			d(ii) = d(ii) + A (localrow, localcol )
		else if (op .eq. 'set') then
			A (localrow, localcol ) = d(ii) + A (localrow, localcol )
		end if
	  end do
	end if




end subroutine

function cbp_singular_value_array_size(dimA, descriptorA, rank, nprocs) result(thesize)
	integer , dimension(4)  :: svddesc, descriptorA
	integer :: ldiag, maxdim , thesize ,rank, nprocs
	integer , dimension(2) :: dimA


	ldiag = ceiling(real(min(dimA(1), dimA(2)))/ real(nprocs))
	svddesc = [ ldiag+2 , ldiag+2 , ldiag*rank+1, ldiag*rank+1 ]

	if (rank .eq. nprocs-1)  then
		svddesc(1:2) = minval(dimA) - ldiag*rank
		if ( dimA(2) .gt. dimA(1) ) then
			svddesc(2) = svddesc(2) +1 ! last bidiagonal element
		end if
	end if

	maxdim = maxval(svddesc(1:2))

	thesize = max(maxval(descriptorA(1:2)),maxdim+1)

	return

end function

!Got l lowerdiagonal ( firs element zero)
!U upperdiagonal
!D diagonal

!cndition are the two x components at 0, and 1 , so

!to get x2 = -l(1) x0 - d(1) x1

!Given rowdim equation and 2 condition i get rowdim+1 solutions
subroutine svd_tridiagonal_ker(l,d,u,sol,rowdim,condition)
	CTYPE, dimension(:) :: l,u, d , sol
	CTYPE, dimension(2) :: condition
	CTYPE, dimension(:), allocatable :: temp
	integer :: rowdim , ii, rank, ierr

	!call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	if (rowdim .gt. 0 ) then
	allocate (temp(rowdim+1))
	temp(1:2) = condition
			  !temp(ii+2) represent the ii+1 element
	do ii = 1, rowdim !! ii is the row
		temp(ii+2) =( -l(ii)*temp(ii)-d(ii)*temp(ii+1) ) / u(ii)
		!print*, temp(ii+1)
		!print*,ii, '|||',abs(temp(ii)*l(ii)+temp(ii+1)*d(ii)+temp(ii+2)*u(ii))
		!print*,ii, " ooo "	,temp(ii+2), temp(ii+1) , temp(ii) , l(ii), d(ii), u(ii)
		!print*,rank, '|||', temp(ii)*l(ii), temp(ii+1)*d(ii), temp(ii+2)*u(ii), '|'
	end do
	sol(:rowdim+1) = temp(2:rowdim+2)
	!print*, "  "
	!debug check solution
	!SOL(0) IS X0 AND SO ON
	do ii = 1, rowdim+1 !! ii is the row
		!print*,"         " , ii, '|', sol(ii)
		!print*, 'heyy2', sol(ii-1)*l(ii) + sol(ii)*d(ii) + sol(ii+1)*u(ii)
	end do

	deallocate(temp)
	end if
end subroutine


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!1D descriptor generator

subroutine cbp_init_1D_descriptor(descriptor, full_dim, rank, nprocs)
  integer, dimension(2) :: descriptor
  integer :: rank, nprocs, full_dim, dstr_dim

  dstr_dim = ceiling( real(full_dim) / real(nprocs) )
  descriptor(1) = dstr_dim
  descriptor(2) = rank*dstr_dim +1
  if ((rank+1)*dstr_dim .gt. full_dim) then
    descriptor(1) = full_dim - rank*dstr_dim
    if (descriptor(1) .le. 0) descriptor = 0
  end if

end subroutine

end module
