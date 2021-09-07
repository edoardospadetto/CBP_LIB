#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define CTYPE1 cmplx(1.0,0.0)
#define MPI_CTYPE MPI_COMPLEX
#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define CTYPE0 0.0
#define CTYPE1 1.0
#define MPI_CTYPE MPI_DOUBLE_PRECISION
#endif

#ifdef CBPZMPLX 
#define CTYPE0 dcmplx(0.0,0.0)
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#endif

#ifdef CBPREAL 
#define CTYPE0 0.0
#define CTYPE1 1.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#endif


module cbp_DistributeMODCUDA
    use cbp_distributeMOD
    use mpi
    use cbpfor_parallel
    use misc_mod
        
    contains 
    
    
subroutine cbp_SVD_GulobReinsch(U0, descriptoru0, A, descriptorA, dimA, V0, descriptorV0, ngpu, rank, nprocs) 
	!use svd
	use kernels
  	use cudafor
  	
	CTYPE , dimension(:,:) :: A, U0,V0
	CTYPE,  dimension(:) , allocatable :: d, u, uu, dd, cosj , sinj
	CTYPE , dimension(:,:), allocatable ::M
	integer, dimension(4) :: descriptorA, descriptorU0, descriptorV0, svddesc
	CTYPE,dimension(:),allocatable, device :: d_d, u_d, s_d, c_d, dd_d, uu_d
	integer :: ii  , it , maxdim, localpivot, localmaxit , localoddeven,jj
	integer :: rank, nprocs, ierr, ngpu, cntr, ldiag
	integer, dimension(2) ::  dimA, descriptord, descriptoru
	integer :: side, oddeven, maxit,pivot,workonpivot
	double precision :: e, sens 
	type(dim3) :: grid, tBlock
	logical :: skip
   

	tBlock = dim3(1,1,1)
	grid = dim3(4,1,1)
  
	e = 1e-11
	sens = 1e-11
	cntr = 0 
	
	allocate(d(maxval(descriptorA(1:2))))
	allocate(u(maxval(descriptorA(1:2))))
	allocate(uu(maxval(descriptorA(1:2))))
	allocate(dd(maxval(descriptorA(1:2))))
	
	descriptord = descriptorA(1:3:2)
	descriptoru = descriptorA(1:3:2)
	
	
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
	
	!device pointers
	allocate(d_d(maxdim+1))
	allocate(dd_d(maxdim+1))
	allocate(uu_d(maxdim+1))
	allocate(u_d(maxdim+1))
	allocate(c_d(maxdim+1))
	allocate(s_d(maxdim+1))
	allocate(M(maxdim,maxdim))
	call cbp_setget_diagonal(A,descriptorA, d, descriptord, 0, 'get')
	call cbp_setget_diagonal(A,[descriptorA(1),descriptorA(2),descriptorA(1),descriptorA(1)-1],& 
				 u, descriptoru, 0, 'get')
	
	call cbp_redispatch_1d(d, descriptord, [maxdim, svddesc(3)],rank,nprocs)
	call cbp_redispatch_1d(u, descriptoru, [maxdim, svddesc(3)],rank,nprocs) 

	d_d(:maxdim) = d(:maxdim)
	u_d(:maxdim) = u(:maxdim) 
	uu_d = CTYPE0
	dd_d = CTYPE0
	
	pivot = 1
	print*, svddesc
	!------------------------------------------------------
	
	do while(.True.)
		
		  maxit = max( (it+2)/2  , 0 )  ! - svddesc(3) works only if the upper left corner of the block is on the global diagonal
		  oddeven =  1-mod(it,4)/2 !first it is odd, second is odd, third even , 4th even , 5 odd.. Global ticking
		  side = mod(it,2) !global ticking
		  it = it+1
		  
		  localoddeven =  (1-oddeven)*(1-mod(svddesc(3),2))+mod(svddesc(3),2)*oddeven  ! 0,0 -> 1 / 1,1->1 /  1,0 -> 0 / 0,1-> 1    
		  localmaxit = max( maxit - svddesc(3) +1 , 0) 
		  
		  if (pivot .ge. svddesc(3) .and. pivot .le. svddesc(3)+maxit -1 ) then 
		  	  
		  	  localpivot = pivot-svddesc(3)+1
		  	  if (judge_pivot_svd(localpivot, maxdim, rank, nprocs) .eq. 1)then 
			  if (  (localpivot .lt. maxdim-1 .and. abs(u(pivot))+abs(uu(pivot))+abs(dd(pivot+1)) .lt. 3*e) .or. & 
			  	(localpivot .eq. maxdim-1 .and. abs(u(pivot))+abs(dd(pivot+1)) .lt. 2*e ))then  
				
				pivot = pivot+1
				localpivot = pivot-svddesc(3)+1
					
				if (localpivot .eq. svddesc(2)) then 
				skip = .TRUE.
				end if
			  end if 
			  end if 
		  else 
		  	localpivot = 0
		  end if 
		  if (.not. skip) then 	 
	          call GulobReinsch_step<<<tBlock,grid>>>(d_d, u_d, dd_d, uu_d, c_d, s_d, & !diagonal, abovediagonal,cosine, since
								     e, sens,&      !precision, sensibility to treat sine and cosine
								     localpivot,& 
								     localmaxit, &  
								     side,&         !multiply left or right, odd indeces in parallel or even indices in parallel
								     localoddeven, &       
								     maxdim, &	   !maximum dimension
				                                    max(0,rank-nprocs+2))! lastproc
		  
		 end if 
		 
		  !to check 
		  if ( localpivot .ne. 0) then 
		  	u(localpivot) = u_d(localpivot) 
		  	uu(localpivot) = uu_d(localpivot) 
		  	dd(localpivot+1) = dd_d(localpivot+1) 
		  end if 
		  
		  !GPU to CPU
		  u(:maxdim) =  u_d(:maxdim)
		  d(:maxdim) = d_d(:maxdim)
		  dd(:maxdim) = dd_d(:maxdim)
		   
		  !GPU to GPU 
		  if(side .eq. 0) then 
		  	if (rank.lt.nprocs-1) then 
		  	call MPI_SEND( d(maxdim-1:maxdim), 2, MPI_CTYPE, rank+1, 2*rank, MPI_COMM_WORLD, ierr )
		  	call MPI_SEND( u(maxdim-1) ,1, MPI_CTYPE, rank+1, 2*rank, MPI_COMM_WORLD, ierr )
		  	call MPI_SEND( dd(maxdim) ,1, MPI_CTYPE, rank+1, 2*rank, MPI_COMM_WORLD, ierr )
		  	!print*, " send", rank
		  	end if 
		  	if (rank.gt.0) then 
		  	call MPI_RECV( d(1:2), 2, MPI_CTYPE, rank-1, 2*(rank-1), MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
		  	call MPI_RECV( u(1) ,1, MPI_CTYPE, rank-1, 2*(rank-1), MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
		  	call MPI_RECV( dd(2) ,1, MPI_CTYPE, rank-1, 2*(rank-1), MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr )
		  	!print*, " recv", rank
		  	end if 
		  else if(side .eq. 1) then 
		  	if (rank.gt.0) then 
		  	call MPI_SEND( d(1:2), 2, MPI_CTYPE, rank-1, 2*rank+1, MPI_COMM_WORLD, ierr )
		  	call MPI_SEND( u(1) ,1, MPI_CTYPE, rank-1, 2*rank+1, MPI_COMM_WORLD, ierr )
		  	call MPI_SEND( dd(2) ,1, MPI_CTYPE, rank-1, 2*rank+1, MPI_COMM_WORLD, ierr )
		  	end if 
		  	if (rank.lt.nprocs-1) then 
		  	call MPI_RECV( d(maxdim-1:maxdim), 2, MPI_CTYPE, rank+1, 2*(rank+1)+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
		  	call MPI_RECV( u(maxdim)  ,1, MPI_CTYPE,         rank+1, 2*(rank+1)+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
		  	call MPI_RECV( dd(maxdim) ,1, MPI_CTYPE,         rank+1, 2*(rank+1)+1, MPI_COMM_WORLD,MPI_STATUS_IGNORE, ierr )
		  	end if 
		  end if 
		 ! print*, "exchanged", rank
		  
		  !CPU to GPU
		  u_d(:maxdim-1) =  u(:maxdim-1)
		  d_d(maxdim-1:maxdim) = d(maxdim-1:maxdim)
		  dd_d(:maxdim) = dd(:maxdim)
		 
		  if (localpivot .ne. 0) then 
		 ! print*, maxit, oddeven,side,localpivot,u(localpivot),uu(localpivot),dd(localpivot+1)
		  end if 
		  !-------------------------------------------------------------------------
		  
			do ii = 1,maxdim
			M (ii,ii) = d_d(ii)
			end do 
			do ii = 2,maxdim
			M (ii,ii-1) = dd_d(ii)
			end do 
			do ii = 1,maxdim-2
			M (ii,ii+2) = uu_d(ii)
			end do 
			do ii = 1,maxdim-1
			M (ii,ii+1) = u_d(ii)
			end do 
			
			do jj = 0, nprocs-1
			if (rank.eq.jj)then 
			print*,"RANK" ,rank
			print*, " " 
			do ii = 1,maxdim 
			
			print "(*('('sf6.2xspf6.2x'i)':x))", M(ii,:)
			!write(*,*), M(ii,:)
			end do 
			
			print*, " " 
			end if
			CALL SLEEP(1)
			call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
			print*, " "
			end do 
			
		  
		  !-------------------------------------------------------------------------
	 if ( it .eq. 10 ) then  
	 	!exit 
	 end if 
  end do 	
  !call saxpy<<<grid, tBlock>>>(d_d,u_d,c_d,s_d,dd_d,uu_d, e,sens,pivot,side,oddeven)
  !y = y_d
  !write(*,*) 'Max error: ', maxval(abs(y-4.0))


end subroutine 


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

    
end module 
