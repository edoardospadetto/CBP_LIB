
  do ii = 0, nprocs-1
  	if (rank.eq.ii ) then 
 	do jj = 1, svddesc(1)
 		!print*, rank, d(jj) , svddesc(3) + jj -1 
 	end do 
 	end if 
  end do 
  if (rank.eq.nprocs-1) then 
  	descriptorD(1) = descriptorD(1)+2 
  end if
  
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  
  	!Obtain angain diagonal and upperdiagonal of A matrix	
  	
	call cbp_setget_diagonal(A,descriptorA, dd, descriptoru, 0, 'get')
	call cbp_setget_diagonal(A,[descriptorA(1),descriptorA(2),descriptorA(3),descriptorA(4)-1],& 
				 uu, descriptoru, 0, 'get')

	descriptoru(2)= max(descriptorA(3), descriptorA(4) -1)
	descriptoru(1) = minval(descriptorA(1:2)+ [descriptorA(3), descriptorA(4) -1] - descriptoru(2) )
	descriptoru(1) = max(0, descriptoru(1))
	
	
	descriptord(2) = max(descriptorA(3), descriptorA(4))
	descriptord(1) = minval(descriptorA(1:2)+ descriptorA(3:4)-descriptord(2) )
	descriptord(1) = max(0, descriptord(1))
	
	call cbp_redispatch_1d(uu, descriptoru, [svddesc(1), svddesc(3)],rank,nprocs) 
	call cbp_redispatch_1d(dd, descriptord, [svddesc(1), svddesc(3)],rank,nprocs)
	
	!compute AA+
	
	do ii = 1, maxdim
		dd(ii) = conjg(dd(ii))*dd(ii)  + conjg(uu(ii))*uu(ii) !diagonal
		uu(ii) = conjg(dd(ii+1))*uu(ii)  !upper diagonal
	end do 
	
	u(2:) = conjg(uu(:maxdim-1)) !lowerdiagonal
	u(1) = CTYPE0
	
	!if (rank.eq.2) then 
	!do kk = 0 , nprocs-1
	!if (rank.eq.kk ) then 
	!print*, rank, "-- u --", u ,"-- u --"
	!print*, rank, "-- uu --",uu,"-- uu --"
	!print*, rank ,"-- dd --",dd,"-- dd --"
	!print*, "  "
	!end if 
	!call sleep(1)
	!end do 
	! dd is diagonal, uu upperdiagonal, u lowerdiagonal
	!end if 
	descriptorD = [min(ldiag,minval(dimA) - rank*ldiag), svddesc(3) ] !descriptor of singular values

	do ii = 0, nprocs-1
		do jj = 1, min(ldiag, minval(dimA) - ldiag*ii) !descriptord(1) 

			if ( rank.eq. ii ) then
				call MPI_BCAST(d(jj),1, MPI_CTYPE, ii ,MPI_COMM_WORLD,ierr )
				temp = d(jj)
				!print*, rank, temp, jj , "---"
			else 
				call MPI_BCAST(temp,1, MPI_CTYPE, ii ,MPI_COMM_WORLD,ierr )
			end if 
			
			!------------ RESOLVE KER -----------------------------------------	
			!solve and pass first val to next process 
			if (rank.eq.0 ) then 
			ker_condition = [CTYPE0,CTYPE1]
			descriptorevec = [maxdim,1]
			!PRINT*,"--ii>>", descriptorevec
			call svd_tridiagonal_ker(u,dd-temp**2,uu,evec,maxdim-1,ker_condition)
			end if 
			
			if(nprocs .gt. 1 ) then 
				if (rank.eq.0 ) then 
				call MPI_SEND(evec(maxdim-1:maxdim), 2 , MPI_CTYPE, rank+1, rank, MPI_COMM_WORLD, ierr)
				else if (rank.eq.1) then 
				call MPI_RECV(ker_condition, 2 , MPI_CTYPE, rank-1, rank-1, MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )
				end if 
				
				do kk = 1, nprocs-2
					if (rank.eq.kk) then 
					evec(1) = ker_condition(1)
					descriptorevec = [maxdim, svddesc(3)]
					call svd_tridiagonal_ker(u(2:),dd(2:)-temp**2,uu(2:),evec(2:),maxdim-2,ker_condition)
					call MPI_SEND(evec(maxdim-2:maxdim-1) , 2 , MPI_CTYPE, rank+1, rank, MPI_COMM_WORLD,  ierr )
					else if (rank.eq.kk+1) then 
					call MPI_RECV(ker_condition , 2 , MPI_CTYPE, rank-1, rank-1, MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )
					end if 
			       end do 
			       
			       if (rank .eq. nprocs-1) then  
			       		descriptorevec = [maxdim, svddesc(3)]
			       		evec(1) = ker_condition(1)
			       		call svd_tridiagonal_ker(u(2:),dd(2:)-temp**2,uu(2:),evec(2:),maxdim-2,ker_condition)
			       end if 
		       end if
		       call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
		       if (rank.eq.0) then 
		       print*, "done", ii*ldiag+jj-1
		       end if
		       !print*, rank, evec(:maxdim-1), descriptorevec !rank, '|', descriptorevec, '!', svddesc(3), svddesc(1)
		       !print*," " 
		       !--------------------------------------------------------------------
		       !test EVEC
		       !do kk = 1, nprocs-1
		      ! if (rank .eq. kk) then 
		       		!print*, "s",  evec(2:maxdim-1), kk , maxdim-2
		       !		call MPI_SEND(evec(2:maxdim-1), maxdim-2, MPI_CTYPE, 0, rank, MPI_COMM_WORLD,  ierr)
		       		
		       !else if (rank.eq.0) then 
		       !		dd_end = [3+ldiag*kk , min( maxval(dimA) , ldiag*(kk+1) +2) ]
		       		!print*, "r", kk ,   dd_end(2)-dd_end(1)+1
		       !		call MPI_RECV(testevect( dd_end(1) : dd_end(2)), &
		       !			                 dd_end(2)-dd_end(1)+1, & 
		       !			                 MPI_CTYPE, kk, kk, MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierr)
		       	!		      
		       !end if
		       !if (rank.eq.0) then 
		       !	testevect(:ldiag+2) = evec(:ldiag+2)
		       !end if 
		      
		       !call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
		       !end do 
		       !PRINT*,"-->>", descriptorevec
		              
		       
		       evec = evec/cbp_l2norm(evec,1,descriptorevec(1),rank,nprocs)
		       !--------------- DEBUG EIGENVECTOORR----------------
		       testevect(:size(evec)) = evec
		       call cbp_redispatch_1d(testevect, descriptorevec, [dimA(1), 1],rank,nprocs) 
		       
		       !if (rank.eq.0) then 
		       ! 	print*,"----------->", testevect
		       !	print*, " " 
		       !end if 
		       
		       call GetBackMatrix(A, descriptorA, M, rowproc, colproc,& 
         	  	idxrow, idxcol,rank) 
         	  	
         	  	!call MPI_BARRIER(MPI_COMM_WORLD, ierr) 	
         	  	 
         	  	if (rank.eq.0) then 
         	  		M = matmul(M, conjg(transpose(M)))
         	  		trace = CTYPE0
         	  		print*, ii*ldiag+jj, abs(dot_product(testevect(:dimA(1)),testevect(:dimA(1))))
         	  		ortho(:,ii*ldiag+jj) = testevect(:dimA(1))
         	  		do kk = 1, dimA(1) 
         	  			trace = trace + M(kk,kk) 
         	  			
         	  			!print "(*('('sf8.2xspf8.2x'i)':x))", M(kk,:)
         	  	 		print*, abs(sum(M(kk,:)*testevect(:dimA(1))) - testevect(kk) * temp**2) ,  testevect(kk)
         	  		end do 
         	  		do kk = 1, ii*ldiag+jj
         	  		!print*, " orthotest " , abs(dot_product(ortho(:,kk),ortho(:,ii*ldiag+jj))), kk, ii*ldiag+jj
         	  	!	do kk = 1, dimA(1) 
         	  			!M(kk,kk) = M(kk,kk) - temp**2
         	  			!print "(*('('sf8.2xspf8.2x'i)':x))", M(kk,:)
         	  			!print*, "--" , M(kk,max(1,kk-1):min(dimA(1),kk+1))
         	  			!print*, testevect(max(1,kk-1):min(dimA(1),kk+1)), "--" ,sum(M(kk,:)*testevect(:dimA(1)))!, * temp**2 , testevect(kk)
         	  			!print*, "  "
         	  			
         	  		end do 
         	  	end if 
		         !----------------------------------
		       !call sleep(100) 
		       !descriptorTtemp = descriptorU0
			 Ttemp = CTYPE0
		     
	!	       print*, "norm" , cbp_l2norm(evec,1,descriptorevec(1),rank,nprocs),& 
	!	        sqrt(dot_product(evec(:descriptorevec(1)),evec(:descriptorevec(1))))
		       !U0(:,jj) = evec(:maxdim) 
		       !print*, "---Z", evec
		        do kk = 0, nprocs-1
		        	if (rank.eq.kk) then 
		        		!print*, "$$----Z>", [descriptorevec(1),1,descriptorevec(2),svddesc(4)+jj-1]
		        	end if 
  				call cbp_setgetprcline(U0, descriptorU0, evec, & 
  				[descriptorevec(1),1,descriptorevec(2),svddesc(3)+jj-1] , kk, rank, nprocs, 'set')
  			end do 
  			   !do kk = 1, dimA(1) 
         	  		!print "(*('('sf8.5xspf8.5x'i)':x))", U0(kk,:)
         	  	!end do 
  			if (+jj-descriptorTtemp(4) .ge. 1 .and. +jj-descriptorTtemp(4) .le. sum(descriptorTtemp(2:4:2))-1 ) then 
  				Ttemp( : , svddesc(4)+jj-descriptorTtemp(4) ) = U0( :, svddesc(4)+jj-descriptorU0(4) ) / temp 
  			end if 
		 
		 end do !end jj  
		 call cbp_adjoint(Ttemp, descriptorttemp)
		 
		 call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
		 
		 if (rank.eq.0) then 
		       print*, "done2", ii*ldiag+jj-1, ii, jj, ldiag
		 end if
		 !From now T temp will be globally rows TTemp = colA, and cols TTemp = rowA
		 !Pad with zeros or cancel columns   
		 !if (descriptorTtemp(3)-dimA(2)+1 .ge. 1 .and. descriptorTtemp(3)-dimA(2)+1 .le. sum(descriptorTtemp(1:3:2))-1 ) then 
		 !	TTemp(descriptorTtemp(3)-dimA(2)+1:,:) = CTYPE0
		 !end if 
		 
		 !call init_process_grid(process_grid, rowproc, colproc,rank,idxrow, idxcol)
		 !tmpdesc = descriptorTtemp
		 !svddesc = descriptorA
		! call cbp_matmul_descriptors(descriptorTtemp, descriptorA, descriptorV0,& 
		!				    [dimA(2),dimA(1)], dimA, [dimA(2),dimA(2)], & 
		!				    ld1, ld2, ld3 , rowproc, & 
		!				    colproc, process_grid,rank) 
	         !call cbp_redispatch(Ttemp,tmpdesc,descriptorTtemp, rank,nprocs)
	         !call cbp_redispatch(A,svddesc,descriptorA, rank,nprocs)
	         
	         !call cbp_matmul(Ttemp,A,V0,descriptorTTemp,descriptorA,descriptorV0, & 
		!			rowproc, colproc, rank, process_grid, context) 
			
		     
		
