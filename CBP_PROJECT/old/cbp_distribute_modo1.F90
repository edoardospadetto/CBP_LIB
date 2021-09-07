#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE0 cmplx(0.0,0.0)
#define MPI_CTYPE MPI_COMPLEX
#endif

#ifdef CBPDOUBLE 
#define CTYPE double precision
#define CTYPE0 0.0

#define MPI_CTYPE MPI_DOUBLE_PRECISION
#endif

#ifdef CBPZMPLX 
#define CTYPE0 dcmplx(0.0,0.0)

#define CTYPE double complex
#define MPI_CTYPE MPI_DOUBLE_COMPLEX
#endif

#ifdef CBPREAL 
#define CTYPE0 0.0
#define CTYPE real
#define MPI_CTYPE MPI_REAL
#endif

include 'misc_module.F90'

module cbp_DistributeMOD
    use mpi
    use cbpfor_parallel
    use misc_mod

    implicit none 
    contains
    
    !Just print a local Matrix on terminal
        SUBROUTINE printmatrix(mat)
          IMPLICIT NONE
          INTEGER :: ii,jj
          INTEGER, DIMENSION(:,:), INTENT(IN) :: mat
          INTEGER, DIMENSION(SIZE(SHAPE(mat))) :: sizem
          sizem=SHAPE(mat)

          PRINT*, ""
          DO ii=1,sizem(1)
        DO jj=1,sizem(2)
          WRITE(*, fmt="(I2, tr2)", advance="no") mat(ii,jj)
        END DO
        PRINT*, ""
          END DO


        END SUBROUTINE printmatrix
    
    
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
        integer , dimension(:,:) :: process_grid
        integer , dimension(:), allocatable :: setofpieces
        integer :: ii, jj, kk, rowproc, colproc, ierr, rank
        integer :: idxrow, idxcol
        
        idxcol= rank/rowproc+1            
        idxrow= mod(rank,rowproc)+1
        
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
    
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
    !Transpose distributed matrix. 
    
    subroutine cbp_transpose(localMat, descriptor)
        
        CTYPE , dimension(:,:,:) :: localMat 
        integer , dimension(4,2) :: descriptor
        integer :: row, col, ii
        CTYPE :: temp
        
        if( size(localMat(:,:,1), dim = 1 ) .eq. size(localMat(:,:,1) , dim= 2 ))  then 
                !transpose matrix
                
                localMat(:,:,2) = transpose(localMat(:,:,1))
                localMat(:,:,1) = localMat(:,:,2)
                !ii = descriptor(1,1) * descriptor(2,1)
               ! print*, descriptor(1,1), descriptor(2,1), (ii-1)/descriptor(2,1) +1, mod(ii-1,descriptor(2,1))+1 
               ! !$omp parallel do private(temp,row,col)
               ! do ii = 1, size(localMat(:,:,1), dim = 1)*size(localMat(:,:,1), dim = 2)
               !          
               !         row = (ii-1)/size(localMat(:,:,1), dim = 1) +1
               !         col = mod(ii-1,size(localMat(:,:,1), dim = 1))+1
               !         if(row .gt. col) then  
               !                 temp = localMat(row,col,1)
               !                 localMat(row,col,1) = localMat(col,row,1)
               !                 localMat(col,row,1) = temp
               !         end if
               ! end do 
               ! !$omp end parallel do
                
                !transpose descriptor dim 
                descriptor(1,2) = descriptor(2,1)
                descriptor(2,2) = descriptor(1,1)
        
                !transpose descriptor starting indexes
                descriptor(3,2) = descriptor(4,1)
                descriptor(4,2) = descriptor(3,1)


                !reset descriptor 
                descriptor(:,1) = descriptor(:,2)
        else  
                print*, "Error, Transposition not safe" 
        end if 
    end subroutine
        
    
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !given new descriptors and old descriptor 
    !redispatch a matrix according to the new descriptors
    subroutine cbp_redispatch(localMat, olddescriptor, newdescriptor, rank, nprocs)
        integer , dimension(4) :: olddescriptor , newdescriptor, other
        CTYPE , dimension(:,:,:) :: localMat 
        integer :: rank, nprocs, ii, jj, kk, ierr
        integer , dimension(4) :: orderr, orderc
        integer , dimension(4) :: row_lims, col_lims
        integer, dimension(MPI_STATUS_SIZE) :: stat
        logical :: colok, rowok
        integer :: tmp
        localMat(:,:,2) = localMat(:,:,1)
        localMat(:,:,1) = 0.0

        do ii = 0, nprocs-1 
                do jj = 0 , nprocs-1 
               
                        if (rank .eq. ii) then 
                                !send my new descriptor and receive its old one.i
                                !print*, ii,jj, newdescriptor, 
                                !print*, ii,jj, olddescriptor
                                if(ii .ne. jj) then 
                                !print*, ii, jj, "firstmet r"
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

       
                                      ! print*,ii, jj , "||"   ,max(newdescriptor(3),other(3)) , &
                                      !           min(newdescriptor(3)+newdescriptor(1),other(3)+other(1))-1,&  
                                      !           max(newdescriptor(4),other(4)), & 
                                      !           min(newdescriptor(4)+newdescriptor(2),other(4))+other(2)-1,"r"
                                         
                                        if (rank .ne. jj ) then
                                        call  MPI_RECV( & 
                                        localMat(max(newdescriptor(3),other(3))-newdescriptor(3)+1 :& 
                                                 min(newdescriptor(3)+newdescriptor(1),other(3)+other(1))-newdescriptor(3), &  
                                                 max(newdescriptor(4),other(4))-newdescriptor(4)+1:& 
                                                 min(newdescriptor(4)+newdescriptor(2),other(4)+other(2))-newdescriptor(4),&
                                                 1) , &
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
                                                 olddescriptor(4)+olddescriptor(2))-newdescriptor(4),&
                                                 1) = & 
                                         localMat(max(olddescriptor(3),newdescriptor(3))-olddescriptor(3)+1:& 
                                                 min(olddescriptor(3)+olddescriptor(1),&
                                                 newdescriptor(3)+newdescriptor(1))-olddescriptor(3), &  
                                                 max(olddescriptor(4),newdescriptor(4))-olddescriptor(4)+1:& 
                                                 min(olddescriptor(4)+olddescriptor(2),&
                                                 newdescriptor(4)+newdescriptor(2))-olddescriptor(4),&
                                                 2)  
                                        
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
                                       !          min(olddescriptor(3)+olddescriptor(1),other(3)+other(1))-1, &  
                                       !          max(olddescriptor(4),other(4)),& 
                                       !          min(olddescriptor(4)+olddescriptor(2),other(4)+other(2))-1,"s"
                                       ! print*, ii, jj , "overlap  send"         
                                        call MPI_SEND( & 
                                        localMat(max(olddescriptor(3),other(3))-olddescriptor(3)+1:& 
                                                 min(olddescriptor(3)+olddescriptor(1),other(3)+other(1))-olddescriptor(3), &  
                                                 max(olddescriptor(4),other(4))-olddescriptor(4)+1:& 
                                                 min(olddescriptor(4)+olddescriptor(2),other(4)+other(2))-olddescriptor(4),&
                                                  2) , &
                                                (min(olddescriptor(3)+olddescriptor(1),other(3)+other(1)) - & 
                                                 max(olddescriptor(3),other(3))) * & 
                                                (min(olddescriptor(4)+olddescriptor(2), other(4)+other(2)) - & 
                                                 max(olddescriptor(4),other(4))) ,& 
                                                 MPI_CTYPE, ii, ii, MPI_COMM_WORLD, ierr)   
                                               
                                
                                end if 

                              
                        end if  
                        

                end do 
        end do 



    end subroutine
        
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    !given row and column in the process grid computes the rank
    function grid2rank(row, col , rowproc, colproc) result(idx)
        implicit none
        integer :: row, col, idx, rowproc, colproc
    
        idx = (col-1)*rowproc + row-1
        
        return
        
    end function
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        
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
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Same as above but take in account more mpi calls to send or receive a 
!local matrix if it is too big
!Not working at the moment, a bug must be solved. 

subroutine MPI_comm_wrap_matmul(LocalMatrix, MatrixDescriptor ,& 
                      merow,mecol,otherrow,othercol,rowproc,colproc,& 
                      mod,tag ,requestM, requestD)
        
        implicit none
        CTYPE, dimension(:,:) :: LocalMatrix
        integer, dimension(4) :: MatrixDescriptor
        
        integer*8 :: dimstmp(2)
        integer :: ii, cut , out, tmp(2) 
        integer :: dims(2) 
        integer :: otherrow,othercol, ierr, tag, merow, mecol, rowproc, colproc
        integer, dimension(:), allocatable :: requestM 
        integer ::requestD
        character*1 :: mod
        integer :: maxint
        
        dimstmp(1) = size(LocalMatrix, dim = 1)
        dimstmp(2) = size(LocalMatrix, dim = 2)
        !print *,"dimension : ",  dims
        dims = dimstmp
        out = 1
        ii = 1
        

                if(dimstmp(1) .gt. huge(maxint)) then 
                        print*, "Error at 415:cbp_distribute.F90 : Matrix too big "
                end if
                cut = dimstmp(2)
                

                do while (cut*dimstmp(1) .gt. huge(maxint)) 
                        cut = ceiling(real(dimstmp(2))/real(ii))
                        out = ii
                        ii = ii +1
                        !print*, ii , cut, out
                end do
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
        
        allocate(requestM(out))

        do  ii = 1, out 
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
    


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
    


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òòòòòò@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine ComputeMatmul(localA, localB, localC, descriptorA, descriptorB , descriptorC, rowproc, colproc, rank, process_grid,& 
                ngpu, cpuratio, handle)
    use iso_c_binding
    CTYPE, dimension(:,:,:):: localA , localB 
    CTYPE, dimension(:,:):: localC
    integer , dimension(4,2) :: descriptorA, descriptorB
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
    
    
    integer , dimension(:,:), allocatable :: stat 
        
        
        
    
    ierr = 0
    allocate (old_process_grid(max(rowproc,colproc), max(rowproc, colproc)))
    
    
    !end here
    nprocs = colproc*rowproc
    pace = max(colproc,rowproc)
    switch = 1
    
    allocate(stat(pace,8))
    
    !Each process understand its identity
    idxcol= rank/rowproc+1            !grid assigned by cols
    idxrow= mod(rank,rowproc)+1
    
    localC = CTYPE0
    if(process_grid(idxrow,idxcol) .gt. 0  ) then 

    call cbpfor_auto_gemm(localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),&
         localB(:descriptorB(1,switch),:descriptorB(2,switch),switch),&
          localC(:descriptorC(1), :descriptorC(2)), 1, ngpu , cpuratio, handle)
        

            
    end if
    
    
    
    
    
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
                LocalB(:,:,mod(switch,2)+1) = LocalB(:,:,switch)
                descriptorB(:,mod(switch,2)+1) = descriptorB(:,switch)
                !print*, "Already here sendB"
            else 
            sb = .true.
            call  MPI_comm_wrap_matmul(LocalB(:,:,switch),descriptorB(:,switch),idxrow,idxcol,& 
                withwhorow, withwhocol,rowproc,colproc,'s',1, requestsendB, requestsendDB)
            !print*,'sb', idxrow, idxcol, withwhorow, withwhocol, ii,  '| b |', descriptorB(:,switch)
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
                LocalA(:,:,mod(switch,2)+1) = LocalA(:,:,switch)
                descriptorA(:,mod(switch,2)+1) = descriptorA(:,switch)
                !print*, "Already Here send A"
            else 
            sa = .true.
            !print*,'sa', idxrow, idxcol, withwhorow, withwhocol, ii , '| a |', descriptorA(:,switch) 
            call  MPI_comm_wrap_matmul(LocalA(:,:,switch),descriptorA(:,switch),idxrow,idxcol,&
                withwhorow, withwhocol,rowproc,colproc,'s',0,requestsendA, requestsendDA) 
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
                
                call MPI_comm_wrap_matmul(LocalB(:,:,switch),descriptorB(:,switch),idxrow,idxcol,& 
                    withwhorow, withwhocol,rowproc,colproc,'r',1, requestrecvB, requestrecvDB)
                
                !print*,'rb',withwhorow, withwhocol, idxrow, idxcol,  ii ,  '| b |'!, descriptorB(:,switch)
        
            end if
        end if
        !############################################################################################################!
        
        
        !############################################################################################################!
        !recv A--------------------------------------------------------------------------------------------------------
        call find_coop(old_process_grid, process_grid, idxrow, idxcol, withwhorow, withwhocol, rowproc, colproc, "ra")
        
        if ((withwhorow .gt. 0 ).and. (withwhocol .gt. 0 )) then
            
            if(.not. ( (withwhorow .eq. idxrow) .and. (withwhocol .eq. idxcol) ) ) then 
            ra = .true.
            
            call MPI_comm_wrap_matmul(LocalA(:,:,switch),descriptorA(:,switch),idxrow,idxcol,& 
                withwhorow, withwhocol,rowproc,colproc,'r',0, requestrecvA, requestrecvDA)
            
            
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
    
           
            call cbpfor_auto_gemm(localA(:descriptorA(1,switch),:descriptorA(2,switch),switch),& 
                                  localB(:descriptorB(1,switch),:descriptorB(2,switch),switch),& 
                                  localC(:descriptorC(1), :descriptorC(2)), 1, ngpu, cpuratio, handle) 
            
        end if
        
        !--------------------------------------------------------------------------
    end do 
    
    if(switch .eq. 2 ) then 
        
            descriptorA(:,1) = descriptorA(:,2)
            descriptorB(:,1) = descriptorB(:,2)
            localA(:,:,1) = localA(:,:,2)
            localB(:,:,1) = localB(:,:,2)

    end if 
    
    deallocate(stat)
    
    deallocate(old_process_grid)
    
    


end subroutine  


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Set entry given global indexes of local matrices and descriptor 
!Notice. localMat has rank2, the input of compute matmul have rank 3.
!same for the descriptors
!If you set before launching matmul use the localinput(:,:,1) and descriptor(:,1) (switch =1)
!the index at 1 is the switch used in the matmul to not overwrite local space. 
!At the start data is in switch = 1, during matmul switch changes between 1 and 2
! if you use this after the matmul with the same process grid you have to know the walue of the switch.

!Not tested yet. 

subroutine cbp_setentry(localMat, descriptor, process_grid, row,col, val)
        
        CTYPE, dimension(:,:) :: localMat
        CTYPE :: val
        integer, dimension(4) :: descriptor 
        integer, dimension(:,:) :: process_grid
        integer :: row, col
        
        if ( (descriptor(3) .le. row) .and. (descriptor(3)+descriptor(1) .gt. row)) then 
                    if ( (descriptor(4) .le. col) .and. (descriptor(4)+descriptor(2).gt. col) ) then 

                        localmat(row-descriptor(3)+1, col-descriptor(4)+1) = val        
                    end if
        end if
        


end subroutine

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Process #0 contains the global Matrix, request the descriptor from all the others, and according to it 
!send back the correct chunk of the global matrix

subroutine dispatch_global_matrix(A, localA, descriptor, rowproc, colproc,& 
         idxrow, idxcol,rank,process_grid)
        use mpi
        implicit none
        !Use the descriptor and the process grid slepy dumb idiot.
        CTYPE, dimension(:,:) :: A
        CTYPE, dimension(:,:) :: localA
        integer, dimension(:,:) :: process_grid
        integer, dimension(4) :: descriptor, extdescriptor
        integer :: rowproc, colproc, idxrow, idxcol, rank
        integer, dimension(MPI_STATUS_SIZE) :: status_
        
        integer :: ierr, ii,jj
        
    
        if(rank .eq. 0) then 
        
        do ii = 1, rowproc
            do jj = 1, colproc
                if(ii*jj .ne. 1) then 
                    
                    !if(process_grid(ii,jj) .gt. 0) then
                    
                    
                    call MPI_RECV(extdescriptor, 4, MPI_INTEGER, & 
                              grid2rank(ii,jj,rowproc,colproc) & 
                             ,3*(ii*(colproc+1)+jj)+1 , MPI_COMM_WORLD, status_, ierr)
                    call mpic_checkerror(ierr)
                    !print* , extdescriptor, ii, jj,  process_grid(ii, jj)
                    !if(process_grid(ii,jj) .gt. 0) then
                    if(sum(extdescriptor) .ne. 0 ) then
                    !print*, "test g", ii , jj,  A(extdescriptor(3), extdescriptor(4)),extdescriptor(1)*extdescriptor(2) !here the problem
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
        
        !if(process_grid(idxrow,idxcol) .gt. 0) then
        
        
        call MPI_SSEND( descriptor , 4,& 
                     MPI_INTEGER, 0 ,& 
                     3*(idxrow*(colproc+1)+idxcol)+1 , MPI_COMM_WORLD, ierr)
        
        call mpic_checkerror(ierr)
        !print* , descriptor, idxrow, idxcol, process_grid(idxrow, idxcol)
        !if(process_grid(ii,jj) .gt. 0) then
        if(sum(descriptor) .ne. 0 ) then
        !if(process_grid(idxrow,idxcol) .gt. 0) then
        call MPI_RECV(localA(:descriptor(1),:descriptor(2)), descriptor(1)*descriptor(2), MPI_CTYPE, 0 & 
            ,idxrow*(colproc+1)+idxcol , MPI_COMM_WORLD, status_, ierr)
        
        call mpic_checkerror(ierr)
        
        
    
        end if

        
        end if
        
        
                
        
    
    
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
                    
                    C(descriptor(3): descriptor(3)+descriptor(1)-1,& 
                           descriptor(4):descriptor(4)+descriptor(2)-1) = localC
                    
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


subroutine cbp_init_local_matrix_descriptor & 
        (descriptorA,dimA,localdimA,&  
        rowproc, colproc, rank) 
        
        implicit none
        
        integer , dimension(4, 2)  :: descriptorA
        integer , dimension(2)  :: dimA 
        integer , dimension(2)  :: localdimA
        integer :: rowproc, colproc
        integer :: idxcol, idxrow, rank
        
        idxcol= rank/rowproc+1            
        idxrow= mod(rank,rowproc)+1
        
        localdimA(1) = ceiling(real(dimA(1))/real(rowproc))
        localdimA(2) = ceiling(real(dimA(2))/real(colproc))
        
        localdimA = maxval(localdimA) 
        
        descriptorA(:,1) =  (/localdimA(1), localdimA(2), (idxrow-1)*localdimA(1)+1, & 
                (idxcol-1)*localdimA(1)+1 /) 
        
        !---LIMIT CASES---
        
        if(idxrow .eq. rowproc) then
                descriptorA(1,1) = dimA(1) - descriptorA(3,1) +1
        end if
        if(idxcol .eq. colproc) then
                descriptorA(2,1) = dimA(2) - descriptorA(4,1) +1
        end if 
        
end subroutine

    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine cbp_matmul_descriptors         & 
        (descriptorA, descriptorB, descriptorC,& 
        dimA, dimB, dimC,            &
        localdimA, localdimB, localdimC,       &  
        rowproc, colproc, process_grid ,rank) 
        
        implicit none
        
        integer, dimension(:,:) :: process_grid
        integer , dimension(4, 2)  :: descriptorA, descriptorB
        integer , dimension(4)     ::descriptorC
        integer , dimension(2)  :: dimA, dimB, dimC 
        integer , dimension(2)  :: localdimA, localdimB, localdimC
        integer :: rowproc, colproc
        integer :: idxcol, idxrow, rank
        
        idxcol= rank/rowproc+1            
        idxrow= mod(rank,rowproc)+1
        
        localdimA(1) = ceiling(real(dimA(1))/real(rowproc))
        localdimB(2) = ceiling(real(dimB(2))/real(colproc))
        localdimC(1) = localdimA(1)
        localdimC(2) = localdimB(2)
        localdimA(2) = ceiling(real(dimA(2))/real(min(colproc,rowproc)))
        localdimB(1) = localdimA(2)
        
        descriptorA(:,1) =  (/localdimA(1), localdimA(2), (idxrow-1)*localdimA(1)+1, & 
                (abs(process_grid(idxrow,idxcol))-1)*localdimA(2)+1 /) 
     
        descriptorB(:,1) =  (/localdimB(1), localdimB(2),& 
                 (abs(process_grid(idxrow,idxcol))-1)*localdimB(1)+1, (idxcol-1)*localdimB(2)+1 /) 
        
        descriptorC =  (/localdimC(1), localdimC(2), (idxrow-1)*localdimC(1)+1,& 
                 (idxcol-1)*localdimC(2)+1/) 
       
        localdimA = maxval(localdimA) 
        localdimB = maxval(localdimB)
        localdimc = maxval(localdimC)


        
        !---LIMIT CASES---
        
        if(idxrow .eq. rowproc) then
                descriptorA(1,1) = dimA(1) - descriptorA(3,1) +1
                descriptorC(1) = dimC(1) - descriptorC(3) +1
        end if
        if(idxcol .eq. colproc) then
                descriptorB(2,1) = dimB(2) - descriptorB(4,1) +1
                descriptorC(2) = dimC(2) -descriptorC(4)+1
        end if 
        if(abs(process_grid(idxrow,idxcol)) .eq. min(colproc,rowproc)) then 
                descriptorA(2,1) = dimA(2) - descriptorA(4,1) +1
                descriptorB(1,1) = dimB(1) - descriptorB(3,1) +1
            if(descriptorA(2,1) - descriptorB(1,1) .ne. 0) then
                    print*, "Erroror Init_matmul_descriptors, & 
                            Getting Wrong Sizes of local matrices blocks."
            end if 
        end if
        
        if( process_grid(idxrow, idxcol) .le. 0 ) then
            if (rowproc .lt. colproc) then
                descriptorA  = 0 
            else if (colproc .lt. rowproc) then
                descriptorB  = 0
            end if
            
        end if
        
end subroutine


!subroutine cbp_unbalanced_matmul(localA, descritptorA, B, localC, descriptorC)
!        integer, dimension(4) :: descriptorA, descriptorC
!
!        if (descriptorA() ) 
!
!end subroutine 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!descriptorv & descriptor res descriptor for unidimensional array 
!localv localres vectors, so matrices with 1column.

subroutine cbp_vec_matmul(localA, descriptorA, localv, descriptorv, localres, descriptorres, idxrow, idxcol, colproc, rowproc)
        
        CTYPE, dimension(:,:) ::localA
        CTYPE, dimension(:,:) ::localv , localres
        
        integer, dimension(2) :: descriptorv, descriptorres !Vector descriptors 
        integer :: idxrow, idxcol, rowproc, colproc

        integer, dimension(4) :: descriptortmpres, descriptortmpv, descriptorA !descriptor matrices
        integer, dimension(:), allocatable :: requestV
        integer :: requestDV
        integer :: rowcomm, ierr
        integer :: ii 
        type(c_ptr) :: handle

        if (size(localv, dim = 2) .ne. 1) then 
                print*, "Error localv do not have a vector shape"
        end if
        if (size(localres, dim = 2) .ne. 1) then 
                print*, "Error localv do not have a vector shape"
        end if        
        if ( descriptorA(2) .ne. descriptorv(1) ) then  !dimension
                print*, "Error , vector matrix multiplication , wrong dim"
        end if
        if ( descriptorA(3) .eq. descriptorv(2) ) then !start
                print*, "Error , vector matrix multiplication , wrong start idx"
        end if

     
        descriptortmpres = (/descriptorv(1) ,1, descriptorv(2) , 1/)   
        descriptortmpv = (/descriptorv(1) ,1, descriptorv(2) , 1/)   
        
        if (idxrow.eq.1) then 
        do ii = 2, rowproc

                call MPI_comm_wrap_matmul(LocalV, descriptortmpv ,&
                              idxrow,idxcol,ii,idxcol,rowproc,colproc,&
                              's', ii ,requestV, requestDV)
        end do
        else  
        call MPI_comm_wrap_matmul(LocalV, descriptortmpv ,&
                      idxrow,idxcol,1,idxcol,rowproc,colproc,&
                      'r', idxrow ,requestV, requestDV)

        end if

        call MPI_WAITALL( size(requestV,dim=1) , & 
                        requestV, MPI_STATUSES_IGNORE, ierr)

        call MPI_WAIT(requestDV, MPI_STATUS_IGNORE, ierr )
          
        descriptorres = (/descriptortmpv(1) , descriptortmpv(3) /)
        descriptorv = (/descriptortmpv(1) , descriptortmpv(3) /)


        call cbpfor_auto_gemm( localA(:descriptorA(1),:descriptorA(2)),&
                               localv(:descriptorv(1), :),&
                               localres(:descriptorres(1),: ) , 1, 1, 0.0, handle)

        deallocate (requestV)

        call MPI_Comm_split(MPI_COMM_WORLD, idxrow, idxcol, rowcomm, ierr)
        call MPI_reduce(localres, localv, descriptorres(1), MPI_CTYPE, MPI_SUM, 1, rowcomm, ierr )


        if (idxrow .eq. 1) then 
                localres(:descriptorres(1),:) = localV(:descriptorv(1),:)
                localv = CTYPE0 
        end if

end subroutine        
        
end module


