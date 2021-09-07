                                !do kk = 1, 4 
                                !        orderr(kk)=kk
                                !        orderc(kk)=kk
                                !end do
                                !print*, orderr, orderc , "b"
                                !do while (.not.rowok)
                                !        rowok = .true.
                                !        do kk = 1, size(row_lims, dim=1)-1
                                !                if(row_lims(kk) .gt. row_lims(kk+1) ) then 
                                !                        rowok = .false.
                                !                        tmp = row_lims(kk) 
                                !                        row_lims(kk) = row_lims(kk+1)
                                !                        row_lims(kk+1) = tmp
                                !
                                !                        tmp = orderr(kk)
                                !                        orderr(kk)=orderr(kk+1)
                                !                        orderr(kk+1)=tmp
                                !                 end if
                                !        end do 
                                !end do

                                !do while (.not. colok) 
                                 !       colok = .true.
                                 !       do kk = 1, size(col_lims, dim=1) -1
                                 !               if(col_lims(kk) .gt. col_lims(kk+1) ) then 
                                 !                       colok = .false.
                                 !                       tmp = col_lims(kk) 
                                 !                       col_lims(kk) = col_lims(kk+1)
                                 !                       col_lims(kk+1) = tmp

                                !                        tmp = orderc(kk)
                                !                        orderc(kk)=orderc(kk+1)
                                !                        orderc(kk+1)=tmp
                                !                 end if
                                !         end do 
                                !end do


                                !call cbp_sort(4,row_lims,orderr)
                                !call cbp_sort(4,col_lims,orderc)

                                !do kk = 1, 4
                                !        orderr(kk) = mod(orderr(kk),2) 
                                !        orderc(kk) = mod(orderc(kk),2)
                                !end do   
                                !print*, orderr,orderc, "a"



                                ! do kk = 1, 4
                               !         orderr(kk)=kk
                               !         orderc(kk)=kk
                               ! end do
                                !print*, orderr, orderc , "b"
                                !do while (.not.rowok)
                                !        rowok = .true.
                                !        do kk = 1, size(row_lims, dim=1)-1
                                !                if(row_lims(kk) .gt. row_lims(kk+1) ) then
                                !                        rowok = .false.
                                !                        tmp = row_lims(kk)
                                !                        row_lims(kk) = row_lims(kk+1)
                                !                        row_lims(kk+1) = tmp
                                !
                               !                         tmp = orderr(kk)
                               !                         orderr(kk)=orderr(kk+1)
                               !                         orderr(kk+1)=tmp
                                !                 end if
                                !        end do
                                !end do

                                !do while (.not. colok)
                                  !      colok = .true.
                                  !      do kk = 1, size(col_lims, dim=1) -1
                                  !              if(col_lims(kk) .gt. col_lims(kk+1) ) then
                                  !                      colok = .false.
                                  !                      tmp = col_lims(kk)
                                  !                      col_lims(kk) = col_lims(kk+1)
                                 !                       col_lims(kk+1) = tmp

                                !                        tmp = orderc(kk)
                                !                        orderc(kk)=orderc(kk+1)
                                !                        orderc(kk+1)=tmp
                                !                 end if
                                !         end do
                                !end do


                                !call cbp_sort(4,row_lims,orderr)
                                !call cbp_sort(4,col_lims,orderr)
                                !do kk = 1, 4
                                !        orderr(kk) = mod(orderr(kk),2)
                                !        orderc(kk) = mod(orderc(kk),2)
                                !end do

                                !if ( xor(orderr(1), orderr(2)) == 1 .and. xor(orderc(1), orderc(2)) == 1 ) then
-- INSERT --
