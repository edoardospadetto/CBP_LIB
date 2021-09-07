module misc_mod
    use, intrinsic :: iso_fortran_env, only: IK=>int32, RK=>real64
    implicit none
contains
    subroutine cbp_sort(n,Array,Index)
        implicit none
        integer(IK), intent(in)  :: n
        integer(IK), intent(in)  :: Array(n)
        integer(IK), intent(out) :: Index(n)
        integer(IK), parameter   :: nn=15, nstack=50
        integer(IK)              :: k,i,j,indext,jstack,l,r
        integer(IK)              :: istack(nstack)
        real(RK)                 :: a
        do j = 1,n
            Index(j) = j
        end do
        jstack=0
        l=1
        r=n
        do
            if (r-l < nn) then
                do j=l+1,r
                    indext=Index(j)
                    a=Array(indext)
                    do i=j-1,l,-1
                        if (Array(Index(i)) <= a) exit
                        Index(i+1)=Index(i)
                    end do
                    Index(i+1)=indext
                end do
                if (jstack == 0) return
                r=istack(jstack)
                l=istack(jstack-1)
                jstack=jstack-2
            else
                k=(l+r)/2
                call swap(Index(k),Index(l+1))
                call exchangeIndex(Index(l),Index(r))
                call exchangeIndex(Index(l+1),Index(r))
                call exchangeIndex(Index(l),Index(l+1))
                i=l+1
                j=r
                indext=Index(l+1)
                a=Array(indext)
                do
                    do
                        i=i+1
                        if (Array(Index(i)) >= a) exit
                    end do
                    do
                        j=j-1
                        if (Array(Index(j)) <= a) exit
                    end do
                    if (j < i) exit
                    call swap(Index(i),Index(j))
                end do
                Index(l+1)=Index(j)
                Index(j)=indext
                jstack=jstack+2
                if (jstack > nstack) then
                    write(*,*) 'NSTACK too small in indexArrayReal()'   ! xxx
                    error stop
                end if
                if (r-i+1 >= j-l) then
                    istack(jstack)=r
                    istack(jstack-1)=i
                    r=j-1
                else
                    istack(jstack)=j-1
                    istack(jstack-1)=l
                    l=i
                end if
            end if
        end do
        contains 

        subroutine exchangeIndex(i,j)
            integer(IK), intent(inout) :: i,j
            integer(IK)                :: swp
            if (Array(j) < Array(i)) then
                swp=i
                i=j
                j=swp
            end if
        end subroutine exchangeIndex
        
        pure elemental subroutine swap(a,b)
            implicit none
            integer(IK), intent(inout) :: a,b
            integer(IK) :: dum
            dum=a
            a=b
            b=dum
        end subroutine swap
    
        end subroutine 
end module misc_mod
