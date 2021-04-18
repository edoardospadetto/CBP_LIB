#include "cbp_parallel.F90"


program test 
    use cbpfor_parallel
    real, dimension(1000,1000) ::e,f
    CTYPE ,dimension(1000,1000) :: a , b , c, d
    integer :: ii , jj 
    !do ii = 1, 100
    !    do jj = 1, 100 
    !        a(ii,jj) = 1.0
    !        b(ii,jj) = 1.0
    !    end do 
    !end do 
    call random_number(e)
    call random_number(f)
    A = e
    call random_number(e)
    call random_number(f)
    B = f
    d= matmul(a,b)
    
    call cbpfor_auto_gemm(A,B,C,0,4)
   
    call random_number(e)
    call random_number(f)
    A = e
    call random_number(e)
    call random_number(f)
    B = f
    d= d+matmul(a,b)
   
    call cbpfor_auto_gemm(A,B,C,1,4)

    print*, "result ",  (sum(c)), sum(d)
      
    
    
        
end program 
