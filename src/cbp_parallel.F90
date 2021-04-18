#ifdef CBPREAL
#define CTYPE   REAL
#define CTYPE1 1.0
#define CTYPE0 0.0
#endif

#ifdef CBPCMPLX
#define CTYPE complex
#define CTYPE1 cmplx(1.0,0.0)
#define CTYPE0 cmplx(0.0,0.0)
#endif

#ifdef CBPZMPLX
#define CTYPE double complex
#define CTYPE1 dcmplx(1.0,0.0)
#define CTYPE0 dcmplx(0.0,0.0)
#endif

#ifdef CBPDOUBLE
#define CTYPE double precision
#define CTYPE1 1.d0
#define CTYPE0 0.d0
#endif



module cbpfor_parallel 
        implicit none 
        contains 
        subroutine cbpfor_auto_gemm(A,B,C,add,ngpu)
            CTYPE, dimension(:,:) :: A,B,C
            integer :: dimA(2), dimB(2), dimC(2)
            integer :: ngpu, grid(2)
            integer :: add
            dimA = shape(A)
            dimB = shape(B)
            dimC = shape(C)
            
            if (dimA(2) .ne. dimB(1)) then
                STOP
            else if (dimA(1) .ne. dimC(1)) then
                STOP
            else if (dimB(2) .ne. dimC(2)) then 
                STOP
            end if
            
            call cbpfor_auto_gemm_wrapped(A,B,C,dimA,dimB,dimC,add,ngpu)
            
        end subroutine 

end module 
