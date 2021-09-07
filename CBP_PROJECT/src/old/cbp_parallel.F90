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
        use iso_c_binding
        implicit none
       
        contains
        
        subroutine destroy_xt_context(h)
                type(c_ptr) :: h

                call destroy_xt_context_wrapped(h)

        end subroutine
        
        subroutine init_xt_context(ngpu, h)
               integer:: ngpu 
               type(c_ptr) :: h

              call init_xt_context_wrapped(ngpu, h) 

        end subroutine
        subroutine cbpfor_auto_gemm(A,B,C,add,ngpu,cpuratio,h)
            CTYPE, dimension(:,:) :: A,B,C
            integer :: dimA(2), dimB(2), dimC(2)
            integer :: ngpu, grid(2)
            integer :: add
            real:: cpuratio
            type(c_ptr):: h
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
            
            call cbpfor_auto_gemm_wrapped(A,B,C,dimA,dimB,dimC,add,ngpu, cpuratio, h)
            
        end subroutine 

end module 
