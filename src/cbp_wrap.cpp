#include <cstddef>
#include <cstdio>

#ifdef CBPREAL
#define WTYPE float
#define PARAMS float alpha =1.0, beta = 0.0;
#endif

#ifdef CBPCMPLX
#define WTYPE cuComplex
#define PARAMS cuComplex alpha{1.0,0.0}, beta {0.0,0.0};
#endif

#ifdef CBPZMPLX
#define WTYPE cuDoubleComplex
#define PARAMS cuDoubleComplex alpha{1.0,0.0}, beta {0.0,0.0};
#endif

#ifdef CBPDOUBLE
#define WTYPE double
#define PARAMS double alpha =1.0, beta = 0.0;
#endif

#include "gemm.cpp"
// notice that the array becomes a line, saved row1-row2...
//integer*4-int
//cmplx-real and real
//zmplx-double float and double float
//real-float
//double-double




extern "C" {
        //test function
        int cbpfor_auto_gemm_wrapped_( WTYPE* a, WTYPE* b, WTYPE* c, int* dimA, int* dimB, int* dimC,int*add, int* ngpu)
        {
		PARAMS
		
		//printf("%d\n", *add);	
		if(*add == 1)
                {
                        beta = alpha;
                }
                else if(*add != 0)
                {
                        printf("Error, Not expected parameter\n");
                }

                
                return cbp_auto_gemm(a,b,c,dimA,dimB,dimC,*ngpu, alpha,beta);

        }
}

