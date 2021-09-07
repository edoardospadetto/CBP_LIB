#ifdef CBPREAL
#define BLASGEMM sgemm_
#define cputype real
#include "blas_gemm_routines/sgemm.c"
#endif
#ifdef CBPCMPLX
#define BLASGEMM cgemm_
#define cputype complex
#include "blas_gemm_routines/cgemm.c"
#endif
#ifdef CBPZMPLX
#define BLASGEMM zgemm_
#define cputype doublecomplex
#include "blas_gemm_routines/zgemm.c"
#endif
#ifdef CBPDOUBLE
#define BLASGEMM dgemm_
#define cputype double 
#include "blas_gemm_routines/dgemm.c"
#endif




int cpu_gemm(char *transa, char *transb, integer *m, integer *
        n, integer *k, cputype *alpha, cputype *a, integer *lda, cputype *b, integer *
        ldb, cputype *beta, cputype *c, integer *ldc)
{
	
	/* n = rows of a 
	 * m = cols of b 
	 * k = rows of b
	*/
	
//	omp_set_nested(1);
	bool nota = (*transa == 'N');
	bool notb = (*transb == 'N');
	
	/*if(nota*notb)
	{
		int t = 1;
		//look columns of columns of B
		//#pragma omp parallel for  
		for (int i =0 ; i<(*n) ; i++)
		{	
				
			BLASGEMM(transa, transb, m,
        		&t, k, alpha, a, lda, b+i*(*k), 
        		ldb, beta, c+i*(*m), ldc);
			
		}
	}
	else*/ 
	{
		BLASGEMM(transa, transb, m,
                        n, k, alpha, a, lda, b,
                        ldb, beta, c, ldc);
	}
			
	return 0;

}


