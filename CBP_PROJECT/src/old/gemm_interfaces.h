 

#include "cpugemms/cpugemm.cpp"

/*AUTO GEMM IMPLEMENTATION TEMPLATE*/

template<typename T>
int cbp_Gemm(T* A, T* B, T* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, T alpha, T beta, float cpuratio )
{
	printf("missing definition of gemm \n");
	return 1;
}
    
template<>
int cbp_Gemm(double* A, double* B, double* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, double alpha, double beta, float cpuratio )
{
     	if(cpuratio >0)
	{
		if (cublasXtSetCpuRatio(xt_,CUBLASXT_GEMM, CUBLASXT_DOUBLE,0.5)!= 0) printf("nogood\n");
		if(cublasXtSetCpuRoutine(xt_,CUBLASXT_GEMM, CUBLASXT_DOUBLE, (void *) cpu_gemm )!= 0 ) printf("nogood\n");    
        }  
        
	if(cublasXtDgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}

template<>
int cbp_Gemm(float* A, float* B, float* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, float alpha, float beta ,float cpuratio)
{
 	
	if(cpuratio > 0)
	{
	  if (cublasXtSetCpuRatio(xt_,CUBLASXT_GEMM, CUBLASXT_FLOAT,cpuratio)!= 0) printf("nogood\n");
	  if(cublasXtSetCpuRoutine(xt_,CUBLASXT_GEMM, CUBLASXT_FLOAT, (void *) cpu_gemm )!= 0 ) printf("nogood\n");    
        }
			
	if(cublasXtSgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
     	//printf("here\n"); 
	
        cudaDeviceSynchronize();

        return 0;
}


template<>
int cbp_Gemm(cuComplex* A, cuComplex* B, cuComplex* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, cuComplex alpha, cuComplex beta, float cpuratio )
{
       if(cpuratio >0)
	{ 
       	if (cublasXtSetCpuRatio(xt_,CUBLASXT_GEMM, CUBLASXT_COMPLEX,cpuratio)!= 0) printf("nogood\n");
	if(cublasXtSetCpuRoutine(xt_,CUBLASXT_GEMM, CUBLASXT_COMPLEX, (void *) cpu_gemm )!= 0 ) printf("nogood\n");    
	} 
	
	if(cublasXtCgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}


template<>
int cbp_Gemm(cuDoubleComplex* A, cuDoubleComplex* B, cuDoubleComplex* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, cuDoubleComplex alpha, cuDoubleComplex beta, float cpuratio )
{
	if(cpuratio>0)
	{
      	 if (cublasXtSetCpuRatio(xt_,CUBLASXT_GEMM, CUBLASXT_DOUBLECOMPLEX,cpuratio)!= 0) printf("nogood\n");
	 if(cublasXtSetCpuRoutine(xt_,CUBLASXT_GEMM, CUBLASXT_DOUBLECOMPLEX, (void *) cpu_gemm )!= 0 ) printf("nogood\n");    
        }	 
        
	if(cublasXtZgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}
