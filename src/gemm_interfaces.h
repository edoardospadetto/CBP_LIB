#include "BLAS_INC/cblas.h" 

template<typename T>
int cbp_Gemm(T* A, T* B, T* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, T alpha, T beta )
{
	printf("missing definition of gemm \n");
	return 1;
}
    
template<>
int cbp_Gemm(double* A, double* B, double* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, double alpha, double beta )
{
        
        
	if(cublasXtDgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}

template<>
int cbp_Gemm(float* A, float* B, float* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, float alpha, float beta )
{
       	if(cublasXtSetCpuRoutine(xt_,CUBLASXT_GEMM, CUBLASXT_FLOAT, nullptr )!= 0 ) printf("nogood\n");    
        if (cublasXtSetCpuRatio(xt_,CUBLASXT_GEMM, CUBLASXT_FLOAT,99.0)!= 0) printf("nogood\n");  
			
	if(cublasXtSgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
     	printf("here\n"); 
	
        cudaDeviceSynchronize();

        return 0;
}


template<>
int cbp_Gemm(cuComplex* A, cuComplex* B, cuComplex* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, cuComplex alpha, cuComplex beta )
{
        
        
	if(cublasXtCgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}


template<>
int cbp_Gemm(cuDoubleComplex* A, cuDoubleComplex* B, cuDoubleComplex* C, int rowa, int colb, int freedim, cublasXtHandle_t &xt_, cuDoubleComplex alpha, cuDoubleComplex beta )
{
        
        
	if(cublasXtZgemm(xt_, CUBLAS_OP_N, CUBLAS_OP_N,
                rowa, colb, freedim, &alpha, A, rowa, B, freedim, &beta, C, rowa) !=0  ) return 1;
        
      

        cudaDeviceSynchronize();

        return 0;
}
