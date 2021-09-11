#include <cstddef>
#include <cstdio>
#include <cuda_runtime.h>
#include <cusolverDn.h>

#ifdef CBPREAL
#define WTYPE float
#define RTYPE float
#define PARAMS float alpha =1.0, beta = 0.0;
#define CUSOLVERSVD cusolverDnSgesvd
#define CUSOLVERSVD_bufferSize cusolverDnSgesvd_bufferSize
#endif

#ifdef CBPCMPLX
#define WTYPE cuComplex
#define RTYPE float
#define PARAMS cuComplex alpha{1.0,0.0}, beta {0.0,0.0};
#define CUSOLVERSVD cusolverDnCgesvd
#define CUSOLVERSVD_bufferSize cusolverDnCgesvd_bufferSize
#endif

#ifdef CBPZMPLX
#define WTYPE cuDoubleComplex
#define RTYPE double
#define PARAMS cuDoubleComplex alpha{1.0,0.0}, beta {0.0,0.0};
#define CUSOLVERSVD cusolverDnZgesvd
#define CUSOLVERSVD_bufferSize cusolverDnZgesvd_bufferSize
#endif

#ifdef CBPDOUBLE
#define WTYPE double
#define RTYPE double
#define PARAMS double alpha =1.0, beta = 0.0;
#define CUSOLVERSVD cusolverDnDgesvd
#define CUSOLVERSVD_bufferSize cusolverDnDgesvd_bufferSize
#endif

#include "gemm.cpp"

// notice that the array becomes a line, saved row1-row2...
//integer*4-int
//cmplx-real and real
//zmplx-double float and double float
//real-float
//double-double




extern "C" {


  static const char *_cudaGetErrorEnum(cublasStatus_t error)
    {
    switch (error)
    {
      case CUBLAS_STATUS_SUCCESS:
          return "CUBLAS_STATUS_SUCCESS";

      case CUBLAS_STATUS_NOT_INITIALIZED:
          return "CUBLAS_STATUS_NOT_INITIALIZED";

      case CUBLAS_STATUS_ALLOC_FAILED:
          return "CUBLAS_STATUS_ALLOC_FAILED";

      case CUBLAS_STATUS_INVALID_VALUE:
          return "CUBLAS_STATUS_INVALID_VALUE";

      case CUBLAS_STATUS_ARCH_MISMATCH:
          return "CUBLAS_STATUS_ARCH_MISMATCH";

      case CUBLAS_STATUS_MAPPING_ERROR:
          return "CUBLAS_STATUS_MAPPING_ERROR";

      case CUBLAS_STATUS_EXECUTION_FAILED:
          return "CUBLAS_STATUS_EXECUTION_FAILED";

      case CUBLAS_STATUS_INTERNAL_ERROR:
          return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
    }

        /*
        int cbpfor_auto_svd_wrapped_ (WTYPE *A, WTYPE *s,  WTYPE* v,  WTYPE *d, int rows, int cols )
        {


		    cusolverDnHandle_t handle;
		    cusolverDnCreate(&handle);

		    int lwork;

		    CUSOLVERSVD_bufferSize(
			handle,
			rows,
			cols,
			&lwork);

		    WTYPE *d_A;
		    cudaMalloc(&d_A, sizeof(WTYPE)*rows*cols);
		    cudaMemcpy(d_A, A, sizeof(WTYPE)*rows*cols, cudaMemcpyHostToDevice);

		    WTYPE *d_U;
		    cudaMalloc(&d_U, sizeof(WTYPE)*rows*rows);

		    WTYPE *d_VT;
		    cudaMalloc(&d_VT, sizeof(WTYPE)*rows*rows);

		    WTYPE *d_work;
		    cudaMalloc(&d_work, sizeof(WTYPE)*lwork);

		    RTYPE *d_S;
		    cudaMalloc(&d_S, sizeof(RTYPE)*rows);

		    RTYPE *d_rwork;
		    cudaMalloc(&d_rwork, sizeof(RTYPE)*(rows - 1));

		    int *devInfo;
		    cudaMalloc(&devInfo, sizeof(int));

		    signed char jobu = 'A';
		    signed char jobvt = 'A';

		    CUSOLVERSVD(
			    handle,
			    jobu,
			    jobvt,
			    rows,
			    cols,
			    d_A,
			    rows,
			    d_S,
			    d_U,
			    rows,
			    d_VT,
			    rows,
			    d_work,
			    lwork,
			    d_rwork,
			    devInfo);


		    cudaFree(d_A);
		    cudaFree(d_rwork);
		    cudaFree(d_S);
		    cudaFree(d_U);
		    cudaFree(d_VT);
		    cudaFree(d_work);

		    return 0;

	}


        */


        int cbpfor_auto_gemm_wrapped_( WTYPE* a, WTYPE* b, WTYPE* c, int* dimA, int* dimB, int* dimC,int*add, int* ngpu, float* cpuratio, void** xt_)
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


                return cbp_auto_gemm(a,b,c,dimA,dimB,dimC,alpha,beta,*cpuratio, (cublasXtHandle_t) *xt_);

        }

	int init_xt_context_wrapped_(int *devcount, void** xt_)
	{

	  if(cublasXtCreate( (cublasXtHandle_t*) xt_) != CUBLAS_STATUS_SUCCESS) {printf("handle create fail\n"); return 1;}

	  if(*devcount == 0)
	  {
 		   cudaGetDeviceCount ( devcount ) ;
	  }

    int* devices = (int*) malloc((*devcount)*sizeof(int));  // add this line

	  for (int i=0; i< (*devcount) ; i++)
	  {
		    devices[i] = i ;
	  }

      int code = cublasXtDeviceSelect((cublasXtHandle_t) *xt_, *devcount, devices);

     if( code != CUBLAS_STATUS_SUCCESS) {printf("-- Set devices failed %s %d: errorcode = %s , devcount : %d \n", __FILE__ , __LINE__,
                                                    _cudaGetErrorEnum( (cublasStatus_t) code), *devcount); return code;}



	 /*cudaDeviceProp deviceProp;
	 for (int i = 0; i < *devcount; i++)
	 {
		cudaGetDeviceProperties(&deviceProp, devices[i]);
    		printf("GPU ID = %d, Name = %s \n", devices[i], deviceProp.name);
	 }*/


	  free(devices);
	  cublasXtSetPinningMemMode( (cublasXtHandle_t) *xt_, CUBLASXT_PINNING_ENABLED );


  	 return 0;



	}

	int destroy_xt_context_wrapped_(void ** xt_)
	{

	return cublasXtDestroy((cublasXtHandle_t) *xt_);

	}

}
