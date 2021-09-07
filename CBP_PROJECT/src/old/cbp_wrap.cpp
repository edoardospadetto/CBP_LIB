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

	 

  	 if(cublasXtDeviceSelect( (cublasXtHandle_t) *xt_, *devcount, devices) != CUBLAS_STATUS_SUCCESS) {printf("set devices fail\n"); return 1;} 
	 
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

