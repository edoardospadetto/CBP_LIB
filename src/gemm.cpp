#include <stdio.h>
#include "cublasXt.h"
#include <curand.h>
#include <typeinfo>
#include "gemm_interfaces.h"

/*Debug function*/

template <typename T>
void MallocUnif(T* &x, long row, long col, T val) {
  x = (T*) malloc(row*col*sizeof(T));
  for (long i = 0; i < col; ++i) {
    for (long j = 0; j < row; ++j) {
      x[i * row + j] = val;
    }
  }
}

/* Get the number of available devices and and calls cublasXtDeviceSelect*/
/* if devcount == 0 the number of devices is set to the max available*/
int cbp_SetDev(cublasXtHandle_t& xt_, int devcount)
{
	  if(devcount == 0)
	  { 
 		cudaGetDeviceCount ( &devcount ) ;
	  }
  	  int* devices = (int*) malloc(devcount*sizeof(int));  // add this line
  	  for (int i=0; i< devcount ; i++)
  	  {
  		devices[i] = i ;
  	  }
  	  if(cublasXtDeviceSelect(xt_, 1, devices) != CUBLAS_STATUS_SUCCESS) {printf("set devices fail\n"); return 1;} 
  	  free(devices);
  	  return 0;

}

void printMat(double * C, int row , int col )
{

	  for (int i = 0; i < col; ++i) {
	    for (int j = 0; j < row; ++j) {
	      printf ("%lf ", C[i *row + j]);
	    }
	    printf ("\n");
  }
}

template <typename T>
int cbp_auto_gemm(T* A, T*B , T*C, int* dimA, int* dimB, int* dimC, int ngpu, T alpha, T beta )
{

  cublasXtHandle_t xt_;
        if(cublasXtCreate(&xt_) != CUBLAS_STATUS_SUCCESS) {printf("handle create fail\n"); return 1;}
        cbp_SetDev(xt_, ngpu);
  
  cbp_Gemm(A,B,C,dimA[0],dimB[1],dimA[1],xt_,alpha,beta);

  return 0;
}


