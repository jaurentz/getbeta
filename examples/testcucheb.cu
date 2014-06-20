#include <stdio.h>
#include <stdlib.h>
#include <cucheb.h>
#include <omp.h>


int main(void){

	// compute variables
	int n;
	double *x, *y, *dx, *dy;
	
	// set size
	n = 10;
	
	// allocate host memory
	x = new double[n];
	y = new double[n];
	
	// print values in host memory
	printf("\nx and y before being set:\n\n");
	for(int ii=0;ii<n;ii++){
		printf("x[%d] = %+1.15e, y[%d] = %+1.15e\n",ii,x[ii],ii,y[ii]);
	}
	printf("\n");
	
	// allocate device memory
	cuchebCheckError(cudaMalloc(&dx,n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dy,n*sizeof(double)),__FILE__,__LINE__);
	
	// set device memory using CUCHEB
	cuchebCheckError(cuchebDinit(n,dx,1,1.0),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(n,dy,1,0.0),__FILE__,__LINE__);

	// create cublas handle
	cublasHandle_t handle;
	cuchebCheckError(cublasCreate(&handle),__FILE__,__LINE__);
	
	// set pointer mode to host
	cuchebCheckError(cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);

	// compute norm of dx
	double nrmx;
	cuchebCheckError(cublasDnrm2(handle, n, dx, 1, &nrmx),__FILE__,__LINE__);

	// print norm of x
	printf("\nnorm of x = %+1.15e\n\n",nrmx);

	// compute norm of dy
	double nrmy;
	cuchebCheckError(cublasDnrm2(handle, n, dy, 1, &nrmy),__FILE__,__LINE__);

	// print norm of y
	printf("\nnorm of y = %+1.15e\n\n",nrmy);

	// destroy cublas handle
	cuchebCheckError(cublasDestroy(handle),__FILE__,__LINE__);
	
	// copy device memory to host
	cuchebCheckError(cudaMemcpy(x,dx,n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(y,dy,n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print values in host memory
	printf("\nx and y after being set:\n\n");
	for(int ii=0;ii<n;ii++){
		printf("x[%d] = %+1.15e, y[%d] = %+1.15e\n",ii,x[ii],ii,y[ii]);
	}
	printf("\n");
	
	// free device memory
	cuchebCheckError(cudaFree(dx),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dy),__FILE__,__LINE__);
	
	// free host memory
	delete[] x;
	delete[] y;
	
	// return
	return 0;

}
