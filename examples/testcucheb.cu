#include <stdio.h>
#include <stdlib.h>
#include <cucheb.h>
#include <omp.h>


int main(void){

	//defining variables
	int n;
	double *x, *y, *dx, *dy, *A;
	
	// set size
	n = 5;
	
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

	cuchebCheckError(cudaMalloc(&A,n*n*sizeof(double)),__FILE__,__LINE__); 		


	// set device memory using CUCHEB
	cuchebCheckError(cuchebDinit(n,dx,1,1.0),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(n,dy,1,1.0),__FILE__,__LINE__);

	cuchebCheckError(cuchebDinit(n*n,A,1,cuRAND),__FILE__,__LINE__);
	//cuchebCheckError(cuchebDinit(n,A,n+1,1.0),__FILE__,__LINE__);
	

	// create cublas handle
	cublasHandle_t handle;
	cuchebCheckError(cublasCreate(&handle),__FILE__,__LINE__);
	
	// set pointer mode to host
	cuchebCheckError(cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);

	// compute norm of dx
	//double nrmx;
	//cuchebCheckError(cublasDnrm2(handle, n, dx, 1, &nrmx),__FILE__,__LINE__);

	// compute sum of dx
	//double sum;
	//cuchebCheckError(cublasDasum(handle, n, dx, 1, &sum),__FILE__,__LINE__);

	// print sum of x
	//printf("\nthe sum of x = %E \n\n", sum);


	//computing alpha times a vector(dx) and storing it in (dy)
	//double alpha = 2.0;
	//int incx = 1; 
	//cuchebCheckError( cublasDaxpy(handle, n, &alpha, dx, incx, dy, incx),__FILE__,__LINE__);


	//compute the dot product between the two vectors 
	//double z;
	//int incx = 1;
	//cuchebCheckError(cublasDdot (handle, n, dx, incx, dy, incx, &z),__FILE__,__LINE__);

	//print the dot product "result"
	//printf("\n this is the dot product between x and y: \n %E", z); 


	//function performs matrix-vector multiplication y = alpha * OP(A) * x + beta * y
	int stride = 1;	
	//int n = 10;	
	double alpha = 2.0;
	double beta = 3.0;
	cuchebCheckError(cublasDgemv(handle, CUBLAS_OP_T, n, n, &alpha, A, n, dx, stride, &beta, dy, stride),__FILE__,__LINE__);	
	
	// compute norm of dy
	//double nrmy;
	//cuchebCheckError(cublasDnrm2(handle, n, dy, 1, &nrmy),__FILE__,__LINE__);

	// print norm of y
	//printf("\nnorm of y = %+1.15e\n\n",nrmy);

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
	cuchebCheckError(cudaFree(A),__FILE__,__LINE__);
	
	// free host memory
	delete[] x;
	delete[] y;
	
	// return
	return 0;

}
