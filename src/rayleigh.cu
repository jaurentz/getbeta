#include <cucheb.h>

cuchebStatus_t rayleigh(int n, cuchebOpMult OPMULT, void *USERDATA, int numeigs, double *rayleigh, double *eigvecs, double *residuals){

	// check n
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check numeigs
	if(numeigs < 1){
		fprintf(stderr,"\nIn %s line: %d, numeigs must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	if(numeigs > n){
		fprintf(stderr,"\nIn %s line: %d, numeigs must be <= n.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// allocate memory for ritz vectors
	double *ritzvecs;
	cuchebCheckError(cudaMalloc(&ritzvecs,numeigs*n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// rayleigh quotients
	double temp;
	int minind;	
	for(int ii=0;ii<numeigs;ii++){
		// normalize eigvecs
		cuchebCheckError(cublasDnrm2(cublas_handle,n,&eigvecs[ii*n],1,&temp),__FILE__,__LINE__);
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,&eigvecs[ii*n],1),__FILE__,__LINE__);
		OPMULT(&eigvecs[ii*n],&ritzvecs[ii*n],USERDATA);
			
		// rayleigh quotient
		cuchebCheckError(cublasDdot(cublas_handle,n,&eigvecs[ii*n],1,&ritzvecs[ii*n],1,&rayleigh[ii]),__FILE__,__LINE__);
	}
	
	// sort in ascending order
	for(int ii=0;ii<numeigs-1;ii++){
		// index smallest rayleigh quotient
		minind = ii;
		for(int jj=ii+1;jj<numeigs;jj++){
			if(rayleigh[jj] < rayleigh[minind]){minind = jj;}
		}
		
		// swap rayleigh quotients
		temp = rayleigh[ii];
		rayleigh[ii] = rayleigh[minind];
		rayleigh[minind] = temp;
		
		// swap eigenvectors
		cuchebCheckError(cublasDswap(cublas_handle,n,&eigvecs[ii*n],1,&eigvecs[minind*n],1),__FILE__,__LINE__);
		
		// swap ritz vectors
		cuchebCheckError(cublasDswap(cublas_handle,n,&ritzvecs[ii*n],1,&ritzvecs[minind*n],1),__FILE__,__LINE__);
	}
	
	
	// residuals
	for(int ii=0;ii<numeigs;ii++){
		
		// residual
		temp = -rayleigh[ii];
		cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&eigvecs[ii*n],1,&ritzvecs[ii*n],1),__FILE__,__LINE__);
		cuchebCheckError(cublasDnrm2(cublas_handle,n,&ritzvecs[ii*n],1,&residuals[ii]),__FILE__,__LINE__);
	}
			
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);
	
	// free memory
	cuchebCheckError(cudaFree(ritzvecs),__FILE__,__LINE__);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}

