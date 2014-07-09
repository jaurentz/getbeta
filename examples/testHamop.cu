#include <omp.h>
#include <getbeta.h>
#include <hamop.h>
#include <cucheb.h>

// spectrans
cuchebStatus_t spectrans(int n, const double *in, int incin, double *out, int incout, void* userdata){

	// set shape parameter
	double *tau = (double*)userdata;
	
	// compute function values
	for(int ii=0;ii<n;ii++){
		out[ii*incout] = exp(-(*tau)*in[ii*incin]*in[ii*incin]);
	}
	
	// return
	return CUCHEB_STATUS_SUCCESS;
}

// hamiltonian
void op(void* x, void* y, void* userdata){
	Hamop* HO = (Hamop*)userdata;
	HO->Mult((double*)x, (double*)y);
}


// driver
int main(){

	// declare hamiltonian
	Hamop HO;
	HO.setDims(3);
	HO.setNx(pow(2,7));
	HO.setNy(pow(2,7));
	HO.setNz(pow(2,7));
	HO.setXminXmax(0.0,1.0);
	HO.setYminYmax(0.0,1.0);
	HO.setZminZmax(0.0,1.0);
	HO.setXpotYpotZpot(10.0,20.0,30.0);
	
	// lanczos handle
	cuchebLanczosHandle LH;
	LH.n = HO.getNx()*HO.getNy()*HO.getNz();
	LH.numeigs = 4;
	LH.runlength = 60;
	LH.restarts = 10;
	LH.tol = 1e-14;
	LH.numconv = 0;
	LH.numrestarts = 0;
	LH.nummatvecs = 0;	
	
	// set n
	printf("\nn = %d\n",LH.n);

	// begin timer
	double begin, end;
	begin = omp_get_wtime();
	
	// compute specrad
	double specrad;
	cuchebCheckError(cuchebDspecrad(LH.n,op,(void*)&HO,&specrad),__FILE__,__LINE__);
	
	// end timer
	end = omp_get_wtime();
	printf("\nTime to compute specrad: %f (secs)\n",end-begin);
	printf("specrad = %e\n",specrad);
	
	// chebpoly memory
	double tau;
	double temp;
	double tol = 1e-2;
	double a, b;
	a = 0.0;
	b = specrad;
	ChebPoly SpecTrans;
	
	// chebop memory
	ChebOp SpecOp;
	
	// allocate memory for eigenvectors
	double *eigvecs;
	cuchebCheckError(cudaMalloc(&eigvecs,LH.numeigs*LH.n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(LH.n,eigvecs,1,1.0),__FILE__,__LINE__);
/*	
	// allocate memory for residuals
	double theta;
	double maxres;
	double *res, *ray;
	cuchebCheckError((void*)(ray = (double*)malloc(numeigs*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((void*)(res = (double*)malloc(numeigs*sizeof(double))),__FILE__,__LINE__);
*/	
	// chebpoly
	tau = 1e5;
	temp = tau/specrad/specrad;
	SpecTrans = ChebPoly(&spectrans,(void*)&temp,&a,&b,&tol);
		
	// set chebop
	SpecOp = ChebOp(LH.n,op,(void*)&HO,&SpecTrans);
	SpecOp.print();
	
	// begin timer
	begin = omp_get_wtime();
	
	// call IRLM
	cuchebCheckError(cuchebDeigs(&LH, &SpecOp, eigvecs),__FILE__,__LINE__);
	
	// end timer
	end = omp_get_wtime();
	printf("\nTime to run IRLM: %f (secs)\n\n",end-begin);
	
	// print Lanczos stats
	printf("\nnumconv = %d\n",LH.numconv);
	printf("numrestarts = %d\n",LH.numrestarts);
	printf("nummatvecs = %d\n\n",LH.nummatvecs);	
/*	
	// compute rayleigh quotients and residuals
	cuchebCheckError(rayleigh(n,op,(void*)&HO,numeigs,ray,eigvecs,res),__FILE__,__LINE__);
		
	// print rayleigh quotients
	printf("\n");
	for(int ii=0;ii<numeigs;ii++){		
		printf("rayleigh[%d] = %+1.2e, res[%d] = %+1.2e\n",ii,ray[ii],ii,res[ii]/specrad);
	}
	printf("\n");
		
	// compute diff
	//maxres = 0.0;
	//for(int ii=0;ii<numeigs;ii++){
	//	theta = M_PI*(ii+1.0)/2.0/(n+1.0);
	//	res[ii] = abs(4.0*sin(theta)*sin(theta) - ray[ii])/abs(4.0*sin(theta)*sin(theta));		
	//	if(res[ii] > maxres){maxres = res[ii];}
	//}

	// print results
	printf("tau = %+1.2e, degree = %d, maxres = %+1.2e, time = %f\n",tau,SpecTrans.getDegree(),maxres,end-begin);	
*/	
	//free memory
	cuchebCheckError(cudaFree(eigvecs),__FILE__,__LINE__);
//	free(ray);
//	free(res);

	// return
	return 0;

}

