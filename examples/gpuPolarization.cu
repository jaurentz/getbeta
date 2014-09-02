#include <omp.h>
#include <hamop.h>
#include <cucheb.h>
extern "C" {
    #include "getbeta.h"
}


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
	HO.setNx(150);
	HO.setNy(150);
	HO.setNz(150);
	HO.setXminXmax(-10.0,10.0);
	HO.setYminYmax(-10.0,10.0);
	HO.setZminZmax(-10.0,10.0);
	HO.setXpotYpotZpot(1.0,2.0,3.0);
	
	// lanczos handle
	cuchebLanczosHandle LH;
	LH.n = HO.getNx()*HO.getNy()*HO.getNz();
	LH.numeigs = 2;
	LH.runlength = 60;
	LH.restarts = 4;
	
	// set n
	printf("\nn = %d\n",LH.n);
	
	// allocate memory for eigenvectors
	double *eigvecs;
	cuchebCheckError(cudaMalloc(&eigvecs,LH.numeigs*LH.n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(LH.n,eigvecs,1,1.0),__FILE__,__LINE__);

	// allocate memory for residuals
	double *res, *ray;
	cuchebCheckError((void*)(ray = (double*)malloc(LH.numeigs*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((void*)(res = (double*)malloc(LH.numeigs*sizeof(double))),__FILE__,__LINE__);
	
	// chebpoly memory
	double tau;
	double temp;
	double tol = 1e-2;
	double a, b;
	a = 0.0;
	ChebPoly SpecTrans;
	
	// chebop memory
	ChebOp SpecOp;

	// timer variables
	double begin, end;
	
	// specrad variable
	double specrad;
	
	// electric field points
	int numpoints = pow(2,3)+1;
	double chebpoints[numpoints];
	cheb_array(numpoints, &chebpoints[0], 0.1, -0.1);
	
	// polarization array
	double polarvals[numpoints];
	
	// loop to compute polarization
	for(int jj=0;jj<numpoints;jj++){
	    // set e_field strength
	    HO.setXfieldYfieldZfield(chebpoints[jj],0.0,0.0);
	    //HO.print();
	
	    // begin timer
	    begin = omp_get_wtime();
	
	    // compute specrad
	    cuchebCheckError(cuchebDspecrad(LH.n,op,(void*)&HO,&specrad),__FILE__,__LINE__);
	    b = specrad;
	
	    // end timer
	    //end = omp_get_wtime();
	    //printf("\nTime to compute specrad: %f (secs)\n",end-begin);
	    //printf("specrad = %e\n",specrad);
	
	    // chebpoly
	    tau = 1e5;
	    temp = tau/specrad/specrad;
	    SpecTrans = ChebPoly(&spectrans,(void*)&temp,&a,&b,&tol);
		
	    // set chebop
	    SpecOp = ChebOp(LH.n,op,(void*)&HO,&SpecTrans);
	    //SpecOp.print();
	
	    // begin timer
	    //begin = omp_get_wtime();
	
	    // call IRLM
	    LH.tol = 1e-14;
	    LH.numconv = 0;
	    LH.numrestarts = 0;
	    LH.nummatvecs = 0;
	    if(jj > 0){LH.runlength = 10;}	
	    cuchebCheckError(cuchebDeigs(&LH, &SpecOp, eigvecs),__FILE__,__LINE__);
	
	    // end timer
	    end = omp_get_wtime();
	    printf("\nTime to run IRLM: %f (secs)\n\n",end-begin);
	
	    // print Lanczos stats
	    printf("\nnumconv = %d\n",LH.numconv);
	    printf("numrestarts = %d\n",LH.numrestarts);
	    //printf("nummatvecs = %d\n\n",LH.nummatvecs);	
	
	    // compute rayleigh quotients and residuals
	    cuchebCheckError(rayleigh(LH.n,op,(void*)&HO,LH.numeigs,ray,eigvecs,res),__FILE__,__LINE__);
		
	    // print rayleigh quotients
	    printf("\n");
	    for(int ii=0;ii<LH.numeigs;ii++){		
		    printf("rayleigh[%d] = %+1.15e, res[%d] = %+1.2e\n",ii,ray[ii],ii,res[ii]/specrad);
	    }
	    printf("\n");
	
	    // compute polarization
	    HO.Polarization(1.0, 0.0, 0.0, eigvecs, &polarvals[jj]);
	    	
	}
	
    // print polarvals
    printf("\n");
    for(int ii=0;ii<numpoints;ii++){		
	    printf("(%+1.15e,%+1.15e)\n",chebpoints[ii],polarvals[ii]);
    }
    printf("\n");
    
    // compute beta
    double beta_abs;
    beta_abs = beta(numpoints, &chebpoints[0], &polarvals[0]);
    printf("beta_abs = %+1.15e\n\n",beta_abs);
	
	//free memory
	cuchebCheckError(cudaFree(eigvecs),__FILE__,__LINE__);
	free(ray);
	free(res);

	// return
	return 0;

}

