//***** this motherfucker will take one vector, normalized afformentioned bastard, and spit out the polarization of that filthy whore. *****

#include <stdlib.h>
#include <math.h>


void normalize_polarize(double * function_1, double left_endpoint, double right_endpoint, int num_points, double * integral){
	
	int ii; 

	double * psi_hat;
	psi_hat = (double*)malloc(num_points*sizeof(double));

	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] * function_1[ii];
	}

	trapz(psi_hat, left_endpoint, right_endpoint, num_points, integral);

	double norm = pow(*integral,0.5); 

	for(ii=0;ii<num_points;ii++){
		psi_hat[ii] = function_1[ii] / norm;
	}

	polarization(&psi_hat[0], left_endpoint, right_endpoint, num_points, integral);


}



